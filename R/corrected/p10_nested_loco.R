# =============================================================================
# R/corrected/p10_nested_loco.R  —  P10: selection-honest transport numbers
#
# THE PROBLEM
# -----------
# build_best_transportable_predictors() (R/transportability_best_model.R) picks
# its predictor set by scoring candidates on leave-one-country-out folds, and
# the selected set is then REPORTED on those same folds. Concretely, inside one
# call:
#
#   screen_bivariate_loco()   ranks every candidate by mean |AUC - 0.5| over all
#                             4 held-out countries
#   select_stepwise_loco()    greedily adds whichever candidate maximises mean
#                             LOCO AUC over all 4 held-out countries
#   add_construct_pca_features()  fits PCA loadings on the pooled matrix of all
#                             4 countries
#
# Every one of those steps has seen the held-out country. The reported LOCO
# number is therefore not out-of-sample with respect to the selection, and it is
# optimistic by an amount nobody has measured for this predictor set.
#
# This repository has already measured that bias once, for a different method:
# sandbox_parsimony/R/34_rescore_spatial_plus_soil.R re-scored the eight locked
# SoilGrids variables of `spatial_plus_soil` by moving selection inside the fold,
# and reported `locked_published` against `selected_in_fold`. This file does the
# same thing for the stepwise predictor set, and adopts that variant structure.
#
# THE FIX
# -------
# Nest the selection inside the outer fold. For each held-out country H:
#
#   1. Restrict BOTH the merged data and the country configs to the training
#      countries. Everything downstream is then a pure function of the training
#      set; there is no code path by which H can be consulted.
#   2. Run the ENTIRE unmodified selection procedure on that subset. Its internal
#      LOCO folds are now leave-one-of-the-training-countries-out.
#   3. Score once on H.
#
# Step 2 deliberately calls build_best_transportable_predictors() itself rather
# than reimplementing the screen, dedup, PCA and stepwise search. Reimplementing
# them would risk the nested and original schemes differing for reasons other
# than the nesting, which is the one thing this comparison must not confound.
#
# SCORING
# -------
# Both schemes are scored by the SAME estimator, .nl_fit_score(), a pooled glm
# on the selected predictors standardized on the training fold only. A glm is
# used rather than the production SuperLearner because the question here is the
# size of the selection bias, not the absolute attainable performance, and the
# comparison is only meaningful if the two schemes are scored identically and
# cheaply enough to run the full grid. The metric definitions match
# run_loco_cv() (R/transportability.R:439) column for column so the output can
# be read against results/tables/transportability_loco_results.csv.
#
# A selected predictor can be absent from the held-out country: the nested
# scheme chooses from the 3-country common set, which is a superset of the
# 4-country common set. Those predictors are dropped at scoring time, which is
# what deployment in a genuinely new country would force, and both `n_selected`
# and `n_scored` are recorded so the reader can see it happen.
#
# INFERENCE
# ---------
# Four countries give four outer folds. Paired per-fold deltas, their mean and
# their range are reported. No p-value is computed; four folds do not support
# one.
# =============================================================================

#' Outcomes the nested comparison runs over. Matches the selection procedure's
#' own joint objective (BEST_MODEL_SHARED_OUTCOMES), because the stepwise search
#' optimises the average across these four and cannot be split.
NESTED_LOCO_OUTCOMES <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

#' Scheme identifiers written into the `scheme` column.
NESTED_LOCO_SCHEMES <- c("original_selection", "nested_selection")

# ── Fold construction ────────────────────────────────────────────────────────

#' Restrict a merged-data list and a config list to a set of countries.
#'
#' Both are restricted together. Restricting only the data would leave
#' get_country_configs() advertising a country that has no rows, and
#' build_pooled_dataset() keys its per-country handling off the config.
#'
#' @param all_merged Named list of per-country merged data.frames
#' @param all_configs get_country_configs() output
#' @param keep Character vector of country names to retain
#' @return list(merged, configs) restricted to `keep`, in `keep` order
.nl_subset <- function(all_merged, all_configs, keep) {
  keep <- intersect(keep, intersect(names(all_merged), names(all_configs)))
  if (length(keep) < 2)
    stop("Nested selection needs at least 2 training countries, got ",
         length(keep), call. = FALSE)
  list(merged = all_merged[keep], configs = all_configs[keep])
}

#' Run the full selection procedure on the training countries only.
#'
#' This is the honest analogue of build_best_transportable_predictors(). It
#' passes a restricted world to that function rather than reimplementing it, so
#' the screen, the construct dedup, the PCA augmentation and the greedy forward
#' search are byte-identical to the reported procedure; only the country set
#' they see changes.
#'
#' @param all_merged Named list of per-country merged data.frames
#' @param all_configs get_country_configs() output
#' @param train_countries Character vector of training country names
#' @param seed Integer seed, recorded and set for reproducibility
#' @param ... Passed through to build_best_transportable_predictors()
#' @return list(selected_vars, path, final_glm_auc, bivariate_screen,
#'   train_countries, n_candidates)
nested_select_features <- function(all_merged, all_configs, train_countries,
                                   seed = 12345L, ...) {
  sub <- .nl_subset(all_merged, all_configs, train_countries)
  set.seed(seed)
  sel <- build_best_transportable_predictors(sub$merged, sub$configs, ...)
  list(
    selected_vars   = sel$selected_vars,
    path            = sel$path,
    final_glm_auc   = sel$final_glm_auc,
    bivariate_screen = sel$bivariate_screen,
    train_countries = train_countries,
    n_candidates    = if (is.null(sel$bivariate_screen)) NA_integer_
                      else nrow(sel$bivariate_screen)
  )
}

# ── Scoring ──────────────────────────────────────────────────────────────────

#' Row weights that give every training country equal total weight.
#'
#' The pooled individual-level fits weight each row equally, so a country with
#' more biomarker rows dominates the fit. This is the equal-country-weight
#' sensitivity: within the training set, each country's rows are scaled to sum
#' to the same total, and the weights are then normalized to mean 1 so the
#' effective sample size is unchanged.
#'
#' @param country Character vector of country labels, one per training row
#' @return Numeric vector of weights, mean 1
.nl_country_weights <- function(country) {
  n_by <- table(country)
  w <- as.numeric(1 / n_by[country])
  w / mean(w)
}

#' Fit a pooled glm on the training countries and score one held-out country.
#'
#' Standardization uses TRAINING means and standard deviations only. Complete
#' cases are taken per fold, matching loco_glm_auc_single()
#' (R/transportability_best_model.R:106).
#'
#' Metric definitions follow run_loco_cv() (R/transportability.R:439):
#'   null_brier        Brier of predicting the held-out country's OWN prevalence,
#'                     an oracle the transported model never sees
#'   null_brier_train  Brier of predicting the TRAINING prevalence on the
#'                     held-out test, the honest transport null
#'
#' @param pooled build_pooled_dataset() output over ALL countries
#' @param held_out Character, the held-out country
#' @param vars Character vector of selected predictor names
#' @param country_weighted Logical, use equal-country weights on the training rows
#' @return One-row data.frame, or a row of NA metrics when the fold is unusable
.nl_fit_score <- function(pooled, held_out, vars, country_weighted = FALSE) {
  d <- pooled$data
  # A predictor selected on the training countries can be absent from the
  # held-out country. Dropping it here is what deployment would force.
  scored <- intersect(vars, names(d))

  na_row <- function(note) data.frame(
    held_out = held_out, n_train = NA_integer_, n_test = NA_integer_,
    prev_test = NA_real_, pred_prev = NA_real_, auc = NA_real_,
    brier = NA_real_, null_brier = NA_real_, null_brier_train = NA_real_,
    calib_intercept = NA_real_, calib_slope = NA_real_, calib_scale = "logit",
    n_selected = length(vars), n_scored = length(scored),
    note = note, stringsAsFactors = FALSE)

  if (length(scored) == 0) return(na_row("no selected predictor present in held-out country"))

  tr <- d[d$country != held_out, c("Y_binary", "country", scored), drop = FALSE]
  te <- d[d$country == held_out, c("Y_binary", "country", scored), drop = FALSE]
  tr <- tr[stats::complete.cases(tr), , drop = FALSE]
  te <- te[stats::complete.cases(te), , drop = FALSE]

  if (nrow(tr) < 30 || nrow(te) < 10 ||
      length(unique(tr$Y_binary)) < 2 || length(unique(te$Y_binary)) < 2)
    return(na_row("insufficient rows or no outcome variation in a fold"))

  for (v in scored) {
    mu <- mean(tr[[v]]); sdv <- stats::sd(tr[[v]])
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    tr[[v]] <- (tr[[v]] - mu) / sdv
    te[[v]] <- (te[[v]] - mu) / sdv
  }

  w <- if (isTRUE(country_weighted)) .nl_country_weights(tr$country) else NULL

  fit <- tryCatch(
    suppressWarnings(stats::glm(
      stats::as.formula(paste("Y_binary ~", paste(scored, collapse = " + "))),
      data = tr, family = stats::binomial(), weights = w)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row("glm did not converge"))

  p <- tryCatch(as.numeric(stats::predict(fit, newdata = te, type = "response")),
                error = function(e) NULL)
  if (is.null(p) || anyNA(p)) return(na_row("prediction failed or produced NA"))

  y <- te$Y_binary
  auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))),
                  error = function(e) NA_real_)
  calib <- tryCatch({
    pc <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    stats::coef(suppressWarnings(
      stats::glm(y ~ stats::qlogis(pc), family = stats::binomial())))
  }, error = function(e) c(NA_real_, NA_real_))

  data.frame(
    held_out         = held_out,
    n_train          = nrow(tr),
    n_test           = nrow(te),
    prev_test        = mean(y),
    pred_prev        = mean(p),
    auc              = auc,
    brier            = mean((y - p)^2),
    null_brier       = mean(y) * (1 - mean(y)),
    null_brier_train = mean((mean(tr$Y_binary) - y)^2),
    calib_intercept  = unname(calib[1]),
    calib_slope      = unname(calib[2]),
    calib_scale      = "logit",
    n_selected       = length(vars),
    n_scored         = length(scored),
    note             = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Fit the production SuperLearner on the training countries and score one fold.
#'
#' The glm scorer above is the workhorse: it is cheap enough to run the full
#' grid and it matches the metric the selection search itself optimises. This
#' function is the comparability path. The published transport numbers were
#' produced by run_loco_cv() with the mlr3 SuperLearner, so a claim about how
#' those numbers change under nested selection has to be scored the same way.
#'
#' Mirrors one iteration of run_loco_cv()'s fold loop (R/transportability.R:290)
#' rather than calling it, because run_loco_cv() refits every country and this
#' needs exactly one held-out country per selected variable set.
#'
#' BART variants are dropped from the library for the same reason run_loco_cv()
#' drops them (R/transportability.R:278): dbarts MCMC has crashed the R process
#' during pooled LOCO fits.
#'
#' @param pooled build_pooled_dataset() output over ALL countries
#' @param held_out Character, the held-out country
#' @param vars Character vector of selected predictor names
#' @param sl_learners setup_mlr3_learners() output
#' @param params get_pipeline_params() output
#' @return One-row data.frame with the same columns as .nl_fit_score()
.nl_fit_score_sl <- function(pooled, held_out, vars, sl_learners, params) {
  d <- pooled$data
  scored <- intersect(vars, names(d))

  na_row <- function(note) data.frame(
    held_out = held_out, n_train = NA_integer_, n_test = NA_integer_,
    prev_test = NA_real_, pred_prev = NA_real_, auc = NA_real_,
    brier = NA_real_, null_brier = NA_real_, null_brier_train = NA_real_,
    calib_intercept = NA_real_, calib_slope = NA_real_, calib_scale = "logit",
    n_selected = length(vars), n_scored = length(scored),
    note = note, stringsAsFactors = FALSE)

  if (length(scored) == 0) return(na_row("no selected predictor present in held-out country"))

  tr <- d[d$country != held_out, , drop = FALSE]
  te <- d[d$country == held_out, , drop = FALSE]
  if (nrow(tr) < 30 || nrow(te) < 10) return(na_row("insufficient rows in a fold"))

  lib <- Filter(function(x) {
    type <- if (is.list(x)) (x[[1]] %||% "") else x
    !grepl("bart", type, ignore.case = TRUE)
  }, sl_learners$library)

  # The full stack's `ranger_low_mtry` hard-codes mtry = 8 (see
  # setup_mlr3_learners(), R/sensitivity/mlr3_fitting.R). ranger raises "mtry
  # can not be larger than number of variables in data" and exits, which takes
  # the whole SuperLearner fit down rather than just that learner. A nested
  # selection routinely returns fewer than 8 predictors, so without this filter
  # 28 of 32 cells return no metrics at all.
  #
  # Dropped rather than clamped: silently rewriting mtry would change what the
  # learner is while keeping its id, and the point of this scorer is to be the
  # published estimator. Filtering mirrors how run_loco_cv() drops BART
  # (R/transportability.R:199).
  #
  # mlr3_SL_clustered() now applies the same rule itself, via
  # sl_filter_learners_for_p(), so this guard is belt-and-braces. Kept because
  # it filters on the NOMINAL predictor count while the one inside filters on
  # the post-recipe count (imputation indicators can push the latter back to 8
  # or above); dropping it would move WS1 numbers mid-comparison.
  n_p <- length(scored)
  dropped_mtry <- Filter(function(x)
    is.list(x) && !is.null(x[["mtry"]]) && x[["mtry"]] > n_p, lib)
  if (length(dropped_mtry)) {
    lib <- Filter(function(x)
      !(is.list(x) && !is.null(x[["mtry"]]) && x[["mtry"]] > n_p), lib)
    cat(sprintf("    [nested_loco] p=%d: dropped %d learner(s) whose mtry exceeds p: %s\n",
                n_p, length(dropped_mtry),
                paste(vapply(dropped_mtry, function(x) x[["id"]] %||% x[[1]],
                             character(1)), collapse = ", ")))
  }

  fit <- tryCatch(
    mlr3_SL_clustered(d = tr, Xvars = scored, outcome = "Y_binary",
                      population = "pooled", id = "pooled_cluster",
                      folds = params$K, mlr3_library = lib,
                      outcome_type = "binomial", prescreen = TRUE),
    error = function(e) { message("  SL fit failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(na_row("SuperLearner fit failed"))

  p <- tryCatch(predict_on_new_data(fit, te, scored), error = function(e) NULL)
  if (is.null(p) || length(p) != nrow(te)) return(na_row("SuperLearner prediction failed"))

  y <- te$Y_binary
  ok <- is.finite(p) & is.finite(y)
  if (sum(ok) < 10 || length(unique(y[ok])) < 2)
    return(na_row("insufficient usable predictions"))
  p <- p[ok]; y <- y[ok]

  auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))),
                  error = function(e) NA_real_)
  calib <- tryCatch({
    pc <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    stats::coef(suppressWarnings(
      stats::glm(y ~ stats::qlogis(pc), family = stats::binomial())))
  }, error = function(e) c(NA_real_, NA_real_))

  data.frame(
    held_out = held_out, n_train = nrow(tr), n_test = length(y),
    prev_test = mean(y), pred_prev = mean(p), auc = auc,
    brier = mean((y - p)^2),
    null_brier = mean(y) * (1 - mean(y)),
    null_brier_train = mean((mean(tr$Y_binary, na.rm = TRUE) - y)^2),
    calib_intercept = unname(calib[1]), calib_slope = unname(calib[2]),
    calib_scale = "logit",
    n_selected = length(vars), n_scored = length(scored),
    note = NA_character_,
    sl_learners_used = length(lib),
    sl_learners_dropped = if (length(dropped_mtry))
      paste(vapply(dropped_mtry, function(x) x[["id"]] %||% x[[1]], character(1)),
            collapse = "|") else NA_character_,
    stringsAsFactors = FALSE
  )
}

# ── Orchestration ────────────────────────────────────────────────────────────

#' Run the nested-versus-original LOCO comparison.
#'
#' The original selection is computed once on all countries, which is exactly
#' what the production pipeline does (`best_model_selection` in _targets.R:1416).
#' The nested selection is recomputed per outer fold.
#'
#' @param all_merged Named list of per-country merged data.frames
#' @param all_configs get_country_configs() output
#' @param outcomes Character vector of outcome tags
#' @param held_out_countries Character vector, defaults to every country
#' @param seed Integer seed
#' @param country_weighted_also Logical, additionally score every cell with
#'   equal-country weights on the training rows. Applies to the glm scorer only;
#'   mlr3_SL_clustered() has no weight role (see docs/PROJECT_STATUS_2026-08_UPDATE.md
#'   section 6), so the SuperLearner rows are equal-row weighted by construction.
#' @param scorers Character vector, any of "glm" and "sl". The glm scorer runs
#'   the full grid cheaply and matches the selection search's own metric. The
#'   "sl" scorer reproduces the estimator the published transport numbers used,
#'   and is the one a claim about those numbers has to be made on.
#' @param sl_learners setup_mlr3_learners() output, required when "sl" is requested
#' @param params get_pipeline_params() output, required when "sl" is requested
#' @param ... Passed through to build_best_transportable_predictors()
#' @return list(metrics, selections) where `metrics` is the long results table
run_nested_loco <- function(all_merged, all_configs,
                            outcomes = NESTED_LOCO_OUTCOMES,
                            held_out_countries = names(all_configs),
                            seed = 12345L,
                            country_weighted_also = TRUE,
                            scorers = "glm",
                            sl_learners = NULL, params = NULL, ...) {

  scorers <- match.arg(scorers, c("glm", "sl"), several.ok = TRUE)
  if ("sl" %in% scorers && (is.null(sl_learners) || is.null(params)))
    stop("The 'sl' scorer needs sl_learners and params.", call. = FALSE)

  countries <- names(all_configs)
  held_out_countries <- intersect(held_out_countries, countries)

  cat(sprintf("[nested_loco] %d outcomes x %d outer folds over %s\n",
              length(outcomes), length(held_out_countries),
              paste(countries, collapse = "/")))

  # Selection under the ORIGINAL scheme: all countries, once.
  cat("[nested_loco] original selection (all countries, sees every held-out fold)\n")
  orig <- nested_select_features(all_merged, all_configs, countries, seed = seed, ...)
  cat(sprintf("[nested_loco]   selected %d: %s\n",
              length(orig$selected_vars), paste(orig$selected_vars, collapse = ", ")))

  # Selection under the NESTED scheme: one per outer fold.
  nested <- list()
  for (ho in held_out_countries) {
    cat(sprintf("[nested_loco] nested selection, holding out %s\n", ho))
    nested[[ho]] <- nested_select_features(
      all_merged, all_configs, setdiff(countries, ho), seed = seed, ...)
    cat(sprintf("[nested_loco]   selected %d: %s\n",
                length(nested[[ho]]$selected_vars),
                paste(nested[[ho]]$selected_vars, collapse = ", ")))
  }

  # Score. The pooled dataset is built once per outcome over ALL countries so
  # both schemes are scored on identical test rows.
  weight_modes <- if (isTRUE(country_weighted_also)) c(FALSE, TRUE) else FALSE
  rows <- list()
  selections <- list()

  for (otag in outcomes) {
    pooled <- tryCatch(build_pooled_dataset(all_merged, all_configs, otag),
                       error = function(e) NULL)
    if (is.null(pooled) || is.null(pooled$data) || !nrow(pooled$data)) {
      cat(sprintf("[nested_loco] %s: no pooled data, skipped\n", otag))
      next
    }
    for (ho in held_out_countries) {
      for (scheme in NESTED_LOCO_SCHEMES) {
        vars <- if (scheme == "original_selection") orig$selected_vars
                else nested[[ho]]$selected_vars

        if ("glm" %in% scorers) for (cw in weight_modes) {
          r <- .nl_fit_score(pooled, ho, vars, country_weighted = cw)
          r$outcome   <- otag
          r$scheme    <- scheme
          r$scorer    <- "glm"
          r$weighting <- if (cw) "equal_country" else "equal_row"
          r$seed      <- seed
          r$selected_vars <- paste(vars, collapse = "|")
          rows[[length(rows) + 1L]] <- r
        }

        if ("sl" %in% scorers) {
          cat(sprintf("[nested_loco]   SL fit: %s / %s / %s\n", otag, ho, scheme))
          set.seed(seed)
          r <- .nl_fit_score_sl(pooled, ho, vars, sl_learners, params)
          r$outcome   <- otag
          r$scheme    <- scheme
          r$scorer    <- "sl"
          r$weighting <- "equal_row"
          r$seed      <- seed
          r$selected_vars <- paste(vars, collapse = "|")
          rows[[length(rows) + 1L]] <- r
        }
      }
    }
    cat(sprintf("[nested_loco] %s scored\n", otag))
  }

  for (ho in held_out_countries)
    selections[[length(selections) + 1L]] <- data.frame(
      scheme = "nested_selection", held_out = ho,
      n_selected = length(nested[[ho]]$selected_vars),
      n_candidates = nested[[ho]]$n_candidates,
      final_glm_auc = nested[[ho]]$final_glm_auc,
      selected_vars = paste(nested[[ho]]$selected_vars, collapse = "|"),
      stringsAsFactors = FALSE)
  selections[[length(selections) + 1L]] <- data.frame(
    scheme = "original_selection", held_out = "ALL",
    n_selected = length(orig$selected_vars),
    n_candidates = orig$n_candidates,
    final_glm_auc = orig$final_glm_auc,
    selected_vars = paste(orig$selected_vars, collapse = "|"),
    stringsAsFactors = FALSE)

  list(metrics = if (length(rows)) dplyr::bind_rows(rows) else NULL,
       selections = dplyr::bind_rows(selections))
}

#' Pair the two schemes within each outer fold and summarise the deltas.
#'
#' Countries are the unit of inference. With four outer folds this reports the
#' paired delta per fold, its mean and its range, and nothing else. No p-value
#' is computed.
#'
#' @param metrics run_nested_loco()$metrics
#' @return list(paired, summary)
summarize_nested_loco <- function(metrics) {
  METRICS <- c("auc", "brier", "null_brier_train", "calib_slope", "n_scored")
  LOWER_BETTER <- c("brier")
  METRICS <- intersect(METRICS, names(metrics))

  paired <- metrics |>
    dplyr::select(dplyr::all_of(c("outcome", "held_out", "scorer", "weighting", "scheme",
                                  METRICS))) |>
    tidyr::pivot_longer(dplyr::all_of(METRICS), names_to = "metric",
                        values_to = "value") |>
    tidyr::pivot_wider(names_from = "scheme", values_from = "value") |>
    dplyr::filter(!is.na(.data$original_selection), !is.na(.data$nested_selection)) |>
    dplyr::mutate(
      delta        = .data$nested_selection - .data$original_selection,
      lower_better = .data$metric %in% LOWER_BETTER,
      nested_better = ifelse(.data$lower_better, .data$delta < 0, .data$delta > 0)
    )

  summary <- paired |>
    dplyr::group_by(.data$scorer, .data$weighting, .data$metric) |>
    dplyr::summarise(
      n_folds          = dplyr::n(),
      folds_nested_better = sum(.data$nested_better, na.rm = TRUE),
      mean_original    = round(mean(.data$original_selection), 4),
      mean_nested      = round(mean(.data$nested_selection), 4),
      mean_delta       = round(mean(.data$delta), 4),
      median_delta     = round(stats::median(.data$delta), 4),
      min_delta        = round(min(.data$delta), 4),
      max_delta        = round(max(.data$delta), 4),
      lower_better     = dplyr::first(.data$lower_better),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$scorer, .data$weighting, .data$metric)

  list(paired = paired, summary = summary)
}

#' Write the nested-LOCO tables.
#'
#' Additive by construction: all three files are new, and nothing under
#' results/tables/ is regenerated.
#'
#' @param res run_nested_loco() output
#' @param dir Output directory
#' @return Character vector of written paths
write_nested_loco_tables <- function(res, dir = here::here("results", "tables", "corrected")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  s <- summarize_nested_loco(res$metrics)
  paths <- c(
    file.path(dir, "loco_nested_selection.csv"),
    file.path(dir, "loco_nested_selection_paired.csv"),
    file.path(dir, "loco_nested_selection_summary.csv"),
    file.path(dir, "loco_nested_selected_vars.csv")
  )
  readr::write_csv(res$metrics,   paths[1])
  readr::write_csv(s$paired,      paths[2])
  readr::write_csv(s$summary,     paths[3])
  readr::write_csv(res$selections, paths[4])
  paths
}
