# =============================================================================
# R/transportability_best_model.R
#
# Builds a small, scientifically-curated, LOCO-optimized predictor set from
# the full common (cross-country-poolable) candidate pool, then evaluates it
# with the real production SL (mlr3_SL_clustered() via run_loco_cv()).
#
# Pipeline: dedup near-duplicate buffer/month/year/depth/stat variants of the
# same construct -> augment with PCA components for large domains (>=20 raw
# columns, via feature_engineering_constructs.R's construct_of()/
# reduce_pca_pooled()) -> rank by LOCO-native bivariate strength (NOT
# within-country strength -- a predictor that looks strong within one country
# can transport poorly, and vice versa) -> greedy forward stepwise selection,
# scored by average glm LOCO AUC across the 4 shared outcomes, stopping when
# additional predictors stop helping.
#
# Why glm for the search but mlr3 SL for final evaluation: the search step
# needs ~20x40 candidate x step combinations evaluated cheaply; a full SL fit
# per candidate would be prohibitively slow. glm on a handful of predictors is
# a fair, fast proxy for "does this candidate help." The FINAL selected set is
# then evaluated with the real learner library via run_loco_cv(), which as of
# 2026-08-27 skips washb_prescreen/step_corr below SMALL_P_SKIP_CORR_THRESHOLD
# predictors (mlr3_fitting.R) -- both are tuned for the large-p regime and
# actively discarded signal from an already-small, deliberately-chosen set.
# As of 2026-08-28 it also drops learners whose hard-coded mtry exceeds the
# predictor count (sl_filter_learners_for_p()); until then `ranger_low_mtry`'s
# mtry = 8 aborted the ENTIRE ensemble for any set smaller than that, which is
# why no loco_best_* target had ever produced metrics.
#
# Context for two specific findings from the investigation that produced this
# file (kept here since they explain design choices below):
#  - A food-price feature naively built with case-sensitive, ungraded-level
#    admin-name matching looked strong (rank #3) but was ~65-100% missing for
#    2 of 4 countries; its apparent strength was an artifact of that
#    missingness acting as a de facto country indicator, not real signal.
#    Properly matched (extract_wfp_price_deviation(), external_data.R) it
#    still contributes only modestly -- included in the candidate pool
#    honestly, not specially boosted.
#  - fsec_ (food security) columns exist in all 4 countries but no single
#    literal column name is common to all 4 (each dropped different all-NA
#    indicators) -- build_pooled_dataset() already logs this via its
#    "contributes 0 pooled predictors" warning. Real limitation, not a bug.
# =============================================================================

#' Prefixes/domains poolable across all 4 countries, for candidate screening
BEST_MODEL_SHARED_OUTCOMES <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

#' Collapse buffer-radius / month / year / depth / summary-stat variants of
#' the same underlying signal into one construct key, keeping the single
#' strongest-ranked representative. Purely a naming-based heuristic (no
#' semantic knowledge of what the columns measure) -- collapses e.g.
#' gee_trmm_Jan_10km / gee_trmm_Feb_25km / gee_trmm_2019 into one "construct"
#' so a stepwise search sees ~1 candidate per real signal instead of the 6+
#' near-duplicate copies most GEE variables come in.
#'
#' @param rank_tbl data.frame with columns `variable`, `strength` (or any
#'   column named `strength_col`) -- one row per candidate, already ranked
#' @param strength_col Name of the column to rank by (higher = keep)
#' @return data.frame, one row per distinct construct, arranged by strength_col desc
dedup_predictor_constructs <- function(rank_tbl, strength_col = "strength") {
  key <- function(v) {
    k <- v
    k <- gsub("_(10|25|50)km", "_RADIUS", k)
    k <- gsub("_(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)_", "_MONTH_", k)
    k <- gsub("_(19|20)[0-9]{2}(_|$)", "_YEAR\\2", k)
    k <- gsub("_(0_20|20_50)", "_DEPTH", k)
    k <- gsub("(mean|stdev|sd|min|max|range|cv)(_|$)", "STAT\\2", k)
    k <- gsub("(0_to_4_years_0_to_4|0_to_99_years_all_ages|15_plus_15|5_to_14_years_5_to_14)", "AGEGRP", k)
    k
  }
  rank_tbl$.construct <- vapply(rank_tbl$variable, key, character(1))
  out <- rank_tbl[order(-rank_tbl[[strength_col]]), ]
  out <- out[!duplicated(out$.construct), ]
  out$.construct <- NULL
  out
}

#' Bivariate LOCO screen: rank candidates by their OWN cross-country
#' transportability, not within-country strength (the two can and do differ —
#' see file header). For each candidate, fits a univariate logistic model
#' pooled across 3 countries and scores it on the 4th, held-out country,
#' averaged across all 4 countries x all 4 shared outcomes.
#'
#' @param pooled_by_outcome Named list (one per outcome in
#'   BEST_MODEL_SHARED_OUTCOMES) of build_pooled_dataset() output
#' @param candidates Character vector of column names to screen
#' @return data.frame: variable, mean_loco_strength (mean |AUC-0.5|), n_valid
screen_bivariate_loco <- function(pooled_by_outcome, candidates) {
  out <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    v <- candidates[i]
    per_outcome <- vapply(names(pooled_by_outcome), function(o) {
      loco_glm_auc_single(pooled_by_outcome[[o]], v)
    }, numeric(1))
    out[[i]] <- data.frame(variable = v,
                            mean_loco_strength = mean(abs(per_outcome - 0.5), na.rm = TRUE),
                            n_valid = sum(is.finite(per_outcome)))
  }
  res <- do.call(rbind, out)
  res[order(-res$mean_loco_strength), ]
}

#' LOCO-averaged glm AUC for one or more predictors (helper shared by the
#' bivariate screen and the stepwise search). Standardizes on TRAINING-fold
#' mean/sd only (no leakage), uses complete-cases per fold.
#'
#' @param pooled build_pooled_dataset() output for one outcome
#' @param vars Character vector of predictor column names (1+)
#' @return Mean AUC across the held-out countries with usable data (NA if none)
loco_glm_auc_single <- function(pooled, vars) {
  d <- pooled$data
  countries <- unique(d$country)
  aucs <- numeric(0)
  for (held in countries) {
    tr <- d[d$country != held, c("Y_binary", vars), drop = FALSE]
    te <- d[d$country == held, c("Y_binary", vars), drop = FALSE]
    tr <- tr[stats::complete.cases(tr), ]
    te <- te[stats::complete.cases(te), ]
    if (nrow(tr) < 30 || nrow(te) < 10 ||
        length(unique(tr$Y_binary)) < 2 || length(unique(te$Y_binary)) < 2) next
    for (v in vars) {
      mu <- mean(tr[[v]]); sdv <- stats::sd(tr[[v]])
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      tr[[v]] <- (tr[[v]] - mu) / sdv
      te[[v]] <- (te[[v]] - mu) / sdv
    }
    fit <- tryCatch(
      suppressWarnings(stats::glm(stats::as.formula(paste("Y_binary ~", paste(vars, collapse = " + "))),
                                  data = tr, family = stats::binomial())),
      error = function(e) NULL)
    if (is.null(fit)) next
    pred <- tryCatch(as.numeric(stats::predict(fit, newdata = te, type = "response")), error = function(e) NULL)
    if (is.null(pred) || anyNA(pred)) next
    auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(te$Y_binary, pred, quiet = TRUE))), error = function(e) NA_real_)
    if (is.finite(auc)) aucs <- c(aucs, auc)
  }
  if (length(aucs) == 0) return(NA_real_)
  mean(aucs)
}

#' Add PCA components for large predictor domains (construct-level, via
#' feature_engineering_constructs.R's construct_of()/reduce_pca_pooled()) as
#' new candidate columns on each outcome's pooled$data.
#'
#' PCA loadings are fit ONCE on the full pooled (all 4 countries) common-
#' predictor matrix -- this is a mild leakage in the SELECTION/SCREENING step
#' only (candidate columns are otherwise treated as static throughout this
#' screen, same as every raw candidate). If the FINAL selected set includes a
#' PC, a production caller wanting a fully LOCO-clean number should refit
#' reduce_pca_pooled() per held-out fold via its own project() closure rather
#' than reuse these globally-fit scores.
#'
#' @param pooled_by_outcome Named list of build_pooled_dataset() output
#' @param common_cols Character vector of common (cross-country) candidate
#'   columns to search for large domains within
#' @param k PCs per construct (default 3)
#' @param min_cols_for_pca Only construct domains with at least this many raw
#'   columns get PCA-reduced (default 20, per the request that motivated this)
#' @return list(pooled_by_outcome = <same list, with PC columns added to each
#'   $data>, pc_names = character vector of the new PC column names)
add_construct_pca_features <- function(pooled_by_outcome, common_cols, k = 3L, min_cols_for_pca = 20L) {
  con_labels <- construct_of(common_cols)
  con_counts <- table(con_labels)
  big_constructs <- names(con_counts[con_counts >= min_cols_for_pca])
  if (length(big_constructs) == 0) return(list(pooled_by_outcome = pooled_by_outcome, pc_names = character(0)))

  safe_impute <- function(df) {
    for (col in colnames(df)) {
      if (anyNA(df[[col]]) || any(!is.finite(df[[col]]))) {
        med <- stats::median(df[[col]], na.rm = TRUE)
        if (!is.finite(med)) med <- 0
        bad <- is.na(df[[col]]) | !is.finite(df[[col]])
        df[[col]][bad] <- med
      }
    }
    df
  }

  first_outcome <- pooled_by_outcome[[1]]
  X_pool <- safe_impute(first_outcome$data[, common_cols, drop = FALSE])
  pca_out <- reduce_pca_pooled(X_pool, k = k, constructs = big_constructs, min_keep = 8L)
  pc_names <- colnames(as.data.frame(pca_out$features))

  for (oname in names(pooled_by_outcome)) {
    Xo <- safe_impute(pooled_by_outcome[[oname]]$data[, common_cols, drop = FALSE])
    proj <- as.data.frame(pca_out$project(Xo))
    for (pc in colnames(proj)) pooled_by_outcome[[oname]]$data[[pc]] <- proj[[pc]]
    pooled_by_outcome[[oname]]$Xvars_common <- union(pooled_by_outcome[[oname]]$Xvars_common, colnames(proj))
  }
  list(pooled_by_outcome = pooled_by_outcome, pc_names = pc_names, pca_fit = pca_out)
}

#' Greedy forward stepwise predictor selection, scored by average glm LOCO
#' AUC across BEST_MODEL_SHARED_OUTCOMES. Stops after two consecutive
#' additions that gain less than `min_gain`.
#'
#' @param pooled_by_outcome Named list of build_pooled_dataset() output
#'   (outcomes to average over = names(pooled_by_outcome))
#' @param candidates Character vector of candidate column names
#' @param max_steps Hard cap on selected predictors (default 20)
#' @param min_gain Minimum AUC gain to keep counting as "still improving" (default 0.003)
#' @return list(selected = character vector, path = data.frame of the search trace)
select_stepwise_loco <- function(pooled_by_outcome, candidates, max_steps = 20L, min_gain = 0.003) {
  selected <- character(0)
  remaining <- candidates
  path <- list()
  best_avg <- 0.5
  patience <- 0L
  for (step in seq_len(max_steps)) {
    if (length(remaining) == 0) break
    scores <- vapply(remaining, function(v) {
      trial_vars <- c(selected, v)
      per_outcome <- vapply(names(pooled_by_outcome), function(o)
        loco_glm_auc_single(pooled_by_outcome[[o]], trial_vars), numeric(1))
      mean(per_outcome, na.rm = TRUE)
    }, numeric(1))
    best_idx <- which.max(scores)
    gain <- scores[best_idx] - best_avg
    path[[step]] <- data.frame(step = step, added = remaining[best_idx],
                               avg_loco_auc = scores[best_idx], gain = gain)
    cat(sprintf("[select_stepwise_loco] step %2d: %-45s avg LOCO AUC = %.4f (gain %+.4f)\n",
                step, remaining[best_idx], scores[best_idx], gain))
    if (gain < min_gain) {
      patience <- patience + 1L
      if (patience >= 2L) { cat("[select_stepwise_loco] stopping: two non-improving steps\n"); break }
    } else {
      patience <- 0L
    }
    selected <- c(selected, remaining[best_idx])
    best_avg <- scores[best_idx]
    remaining <- setdiff(remaining, remaining[best_idx])
  }
  list(selected = selected, path = do.call(rbind, path), final_glm_auc = best_avg)
}

#' Orchestrator: build the best transportable predictor set from scratch.
#'
#' @param all_merged Named list of per-country merged data.frames
#' @param all_configs get_country_configs() output
#' @param n_top_by_strength How many top (deduplicated) raw candidates to
#'   seed the stepwise search with, before adding ruralness/access + PCA
#'   (default 40)
#' @param ruralness_access_vars Character vector of always-included
#'   ruralness/healthcare-access candidates regardless of their individual
#'   rank (default: population density, urbanization, travel-time
#'   accessibility, elevation, demographic structure -- confirmed common
#'   across all 4 countries but individually weak, so a pure top-N-by-
#'   strength cutoff would never give them a chance)
#' @return list(selected_vars, path, pooled_by_outcome (with PCA columns),
#'   bivariate_screen, pca_fit)
build_best_transportable_predictors <- function(all_merged, all_configs,
    n_top_by_strength = 40L,
    ruralness_access_vars = c("gee_smod_50km", "gee_smod_25km", "gee_smod_10km",
                              "gee_ghsl_smod_2015", "gee_elevation_2000",
                              "gee_gpw_demographic_2010_annual_mean",
                              "gee_popdensity_2010", "gee_accessibility_2019")) {

  pooled_by_outcome <- setNames(
    lapply(BEST_MODEL_SHARED_OUTCOMES, function(o) build_pooled_dataset(all_merged, all_configs, o)),
    BEST_MODEL_SHARED_OUTCOMES)

  common_cols <- pooled_by_outcome[[1]]$Xvars_common
  is_num <- vapply(common_cols, function(v) is.numeric(pooled_by_outcome[[1]]$data[[v]]), logical(1))
  common_cols <- common_cols[is_num]

  cat(sprintf("[build_best_transportable_predictors] %d numeric common candidates\n", length(common_cols)))

  biv <- screen_bivariate_loco(pooled_by_outcome, common_cols)
  dedup <- dedup_predictor_constructs(biv, "mean_loco_strength")
  top_n <- head(dedup$variable, n_top_by_strength)

  pca <- add_construct_pca_features(pooled_by_outcome, common_cols)
  pooled_by_outcome <- pca$pooled_by_outcome

  wfp_dev_vars <- grep("^wfp_dev_", common_cols, value = TRUE)
  ruralness_access_vars <- intersect(ruralness_access_vars, common_cols)

  stepwise_pool <- unique(c(top_n, ruralness_access_vars, wfp_dev_vars, pca$pc_names))
  cat(sprintf(
    "[build_best_transportable_predictors] stepwise candidate pool: %d unique vars (%d top-strength, %d ruralness/access, %d wfp_dev, %d PCA)\n",
    length(stepwise_pool), length(top_n), length(ruralness_access_vars),
    length(wfp_dev_vars), length(pca$pc_names)))

  sel <- select_stepwise_loco(pooled_by_outcome, stepwise_pool)

  list(selected_vars = sel$selected, path = sel$path, final_glm_auc = sel$final_glm_auc,
       pooled_by_outcome = pooled_by_outcome, bivariate_screen = biv, pca_fit = pca$pca_fit)
}

#' Build a pooled dataset restricted to a pre-selected predictor set
#' (mirrors build_pooled_gee_only()'s pattern in transportability.R).
#'
#' @param all_merged Named list of per-country merged data.frames
#' @param all_configs get_country_configs() output
#' @param outcome_tag Character, e.g. "child_vitA"
#' @param selected_vars Character vector of predictor names (from
#'   build_best_transportable_predictors()$selected_vars)
#' @return build_pooled_dataset() output with Xvars_common restricted to selected_vars
build_pooled_best_model <- function(all_merged, all_configs, outcome_tag, selected_vars) {
  pooled <- build_pooled_dataset(all_merged, all_configs, outcome_tag)
  if (is.null(pooled) || nrow(pooled$data) == 0) return(pooled)
  use_vars <- intersect(selected_vars, pooled$Xvars_common)
  if (length(use_vars) == 0) {
    warning("None of the selected best-model predictors are present in this outcome's pooled data")
    return(NULL)
  }
  pooled$Xvars_common <- use_vars
  pooled
}

#' Run leave-one-country-out CV using the best-model predictor set.
#' Identical to run_loco_gee_only() -- delegates entirely to run_loco_cv(),
#' which as of 2026-08-27 handles small predictor counts correctly (see
#' mlr3_fitting.R's SMALL_P_SKIP_CORR_THRESHOLD).
#'
#' @param pooled_best build_pooled_best_model() output
#' @param sl_learners setup_mlr3_learners() output
#' @param params get_pipeline_params() output
#' @return run_loco_cv() output data.frame
run_loco_best_model <- function(pooled_best, sl_learners, params) {
  if (is.null(pooled_best) || nrow(pooled_best$data) == 0) {
    warning("No pooled best-model data available for LOCO")
    return(data.frame())
  }
  cat(sprintf("\n[LOCO-BEST] Using %d selected predictors\n", length(pooled_best$Xvars_common)))
  run_loco_cv(pooled_best, sl_learners, params)
}
