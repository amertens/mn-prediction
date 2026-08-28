# =============================================================================
# R/area_level_comparison.R
#
# Area-level SuperLearner models for Admin-2 prevalence prediction.
# Compares individual-level SL (aggregated to Admin-2) with area-level SL
# using both MSE (continuous) and NLL (binomial/beta-inspired) loss functions.
#
# Also provides LOCO cross-country area-level validation.
# =============================================================================

#' Fit an area-level SuperLearner on Admin-2 survey prevalence
#'
#' Fits mlr3superlearner on Admin-2 level: Y = survey prevalence, X = GEE covariates.
#' Uses LOO-CV (appropriate for small n ~ 20-90 Admin-2 units per country).
#'
#' @param svy_admin2 Data.frame with columns Admin2, svy_prev, n_svy
#' @param gee_admin2 Data.frame with columns Admin2 + gee_* variables
#' @param outcome_type "continuous" (MSE loss) or "binomial" (NLL loss)
#' @param sl_library List of mlr3 learner specs
#' @param n_folds Number of CV folds (default: LOO via nrow)
#' @return List with area_preds, loo_preds, metrics, model_type
fit_area_sl <- function(svy_admin2, gee_admin2, outcome_type = "binomial",
                        sl_library = NULL, n_folds = NULL) {

  if (!requireNamespace("mlr3superlearner", quietly = TRUE))
    stop("mlr3superlearner required")
  library(mlr3extralearners, quietly = TRUE)

  # Default library: compact for small-n area-level data (p >> n typical)
  # Avoid BART and GP which fail with NaN/sigma errors when p >> n
  if (is.null(sl_library)) {
    sl_library <- list(
      list("glmnet", alpha = 1, id = "lasso"),
      list("glmnet", alpha = 0, id = "ridge"),
      list("glmnet", alpha = 0.5, id = "elastic_net"),
      list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
      list("xgboost", max_depth = 2, eta = 0.05, nrounds = 200,
           min_child_weight = 3, subsample = 0.8, id = "xgb")
    )
  }

  # Merge survey prevalence with GEE covariates
  train_df <- dplyr::inner_join(
    svy_admin2 |> dplyr::select(Admin2, svy_prev, dplyr::any_of("n_svy")),
    gee_admin2,
    by = "Admin2"
  ) |>
    dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))

  if (nrow(train_df) < 5) {
    warning("Too few Admin-2 areas for area-level SL (n = ", nrow(train_df), ")")
    return(list(area_preds = NULL, loo_preds = NULL,
                metrics = NULL, model_type = outcome_type))
  }

  # Identify valid GEE predictors (non-constant, non-NA)
  gee_vars <- setdiff(names(gee_admin2), c("Admin1", "Admin2"))
  valid_vars <- gee_vars[sapply(gee_vars, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  if (length(valid_vars) < 2) {
    warning("Too few valid GEE predictors for area-level SL")
    return(list(area_preds = NULL, loo_preds = NULL,
                metrics = NULL, model_type = outcome_type))
  }

  cat(sprintf("  [area_sl] n=%d Admin-2 areas, p=%d GEE vars, loss=%s\n",
              nrow(train_df), length(valid_vars), outcome_type))

  # Prepare Y and X
  Y <- train_df$svy_prev
  logit_transform <- FALSE
  if (outcome_type == "binomial") {
    # mlr3superlearner's "binomial" expects 0/1, not proportions.
    # Instead: logit-transform prevalence, fit as continuous on logit scale,
    # then back-transform predictions via plogis(). This gives logistic loss
    # on the original scale — equivalent to beta regression's link function.
    Y <- pmin(pmax(Y, 0.001), 0.999)
    Y <- qlogis(Y)  # logit transform
    outcome_type <- "continuous"  # fit on logit scale
    logit_transform <- TRUE
    cat("  [area_sl] Using logit-transformed prevalence (beta-inspired)\n")
  }

  X <- as.data.frame(train_df[, valid_vars, drop = FALSE])

  # Impute NAs with column medians
  for (col in names(X)) {
    if (inherits(X[[col]], "haven_labelled")) X[[col]] <- as.double(unclass(X[[col]]))
    na_idx <- is.na(X[[col]]) | !is.finite(X[[col]])
    if (any(na_idx)) X[[col]][na_idx] <- median(X[[col]][!na_idx], na.rm = TRUE)
  }

  # Determine folds (LOO for small n, otherwise 5-fold)
  if (is.null(n_folds)) {
    n_folds <- if (nrow(X) <= 50) nrow(X) else 5L
  }

  # Monkey-patch NLL to clip probabilities
  env <- asNamespace("mlr3superlearner")
  if (exists("loss_nll", envir = env)) {
    clipped_nll <- function(x, y) {
      x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
      -mean(y * log(x) + (1 - y) * log(1 - x))
    }
    tryCatch(assignInNamespace("loss_nll", clipped_nll, "mlr3superlearner"),
             error = function(e) NULL)
  }

  # Fit SL
  fit_df <- data.frame(Y = Y, X)
  sl_fit <- tryCatch(
    mlr3superlearner::mlr3superlearner(
      data = fit_df, target = "Y", library = sl_library,
      outcome_type = outcome_type, folds = n_folds
    ),
    error = function(e) {
      cat(sprintf("  [area_sl] SL fit failed: %s\n", e$message))
      NULL
    }
  )

  if (is.null(sl_fit)) {
    return(list(area_preds = NULL, loo_preds = NULL,
                metrics = NULL, model_type = outcome_type))
  }

  # Full-data (in-sample) predictions. Used ONLY as the all-areas point estimate
  # (area_pred_full); NEVER for performance metrics.
  preds <- tryCatch(as.numeric(predict(sl_fit, fit_df)), error = function(e) NULL)
  if (is.null(preds)) preds <- rep(NA_real_, nrow(fit_df))

  # ── Genuine out-of-fold (OOF) predictions for honest metrics ────────────────
  # 2026-06-17 fix: mlr3superlearner's internal CV only tunes the ensemble
  # WEIGHTS; predict(sl_fit, fit_df) returns IN-SAMPLE fits, which inflate
  # within-country r to ~0.97 (even 1.00). We therefore refit the whole SL in an
  # OUTER CV loop and predict each held-out fold, so loo_preds and `metrics` are
  # genuinely out-of-sample.
  set.seed(20260617L)
  n_obs   <- nrow(fit_df)
  # Outer CV for honest OOF — capped at 10 folds. LOO outer-CV would be O(n^2)
  # SL refits (and is unnecessary); 10-fold gives honest OOF cheaply. The inner
  # SL weight-tuning CV is also capped (5-fold) so the outer loop stays fast.
  n_oof   <- min(10L, n_obs)
  oof_fold <- if (n_oof >= n_obs) seq_len(n_obs)
              else sample(rep_len(seq_len(n_oof), n_obs))
  oof_raw <- rep(NA_real_, n_obs)
  for (k in unique(oof_fold)) {
    tr <- which(oof_fold != k); te <- which(oof_fold == k)
    if (length(tr) < 3) next
    f_k <- tryCatch(
      mlr3superlearner::mlr3superlearner(
        data = fit_df[tr, , drop = FALSE], target = "Y",
        library = sl_library, outcome_type = outcome_type,
        folds = min(5L, length(tr))),
      error = function(e) NULL)
    if (is.null(f_k)) next
    p_k <- tryCatch(as.numeric(predict(f_k, fit_df[te, , drop = FALSE])),
                    error = function(e) NULL)
    if (!is.null(p_k) && length(p_k) == length(te)) oof_raw[te] <- p_k
  }

  # Back-transform both to the [0,1] prevalence scale.
  if (logit_transform) {
    preds   <- plogis(preds)
    oof_raw <- plogis(oof_raw)
    Y <- train_df$svy_prev  # original-scale outcome for metrics
  }
  preds     <- pmin(pmax(preds, 0), 1)
  oof_preds <- pmin(pmax(oof_raw, 0), 1)

  # loo_preds are genuine OUT-OF-FOLD predictions (not in-sample).
  loo_preds <- oof_preds

  # Metrics computed on OOF predictions only.
  valid <- is.finite(oof_preds) & is.finite(Y)
  metrics <- if (sum(valid) >= 3) {
    data.frame(
      model_type = outcome_type,
      n = sum(valid),
      mae_pp = round(mean(abs(Y[valid] - oof_preds[valid])) * 100, 2),
      rmse_pp = round(sqrt(mean((Y[valid] - oof_preds[valid])^2)) * 100, 2),
      pearson_r = round(cor(Y[valid], oof_preds[valid]), 3),
      mean_bias_pp = round(mean(oof_preds[valid] - Y[valid]) * 100, 2),
      stringsAsFactors = FALSE
    )
  } else NULL

  # Get winner
  winner <- tryCatch({
    wts <- sl_fit$weights
    names(which.max(wts))
  }, error = function(e) "unknown")

  if (!is.null(metrics)) {
    metrics$winner <- winner
    cat(sprintf("  [area_sl] %s: MAE=%.1f pp, r=%.3f, winner=%s\n",
                outcome_type, metrics$mae_pp, metrics$pearson_r, winner))
  }

  # Build prediction set (with Admin2 names). `area_pred` is the genuine OOF
  # prediction (so compare_admin2_approaches reports honest within-country
  # metrics); `area_pred_full` keeps the in-sample full-data fit for reference.
  area_preds <- data.frame(
    Admin2 = train_df$Admin2,
    area_pred = oof_preds,
    area_pred_full = preds,
    svy_prev = train_df$svy_prev,
    stringsAsFactors = FALSE
  )

  list(
    area_preds = area_preds,
    loo_preds = loo_preds,
    metrics = metrics,
    model_type = outcome_type,
    valid_vars = valid_vars,
    sl_fit = sl_fit
  )
}


#' Compare individual-level SL, area MSE, and area NLL predictions at Admin-2
#'
#' @param sl_admin2 Data.frame from aggregate_admin2_sl() — individual SL agg to Admin-2
#' @param area_mse List from fit_area_sl(outcome_type="continuous")
#' @param area_nll List from fit_area_sl(outcome_type="binomial")
#' @param svy_admin2 Data.frame from compute_svy_admin2() — survey-weighted prevalence
#' @param country Country name
#' @param outcome Outcome tag
#' @return Data.frame with error metrics for each approach
compare_admin2_approaches <- function(sl_admin2, area_mse, area_nll,
                                       svy_admin2, country, outcome, cc = NULL) {
  rows <- list()

  # Helper to compute metrics
  compute_err <- function(pred_df, pred_col, label) {
    merged <- dplyr::inner_join(
      svy_admin2 |> dplyr::select(Admin2, svy_prev),
      pred_df |> dplyr::select(Admin2, pred = dplyr::all_of(pred_col)),
      by = "Admin2"
    ) |> dplyr::filter(!is.na(svy_prev), !is.na(pred))

    if (nrow(merged) < 3) return(NULL)

    data.frame(
      country = country,
      outcome = outcome,
      approach = label,
      n_admin2 = nrow(merged),
      mae_pp = round(mean(abs(merged$svy_prev - merged$pred)) * 100, 2),
      rmse_pp = round(sqrt(mean((merged$svy_prev - merged$pred)^2)) * 100, 2),
      # constant predictions (the null row) legitimately give NA, not a warning
      pearson_r = round(suppressWarnings(cor(merged$svy_prev, merged$pred)), 3),
      mean_bias_pp = round(mean(merged$pred - merged$svy_prev) * 100, 2),
      stringsAsFactors = FALSE
    )
  }

  # Individual-level SL aggregated to Admin-2
  if (!is.null(sl_admin2) && nrow(sl_admin2) > 0) {
    rows[["individual_sl"]] <- compute_err(sl_admin2, "sl_prev", "Individual SL")
  }

  # Area-level SL with MSE loss
  if (!is.null(area_mse) && !is.null(area_mse$area_preds)) {
    rows[["area_mse"]] <- compute_err(area_mse$area_preds, "area_pred", "Area SL (MSE)")
  }

  # Area-level SL with NLL loss
  if (!is.null(area_nll) && !is.null(area_nll$area_preds)) {
    rows[["area_nll"]] <- compute_err(area_nll$area_preds, "area_pred", "Area SL (NLL)")
  }

  # ── Benchmarks the table needs in order to be readable ───────────────────
  # Without these two rows there is no way to tell whether the GEE block is
  # contributing anything. In the sandbox bake-off a covariate-free lon/lat
  # smoother OUT-SCORED the production recipe in several country x outcome
  # cells, which only became visible once these rows existed.
  svy_ok <- svy_admin2[is.finite(svy_admin2$svy_prev), , drop = FALSE]
  if (nrow(svy_ok) >= 5) {
    # (1) national mean painted on every district -- what a programme has
    #     without any model at all
    nat <- stats::weighted.mean(svy_ok$svy_prev,
                                pmax(svy_ok$n_svy %||% rep(1, nrow(svy_ok)), 1))
    rows[["null_mean"]] <- compute_err(
      data.frame(Admin2 = svy_ok$Admin2, pred = nat, stringsAsFactors = FALSE),
      "pred", "National mean (null)")

    # (2) leave-one-out spatial smoother on district centroids, no covariates
    if (!is.null(cc) && !is.null(cc$gadm_code) && requireNamespace("mgcv", quietly = TRUE)) {
      ctr <- tryCatch(
        sf::st_drop_geometry(load_admin2_centroids(cc$gadm_code))[, c("Admin2", "lon", "lat")],
        error = function(e) NULL)
      if (!is.null(ctr)) {
        sp <- dplyr::inner_join(svy_ok, ctr, by = "Admin2")
        sp <- sp[is.finite(sp$lon) & is.finite(sp$lat), , drop = FALSE]
        if (nrow(sp) >= 20) {
          y <- stats::qlogis(pmin(pmax(sp$svy_prev, 0.005), 0.995))
          oof <- rep(NA_real_, nrow(sp))
          kk <- min(20L, max(5L, floor(nrow(sp) / 4)))
          for (i in seq_len(nrow(sp))) {
            g <- tryCatch(mgcv::gam(y[-i] ~ s(lon, lat, k = kk),
                                    data = sp[-i, , drop = FALSE],
                                    weights = pmax(sp$n_svy[-i], 1), method = "REML"),
                          error = function(e) NULL)
            if (!is.null(g))
              oof[i] <- stats::plogis(as.numeric(
                stats::predict(g, newdata = sp[i, , drop = FALSE])))
          }
          if (sum(is.finite(oof)) >= 10)
            rows[["spatial_only"]] <- compute_err(
              data.frame(Admin2 = sp$Admin2, pred = oof, stringsAsFactors = FALSE),
              "pred", "Spatial only (no covariates)")
        }
      }
    }
  }

  if (length(rows) == 0) return(NULL)
  out <- dplyr::bind_rows(rows)
  # r_max / r_share: is a Pearson r of 0.30 near the ceiling or far from it?
  if (exists("add_reliability_columns")) out <- add_reliability_columns(out, svy_admin2)

  # ── Report the LEVEL and the PATTERN as separate columns ─────────────────
  # A single RMSE hides which of the two failed. On these data a 40 pp national
  # level miss can sit alongside a perfectly usable district ranking, and the
  # benchmarking step (R/benchmark_area.R) fixes the first without touching the
  # second -- so they have to be visible separately for anyone to see that.
  nd <- attr(svy_admin2, "national_design_based")
  if (!is.null(nd) && is.finite(nd$prev)) {
    out$national_design_pp <- round(100 * nd$prev, 2)
    # national error of each approach, weighted by survey n as a population
    # proxy (a true population weight is better; see .area_population_weights)
    lvl <- vapply(out$approach, function(ap) {
      src <- switch(ap,
        "Individual SL" = if (!is.null(sl_admin2)) sl_admin2 else NULL,
        "Area SL (MSE)" = if (!is.null(area_mse)) area_mse$area_preds else NULL,
        "Area SL (NLL)" = if (!is.null(area_nll)) area_nll$area_preds else NULL,
        NULL)
      if (is.null(src)) return(NA_real_)
      col <- if ("sl_prev" %in% names(src)) "sl_prev" else "area_pred"
      m <- merge(src[, c("Admin2", col)],
                 svy_admin2[, c("Admin2", "n_svy")], by = "Admin2")
      if (!nrow(m)) return(NA_real_)
      100 * (stats::weighted.mean(m[[col]], pmax(m$n_svy, 1), na.rm = TRUE) - nd$prev)
    }, numeric(1))
    out$national_err_pp <- round(lvl, 2)
  }
  out
}


#' Build pooled area-level dataset for LOCO cross-validation
#'
#' Pools Admin-2 survey prevalences and GEE covariates across countries,
#' keeping only common GEE variables.
#'
#' @param svy_admin2_list Named list of svy_admin2 data.frames (one per country)
#' @param gee_admin2_list Named list of gee_admin2 data.frames (one per country)
#' @return List with pooled_data, common_gee_vars, country_names
build_area_loco_dataset <- function(svy_admin2_list, gee_admin2_list) {

  # Find common GEE variables across all countries
  gee_var_lists <- lapply(gee_admin2_list, function(g) {
    setdiff(names(g), c("Admin1", "Admin2"))
  })
  common_vars <- Reduce(intersect, gee_var_lists)
  cat(sprintf("[area_loco] %d common GEE variables across %d countries\n",
              length(common_vars), length(gee_admin2_list)))

  # Optional covariate hygiene: drop cross-band _annual_* summaries over
  # non-commensurable bands and static layers' per-year duplicates. Applied
  # AFTER the intersection so the pruning is identical for every country.
  # Active by default; set GEE_COVARIATE_HYGIENE=false to keep the unpruned set.
  # See R/gee_band_semantics.R and docs/findings/WS2b_hygiene_flip.md.
  hy_drop <- prune_gee_covariates(common_vars)
  if (length(hy_drop)) {
    common_vars <- setdiff(common_vars, hy_drop)
    cat(sprintf("[area_loco] %d GEE variables after covariate hygiene\n",
                length(common_vars)))
  }

  # Per-country Admin-2 centroid (lon/lat) lookup so spatial comparators
  # (spatial_coords / spatial_plus_soil) have coordinates to work with. Pulled
  # from the shared GADM cache via the same helper GWR uses; degrades
  # gracefully (lon/lat = NA) if a code/cache is missing, in which case spatial
  # methods simply no-op as before.
  configs <- tryCatch(get_country_configs(), error = function(e) NULL)
  centroids_for <- function(ctry) {
    if (is.null(configs) || is.null(configs[[ctry]]$gadm_code)) return(NULL)
    tryCatch(
      sf::st_drop_geometry(
        load_admin2_centroids(configs[[ctry]]$gadm_code)
      )[, c("Admin2", "Admin1", "lon", "lat")],
      error = function(e) NULL)
  }

  # Pool data
  pooled <- list()
  for (ctry in names(svy_admin2_list)) {
    svy <- svy_admin2_list[[ctry]]
    gee <- gee_admin2_list[[ctry]]
    if (is.null(svy) || is.null(gee) || nrow(svy) == 0) next

    # Join on (Admin1, Admin2) where both sides carry Admin1, on the name alone
    # otherwise. GADM Admin-2 names are not unique -- Malawi has four same-named
    # district pairs in different regions -- so a name-only join against a GADM
    # -derived table multiplies rows silently. See R/admin2_key_hygiene.R.
    svy_sel <- svy |> dplyr::select(dplyr::any_of("Admin1"), Admin2, svy_prev,
                                    dplyr::any_of("n_svy"))
    gee_sel <- gee |> dplyr::select(dplyr::any_of("Admin1"), Admin2,
                                    dplyr::all_of(common_vars))
    by_sg <- admin2_join_by(svy_sel, gee_sel)

    n_before <- nrow(svy_sel)
    report_pair_join_losses(svy_sel, gee_sel, by_sg,
                            sprintf("%s survey-to-GEE join", ctry))
    merged <- dplyr::inner_join(svy_sel, gee_sel, by = by_sg) |>
      dplyr::filter(!is.na(svy_prev), is.finite(svy_prev)) |>
      dplyr::mutate(country = ctry)
    merged <- warn_if_join_multiplied(
      n_before, merged, sprintf("%s survey-to-GEE join", ctry))

    ctr <- centroids_for(ctry)
    if (!is.null(ctr)) {
      by_ctr <- admin2_join_by(merged, ctr)
      # When the pair key is unavailable the centroid table must be collapsed
      # first: it comes straight from GADM and carries the duplicates.
      if (identical(by_ctr, "Admin2"))
        ctr <- dedupe_admin2_key(ctr, sprintf("%s centroids", ctry))
      n_before <- nrow(merged)
      merged <- dplyr::left_join(merged, ctr, by = by_ctr)
      merged <- warn_if_join_multiplied(
        n_before, merged, sprintf("%s centroid join", ctry))
    }

    pooled[[ctry]] <- merged
  }

  pooled_df <- dplyr::bind_rows(pooled)
  cat(sprintf("[area_loco] Pooled %d Admin-2 areas from %d countries\n",
              nrow(pooled_df), length(pooled)))

  list(
    pooled_data = pooled_df,
    common_gee_vars = common_vars,
    country_names = names(pooled)
  )
}


#' Run area-level LOCO cross-validation
#'
#' Train area-level SL on N-1 countries' Admin-2 data, predict on held-out
#' country. Tests both MSE and NLL loss functions.
#'
#' @param pooled Output from build_area_loco_dataset()
#' @param sl_library List of mlr3 learner specs (or NULL for default)
#' @return Data.frame with error metrics per held-out country and loss function
run_area_loco <- function(pooled, sl_library = NULL) {

  if (!requireNamespace("mlr3superlearner", quietly = TRUE))
    stop("mlr3superlearner required")
  library(mlr3extralearners, quietly = TRUE)

  data <- pooled$pooled_data
  gee_vars <- pooled$common_gee_vars
  countries <- pooled$country_names

  if (is.null(sl_library)) {
    sl_library <- list(
      list("glmnet", alpha = 1, id = "lasso"),
      list("glmnet", alpha = 0, id = "ridge"),
      list("glmnet", alpha = 0.5, id = "elastic_net"),
      list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
      list("xgboost", max_depth = 2, eta = 0.05, nrounds = 200,
           min_child_weight = 3, subsample = 0.8, id = "xgb")
    )
  }

  # Monkey-patch NLL
  env <- asNamespace("mlr3superlearner")
  if (exists("loss_nll", envir = env)) {
    clipped_nll <- function(x, y) {
      x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
      -mean(y * log(x) + (1 - y) * log(1 - x))
    }
    tryCatch(assignInNamespace("loss_nll", clipped_nll, "mlr3superlearner"),
             error = function(e) NULL)
  }

  results <- list()

  for (held_out in countries) {
    train <- data |> dplyr::filter(country != held_out)
    test  <- data |> dplyr::filter(country == held_out)

    if (nrow(train) < 5 || nrow(test) < 3) {
      cat(sprintf("  [area_loco] Skipping %s (train=%d, test=%d)\n",
                  held_out, nrow(train), nrow(test)))
      next
    }

    cat(sprintf("\n  [area_loco] Held out: %s (train=%d, test=%d)\n",
                held_out, nrow(train), nrow(test)))

    # Prepare data
    valid_vars <- gee_vars[sapply(gee_vars, function(v) {
      x <- train[[v]]
      sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
    })]

    X_train <- as.data.frame(train[, valid_vars, drop = FALSE])
    X_test  <- as.data.frame(test[, valid_vars, drop = FALSE])

    # Impute NAs
    col_medians <- sapply(X_train, median, na.rm = TRUE)
    for (col in names(X_train)) {
      na_idx <- is.na(X_train[[col]]) | !is.finite(X_train[[col]])
      if (any(na_idx)) X_train[[col]][na_idx] <- col_medians[col]
      na_idx <- is.na(X_test[[col]]) | !is.finite(X_test[[col]])
      if (any(na_idx)) X_test[[col]][na_idx] <- col_medians[col]
    }

    for (otype in c("continuous", "logit")) {
      Y_train_raw <- train$svy_prev
      Y_test_raw  <- test$svy_prev

      if (otype == "logit") {
        # Logit-transform for beta-inspired loss
        Y_train <- qlogis(pmin(pmax(Y_train_raw, 0.001), 0.999))
        fit_otype <- "continuous"
      } else {
        Y_train <- Y_train_raw
        fit_otype <- "continuous"
      }

      fit_df <- data.frame(Y = Y_train, X_train)
      pred_df <- data.frame(Y = 0, X_test)  # dummy Y for prediction

      n_folds <- min(nrow(fit_df), 10L)

      sl_fit <- tryCatch(
        mlr3superlearner::mlr3superlearner(
          data = fit_df, target = "Y", library = sl_library,
          outcome_type = fit_otype, folds = n_folds
        ),
        error = function(e) {
          cat(sprintf("    %s SL failed: %s\n", otype, e$message))
          NULL
        }
      )

      if (is.null(sl_fit)) next

      preds <- tryCatch(as.numeric(predict(sl_fit, pred_df)),
                          error = function(e) rep(NA_real_, nrow(pred_df)))
      # Back-transform if logit
      if (otype == "logit") preds <- plogis(preds)
      preds <- pmin(pmax(preds, 0), 1)

      valid <- !is.na(preds) & !is.na(Y_test_raw)
      if (sum(valid) < 3) next

      winner <- tryCatch(names(which.max(sl_fit$weights)), error = function(e) "?")

      results[[paste0(held_out, "_", otype)]] <- data.frame(
        held_out = held_out,
        model_type = otype,
        n_train = nrow(train),
        n_test = nrow(test),
        n_valid_vars = length(valid_vars),
        obs_prev = mean(Y_test_raw, na.rm = TRUE),
        pred_prev = mean(preds, na.rm = TRUE),
        mae_pp = round(mean(abs(Y_test_raw[valid] - preds[valid])) * 100, 2),
        rmse_pp = round(sqrt(mean((Y_test_raw[valid] - preds[valid])^2)) * 100, 2),
        pearson_r = round(cor(Y_test_raw[valid], preds[valid]), 3),
        mean_bias_pp = round(mean(preds[valid] - Y_test_raw[valid]) * 100, 2),
        winner = winner,
        stringsAsFactors = FALSE
      )

      cat(sprintf("    %s: MAE=%.1f pp, r=%.3f, pred_nat=%.1f%%, winner=%s\n",
                  otype,
                  results[[paste0(held_out, "_", otype)]]$mae_pp,
                  results[[paste0(held_out, "_", otype)]]$pearson_r,
                  results[[paste0(held_out, "_", otype)]]$pred_prev * 100,
                  winner))
    }
  }

  if (length(results) == 0) return(data.frame())
  dplyr::bind_rows(results)
}


#' Run within-country area-level comparison for one outcome
#'
#' Fits area SL (MSE + NLL), compares with individual SL aggregation.
#' Convenience wrapper for use in targets pipeline.
#'
#' @param svy_admin2 Survey-weighted Admin-2 prevalence
#' @param sl_admin2 Individual-level SL aggregated to Admin-2
#' @param gee_admin2 GEE Admin-2 covariates
#' @param cc Country config
#' @param oc Outcome config
#' @return List with area_mse, area_nll, comparison
run_area_comparison <- function(svy_admin2, sl_admin2, gee_admin2, cc, oc) {

  cat(sprintf("\n[area_comparison] %s — %s\n", cc$country, oc$tag))

  # Fit area-level SL with MSE loss
  area_mse <- tryCatch(
    fit_area_sl(svy_admin2, gee_admin2, outcome_type = "continuous"),
    error = function(e) {
      cat(sprintf("  MSE fit failed: %s\n", e$message))
      list(area_preds = NULL, metrics = NULL, model_type = "continuous")
    }
  )

  # Fit area-level SL with NLL loss (beta-inspired)
  area_nll <- tryCatch(
    fit_area_sl(svy_admin2, gee_admin2, outcome_type = "binomial"),
    error = function(e) {
      cat(sprintf("  NLL fit failed: %s\n", e$message))
      list(area_preds = NULL, metrics = NULL, model_type = "binomial")
    }
  )

  # Compare all three approaches
  comparison <- compare_admin2_approaches(
    sl_admin2, area_mse, area_nll, svy_admin2, cc$country, oc$tag, cc = cc
  )

  list(
    area_mse = area_mse,
    area_nll = area_nll,
    comparison = comparison
  )
}
