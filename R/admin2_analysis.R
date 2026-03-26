# =============================================================================
# R/admin2_analysis.R
#
# Functions for Admin2-level analysis:
#   - Aggregate individual SL predictions to Admin2
#   - Survey-weighted Admin2 prevalence (srvyr)
#   - Area-level model: GEE rasters -> Admin2 prevalence (for unsurveyed areas)
#   - Error analysis
# =============================================================================

#' Aggregate individual-level SL predictions to Admin2 prevalence
#'
#' @param sl_fit Output from fit_sl_models()
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with Admin2, n, sl_prev
aggregate_admin2_sl <- function(sl_fit, outcome_data, cc, oc) {
  if (is.null(sl_fit)) return(NULL)

  fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
  if (is.null(fit)) return(NULL)

  d   <- outcome_data$data
  res <- fit$res

  # Predicted prevalence: binary model gives P(deficient), continuous uses threshold
  if (!is.null(sl_fit$bin_fit)) {
    d$deficient <- res$yhat_full
  } else {
    d$deficient <- apply_threshold(res$yhat_full, oc$cutoff, oc$cutoff_dir)
  }

  admin2_col <- cc$admin2_col

  d %>%
    dplyr::filter(!is.na(.data[[admin2_col]])) %>%
    dplyr::group_by(.data[[admin2_col]]) %>%
    dplyr::summarise(
      n_sl    = dplyr::n(),
      sl_prev = mean(deficient, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(Admin2 = dplyr::all_of(admin2_col))
}


#' Compute survey-weighted Admin2 prevalence
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with Admin2, svy_prev, svy_prev_se, n_svy, etc.
compute_svy_admin2 <- function(outcome_data, cc, oc) {
  d <- outcome_data$data

  svy_cols <- c(cc$psu_col, cc$weight_col, cc$admin2_col,
                cc$admin1_col, oc$binary)
  if (!is.null(cc$strata_col)) svy_cols <- c(svy_cols, cc$strata_col)
  svy_cols <- unique(svy_cols[svy_cols %in% colnames(d)])

  svy_df <- d %>%
    dplyr::select(dplyr::all_of(svy_cols)) %>%
    dplyr::filter(!is.na(.data[[oc$binary]]),
                  !is.na(.data[[cc$admin2_col]]),
                  !is.na(.data[[cc$weight_col]])) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(cc$psu_col), as.factor),
      !!oc$binary := as.numeric(.data[[oc$binary]]),
      !!cc$weight_col := as.numeric(.data[[cc$weight_col]])
    )

  if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(svy_df))
    svy_df <- svy_df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(cc$strata_col), as.factor))

  options(survey.lonely.psu = "adjust")

  if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(svy_df)) {
    svy_des <- srvyr::as_survey_design(
      svy_df,
      ids     = !!rlang::sym(cc$psu_col),
      strata  = !!rlang::sym(cc$strata_col),
      weights = !!rlang::sym(cc$weight_col),
      nest    = TRUE
    )
  } else {
    svy_des <- srvyr::as_survey_design(
      svy_df,
      ids     = !!rlang::sym(cc$psu_col),
      weights = !!rlang::sym(cc$weight_col)
    )
  }

  # Some Admin2 areas have too few clusters for survey_mean to compute SEs.
  # Compute per-group safely: first get simple counts, then try survey_mean.
  svy_admin2 <- tryCatch({
    svy_des %>%
      dplyr::group_by(.data[[cc$admin2_col]]) %>%
      dplyr::summarise(
        svy_prev = srvyr::survey_mean(!!rlang::sym(oc$binary),
                                       vartype = c("se", "ci"), na.rm = TRUE),
        n_svy    = dplyr::n(),
        .groups  = "drop"
      ) %>%
      dplyr::rename(Admin2 = dplyr::all_of(cc$admin2_col)) %>%
      dplyr::mutate(
        svy_cv       = dplyr::if_else(svy_prev > 0, svy_prev_se / svy_prev, NA_real_),
        svy_ci_width = svy_prev_upp - svy_prev_low
      )
  }, error = function(e) {
    # Fallback: compute per Admin2 area individually, skipping failures
    cat(sprintf("  [svy_admin2] Grouped survey_mean failed (%s), computing per-area...\n",
                conditionMessage(e)))
    admin2_vals <- unique(svy_df[[cc$admin2_col]])
    admin2_vals <- admin2_vals[!is.na(admin2_vals)]

    results <- lapply(admin2_vals, function(a2) {
      sub <- svy_df[svy_df[[cc$admin2_col]] == a2, ]
      n_clusters <- length(unique(sub[[cc$psu_col]]))
      n_obs <- nrow(sub)

      if (n_obs < 2) return(NULL)

      res <- tryCatch({
        if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(sub)) {
          des <- srvyr::as_survey_design(sub,
            ids = !!rlang::sym(cc$psu_col),
            strata = !!rlang::sym(cc$strata_col),
            weights = !!rlang::sym(cc$weight_col), nest = TRUE)
        } else {
          des <- srvyr::as_survey_design(sub,
            ids = !!rlang::sym(cc$psu_col),
            weights = !!rlang::sym(cc$weight_col))
        }
        sm <- des %>%
          dplyr::summarise(
            svy_prev = srvyr::survey_mean(!!rlang::sym(oc$binary),
                                           vartype = c("se", "ci"), na.rm = TRUE),
            .groups = "drop")
        data.frame(Admin2 = a2, svy_prev = sm$svy_prev,
                   svy_prev_se = sm$svy_prev_se,
                   svy_prev_low = sm$svy_prev_low,
                   svy_prev_upp = sm$svy_prev_upp,
                   n_svy = n_obs, stringsAsFactors = FALSE)
      }, error = function(e2) {
        # Last resort: unweighted mean, no SE
        prev <- mean(sub[[oc$binary]], na.rm = TRUE)
        data.frame(Admin2 = a2, svy_prev = prev,
                   svy_prev_se = NA_real_,
                   svy_prev_low = NA_real_,
                   svy_prev_upp = NA_real_,
                   n_svy = n_obs, stringsAsFactors = FALSE)
      })
      res
    })

    valid_results <- results[!sapply(results, is.null)]
    if (length(valid_results) == 0) {
      return(data.frame(Admin2 = character(), svy_prev = numeric(),
                        svy_prev_se = numeric(), svy_prev_low = numeric(),
                        svy_prev_upp = numeric(), n_svy = integer(),
                        svy_cv = numeric(), svy_ci_width = numeric(),
                        stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, valid_results)
    out$svy_cv <- ifelse(out$svy_prev > 0 & !is.na(out$svy_prev_se),
                         out$svy_prev_se / out$svy_prev, NA_real_)
    out$svy_ci_width <- out$svy_prev_upp - out$svy_prev_low
    out
  })

  svy_admin2
}


#' Fit area-level model (GEE rasters -> Admin2 prevalence) and predict everywhere
#'
#' This trains an elastic net on surveyed Admin2 units and predicts to ALL
#' polygons (including unsurveyed areas with no biomarker data).
#'
#' @param svy_admin2 Survey-weighted Admin2 prevalence (from compute_svy_admin2)
#' @param cc Country config
#' @param oc Outcome config
#' @param params Pipeline parameters
#' @return list with area_preds, loo_summary, fit_area (glmnet object)
fit_area_level_model <- function(svy_admin2, cc, oc, params) {
  if (is.null(svy_admin2) || !is.data.frame(svy_admin2) ||
      nrow(svy_admin2) < 3 || all(is.na(svy_admin2$svy_prev))) {
    warning("Too few survey Admin2 estimates for area-level model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = NULL))
  }

  # Download GADM boundaries (with retry and error handling)
  gadm_raw <- tryCatch(
    geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
    error = function(e) {
      # Retry once — GADM server can be flaky
      tryCatch(
        geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(gadm_raw)) {
    warning(sprintf("GADM download failed for %s — skipping area model", cc$gadm_code))
    return(list(area_preds = NULL, loo_summary = NULL, polygons = NULL))
  }
  gadm <- sf::st_as_sf(gadm_raw)
  gadm$Admin2 <- gadm$NAME_2
  all_polys <- gadm

  # Identify surveyed/unsurveyed
  surveyed_names <- svy_admin2$Admin2
  n_surveyed     <- length(surveyed_names)
  n_unsurveyed   <- sum(!all_polys$Admin2 %in% surveyed_names)

  cat(sprintf("[area_model] %d surveyed, %d unsurveyed Admin2 units\n",
              n_surveyed, n_unsurveyed))

  # Find raster directory (try both naming conventions)
  raster_dir <- cc$raster_dir
  if (!dir.exists(raster_dir)) {
    # Try alternate name with space
    alt <- gsub("_GEE_rasters$", "_GEE rasters", raster_dir)
    if (dir.exists(alt)) raster_dir <- alt
  }

  if (!dir.exists(raster_dir)) {
    warning("GEE raster directory not found: ", raster_dir)
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }

  # Extract zonal means for all polygons
  tif_files <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  cat(sprintf("  %d .tif files found\n", length(tif_files)))

  gee_admin2 <- data.frame(Admin2 = all_polys$Admin2)
  for (tif in tif_files) {
    base <- tools::file_path_sans_ext(basename(tif))
    varname <- paste0("gee_", tolower(gsub("[^A-Za-z0-9]+", "_", base)))
    varname <- sub("^_|_$", "", varname)

    r <- tryCatch(terra::rast(tif), error = function(e) NULL)
    if (is.null(r)) next
    if (terra::nlyr(r) > 1) r <- r[[1]]

    vals <- tryCatch(
      exactextractr::exact_extract(r, all_polys, fun = "mean"),
      error = function(e) rep(NA_real_, nrow(all_polys))
    )
    gee_admin2[[varname]] <- vals
  }

  gee_vars <- setdiff(colnames(gee_admin2), "Admin2")
  cat(sprintf("  Extracted %d GEE variables\n", length(gee_vars)))

  # Training set: surveyed areas with GEE covariates

  train_df <- gee_admin2 %>%
    dplyr::inner_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    ) %>%
    dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))

  if (nrow(train_df) < 3) {
    warning("Too few surveyed Admin2 areas with valid prevalence for area model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }

  # Drop zero-variance or all-NA predictors
  valid_vars <- gee_vars[sapply(gee_vars, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  if (length(valid_vars) < 2) {
    warning("Too few valid GEE predictors for area-level model")
    return(list(area_preds = NULL, loo_summary = NULL, polygons = all_polys))
  }

  X_train <- as.matrix(train_df[, valid_vars])
  Y_train <- train_df$svy_prev
  col_medians <- apply(X_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_train))) {
    na_idx <- is.na(X_train[, j])
    if (any(na_idx)) X_train[na_idx, j] <- col_medians[j]
  }

  X_all <- as.matrix(gee_admin2[, valid_vars])
  for (j in seq_len(ncol(X_all))) {
    na_idx <- is.na(X_all[, j])
    if (any(na_idx)) X_all[na_idx, j] <- col_medians[j]
  }

  # LOO-CV
  cat("  Running LOO-CV...\n")
  loo_preds <- rep(NA_real_, nrow(X_train))
  n_loo_folds <- max(min(nrow(X_train) - 1, 10), 3)
  for (i in seq_len(nrow(X_train))) {
    n_remain <- nrow(X_train) - 1
    if (n_remain < 3) next
    fit_loo <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train[-i, , drop = FALSE], y = Y_train[-i],
        alpha = 0.5, nfolds = min(n_remain, 10),
        type.measure = "mse"
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_loo))
      loo_preds[i] <- as.numeric(predict(fit_loo, X_train[i, , drop = FALSE],
                                          s = "lambda.min"))
  }
  loo_preds <- pmin(pmax(loo_preds, 0), 1)

  loo_valid <- !is.na(loo_preds)
  loo_summary <- if (sum(loo_valid) >= 3) {
    data.frame(
      outcome      = oc$tag,
      country      = cc$country,
      n_train      = sum(loo_valid),
      n_predictors = length(valid_vars),
      loo_mae_pp   = round(mean(abs(Y_train[loo_valid] - loo_preds[loo_valid])) * 100, 2),
      loo_rmse_pp  = round(sqrt(mean((Y_train[loo_valid] - loo_preds[loo_valid])^2)) * 100, 2),
      loo_r        = round(cor(Y_train[loo_valid], loo_preds[loo_valid]), 3)
    )
  } else NULL

  if (!is.null(loo_summary))
    cat(sprintf("  LOO-CV: MAE = %.1f pp, r = %.3f\n",
                loo_summary$loo_mae_pp, loo_summary$loo_r))

  # Final model
  set.seed(params$seed)
  fit_area <- tryCatch(
    glmnet::cv.glmnet(
      x = X_train, y = Y_train, alpha = 0.5,
      nfolds = max(min(nrow(X_train), 10), 3), type.measure = "mse"
    ),
    error = function(e) {
      warning("Area-level glmnet failed: ", e$message)
      NULL
    }
  )

  if (is.null(fit_area)) {
    return(list(area_preds = NULL, loo_summary = loo_summary, polygons = all_polys))
  }

  yhat_all <- pmin(pmax(as.numeric(predict(fit_area, X_all, s = "lambda.min")), 0), 1)

  area_preds <- data.frame(
    Admin2         = gee_admin2$Admin2,
    area_pred_prev = yhat_all,
    has_survey     = gee_admin2$Admin2 %in% surveyed_names
  ) %>%
    dplyr::left_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    )

  # Bootstrap CIs
  B_AREA <- params$B_area
  cat(sprintf("  Bootstrapping area model (B = %d)...\n", B_AREA))
  set.seed(params$seed)
  boot_all <- matrix(NA_real_, nrow = nrow(X_all), ncol = B_AREA)

  for (b in seq_len(B_AREA)) {
    idx_b <- sample(nrow(X_train), replace = TRUE)
    fit_b <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train[idx_b, , drop = FALSE], y = Y_train[idx_b],
        alpha = 0.5, nfolds = min(length(unique(idx_b)), 10),
        type.measure = "mse"
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b))
      boot_all[, b] <- pmin(pmax(
        as.numeric(predict(fit_b, X_all, s = "lambda.min")), 0), 1)
  }

  valid_boots <- sum(!is.na(boot_all[1, ]))
  if (valid_boots >= 3) {
    area_preds$boot_mean <- apply(boot_all, 1, mean, na.rm = TRUE)
    area_preds$ci_lo     <- apply(boot_all, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
    area_preds$ci_hi     <- apply(boot_all, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
    area_preds$ci_width  <- area_preds$ci_hi - area_preds$ci_lo
  }

  list(
    area_preds  = area_preds,
    loo_summary = loo_summary,
    polygons    = all_polys,
    fit_area    = fit_area,
    gee_vars    = valid_vars
  )
}


#' Compute error metrics between SL predictions and survey estimates
#'
#' @param sl_admin2 Admin2 SL predictions (from aggregate_admin2_sl)
#' @param svy_admin2 Survey-weighted Admin2 prevalence (from compute_svy_admin2)
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with error metrics, or NULL if inputs are missing
compute_admin2_error <- function(sl_admin2, svy_admin2, cc, oc) {
  if (is.null(sl_admin2) || is.null(svy_admin2)) return(NULL)

  merged <- sl_admin2 %>%
    dplyr::inner_join(
      svy_admin2 %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      by = "Admin2"
    ) %>%
    dplyr::filter(!is.na(sl_prev), !is.na(svy_prev))

  if (nrow(merged) < 2) return(NULL)

  data.frame(
    outcome      = oc$tag,
    country      = cc$country,
    n_admin2     = nrow(merged),
    mae_pp       = round(mean(abs(merged$sl_prev - merged$svy_prev)) * 100, 2),
    rmse_pp      = round(sqrt(mean((merged$sl_prev - merged$svy_prev)^2)) * 100, 2),
    pearson_r    = round(cor(merged$sl_prev, merged$svy_prev), 3),
    mean_bias_pp = round(mean(merged$sl_prev - merged$svy_prev) * 100, 2)
  )
}
