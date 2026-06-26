# =============================================================================
# R/oos_prediction.R
#
# Out-of-sample prediction for unsurveyed countries using GEE rasters.
# Takes an area-level model trained on surveyed countries and applies it
# to GADM Admin-2 polygons in a new country, extracting the same GEE
# covariates from locally stored rasters.
# =============================================================================

#' Normalize a raster filename into a canonical GEE variable name
#'
#' Strips country names, standardizes separators, and creates names that
#' match across countries. E.g.:
#'   "CCNL_2013_Gambia.tif"         -> "gee_ccnl_2013"
#'   "CCNL_2013_Cote_dIvoire.tif"   -> "gee_ccnl_2013"
#'   "NDVI_Ghana_2017.tif"          -> "gee_ndvi_2017"
#'   "FLDAS_2017_Annual_Mean_Gambia" -> "gee_fldas_2017_annual_mean"
#'   "FLDAS_2017_Monthly_Cote_dIvoire" -> "gee_fldas_2017_monthly"
#'
#' @param filename Basename of .tif file (without extension)
#' @return Canonical variable name starting with "gee_"
canonicalize_gee_varname <- function(filename) {
  # Remove file extension
  base <- tools::file_path_sans_ext(basename(filename))

  # Country names to strip (order matters — longer first)
  country_patterns <- c(
    "Cote_dIvoire", "Cote_d_Ivoire", "CotedIvoire",
    "Sierra_Leone", "Sierra Leone", "SierraLeone",
    "Gambia", "Ghana", "Malawi", "Nigeria", "Senegal",
    "Burkina_Faso", "BurkinaFaso", "Guinea", "Mali", "Niger",
    "Togo", "Benin", "Liberia"
  )

  cleaned <- base
  for (pat in country_patterns) {
    # Remove country name and any surrounding underscores/spaces
    cleaned <- gsub(paste0("_?", pat, "_?"), "_", cleaned, ignore.case = TRUE)
  }

  # Clean up: lowercase, replace non-alphanumeric with _, collapse multiple _
  cleaned <- tolower(cleaned)
  cleaned <- gsub("[^a-z0-9]+", "_", cleaned)
  cleaned <- gsub("_+", "_", cleaned)
  cleaned <- sub("^_|_$", "", cleaned)

  paste0("gee_", cleaned)
}


#' Extract GEE covariates for a new country's Admin-2 polygons
#'
#' @param gadm_code ISO3 code for geodata::gadm()
#' @param raster_dir Path to folder of .tif rasters
#' @param gadm_path Path for caching GADM downloads
#' @return data.frame with Admin1, Admin2, and gee_* columns
extract_gee_for_country <- function(gadm_code, raster_dir, gadm_path = here::here("data", "gadm")) {

  gadm_raw <- geodata::gadm(gadm_code, level = 2, path = gadm_path)
  gadm <- sf::st_as_sf(gadm_raw)
  gadm$Admin2 <- gadm$NAME_2
  gadm$Admin1 <- gadm$NAME_1

  tif_files <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  cat(sprintf("[oos] %s: %d Admin2 polygons, %d .tif files\n",
              gadm_code, nrow(gadm), length(tif_files)))

  gee_df <- data.frame(Admin1 = gadm$Admin1, Admin2 = gadm$Admin2)

  for (tif in tif_files) {
    varname <- canonicalize_gee_varname(tif)

    r <- tryCatch(terra::rast(tif), error = function(e) NULL)
    if (is.null(r)) next
    if (terra::nlyr(r) > 1) r <- r[[1]]

    vals <- tryCatch(
      exactextractr::exact_extract(r, gadm, fun = "mean"),
      error = function(e) rep(NA_real_, nrow(gadm))
    )

    # If duplicate canonical name (e.g., multiple years), keep first
    if (!varname %in% colnames(gee_df)) {
      gee_df[[varname]] <- vals
    }
  }

  gee_vars <- setdiff(colnames(gee_df), c("Admin1", "Admin2"))
  cat(sprintf("[oos] Extracted %d GEE variables for %s\n", length(gee_vars), gadm_code))

  attr(gee_df, "polygons") <- gadm
  gee_df
}


#' CORAL domain alignment for out-of-sample prediction
#'
#' Aligns the training (source) covariate covariance to the target country's
#' (unsurveyed) covariate covariance, then refits an elastic net on the aligned
#' features. Transductive but label-honest: only the target's COVARIATES are
#' used to estimate its covariance — no target outcome is ever touched. This
#' corrects the covariate shift between surveyed and unsurveyed areas, which is
#' the setting where CORAL is genuinely justified (unlike within-country k-fold).
#'
#' @param Xs n_train x p raw (imputed) training covariate matrix
#' @param Ys length-n_train training prevalence
#' @param Xt n_target x p raw (imputed) target covariate matrix (same columns)
#' @param ridge covariance ridge for the matrix square roots
#' @return list(fit, X_source, X_target, loo) — fit + both matrices in the
#'   aligned feature space, and LOO residuals for conformal calibration
.coral_oos_refit <- function(Xs, Ys, Xt, ridge = 1e-3) {
  mu  <- colMeans(Xs); sdv <- apply(Xs, 2, stats::sd); sdv[sdv < 1e-8] <- 1
  Zs  <- sweep(sweep(Xs, 2, mu, "-"), 2, sdv, "/")
  Zt  <- sweep(sweep(Xt, 2, mu, "-"), 2, sdv, "/")
  p <- ncol(Zs); Ip <- diag(p)
  msqrt <- function(M, pw) {
    e <- eigen(M, symmetric = TRUE); v <- pmax(e$values, 1e-8)
    e$vectors %*% diag(v^pw, p) %*% t(e$vectors)
  }
  A <- msqrt(stats::cov(Zs) + ridge * Ip, -0.5) %*%
       msqrt(stats::cov(Zt) + ridge * Ip,  0.5)
  Xs_a <- Zs %*% A          # align source onto target covariance
  Xt_a <- Zt                # target stays on its own standardized scale
  colnames(Xs_a) <- colnames(Xt_a) <- colnames(Xs)
  # Reuse the deterministic-fold cv.glmnet wrapper when available.
  cvg <- if (exists(".cv_glmnet")) .cv_glmnet else
           function(x, y, ...) glmnet::cv.glmnet(x, y, ...)
  fit <- cvg(Xs_a, Ys, alpha = 0.5, nfolds = max(min(nrow(Xs_a), 10), 3))
  loo <- rep(NA_real_, nrow(Xs_a))
  for (i in seq_len(nrow(Xs_a))) {
    lf <- tryCatch(cvg(Xs_a[-i, , drop = FALSE], Ys[-i], alpha = 0.5,
                       nfolds = max(min(nrow(Xs_a) - 1, 10), 3)),
                   error = function(e) NULL)
    if (!is.null(lf))
      loo[i] <- abs(Ys[i] -
        as.numeric(predict(lf, Xs_a[i, , drop = FALSE], s = "lambda.min")))
  }
  list(fit = fit, X_source = Xs_a, X_target = Xt_a, loo = loo[!is.na(loo)])
}


#' Predict micronutrient deficiency in an unsurveyed country
#'
#' Uses an area-level glmnet model (from fit_area_level_model) trained on
#' surveyed countries' Admin-2 prevalences, and applies it to a new country's
#' Admin-2 GEE covariates. Includes weighted conformal prediction intervals.
#'
#' @param area_model Output from predict_oos_pooled training step
#'   (must contain fit_area, gee_vars, X_train, loo_residuals)
#' @param oos_gee Output from extract_gee_for_country() for the target country
#' @param outcome_tag Character label for the outcome
#' @param alpha Nominal miscoverage rate (default 0.1 for 90% CIs)
#' @return list with predictions data.frame and polygons
predict_oos_country <- function(area_model, oos_gee, outcome_tag, alpha = 0.1) {

  if (is.null(area_model$fit_area) || is.null(area_model$gee_vars)) {
    warning("Area model has no fit or gee_vars — cannot predict OOS")
    return(NULL)
  }

  model_vars <- area_model$gee_vars
  fit        <- area_model$fit_area
  polygons   <- attr(oos_gee, "polygons")

  # Build predictor matrix, padding missing columns with 0
  available <- intersect(model_vars, colnames(oos_gee))
  missing   <- setdiff(model_vars, colnames(oos_gee))
  cat(sprintf("[oos_predict] %d/%d model variables found in target country (%d missing)\n",
              length(available), length(model_vars), length(missing)))

  if (length(missing) > 0) {
    cat(sprintf("  Missing vars: %s\n", paste(head(missing, 10), collapse = ", ")))
  }

  X_oos <- as.matrix(oos_gee[, available, drop = FALSE])
  # Pad missing columns with 0
  if (length(missing) > 0) {
    pad <- matrix(0, nrow = nrow(X_oos), ncol = length(missing))
    colnames(pad) <- missing
    X_oos <- cbind(X_oos, pad)
  }
  # Reorder to match model
  X_oos <- X_oos[, model_vars, drop = FALSE]

  # Impute any NAs with column medians from training
  col_medians <- if (!is.null(area_model$X_train)) {
    apply(area_model$X_train, 2, median, na.rm = TRUE)
  } else {
    apply(X_oos, 2, median, na.rm = TRUE)
  }
  for (j in seq_len(ncol(X_oos))) {
    na_idx <- is.na(X_oos[, j])
    if (any(na_idx)) X_oos[na_idx, j] <- col_medians[j]
  }

  # ── Optional CORAL domain alignment (covariate-shift correction) ──────────
  # When area_model$align == "coral", realign the source covariance onto this
  # target country's covariance and refit before predicting. Uses target
  # covariates only (no labels), so it stays out-of-sample-honest.
  Xtr_used <- area_model$X_train
  if (identical(area_model$align, "coral") &&
      !is.null(area_model$Y_train) && !is.null(area_model$X_train)) {
    al <- tryCatch(.coral_oos_refit(area_model$X_train, area_model$Y_train, X_oos),
                   error = function(e) {
                     cat(sprintf("[oos_predict] CORAL align failed (%s) — using raw model\n",
                                 conditionMessage(e))); NULL })
    if (!is.null(al)) {
      fit      <- al$fit
      X_oos    <- al$X_target
      Xtr_used <- al$X_source
      area_model$loo_residuals <- al$loo
      cat(sprintf("[oos_predict] CORAL alignment applied (%d source, %d target areas)\n",
                  nrow(Xtr_used), nrow(X_oos)))
    }
  }

  # Predict
  yhat <- pmin(pmax(
    as.numeric(predict(fit, X_oos, s = "lambda.min")),
    0), 1)

  # ── Conformal prediction intervals ───────────────────────────────────────
  # Use weighted conformal with covariate-shift density ratio.
  # The density ratio p_target(x)/p_train(x) is estimated via logistic
  # regression: P(target=1 | x). The importance weight for each training
  # calibration point is w_i = p_hat(x_i) / (1 - p_hat(x_i)).
  ci_lo <- rep(NA_real_, nrow(X_oos))
  ci_hi <- rep(NA_real_, nrow(X_oos))
  ci_method <- "none"

  loo_res <- area_model$loo_residuals
  X_train <- Xtr_used   # aligned source matrix when CORAL is active, else raw

  if (!is.null(loo_res) && length(loo_res) >= 5 && !is.null(X_train)) {
    n_cal <- length(loo_res)

    # ── Estimate density ratio via logistic regression ──
    # Combine training and target covariate matrices
    X_combined <- rbind(X_train[seq_len(n_cal), , drop = FALSE], X_oos)
    label <- c(rep(0L, n_cal), rep(1L, nrow(X_oos)))  # 0=train, 1=target

    # Fit logistic regression (regularized to handle collinearity)
    density_fit <- tryCatch({
      glmnet::cv.glmnet(x = X_combined, y = label, family = "binomial",
                          alpha = 0.5,
                          nfolds = min(10, max(3, floor(nrow(X_combined) / 5))),
                          type.measure = "deviance")
    }, error = function(e) NULL)

    if (!is.null(density_fit)) {
      # P(target | x) for training points
      p_train <- as.numeric(predict(density_fit, X_train[seq_len(n_cal), , drop = FALSE],
                                     s = "lambda.min", type = "response"))
      # Clip to avoid division by zero
      p_train <- pmin(pmax(p_train, 0.01), 0.99)
      # Importance weights: p/(1-p) proportional to density ratio
      weights_cal <- p_train / (1 - p_train)

      # Per-target-point weighted conformal quantile
      # For each target point j, the weighted conformal quantile is:
      # q_alpha such that sum(w_i * I(r_i <= q)) / (sum(w_i) + w_j) >= 1 - alpha
      # where w_j = p_j / (1 - p_j)
      p_target <- as.numeric(predict(density_fit, X_oos,
                                      s = "lambda.min", type = "response"))
      p_target <- pmin(pmax(p_target, 0.01), 0.99)
      w_target <- p_target / (1 - p_target)

      for (j in seq_len(nrow(X_oos))) {
        # Weighted quantile
        total_w <- sum(weights_cal) + w_target[j]
        # Sort residuals by size, compute cumulative weighted proportion
        ord <- order(loo_res)
        sorted_res <- loo_res[ord]
        sorted_w   <- weights_cal[ord]
        cum_w <- cumsum(sorted_w) / total_w
        # Find smallest residual threshold where cumulative weight >= 1-alpha
        q_idx <- which(cum_w >= (1 - alpha))[1]
        if (is.na(q_idx)) q_idx <- n_cal
        q_val <- sorted_res[q_idx]

        ci_lo[j] <- pmax(yhat[j] - q_val, 0)
        ci_hi[j] <- pmin(yhat[j] + q_val, 1)
      }
      ci_method <- "weighted_conformal"
      cat(sprintf("[oos_predict] Weighted conformal %d%% CIs computed (density ratio via logistic glmnet)\n",
                  round((1 - alpha) * 100)))
    } else {
      # Fallback: simple (unweighted) conformal
      q_val <- quantile(loo_res, probs = 1 - alpha, na.rm = TRUE)
      ci_lo <- pmax(yhat - q_val, 0)
      ci_hi <- pmin(yhat + q_val, 1)
      ci_method <- "simple_conformal"
      cat(sprintf("[oos_predict] Simple conformal %d%% CIs (density ratio fit failed)\n",
                  round((1 - alpha) * 100)))
    }
  } else if (!is.null(loo_res) && length(loo_res) >= 5) {
    # No X_train available — use simple conformal
    q_val <- quantile(loo_res, probs = 1 - alpha, na.rm = TRUE)
    ci_lo <- pmax(yhat - q_val, 0)
    ci_hi <- pmin(yhat + q_val, 1)
    ci_method <- "simple_conformal"
    cat(sprintf("[oos_predict] Simple conformal %d%% CIs (no covariate shift adjustment)\n",
                round((1 - alpha) * 100)))
  } else {
    cat("[oos_predict] No LOO residuals available — CIs not computed\n")
  }

  preds <- data.frame(
    Admin1    = oos_gee$Admin1,
    Admin2    = oos_gee$Admin2,
    pred_prev = yhat,
    ci_lo     = ci_lo,
    ci_hi     = ci_hi,
    ci_width  = ci_hi - ci_lo,
    stringsAsFactors = FALSE
  )

  cat(sprintf("[oos_predict] %s: %d Admin2 predictions, range [%.1f%%, %.1f%%], median CI width = %.1f pp\n",
              outcome_tag, nrow(preds),
              100 * min(yhat, na.rm = TRUE), 100 * max(yhat, na.rm = TRUE),
              100 * median(preds$ci_width, na.rm = TRUE)))

  # Admin1 summaries
  admin1_preds <- preds %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n_admin2  = dplyr::n(),
      pred_prev = mean(pred_prev, na.rm = TRUE),
      ci_lo     = mean(ci_lo, na.rm = TRUE),
      ci_hi     = mean(ci_hi, na.rm = TRUE),
      ci_width  = mean(ci_width, na.rm = TRUE),
      .groups   = "drop"
    )
  cat(sprintf("[oos_predict] Admin1 summary: %d regions\n", nrow(admin1_preds)))

  list(
    admin2_preds = preds,
    admin1_preds = admin1_preds,
    polygons     = polygons,
    outcome      = outcome_tag,
    n_model_vars = length(model_vars),
    n_matched    = length(available),
    n_missing    = length(missing),
    ci_method    = ci_method,
    ci_alpha     = alpha,
    n_calibration = length(loo_res)
  )
}


#' Run out-of-sample predictions for all outcomes using pooled area models
#'
#' Trains area-level models on ALL surveyed countries' Admin-2 data combined,
#' then predicts for the target unsurveyed country.
#'
#' @param svy_admin2_list Named list of svy_admin2 data.frames (one per country)
#' @param country_configs Country configs for surveyed countries
#' @param oos_gadm_code GADM ISO3 code for target country
#' @param oos_raster_dir Path to GEE rasters for target country
#' @param oc Outcome config
#' @param params Pipeline parameters
#' @param ext_cache_dir Path to external predictor cache directory (default: data/external_cache)
#' @param oos_country_name Country name for loading cached external predictors (e.g., "Cote d'Ivoire")
#' @return Output from predict_oos_country()
predict_oos_pooled <- function(svy_admin2_list, country_configs,
                                oos_gadm_code, oos_raster_dir,
                                oc, params,
                                ext_cache_dir = here::here("data", "external_cache"),
                                oos_country_name = NULL,
                                align = c("none", "coral")) {
  align <- match.arg(align)

  # 1. Pool all surveyed countries' Admin-2 data
  pooled_svy <- dplyr::bind_rows(lapply(names(svy_admin2_list), function(cn) {
    svy <- svy_admin2_list[[cn]]
    if (is.null(svy) || !is.data.frame(svy) || nrow(svy) == 0) return(NULL)
    svy$source_country <- cn
    svy
  }))

  if (nrow(pooled_svy) < 5) {
    warning("Too few pooled Admin2 observations for OOS prediction")
    return(NULL)
  }

  cat(sprintf("[oos_pooled] Pooled %d Admin2 areas from %d countries\n",
              nrow(pooled_svy), length(unique(pooled_svy$source_country))))

  # 2. Extract GEE for all surveyed countries' polygons
  all_gee <- list()
  for (cn in names(country_configs)) {
    cc <- country_configs[[cn]]
    rdir <- cc$raster_dir
    if (!dir.exists(rdir)) {
      alt <- gsub("_GEE_rasters$", "_GEE rasters", rdir)
      if (dir.exists(alt)) rdir <- alt
    }
    if (!dir.exists(rdir)) next

    gadm_raw <- tryCatch(
      geodata::gadm(cc$gadm_code, level = 2, path = here::here("data", "gadm")),
      error = function(e) NULL
    )
    if (is.null(gadm_raw)) next
    gadm <- sf::st_as_sf(gadm_raw)
    gadm$Admin2 <- gadm$NAME_2

    tif_files <- sort(list.files(rdir, pattern = "\\.tif$", full.names = TRUE))
    gee_df <- data.frame(Admin2 = gadm$Admin2)
    for (tif in tif_files) {
      varname <- canonicalize_gee_varname(tif)
      r <- tryCatch(terra::rast(tif), error = function(e) NULL)
      if (is.null(r)) next
      if (terra::nlyr(r) > 1) r <- r[[1]]
      vals <- tryCatch(exactextractr::exact_extract(r, gadm, fun = "mean"),
                        error = function(e) rep(NA_real_, nrow(gadm)))
      if (!varname %in% colnames(gee_df)) {
        gee_df[[varname]] <- vals
      }
    }
    all_gee[[cn]] <- gee_df
  }

  # Find common GEE variables across all countries
  gee_var_lists <- lapply(all_gee, function(df) setdiff(colnames(df), "Admin2"))
  common_gee <- Reduce(intersect, gee_var_lists)
  cat(sprintf("[oos_pooled] %d common GEE variables across %d countries\n",
              length(common_gee), length(all_gee)))

  # 3. Build training data: merge survey prevalences with GEE + external covariates
  train_rows <- list()
  for (cn in names(all_gee)) {
    gee <- all_gee[[cn]]
    svy <- svy_admin2_list[[cn]]
    if (is.null(svy) || nrow(svy) == 0) next

    # Also load cached external predictors for this training country
    cc_train <- country_configs[[cn]]
    if (!is.null(cc_train) && !is.null(ext_cache_dir)) {
      ext_file <- file.path(ext_cache_dir,
                            paste0(tolower(gsub(" ", "_", cn)), "_external_predictors.rds"))
      if (file.exists(ext_file)) {
        ext_df <- readRDS(ext_file)
        # Deduplicate Admin2 in external data (same fix as merge_external.R)
        ext_dup <- ext_df$Admin2[duplicated(ext_df$Admin2)]
        if (length(ext_dup) > 0) {
          ext_num <- names(ext_df)[sapply(ext_df, is.numeric)]
          ext_df <- ext_df |>
            dplyr::group_by(Admin2) |>
            dplyr::summarise(dplyr::across(dplyr::all_of(ext_num),
                                            ~ mean(.x, na.rm = TRUE)),
                             .groups = "drop") |>
            as.data.frame()
        }
        ext_pred_cols <- setdiff(colnames(ext_df), c("Admin1", "Admin2"))
        gee <- dplyr::left_join(gee, ext_df[, c("Admin2", ext_pred_cols)], by = "Admin2")
        cat(sprintf("[oos_pooled] Merged %d external predictors for %s\n",
                    length(ext_pred_cols), cn))
      }
    }

    merged <- dplyr::inner_join(
      svy %>% dplyr::select(Admin2, svy_prev, svy_prev_se, n_svy),
      gee,
      by = "Admin2"
    ) %>% dplyr::filter(!is.na(svy_prev), is.finite(svy_prev))
    if (nrow(merged) > 0) train_rows[[cn]] <- merged
  }
  train_df <- dplyr::bind_rows(train_rows)
  cat(sprintf("[oos_pooled] Training on %d Admin2 areas with survey data\n", nrow(train_df)))

  if (nrow(train_df) < 5) {
    warning("Too few training areas for pooled OOS model")
    return(NULL)
  }

  # 4. Fit pooled glmnet
  valid_vars <- common_gee[sapply(common_gee, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  X_train <- as.matrix(train_df[, valid_vars])
  Y_train <- train_df$svy_prev
  col_medians <- apply(X_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_train))) {
    na_idx <- is.na(X_train[, j])
    if (any(na_idx)) X_train[na_idx, j] <- col_medians[j]
  }

  set.seed(params$seed)
  fit <- tryCatch(
    glmnet::cv.glmnet(x = X_train, y = Y_train, alpha = 0.5,
                       nfolds = max(min(nrow(X_train), 10), 3),
                       type.measure = "mse"),
    error = function(e) { warning("Pooled glmnet failed: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  # ── 4b. Compute LOO residuals for conformal calibration ──────────────────
  # Leave-one-out: for each training area, predict using model fit on the
  # remaining N-1 areas. The absolute residuals form the conformal score set.
  cat("[oos_pooled] Computing LOO residuals for conformal calibration...\n")
  loo_residuals <- rep(NA_real_, nrow(X_train))
  for (i in seq_len(nrow(X_train))) {
    loo_fit <- tryCatch(
      glmnet::cv.glmnet(x = X_train[-i, , drop = FALSE], y = Y_train[-i],
                          alpha = 0.5,
                          nfolds = max(min(nrow(X_train) - 1, 10), 3),
                          type.measure = "mse"),
      error = function(e) NULL
    )
    if (!is.null(loo_fit)) {
      yhat_i <- as.numeric(predict(loo_fit, X_train[i, , drop = FALSE], s = "lambda.min"))
      loo_residuals[i] <- abs(Y_train[i] - yhat_i)
    }
  }
  loo_residuals <- loo_residuals[!is.na(loo_residuals)]
  cat(sprintf("[oos_pooled] %d LOO residuals computed (median = %.4f)\n",
              length(loo_residuals), median(loo_residuals)))

  area_model <- list(fit_area = fit, gee_vars = valid_vars,
                     X_train = X_train, Y_train = Y_train,
                     loo_residuals = loo_residuals,
                     align = align)   # "coral" triggers covariate-shift alignment in predict_oos_country()

  # 5. Extract GEE for target country and merge external predictors
  oos_gee <- extract_gee_for_country(oos_gadm_code, oos_raster_dir)

  # Also load cached external predictors for the target country
  if (!is.null(oos_country_name) && !is.null(ext_cache_dir)) {
    oos_ext_file <- file.path(ext_cache_dir,
                               paste0(tolower(gsub("[' ]", "_", gsub("'", "", oos_country_name))),
                                      "_external_predictors.rds"))
    if (file.exists(oos_ext_file)) {
      oos_ext <- readRDS(oos_ext_file)
      # Deduplicate
      oos_ext_dup <- oos_ext$Admin2[duplicated(oos_ext$Admin2)]
      if (length(oos_ext_dup) > 0) {
        ext_num <- names(oos_ext)[sapply(oos_ext, is.numeric)]
        oos_ext <- oos_ext |>
          dplyr::group_by(Admin2) |>
          dplyr::summarise(dplyr::across(dplyr::all_of(ext_num),
                                          ~ mean(.x, na.rm = TRUE)),
                           .groups = "drop") |>
          as.data.frame()
      }
      ext_pred_cols <- setdiff(colnames(oos_ext), c("Admin1", "Admin2"))
      # Merge external predictors into OOS GEE data
      if (is.data.frame(oos_gee)) {
        oos_gee <- dplyr::left_join(oos_gee, oos_ext[, c("Admin2", ext_pred_cols)], by = "Admin2")
      } else if (inherits(oos_gee, "sf")) {
        oos_gee <- dplyr::left_join(oos_gee, oos_ext[, c("Admin2", ext_pred_cols)], by = "Admin2")
      }
      cat(sprintf("[oos_pooled] Merged %d external predictors for target country %s\n",
                  length(ext_pred_cols), oos_country_name))
    } else {
      cat(sprintf("[oos_pooled] No cached external predictors found for %s at %s\n",
                  oos_country_name, oos_ext_file))
    }
  }

  predict_oos_country(area_model, oos_gee, oc$tag)
}
