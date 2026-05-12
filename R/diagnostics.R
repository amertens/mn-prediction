# =============================================================================
# R/diagnostics.R
#
# Post-hoc model diagnostics from existing cross-validated predictions.
# Does NOT refit any models — reads only from sl_fit$res (Y, yhat_full).
# =============================================================================

#' Extract individual-level predictions from an sl_fit target
#'
#' @param sl_fit Output from fit_sl_models()
#' @param cc Country config
#' @param oc Outcome config
#' @return list with $continuous and $binary data.frames (Y, yhat, country, outcome)
extract_predictions <- function(sl_fit, cc, oc) {
  out <- list()

  if (!is.null(sl_fit$cont_fit)) {
    res <- sl_fit$cont_fit$res
    out$continuous <- data.frame(
      Y       = res$Y,
      yhat    = res$yhat_full,
      country = cc$country,
      outcome = oc$tag,
      type    = "continuous",
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(sl_fit$bin_fit)) {
    res <- sl_fit$bin_fit$res
    out$binary <- data.frame(
      Y       = res$Y,
      yhat    = res$yhat_full,
      country = cc$country,
      outcome = oc$tag,
      type    = "binary",
      stringsAsFactors = FALSE
    )
  }

  out
}


#' Compute extended diagnostics from cross-validated predictions
#'
#' Produces PR-AUC, calibration slope/intercept, Brier skill score,
#' prevalence, event counts — all from existing predictions.
#'
#' @param preds data.frame with Y (0/1) and yhat (predicted probability)
#' @param n_bins Number of bins for calibration table (default 10)
#' @return list with $metrics (1-row df), $calibration_table, $pr_curve
compute_binary_diagnostics <- function(preds, n_bins = 10) {
  Y    <- preds$Y
  yhat <- preds$yhat

  # Clamp predictions to (0,1)
  yhat <- pmin(pmax(yhat, 1e-6), 1 - 1e-6)

  n        <- length(Y)
  n_events <- sum(Y == 1, na.rm = TRUE)
  prev     <- mean(Y, na.rm = TRUE)

  # ── ROC-AUC ──
  roc_auc <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(Y, yhat, quiet = TRUE))),
    error = function(e) NA_real_
  )

  # ── PR-AUC (using ROCR) ──
  pr_auc <- tryCatch({
    pred_obj <- ROCR::prediction(yhat, Y)
    perf_pr  <- ROCR::performance(pred_obj, "prec", "rec")
    # Trapezoidal integration (remove NaN from precision at recall=0)
    rec  <- perf_pr@x.values[[1]]
    prec <- perf_pr@y.values[[1]]
    valid <- !is.na(rec) & !is.na(prec)
    rec  <- rec[valid]
    prec <- prec[valid]
    if (length(rec) > 1) {
      sum(diff(rec) * (prec[-1] + prec[-length(prec)]) / 2)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # ── Brier score and Brier skill score ──
  brier      <- mean((Y - yhat)^2, na.rm = TRUE)
  brier_null <- mean((Y - prev)^2, na.rm = TRUE)  # always-predict-prevalence

  brier_skill <- 1 - brier / brier_null

  # ── Calibration slope and intercept (logistic recalibration) ──
  calib_int   <- NA_real_
  calib_slope <- NA_real_
  calib_slope_lo <- NA_real_
  calib_slope_hi <- NA_real_
  calib_int_lo   <- NA_real_
  calib_int_hi   <- NA_real_
  tryCatch({
    logit_yhat <- log(yhat / (1 - yhat))
    fit <- glm(Y ~ logit_yhat, family = binomial)
    calib_int   <- coef(fit)[1]
    calib_slope <- coef(fit)[2]
    ci <- suppressMessages(stats::confint.default(fit))
    if (is.matrix(ci) && all(c("(Intercept)", "logit_yhat") %in% rownames(ci))) {
      calib_int_lo   <- ci["(Intercept)", 1]
      calib_int_hi   <- ci["(Intercept)", 2]
      calib_slope_lo <- ci["logit_yhat",  1]
      calib_slope_hi <- ci["logit_yhat",  2]
    }
  }, error = function(e) NULL)

  # ── Expected Calibration Error (ECE) ──
  # |mean(yhat) - mean(Y)| within each equal-width bin, weighted by bin size.
  # Uses the same n_bins as the calibration table for consistency.
  ece <- tryCatch({
    bins <- cut(yhat, breaks = n_bins, include.lowest = TRUE)
    n_per <- as.integer(table(bins))
    mp    <- tapply(yhat, bins, mean, na.rm = TRUE)
    mo    <- tapply(Y,    bins, mean, na.rm = TRUE)
    ok    <- !is.na(mp) & !is.na(mo)
    if (sum(n_per[ok]) == 0) NA_real_
    else sum(n_per[ok] * abs(mp[ok] - mo[ok])) / sum(n_per[ok])
  }, error = function(e) NA_real_)

  # ── Maximum Calibration Error (MCE) — worst-case bin gap ──
  mce <- tryCatch({
    bins <- cut(yhat, breaks = n_bins, include.lowest = TRUE)
    mp <- tapply(yhat, bins, mean, na.rm = TRUE)
    mo <- tapply(Y,    bins, mean, na.rm = TRUE)
    ok <- !is.na(mp) & !is.na(mo)
    if (!any(ok)) NA_real_ else max(abs(mp[ok] - mo[ok]))
  }, error = function(e) NA_real_)

  # ── Binned calibration table ──
  cal_table <- tryCatch({
    bins <- cut(yhat, breaks = n_bins, include.lowest = TRUE)
    cal <- data.frame(
      bin          = levels(bins),
      n            = as.integer(table(bins)),
      mean_pred    = tapply(yhat, bins, mean, na.rm = TRUE),
      mean_obs     = tapply(Y, bins, mean, na.rm = TRUE),
      n_events     = tapply(Y, bins, sum, na.rm = TRUE)
    )
    cal$bin_idx <- seq_len(nrow(cal))
    cal
  }, error = function(e) NULL)

  # ── PR curve data ──
  pr_curve <- tryCatch({
    pred_obj <- ROCR::prediction(yhat, Y)
    perf_pr  <- ROCR::performance(pred_obj, "prec", "rec")
    data.frame(
      recall    = perf_pr@x.values[[1]],
      precision = perf_pr@y.values[[1]]
    )
  }, error = function(e) NULL)

  # ── ROC curve data ──
  roc_curve <- tryCatch({
    roc_obj <- pROC::roc(Y, yhat, quiet = TRUE)
    data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    )
  }, error = function(e) NULL)

  metrics <- data.frame(
    n              = n,
    n_events       = n_events,
    prevalence     = round(prev, 4),
    roc_auc        = round(roc_auc, 4),
    pr_auc         = round(pr_auc, 4),
    brier          = round(brier, 4),
    brier_skill    = round(brier_skill, 4),
    calib_int      = round(calib_int, 4),
    calib_int_lo   = round(calib_int_lo, 4),
    calib_int_hi   = round(calib_int_hi, 4),
    calib_slope    = round(calib_slope, 4),
    calib_slope_lo = round(calib_slope_lo, 4),
    calib_slope_hi = round(calib_slope_hi, 4),
    ece            = round(ece, 4),
    mce            = round(mce, 4),
    stringsAsFactors = FALSE
  )

  list(
    metrics           = metrics,
    calibration_table = cal_table,
    pr_curve          = pr_curve,
    roc_curve         = roc_curve
  )
}


#' Compute continuous model diagnostics
#'
#' @param preds data.frame with Y and yhat (continuous scale)
#' @param oc Outcome config (for cutoff/cutoff_dir to derive binary prevalence)
#' @return list with $metrics (1-row df)
compute_continuous_diagnostics <- function(preds, oc = NULL) {
  Y    <- preds$Y
  yhat <- preds$yhat
  n    <- length(Y)

  rmse <- sqrt(mean((Y - yhat)^2, na.rm = TRUE))
  mae  <- mean(abs(Y - yhat), na.rm = TRUE)
  r2   <- 1 - sum((Y - yhat)^2, na.rm = TRUE) / sum((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE)
  r    <- cor(Y, yhat, use = "complete.obs")

  metrics <- data.frame(
    n    = n,
    rmse = round(rmse, 4),
    mae  = round(mae, 4),
    r2   = round(r2, 4),
    r    = round(r, 4),
    stringsAsFactors = FALSE
  )

  # If we have cutoff info, compute prevalence recovery
  if (!is.null(oc) && !is.null(oc$cutoff)) {
    dir <- oc$cutoff_dir %||% "<"
    obs_bin  <- if (dir == "<") as.integer(Y < oc$cutoff) else as.integer(Y >= oc$cutoff)
    pred_bin <- if (dir == "<") as.integer(yhat < oc$cutoff) else as.integer(yhat >= oc$cutoff)
    metrics$prev_obs  <- round(mean(obs_bin, na.rm = TRUE), 4)
    metrics$prev_pred <- round(mean(pred_bin, na.rm = TRUE), 4)
    metrics$prev_bias <- round(metrics$prev_pred - metrics$prev_obs, 4)
  }

  list(metrics = metrics)
}


#' Run all diagnostics for one country x outcome
#'
#' @param sl_fit Output from fit_sl_models()
#' @param cc Country config
#' @param oc Outcome config
#' @return list with combined metrics df, calibration table, PR curve data
run_diagnostics <- function(sl_fit, cc, oc) {
  preds <- extract_predictions(sl_fit, cc, oc)

  results <- list()

  # Continuous diagnostics
  if (!is.null(preds$continuous)) {
    cont_diag <- compute_continuous_diagnostics(preds$continuous, oc)
    cont_row <- cont_diag$metrics
    cont_row$model_type <- "continuous"
    cont_row$country    <- cc$country
    cont_row$outcome    <- oc$tag
    results$continuous_metrics <- cont_row
  }

  # Binary diagnostics
  if (!is.null(preds$binary)) {
    bin_diag <- compute_binary_diagnostics(preds$binary)
    bin_row <- bin_diag$metrics
    bin_row$model_type <- "binary"
    bin_row$country    <- cc$country
    bin_row$outcome    <- oc$tag
    results$binary_metrics       <- bin_row
    results$calibration_table    <- bin_diag$calibration_table
    results$pr_curve             <- bin_diag$pr_curve
    results$roc_curve            <- bin_diag$roc_curve
    if (!is.null(bin_diag$calibration_table)) {
      results$calibration_table$country <- cc$country
      results$calibration_table$outcome <- oc$tag
    }
    if (!is.null(bin_diag$pr_curve)) {
      results$pr_curve$country <- cc$country
      results$pr_curve$outcome <- oc$tag
    }
    if (!is.null(bin_diag$roc_curve)) {
      results$roc_curve$country <- cc$country
      results$roc_curve$outcome <- oc$tag
    }
  }

  results
}


# ── Plot functions ──────────────────────────────────────────────────────────

#' Reliability (calibration) plot
#'
#' @param cal_table Output from compute_binary_diagnostics()$calibration_table
#' @param title Plot title
#' @return ggplot
plot_reliability <- function(cal_table, title = "Calibration plot") {
  if (is.null(cal_table) || nrow(cal_table) == 0) return(NULL)

  ggplot2::ggplot(cal_table, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.7) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(
      title = title,
      x = "Mean predicted probability",
      y = "Observed proportion",
      size = "n"
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_minimal()
}


#' PR curve plot
#'
#' @param pr_data Output from compute_binary_diagnostics()$pr_curve
#' @param prevalence Baseline prevalence (for no-skill line)
#' @param title Plot title
#' @return ggplot
plot_pr_curve <- function(pr_data, prevalence = NULL, title = "Precision-Recall curve") {
  if (is.null(pr_data) || nrow(pr_data) == 0) return(NULL)

  p <- ggplot2::ggplot(pr_data, ggplot2::aes(x = recall, y = precision)) +
    ggplot2::geom_line(color = "#1b9e77", linewidth = 1) +
    ggplot2::labs(title = title, x = "Recall", y = "Precision") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_minimal()

  if (!is.null(prevalence)) {
    p <- p + ggplot2::geom_hline(
      yintercept = prevalence, linetype = "dashed", color = "grey50"
    ) +
      ggplot2::annotate("text", x = 0.5, y = prevalence + 0.03,
                        label = sprintf("No skill (prev = %.1f%%)", prevalence * 100),
                        size = 3, color = "grey40")
  }

  p
}


#' Moran's I on Admin2 SL residuals
#'
#' Tests whether the global SuperLearner leaves spatial structure unexplained
#' at Admin2. Builds a row-standardised k-nearest-neighbour weights matrix on
#' Admin2 polygon centroids (GADM level 2) and runs spdep::moran.test on
#' (svy_prev - sl_prev). Returns NA fields if spatial packages are missing or
#' if too few Admin2s have both survey and SL predictions.
#'
#' @param admin2_df data.frame with Admin1, Admin2, svy_prev, sl_prev
#' @param gadm_code 3-letter ISO code for GADM (e.g. "GMB")
#' @param k_neighbors integer (default 5)
#' @return 1-row data.frame: n, morans_i, expectation, variance, z, p_value
spatial_residual_diagnostics <- function(admin2_df, gadm_code, k_neighbors = 5) {
  if (!requireNamespace("spdep", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("geodata", quietly = TRUE)) {
    return(data.frame(n = NA_integer_, morans_i = NA_real_,
                      expectation = NA_real_, variance = NA_real_,
                      z = NA_real_, p_value = NA_real_,
                      note = "spdep/sf/geodata not installed"))
  }
  if (!exists("morans_i_residuals", mode = "function")) {
    src <- here::here("R", "gwr_analysis.R")
    if (file.exists(src)) source(src)
  }
  polys <- tryCatch(load_admin2_centroids(gadm_code),
                    error = function(e) NULL)
  if (is.null(polys)) {
    return(data.frame(n = NA_integer_, morans_i = NA_real_,
                      expectation = NA_real_, variance = NA_real_,
                      z = NA_real_, p_value = NA_real_,
                      note = "GADM load failed"))
  }
  morans_i_residuals(admin2_df, polys, k_neighbors = k_neighbors)$morans_i
}
