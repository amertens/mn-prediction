# =============================================================================
# R/utils.R
#
# Utility functions shared across all tutorial scripts.
# All functions are pure (no side effects beyond returned value).
# =============================================================================

# ── Performance metrics ───────────────────────────────────────────────────────

#' Compute performance metrics for continuous or binary predictions.
#'
#' @param obs    Numeric vector of observed outcomes.
#' @param pred   Numeric vector of predictions (probabilities for binary).
#' @param type   "continuous" or "binary".
#' @param weights Optional numeric weights (e.g., survey weights).
#' @return Named list of metrics.
compute_metrics <- function(obs, pred, type = c("continuous", "binary"),
                             weights = NULL) {
  type <- match.arg(type)
  keep <- !is.na(obs) & !is.na(pred)
  obs  <- obs[keep];  pred <- pred[keep]
  w    <- if (!is.null(weights)) weights[keep] else rep(1, sum(keep))

  if (type == "continuous") {
    res <- mean(w * (obs - pred)^2) / mean(w)   # weighted MSE
    list(
      rmse        = sqrt(res),
      mae         = weighted.mean(abs(obs - pred), w),
      r2          = 1 - sum(w * (obs - pred)^2) / sum(w * (obs - mean(obs))^2),
      pearson_r   = cor(obs, pred)
    )
  } else {
    # Brier score (weighted)
    brier <- weighted.mean((obs - pred)^2, w)

    # AUC (unweighted; pROC)
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    # PR-AUC
    pr_obj   <- PRROC::pr.curve(scores.class0 = pred[obs == 1],
                                 scores.class1 = pred[obs == 0],
                                 curve = FALSE)
    prauc_val <- pr_obj$auc.integral

    # Calibration slope & intercept (logistic recalibration)
    cal <- tryCatch({
      lp   <- qlogis(pmax(1e-6, pmin(1 - 1e-6, pred)))
      m    <- glm(obs ~ lp, family = binomial(), weights = w)
      list(slope = coef(m)[2], intercept = coef(m)[1])
    }, error = function(e) list(slope = NA, intercept = NA))

    list(
      brier       = brier,
      auc         = auc_val,
      pr_auc      = prauc_val,
      cal_slope   = cal$slope,
      cal_int     = cal$intercept
    )
  }
}

# ── Threshold helper ──────────────────────────────────────────────────────────

#' Convert continuous prediction to binary deficiency indicator.
#'
#' @param x         Numeric vector of predictions.
#' @param cutoff    Scalar threshold.
#' @param direction "less" (deficient if x < cutoff) or "greater".
#' @return Integer 0/1 vector.
apply_threshold <- function(x, cutoff, direction = "less") {
  if (direction == "less") as.integer(x < cutoff) else as.integer(x > cutoff)
}

# ── Calibration ───────────────────────────────────────────────────────────────

#' Logistic recalibration of predicted probabilities.
#'
#' Fits logit(P(D=1)) = beta_0 + beta_1 * logit(p_hat) on a calibration
#' dataset and returns recalibrated probabilities for new data.
#'
#' @param pred_cal  Predicted probabilities on calibration set.
#' @param obs_cal   Observed binary outcomes on calibration set.
#' @param pred_new  Predicted probabilities to recalibrate.
#' @return Recalibrated probability vector.
recalibrate_probs <- function(pred_cal, obs_cal, pred_new) {
  lp_cal <- qlogis(pmax(1e-6, pmin(1 - 1e-6, pred_cal)))
  m      <- glm(obs_cal ~ lp_cal, family = binomial())
  lp_new <- qlogis(pmax(1e-6, pmin(1 - 1e-6, pred_new)))
  plogis(coef(m)[1] + coef(m)[2] * lp_new)
}

# ── Loss functions ────────────────────────────────────────────────────────────

#' Weighted log loss (cross-entropy).
#' @param y       Binary outcome vector.
#' @param p       Predicted probability vector.
#' @param weights Optional weight vector.
weighted_log_loss <- function(y, p, weights = NULL) {
  p <- pmax(1e-9, pmin(1 - 1e-9, p))
  w <- if (!is.null(weights)) weights / mean(weights) else rep(1, length(y))
  -mean(w * (y * log(p) + (1 - y) * log(1 - p)))
}

#' Focal loss (Lin et al., 2017).
#' Downweights easy examples to focus learning on hard cases.
#'
#' @param y     Binary outcome vector.
#' @param p     Predicted probability vector.
#' @param gamma Focusing parameter (>= 0; 0 = standard log loss).
#' @param alpha Class-balance weight for positive class.
focal_loss <- function(y, p, gamma = 2, alpha = 0.5) {
  p  <- pmax(1e-9, pmin(1 - 1e-9, p))
  pt <- ifelse(y == 1, p, 1 - p)
  at <- ifelse(y == 1, alpha, 1 - alpha)
  mean(-at * (1 - pt)^gamma * log(pt))
}

#' Tversky loss (generalized Dice, controls FP vs FN trade-off).
#'
#' @param y     Binary outcome vector.
#' @param p     Predicted probability vector.
#' @param alpha FP penalty weight.
#' @param beta  FN penalty weight.
tversky_loss <- function(y, p, alpha = 0.5, beta = 0.5) {
  eps <- 1e-6
  tp  <- sum(y * p)
  fp  <- sum((1 - y) * p)
  fn  <- sum(y * (1 - p))
  1 - (tp + eps) / (tp + alpha * fp + beta * fn + eps)
}

# ── Prevalence helpers ────────────────────────────────────────────────────────

#' Summarise Admin1 prevalence (unweighted mean of binary outcome).
#' @param df         Data frame with outcome and admin1 columns.
#' @param admin1_col Column name for Admin1 identifier.
#' @param outcome    Column name for binary or continuous outcome.
#' @return Data frame: admin1, n, mean_prev.
summarise_admin1 <- function(df, admin1_col, outcome) {
  df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(
      n         = dplyr::n(),
      mean_prev = mean(.data[[outcome]], na.rm = TRUE),
      .groups   = "drop"
    )
}

#' Compute Admin1 bias, RMSE, and coverage for a vector of estimates.
#' @param true_prev  Named numeric vector of true Admin1 prevalences.
#' @param est_prev   Named numeric vector of estimated Admin1 prevalences.
#' @param ci_lo      Optional lower 95% CI bound (named numeric).
#' @param ci_hi      Optional upper 95% CI bound (named numeric).
#' @return Data frame with columns: admin1, true, est, bias, abs_error.
eval_prevalence <- function(true_prev, est_prev, ci_lo = NULL, ci_hi = NULL) {
  admins <- intersect(names(true_prev), names(est_prev))
  df <- data.frame(
    admin1    = admins,
    true      = true_prev[admins],
    est       = est_prev[admins],
    bias      = est_prev[admins] - true_prev[admins],
    abs_error = abs(est_prev[admins] - true_prev[admins])
  )
  if (!is.null(ci_lo) && !is.null(ci_hi)) {
    df$ci_lo   <- ci_lo[admins]
    df$ci_hi   <- ci_hi[admins]
    df$covered <- true_prev[admins] >= ci_lo[admins] &
                  true_prev[admins] <= ci_hi[admins]
  }
  df
}

# ── Plotting helpers ──────────────────────────────────────────────────────────

#' Scatter plot of observed vs predicted values.
#' @param obs   Numeric vector of observed outcomes.
#' @param pred  Numeric vector of predictions.
#' @param label Character label for the plot subtitle.
#' @return ggplot object.
make_scatter <- function(obs, pred, label = "") {
  df <- data.frame(obs = obs, pred = pred)
  r2 <- round(1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2), 3)
  ggplot2::ggplot(df, ggplot2::aes(x = obs, y = pred)) +
    ggplot2::geom_point(alpha = 0.4, colour = "steelblue", size = 1.5) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                          colour = "grey40") +
    ggplot2::annotate("text", x = min(obs), y = max(pred),
                       label = sprintf("R² = %.3f", r2),
                       hjust = 0, vjust = 1, size = 3.5) +
    ggplot2::labs(x = "Observed", y = "Predicted", subtitle = label) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Calibration plot for binary predictions.
#' @param obs   Binary observed outcomes.
#' @param pred  Predicted probabilities.
#' @param n_bins Number of quantile bins.
#' @return ggplot object.
calibration_plot <- function(obs, pred, n_bins = 10) {
  breaks <- quantile(pred, probs = seq(0, 1, length.out = n_bins + 1),
                     na.rm = TRUE)
  breaks <- unique(breaks)
  bin    <- cut(pred, breaks = breaks, include.lowest = TRUE)
  cal_df <- data.frame(obs = obs, pred = pred, bin = bin) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(pred),
      mean_obs  = mean(obs),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  ggplot2::ggplot(cal_df, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                          colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = "steelblue", alpha = 0.8) +
    ggplot2::geom_line(colour = "steelblue", linewidth = 0.8) +
    ggplot2::scale_x_continuous("Mean predicted probability", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("Observed fraction",          limits = c(0, 1)) +
    ggplot2::labs(title = "Calibration plot",
                   subtitle = "Perfect calibration: points on 45° line") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Pretty metrics table (formatted data frame for knitr::kable).
#' @param metrics_list Named list of named-list metrics (one entry per model).
#' @return Data frame.
make_perf_table <- function(metrics_list) {
  do.call(rbind, lapply(names(metrics_list), function(nm) {
    m  <- metrics_list[[nm]]
    df <- as.data.frame(t(unlist(m)))
    df$model <- nm
    df
  })) %>%
    dplyr::select(model, dplyr::everything())
}

# ── Miscellaneous ─────────────────────────────────────────────────────────────

#' Soft-clip a vector to [lo, hi] to avoid log(0) issues.
clip_probs <- function(p, lo = 1e-6, hi = 1 - 1e-6) {
  pmax(lo, pmin(hi, p))
}

#' Exponential spatial covariance matrix.
#' @param dists  Distance matrix (numeric).
#' @param sigma2 Marginal variance.
#' @param rho    Range parameter (larger = longer-range correlation).
exp_cov <- function(dists, sigma2 = 1, rho = 0.3) {
  sigma2 * exp(-dists / rho)
}

#' Print a named summary vector cleanly.
print_summary <- function(x, digits = 4) {
  cat(paste(sprintf("  %-20s %s", names(x),
                    round(unlist(x), digits)), collapse = "\n"), "\n")
}
