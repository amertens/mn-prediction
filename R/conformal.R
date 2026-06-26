# =============================================================================
# R/conformal.R
#
# Conformal prediction intervals for SuperLearner predictions.
#
# Uses the cross-validated residuals from the SL fit — no refitting needed.
# This replaces the expensive bootstrap (B × SL refits) with a method that
# is essentially free to compute and provides distribution-free coverage
# guarantees under exchangeability.
#
# Two approaches:
#   1. Split conformal: constant-width intervals based on quantile of |residuals|
#   2. Locally-adaptive conformal: width varies with predicted probability
#      (wider intervals where the model is less certain)
#
# For aggregated (Admin-1/national) estimates, we propagate individual
# intervals to regional means using the delta method.
#
# Reference: Barber, Candes, Ramdas, Tibshirani (2021).
#            "Predictive Inference with the Jackknife+"
# =============================================================================


#' Compute conformal prediction intervals from cross-validated SL predictions
#'
#' Uses the existing CV predictions (each observation predicted by a model
#' that didn't see it) to construct valid prediction intervals.
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param sl_fit Output from fit_sl_models()
#' @param cc Country config
#' @param oc Outcome config
#' @param alpha Miscoverage rate (default 0.05 for 95% intervals)
#' @param rho Assumed intra-area correlation of prediction errors used when
#'   aggregating individual interval widths to an Admin-1/national MEAN CI
#'   (2026-06-23 H2 fix). The previous code implicitly used rho = 0 (errors
#'   independent), which shrank the aggregated CI by ~1/sqrt(m) and ignored
#'   spatial/sampling error correlation, badly understating mean uncertainty.
#'   Default 0.5 (moderate, conservative); rho = 0 reproduces the old behaviour,
#'   rho = 1 applies no variance reduction. NOTE: the aggregated interval is a
#'   CI for the AREA MEAN prevalence, not a predictive interval for new units.
#' @return list with admin1_ci, national_ci (same structure as bootstrap output)
compute_conformal_ci <- function(outcome_data, sl_fit, cc, oc, alpha = 0.05,
                                  rho = 0.5) {

  fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
  if (is.null(fit)) {
    warning("No SL fit available for conformal intervals")
    return(list(admin1_ci = NULL, national_ci = NULL, method = "conformal"))
  }

  d   <- outcome_data$data
  res <- fit$res

  # Extract CV predictions and observed outcomes
  Y    <- as.double(unclass(res$Y))
  yhat <- as.double(unclass(res$yhat_full))
  n    <- length(Y)

  # Survey weights
  wt_col <- cc$weight_col
  if (!is.null(wt_col) && wt_col %in% colnames(d)) {
    wts <- as.numeric(d[[wt_col]])
    wts[is.na(wts) | wts <= 0] <- 1
  } else {
    wts <- rep(1, nrow(d))
  }

  admin1_col <- cc$admin1_col
  admin1_vec <- if (admin1_col %in% colnames(d)) d[[admin1_col]] else rep("National", n)

  use_binary <- !is.null(sl_fit$bin_fit)
  model_type <- if (use_binary) "binary" else "continuous"

  cat(sprintf("\n[conformal] %s — %s | n = %d | model = %s\n",
              cc$country, oc$tag, n, model_type))

  # ── Step 1: Compute conformity scores (CV residuals) ─────────────────
  # For binary outcomes, use absolute residuals |Y - p_hat|
  # For continuous outcomes, use absolute residuals |Y - y_hat|
  residuals <- abs(Y - yhat)

  # ── Step 2: Standard conformal interval ──────────────────────────────
  # The (1-alpha) quantile of |residuals| gives the constant half-width
  # that guarantees marginal coverage >= 1-alpha under exchangeability.
  q_level <- ceiling((1 - alpha) * (n + 1)) / n
  q_level <- min(q_level, 1)  # cap at 1
  conformal_width <- as.numeric(quantile(residuals, q_level))

  cat(sprintf("  Conformal half-width (constant): %.4f (%.1f%%)\n",
              conformal_width, conformal_width * 100))

  # Individual-level intervals
  ci_lo_ind <- pmax(yhat - conformal_width, 0)
  ci_hi_ind <- pmin(yhat + conformal_width, 1)

  # Coverage check. NOTE (2026-06-23, H1): this coverage is computed on the SAME
  # OOF residuals used to set `conformal_width`, so it is ~self-fulfilling
  # (mechanically >= 1-alpha) — a sanity check that the quantile was applied
  # correctly, NOT an independent validity test. Treat as in-sample coverage.
  coverage <- mean(Y >= ci_lo_ind & Y <= ci_hi_ind)
  cat(sprintf("  Individual coverage (in-sample, not held-out): %.1f%% (target: %.1f%%)\n",
              coverage * 100, (1 - alpha) * 100))

  # ── Step 3: Locally-adaptive conformal intervals ─────────────────────
  # Width varies with predicted probability: wider where predictions are
  # more uncertain (near 0.5 for binary, or in sparse regions for continuous).
  # Scale residuals by a local measure of difficulty.
  if (use_binary) {
    # For binary: scale by predicted variance p*(1-p)
    pred_var <- pmax(yhat * (1 - yhat), 0.001)
    normalized_residuals <- residuals / sqrt(pred_var)
  } else {
    # For continuous: use MAD of predictions in local neighborhood
    normalized_residuals <- residuals / pmax(median(residuals), 0.001)
  }

  q_adaptive <- as.numeric(quantile(normalized_residuals, q_level))

  if (use_binary) {
    adaptive_width <- q_adaptive * sqrt(pmax(yhat * (1 - yhat), 0.001))
  } else {
    adaptive_width <- q_adaptive * pmax(median(residuals), 0.001)
  }

  ci_lo_adaptive <- pmax(yhat - adaptive_width, 0)
  ci_hi_adaptive <- pmin(yhat + adaptive_width, 1)

  adaptive_coverage <- mean(Y >= ci_lo_adaptive & Y <= ci_hi_adaptive)
  cat(sprintf("  Adaptive coverage (in-sample, not held-out): %.1f%% (target: %.1f%%)\n",
              adaptive_coverage * 100, (1 - alpha) * 100))

  # ── Step 4: Aggregate to Admin-1 ────────────────────────────────────
  # For regional prevalence = weighted mean of individual predictions,
  # the CI is derived from propagating individual intervals.
  #
  # Two approaches:
  #   a) Simple: average of individual CI bounds (conservative)
  #   b) Delta method: SE(mean) = sqrt(sum(w_i^2 * var_i)) / sum(w_i)
  #      where var_i is estimated from the conformal interval width

  admin1_df <- data.frame(
    Admin1   = admin1_vec,
    yhat     = yhat,
    Y        = Y,
    wt       = wts,
    ci_lo    = ci_lo_adaptive,
    ci_hi    = ci_hi_adaptive,
    half_w   = adaptive_width,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Admin1))

  admin1_ci <- admin1_df %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n_obs     = dplyr::n(),
      # Point estimate: survey-weighted mean of predictions
      conf_mean = stats::weighted.mean(yhat, wt, na.rm = TRUE),
      # Observed prevalence (for comparison)
      obs_prev  = stats::weighted.mean(Y, wt, na.rm = TRUE),
      # Delta method SE for weighted mean
      # SE = sqrt(sum(w_i^2 * sigma_i^2)) / sum(w_i)
      # where sigma_i is approximated as half_w_i / qnorm(1 - alpha/2)
      conf_se   = {
        w <- wt / sum(wt)
        sigma_i <- half_w / stats::qnorm(1 - alpha / 2)
        # equicorrelation variance of the weighted mean (H2 fix): rho=0 -> old
        # independence value sum(w^2 sigma^2); rho=1 -> (sum w*sigma)^2.
        sqrt((1 - rho) * sum(w^2 * sigma_i^2) + rho * (sum(w * sigma_i))^2)
      },
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      ci_lo    = pmax(conf_mean - stats::qnorm(1 - alpha / 2) * conf_se, 0),
      ci_hi    = pmin(conf_mean + stats::qnorm(1 - alpha / 2) * conf_se, 1),
      ci_width = ci_hi - ci_lo
    )

  # Rename to match bootstrap output format (boot_mean → conf_mean)
  admin1_ci <- admin1_ci %>%
    dplyr::rename(boot_mean = conf_mean) %>%
    dplyr::select(Admin1, n_boot = n_obs, boot_mean, ci_lo, ci_hi, ci_width)

  cat(sprintf("  Admin-1 CI widths: median = %.1f pp, range = [%.1f, %.1f] pp\n",
              100 * median(admin1_ci$ci_width),
              100 * min(admin1_ci$ci_width),
              100 * max(admin1_ci$ci_width)))

  # ── Step 5: National CI ──────────────────────────────────────────────
  nat_mean <- stats::weighted.mean(yhat, wts)
  nat_se <- {
    w <- wts / sum(wts)
    sigma_i <- adaptive_width / stats::qnorm(1 - alpha / 2)
    # equicorrelation variance (H2 fix); see `rho` in the function header.
    sqrt((1 - rho) * sum(w^2 * sigma_i^2) + rho * (sum(w * sigma_i))^2)
  }

  national_ci <- data.frame(
    outcome   = oc$tag,
    country   = cc$country,
    n_reps    = n,  # not bootstrap reps, but n observations used
    boot_mean = nat_mean,
    ci_lo     = pmax(nat_mean - stats::qnorm(1 - alpha / 2) * nat_se, 0),
    ci_hi     = pmin(nat_mean + stats::qnorm(1 - alpha / 2) * nat_se, 1),
    method    = "conformal_adaptive"
  )

  cat(sprintf("  National: %.1f%% (%.1f%%, %.1f%%)\n",
              100 * nat_mean, 100 * national_ci$ci_lo, 100 * national_ci$ci_hi))

  # ── Step 6: Admin-2 conformal CIs ──────────────────────────────────────
  # Aggregate individual conformal half-widths to Admin-2 level using the
  # same delta method as Admin-1 above. This produces district-specific
  # uncertainty rather than broadcasting a single Admin-1 CI to all
  # constituent districts.
  admin2_ci <- NULL
  admin2_col <- cc$admin2_col %||% "Admin2"
  admin2_vec <- if (admin2_col %in% colnames(d)) d[[admin2_col]]
                else if ("Admin2" %in% colnames(d)) d[["Admin2"]] else NULL
  if (!is.null(admin2_vec)) {
    admin2_df <- data.frame(
      Admin2 = admin2_vec,
      yhat   = yhat,
      Y      = Y,
      wt     = wts,
      ci_lo  = ci_lo_adaptive,
      ci_hi  = ci_hi_adaptive,
      half_w = adaptive_width,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(!is.na(Admin2))

    admin2_ci <- admin2_df %>%
      dplyr::group_by(Admin2) %>%
      dplyr::summarise(
        n_obs     = dplyr::n(),
        conf_mean = stats::weighted.mean(yhat, wt, na.rm = TRUE),
        obs_prev  = stats::weighted.mean(Y, wt, na.rm = TRUE),
        conf_se   = {
          w <- wt / sum(wt)
          sigma_i <- half_w / stats::qnorm(1 - alpha / 2)
          # For very small Admin-2 samples, the CI can collapse — apply
          # a minimum SE floor based on the global conformal width.
          se_raw <- sqrt(sum(w^2 * sigma_i^2))
          pmax(se_raw, conformal_width / (4 * stats::qnorm(1 - alpha / 2)))
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        ci_lo    = pmax(conf_mean - stats::qnorm(1 - alpha / 2) * conf_se, 0),
        ci_hi    = pmin(conf_mean + stats::qnorm(1 - alpha / 2) * conf_se, 1),
        ci_width = ci_hi - ci_lo
      ) %>%
      dplyr::rename(boot_mean = conf_mean) %>%
      dplyr::select(Admin2, n_boot = n_obs, boot_mean, ci_lo, ci_hi, ci_width)

    cat(sprintf("  Admin-2 CI widths: median = %.1f pp, range = [%.1f, %.1f] pp (%d districts)\n",
                100 * median(admin2_ci$ci_width),
                100 * min(admin2_ci$ci_width),
                100 * max(admin2_ci$ci_width),
                nrow(admin2_ci)))
  }

  list(
    admin1_ci   = admin1_ci,
    admin2_ci   = admin2_ci,
    national_ci = national_ci,
    method      = "conformal",
    # Extra diagnostics (not used by downstream targets but useful for reports)
    diagnostics = list(
      constant_width    = conformal_width,
      constant_coverage = coverage,
      adaptive_coverage = adaptive_coverage,
      individual_ci     = data.frame(
        yhat = yhat, Y = Y,
        ci_lo = ci_lo_adaptive, ci_hi = ci_hi_adaptive,
        width = ci_hi_adaptive - ci_lo_adaptive
      )
    )
  )
}
