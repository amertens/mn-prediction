# =============================================================================
# R/national_estimates.R
#
# Compute survey-weighted national prevalence estimates (observed) and
# model-predicted national estimates for comparison.
# =============================================================================

#' Compute national prevalence estimates: survey-weighted observed vs predicted
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param sl_fit Output from fit_sl_models()
#' @param cc Country config
#' @param oc Outcome config
#' @return data.frame with one row: country, outcome, observed/predicted national estimates
compute_national_estimates <- function(outcome_data, sl_fit, cc, oc) {

  d <- outcome_data$data
  if (nrow(d) == 0) return(NULL)

  # Survey weights
  wt_col <- cc$weight_col
  has_weights <- !is.null(wt_col) && wt_col %in% colnames(d) &&
                 any(!is.na(d[[wt_col]]) & d[[wt_col]] > 0)
  if (has_weights) {
    wts <- as.numeric(d[[wt_col]])
    wts[is.na(wts) | wts <= 0] <- 1
  } else {
    wts <- rep(1, nrow(d))
  }

  result <- data.frame(
    country  = cc$country,
    outcome  = oc$tag,
    n        = nrow(d),
    stringsAsFactors = FALSE
  )

  # ── Survey-weighted observed prevalence ──
  # Binary outcome
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    y_bin <- as.numeric(d[[oc$binary]])
    valid <- !is.na(y_bin)
    result$obs_prev <- stats::weighted.mean(y_bin[valid], wts[valid])
    # 2026-06-23 (M5): observed COUNT of deficient individuals (raw, unweighted).
    # Was weighted-prevalence * n = EXPECTED events, which mislabels a model-style
    # quantity as an observed count. `result$exp_n_events` keeps the weighted estimate.
    result$obs_n_events <- sum(y_bin[valid] == 1, na.rm = TRUE)
    result$exp_n_events <- stats::weighted.mean(y_bin[valid], wts[valid]) * sum(valid)

    # Survey design SE (using srvyr if possible)
    svy_se <- tryCatch({
      options(survey.lonely.psu = "adjust")
      svy_df <- d
      svy_df$`.y` <- y_bin
      svy_df$`.wt` <- wts

      if (!is.null(cc$strata_col) && cc$strata_col %in% colnames(svy_df)) {
        des <- srvyr::as_survey_design(
          svy_df,
          ids    = !!rlang::sym(cc$psu_col),
          strata = !!rlang::sym(cc$strata_col),
          weights = `.wt`,
          nest = TRUE
        )
      } else {
        des <- srvyr::as_survey_design(
          svy_df,
          ids     = !!rlang::sym(cc$psu_col),
          weights = `.wt`
        )
      }

      sm <- des %>%
        dplyr::summarise(
          p = srvyr::survey_mean(`.y`, vartype = c("se", "ci"), na.rm = TRUE)
        )
      list(se = sm$p_se, lo = sm$p_low, hi = sm$p_upp)
    }, error = function(e) list(se = NA_real_, lo = NA_real_, hi = NA_real_))

    result$obs_se  <- svy_se$se
    result$obs_lo  <- svy_se$lo
    result$obs_hi  <- svy_se$hi
  }

  # ── Model-predicted national prevalence ──
  fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
  if (!is.null(fit)) {
    preds <- as.numeric(fit$res$yhat_full)

    if (!is.null(sl_fit$bin_fit)) {
      # Binary model: predictions are P(deficient)
      result$pred_prev <- stats::weighted.mean(preds, wts)
    } else {
      # Continuous model: 2026-06-23 (H3) integrate over the empirical OOF
      # residual distribution (survey-weighted) instead of thresholding the
      # point predictions, which is a Jensen-biased plug-in. preds and the
      # continuous outcome are on the same modelled scale (so is oc$cutoff).
      # Falls back to the plug-in if residuals cannot be aligned.
      dir <- oc$cutoff_dir %||% "<"; cut <- oc$cutoff
      ycont <- suppressWarnings(as.numeric(d[[oc$continuous]]))
      if (length(ycont) == length(preds) &&
          sum(is.finite(ycont) & is.finite(preds)) >= 10) {
        okc   <- is.finite(ycont) & is.finite(preds)
        resid <- ycont[okc] - preds[okc]
        p_lt  <- stats::ecdf(resid)(cut - preds)         # P(resid < cut - yhat_i)
        pcont <- if (dir == "<") p_lt else 1 - p_lt
        result$pred_prev <- stats::weighted.mean(pcont, wts, na.rm = TRUE)
      } else {
        pred_bin <- as.numeric(apply_threshold(preds, oc$cutoff, oc$cutoff_dir))
        result$pred_prev <- stats::weighted.mean(pred_bin, wts)
      }
    }

    result$model_type <- if (!is.null(sl_fit$bin_fit)) "binary" else "continuous"
  }

  # Difference
  if ("obs_prev" %in% names(result) && "pred_prev" %in% names(result)) {
    result$abs_diff <- abs(result$pred_prev - result$obs_prev)
    result$pct_diff <- ifelse(result$obs_prev > 0.01,
                               100 * abs(result$pred_prev - result$obs_prev) / result$obs_prev,
                               NA_real_)
  }

  cat(sprintf("[national] %s %s: observed = %.1f%% (%.1f-%.1f%%), predicted = %.1f%%, diff = %.1f pp\n",
              cc$country, oc$tag,
              100 * (result$obs_prev %||% NA),
              100 * (result$obs_lo %||% NA),
              100 * (result$obs_hi %||% NA),
              100 * (result$pred_prev %||% NA),
              100 * (result$abs_diff %||% NA)))

  result
}
