# =============================================================================
# R_corrected/p7_area.R — P7: design-aware partial-pooling area estimator
#
# A Fay-Herriot area-level model (proper per-area sampling variances) used as a
# FIRST-CLASS estimator of district prevalence, reusing R/benchmark_models.R::
# fit_predict_fh. Falls back to empirical-Bayes shrinkage toward the global mean
# (weight n/(n+k)) if {sae} is unavailable or FH fails — so the target always
# yields a smoothed, design-aware estimate. Output places the partial-pooled
# estimate beside the raw direct survey estimate and the corrected ML estimate.
# =============================================================================

area_partial_pooling_corrected <- function(slfit, svy_admin2, gee_admin2,
                                           n_var = 12L) {
  if (!is.null(slfit$error)) return(NULL)
  if (is.null(svy_admin2) || !"svy_prev" %in% colnames(svy_admin2)) return(NULL)
  sv <- svy_admin2
  if (!"n_svy" %in% colnames(sv)) return(NULL)

  # join area covariates on Admin2 (single country; key conceptually country+Admin2)
  gee_cols <- if (!is.null(gee_admin2)) grep("^gee_", colnames(gee_admin2),
                                             value = TRUE) else character(0)
  area <- sv
  if (length(gee_cols) > 0) {
    g <- gee_admin2[, c("Admin2", gee_cols), drop = FALSE]
    area <- merge(sv, g, by = "Admin2", all.x = TRUE, sort = FALSE)
  }
  # rank covariates by |corr with svy_prev| (area-model in-sample selection)
  vars <- character(0)
  if (length(gee_cols) > 0) {
    cc <- vapply(gee_cols, function(v) {
      r <- tryCatch(suppressWarnings(
        stats::cor(area[[v]], area$svy_prev, use = "pairwise.complete.obs")),
        error = function(e) NA_real_)
      if (is.na(r)) 0 else abs(r)
    }, numeric(1))
    vars <- names(sort(cc, decreasing = TRUE))[seq_len(min(n_var, length(cc)))]
  }

  direct <- area$svy_prev
  direct_se <- if ("svy_prev_se" %in% colnames(area)) area$svy_prev_se else NA_real_
  est <- rep(NA_real_, nrow(area)); method <- "none"

  fh <- tryCatch({
    if (exists("fit_predict_fh") && length(vars) >= 1) {
      tr <- area; tr$svy_prev <- direct
      fit_predict_fh(tr, tr, vars, model_type = "continuous")
    } else NULL
  }, error = function(e) NULL)

  if (!is.null(fh) && !is.null(fh$train_pred) &&
      length(fh$train_pred) == nrow(area)) {
    est <- fh$train_pred; method <- "Fay-Herriot (REML)"
  } else {
    # Empirical-Bayes shrinkage fallback (design-aware via n_svy)
    gm <- stats::weighted.mean(direct, area$n_svy, na.rm = TRUE)
    vt <- direct * (1 - direct) / pmax(area$n_svy, 1)          # sampling var
    bv <- max(stats::var(direct, na.rm = TRUE) - mean(vt, na.rm = TRUE), 1e-6)
    shrink <- bv / (bv + vt)                                    # n/(n+k)-type
    est <- shrink * direct + (1 - shrink) * gm
    method <- "Empirical-Bayes shrinkage (fallback)"
  }

  # attach corrected ML estimate for side-by-side
  ml <- corrected_admin2(slfit, sv)[, c("Admin2", "pred_prev")]
  out <- data.frame(country = slfit$country, outcome = slfit$outcome,
                    Admin2 = area$Admin2,
                    direct_prev = round(direct, 4),
                    direct_se = round(direct_se, 4),
                    partial_pooled_prev = round(est, 4),
                    method = method, stringsAsFactors = FALSE)
  out <- merge(out, ml, by = "Admin2", all.x = TRUE, sort = FALSE)
  names(out)[names(out) == "pred_prev"] <- "ml_prev"
  attr(out, "method") <- method
  out
}
