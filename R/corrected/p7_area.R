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

# In-fold covariate ranking by |corr with svy_prev| (uses ONLY the rows in df).
.fh_rank_vars <- function(df, cols, k) {
  if (length(cols) == 0) return(character(0))
  cc <- vapply(cols, function(v) {
    r <- tryCatch(suppressWarnings(
      stats::cor(df[[v]], df$svy_prev, use = "pairwise.complete.obs")),
      error = function(e) NA_real_)
    if (is.na(r)) 0 else abs(r)
  }, numeric(1))
  names(sort(cc, decreasing = TRUE))[seq_len(min(k, length(cc)))]
}

area_partial_pooling_corrected <- function(slfit, svy_admin2, gee_admin2,
                                           n_var = 12L) {
  if (!is.null(slfit$error)) return(NULL)
  if (is.null(svy_admin2) || !"svy_prev" %in% colnames(svy_admin2)) return(NULL)
  sv <- svy_admin2
  if (!"n_svy" %in% colnames(sv)) return(NULL)

  # join area covariates on Admin2 (single country; key conceptually country+Admin2)
  gee_cols <- if (!is.null(gee_admin2)) area_covariate_cols(gee_admin2) else character(0)
  area <- sv
  if (length(gee_cols) > 0) {
    g <- gee_admin2[, c("Admin2", gee_cols), drop = FALSE]
    area <- merge(sv, g, by = "Admin2", all.x = TRUE, sort = FALSE)
  }
  direct    <- area$svy_prev
  direct_se <- if ("svy_prev_se" %in% colnames(area)) area$svy_prev_se else NA_real_
  nA        <- nrow(area)

  # ── (reference, NOT predictive) in-sample EBLUP ceiling ────────────────────
  # Covariates selected on the FULL surveyed set and FH fit with train == test.
  # Retained only as a LABELLED ceiling so the in-sample -> out-of-sample optimism
  # gap is visible (cf. p9 `sae_glmm_insample`). Its correlation with `direct`
  # must NEVER be reported as predictive skill (it selects + fits + evaluates on
  # the same areas). Issue 1 of CORRECTED_PIPELINE_PATCH_PLAN.md.
  vars_all <- .fh_rank_vars(area, gee_cols, n_var)
  est_in <- rep(NA_real_, nA)
  fh_in <- tryCatch({
    if (exists("fit_predict_fh") && length(vars_all) >= 1)
      fit_predict_fh(area, area, vars_all, model_type = "continuous") else NULL
  }, error = function(e) NULL)
  if (!is.null(fh_in) && !is.null(fh_in$train_pred) &&
      length(fh_in$train_pred) == nA) est_in <- fh_in$train_pred

  # ── (headline) honest leave-one-AREA-out FH ────────────────────────────────
  # For each held-out Admin2 i: rank covariates and fit FH on the OTHER areas
  # only, then predict area i from its covariates (synthetic-only `fh$pred`). No
  # row uses its own svy_prev, so cor(partial_pooled_prev, direct) is honest
  # out-of-sample skill. Per-fold fallback (FH NULL or degenerate — see Issue 4
  # guard) is the training grand mean, the best covariate-free out-of-sample guess.
  est     <- rep(NA_real_, nA)
  used_fh <- logical(nA)
  for (i in seq_len(nA)) {
    tr <- area[-i, , drop = FALSE]
    te <- area[ i, , drop = FALSE]
    vars_i <- .fh_rank_vars(tr, gee_cols, n_var)
    fh <- tryCatch({
      if (exists("fit_predict_fh") && length(vars_i) >= 1)
        fit_predict_fh(tr, te, vars_i, model_type = "continuous") else NULL
    }, error = function(e) NULL)
    if (!is.null(fh) && !is.null(fh$pred) && length(fh$pred) == 1 &&
        is.finite(fh$pred)) {
      est[i] <- fh$pred; used_fh[i] <- TRUE
    } else {
      est[i] <- stats::weighted.mean(tr$svy_prev, tr$n_svy, na.rm = TRUE)
    }
  }
  method <- if (all(used_fh)) "leave-area-out FH (out-of-sample)"
            else if (any(used_fh))
              sprintf("leave-area-out FH (out-of-sample); %d/%d folds grand-mean fallback",
                      sum(!used_fh), nA)
            else "leave-area-out grand-mean (FH unavailable)"

  # attach corrected ML estimate for side-by-side
  ml <- corrected_admin2(slfit, sv)[, c("Admin2", "pred_prev")]
  out <- data.frame(country = slfit$country, outcome = slfit$outcome,
                    Admin2 = area$Admin2,
                    direct_prev = round(direct, 4),
                    direct_se = round(direct_se, 4),
                    partial_pooled_prev  = round(est, 4),       # honest OOS (headline)
                    insample_pooled_prev = round(est_in, 4),    # labelled ceiling — NOT skill
                    loao_used_fh = used_fh,
                    estimate_type = "leave-area-out FH (out-of-sample)",
                    method = method, stringsAsFactors = FALSE)
  out <- merge(out, ml, by = "Admin2", all.x = TRUE, sort = FALSE)
  names(out)[names(out) == "pred_prev"] <- "ml_prev"
  attr(out, "method") <- method
  out
}
