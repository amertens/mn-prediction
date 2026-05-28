# =============================================================================
# sandbox/11_bootstrap_h6.R
#
# Methodological: 80% and 95% bootstrap confidence intervals on LOCO
# Pearson r for the winning H6 combined_filter method.  Resamples
# held-out admin-2 polygons with replacement (B = 500 by default),
# refits the regression on the resampled set, recomputes r.
#
# Honest reporting: with n_test as small as 14 admin-2 polygons
# (Sierra Leone), a point estimate of r = 0.36 carries a wide CI.
# Manuscript should report point + bootstrap CI, not point alone.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

pooled_all <- load_all_pooled()
cross_pooled <- setNames(lapply(pooled_all, function(x) x$pooled_data),
                          names(pooled_all))

B <- 500L

ci_per_outcome <- function(this_outcome) {
  p <- pooled_all[[this_outcome]]
  pd <- p$pooled_data
  rows <- list()
  for (held_out in p$country_names) {
    train <- pd[pd$country != held_out, , drop = FALSE]
    test  <- pd[pd$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    out <- fit_predict_combined_filter(
      train, test, p$common_gee_vars,
      cross_outcome_pooled = cross_pooled, this_outcome = this_outcome,
      min_within_ratio = 0.70, top_k = 12L)
    if (is.null(out) || is.null(out$pred)) next
    obs <- test$svy_prev; pred <- out$pred
    ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
    if (length(obs) < 3) next
    r_point <- suppressWarnings(cor(obs, pred))
    # Bootstrap CI by resampling held-out polygons.
    rs <- numeric(B)
    n  <- length(obs)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      rs[b] <- suppressWarnings(cor(obs[idx], pred[idx]))
    }
    rs <- rs[is.finite(rs)]
    ci80 <- stats::quantile(rs, c(0.1,  0.9 ), na.rm = TRUE)
    ci95 <- stats::quantile(rs, c(0.025,0.975), na.rm = TRUE)
    rows[[held_out]] <- data.frame(
      outcome = this_outcome, held_out = held_out, n_test = n,
      r_point = round(r_point, 3),
      r_lo_80 = round(ci80[1], 3), r_hi_80 = round(ci80[2], 3),
      r_lo_95 = round(ci95[1], 3), r_hi_95 = round(ci95[2], 3),
      B = B, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}

all_ci <- list()
for (oc in OUTCOMES) {
  cat(sprintf("\n[bootstrap] %s ...\n", oc))
  all_ci[[oc]] <- ci_per_outcome(oc)
}
out <- dplyr::bind_rows(all_ci)
print(out, row.names = FALSE)
cat(sprintf("\nmean r = %.3f (median %.3f); 95%% CI half-widths range %.2f - %.2f\n",
            mean(out$r_point, na.rm = TRUE), median(out$r_point, na.rm = TRUE),
            min((out$r_hi_95 - out$r_lo_95)/2, na.rm = TRUE),
            max((out$r_hi_95 - out$r_lo_95)/2, na.rm = TRUE)))

readr::write_csv(out, file.path(SANDBOX_RESULTS, "11_bootstrap_h6.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "11_bootstrap_h6.csv")))
