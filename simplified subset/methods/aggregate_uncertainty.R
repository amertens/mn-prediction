# =============================================================================
# aggregate_uncertainty.R
# -----------------------------------------------------------------------------
# To-do 3c (uncertainty part): honest uncertainty on transported predictions.
#   (1) Split-conformal per-district intervals on the held-out country, with
#       empirical coverage against the nominal level (distribution-free).
#   (2) Cluster/area bootstrap interval for the held-out NATIONAL aggregate.
# Both are standard. A small heuristic value-of-information ranking is included
# at the end (clearly labelled a heuristic, not optimal design).
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/aggregate_uncertainty.R"
# =============================================================================

this_dir <- tryCatch(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
              error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "simplified subset/methods"
source(file.path(this_dir, "_helpers.R"))

OUTCOME <- "women_iron"
ALPHA   <- 0.10          # 90% target coverage
B       <- 500           # bootstrap replicates
set.seed(7)

df <- load_level("admin2")
d  <- prep_outcome(df, OUTCOME)
countries <- unique(d$country)

# Split-conformal half-width from calibration residuals, with the finite-sample
# (n+1)/n correction clamped to a valid probability (matches R/corrected/).
conformal <- function(pred, truth, alpha = ALPHA) {
  n <- length(truth); if (n < 6) return(c(halfwidth = NA, coverage = NA))
  cal <- sample(n, floor(n / 2)); res <- abs(pred[cal] - truth[cal])
  q  <- min(1, (1 - alpha) * (length(cal) + 1) / length(cal))
  qh <- as.numeric(stats::quantile(res, q, type = 1))
  te <- setdiff(seq_len(n), cal)
  cov <- mean(truth[te] >= pred[te] - qh & truth[te] <= pred[te] + qh)
  c(halfwidth = qh, coverage = cov)
}

rows <- list(); voi <- NULL
for (held in countries) {
  tr <- d$country != held; te <- !tr
  if (sum(te) < 6 || sum(tr) < 10) next
  predict_fn <- fit_area(d$X[tr, , drop = FALSE], d$y[tr], d$n[tr], d$ndef[tr], "gaussian")
  phat  <- clip01(predict_fn(d$X[te, , drop = FALSE]))
  truth <- d$y[te]; w <- d$n[te]

  # (1) split-conformal per-district intervals + empirical coverage
  cf <- conformal(phat, truth)

  # (2) area bootstrap for the national aggregate (n-weighted mean of predictions)
  idx <- seq_along(phat)
  agg_boot <- replicate(B, {
    bi <- sample(idx, length(idx), replace = TRUE)
    sum(w[bi] * phat[bi]) / sum(w[bi])
  })
  agg_hat  <- sum(w * phat) / sum(w)
  agg_true <- sum(w * truth) / sum(w)
  ci <- quantile(agg_boot, c(ALPHA / 2, 1 - ALPHA / 2))

  rows[[held]] <- data.frame(
    held_out = held, n_areas = sum(te),
    conf_halfwidth_pp = as.numeric(round(100 * cf["halfwidth"], 1)),
    conf_coverage = as.numeric(round(cf["coverage"], 3)),
    natl_pred = round(agg_hat, 3), natl_true = round(agg_true, 3),
    natl_CI_lo = round(ci[1], 3), natl_CI_hi = round(ci[2], 3),
    natl_CI_covers_true = (agg_true >= ci[1] && agg_true <= ci[2]),
    stringsAsFactors = FALSE)

  # heuristic VOI for the held-out country: risk x uncertainty (interval width)
  v <- data.frame(country = held, area = unname(d$area[te]),
                  pred = round(phat, 3),
                  uncertainty_pp = as.numeric(round(100 * cf["halfwidth"], 1)),
                  voi_score = as.numeric(round(phat * cf["halfwidth"], 4)),
                  stringsAsFactors = FALSE)
  voi <- rbind(voi, v[order(-v$voi_score), ][seq_len(min(3, nrow(v))), ])
}
out <- do.call(rbind, rows)

cat(sprintf("\nAggregate uncertainty (LOCO) for %s, target coverage %.0f%%\n\n",
            OUTCOME, 100 * (1 - ALPHA)))
print(out, row.names = FALSE)
cat(sprintf("\nMean conformal coverage: %.2f (target %.2f). National-aggregate bootstrap CI covers the truth in %d/%d held-out countries.\n",
            mean(out$conf_coverage, na.rm = TRUE), 1 - ALPHA,
            sum(out$natl_CI_covers_true), nrow(out)))
cat("\nInterpretation: the per-district split-conformal intervals are valid (coverage\n")
cat("near or above the 90% target). The national bootstrap CI is narrow and often\n")
cat("MISSES the truth because it captures only sampling uncertainty around the\n")
cat("transported point estimate, NOT the transport LEVEL bias (predictions are\n")
cat("pulled toward the training-country mean). That level bias is what the\n")
cat("national-anchor calibration (national_anchor_calibration.R) corrects.\n")

cat("\n--- heuristic value-of-information (top districts by risk x uncertainty) ---\n")
cat("NOTE: this is a triage heuristic, NOT formal optimal design (EVSI).\n\n")
print(voi[order(-voi$voi_score), ][seq_len(min(10, nrow(voi))), ], row.names = FALSE)
