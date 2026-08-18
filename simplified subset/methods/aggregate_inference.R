# =============================================================================
# aggregate_inference.R
# -----------------------------------------------------------------------------
# The cheap, defensible pieces of "loss-based inference for the aggregate target"
# (see docs/METHODS_TODO_IMPLEMENTATION_PLAN.md, To-do 3c deep dive):
#
#   Option A. Cross-validated aggregate error with bootstrap CIs:
#             (i)  per-area CV error (RMSE pp), area-resampling bootstrap stratified
#                  by country (stable, ~200 units);
#             (ii) between-country national-level error, country-block bootstrap
#                  (honest but ~4 effective units, so wide).
#
#   Option D. Partial-identification band for an unsurveyed country's national
#             prevalence, from the across-country range of the logit level-shift.
#             For each held-out country the plausible shift range is estimated from
#             the OTHER held-out countries (leave-country-out), so the band is
#             genuinely out-of-sample.
#
# Both are honest about the binding constraint: with four countries, between-
# country uncertainty is large and barely estimable.
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/aggregate_inference.R"
# =============================================================================

this_dir <- tryCatch(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
              error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "simplified subset/methods"
source(file.path(this_dir, "_helpers.R"))

OUTCOME <- "women_iron"
B       <- 2000          # bootstrap replicates
ALPHA   <- 0.10          # 90% intervals
set.seed(11)

df <- load_level("admin2")
d  <- prep_outcome(df, OUTCOME)
countries <- unique(d$country)

# Solve the logit shift delta so the weighted mean of shifted predictions hits a
# target national prevalence (NA if the target is unreachable in the bracket).
solve_delta <- function(phat, w, target) {
  phat <- clip01(phat)
  f <- function(delta) sum(w * plogis(qlogis(phat) + delta)) / sum(w) - target
  if (is.na(f(-12)) || is.na(f(12)) || f(-12) * f(12) > 0) return(NA_real_)
  uniroot(f, c(-12, 12))$root
}
shift_mean <- function(phat, w, delta)
  sum(w * plogis(qlogis(clip01(phat)) + delta)) / sum(w)

# ── LOCO pass: per held-out country store predictions, truth, weights, shift ──
loco <- list()
for (held in countries) {
  tr <- d$country != held; te <- !tr
  if (sum(te) < 3 || sum(tr) < 10) next
  predict_fn <- fit_area(d$X[tr, , drop = FALSE], d$y[tr], d$n[tr], d$ndef[tr], "gaussian")
  phat  <- clip01(predict_fn(d$X[te, , drop = FALSE]))
  truth <- d$y[te]; w <- d$n[te]
  natl_pred <- sum(w * phat) / sum(w); natl_true <- sum(w * truth) / sum(w)
  loco[[held]] <- list(
    country = held, phat = phat, truth = truth, w = w,
    natl_pred = natl_pred, natl_true = natl_true,
    natl_abs_err = abs(natl_pred - natl_true),
    delta = solve_delta(phat, w, natl_true))   # shift that would match the truth
}

# =========================================================================
# Option A(i): per-area CV error (pooled RMSE pp) + country-stratified
#              area-resampling bootstrap CI (stable).
# =========================================================================
all_err <- unlist(lapply(loco, function(s) (s$phat - s$truth)))
all_ctry <- unlist(lapply(loco, function(s) rep(s$country, length(s$truth))))
pooled_rmse <- 100 * sqrt(mean(all_err^2))
boot_rmse <- replicate(B, {
  # resample areas WITHIN each country (stratified), preserving country sizes
  idx <- unlist(lapply(split(seq_along(all_err), all_ctry),
                       function(ii) sample(ii, length(ii), replace = TRUE)))
  100 * sqrt(mean(all_err[idx]^2))
})
ci_area <- quantile(boot_rmse, c(ALPHA/2, 1 - ALPHA/2))

# =========================================================================
# Option A(ii): between-country national error + country-block bootstrap
#               (~4 effective units; honest and wide).
# =========================================================================
natl_errs <- vapply(loco, function(s) 100 * s$natl_abs_err, numeric(1))
mean_natl_err <- mean(natl_errs)
boot_block <- replicate(B, mean(sample(natl_errs, length(natl_errs), replace = TRUE)))
ci_block <- quantile(boot_block, c(ALPHA/2, 1 - ALPHA/2))

# =========================================================================
# Option D: partial-identification band for each held-out country's national
#           prevalence, from the leave-country-out range of observed shifts.
# =========================================================================
all_deltas <- vapply(loco, function(s) s$delta, numeric(1))
band_rows <- list()
for (s in loco) {
  others <- all_deltas[names(all_deltas) != s$country]
  others <- others[is.finite(others)]
  if (length(others) < 2) next
  d_lo <- min(others); d_hi <- max(others)
  lo <- shift_mean(s$phat, s$w, d_lo); hi <- shift_mean(s$phat, s$w, d_hi)
  band <- sort(c(lo, hi))
  band_rows[[s$country]] <- data.frame(
    held_out = s$country,
    natl_pred_raw = round(s$natl_pred, 3),
    band_lo = round(band[1], 3), band_hi = round(band[2], 3),
    band_width_pp = round(100 * (band[2] - band[1]), 1),
    natl_true = round(s$natl_true, 3),
    covered = (s$natl_true >= band[1] && s$natl_true <= band[2]),
    stringsAsFactors = FALSE)
}
band_tab <- do.call(rbind, band_rows)

# ── Report ───────────────────────────────────────────────────────────────────
cat(sprintf("\nLoss-based aggregate inference (LOCO) for %s\n", OUTCOME))

cat("\n=== Option A(i): cross-validated per-area error (pooled RMSE) ===\n")
cat(sprintf("  RMSE = %.2f pp   90%% CI [%.2f, %.2f]  (area resampling, country-stratified; ~%d areas)\n",
            pooled_rmse, ci_area[1], ci_area[2], length(all_err)))

cat("\n=== Option A(ii): between-country national-level error ===\n")
cat("  per-country |national pred - national truth| (pp): ",
    paste(sprintf("%s=%.1f", names(natl_errs), natl_errs), collapse = "  "), "\n")
cat(sprintf("  mean = %.2f pp   90%% CI [%.2f, %.2f]  (country-block bootstrap, ~%d units -> wide/unstable)\n",
            mean_natl_err, ci_block[1], ci_block[2], length(natl_errs)))

cat("\n=== Option D: partial-identification band for the national aggregate ===\n")
cat("  Band = apply the leave-country-out range of fitted level-shifts to the raw\n")
cat("  transported prediction (the honest interval when no national anchor exists).\n\n")
print(band_tab, row.names = FALSE)
cat(sprintf("\n  Band covers the truth in %d/%d held-out countries (compare: the naive\n",
            sum(band_tab$covered), nrow(band_tab)))
cat("  bootstrap CI in aggregate_uncertainty.R covered 1/4). Mean band width = ")
cat(sprintf("%.1f pp.\n", mean(band_tab$band_width_pp)))
cat("\n  Caveat: the band is wide because four countries give little information\n")
cat("  about between-country level shifts. With a national anchor (VMNIS/BRINDA),\n")
cat("  the mean is identified and this band collapses to within-country ranking.\n")
