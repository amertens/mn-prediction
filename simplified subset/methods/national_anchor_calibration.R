# =============================================================================
# national_anchor_calibration.R
# -----------------------------------------------------------------------------
# To-do 2: shift transported (leave-one-country-out) predictions so their
# population-weighted mean matches a known national prevalence. This corrects the
# documented transport LEVEL bias while leaving the district RANKING unchanged.
#
# Method: a one-parameter logit shift solved with uniroot (rank-preserving
# benchmarking of small-area estimates to a known aggregate). Standard SAE /
# survey-calibration practice; low methodological novelty.
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/national_anchor_calibration.R"
# =============================================================================

this_dir <- tryCatch(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
              error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "simplified subset/methods"
source(file.path(this_dir, "_helpers.R"))

OUTCOME <- "women_iron"   # any outcome present in mn_admin2.csv
ALPHA   <- NULL           # unused; kept for symmetry with other scripts

df <- load_level("admin2")
d  <- prep_outcome(df, OUTCOME)
countries <- unique(d$country)

# Logit shift: find delta so the population-weighted mean of the shifted
# predictions equals `target`. Monotone in delta, so district ranking is fixed.
logit_shift <- function(phat, w, target) {
  phat <- clip01(phat)
  f <- function(delta) sum(w * plogis(qlogis(phat) + delta)) / sum(w) - target
  lo <- f(-12); hi <- f(12)
  if (is.na(lo) || is.na(hi) || lo * hi > 0) return(phat)  # target unreachable
  delta <- uniroot(f, c(-12, 12))$root
  plogis(qlogis(phat) + delta)
}

rows <- list()
for (held in countries) {
  tr <- d$country != held; te <- !tr
  if (sum(te) < 3 || sum(tr) < 10) next
  predict_fn <- fit_area(d$X[tr, , drop = FALSE], d$y[tr], d$n[tr], d$ndef[tr],
                         family = "gaussian")
  phat  <- clip01(predict_fn(d$X[te, , drop = FALSE]))
  truth <- d$y[te]; w <- d$n[te]
  target <- sum(w * truth) / sum(w)            # stand-in anchor = true national mean
  pcal   <- logit_shift(phat, w, target)
  rows[[held]] <- data.frame(
    held_out = held, n_areas = sum(te),
    target_natl = round(target, 3),
    MAE_before = mae_pp(phat, truth),  MAE_after = mae_pp(pcal, truth),
    bias_before = bias_pp(phat, truth), bias_after = bias_pp(pcal, truth),
    spearman_before = spear(phat, truth), spearman_after = spear(pcal, truth),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)

cat(sprintf("\nNational-anchor calibration (LOCO) for %s\n", OUTCOME))
cat("Anchor = held-out country's true national prevalence (stand-in for VMNIS/BRINDA).\n")
cat("Expect: bias_after ~ 0 by construction; MAE usually drops; Spearman unchanged.\n\n")
print(out, row.names = FALSE)
cat(sprintf("\nMean MAE (pp): before %.2f  ->  after %.2f\n",
            mean(out$MAE_before), mean(out$MAE_after)))
cat(sprintf("Mean |bias| (pp): before %.2f  ->  after %.2f\n",
            mean(abs(out$bias_before)), mean(abs(out$bias_after))))
cat("Ranking identical (spearman_before == spearman_after):",
    isTRUE(all.equal(out$spearman_before, out$spearman_after)), "\n")
