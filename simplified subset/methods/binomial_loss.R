# =============================================================================
# binomial_loss.R
# -----------------------------------------------------------------------------
# To-do 3b: model area prevalence with a proper binomial likelihood (on the
# deficient / not-deficient counts) instead of a Gaussian link, so predictions
# stay in [0,1] and the mean-variance relationship of a proportion is respected.
#
# Compares, by within-country 5-fold CV, a weighted-Gaussian elastic net against
# a binomial-counts elastic net. Standard modelling; low novelty.
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/binomial_loss.R"
# =============================================================================

this_dir <- tryCatch(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
              error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "simplified subset/methods"
source(file.path(this_dir, "_helpers.R"))

OUTCOME <- "women_iron"
set.seed(1)

df <- load_level("admin2")
d  <- prep_outcome(df, OUTCOME)

cv_one_country <- function(ctry) {
  idx <- which(d$country == ctry)
  if (length(idx) < 12) return(NULL)        # too few areas for a stable 5-fold CV
  X <- d$X[idx, , drop = FALSE]; y <- d$y[idx]; n <- d$n[idx]; ndef <- d$ndef[idx]
  folds <- sample(rep(1:5, length.out = length(idx)))
  pg <- pb <- rep(NA_real_, length(idx))
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    if (sum(tr) < 8) next
    gfit <- tryCatch(fit_area(X[tr, , drop=FALSE], y[tr], n[tr], ndef[tr], "gaussian"),
                     error = function(e) NULL)
    bfit <- tryCatch(fit_area(X[tr, , drop=FALSE], y[tr], n[tr], ndef[tr], "binomial"),
                     error = function(e) NULL)
    if (!is.null(gfit)) pg[te] <- gfit(X[te, , drop = FALSE])  # raw (NOT clipped)
    if (!is.null(bfit)) pb[te] <- bfit(X[te, , drop = FALSE])
  }
  ok <- is.finite(pg) & is.finite(pb)
  if (sum(ok) < 5) return(NULL)
  data.frame(
    country = ctry, n_areas = length(idx),
    RMSE_gauss = rmse_pp(pg[ok], y[ok]),  RMSE_binom = rmse_pp(pb[ok], y[ok]),
    # fraction of predictions that fall OUTSIDE [0,1] (Gaussian can; binomial can't)
    out_of_range_gauss = round(mean(pg[ok] < 0 | pg[ok] > 1), 3),
    out_of_range_binom = round(mean(pb[ok] < 0 | pb[ok] > 1), 3),
    stringsAsFactors = FALSE)
}

rows <- Filter(Negate(is.null), lapply(unique(d$country), cv_one_country))
out  <- do.call(rbind, rows)

cat(sprintf("\nGaussian vs binomial area model, within-country 5-fold CV, %s\n\n", OUTCOME))
print(out, row.names = FALSE)
cat(sprintf("\nMean RMSE (pp): gaussian %.2f  |  binomial %.2f\n",
            mean(out$RMSE_gauss), mean(out$RMSE_binom)))
cat("Binomial keeps every prediction in [0,1] by construction; the Gaussian link\n")
cat("can produce out-of-range values (see out_of_range_gauss). Use binomial /\n")
cat("quasi-binomial as the default for the area-prevalence learners.\n")
