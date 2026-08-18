# =============================================================================
# 14_outcome_resolution_calibration.R
# -----------------------------------------------------------------------------
# Does the transport / calibration EVALUATION change if we model and score at
# CLUSTER level instead of ADMIN-2 aggregated outcomes?
#
# Why this matters. The predictors are area-resolved (~100% of predictor variance
# is between-area), so the effective sample size is the number of AREAS, not
# clusters or individuals. Moving the OUTCOME to a finer resolution does not add
# predictor signal; it mostly adds binomial sampling noise (cluster n is ~9-30
# people, so a cluster "prevalence" moves in steps of 1/n). This script shows
# what that does to (i) ranking, (ii) the transport LEVEL bias, and (iii) the
# calibration that fixes the level.
#
# Three pipelines, each self-contained, per outcome, under LOCO:
#   A. admin2 -> admin2      model and score at admin-2 (the analysis of record)
#   B. cluster -> cluster    model and score at cluster level (noisy target)
#   C. cluster -> admin2     model at cluster, AGGREGATE predictions to admin-2,
#                            score at admin-2 (separates model resolution from
#                            evaluation resolution)
#
# For each we report the raw transport metrics and the metrics after the ORACLE
# national-anchor logit shift (rung 1), so the LEVEL bias and its correction are
# comparable across resolutions.
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "sensitivity/14_outcome_resolution_calibration.R"
# =============================================================================

find_helpers <- function() {
  cands <- c(file.path("simplified subset", "methods", "_helpers.R"),
             file.path("..", "simplified subset", "methods", "_helpers.R"),
             file.path("methods", "_helpers.R"))
  p <- cands[file.exists(cands)][1]
  if (is.na(p)) stop("Could not find _helpers.R from ", getwd())
  p
}
source(find_helpers())
set.seed(1)

OUTCOMES <- c("women_iron", "women_vitA", "child_iron")

logit_shift_to <- function(phat, w, target) {         # rung-1 oracle level shift
  phat <- clip01(phat)
  f <- function(delta) sum(w * plogis(qlogis(phat) + delta)) / sum(w) - target
  if (!is.finite(f(-12)) || !is.finite(f(12)) || f(-12) * f(12) > 0) return(phat)
  plogis(qlogis(phat) + uniroot(f, c(-12, 12))$root)
}
metrics <- function(phat, truth, w) {
  natl  <- sum(w * truth) / sum(w)
  pcal  <- logit_shift_to(phat, w, natl)
  data.frame(
    MAE_raw  = mae_pp(phat,  truth),  MAE_anchored  = mae_pp(pcal, truth),
    bias_raw = bias_pp(phat, truth),  bias_anchored = bias_pp(pcal, truth),
    spearman = spear(phat, truth))
}

run_level <- function(outcome) {
  dc <- prep_outcome(load_level("cluster"), outcome)   # cluster rows
  da <- prep_outcome(load_level("admin2"),  outcome)   # admin-2 rows
  countries <- intersect(unique(dc$country), unique(da$country))
  A <- B <- C <- list()

  for (held in countries) {
    # ---- A. admin2 -> admin2 ------------------------------------------------
    tr <- da$country != held; te <- !tr
    if (sum(te) >= 3 && sum(tr) >= 10) {
      fa <- fit_area(da$X[tr, , drop=FALSE], da$y[tr], da$n[tr], da$ndef[tr], "gaussian")
      A[[held]] <- cbind(held_out = held, n_units = sum(te),
                         metrics(clip01(fa(da$X[te, , drop=FALSE])), da$y[te], da$n[te]))
    }
    # ---- cluster model (shared by B and C) ----------------------------------
    trc <- dc$country != held; tec <- !trc
    if (sum(tec) < 6 || sum(trc) < 10) next
    fc   <- fit_area(dc$X[trc, , drop=FALSE], dc$y[trc], dc$n[trc], dc$ndef[trc], "gaussian")
    phc  <- clip01(fc(dc$X[tec, , drop=FALSE]))
    yc   <- dc$y[tec]; wc <- dc$n[tec]; areac <- dc$area[tec]

    # ---- B. cluster -> cluster ---------------------------------------------
    B[[held]] <- cbind(held_out = held, n_units = sum(tec), metrics(phc, yc, wc))

    # ---- C. cluster -> admin2 (aggregate predictions and truth by weight) ---
    ph_a <- tapply(phc * wc, areac, sum) / tapply(wc, areac, sum)
    y_a  <- tapply(yc  * wc, areac, sum) / tapply(wc, areac, sum)
    w_a  <- tapply(wc, areac, sum)
    C[[held]] <- cbind(held_out = held, n_units = length(ph_a),
                       metrics(as.numeric(ph_a), as.numeric(y_a), as.numeric(w_a)))
  }
  list(A = do.call(rbind, A), B = do.call(rbind, B), C = do.call(rbind, C))
}

summ <- function(tab, label) {
  if (is.null(tab)) return(NULL)
  data.frame(pipeline = label,
             n_units_med  = round(median(tab$n_units)),
             MAE_raw      = round(mean(tab$MAE_raw), 2),
             MAE_anchored = round(mean(tab$MAE_anchored), 2),   # level removed
             abs_bias_raw = round(mean(abs(tab$bias_raw)), 2),  # |level bias|, no cancelling
             spearman     = round(mean(tab$spearman, na.rm = TRUE), 3))
}

for (o in OUTCOMES) {
  r <- run_level(o)
  cat("\n=====================================================================\n")
  cat(sprintf(" %s  — model/score resolution vs transport & calibration\n", o))
  cat("=====================================================================\n")
  tab <- rbind(summ(r$A, "A. admin2 -> admin2"),
               summ(r$B, "B. cluster -> cluster"),
               summ(r$C, "C. cluster -> admin2 (aggregated)"))
  print(tab, row.names = FALSE)
}

cat("\n---------------------------------------------------------------------\n")
cat("How to read it (means over the 4 held-out countries):\n")
cat(" * spearman   = district/cluster RANKING under transport.\n")
cat(" * bias_raw   = the transport LEVEL bias (pp). The thing calibration fixes.\n")
cat(" * MAE_raw    = level error incl. sampling noise; MAE_anchored removes the\n")
cat("               level with the oracle national anchor -> what's LEFT is the\n")
cat("               resolution's irreducible noise floor + ranking error.\n")
cat("Expect: B (cluster) has attenuated spearman and a high MAE_anchored noise\n")
cat("floor; C (aggregate cluster preds to admin-2) largely recovers A; the LEVEL\n")
cat("bias (bias_raw) is similar across resolutions because it is an aggregate\n")
cat("property, so the calibration ladder is unchanged by outcome resolution.\n")
