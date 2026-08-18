# =============================================================================
# qa_reporting.R  — standard output QA for transported prevalence surfaces
# -----------------------------------------------------------------------------
# Two reusable checks, prompted by the WFP / MIMI partner materials, that attach
# to ANY predicted surface (within-country or transported):
#
#   1. skill_over_baseline()   -- report every model as improvement over a naive
#      baseline, the way WFP report "normalized difference vs a dummy model".
#      Makes the effective-n story quantitative: if the transported model barely
#      beats "predict the training-country mean", it has learned little that
#      transports, and the skill score says so in one number.
#
#   2. gradient_sanity_check() -- WFP's "sanity check" generalised into an output
#      QA gate: does the predicted surface reproduce an EXPECTED structural
#      gradient (deficiency higher in less-developed / poorer / more-rural areas)?
#      In a data-constrained deployment with no ground truth, this is often the
#      ONLY validation available. Returns a per-country pass/fail.
#
# Both are self-contained (depend only on _helpers.R + glmnet) and are written so
# they can be lifted into R/corrected/ as post-prediction QA. Running this file
# demonstrates them on the simplified subset under leave-one-country-out (LOCO).
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/methods/qa_reporting.R"
# =============================================================================

this_dir <- tryCatch(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
              error = function(e) ".")
if (is.na(this_dir) || this_dir == "") this_dir <- "simplified subset/methods"
source(file.path(this_dir, "_helpers.R"))
set.seed(1)

# =============================================================================
# 1. SKILL OVER BASELINE
# -----------------------------------------------------------------------------
# Skill = 1 - error(model) / error(baseline). Positive = better than baseline,
# 0 = no better, negative = worse. Reported for RMSE and MAE, as a fraction and
# a percent, alongside the raw errors so nothing is hidden.
#
# Two baselines make the transport story legible:
#   * TRANSPORT skill  (baseline = pooled TRAINING-country weighted mean): the
#     honest, deployable "did we learn anything that transports beyond the mean?"
#   * RANKING skill    (baseline = the held-out country's OWN weighted mean,
#     model first mean-anchored to remove the level bias): "within the country,
#     does the surface rank areas better than a flat map?" This isolates the part
#     of the signal that survives the documented level bias.
# =============================================================================

# One-parameter logit shift so the weighted mean of `p` equals `target` (reused
# from the calibration work; rank-preserving, used only for the RANKING baseline).
.shift_to_mean <- function(p, w, target) {
  p <- clip01(p)
  f <- function(d) sum(w * plogis(qlogis(p) + d)) / sum(w) - target
  if (!is.finite(f(-12)) || !is.finite(f(12)) || f(-12) * f(12) > 0) return(p)
  plogis(qlogis(p) + uniroot(f, c(-12, 12))$root)
}

skill_over_baseline <- function(pred, truth, weight = NULL,
                                baseline = NULL) {
  n <- length(truth)
  if (is.null(weight)) weight <- rep(1, n)
  # default baseline = weighted mean of the truth passed in (flat map)
  if (is.null(baseline)) baseline <- rep(sum(weight * truth) / sum(weight), n)
  if (length(baseline) == 1L) baseline <- rep(baseline, n)

  wrmse <- function(a) sqrt(sum(weight * (a - truth)^2) / sum(weight))
  wmae  <- function(a) sum(weight * abs(a - truth)) / sum(weight)
  rm_m <- wrmse(pred); rm_b <- wrmse(baseline)
  ma_m <- wmae(pred);  ma_b <- wmae(baseline)
  data.frame(
    RMSE_model = round(100 * rm_m, 2), RMSE_base = round(100 * rm_b, 2),
    skill_RMSE = round(1 - rm_m / rm_b, 3),
    MAE_model  = round(100 * ma_m, 2), MAE_base  = round(100 * ma_b, 2),
    skill_MAE  = round(1 - ma_m / ma_b, 3))
}

# =============================================================================
# 2. GRADIENT SANITY CHECK  (output QA gate)
# -----------------------------------------------------------------------------
# For each country, split areas into the bottom vs top tertile of a structural
# gradient variable (development / urbanicity / wealth proxy) and measure the
# PREDICTED contrast = mean(pred | bottom) - mean(pred | top). Compare its SIGN
# to the expected direction. `expected_sign = +1` encodes "deficiency higher in
# the LESS-developed (bottom) tertile", the usual public-health expectation.
#
# pass = the predicted gradient runs in the expected direction. A FALSE is a red
# flag that the transported surface is not reproducing known structure and should
# not be trusted for targeting in that country.
# =============================================================================

gradient_sanity_check <- function(pred, gradient, country = NULL,
                                  weight = NULL, expected_sign = +1,
                                  min_per_stratum = 2) {
  n <- length(pred)
  if (is.null(weight))  weight  <- rep(1, n)
  if (is.null(country)) country <- rep("all", n)
  wm <- function(v, w) if (sum(w) > 0) sum(w * v) / sum(w) else NA_real_

  rows <- lapply(unique(country), function(c) {
    ic <- country == c
    g  <- gradient[ic]; p <- pred[ic]; w <- weight[ic]
    lo <- suppressWarnings(stats::quantile(g, 1/3, na.rm = TRUE, names = FALSE))
    hi <- suppressWarnings(stats::quantile(g, 2/3, na.rm = TRUE, names = FALSE))
    bot <- is.finite(g) & g <= lo
    top <- is.finite(g) & g >= hi
    if (sum(bot) < min_per_stratum || sum(top) < min_per_stratum)
      return(data.frame(country = c, n_areas = sum(ic),
                        pred_bottom = NA, pred_top = NA,
                        contrast = NA, expected_sign = expected_sign,
                        pass = NA))
    pb <- wm(p[bot], w[bot]); pt <- wm(p[top], w[top]); contr <- pb - pt
    data.frame(country = c, n_areas = sum(ic),
               pred_bottom = round(pb, 3), pred_top = round(pt, 3),
               contrast = round(contr, 3), expected_sign = expected_sign,
               pass = sign(contr) == sign(expected_sign))
  })
  out <- do.call(rbind, rows)
  attr(out, "pass_rate") <- mean(out$pass, na.rm = TRUE)
  out
}

# =============================================================================
# DEMO: both checks on the simplified subset under LOCO
# =============================================================================

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  df       <- load_level("admin2")
  OUTCOMES <- c("women_iron", "women_vitA", "child_iron")
  GRAD_VAR <- "human_modification_index"   # development/urbanicity gradient

  cat("=====================================================================\n")
  cat(" QA report: skill-over-baseline + gradient sanity check (LOCO)\n")
  cat(" Gradient variable:", GRAD_VAR, "(expected: deficiency higher in the\n")
  cat("   LESS-developed / bottom tertile -> expected_sign = +1)\n")
  cat("=====================================================================\n")

  for (outcome in OUTCOMES) {
    d  <- prep_outcome(df, outcome)
    countries <- unique(d$country)
    skill_rows <- list(); grad_all <- list()

    for (held in countries) {
      tr <- d$country != held; te <- !tr
      if (sum(te) < 6 || sum(tr) < 10) next
      predict_fn <- fit_area(d$X[tr, , drop = FALSE], d$y[tr], d$n[tr],
                             d$ndef[tr], family = "gaussian")
      phat  <- clip01(predict_fn(d$X[te, , drop = FALSE]))
      truth <- d$y[te]; w <- d$n[te]

      # --- skill: transport (vs training mean) and ranking (vs own mean) -----
      train_mean <- sum(d$n[tr] * d$y[tr]) / sum(d$n[tr])
      own_mean   <- sum(w * truth) / sum(w)
      transport  <- skill_over_baseline(phat, truth, w,
                                        baseline = rep(train_mean, length(phat)))
      p_anchored <- .shift_to_mean(phat, w, own_mean)   # remove level bias
      ranking    <- skill_over_baseline(p_anchored, truth, w,
                                        baseline = rep(own_mean, length(phat)))
      skill_rows[[held]] <- data.frame(
        held_out = held, n_areas = sum(te),
        RMSE_model = transport$RMSE_model,
        skill_vs_trainmean = transport$skill_RMSE,   # deployable transport skill
        skill_vs_ownmean   = ranking$skill_RMSE)     # pure ranking skill

      # --- gradient sanity check on the transported surface -----------------
      grad_all[[held]] <- gradient_sanity_check(
        phat, d$X[te, GRAD_VAR], country = rep(held, length(phat)),
        weight = w, expected_sign = +1)
    }

    cat(sprintf("\n----- %s -----\n", outcome))
    cat("Skill over baseline (RMSE-based; 1 = perfect, 0 = no better, <0 = worse):\n")
    print(do.call(rbind, skill_rows), row.names = FALSE)
    gtab <- do.call(rbind, grad_all)
    cat(sprintf("Gradient sanity check (pass rate %.0f%%):\n",
                100 * mean(gtab$pass, na.rm = TRUE)))
    print(gtab[, c("country","pred_bottom","pred_top","contrast","pass")],
          row.names = FALSE)
  }

  cat("\n---------------------------------------------------------------------\n")
  cat("Reading it:\n")
  cat(" * skill_vs_trainmean <= ~0  => the transported model barely beats guessing\n")
  cat("   the training mean: little transportable signal (the effective-n story).\n")
  cat(" * skill_vs_ownmean is the RANKING skill left after removing the level bias;\n")
  cat("   positive here but near-zero for trainmean = 'ranks a bit, level is wrong'.\n")
  cat(" * gradient pass = FALSE flags a surface that does NOT reproduce the expected\n")
  cat("   development->deficiency gradient and should not be trusted for targeting.\n")
}
