# =============================================================================
# 13_gradient_anchored_calibration.R
# -----------------------------------------------------------------------------
# Calibration ladder, RUNG 2: correct the transport LEVEL bias when NO national
# prevalence is available for the target country.
#
# Context. Under leave-one-country-out (LOCO) transport, predictions keep the
# district RANKING but are pulled toward the training-country mean (documented
# LEVEL bias). Rung 1 (national_anchor_calibration.R) fixes this with a one-
# parameter logit shift to a KNOWN national prevalence. But the operational use
# case is a country with no biomarker survey at all -- so no national anchor.
#
# Idea (adapted from WFP / MIMI, Voukelatou et al. 2025, Sci Rep). There they
# calibrate a transported model with no ground truth by tuning one parameter so
# the output reproduces an expected STRUCTURAL GRADIENT (deficiency falls as
# wealth / urbanicity rises). We have no wealth or urban/rural field in the
# simplified subset, so we use the available DEVELOPMENT / URBANICITY gradient
# (Global Human Modification index, built-up surface, cropland population
# density) as the structural axis, and we anchor the LEVEL to a REFERENCE
# STRATUM rather than the whole distribution:
#
#   * Reference stratum = the most-developed tertile of areas. Its deficiency
#     level is the most cross-country-comparable part of the distribution
#     (development converges at the top; the national mean is dominated by each
#     country's idiosyncratic rural burden). So we BORROW the top-tertile level
#     from the surveyed (training) countries and shift the transported target
#     predictions to match it.
#   * The shift is a single logit delta -> strictly rank-preserving, exactly
#     like rung 1. Only the anchor changes: a borrowed subgroup level instead of
#     a known national mean.
#
# Deployability. Unlike rung 1's oracle (which here uses the TARGET country's own
# true mean), the gradient anchor borrows only from TRAINING countries, so it is
# usable with zero target-country biomarker data. We report all three -- raw,
# gradient-anchored (deployable), national-anchored (oracle upper bound) -- so
# the question is: how much of the oracle's level correction does the deployable
# gradient anchor recover?
#
# HONESTY FLAGS.
#   * The logit shift / benchmarking is standard, low-novelty SAE.
#   * NOVEL-IN-THIS-SETTING and NOT YET A MANUSCRIPT CLAIM (take to Alan):
#     anchoring the transported LEVEL to a borrowed structural-gradient stratum
#     is a non-standard identification assumption. It is only valid if (a) the
#     top-tertile deficiency level is genuinely transportable, and (b) the
#     development->deficiency gradient has a STABLE SIGN across countries. Both
#     are testable here and the script prints the diagnostics. The known
#     counter-example is Malawi (the non-West-African country that anti-
#     correlates on transport): if its internal gradient sign disagrees with the
#     training countries', the anchor will MISLEAD, and the script flags it.
#
# Run from the repo root:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "sensitivity/13_gradient_anchored_calibration.R"
# =============================================================================

# --- locate and reuse the simplified-subset helpers --------------------------
find_helpers <- function() {
  cands <- c(
    file.path("simplified subset", "methods", "_helpers.R"),
    file.path("..", "simplified subset", "methods", "_helpers.R"),
    file.path("methods", "_helpers.R")
  )
  p <- cands[file.exists(cands)][1]
  if (is.na(p)) stop("Could not find simplified subset/methods/_helpers.R from ", getwd())
  p
}
source(find_helpers())

set.seed(1)

OUTCOMES <- c("women_iron", "women_vitA", "child_iron")   # illustrative panel
DEV_VARS <- c("human_modification_index",                 # primary + robustness
              "built_surface_nonres_10km",
              "pop_density_cropland_10km")
PRIMARY_DEV <- DEV_VARS[1]

# --- calibration primitives --------------------------------------------------

# One global logit shift chosen so the weighted mean of the shifted predictions
# over `idx` equals `target`. Applied to ALL predictions -> ranking unchanged.
logit_shift_to <- function(phat, w, idx, target) {
  phat <- clip01(phat)
  f <- function(delta) {
    p <- plogis(qlogis(phat) + delta)
    sum(w[idx] * p[idx]) / sum(w[idx]) - target
  }
  lo <- f(-12); hi <- f(12)
  if (!is.finite(lo) || !is.finite(hi) || lo * hi > 0) return(phat)  # unreachable
  plogis(qlogis(phat) + uniroot(f, c(-12, 12))$root)
}

# Top-tertile membership (most-developed third) of a development variable.
top_tertile <- function(dev) {
  thr <- suppressWarnings(stats::quantile(dev, 2/3, na.rm = TRUE, names = FALSE))
  is.finite(dev) & dev >= thr
}
bot_tertile <- function(dev) {
  thr <- suppressWarnings(stats::quantile(dev, 1/3, na.rm = TRUE, names = FALSE))
  is.finite(dev) & dev <= thr
}

# Pooled pop-weighted TRUTH level in the top development tertile of the TRAINING
# countries (tertiles computed within each training country, then pooled).
borrow_top_tertile_level <- function(dev_tr, y_tr, w_tr, country_tr) {
  num <- den <- 0
  for (c in unique(country_tr)) {
    ic  <- country_tr == c
    top <- top_tertile(dev_tr[ic])
    if (sum(top) == 0) next
    num <- num + sum(w_tr[ic][top] * y_tr[ic][top])
    den <- den + sum(w_tr[ic][top])
  }
  if (den == 0) return(NA_real_)
  num / den
}

# Sign of the development->deficiency gradient (bottom-tertile minus top-tertile
# mean prevalence), pooled across whatever countries are passed. Positive means
# "deficiency higher in LESS-developed areas" (the expected direction).
gradient_contrast <- function(dev, y, w, country) {
  topn <- topd <- botn <- botd <- 0
  for (c in unique(country)) {
    ic <- country == c
    tp <- top_tertile(dev[ic]); bt <- bot_tertile(dev[ic])
    if (sum(tp)) { topn <- topn + sum(w[ic][tp]*y[ic][tp]); topd <- topd + sum(w[ic][tp]) }
    if (sum(bt)) { botn <- botn + sum(w[ic][bt]*y[ic][bt]); botd <- botd + sum(w[ic][bt]) }
  }
  if (topd == 0 || botd == 0) return(NA_real_)
  (botn / botd) - (topn / topd)
}

# --- main LOCO loop ----------------------------------------------------------

run_outcome <- function(df, outcome, dev_var) {
  d <- prep_outcome(df, outcome)
  if (!(dev_var %in% d$predictors))
    stop("Development variable ", dev_var, " not among predictors.")
  dev_all <- d$X[, dev_var]
  countries <- unique(d$country)
  rows <- list()

  for (held in countries) {
    tr <- d$country != held; te <- !tr
    if (sum(te) < 6 || sum(tr) < 10) next

    predict_fn <- fit_area(d$X[tr, , drop = FALSE], d$y[tr], d$n[tr], d$ndef[tr],
                           family = "gaussian")
    phat  <- clip01(predict_fn(d$X[te, , drop = FALSE]))
    truth <- d$y[te]
    w     <- d$n[te]                 # population weight proxy (denominator)
    dev_t <- dev_all[te]

    # Rung 1 (oracle): shift to the target's OWN true national mean.
    natl_target <- sum(w * truth) / sum(w)
    p_natl <- logit_shift_to(phat, w, rep(TRUE, length(phat)), natl_target)

    # Rung 2 (deployable): borrow the top-tertile level from TRAINING countries,
    # shift so the target's top-tertile mean prediction matches it.
    borrowed <- borrow_top_tertile_level(dev_all[tr], d$y[tr], d$n[tr], d$country[tr])
    top_t    <- top_tertile(dev_t)
    p_grad   <- if (is.finite(borrowed) && sum(top_t) >= 2)
                  logit_shift_to(phat, w, top_t, borrowed) else phat

    # Diagnostics for the identifying assumptions ---------------------------
    grad_train  <- gradient_contrast(dev_all[tr], d$y[tr], d$n[tr], d$country[tr])
    grad_target <- gradient_contrast(dev_t, truth, w, rep(held, length(truth)))
    sign_ok <- is.finite(grad_train) && is.finite(grad_target) &&
               sign(grad_train) == sign(grad_target)
    # Is the borrowed top-tertile level close to the target's actual one?
    true_top <- if (sum(top_t)) sum(w[top_t]*truth[top_t]) / sum(w[top_t]) else NA_real_

    rows[[held]] <- data.frame(
      outcome = outcome, dev = dev_var, held_out = held, n_areas = sum(te),
      MAE_raw   = mae_pp(phat,   truth),
      MAE_grad  = mae_pp(p_grad, truth),
      MAE_natl  = mae_pp(p_natl, truth),
      bias_raw  = bias_pp(phat,   truth),
      bias_grad = bias_pp(p_grad, truth),
      bias_natl = bias_pp(p_natl, truth),
      spearman  = spear(phat, truth),                 # identical across all three
      borrow_top = round(borrowed, 3),
      true_top   = round(true_top, 3),
      grad_train = round(grad_train, 3),
      grad_target= round(grad_target, 3),
      sign_ok    = sign_ok,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# --- 1) primary panel: three outcomes, primary development gradient ----------
df <- load_level("admin2")

cat("=====================================================================\n")
cat(" RUNG 2  Gradient-anchored calibration (deployable, no national total)\n")
cat(" Anchor = top-development-tertile level borrowed from TRAINING countries\n")
cat(" Development axis:", PRIMARY_DEV, "\n")
cat("=====================================================================\n")

primary <- do.call(rbind, lapply(OUTCOMES, function(o) run_outcome(df, o, PRIMARY_DEV)))
show <- primary[, c("outcome","held_out","n_areas",
                    "MAE_raw","MAE_grad","MAE_natl",
                    "bias_raw","bias_grad","bias_natl","sign_ok")]
print(show, row.names = FALSE)

cat("\n--- Mean over held-out countries, by outcome (MAE in pp) ---\n")
agg <- aggregate(cbind(MAE_raw, MAE_grad, MAE_natl) ~ outcome, primary, mean)
agg[, -1] <- round(agg[, -1], 2)
print(agg, row.names = FALSE)

cat("\nReading it: MAE_natl is the ORACLE (uses the target's true mean).\n")
cat("MAE_grad is DEPLOYABLE (borrows only from surveyed countries).\n")
cat("Recovery = how far grad closes the raw->oracle gap:\n")
rec <- with(agg, round(100 * (MAE_raw - MAE_grad) /
                       pmax(MAE_raw - MAE_natl, 1e-6)))
print(data.frame(outcome = agg$outcome, recovery_pct = rec), row.names = FALSE)

# --- 2) the honesty diagnostic: when does the anchor MISLEAD? ----------------
cat("\n--- Gradient-sign check (identifying assumption) ---\n")
cat("sign_ok = FALSE means the target country's development->deficiency gradient\n")
cat("runs OPPOSITE to the training countries', so the borrowed anchor is invalid\n")
cat("and the gradient calibration should NOT be trusted for that country.\n\n")
diag <- primary[, c("outcome","held_out","grad_train","grad_target","sign_ok",
                    "borrow_top","true_top","MAE_grad","MAE_natl")]
print(diag, row.names = FALSE)

viol <- primary[!primary$sign_ok, c("outcome","held_out")]
if (nrow(viol)) {
  cat("\nFlagged (gradient assumption violated) -- gradient anchor unreliable here:\n")
  print(unique(viol), row.names = FALSE)
} else {
  cat("\nNo gradient-sign violations in this panel.\n")
}

# --- 3) robustness across alternative development axes (women_iron) ----------
cat("\n--- Robustness of the gradient anchor to the development axis (women_iron) ---\n")
rob <- do.call(rbind, lapply(DEV_VARS, function(v) {
  r <- run_outcome(df, "women_iron", v)
  if (is.null(r)) return(NULL)
  data.frame(dev = v,
             MAE_raw  = round(mean(r$MAE_raw),  2),
             MAE_grad = round(mean(r$MAE_grad), 2),
             MAE_natl = round(mean(r$MAE_natl), 2),
             n_sign_ok = sum(r$sign_ok), n_countries = nrow(r))
}))
print(rob, row.names = FALSE)

# --- 4) persist ---------------------------------------------------------------
out_dir <- if (dir.exists("results/sensitivity")) "results/sensitivity" else "."
out_csv <- file.path(out_dir, "gradient_anchored_calibration.csv")
utils::write.csv(primary, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))

cat("\nTakeaways to check against the manuscript:\n")
cat(" * Does the deployable gradient anchor recover a useful share of the oracle\n")
cat("   national-anchor's level correction, WITHOUT any target-country data?\n")
cat(" * Where sign_ok is FALSE (expected: Malawi on some outcomes), the anchor\n")
cat("   is invalid -- this is the honest limit of rung 2 and belongs in the text.\n")
