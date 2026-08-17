# Task 2: anchor LOCO-transported predictions to a known national prevalence.
#
# Problem (from the plan): a model trained on 3 countries and applied to a
# 4th tends to get district RANKING right but the overall LEVEL drifts toward
# the training countries' average. Many un-surveyed countries still have a
# single external national prevalence number (VMNIS, BRINDA, a past DHS
# anaemia survey). This task uses that one number to correct the level
# without disturbing the ranking.
#
# Self-contained (not sourcing task1_explore.R) -- reuses the same data-
# loading/imputation pattern but keeps the model to plain lm, since Task 2 is
# about the anchoring step, not model choice.

library(dplyr)

admin2 <- read.csv("simplified subset/data/mn_admin2.csv")
dict   <- read.csv("simplified subset/data_dictionary.csv")
predictors <- dict$variable[dict$role == "predictor"]

impute_predictors <- function(data, cols, method = c("mean", "median")) {
  method <- match.arg(method)
  fn <- switch(method, mean = mean, median = median)
  for (col in cols) {
    data[[col]][is.na(data[[col]])] <- fn(data[[col]], na.rm = TRUE)
  }
  data
}

model_data <- admin2 %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors))
data_imputed <- impute_predictors(model_data, predictors, method = "mean")

baseline_formula <- reformulate(predictors, response = "prev_women_iron")

# ── Step 1: LOCO predictions + each held-out country's "known" national
# prevalence (standing in for an external number like VMNIS/BRINDA/DHS) ─────
run_loco <- function(data) {
  bind_rows(lapply(unique(data$country), function(held_out) {
    train <- data[data$country != held_out, ]
    test  <- data[data$country == held_out, ]
    fit <- lm(baseline_formula, data = train, weights = n_women_iron)
    test$pred_loco <- predict(fit, newdata = test)
    test
  }))
}

loco_results <- run_loco(data_imputed)

# Sample-size-weighted average of the TRUE prevalence within the held-out
# country -- this is what we're pretending is externally known, and it's
# the value every held-out country's predictions get anchored to.
national_prevalence <- loco_results %>%
  group_by(country) %>%
  summarise(
    national_prev_true = weighted.mean(prev_women_iron, n_women_iron),
    national_prev_pred = weighted.mean(pred_loco,       n_women_iron),
    .groups = "drop"
  )

cat("True vs. LOCO-predicted national prevalence, by held-out country:\n")
print(national_prevalence)
cat("\n(The gap between these two columns is exactly the level drift Task 2 is meant to fix.)\n")

# =============================================================================
# Step 2: extrapolation diagnostic -- random Tukey (halfspace) depth
#
# Sierra Leone's LOCO predictions came back LOWER than its true prevalence,
# not pulled toward the training mean as naive shrinkage logic would predict
# -- because several of its predictors sit outside the training countries'
# joint covariate range, so the linear model is extrapolating, not
# interpolating (4 of 14 SL districts got impossible negative predictions).
#
# Mahalanobis distance would catch this but assumes the training covariates
# form one elliptical cloud -- not a safe assumption with 4 separate
# countries. Tukey (halfspace) depth makes no shape assumption: a point's
# depth is the SMALLEST fraction of training points on one side of any
# hyperplane through it, so a point outside the training data's convex hull
# gets depth 0 regardless of the cloud's shape.
#
# TRIED FIRST, DOESN'T WORK HERE: mrfDepth::hdepth() (halfspace depth
# directly). With p=16 predictors and n~180 training districts, halfspace
# depth saturates near 0 for almost every point, not just extrapolated
# ones -- verified with a sanity check on a well-behaved 180x16 Gaussian
# cloud evaluated against itself, where even the most central point only
# reached depth ~0.03 (typical point ~0.006, i.e. ~1 of 180). This is a
# known curse-of-dimensionality property of rank-based depth at this n/p
# ratio, not specific to our data -- confirmed when EVERY district in ALL
# FOUR countries came back with LOCO halfspace depth exactly 0, which
# carries no information (see git history for that version).
#
# USED INSTEAD: mrfDepth::outlyingness() -- Stahel-Donoho outlyingness. Same
# worst-case-direction search as hdepth() (affine-invariant, directions
# orthogonal to hyperplanes spanned by training points), but instead of a
# rank-based depth that floors at 0, it reports a robustly standardized
# DISTANCE along the worst direction -- continuous and still informative
# past the point where rank-based depth collapses.
#
# Computed purely from covariates (predictors), never the outcome -- exactly
# what's needed for a genuinely unsurveyed country where prev_women_iron
# isn't available to check against.
# =============================================================================

library(mrfDepth)

# ── Step 2a: outlyingness for every district, computed against its own LOCO
# training set (the other 3 countries) -- same train/test split as Step 1.
# Higher = further outside the training cloud in its worst direction. ──────
depth_results <- bind_rows(lapply(unique(data_imputed$country), function(held_out) {
  train <- data_imputed[data_imputed$country != held_out, ]
  test  <- data_imputed[data_imputed$country == held_out, ]
  out <- outlyingness(x = as.matrix(train[predictors]), z = as.matrix(test[predictors]))
  test$outlyingness <- out$outlyingnessZ
  test %>% select(country, admin1, admin2, outlyingness)
}))

loco_results <- loco_results %>%
  left_join(depth_results, by = c("country", "admin1", "admin2"))

# ── Step 2b: does high outlyingness line up with the extrapolated (negative)
# predictions we found in Sierra Leone? ─────────────────────────────────────
cat("\nSierra Leone: outlyingness vs. LOCO prediction, sorted high to low (most extrapolated first):\n")
loco_results %>%
  filter(country == "Sierra Leone") %>%
  select(admin2, outlyingness, prev_women_iron, pred_loco) %>%
  arrange(desc(outlyingness)) %>%
  print(row.names = FALSE)

cat("\nOutlyingness summary by country (higher = more extrapolated relative to the other 3):\n")
loco_results %>%
  group_by(country) %>%
  summarise(min = min(outlyingness), median = median(outlyingness), max = max(outlyingness)) %>%
  print()

# ── Step 2c: Mahalanobis distance, for direct comparison -- the classical
# elliptical-cloud distance outlyingness was chosen instead of.
#
# The 16 predictors are severely multicollinear (multiple NDVI/rainfall
# variables at different radii/months capturing overlapping signal) -- the
# training covariance's condition number runs 2x10^5 to 1.8x10^6 across the
# four LOCO training sets (checked directly; well-conditioned would be
# closer to 1-100). Naive Mahalanobis distance inverts that covariance, so
# it gets dominated by whichever near-zero-variance direction the inversion
# amplifies most -- numerically real, but not obviously meaningful. Computed
# two ways so the instability is visible rather than hidden: naive (raw
# covariance) and ridge-regularized (shrunk toward the diagonal, damping the
# near-zero eigenvalues before inverting).
# =============================================================================

ridge_shrink <- 0.1  # 0 = no shrinkage (naive); 1 = fully diagonal (ignores all correlation)

mahalanobis_results <- bind_rows(lapply(unique(data_imputed$country), function(held_out) {
  train <- data_imputed[data_imputed$country != held_out, ]
  test  <- data_imputed[data_imputed$country == held_out, ]
  X_train <- as.matrix(train[predictors])
  center  <- colMeans(X_train)
  cov_raw <- cov(X_train)
  cov_shrunk <- (1 - ridge_shrink) * cov_raw + ridge_shrink * diag(diag(cov_raw))

  test$mahal_naive   <- mahalanobis(as.matrix(test[predictors]), center, cov_raw)
  test$mahal_shrunk  <- mahalanobis(as.matrix(test[predictors]), center, cov_shrunk)
  test %>% select(country, admin1, admin2, mahal_naive, mahal_shrunk)
}))

loco_results <- loco_results %>%
  left_join(mahalanobis_results, by = c("country", "admin1", "admin2"))

cat("\nSierra Leone: all three diagnostics vs. LOCO prediction, sorted by shrunk Mahalanobis:\n")
loco_results %>%
  filter(country == "Sierra Leone") %>%
  select(admin2, outlyingness, mahal_naive, mahal_shrunk, prev_women_iron, pred_loco) %>%
  arrange(desc(mahal_shrunk)) %>%
  print(row.names = FALSE)

cat("\nMahalanobis summary by country (naive vs. shrunk; higher = more extrapolated):\n")
loco_results %>%
  group_by(country) %>%
  summarise(
    naive_median  = median(mahal_naive),  naive_max  = max(mahal_naive),
    shrunk_median = median(mahal_shrunk), shrunk_max = max(mahal_shrunk)
  ) %>%
  print()

cat("\nRank correlation (Spearman) between outlyingness and each Mahalanobis version, pooled across countries:\n")
cat("  outlyingness vs. naive Mahalanobis: ",
    cor(loco_results$outlyingness, loco_results$mahal_naive,  method = "spearman"), "\n")
cat("  outlyingness vs. shrunk Mahalanobis:",
    cor(loco_results$outlyingness, loco_results$mahal_shrunk, method = "spearman"), "\n")

# =============================================================================
# Step 3: anchor correction
#
# Shift each held-out country's LOCO predictions on the log-odds scale so
# their n_women_iron-weighted average matches Step 1's "known" national
# prevalence -- a single number per country, applied uniformly, so within-
# country ORDERING is untouched (a log-odds shift is monotonic) while the
# LEVEL is corrected. Reports both: ranking metrics (should not move) and
# level metrics -- MAE and bias -- before vs. after (should improve, since
# Step 1 showed the drift is large and systematic, not just noise).
# =============================================================================

# Raw lm predictions aren't bounded to (0,1) -- 14 of 206 are <= 0 (Step 2's
# extrapolation story, worst for Sierra Leone). Clip before the logit
# transform; this only affects how the most extreme extrapolated districts
# enter the anchor calculation, not what the anchor target is.
clip_prob <- function(p, eps = 1e-3) pmin(pmax(p, eps), 1 - eps)

# The single log-odds shift `delta` for one country such that the weighted
# average of plogis(qlogis(clipped pred) + delta) equals the target national
# prevalence. Solved numerically (weighted mean of a back-transformed shift
# has no closed form).
solve_anchor_shift <- function(pred, weights, target) {
  logit <- qlogis(clip_prob(pred))
  objective <- function(delta) weighted.mean(plogis(logit + delta), weights) - target
  uniroot(objective, interval = c(-30, 30))$root
}

loco_results <- loco_results %>%
  left_join(national_prevalence %>% select(country, national_prev_true), by = "country")

anchor_shifts <- loco_results %>%
  group_by(country) %>%
  summarise(shift = solve_anchor_shift(pred_loco, n_women_iron, unique(national_prev_true)), .groups = "drop")

cat("\nLog-odds anchor shift, per held-out country:\n")
print(anchor_shifts, row.names = FALSE)

loco_results <- loco_results %>%
  left_join(anchor_shifts, by = "country") %>%
  mutate(pred_anchored = plogis(qlogis(clip_prob(pred_loco)) + shift))

# ── Step 3a: sanity check -- did the anchor actually hit its target? ────────
cat("\nSanity check -- weighted mean of anchored predictions vs. true national prevalence:\n")
loco_results %>%
  group_by(country) %>%
  summarise(anchored_mean = weighted.mean(pred_anchored, n_women_iron), true_mean = unique(national_prev_true)) %>%
  print()

# ── Step 3b: ranking metrics -- should not change (plan's explicit ask) ─────
cat("\nRanking check -- Spearman correlation between raw and anchored predictions, within country (expect exactly 1.0):\n")
loco_results %>%
  group_by(country) %>%
  summarise(rank_corr = cor(pred_loco, pred_anchored, method = "spearman")) %>%
  print()

cat("\nDecision-relevant check -- Task 1's top-20%-within-country flag, identical before/after anchoring?\n")
loco_results %>%
  group_by(country) %>%
  mutate(
    high_burden_before = pred_loco     >= quantile(pred_loco,     0.80),
    high_burden_after  = pred_anchored >= quantile(pred_anchored, 0.80)
  ) %>%
  ungroup() %>%
  summarise(pct_identical = mean(high_burden_before == high_burden_after)) %>%
  print()

# ── Step 3c: level metrics -- MAE and bias, before vs. after ────────────────
cat("\nLevel metrics by country -- MAE and bias, before vs. after anchoring:\n")
loco_results %>%
  group_by(country) %>%
  summarise(
    mae_before  = mean(abs(pred_loco     - prev_women_iron)),
    mae_after   = mean(abs(pred_anchored - prev_women_iron)),
    bias_before = mean(pred_loco     - prev_women_iron),
    bias_after  = mean(pred_anchored - prev_women_iron)
  ) %>%
  print()

cat("\nLevel metrics pooled (all 4 held-out countries):\n")
loco_results %>%
  summarise(
    mae_before  = mean(abs(pred_loco     - prev_women_iron)),
    mae_after   = mean(abs(pred_anchored - prev_women_iron)),
    bias_before = mean(pred_loco     - prev_women_iron),
    bias_after  = mean(pred_anchored - prev_women_iron)
  ) %>%
  print()

# =============================================================================
# Step 4: isotonic-regression recalibration -- an alternative to Step 3's
# rigid log-odds shift. The log-odds shift assumes ONE fixed offset applies
# uniformly across the whole prediction range; isotonic regression instead
# fits a flexible MONOTONIC step function g() that can correct low and high
# predictions differently -- at the cost of needing real calibration data to
# fit g() against, rather than just the one known national number.
#
# Two variants:
#  (a) PURE recalibration -- fit g() only from the training countries' own
#      out-of-fold CV predictions, apply it to the held-out country, and see
#      how close its own weighted average lands to the true national value
#      ON ITS OWN, no forced correction. Tests whether a flexible-shape but
#      externally-blind recalibration can match a rigid shift that KNOWS the
#      answer.
#  (b) ANCHORED -- take (a)'s isotonic output and layer Step 3's log-odds
#      shift on top (reusing solve_anchor_shift unchanged -- composing two
#      monotonic maps is still monotonic, so ranking stays intact), so it
#      also hits the known national value exactly. This is the practical
#      version for when the national number actually is available, which is
#      what the plan's Task 2 setup assumes.
#
# Calibration data must be OUT-OF-FOLD (not in-sample -- fitting g() against
# the same lm's in-sample predictions would be circular) and can't include
# Sierra Leone as a contributor -- same 14-districts/17-parameters rank-
# deficiency Task 1 hit, so SL is excluded from calibration pools exactly as
# it was excluded from Task 1's within-country CV.
# =============================================================================

# ── Step 4a: within-country CV, self-contained here (not sourcing
# task1_explore.R) -- same approach Task 1 used. ─────────────────────────────
assign_folds <- function(n, k, seed = 1) {
  set.seed(seed)
  sample(rep(1:k, length.out = n))
}

run_country_cv <- function(data, k = 5) {
  data$fold <- assign_folds(nrow(data), k)
  data$pred_cv <- NA_real_
  for (f in seq_len(k)) {
    train <- data[data$fold != f, ]
    test  <- data[data$fold == f, ]
    fit <- lm(baseline_formula, data = train, weights = n_women_iron)
    data$pred_cv[data$fold == f] <- predict(fit, newdata = test)
  }
  data
}

# Sierra Leone (14 districts, 17 lm parameters) is rank-deficient even on
# its full sample, let alone a CV fold -- excluded from calibration pools.
cv_viable_countries <- c("Gambia", "Ghana", "Malawi")

# ── Step 4b: fit an isotonic calibrator per LOCO round, from that round's
# training countries' own out-of-fold predictions (excluding SL whenever
# it's one of them -- so the SL-held-out round is the only one with all 3
# training countries contributing). ──────────────────────────────────────
fit_isotonic_calibrator <- function(held_out) {
  calibration_countries <- intersect(setdiff(unique(data_imputed$country), held_out), cv_viable_countries)

  calib_data <- data_imputed %>%
    filter(country %in% calibration_countries) %>%
    group_by(country) %>%
    group_modify(~ run_country_cv(.x)) %>%
    ungroup()

  iso <- isoreg(x = calib_data$pred_cv, y = calib_data$prev_women_iron)
  as.stepfun(iso)
}

isotonic_results <- bind_rows(lapply(unique(data_imputed$country), function(held_out) {
  g <- fit_isotonic_calibrator(held_out)
  test <- loco_results[loco_results$country == held_out, ]
  test$pred_isotonic <- g(test$pred_loco)
  test %>% select(country, admin1, admin2, pred_isotonic)
}))

loco_results <- loco_results %>%
  left_join(isotonic_results, by = c("country", "admin1", "admin2"))

# ── Step 4c: Option 1 -- how close does the PURE isotonic recalibration land
# to the true national value, with no forced correction? ────────────────────
cat("\nOption 1 (pure isotonic, no anchor) -- weighted mean vs. true national prevalence:\n")
loco_results %>%
  group_by(country) %>%
  summarise(isotonic_mean = weighted.mean(pred_isotonic, n_women_iron), true_mean = unique(national_prev_true)) %>%
  print()

cat("\nOption 1 level metrics -- MAE/bias: raw vs. Step 3's log-odds anchor vs. pure isotonic:\n")
loco_results %>%
  group_by(country) %>%
  summarise(
    mae_raw    = mean(abs(pred_loco     - prev_women_iron)),
    mae_anchor = mean(abs(pred_anchored - prev_women_iron)),
    mae_iso    = mean(abs(pred_isotonic - prev_women_iron)),
    bias_raw    = mean(pred_loco     - prev_women_iron),
    bias_anchor = mean(pred_anchored - prev_women_iron),
    bias_iso    = mean(pred_isotonic - prev_women_iron)
  ) %>%
  print()

# ── Step 4d: Option 2 -- layer Step 3's log-odds shift on top of the
# isotonic output, so it also hits the known national anchor exactly ────────
isotonic_anchor_shifts <- loco_results %>%
  group_by(country) %>%
  summarise(shift_iso = solve_anchor_shift(pred_isotonic, n_women_iron, unique(national_prev_true)), .groups = "drop")

loco_results <- loco_results %>%
  left_join(isotonic_anchor_shifts, by = "country") %>%
  mutate(pred_isotonic_anchored = plogis(qlogis(clip_prob(pred_isotonic)) + shift_iso))

cat("\nOption 2 sanity check -- weighted mean after anchoring the isotonic output:\n")
loco_results %>%
  group_by(country) %>%
  summarise(mean = weighted.mean(pred_isotonic_anchored, n_women_iron), true_mean = unique(national_prev_true)) %>%
  print()

cat("\nOption 2 level metrics -- MAE/bias: Step 3's log-odds anchor alone vs. isotonic+anchor:\n")
loco_results %>%
  group_by(country) %>%
  summarise(
    mae_anchor_only  = mean(abs(pred_anchored          - prev_women_iron)),
    mae_iso_anchor   = mean(abs(pred_isotonic_anchored - prev_women_iron)),
    bias_anchor_only = mean(pred_anchored          - prev_women_iron),
    bias_iso_anchor  = mean(pred_isotonic_anchored - prev_women_iron)
  ) %>%
  print()

cat("\nDecision-relevant check -- Task 1's top-20%-within-country flag, raw vs. isotonic+anchor:\n")
loco_results %>%
  group_by(country) %>%
  mutate(
    high_burden_before = pred_loco              >= quantile(pred_loco,              0.80),
    high_burden_after   = pred_isotonic_anchored >= quantile(pred_isotonic_anchored, 0.80)
  ) %>%
  ungroup() %>%
  group_by(country) %>%
  summarise(pct_identical = mean(high_burden_before == high_burden_after)) %>%
  print()

cat("\nRanking check -- Spearman correlation, raw vs. isotonic+anchor predictions, within country:\n")
loco_results %>%
  group_by(country) %>%
  summarise(rank_corr = cor(pred_loco, pred_isotonic_anchored, method = "spearman")) %>%
  print()

# =============================================================================
# Step 5: summary table -- Task 2 at a glance
#
# One table: MAE and bias by held-out country, across all four prediction
# variants built above (raw LOCO, log-odds anchor, pure isotonic, isotonic+
# anchor). A second, small table: whether Task 1's top-20%-within-country
# flag survived each of the two ANCHORED corrections, relative to raw LOCO
# (pure isotonic alone was never claimed to preserve ranking, so it's not
# included in that check -- only the weighted-mean-matches-target property
# was ever guaranteed for it).
# =============================================================================

task2_summary <- bind_rows(
  loco_results %>% mutate(pred = pred_loco,             method = "Raw LOCO"),
  loco_results %>% mutate(pred = pred_anchored,          method = "Log-odds anchor"),
  loco_results %>% mutate(pred = pred_isotonic,          method = "Pure isotonic (no anchor)"),
  loco_results %>% mutate(pred = pred_isotonic_anchored, method = "Isotonic + anchor")
) %>%
  group_by(country, method) %>%
  summarise(mae = mean(abs(pred - prev_women_iron)), bias = mean(pred - prev_women_iron), .groups = "drop")

cat("\n===== Task 2 summary: MAE and bias by method =====\n")
print(task2_summary, n = Inf)

ranking_summary <- loco_results %>%
  group_by(country) %>%
  mutate(
    flag_raw        = pred_loco               >= quantile(pred_loco, 0.80),
    flag_anchor     = pred_anchored            >= quantile(pred_anchored, 0.80),
    flag_iso_anchor = pred_isotonic_anchored   >= quantile(pred_isotonic_anchored, 0.80)
  ) %>%
  summarise(
    top20_match_anchor     = mean(flag_raw == flag_anchor),
    top20_match_iso_anchor = mean(flag_raw == flag_iso_anchor),
    .groups = "drop"
  )

cat("\n===== Task 2 summary: top-20%-within-country flag preserved vs. raw LOCO =====\n")
print(ranking_summary, n = Inf)

# =============================================================================
# Step 6: conclusions and takeaways
#
# Task 1 was about RANKING -- which districts are highest-burden. Task 2
# picks up where it left off: LOCO-transported predictions get the ranking
# roughly right but drift in LEVEL, so the plan's ask was a one-number
# correction to fix that drift using a known national anchor.
#
# Investigating WHY the drift happens (Sierra Leone's impossible negative
# predictions) raised a question the plan didn't ask but that any anchoring
# scheme should care about: can we tell, in advance, how much to trust a
# given country's transported predictions? That detour produced a real
# split -- a MODEL-AGNOSTIC distance (outlyingness) that would survive
# swapping lm for a different learner, and a MODEL-SPECIFIC one (Mahalanobis,
# ~= leverage) that explains lm's own failures well but has no reason to
# transfer elsewhere -- and they disagree with each other (Spearman ~0.6)
# because they're answering different questions, not because one is wrong.
# Plain halfspace/Tukey depth was a dead end at this n/p ratio (16 collinear
# covariates, ~180 training districts) -- worth retrying only once n grows
# substantially relative to p; outlyingness isn't a fallback for that, it's
# the right tool at this scale.
#
# Back to the actual correction: anchoring reliably fixes LEVEL (bias ~0,
# MAE -38% pooled, every country) but is MATHEMATICALLY INCAPABLE of
# improving RANKING -- a monotonic log-odds shift can't reorder anything, so
# Step 3's 100%-identical top-20% flag was guaranteed, not measured. Good
# targeting is entirely inherited from whichever model produced the raw
# predictions; that's still Task 1's problem, not something a correction
# layered on top can fix. Isotonic recalibration relaxes the "monotonic"
# constraint to "flexible-shaped but still monotonic," and that flexibility
# is real -- best of every method tried when the calibration pool covered
# the target's true range (Sierra Leone), worse than doing nothing when it
# didn't (Gambia). The one time it touched ranking, sparse-data ties made
# it worse, never better.
#
# Two natural revisits when more data or tools are available: (1) more
# countries/districts should widen isotonic calibration coverage and make
# its Sierra-Leone-vs-Gambia split less of a coin flip; (2) swapping in a
# different base learner (Task 1's blend, or ranger) changes which
# diagnostic is the model-specific one to pair with outlyingness -- ranger's
# analogue is proximity/leaf co-occurrence, not Mahalanobis distance.
# =============================================================================
