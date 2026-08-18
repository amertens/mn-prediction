# Task 3: a proper model for proportions, with sampling weights.
#
# Task 1/2 used plain lm() throughout -- weighted, but with no constraint
# that predictions stay in [0,1] (Task 2's whole extrapolation storyline
# started from exactly that gap: lm predicted negative prevalence for
# several Sierra Leone districts). This task refits using a binomial/
# quasi-binomial GLM (logit link), which bounds predictions to (0,1) by
# construction, and compares calibration/stability against the plain-lm
# baseline.
#
# Self-contained (not sourcing task1/task2_explore.R), same data-loading
# pattern as both.

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
  select(country, admin1, admin2, prev_women_iron, n_women_iron, var_women_iron, all_of(predictors))
data_imputed <- impute_predictors(model_data, predictors, method = "mean")

baseline_formula <- reformulate(predictors, response = "prev_women_iron")

# ── Step 1: in-sample comparison -- does the quasi-binomial GLM actually
# stay bounded to (0,1) where lm didn't? ─────────────────────────────────────
#
# n_women_iron weighting only for now (matches Task 1/2's weighting choice).
# var_women_iron -- the plan's other suggested weight, as 1/var -- is exactly
# 0 for 71 of 206 districts (all zero-deficient-case districts, small n:
# median n=7), so 1/var is undefined for a third of the data. Needs a fix
# (floor, continuity correction) before it's usable -- not attempted yet,
# flagged for a later step.
fit_lm <- lm(baseline_formula, data = data_imputed, weights = n_women_iron)
fit_qb <- glm(baseline_formula, data = data_imputed, weights = n_women_iron, family = quasibinomial(link = "logit"))

data_imputed$pred_lm <- predict(fit_lm)
data_imputed$pred_qb <- predict(fit_qb, type = "response")

cat("In-sample predicted range, plain lm:\n")
cat("  range:", range(data_imputed$pred_lm), " -- n outside [0,1]:", sum(data_imputed$pred_lm < 0 | data_imputed$pred_lm > 1), "\n")

cat("\nIn-sample predicted range, quasi-binomial GLM (logit link):\n")
cat("  range:", range(data_imputed$pred_qb), " -- n outside [0,1]:", sum(data_imputed$pred_qb < 0 | data_imputed$pred_qb > 1), "\n")

cat("\nIn-sample fit comparison -- MAE (both weighted by n_women_iron):\n")
weighted_mae <- function(pred, actual, w) sum(w * abs(pred - actual)) / sum(w)
cat("  lm: ", weighted_mae(data_imputed$pred_lm, data_imputed$prev_women_iron, data_imputed$n_women_iron), "\n")
cat("  qb: ", weighted_mae(data_imputed$pred_qb, data_imputed$prev_women_iron, data_imputed$n_women_iron), "\n")

# ── Step 1a: precision weighting (1/var_women_iron) -- guarded, not fixed ───
# 71 of 206 districts have var_women_iron == 0 (zero deficient cases
# observed, all small-n), so 1/var is undefined for a third of the data.
# Rather than silently producing Inf weights, this stops loudly. A real fix
# (variance floor, continuity correction e.g. Wilson/Agresti-Coull) is left
# for later -- this dataset happens to hit the degenerate case for this
# outcome, but other outcomes/countries may not, so the function itself
# stays general rather than baking in a workaround for THIS data's problem.
precision_weights <- function(var) {
  if (any(var == 0)) {
    stop("precision_weights(): ", sum(var == 0), " of ", length(var),
         " observations have var == 0 (1/var undefined) -- not handled yet.")
  }
  1 / var
}

cat("\nConfirming the guard fires on this outcome (expected -- do not proceed past this for prev_women_iron):\n")
result <- tryCatch(precision_weights(data_imputed$var_women_iron), error = function(e) conditionMessage(e))
cat(" ", result, "\n")

# =============================================================================
# Step 2: within-country CV and LOCO for both models -- Task 1's Steps 4-5
# machinery, applied to lm vs. quasi-binomial GLM. Deliberately NOT Task 1's
# multi-learner blend or nested-CV honesty layer (that answers a different
# question -- how to combine several learners -- not "does this one
# functional-form choice change decisions," which is what Task 3 is asking).
# =============================================================================

# ── Step 2a: both models as fit/predict pairs, same shape as Task 1's Step 7a
# learner list, so CV/LOCO code is written once and shared. ─────────────────
learner_lm <- list(
  fit     = function(train) lm(baseline_formula, data = train, weights = n_women_iron),
  predict = function(model, test) predict(model, newdata = test)
)

learner_qb <- list(
  fit     = function(train) glm(baseline_formula, data = train, weights = n_women_iron, family = quasibinomial(link = "logit")),
  predict = function(model, test) predict(model, newdata = test, type = "response")
)

learners <- list(lm = learner_lm, qb = learner_qb)

# ── Step 2b: within-country CV -- Gambia/Ghana/Malawi only. Sierra Leone (14
# districts, 17 lm parameters) is rank-deficient, same exclusion as Task 1. ──
assign_folds <- function(n, k, seed = 1) {
  set.seed(seed)
  sample(rep(1:k, length.out = n))
}

run_country_cv_generic <- function(data, learner, k = 5) {
  data$fold <- assign_folds(nrow(data), k)
  data$pred <- NA_real_
  for (f in seq_len(k)) {
    train <- data[data$fold != f, ]
    test  <- data[data$fold == f, ]
    model <- learner$fit(train)
    data$pred[data$fold == f] <- learner$predict(model, test)
  }
  data
}

cv_countries <- c("Gambia", "Ghana", "Malawi")

cv_results <- data_imputed %>%
  filter(country %in% cv_countries) %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron)

for (nm in names(learners)) {
  res <- data_imputed %>%
    filter(country %in% cv_countries) %>%
    group_by(country) %>%
    group_modify(~ run_country_cv_generic(.x, learners[[nm]])) %>%
    ungroup() %>%
    select(country, admin1, admin2, pred)
  names(res)[names(res) == "pred"] <- paste0("pred_", nm)
  cv_results <- left_join(cv_results, res, by = c("country", "admin1", "admin2"))
}

# ── Step 2c: LOCO -- all 4 countries. ────────────────────────────────────────
run_loco_generic <- function(data, learner) {
  bind_rows(lapply(unique(data$country), function(held_out) {
    train <- data[data$country != held_out, ]
    test  <- data[data$country == held_out, ]
    model <- learner$fit(train)
    test$pred <- learner$predict(model, test)
    test
  }))
}

loco_results <- data_imputed %>% select(country, admin1, admin2, prev_women_iron, n_women_iron)

for (nm in names(learners)) {
  res <- run_loco_generic(data_imputed, learners[[nm]]) %>%
    select(country, admin1, admin2, pred)
  names(res)[names(res) == "pred"] <- paste0("pred_", nm)
  loco_results <- left_join(loco_results, res, by = c("country", "admin1", "admin2"))
}

# ── Step 2d: does qb avoid the Task 2 negative-prediction pathology
# specifically under LOCO transport, not just in-sample? ────────────────────
cat("\nLOCO predicted range and out-of-[0,1] counts:\n")
cat("  lm:", range(loco_results$pred_lm), " -- n outside [0,1]:", sum(loco_results$pred_lm < 0 | loco_results$pred_lm > 1), "\n")
cat("  qb:", range(loco_results$pred_qb), " -- n outside [0,1]:", sum(loco_results$pred_qb < 0 | loco_results$pred_qb > 1), "\n")

# =============================================================================
# Step 3: threshold sweep under a false-alarm cap -- Task 1's Step 6, applied
# to both models under both evaluation setups. This is the decision-relevant
# comparison: does the choice of lm vs. quasi-binomial change WHICH districts
# get flagged high-burden, not just the raw predicted values.
# =============================================================================

evaluate_threshold <- function(data, pred_col, actual_col, frac) {
  preds  <- data[[pred_col]]
  actual <- data[[actual_col]]
  cutoff <- ave(preds, data$country, FUN = function(x) quantile(x, 1 - frac))
  predicted_flag <- preds >= cutoff
  data.frame(
    frac        = frac,
    sensitivity = sum(actual & predicted_flag)  / sum(actual),
    false_alarm = sum(!actual & predicted_flag) / sum(!actual)
  )
}

best_threshold <- function(sweep, cap) {
  under_cap <- sweep[sweep$false_alarm <= cap, ]
  if (nrow(under_cap) == 0) return(data.frame(frac = NA_real_, sensitivity = -Inf, false_alarm = NA_real_))
  under_cap[which.max(under_cap$sensitivity), ]
}

fractions <- seq(0.05, 0.50, by = 0.05)
false_alarm_cap <- 0.20
high_burden_quantile <- 0.80

# Same top-20%-within-country target as Task 1, computed once so both models
# are scored against the identical target.
cv_results <- cv_results %>%
  group_by(country) %>%
  mutate(actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

loco_results <- loco_results %>%
  group_by(country) %>%
  mutate(actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

cat("\nBest flagged fraction under a", false_alarm_cap * 100, "% false-alarm cap -- within-country CV:\n")
for (nm in names(learners)) {
  sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = cv_results, pred_col = paste0("pred_", nm), actual_col = "actual_high_burden"))
  cat(" ", nm, ":\n")
  print(best_threshold(sweep, false_alarm_cap), row.names = FALSE)
}

cat("\nBest flagged fraction under a", false_alarm_cap * 100, "% false-alarm cap -- LOCO:\n")
for (nm in names(learners)) {
  sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = loco_results, pred_col = paste0("pred_", nm), actual_col = "actual_high_burden"))
  cat(" ", nm, ":\n")
  print(best_threshold(sweep, false_alarm_cap), row.names = FALSE)
}

# ── Step 3a: top-20%-within-country flag match rate, lm vs. qb (each
# model's OWN quantile, not the sweep's tuned fraction) -- the simplest,
# most direct "do decisions change" check, both evaluation setups. ──────────
cat("\nTop-20%-within-country flag match rate, lm vs. qb (each model's own quantile):\n")
cat("Within-country CV:\n")
cv_results %>%
  group_by(country) %>%
  summarise(pct_match = mean(
    (pred_lm >= quantile(pred_lm, high_burden_quantile)) ==
    (pred_qb >= quantile(pred_qb, high_burden_quantile))
  )) %>%
  print()

cat("LOCO:\n")
loco_results %>%
  group_by(country) %>%
  summarise(pct_match = mean(
    (pred_lm >= quantile(pred_lm, high_burden_quantile)) ==
    (pred_qb >= quantile(pred_qb, high_burden_quantile))
  )) %>%
  print()

# =============================================================================
# Step 4: summary table -- Task 3 at a glance
# =============================================================================

task3_summary <- bind_rows(lapply(names(learners), function(nm) {
  cv_sweep   <- bind_rows(lapply(fractions, evaluate_threshold, data = cv_results,   pred_col = paste0("pred_", nm), actual_col = "actual_high_burden"))
  loco_sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = loco_results, pred_col = paste0("pred_", nm), actual_col = "actual_high_burden"))
  bind_rows(
    cbind(model = nm, evaluation = "Within-country CV", best_threshold(cv_sweep, false_alarm_cap)),
    cbind(model = nm, evaluation = "LOCO",              best_threshold(loco_sweep, false_alarm_cap))
  )
}))

cat("\n===== Task 3 summary: sensitivity/false-alarm under a", false_alarm_cap * 100, "% cap, lm vs. qb =====\n")
print(task3_summary, row.names = FALSE)

bounded_support_summary <- data.frame(
  model      = rep(c("lm", "qb"), 2),
  evaluation = rep(c("In-sample", "LOCO"), each = 2),
  n_out_of_range = c(
    sum(data_imputed$pred_lm < 0 | data_imputed$pred_lm > 1),
    sum(data_imputed$pred_qb < 0 | data_imputed$pred_qb > 1),
    sum(loco_results$pred_lm < 0 | loco_results$pred_lm > 1),
    sum(loco_results$pred_qb < 0 | loco_results$pred_qb > 1)
  )
)

cat("\n===== Task 3 summary: out-of-[0,1] predictions, lm vs. qb =====\n")
print(bounded_support_summary, row.names = FALSE)

flag_match_summary <- bind_rows(
  cv_results %>% mutate(evaluation = "Within-country CV"),
  loco_results %>% mutate(evaluation = "LOCO")
) %>%
  group_by(evaluation, country) %>%
  summarise(top20_flag_match = mean(
    (pred_lm >= quantile(pred_lm, high_burden_quantile)) ==
    (pred_qb >= quantile(pred_qb, high_burden_quantile))
  ), .groups = "drop")

cat("\n===== Task 3 summary: top-20% flag match rate, lm vs. qb, by country =====\n")
print(flag_match_summary, n = Inf)

# =============================================================================
# Step 5: conclusions and takeaways
#
# CALIBRATION: qb wins cleanly. lm produces 7 impossible predictions
# in-sample, 14 under LOCO; qb produces zero in either setting, by
# construction of the logit link -- the exact failure Task 2's Sierra Leone
# investigation was built around, closed off with no correction step needed.
#
# DECISION PERFORMANCE: qb does not win -- within-country CV it's clearly
# worse (34.1% vs. lm's 48.8% sensitivity at a comparable false-alarm rate),
# roughly a wash under LOCO. Why: lm's impossible predictions sat in LOW-
# prevalence districts, never the ones competing for the top-20% flag, so
# fixing them didn't help targeting -- while the logit link's bounded shape
# costs some of the fitting flexibility that lm had at the HIGH end, which
# is what sensitivity-under-a-cap rewards. Flags disagree 5-14% of the time
# (worst for Sierra Leone, unsurprisingly).
#
# This is Task 2's conclusion again, by a different mechanism: a fix for
# calibration is not a fix for ranking, because they're governed by
# different parts of the fitted curve. That lever is still the learner and
# how it's fit -- Task 1's territory.
#
# UNFINISHED, both flagged for a possible revisit rather than resolved here:
#  - Precision weighting (1/var_women_iron) was guarded, not fixed -- 71 of
#    206 districts have var=0 for this outcome.
#  - This comparison used single lm vs. single qb under plain CV/LOCO,
#    deliberately skipping Task 1's multi-learner blend + nested-CV honesty
#    layer (scoped out earlier as answering a different question). Worth
#    redoing with qb as one learner in that ensemble, sensitivity-tuned
#    weights, and a nested-CV check -- the within-country sensitivity gap
#    found here might shrink, grow, or reverse once qb isn't evaluated
#    alone against lm alone.
# =============================================================================
