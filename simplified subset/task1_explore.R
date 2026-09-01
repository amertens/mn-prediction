# Task 1, step 1: look at the women_iron columns in mn_admin2.csv before
# doing anything else.

library(dplyr)

admin2 <- read.csv("simplified subset/data/mn_admin2.csv")

women_iron <- admin2 %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron)

cat("Rows total:", nrow(admin2), "\n")
cat("Rows with non-missing prev_women_iron:", sum(!is.na(women_iron$prev_women_iron)), "\n\n")

cat("By country:\n")
women_iron %>%
  group_by(country) %>%
  summarise(
    n_districts = n(),
    n_with_data = sum(!is.na(prev_women_iron)),
    mean_prev   = mean(prev_women_iron, na.rm = TRUE),
    min_prev    = min(prev_women_iron, na.rm = TRUE),
    max_prev    = max(prev_women_iron, na.rm = TRUE)
  ) %>%
  print()

cat("\nFirst few rows:\n")
print(head(women_iron, 10))

# ── Step 2: define the target — top 20% highest-prevalence districts,
# within each country separately ────────────────────────────────────────────
# Single source of truth for "high-burden" -- every later step (predicted-
# flag defaults, actual_high_burden, fold stratification) references this
# constant instead of repeating 0.80, so changing it only requires an edit
# here.
high_burden_quantile <- 0.80

women_iron <- women_iron %>%
  group_by(country) %>%
  mutate(high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

cat("\nHigh-burden districts flagged, by country:\n")
women_iron %>%
  group_by(country) %>%
  summarise(
    n_districts   = n(),
    n_high_burden = sum(high_burden),
    cutoff_prev   = quantile(prev_women_iron, high_burden_quantile)
  ) %>%
  print()

# ── Step 3: fit a plain baseline model — linear regression of
# prev_women_iron on the 16 proxy predictors, weighted by sample size ───────
dict <- read.csv("simplified subset/data_dictionary.csv")
predictors <- dict$variable[dict$role == "predictor"]

model_data <- admin2 %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors))

baseline_formula <- reformulate(predictors, response = "prev_women_iron")

# Missingness: 10 of 206 districts have 1+ missing predictor value (see
# na_counts check). Handle it two ways and compare: drop those rows, or
# mean-impute them. `impute_predictors()` takes a `method` so we can swap in
# a different strategy later without touching the rest of the pipeline.
impute_predictors <- function(data, cols, method = c("mean", "median")) {
  method <- match.arg(method)
  fn <- switch(method, mean = mean, median = median)
  for (col in cols) {
    data[[col]][is.na(data[[col]])] <- fn(data[[col]], na.rm = TRUE)
  }
  data
}

# Fits the baseline model on `data` and returns it with predictions + the
# top-20%-within-country high-burden flags (actual and predicted) attached.
fit_baseline <- function(data) {
  fit <- lm(baseline_formula, data = data, weights = n_women_iron)
  data %>%
    mutate(pred_baseline = predict(fit, newdata = data)) %>%
    group_by(country) %>%
    mutate(
      actual_high_burden    = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile),
      predicted_high_burden = pred_baseline   >= quantile(pred_baseline, high_burden_quantile)
    ) %>%
    ungroup()
}

data_complete <- model_data[complete.cases(model_data[predictors]), ]
data_imputed  <- impute_predictors(model_data, predictors, method = "mean")

cat("\nDropped-rows approach (", nrow(data_complete), "of", nrow(model_data), "districts):\n")
result_complete <- fit_baseline(data_complete)
print(table(actual = result_complete$actual_high_burden,
            predicted = result_complete$predicted_high_burden))

cat("\nMean-imputed approach (all", nrow(data_imputed), "districts):\n")
result_imputed <- fit_baseline(data_imputed)
print(table(actual = result_imputed$actual_high_burden,
            predicted = result_imputed$predicted_high_burden))

# =============================================================================
# Step 4: within-country cross-validation
#
# Steps 1-3 evaluated the model on the same data it was fit on (in-sample) --
# too optimistic. Now: for each country, split its districts into folds,
# fit on the other folds, predict the held-out fold, rotate. This tests how
# well the model does inside a country that already has survey data.
#
# Sierra Leone has only 14 districts, fewer than the 17 parameters (16
# predictors + intercept) a plain linear model needs -- see the CV/rank-
# deficiency discussion. So: CV runs for Gambia/Ghana/Malawi (4a-4c), and
# Sierra Leone is shown separately (4d) as a concrete illustration of the
# degenerate fit, not included in the CV results.
#
# Uses data_imputed (all 206 districts, no dropped rows) so country group
# sizes stay at their full count going into the folds.
# =============================================================================

# ── Step 4a: assign folds within a country ──────────────────────────────────
# strata (optional): a vector the same length as n (e.g. TRUE/FALSE for
# high-burden). When given, folds are balanced within each stratum
# separately, so every fold gets its fair share of positives instead of
# risking a fold with zero (see Step 7's stratified-CV discussion). Default
# (strata = NULL) is plain random assignment, unchanged from Step 4/5.
assign_folds <- function(n, k, seed = 1, strata = NULL) {
  set.seed(seed)
  if (is.null(strata)) return(sample(rep(1:k, length.out = n)))
  folds <- integer(n)
  for (s in unique(strata)) {
    idx <- which(strata == s)
    folds[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  folds
}

# ── Step 4b: run k-fold CV for one country, return out-of-fold predictions ──
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

cv_countries <- c("Gambia", "Ghana", "Malawi")
cv_results <- data_imputed %>%
  filter(country %in% cv_countries) %>%
  group_by(country) %>%
  group_modify(~ run_country_cv(.x)) %>%
  ungroup()

# ── Step 4c: evaluate -- confusion matrix on out-of-fold predictions ────────
cv_results <- cv_results %>%
  group_by(country) %>%
  mutate(
    actual_high_burden    = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile),
    predicted_high_burden = pred_cv         >= quantile(pred_cv, high_burden_quantile)
  ) %>%
  ungroup()

cat("\nWithin-country CV (Gambia, Ghana, Malawi), out-of-fold confusion matrix:\n")
print(table(actual = cv_results$actual_high_burden,
            predicted = cv_results$predicted_high_burden))

# ── Step 4d: Sierra Leone -- illustrate the rank-deficient fit ──────────────
sl_data <- data_imputed %>% filter(country == "Sierra Leone")
sl_fit  <- lm(baseline_formula, data = sl_data, weights = n_women_iron)

cat("\nSierra Leone (14 districts, 17 parameters) -- full-sample fit, no CV:\n")
cat("Residual degrees of freedom:", sl_fit$df.residual, "\n")
cat("Coefficients R could not estimate (aliased, shown as NA):\n")
print(names(coef(sl_fit))[is.na(coef(sl_fit))])

# =============================================================================
# Step 5: leave-one-country-out (LOCO) evaluation
#
# Train on 3 countries, predict the 4th (held out entirely), rotate through
# all 4. Unlike Step 4, Sierra Leone is included here as a held-out test set
# -- the training data is the other 3 countries pooled (~176-192 rows), so
# there's no rank-deficiency problem; SL is only ever predicted on, never
# fit on alone.
# =============================================================================

# ── Step 5a: fit on 3 countries, predict the 4th, for every choice of held-out
run_loco <- function(data) {
  all_countries <- unique(data$country)
  results <- lapply(all_countries, function(held_out) {
    train <- data[data$country != held_out, ]
    test  <- data[data$country == held_out, ]
    fit <- lm(baseline_formula, data = train, weights = n_women_iron)
    test$pred_loco <- predict(fit, newdata = test)
    test
  })
  bind_rows(results)
}

loco_results <- run_loco(data_imputed)

# ── Step 5b: evaluate -- confusion matrix on held-out-country predictions ───
loco_results <- loco_results %>%
  group_by(country) %>%
  mutate(
    actual_high_burden    = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile),
    predicted_high_burden = pred_loco       >= quantile(pred_loco, high_burden_quantile)
  ) %>%
  ungroup()

cat("\nLOCO (all 4 countries), out-of-country confusion matrix:\n")
print(table(actual = loco_results$actual_high_burden,
            predicted = loco_results$predicted_high_burden))

cat("\nLOCO confusion matrix, broken out by held-out country:\n")
for (c in unique(loco_results$country)) {
  cat("\n", c, ":\n", sep = "")
  print(table(actual = loco_results$actual_high_burden[loco_results$country == c],
              predicted = loco_results$predicted_high_burden[loco_results$country == c]))
}

# =============================================================================
# Step 6: threshold sweep
#
# So far "predicted high-burden" always meant "predicted top 20% within
# country" -- matching the flagged rate to the target's own base rate. That
# was one arbitrary choice among many. Here: sweep the flagged fraction (5%
# to 50% of districts, within country), trace out sensitivity vs.
# false-alarm rate at each, and pick the fraction that maximizes sensitivity
# while keeping the false-alarm rate under a cap (20% here, arbitrary --
# easy to change). Uses the CV/LOCO predictions already computed in Steps
# 4-5, no refitting needed.
# =============================================================================

# ── Step 6a: pooled sensitivity/false-alarm at one flagged fraction ─────────
evaluate_threshold <- function(data, pred_col, frac) {
  preds  <- data[[pred_col]]
  cutoff <- ave(preds, data$country, FUN = function(x) quantile(x, 1 - frac))
  predicted_flag <- preds >= cutoff
  actual <- data$actual_high_burden

  data.frame(
    frac        = frac,
    sensitivity = sum(actual & predicted_flag)  / sum(actual),
    false_alarm = sum(!actual & predicted_flag) / sum(!actual)
  )
}

fractions <- seq(0.05, 0.50, by = 0.05)

# ── Step 6b: sweep for within-country CV predictions (Gambia/Ghana/Malawi) ──
cv_sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = cv_results, pred_col = "pred_cv"))
cat("\nWithin-country CV: sensitivity / false-alarm by flagged fraction:\n")
print(cv_sweep, row.names = FALSE)

# ── Step 6c: sweep for LOCO predictions (all 4 countries) ───────────────────
loco_sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = loco_results, pred_col = "pred_loco"))
cat("\nLOCO: sensitivity / false-alarm by flagged fraction:\n")
print(loco_sweep, row.names = FALSE)

# ── Step 6d: pick the best fraction under a false-alarm cap ─────────────────
false_alarm_cap <- 0.20

best_threshold <- function(sweep, cap) {
  under_cap <- sweep[sweep$false_alarm <= cap, ]
  if (nrow(under_cap) == 0) {
    # No fraction in the grid satisfies the cap -- e.g. a degenerate score
    # (extreme weight combo, ties) on a small nested-inner subset. Return a
    # sensitivity of -Inf so this combo is never picked by which.max()
    # upstream, rather than crashing.
    return(data.frame(frac = NA_real_, sensitivity = -Inf, false_alarm = NA_real_))
  }
  under_cap[which.max(under_cap$sensitivity), ]
}

cv_best_lm   <- best_threshold(cv_sweep, false_alarm_cap)
loco_best_lm <- best_threshold(loco_sweep, false_alarm_cap)

cat("\nBest flagged fraction under a", false_alarm_cap * 100, "% false-alarm cap:\n")
cat("Within-country CV:\n")
print(cv_best_lm, row.names = FALSE)
cat("LOCO:\n")
print(loco_best_lm, row.names = FALSE)

# =============================================================================
# Step 7: combined learner -- three base learners blended with sensitivity-
# tuned combination weights, following Zheng et al. 2018 (Stat Med) "Constrained
# binary classification using ensemble learning" -- see memory reference
# reference_zheng_sensitivity_constrained_sl.
#
# Adds glmnet (elastic net) and ranger (random forest) alongside the linear
# model. Steps 4/5's lm-only CV/LOCO code is left as-is; this section adds
# its own generic versions (parameterized by learner) so all three learners
# reuse the same CV/LOCO logic instead of three copy-pasted loops.
#
# Two deliberate departures from Zheng et al.'s exact recipe, both driven by
# our small per-country sample sizes:
#   1. Constraint: we bound the FALSE-ALARM RATE (P(flagged | actually low-
#      burden)), matching the plan's wording ("false-alarm rate below a set
#      limit") and their section 3.2 general/Neyman-Pearson case -- not their
#      worked example's RPP (P(flagged) unconditionally), which is a
#      different quantity.
#   2. Threshold timing: Zheng et al.'s CV objective picks a threshold
#      SEPARATELY on each fold's validation data, then averages the objective
#      across folds. Our CV folds are tiny (a Gambia test fold has ~6
#      districts, ~1 of them high-burden) -- a per-fold sensitivity threshold
#      from 1-2 positives would be too noisy to mean anything. Instead we
#      pool all out-of-fold predictions (across folds/countries) first, then
#      sweep ONE threshold against the pooled set. This is an approximation
#      of their CV objective, not their literal procedure.
#
# The three learners' out-of-fold predictions are fit ONCE (7a-7c), then
# combined via combined = w_lm*pred_lm + w_glmnet*pred_glmnet + w_ranger*pred_ranger
# (weights >= 0, summing to 1) -- Zheng et al.'s Psi_alpha (their eq. in
# section 2.2.2). We search the weight grid + threshold jointly (their eq. 9,
# with a grid search standing in for their nloptr optimizer) to find the
# combination that maximizes sensitivity under the false-alarm cap.
# =============================================================================

# ── Step 7a: define each learner as a fit/predict pair ──────────────────────
learner_lm <- list(
  fit     = function(train) lm(baseline_formula, data = train, weights = n_women_iron),
  predict = function(model, test) predict(model, newdata = test)
)

learner_glmnet <- list(
  fit = function(train) {
    x <- as.matrix(train[predictors])
    nfolds <- max(3, min(10, floor(nrow(train) / 3)))  # small training sets need fewer internal folds
    glmnet::cv.glmnet(x, train$prev_women_iron, weights = train$n_women_iron,
                       alpha = 0.5, nfolds = nfolds)
  },
  predict = function(model, test) {
    as.numeric(predict(model, newx = as.matrix(test[predictors]), s = "lambda.min"))
  }
)

learner_ranger <- list(
  fit     = function(train) ranger::ranger(baseline_formula, data = train, case.weights = train$n_women_iron),
  predict = function(model, test) predict(model, data = test)$predictions
)

learners <- list(lm = learner_lm, glmnet = learner_glmnet, ranger = learner_ranger)

# Toggle: TRUE balances high-burden districts across folds (avoids a fold
# like Gambia's fold 4 having zero positives); FALSE is plain random
# assignment. Easy to flip and re-run.
stratify_cv_folds <- TRUE

# ── Step 7b: generic CV/LOCO runners, parameterized by learner ──────────────
# folds (optional): a pre-computed fold-label vector to reuse instead of
# assigning fresh ones -- used by Step 8's nested CV to run this same
# function as the INNER rotation over an outer-train subset's own labels.
run_country_cv_generic <- function(data, learner, k = 5, stratify = stratify_cv_folds, folds = NULL) {
  if (is.null(folds)) {
    strata <- if (stratify) data$prev_women_iron >= quantile(data$prev_women_iron, high_burden_quantile) else NULL
    folds <- assign_folds(nrow(data), k, strata = strata)
  }
  data$fold <- folds
  data$pred <- NA_real_
  for (f in sort(unique(data$fold))) {
    train <- data[data$fold != f, ]
    test  <- data[data$fold == f, ]
    model <- learner$fit(train)
    data$pred[data$fold == f] <- learner$predict(model, test)
  }
  data
}

run_loco_generic <- function(data, learner) {
  results <- lapply(unique(data$country), function(held_out) {
    train <- data[data$country != held_out, ]
    test  <- data[data$country == held_out, ]
    model <- learner$fit(train)
    test$pred <- learner$predict(model, test)
    test
  })
  bind_rows(results)
}

# ── Step 7c: run all three learners once, collect predictions side by side ──
cv_multi <- data_imputed %>%
  filter(country %in% cv_countries) %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron)

loco_multi <- data_imputed %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron)

for (nm in names(learners)) {
  cv_res <- data_imputed %>%
    filter(country %in% cv_countries) %>%
    group_by(country) %>%
    group_modify(~ run_country_cv_generic(.x, learners[[nm]])) %>%
    ungroup() %>%
    select(country, admin1, admin2, pred)
  names(cv_res)[names(cv_res) == "pred"] <- paste0("pred_", nm)
  cv_multi <- left_join(cv_multi, cv_res, by = c("country", "admin1", "admin2"))

  loco_res <- run_loco_generic(data_imputed, learners[[nm]]) %>%
    select(country, admin1, admin2, pred)
  names(loco_res)[names(loco_res) == "pred"] <- paste0("pred_", nm)
  loco_multi <- left_join(loco_multi, loco_res, by = c("country", "admin1", "admin2"))
}

cv_multi <- cv_multi %>%
  group_by(country) %>%
  mutate(actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

loco_multi <- loco_multi %>%
  group_by(country) %>%
  mutate(actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

# ── Step 7d: search combination weights + threshold jointly ─────────────────
# Grid over weights (must be >= 0 and sum to 1); for each weight combo, sweep
# thresholds and keep the best one under the false-alarm cap; then take the
# best weight combo overall.
simplex_grid <- function(step = 0.1) {
  vals <- seq(0, 1, by = step)
  grid <- expand.grid(w_lm = vals, w_glmnet = vals, w_ranger = vals)
  grid[abs(rowSums(grid) - 1) < 1e-8, ]
}

search_combination_weights <- function(data, pred_cols, grid, fractions, cap) {
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    w <- grid[i, ]
    data$combined <- w$w_lm * data[[pred_cols[1]]] +
      w$w_glmnet * data[[pred_cols[2]]] +
      w$w_ranger * data[[pred_cols[3]]]
    sweep <- bind_rows(lapply(fractions, evaluate_threshold, data = data, pred_col = "combined"))
    cbind(w, best_threshold(sweep, cap))
  })
  bind_rows(rows)
}

pred_cols <- c("pred_lm", "pred_glmnet", "pred_ranger")
weight_grid <- simplex_grid(step = 0.1)

cv_weight_search   <- search_combination_weights(cv_multi,   pred_cols, weight_grid, fractions, false_alarm_cap)
loco_weight_search <- search_combination_weights(loco_multi, pred_cols, weight_grid, fractions, false_alarm_cap)

cv_best_combo   <- cv_weight_search[which.max(cv_weight_search$sensitivity), ]
loco_best_combo <- loco_weight_search[which.max(loco_weight_search$sensitivity), ]

cat("\nBest combination weights -- within-country CV (vs. lm-alone: 48.8% sensitivity at frac 0.25):\n")
print(cv_best_combo, row.names = FALSE)

cat("\nBest combination weights -- LOCO (vs. lm-alone: 29.5% sensitivity at frac 0.20):\n")
print(loco_best_combo, row.names = FALSE)

# =============================================================================
# Step 8: nested (outer) cross-validation
#
# Step 7's weight + threshold search is chosen to maximize sensitivity on the
# same pooled out-of-fold predictions we then report sensitivity from -- an
# honest first layer (model fitting) with an optimistic second layer (weight
# and threshold selection) stacked on top.
#
# Nested CV removes the optimism: for each outer fold/country, hold it out
# completely; do Step 7's ENTIRE process (inner CV + weight/threshold search)
# using only the remaining data; apply the resulting FIXED (weights,
# threshold) -- decided without ever seeing the outer fold -- to score it
# once. Rotate. The outer fold's own labels never influence what gets
# applied to it.
#
# Implementation reuses Step 7's machinery for the inner rotation: for CV,
# the same stratified 5-way `outer_fold` labels are reused as the inner
# folds too (run_country_cv_generic's `folds` argument, restricted to
# outer-train); for LOCO, the inner rotation is just run_loco_generic()
# called on the 3 outer-train countries.
# =============================================================================

# ── Step 8a: persist one fixed outer_fold assignment + actual_high_burden,
# computed ONCE over each country's full data (never recomputed on a subset,
# which would silently shift the target definition) ─────────────────────────
cv_nested_data <- data_imputed %>%
  filter(country %in% cv_countries) %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors)) %>%
  group_by(country) %>%
  mutate(
    actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile),
    outer_fold = assign_folds(n(), 5, strata = if (stratify_cv_folds) actual_high_burden else NULL)
  ) %>%
  ungroup()

loco_nested_data <- data_imputed %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors)) %>%
  group_by(country) %>%
  mutate(actual_high_burden = prev_women_iron >= quantile(prev_women_iron, high_burden_quantile)) %>%
  ungroup()

# ── Step 8b: final-refit helper for CV (per-country, like Step 7) ───────────
fit_predict_by_country <- function(train, test, learner) {
  bind_rows(lapply(unique(test$country), function(c) {
    tr <- train[train$country == c, ]
    te <- test[test$country == c, ]
    model <- learner$fit(tr)
    te$pred <- learner$predict(model, te)
    te
  }))
}

# ── Step 8c: nested within-country CV ────────────────────────────────────────
nested_within_country_cv <- function(data, learners, pred_cols, grid, fractions, cap) {
  outer_rows <- list()
  round_choices <- list()

  for (o in sort(unique(data$outer_fold))) {
    outer_train <- data[data$outer_fold != o, ]
    outer_test  <- data[data$outer_fold == o, ]

    # Inner CV per learner, reusing the SAME outer_fold labels restricted to
    # outer_train -- this is the "k-1 folds get their own rotation" step.
    inner_preds <- outer_train %>% select(country, admin1, admin2, prev_women_iron, n_women_iron, actual_high_burden)
    for (nm in names(learners)) {
      res <- outer_train %>%
        group_by(country) %>%
        group_modify(~ run_country_cv_generic(.x, learners[[nm]], folds = .x$outer_fold)) %>%
        ungroup() %>%
        select(country, admin1, admin2, pred)
      names(res)[names(res) == "pred"] <- paste0("pred_", nm)
      inner_preds <- left_join(inner_preds, res, by = c("country", "admin1", "admin2"))
    }

    # Honest weight + threshold search -- never sees outer_test.
    search <- search_combination_weights(inner_preds, pred_cols, grid, fractions, cap)
    best <- search[which.max(search$sensitivity), ]
    round_choices[[as.character(o)]] <- best

    # Final refit on ALL of outer_train, predict the untouched outer_test once.
    outer_test_preds <- outer_test %>% select(country, admin1, admin2, prev_women_iron, n_women_iron, actual_high_burden)
    for (nm in names(learners)) {
      res <- fit_predict_by_country(outer_train, outer_test, learners[[nm]]) %>%
        select(country, admin1, admin2, pred)
      names(res)[names(res) == "pred"] <- paste0("pred_", nm)
      outer_test_preds <- left_join(outer_test_preds, res, by = c("country", "admin1", "admin2"))
    }

    outer_test_preds$combined <- best$w_lm     * outer_test_preds[[pred_cols[1]]] +
                                  best$w_glmnet * outer_test_preds[[pred_cols[2]]] +
                                  best$w_ranger * outer_test_preds[[pred_cols[3]]]
    outer_test_preds <- outer_test_preds %>%
      group_by(country) %>%
      mutate(predicted_high_burden = combined >= quantile(combined, 1 - best$frac)) %>%
      ungroup()

    outer_rows[[as.character(o)]] <- outer_test_preds
  }

  all_outer <- bind_rows(outer_rows)
  list(
    sensitivity   = sum(all_outer$actual_high_burden & all_outer$predicted_high_burden) / sum(all_outer$actual_high_burden),
    false_alarm   = sum(!all_outer$actual_high_burden & all_outer$predicted_high_burden) / sum(!all_outer$actual_high_burden),
    round_choices = bind_rows(round_choices, .id = "outer_fold")
  )
}

# ── Step 8d: nested LOCO ─────────────────────────────────────────────────────
nested_loco <- function(data, learners, pred_cols, grid, fractions, cap) {
  outer_rows <- list()
  round_choices <- list()

  for (held_out in unique(data$country)) {
    outer_train <- data[data$country != held_out, ]
    outer_test  <- data[data$country == held_out, ]

    # Inner LOCO among just the outer-train countries.
    inner_preds <- outer_train %>% select(country, admin1, admin2, prev_women_iron, n_women_iron, actual_high_burden)
    for (nm in names(learners)) {
      res <- run_loco_generic(outer_train, learners[[nm]]) %>%
        select(country, admin1, admin2, pred)
      names(res)[names(res) == "pred"] <- paste0("pred_", nm)
      inner_preds <- left_join(inner_preds, res, by = c("country", "admin1", "admin2"))
    }

    search <- search_combination_weights(inner_preds, pred_cols, grid, fractions, cap)
    best <- search[which.max(search$sensitivity), ]
    round_choices[[held_out]] <- best

    # Final refit on all of outer_train (pooled across its 3 countries),
    # predict the untouched held-out country once.
    outer_test_preds <- outer_test %>% select(country, admin1, admin2, prev_women_iron, n_women_iron, actual_high_burden)
    for (nm in names(learners)) {
      model <- learners[[nm]]$fit(outer_train)
      outer_test_preds[[paste0("pred_", nm)]] <- learners[[nm]]$predict(model, outer_test)
    }

    outer_test_preds$combined <- best$w_lm     * outer_test_preds[[pred_cols[1]]] +
                                  best$w_glmnet * outer_test_preds[[pred_cols[2]]] +
                                  best$w_ranger * outer_test_preds[[pred_cols[3]]]
    outer_test_preds$predicted_high_burden <-
      outer_test_preds$combined >= quantile(outer_test_preds$combined, 1 - best$frac)

    outer_rows[[held_out]] <- outer_test_preds
  }

  all_outer <- bind_rows(outer_rows)
  list(
    sensitivity   = sum(all_outer$actual_high_burden & all_outer$predicted_high_burden) / sum(all_outer$actual_high_burden),
    false_alarm   = sum(!all_outer$actual_high_burden & all_outer$predicted_high_burden) / sum(!all_outer$actual_high_burden),
    round_choices = bind_rows(round_choices, .id = "held_out_country")
  )
}

# ── Step 8e: run both, compare to Step 7's (optimistic) pooled numbers ──────
nested_cv_result   <- nested_within_country_cv(cv_nested_data,   learners, pred_cols, weight_grid, fractions, false_alarm_cap)
nested_loco_result <- nested_loco(loco_nested_data, learners, pred_cols, weight_grid, fractions, false_alarm_cap)

cat("\nNested within-country CV -- honest sensitivity:", nested_cv_result$sensitivity,
    " false_alarm:", nested_cv_result$false_alarm,
    " (Step 7 pooled estimate: 51.2% / 18.5%)\n")
cat("Per-outer-fold weight/threshold choices (stability check):\n")
print(nested_cv_result$round_choices, row.names = FALSE)

cat("\nNested LOCO -- honest sensitivity:", nested_loco_result$sensitivity,
    " false_alarm:", nested_loco_result$false_alarm,
    " (Step 7 pooled estimate: 38.6% / 15.4%)\n")
cat("Per-outer-country weight/threshold choices (stability check):\n")
print(nested_loco_result$round_choices, row.names = FALSE)

# =============================================================================
# Step 9: a second target definition -- WHO public-health threshold
#
# The plan asks for TWO ways to define "high-burden" (its wording: districts
# "above a WHO public-health threshold, you can pick an arbitrary cutoff at
# this point"). Steps 2-8 used top-20%-within-country, which always flags
# exactly 20% of every country's districts regardless of that country's
# actual level. This is the alternative: a single FIXED cutoff (20%, chosen
# as the analog to WHO's anemia "moderate public health significance" band --
# see the false_alarm_cap discussion for why 20% and not 40%) applied the
# same way everywhere, so a low-burden country can come back with few or zero
# flagged districts and a high-burden one (Gambia: 24/30 districts >= 20%)
# can come back almost entirely flagged.
#
# Reuses every function from Steps 6-8 unchanged (evaluate_threshold,
# best_threshold, search_combination_weights, nested_within_country_cv,
# nested_loco) -- they all take `data` as an argument and read
# data$actual_high_burden off it, so relabeling that one column on copies of
# Step 7's already-fitted prediction sets is the only new work needed; no
# model is refit for this section.
# =============================================================================

who_threshold <- 0.20

label_who_high_burden <- function(data) {
  data$actual_high_burden <- data$prev_women_iron >= who_threshold
  data
}

cat("\nWHO-threshold (>=", who_threshold * 100, "%) districts flagged, by country:\n")
data_imputed %>%
  group_by(country) %>%
  summarise(n_districts = n(), n_high_burden = sum(prev_women_iron >= who_threshold)) %>%
  print()

# ── Step 9a: relabel Step 7's already-fitted CV/LOCO prediction sets ────────
cv_multi_who   <- label_who_high_burden(cv_multi)
loco_multi_who <- label_who_high_burden(loco_multi)

# Plain-lm comparator under the WHO definition, same shape as Step 6d.
cv_sweep_lm_who   <- bind_rows(lapply(fractions, evaluate_threshold, data = cv_multi_who,   pred_col = "pred_lm"))
loco_sweep_lm_who <- bind_rows(lapply(fractions, evaluate_threshold, data = loco_multi_who, pred_col = "pred_lm"))
cv_best_lm_who   <- best_threshold(cv_sweep_lm_who,   false_alarm_cap)
loco_best_lm_who <- best_threshold(loco_sweep_lm_who, false_alarm_cap)

# Blended-ensemble comparator, same weight/threshold search as Step 7d.
cv_weight_search_who   <- search_combination_weights(cv_multi_who,   pred_cols, weight_grid, fractions, false_alarm_cap)
loco_weight_search_who <- search_combination_weights(loco_multi_who, pred_cols, weight_grid, fractions, false_alarm_cap)
cv_best_combo_who   <- cv_weight_search_who[which.max(cv_weight_search_who$sensitivity), ]
loco_best_combo_who <- loco_weight_search_who[which.max(loco_weight_search_who$sensitivity), ]

cat("\n[WHO threshold] Plain lm -- within-country CV:\n")
print(cv_best_lm_who, row.names = FALSE)
cat("[WHO threshold] Plain lm -- LOCO:\n")
print(loco_best_lm_who, row.names = FALSE)
cat("[WHO threshold] Best combination weights -- within-country CV:\n")
print(cv_best_combo_who, row.names = FALSE)
cat("[WHO threshold] Best combination weights -- LOCO:\n")
print(loco_best_combo_who, row.names = FALSE)

# ── Step 9b: nested (honest) versions, mirroring Step 8 ─────────────────────
# Fresh outer_fold labels stratified on the WHO flag (not reused from Step
# 8a's quantile-stratified folds) -- the two target definitions disagree on
# which districts are positive, so a fold balanced for one is not balanced
# for the other.
cv_nested_data_who <- data_imputed %>%
  filter(country %in% cv_countries) %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors)) %>%
  label_who_high_burden() %>%
  group_by(country) %>%
  mutate(outer_fold = assign_folds(n(), 5, strata = if (stratify_cv_folds) actual_high_burden else NULL)) %>%
  ungroup()

loco_nested_data_who <- data_imputed %>%
  select(country, admin1, admin2, prev_women_iron, n_women_iron, all_of(predictors)) %>%
  label_who_high_burden()

nested_cv_result_who   <- nested_within_country_cv(cv_nested_data_who,   learners, pred_cols, weight_grid, fractions, false_alarm_cap)
nested_loco_result_who <- nested_loco(loco_nested_data_who, learners, pred_cols, weight_grid, fractions, false_alarm_cap)

cat("\n[WHO threshold] Nested within-country CV -- honest sensitivity:", nested_cv_result_who$sensitivity,
    " false_alarm:", nested_cv_result_who$false_alarm, "\n")
cat("[WHO threshold] Nested LOCO -- honest sensitivity:", nested_loco_result_who$sensitivity,
    " false_alarm:", nested_loco_result_who$false_alarm, "\n")

# =============================================================================
# Step 10: summary table -- Task 1's stated deliverable
#
# "Output could be a table or curve showing the share of high-burden
# districts correctly flagged, within a country and when transported to a
# new one, compared with an ordinary error-minimising model." One table,
# both target definitions, both evaluation setups, plain lm vs. the
# sensitivity-tuned blend. All numbers are out-of-fold/held-out (never fit on
# what they're scored against) but at "pooled" rigor -- honest about model
# fit, not about the threshold/weight-selection step layered on top. The
# nested (Step 8/9b) numbers, which correct for that remaining optimism, are
# reported separately below since they only exist for the blended model
# (plain lm has no selection step to be optimistic about).
# =============================================================================

summary_table <- tibble::tribble(
  ~target_definition,        ~evaluation,          ~model,     ~sensitivity,                     ~false_alarm,
  "Top 20% within country",  "Within-country CV",  "Plain lm", cv_best_lm$sensitivity,            cv_best_lm$false_alarm,
  "Top 20% within country",  "Within-country CV",  "Blended",  cv_best_combo$sensitivity,         cv_best_combo$false_alarm,
  "Top 20% within country",  "LOCO",               "Plain lm", loco_best_lm$sensitivity,          loco_best_lm$false_alarm,
  "Top 20% within country",  "LOCO",               "Blended",  loco_best_combo$sensitivity,       loco_best_combo$false_alarm,
  "WHO threshold (20%)",     "Within-country CV",  "Plain lm", cv_best_lm_who$sensitivity,        cv_best_lm_who$false_alarm,
  "WHO threshold (20%)",     "Within-country CV",  "Blended",  cv_best_combo_who$sensitivity,     cv_best_combo_who$false_alarm,
  "WHO threshold (20%)",     "LOCO",               "Plain lm", loco_best_lm_who$sensitivity,      loco_best_lm_who$false_alarm,
  "WHO threshold (20%)",     "LOCO",               "Blended",  loco_best_combo_who$sensitivity,   loco_best_combo_who$false_alarm
)

cat("\n===== Task 1 summary: share of high-burden districts caught =====\n")
print(summary_table, n = Inf)

nested_summary <- tibble::tribble(
  ~target_definition,       ~evaluation,          ~sensitivity,                       ~false_alarm,
  "Top 20% within country", "Within-country CV",  nested_cv_result$sensitivity,       nested_cv_result$false_alarm,
  "Top 20% within country", "LOCO",               nested_loco_result$sensitivity,     nested_loco_result$false_alarm,
  "WHO threshold (20%)",    "Within-country CV",  nested_cv_result_who$sensitivity,   nested_cv_result_who$false_alarm,
  "WHO threshold (20%)",    "LOCO",               nested_loco_result_who$sensitivity, nested_loco_result_who$false_alarm
)

cat("\n===== Honesty check: nested (outer) CV, blended model only =====\n")
print(nested_summary, n = Inf)

# =============================================================================
# Step 11: conclusions and takeaways
#
# The plan's framing for Task 1 is decision-focused targeting: not "predict
# prevalence well on average" but "catch the highest-burden districts without
# crying wolf too often." Getting there took four separate refinements, each
# stripping out one form of over-optimism left by the one before it. None of
# them is optional -- skip any one and the reported numbers overstate what
# the model would actually do on a new survey round.
#
#  1. WITHIN-COUNTRY CV vs. LOCO (Steps 4-5, 7, 9a). Same target, same
#     models, two honesty regimes. Within-country CV asks "can this model
#     find the worst districts in a country it has already seen data from?"
#     LOCO asks the harder, more decision-relevant question: "can it find
#     them in a country it has NEVER seen?" -- the actual use case for
#     countries without district-level survey data. LOCO sensitivity comes in
#     lower every time it's compared to within-country CV in this script; the
#     gap between the two IS the transportability cost, and reporting only
#     the within-country number would hide it.
#
#  2. SENSITIVITY UNDER A FALSE-ALARM CONSTRAINT, not raw accuracy (Steps
#     6-7). A model that minimizes overall prediction error is not the same
#     as one that's good at flagging the top of the distribution -- Step 3's
#     in-sample confusion matrices make that gap visible immediately. Once
#     the objective is reframed as "maximize sensitivity subject to
#     false_alarm <= cap" (Zheng et al.'s constrained-classification setup),
#     blending three learners with tuned combination weights beats the
#     single plain-lm model under that objective in the quantile-target
#     comparisons above -- but "better at this specific decision task" only
#     means something once the objective matches the actual decision.
#
#  3. NESTED CV so the threshold/weight search doesn't double-dip (Step 8,
#     9b). Steps 6-7's weight and threshold choices are themselves fit to
#     data -- picked to maximize sensitivity on the same pooled out-of-fold
#     predictions the sensitivity is then reported from. That's a second,
#     easy-to-miss layer of leakage on top of ordinary CV: honest about the
#     MODEL fit, optimistic about the SELECTION on top of it. Nesting an
#     outer hold-out around the entire inner process (fit + select) removes
#     that -- and the size of the correction is itself informative: it was
#     modest and stable for the quantile target, but for the WHO-threshold
#     LOCO comparison the nested honesty check strips away most of the
#     pooled estimate's apparent sensitivity, because with only 3 countries
#     to train on and wildly uneven WHO-flag base rates across them (Gambia
#     ~80% flagged vs. Ghana ~8%), the "best" weights swing hard depending on
#     which 3 countries happen to be in the training set. A model that looks
#     good un-nested here would not have looked good on the next country.
#
#  4. TWO TARGET DEFINITIONS, not one (Steps 2, 9). "High-burden" is a
#     choice, not a fact -- the plan explicitly asks for both a within-
#     country relative definition (top 20%) and a fixed public-health-style
#     cutoff (WHO threshold). They tell different stories: the quantile
#     target is stable across the honesty layers above (its LOCO sensitivity
#     barely moves once nested), while the WHO threshold's LOCO performance
#     is fragile -- exactly because a fixed cutoff, unlike a per-country
#     quantile, can and does produce lopsided class balance across countries.
#     Reporting a single number under a single target definition would have
#     hidden that fragility entirely.
#
# Put together, the throughline is that a "fit one model, report one
# sensitivity number" pass is never the full answer for a task like this --
# every added layer of rigor (transport, objective, honesty, target
# definition) moved the headline numbers, sometimes by a little (quantile
# target, within-country vs. LOCO) and sometimes drastically (WHO threshold,
# LOCO, nested vs. pooled). The four comparisons above are the model of how
# much confidence to place in any single reported number: check whether it
# survives transport, whether it's optimized for the actual decision, whether
# it survives an honest outer hold-out, and whether it survives a different
# reasonable choice of target definition. A number that survives all four is
# a much stronger claim than one that was only ever checked in-sample.
# =============================================================================
