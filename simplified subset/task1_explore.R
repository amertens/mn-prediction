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
women_iron <- women_iron %>%
  group_by(country) %>%
  mutate(high_burden = prev_women_iron >= quantile(prev_women_iron, 0.80)) %>%
  ungroup()

cat("\nHigh-burden districts flagged, by country:\n")
women_iron %>%
  group_by(country) %>%
  summarise(
    n_districts   = n(),
    n_high_burden = sum(high_burden),
    cutoff_prev   = quantile(prev_women_iron, 0.80)
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
      actual_high_burden    = prev_women_iron >= quantile(prev_women_iron, 0.80),
      predicted_high_burden = pred_baseline   >= quantile(pred_baseline, 0.80)
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
assign_folds <- function(n, k, seed = 1) {
  set.seed(seed)
  sample(rep(1:k, length.out = n))
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
    actual_high_burden    = prev_women_iron >= quantile(prev_women_iron, 0.80),
    predicted_high_burden = pred_cv         >= quantile(pred_cv, 0.80)
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
