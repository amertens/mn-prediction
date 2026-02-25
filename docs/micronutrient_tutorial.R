#############################################################
# Tutorial: Predicting Micronutrient Deficiency and Mapping
#
# This script demonstrates how to fit a predictive model for a
# binary micronutrient deficiency outcome using simulated
# individual‑level data.  It then shows how to aggregate
# predictions to the administrative level, estimate uncertainty
# for aggregated predictions, and transport the model to a new
# geographic area.  The analyses are illustrative and use
# simulated data to mimic vitamin A or iron deficiency
# studies.  Adjust the model specification, covariates, and
# aggregation functions as needed for real datasets.

## ------------------------------------------------------------------
## 1.  Generate simulated individual‑level data
##
## Each record represents a child with demographic covariates,
## membership in a Ghanaian administrative level (admin2) and a
## binary deficiency outcome.  The covariate distributions and
## coefficients used to generate the outcome are arbitrary and
## chosen for illustration only.

set.seed(123)
n <- 1000  # number of individuals

# Simulate the admin2 code (ten hypothetical districts)
admin2_levels <- paste0("admin2_", sprintf("%02d", 1:10))
admin2 <- sample(admin2_levels, size = n, replace = TRUE)

# Covariates: age in months, sex, rural residence, maternal education level
age_months     <- rpois(n, lambda = 60)
sex            <- rbinom(n, size = 1, prob = 0.5)       # 1 = male, 0 = female
rural          <- rbinom(n, size = 1, prob = 0.4)       # 1 = rural, 0 = urban
maternal_edu   <- sample(0:3, size = n, replace = TRUE) # 0 = none, 3 = tertiary

# Generate the log odds of deficiency using a simple linear model
logit_p <- -1.0 + 0.01 * age_months - 0.4 * sex - 0.6 * rural +
           0.3 * maternal_edu + rnorm(n, mean = 0, sd = 0.2)
prob_def <- 1 / (1 + exp(-logit_p))
deficiency <- rbinom(n, size = 1, prob = prob_def)

# Assemble the data frame
df <- data.frame(
  deficiency    = as.factor(deficiency),
  age_months    = age_months,
  sex           = factor(sex, levels = c(0, 1), labels = c("female", "male")),
  rural         = factor(rural, levels = c(0, 1), labels = c("urban", "rural")),
  maternal_edu  = factor(maternal_edu),
  admin2        = factor(admin2)
)

## ------------------------------------------------------------------
## 2.  Fit a predictive model for deficiency
##
## Use the caret package to train a logistic regression model with
## cross‑validation.  The outcome is binary; caret will internally
## convert factors to dummy variables.  Other algorithms (e.g.,
## random forest, gradient boosting) can be substituted by
## changing the ‘method’ argument.  Cross‑validation results can
## be inspected to assess model performance and calibrate
## hyperparameters.

library(caret)

# Define cross‑validation settings
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Recode the outcome to have descriptive class names (required by
# caret’s twoClassSummary)
df$deficiency <- relevel(df$deficiency, ref = "0")
levels(df$deficiency) <- c("no", "yes")

# Fit a logistic regression model
set.seed(20260224)
model_fit <- train(
  deficiency ~ age_months + sex + rural + maternal_edu,
  data = df,
  method = "glm",
  family = binomial(link = "logit"),
  trControl = train_control,
  metric = "ROC"
)

## Examine model performance
print(model_fit)

# Predict individual‑level probabilities
df$pred_prob <- predict(model_fit, newdata = df, type = "prob")[, "yes"]

## ------------------------------------------------------------------
## 3.  Aggregate predictions to administrative level 2 (admin2)
##
## To obtain district‑level prevalence estimates, compute the mean
## predicted probability for each admin2 region.  Estimate
## uncertainty via bootstrap sampling: repeatedly draw samples of
## individuals within each admin2 and recompute the mean predicted
## probability.  The bootstrap standard error and confidence
## interval provide a measure of uncertainty for the aggregated
## estimate.

library(dplyr)

aggregate_predictions <- function(data, n_boot = 500, conf_level = 0.95) {
  # Summarise predictions by admin2
  summary_df <- data %>%
    group_by(admin2) %>%
    summarise(
      pred_mean = mean(pred_prob, na.rm = TRUE),
      .groups = "drop"
    )

  # Bootstrap standard error and confidence intervals
  boot_est <- lapply(split(data, data$admin2), function(sub_df) {
    preds <- sub_df$pred_prob
    n     <- length(preds)
    boot_means <- replicate(n_boot, {
      idx <- sample(seq_len(n), size = n, replace = TRUE)
      mean(preds[idx], na.rm = TRUE)
    })
    # Compute standard error and CI
    se  <- sd(boot_means)
    alpha <- 1 - conf_level
    lower <- quantile(boot_means, probs = alpha / 2)
    upper <- quantile(boot_means, probs = 1 - alpha / 2)
    list(se = se, lower = lower, upper = upper)
  })

  summary_df$se    <- sapply(boot_est, `[[`, "se")
  summary_df$lower <- sapply(boot_est, `[[`, "lower")
  summary_df$upper <- sapply(boot_est, `[[`, "upper")
  summary_df
}

admin2_summary <- aggregate_predictions(df, n_boot = 200)
print(admin2_summary)

## ------------------------------------------------------------------
## 4.  Transport predictions to a new geographic area
##
## To apply the model to another area (e.g., a neighbouring country
## or a new survey), construct a new dataset with the same
## covariates.  Here we simulate a new area with different
## covariate distributions (e.g., older children, higher rural
## proportion).  The trained model is used to predict deficiency
## probabilities.  Aggregation and uncertainty estimation are
## performed as above.

# Simulate a new area with different covariate distributions
m <- 800  # number of individuals in the new area
admin2_new <- sample(paste0("new_admin2_", sprintf("%02d", 1:5)), size = m, replace = TRUE)
age_months_new   <- rpois(m, lambda = 48)       # younger children
sex_new          <- rbinom(m, size = 1, prob = 0.6)
rural_new        <- rbinom(m, size = 1, prob = 0.7)
maternal_edu_new <- sample(0:3, size = m, replace = TRUE)

new_df <- data.frame(
  age_months    = age_months_new,
  sex           = factor(sex_new, levels = c(0, 1), labels = c("female", "male")),
  rural         = factor(rural_new, levels = c(0, 1), labels = c("urban", "rural")),
  maternal_edu  = factor(maternal_edu_new),
  admin2        = factor(admin2_new)
)

# Predict deficiency probabilities in the new area
new_df$pred_prob <- predict(model_fit, newdata = new_df, type = "prob")[, "yes"]

# Aggregate predictions and estimate uncertainty in the new area
new_admin2_summary <- aggregate_predictions(new_df, n_boot = 200)
print(new_admin2_summary)

## ------------------------------------------------------------------
# End of tutorial
#
# This script generates a synthetic dataset, fits a predictive
# model for micronutrient deficiency, aggregates individual
# predictions to an administrative level, quantifies uncertainty
# using bootstrap methods, and transports predictions to a new
# geographic area.  Replace the simulated data with real survey or
# biomarker data and adjust the model specification and
# aggregation functions as required by your research question.