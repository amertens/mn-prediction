###############################################################################
# Admin2 prevalence modeling with design-based outcomes + Spatial SuperLearner
# Rows for ML: admin2 units (not individuals)
# Training outcomes: design-based admin2 prevalence estimates from biomarker surveys
# Prediction targets: admin2 units in countries without biomarker data
###############################################################################

library(dplyr)
library(tidyr)
library(survey)
library(SuperLearner)

set.seed(1)
options(survey.lonely.psu = "adjust")

###############################################################################
# 0) USER-SUPPLIED INPUTS (you must adapt variable names)
###############################################################################
# (A) Individual-level biomarker survey microdata (stacked across countries)
# Required columns:
#   country      : country code/name
#   admin2_id    : admin2 identifier (consistent with admin2_covars)
#   y_bin        : deficiency indicator (0/1) defined from biomarkers
#   w_svy        : survey weight
#   psu          : PSU/cluster id
#   strata       : strata id
#
# df_ind <- readRDS("biomarker_microdata_4countries.rds")

# (B) Admin2 covariates for ALL admin2 in ALL countries you want to predict
# Required columns:
#   country, admin2_id
#   lon, lat       : admin2 centroid (for spatial blocks)
#   admin1_id      : admin1 id (useful for spatial CV)
#   X_* predictors : aggregated predictors at admin2 level
#

set.seed(123)

library(dplyr)
library(tidyr)

###############################################################################
# 1) GLOBAL STRUCTURE
###############################################################################

countries <- c("A", "B", "C", "D")
n_admin1_per_country <- 5
n_admin2_per_admin1  <- 6
n_psu_per_admin2     <- 4
n_ind_per_psu        <- 25

###############################################################################
# 2) SIMULATE ADMIN STRUCTURE + COVARIATES
###############################################################################

admin2_covars <- expand.grid(
  country = countries,
  admin1  = paste0("adm1_", 1:n_admin1_per_country),
  admin2  = paste0("adm2_", 1:n_admin2_per_admin1)
) %>%
  mutate(
    admin1_id = paste(country, admin1, sep = "_"),
    admin2_id = paste(country, admin1, admin2, sep = "_")
  ) %>%
  select(country, admin1_id, admin2_id)

# Assign spatial coordinates (clustered by country)
admin2_covars <- admin2_covars %>%
  group_by(country) %>%
  mutate(
    lon = rnorm(n(), mean = match(country, countries) * 5, sd = 1),
    lat = rnorm(n(), mean = match(country, countries) * 3, sd = 1)
  ) %>%
  ungroup()

# Generate admin2-level predictors
admin2_covars <- admin2_covars %>%
  mutate(
    X_ndvi = runif(n(), 0, 1),
    X_poverty = runif(n(), 0, 1),
    X_malaria = runif(n(), 0, 1),
    X_popdens = runif(n(), 0, 1)
  )

###############################################################################
# 3) TRUE PREVALENCE MODEL (ADMIN2 LEVEL)
###############################################################################

# Spatial + covariate-driven true prevalence
admin2_covars <- admin2_covars %>%
  mutate(
    linear_pred =
      -1.5 +
      1.2 * X_poverty +
      0.8 * X_malaria -
      0.7 * X_ndvi +
      0.3 * X_popdens +
      0.2 * scale(lon)[,1] +
      0.2 * scale(lat)[,1] +
      rnorm(n(), 0, 0.3),
    true_prev = plogis(linear_pred)
  )

###############################################################################
# 4) SIMULATE INDIVIDUAL-LEVEL SURVEY DATA
###############################################################################

df_ind <- admin2_covars %>%
  select(country, admin1_id, admin2_id, true_prev) %>%
  slice(rep(1:n(), each = n_psu_per_admin2)) %>%
  group_by(admin2_id) %>%
  mutate(psu = paste0(admin2_id, "_psu_", row_number())) %>%
  ungroup() %>%
  slice(rep(1:n(), each = n_ind_per_psu)) %>%
  group_by(psu) %>%
  mutate(
    strata = paste0(country, "_strata_", sample(1:3, 1)),
    # Slight PSU-level random effect
    psu_effect = rnorm(1, 0, 0.2)
  ) %>%
  ungroup()

# Generate individual binary biomarker outcome
df_ind <- df_ind %>%
  mutate(
    p_ind = plogis(qlogis(true_prev) + psu_effect),
    y_bin = rbinom(n(), 1, p_ind)
  )

###############################################################################
# 5) SIMULATE COMPLEX SURVEY WEIGHTS
###############################################################################

# Unequal probability sampling by poverty (oversample high poverty areas)
df_ind <- df_ind %>%
  left_join(admin2_covars %>% select(admin2_id, X_poverty), by = "admin2_id") %>%
  mutate(
    inclusion_prob = 0.02 + 0.05 * X_poverty,  # oversampling poorer areas
    w_svy = 1 / inclusion_prob,
    w_svy = w_svy * runif(n(), 0.8, 1.2)  # noise
  ) %>%
  select(-X_poverty)

###############################################################################
# 6) FINAL DATASETS
###############################################################################

# df_ind is now ready for survey-based direct admin2 prevalence estimation
print(head(df_ind))

# admin2_covars includes predictors and true_prev (for validation only)
print(head(admin2_covars))



###############################################################################
# 1) FUNCTION: DESIGN-BASED ADMIN2 PREVALENCE FROM MICRODATA
###############################################################################
direct_admin2_from_survey <- function(df_country,
                                      y = "y_bin",
                                      w = "w_svy",
                                      psu = "psu",
                                      strata = "strata",
                                      admin2 = "admin2_id") {

  des <- svydesign(
    ids = as.formula(paste0("~", psu)),
    strata = as.formula(paste0("~", strata)),
    weights = as.formula(paste0("~", w)),
    data = df_country,
    nest = TRUE
  )

  est <- svyby(
    formula = as.formula(paste0("~", y)),
    by = as.formula(paste0("~", admin2)),
    design = des,
    FUN = svymean,
    vartype = c("se"),
    na.rm = TRUE
  ) %>% as_tibble()

  out <- est %>%
    transmute(
      admin2_id = .data[[admin2]],
      y_hat = .data[[y]],
      se = se
    ) %>%
    filter(!is.na(y_hat), !is.na(se), se > 0)

  out
}

###############################################################################
# 2) BUILD TRAINING DATASET: ADMIN2 OUTCOMES + ADMIN2 COVARIATES
###############################################################################
# Direct estimates per country, then stack
direct_admin2 <- df_ind %>%
  group_by(country) %>%
  group_modify(~ direct_admin2_from_survey(.x)) %>%
  ungroup()

# Merge direct outcomes with admin2 covariates (admin2_covars must include country+admin2_id)
train_admin2 <- admin2_covars %>%
  inner_join(
    direct_admin2 %>% mutate(country = NULL),  # admin2_id should be globally unique, or keep country key
    by = "admin2_id"
  )

# If admin2_id is not globally unique, join by (country, admin2_id) instead:
# train_admin2 <- admin2_covars %>%
#   inner_join(direct_admin2, by = c("country","admin2_id"))

# Define ML inputs
Y <- train_admin2$y_hat
X <- train_admin2 %>% select(starts_with("X_")) %>% as.data.frame()

# Precision weights: give more influence to admin2 with more precise design-based estimates
w_admin2 <- 1 / (train_admin2$se^2)
w_admin2 <- pmin(w_admin2, quantile(w_admin2, 0.99, na.rm = TRUE))  # cap extremes

###############################################################################
# 3) SPATIAL / TRANSPORT CV FOLDS
###############################################################################
# You should evaluate extrapolation. Two recommended schemes:
#   (i) Leave-one-country-out (LOCO): strict transportability
#   (ii) Leave-one-admin1-out (LOAO) within countries: geographic extrapolation
#
# SuperLearner wants validRows = list of indices held out for each fold.

make_folds_by_group <- function(group_vec) {
  u <- unique(group_vec)
  lapply(u, function(g) which(group_vec == g))
}

# (i) Leave-one-country-out folds (best for “countries without biomarker data” claim)
validRows_loco <- make_folds_by_group(train_admin2$country)

# (ii) Leave-one-admin1-out folds (within-country spatial transport)
# validRows_loao <- make_folds_by_group(paste(train_admin2$country, train_admin2$admin1_id, sep="__"))

# Choose one CV scheme for fitting / tuning
cv_ctrl <- list(V = length(validRows_loco), validRows = validRows_loco)

###############################################################################
# 4) FIT SURVEY-WEIGHTED SUPERLEARNER (ADMIN2-LEVEL OUTCOME)
###############################################################################
# Outcome is continuous prevalence estimate -> gaussian family.
# Use a modest library; expand later.
SL_lib <- c(
  "SL.glm",
  "SL.glmnet",
  "SL.gam",
  "SL.randomForest"
  # add "SL.xgboost" or "SL.ranger" if available and stable in your environment
)

sl_fit <- SuperLearner(
  Y = Y,
  X = X,
  family = gaussian(),
  SL.library = SL_lib,
  method = "method.NNLS",
  obsWeights = w_admin2,
  cvControl = cv_ctrl
)

print(sl_fit)

###############################################################################
# 5) PREDICT TO ALL ADMIN2 IN ALL COUNTRIES (INCLUDING NO-BIOMARKER COUNTRIES)
###############################################################################
newX_all <- admin2_covars %>% select(starts_with("X_")) %>% as.data.frame()
pred_all <- predict(sl_fit, newdata = newX_all)$pred
pred_all <- pmin(pmax(pred_all, 0), 1)

admin2_pred <- admin2_covars %>%
  transmute(country, admin2_id,
            pred_prev = pred_all)

###############################################################################
# 6) OPTIONAL: REPORT TRANSPORTABILITY PERFORMANCE FROM LOCO CV
###############################################################################
# SuperLearner stores cross-validated risks; to get per-fold diagnostics,
# do explicit LOCO evaluation (recommended for papers and talks).

loco_eval <- function(train_df, SL_lib, cap_quantile = 0.99) {

  countries <- unique(train_df$country)
  res <- vector("list", length(countries))

  for (k in seq_along(countries)) {
    c_hold <- countries[k]
    tr <- train_df %>% filter(country != c_hold)
    te <- train_df %>% filter(country == c_hold)

    Ytr <- tr$y_hat
    Xtr <- tr %>% select(starts_with("X_")) %>% as.data.frame()
    wtr <- 1 / (tr$se^2)
    wtr <- pmin(wtr, quantile(wtr, cap_quantile, na.rm = TRUE))

    fit <- SuperLearner(
      Y = Ytr,
      X = Xtr,
      family = gaussian(),
      SL.library = SL_lib,
      method = "method.NNLS",
      obsWeights = wtr
    )

    Xte <- te %>% select(starts_with("X_")) %>% as.data.frame()
    pred <- pmin(pmax(predict(fit, newdata = Xte)$pred, 0), 1)

    res[[k]] <- tibble(
      heldout_country = c_hold,
      n_admin2 = nrow(te),
      rmse = sqrt(mean((pred - te$y_hat)^2, na.rm = TRUE)),
      mae  = mean(abs(pred - te$y_hat), na.rm = TRUE),
      corr = suppressWarnings(cor(pred, te$y_hat, use = "complete.obs"))
    )
  }

  bind_rows(res)
}

loco_summary <- loco_eval(train_admin2, SL_lib)
print(loco_summary)

###############################################################################
# OUTPUT
###############################################################################
# admin2_pred: predictions for every admin2 in admin2_covars (including countries with no biomarkers)
# loco_summary: strict transportability diagnostics
###############################################################################
