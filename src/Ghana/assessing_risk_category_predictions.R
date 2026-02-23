


rm(list=ls())
library(dplyr)
library(tidyverse)
library(haven)
library(here)
library(purrr)
library(labelled)
library(tlverse)
library(sl3)
library(sf)
library(survey)
library(srvyr)     # tidy interface to 'survey'
library(pROC)

source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#read clean and merged data
df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))
df$gw_wVAD <- as.numeric(df$gw_wVAD)
df_child_VitA <- df %>% select(gw_cVAD, gw_cnum, gw_Strata, gw_sWeight) %>% as.data.frame() %>% filter(!is.na(df$gw_childid), !is.na(gw_cVAD))


adm2_prev = readRDS(here("data", "IPD", "Ghana", "Ghana_GMS_admin2_prevalence.rds"))

true_admin2_prev <- adm2_prev %>%
  filter(outcome == "Child VAD") %>%
  select(Admin2, prev, prev_upp ) %>%
  #mild  ≥2%–<10% moderate  ≥10%–<20%b severe  ≥20%
  mutate(true_cat = case_when(
    prev < 0.02 ~ "None",
    prev >= 0.02 & prev < 0.10 ~ "Mild",
    prev >= 0.10 & prev < 0.20 ~ "Moderate",
    prev >= 0.20 ~ "Severe",
    TRUE ~ NA_character_
  ),
  true_cat_ub = case_when(
    prev_upp  < 0.02 ~ "None",
    prev_upp  >= 0.02 & prev_upp  < 0.10 ~ "Mild",
    prev_upp  >= 0.10 & prev_upp  < 0.20 ~ "Moderate",
    prev_upp  >= 0.20 ~ "Severe",
    TRUE ~ NA_character_
  ))
table(true_admin2_prev$true_cat)
table(true_admin2_prev$true_cat_ub)

res1 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA_full.rds"))
res1 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA.rds"))


roc_obj <- pROC::roc(response = df_child_VitA$gw_cVAD, predictor = res1$yhat_full, quiet = T, direction = "<")
roc_auc <- as.numeric(pROC::auc(roc_obj))
roc_auc

best_threshold <- coords(roc_obj, "best", ret = "threshold")
print(paste("Optimal threshold (Youden):", best_threshold))



res1$Y
df_child_VitA$pred_Y <- ifelse(res1$yhat_full < as.numeric(best_threshold), 1, 0) #How to pick cutoff?? AUC prioritized predicting 0- pROC?
hist(res1$yhat_full)

mean(df_child_VitA$gw_cVAD)
mean(df_child_VitA$pred_Y)
table(df_child_VitA$gw_cVAD==1, df_child_VitA$pred_Y)
table(res1$Y==1, df_child_VitA$pred_Y)





gps = read.csv(file=here("data", "IPD", "Ghana", "Ghana_GMS_GPS_cleaned.csv"))
colnames(gps)[2] <- "gw_cnum"
d<-left_join(df_child_VitA,gps, by="gw_cnum")



summary(d$latitude)
summary(d$longitude)
d$lat= as.numeric(d$latitude)
d$lon= as.numeric(d$longitude)

poly.adm <- geodata::gadm(country="GH", level=2, path=tempdir())
poly.adm <- sf::st_as_sf(poly.adm) %>% select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
d_sf <- st_as_sf(d, coords = c("longitude","latitude"), crs = 4326)
poly.adm <- st_transform(poly.adm, crs = 4326)
#df <- as.data.frame(st_join(d_sf, poly.adm, join = st_within))
df <- (st_join(d_sf, poly.adm, join = st_within))


# Drop geometry
d_nog <- sf::st_drop_geometry(df)

# --- Set PSU and strata identifiers manually ---
# Replace with the correct names in your dataset
psu_var    <- "EANAME"   # e.g., cluster / EA / PSU column
strata_var <- "gw_Strata"      # e.g., survey strata column

# Make sure they're factors/numeric as needed
d_nog[[psu_var]]    <- as.factor(d_nog[[psu_var]])
d_nog[[strata_var]] <- as.factor(d_nog[[strata_var]])
d_nog$gw_sWeight       <- as.numeric(d_nog$gw_sWeight)

#set options for survey design
options(survey.lonely.psu = "adjust")

# --- Helper: compute admin-2 prevalence for outcomes ---
compute_adm2_prev <- function(data, outcomes) {
  dat <- data[rowSums(sapply(outcomes, function(v) !is.na(data[[v]]))) > 0, , drop = FALSE]

  des <- srvyr::as_survey_design(
    dat,
    ids     = !!rlang::sym(psu_var),
    strata  = !!rlang::sym(strata_var),
    weights = !!rlang::sym("gw_sWeight"),
    nest    = TRUE
  )

  purrr::map_dfr(outcomes, function(v) {
    des %>%
      dplyr::group_by(Admin2) %>%
      dplyr::summarise(prev = survey_mean(!!rlang::sym(v), vartype = c("se","ci"), na.rm = TRUE)) %>%
      dplyr::mutate(outcome = v)
  })
}

# predict outcomes

N_per_cluster =d_nog%>% group_by(Admin2) %>% summarise(N=n())

pred_adm2_prev <- compute_adm2_prev(d_nog, "pred_Y")
colnames(pred_adm2_prev) <- c("Admin2","pred_prev",  "pred_prev_se", "pred_prev_low", "pred_prev_upp", "outcome")
pred_adm2_prev <- pred_adm2_prev %>%
  select(Admin2, pred_prev, pred_prev_upp ) %>%
  #mild  ≥2%–<10% moderate  ≥10%–<20%b severe  ≥20%
  mutate(pred_cat = case_when(
    pred_prev < 0.02 ~ "None",
    pred_prev >= 0.02 & pred_prev < 0.10 ~ "Mild",
    pred_prev >= 0.10 & pred_prev < 0.20 ~ "Moderate",
    pred_prev >= 0.20 ~ "Severe",
    TRUE ~ NA_character_
  ),
  pred_cat_ub = case_when(
    pred_prev_upp  < 0.02 ~ "None",
    pred_prev_upp  >= 0.02 & pred_prev_upp  < 0.10 ~ "Mild",
    pred_prev_upp  >= 0.10 & pred_prev_upp  < 0.20 ~ "Moderate",
    pred_prev_upp  >= 0.20 ~ "Severe",
    TRUE ~ NA_character_
  ))

table(true_admin2_prev$true_cat )
table(pred_adm2_prev$pred_cat )

res <- left_join(true_admin2_prev, pred_adm2_prev, by="Admin2")
res <- left_join(res, N_per_cluster, by="Admin2")


head(res)

#drop admn areas based on limited data
summary(res$N)
res %>% filter(true_cat=="Severe" & pred_cat=="None")

# Required packages
library(caret)
library(weights)
library(psych)

# Convert to ordered factors (if not already)
true_categories <- factor(res$true_cat, levels = c("None", "Mild", "Moderate", "Severe"), ordered = TRUE)
pred_categories <- factor(res$pred_cat, levels = c("None", "Mild", "Moderate", "Severe"), ordered = TRUE)

levels(true_categories)
levels(pred_categories)

table(true_categories, pred_categories)

# Convert to numeric for ordinal calculations (1, 2, 3, 4)
true_numeric <- as.numeric(true_categories)
pred_numeric <- as.numeric(pred_categories)

print("Data summary:")
print(paste("Levels:", paste(levels(true_categories), collapse = ", ")))
print(paste("Sample size:", length(true_categories)))

# 1. CORE ORDINAL METRICS
#=====================================

# Mean Absolute Error (MAE) - primary metric for ordinal data
mae <- mean(abs(true_numeric - pred_numeric))
print(paste("Mean Absolute Error (MAE):", round(mae, 4)))

# Root Mean Squared Error (RMSE)
rmse <- sqrt(mean((true_numeric - pred_numeric)^2))
print(paste("Root Mean Squared Error (RMSE):", round(rmse, 4)))

# Mean Squared Error (MSE)
mse <- mean((true_numeric - pred_numeric)^2)
print(paste("Mean Squared Error (MSE):", round(mse, 4)))

# 2. CORRELATION MEASURES
#=====================================

# Spearman's rank correlation (primary for ordinal)
spearman_cor <- cor(true_numeric, pred_numeric, method = "spearman")
print(paste("Spearman Correlation:", round(spearman_cor, 4)))

# Kendall's tau (alternative rank correlation)
kendall_tau <- cor(true_numeric, pred_numeric, method = "kendall")
print(paste("Kendall's Tau:", round(kendall_tau, 4)))

# Pearson correlation (for comparison)
pearson_cor <- cor(true_numeric, pred_numeric, method = "pearson")
print(paste("Pearson Correlation:", round(pearson_cor, 4)))

# 3. EXACT AND ADJACENT ACCURACY
#=====================================

# Exact match accuracy
exact_accuracy <- mean(true_numeric == pred_numeric)
print(paste("Exact Accuracy:", round(exact_accuracy, 4)))

# Adjacent accuracy (within 1 category)
adjacent_accuracy <- mean(abs(true_numeric - pred_numeric) <= 1)
print(paste("Adjacent Accuracy (±1 level):", round(adjacent_accuracy, 4)))

# Within 2 categories accuracy (for 4-level outcome)
within_2_accuracy <- mean(abs(true_numeric - pred_numeric) <= 2)
print(paste("Within 2 Levels Accuracy:", round(within_2_accuracy, 4)))

# 4. ORDINAL-SPECIFIC PERFORMANCE METRICS
#=====================================

# Directional accuracy (same direction from middle)
# Assuming levels 1,2 are "lower" and 3,4 are "higher"
middle_point <- (max(true_numeric) + min(true_numeric)) / 2

true_direction <- ifelse(true_numeric > middle_point, 1, 0)
pred_direction <- ifelse(pred_numeric > middle_point, 1, 0)
directional_accuracy <- mean(true_direction == pred_direction)
print(paste("Directional Accuracy:", round(directional_accuracy, 4)))

# 5. WEIGHTED METRICS (heavier penalty for larger errors)
#=====================================

# Weighted MAE (quadratic penalty for distance)
weights_quadratic <- (abs(true_numeric - pred_numeric))^2
weighted_mae_quad <- mean(weights_quadratic)
print(paste("Quadratically Weighted MAE:", round(weighted_mae_quad, 4)))

# Linear weighted MAE (already calculated above as regular MAE)
print(paste("Linearly Weighted MAE:", round(mae, 4)))

# 6. CONFUSION MATRIX FOR ORDINAL DATA
#=====================================

# Standard confusion matrix
conf_matrix <- table(Predicted = pred_categories, Actual = true_categories)
print("Confusion Matrix:")
print(conf_matrix)

# Calculate off-diagonal penalties
n_levels <- length(levels(true_categories))
penalty_matrix <- abs(outer(1:n_levels, 1:n_levels, "-"))

# Weighted accuracy using distance penalties
total_penalty <- sum(conf_matrix * penalty_matrix)
max_possible_penalty <- sum(conf_matrix) * max(penalty_matrix)
ordinal_accuracy <- 1 - (total_penalty / max_possible_penalty)
print(paste("Ordinal Accuracy (distance-weighted):", round(ordinal_accuracy, 4)))

# 7. DISTRIBUTION COMPARISON
#=====================================

# Compare distributions
true_dist <- table(true_categories) / length(true_categories)
pred_dist <- table(pred_categories) / length(pred_categories)

print("Distribution Comparison:")
comparison_df <- data.frame(
  Level = names(true_dist),
  True_Prop = as.numeric(true_dist),
  Pred_Prop = as.numeric(pred_dist),
  Difference = as.numeric(pred_dist - true_dist)
)
print(round(comparison_df, 4))

# Kolmogorov-Smirnov test for distribution difference
ks_test <- ks.test(true_numeric, pred_numeric)
print(paste("KS Test p-value:", round(ks_test$p.value, 4)))

# 8. ERROR ANALYSIS
#=====================================

# Error distribution
errors <- pred_numeric - true_numeric
error_dist <- table(errors)
print("Error Distribution:")
print(error_dist)

# Mean error (bias)
mean_error <- mean(errors)
print(paste("Mean Error (bias):", round(mean_error, 4)))

# Error by true category
error_by_category <- aggregate(errors, by = list(true_categories), FUN = mean)
names(error_by_category) <- c("True_Category", "Mean_Error")
print("Mean Error by True Category:")
print(round(error_by_category, 4))

# 9. VISUALIZATION
#=====================================

# Install if needed: install.packages(c("ggplot2", "corrplot"))
library(ggplot2)

# Error distribution plot
error_df <- data.frame(
  True = true_numeric,
  Predicted = pred_numeric,
  Error = errors
)

# Scatter plot with perfect prediction line
ggplot(error_df, aes(x = True, y = Predicted)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(breaks = 1:4) +
  labs(title = "Predicted vs True Categories",
       x = "True Category",
       y = "Predicted Category") +
  theme_minimal()

# 10. COMPREHENSIVE SUMMARY FUNCTION
#=====================================

ordinal_performance_summary <- function(true_ord, pred_ord) {
  # Convert to numeric
  true_num <- as.numeric(as.ordered(true_ord))
  pred_num <- as.numeric(as.ordered(pred_ord))

  # Calculate all metrics
  results <- list(
    mae = mean(abs(true_num - pred_num)),
    rmse = sqrt(mean((true_num - pred_num)^2)),
    spearman = cor(true_num, pred_num, method = "spearman"),
    kendall = cor(true_num, pred_num, method = "kendall"),
    exact_accuracy = mean(true_num == pred_num),
    adjacent_accuracy = mean(abs(true_num - pred_num) <= 1),
    mean_error = mean(pred_num - true_num)
  )

  return(results)
}

# Use the summary function
summary_results <- ordinal_performance_summary(df_child_VitA$gw_cVAD,
                                               df_child_VitA$pred_Y)

print("\n=== ORDINAL PERFORMANCE SUMMARY ===")
for(i in 1:length(summary_results)) {
  print(paste(names(summary_results)[i], ":",
              round(summary_results[[i]], 4)))
}
