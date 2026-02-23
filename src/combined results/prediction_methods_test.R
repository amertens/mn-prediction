############################################################
# Interpretable evaluation of micronutrient prediction models
# Simulated example with survey weights
############################################################

rm(list = ls())

# ---- Packages ----
library(tidyverse)
library(yardstick)   # weighted performance metrics
library(survey)      # survey-weighted prevalence
library(sf)          # optional, for maps if desired
library(scales)
library(here)



load(here("data/IPD/Ghana/clean_mn_data_Ghana.rdata"))
#df_child_VitA, df_mom_VitA, df_child_iron, df_mom_iron, df_mom_b12, df_mom_folate,
head(df_mom_folate)
df_admin <- df_mom_folate %>% select(Admin1, Admin2, gw_PSUStrat_weight)
res_bin=readRDS(file = here("results/models/res_full_bin_GW_Ghana_SL_women_folate.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_women_folate.rds"))
names(res)
res$Y
res$yhat_full
dim(df_mom_folate)
length(res$Y)


dat_individ = data.frame(df_admin, true_def=as.numeric(unclass(res_bin$Y)),true_conc=as.numeric(unclass(res$Y)), pred_prob=res_bin$yhat_full, pred_conc=res$yhat_full) %>%
  rename(weight=gw_PSUStrat_weight)
head(dat_individ)
table(dat_individ$true_def)

range(dat_individ$pred_prob, na.rm = TRUE)
mean(dat_individ$pred_prob >= 0.56, na.rm = TRUE)

def_cutoff=0.56




#uncertainty estimates
set.seed(123)
B <- 1000

boot_nat <- replicate(B, {
  idx <- sample(seq_len(nrow(dat_individ)), replace = TRUE)
  with(dat_individ[idx, ],
       sum(weight * pred_prob) / sum(weight))
})

quantile(boot_nat, c(0.025, 0.5, 0.975))


#performance evaluation
dat <- dat_individ %>%
  group_by(Admin1, Admin2) %>%
  summarise(
    # population size represented
    total_weight = sum(weight, na.rm = TRUE),

    # observed prevalence (survey-weighted)
    true_def = sum(weight * true_def, na.rm = TRUE) / total_weight,

    # predicted prevalence (probability-based, recommended)
    pred_prev_prob = sum(weight * pred_prob, na.rm = TRUE) / total_weight,

    # screening-based prevalence (decision-specific)
    pred_def = sum(weight * (pred_prob >= def_cutoff), na.rm = TRUE) / total_weight,

    # counts for diagnostics
    n_indiv = n(),
    weight=total_weight,
    .groups = "drop"
  )



thr <- 0.56  # set policy threshold explicitly

dat <- dat_individ %>%
  transmute(
    Admin1,
    Admin2,
    weight = as.numeric(weight),
    pred_prob = as.numeric(pred_prob),

    # ensure true_def is strictly 0/1
    true_def = as.integer(true_def == 1),
    true_level = true_conc,
    pred_level = pred_conc,

    # class prediction derived from probabilities
    pred_def = as.integer(pred_prob >= thr)
  ) %>%
  mutate(
    # yardstick expects factors for class metrics
    true_def = factor(true_def, levels = c(0, 1)),
    pred_def = factor(pred_def, levels = c(0, 1))
  ) %>%
  # optional: drop missing values
  filter(!is.na(weight), !is.na(true_def), !is.na(pred_def))


# ---- 2. Individual-level classification performance (weighted) ----

class_metrics <- metric_set(
  sens, spec, ppv, npv, accuracy
)

indiv_perf <- class_metrics(
  data = dat,
  truth = true_def,
  estimate = pred_def,
  case_weights = weight
)

print(indiv_perf)

# ---- 3. Confusion matrix translated to counts (weighted) ----

confusion_weighted <- dat %>%
  mutate(
    cell = case_when(
      true_def == 1 & pred_def == 1 ~ "True positive",
      true_def == 1 & pred_def == 0 ~ "False negative",
      true_def == 0 & pred_def == 1 ~ "False positive",
      true_def == 0 & pred_def == 0 ~ "True negative"
    )
  ) %>%
  group_by(cell) %>%
  summarise(
    weighted_population = sum(weight),
    .groups = "drop"
  ) %>%
  mutate(
    pct = weighted_population / sum(weighted_population)
  )

print(confusion_weighted)

# ---- 4. Continuous outcome performance ----

cont_perf <- dat %>%
  summarise(
    weighted_MAE = weighted.mean(abs(pred_level - true_level), weight),
    weighted_bias = weighted.mean(pred_level - true_level, weight),
    correlation = cor(true_level, pred_level)
  )

print(cont_perf)

# ---- 5. Scatter plot: predicted vs true continuous values ----

ggplot(dat, aes(true_level, pred_level)) +
  geom_point(alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Measured biomarker level",
    y = "Predicted biomarker level",
    title = "Agreement between predicted and measured biomarker values",
    subtitle = "Dashed line indicates perfect agreement"
  ) +
  theme_minimal()

# ---- 6. Survey-weighted prevalence (national) ----

design <- svydesign(
  ids = ~1,
  weights = ~weight,
  data = dat
)

true_prev_nat <- svymean(~true_def, design)
pred_prev_nat <- svymean(~pred_def, design)

rbind(
  true = coef(true_prev_nat),
  predicted = coef(pred_prev_nat)
)

# ---- 7. Prevalence by Admin-1 (weighted) ----

prev_admin1 <- dat %>%
  group_by(admin1) %>%
  summarise(
    true_prev = sum(weight * as.numeric(true_def)-1, na.rm=T) / sum(weight, na.rm=T),
    pred_prev = sum(weight * as.numeric(pred_def)-1, na.rm=T) / sum(weight, na.rm=T),
    .groups = "drop"
  ) %>%
  mutate(
    abs_error = abs(pred_prev - true_prev)
  )

print(prev_admin1)

# ---- 8. Aggregate performance summaries (Admin-1) ----

agg_perf <- prev_admin1 %>%
  summarise(
    correlation = cor(true_prev, pred_prev),
    mean_abs_error = mean(abs_error),
    pct_within_5pp = mean(abs_error <= 0.05)
  )

print(agg_perf)

# ---- 9. Scatter plot: predicted vs observed prevalence ----

ggplot(prev_admin1, aes(true_prev, pred_prev, label = admin1)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(nudge_y = 0.01, size = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Observed deficiency prevalence",
    y = "Predicted deficiency prevalence",
    title = "Model agreement with survey prevalence by region",
    subtitle = "Each point represents an Admin-1 unit"
  ) +
  theme_minimal()

# ---- 10. Calibration plot (probability calibration) ----

calibration_df <- dat %>%
  mutate(
    risk_bin = ntile(pred_prob, 10)
  ) %>%
  group_by(risk_bin) %>%
  summarise(
    mean_pred = weighted.mean(as.numeric(true_def)-1, weight),
    obs_prev = weighted.mean(as.numeric(true_def)-1, weight),
    .groups = "drop"
  )

ggplot(calibration_df, aes(mean_pred, obs_prev)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Average predicted probability",
    y = "Observed deficiency prevalence",
    title = "Calibration of deficiency risk predictions",
    subtitle = "Points near the dashed line indicate good calibration"
  ) +
  theme_minimal()

############################################################
# End of script
############################################################




############################################################
# One-page model performance dashboard (slide-ready)
############################################################

library(tidyverse)
library(yardstick)
library(patchwork)
library(scales)

# ---- 1. Key weighted classification metrics (top-left panel) ----

metrics_tbl <- metric_set(sens, spec, ppv, npv)(
  data = dat,
  truth = true_def,
  estimate = pred_def,
  case_weights = weight
) %>%
  mutate(
    metric_label = recode(
      .metric,
      sens = "Sensitivity (deficient cases identified)",
      spec = "Specificity (non-deficient cleared)",
      ppv  = "Positive predictive value",
      npv  = "Negative predictive value"
    ),
    value = percent(.estimate, accuracy = 1)
  ) %>%
  select(metric_label, value)

metrics_plot <- ggplot(metrics_tbl, aes(metric_label, value)) +
  geom_text(aes(label = value), size = 6) +
  coord_flip() +
  labs(
    title = "Does the model find deficient individuals?",
    subtitle = "Survey-weighted individual-level performance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

# ---- 2. Calibration plot (top-right panel) ----

cal_df <- dat %>%
  mutate(risk_bin = ntile(pred_prob, 10)) %>%
  group_by(risk_bin) %>%
  summarise(
    mean_pred = weighted.mean(pred_prob, weight),
    obs_prev  = weighted.mean(as.numeric(as.character(true_def)), weight),
    .groups = "drop"
  )

cal_plot <- ggplot(cal_df, aes(mean_pred, obs_prev)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Are predicted risks well calibrated?",
    subtitle = "Observed vs predicted deficiency prevalence",
    x = "Average predicted risk",
    y = "Observed prevalence"
  ) +
  theme_minimal(base_size = 14)

# ---- 3. Admin-1 prevalence agreement (bottom-left panel) ----

Admin2_df <- dat %>%
  group_by(Admin2) %>%
  summarise(
    true_prev = sum(weight * as.numeric(as.character(true_def))) / sum(weight),
    pred_prev = sum(weight * pred_prob) / sum(weight),
    .groups = "drop"
  )

prev_plot <- ggplot(Admin2_df, aes(true_prev, pred_prev, label = Admin2)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(nudge_y = 0.015, size = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Does the model reproduce regional prevalence?",
    subtitle = "Admin-1 observed vs predicted prevalence",
    x = "Observed prevalence",
    y = "Predicted prevalence"
  ) +
  theme_minimal(base_size = 14)

# ---- 4. Continuous biomarker agreement (bottom-right panel) ----

bio_plot <- ggplot(dat, aes(true_level, pred_level)) +
  geom_point(alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "How close are biomarker predictions?",
    subtitle = "Predicted vs measured concentrations",
    x = "Measured biomarker",
    y = "Predicted biomarker"
  ) +
  theme_minimal(base_size = 14)

# ---- 5. Assemble dashboard ----

dashboard <- (metrics_plot | cal_plot) /
  (prev_plot    | bio_plot) +
  plot_annotation(
    title = "Model performance summary for policy and program use",
    subtitle = "Screening accuracy, calibration, geographic validity, and biological agreement",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13)
    )
  )

# ---- 6. Export (slide-ready) ----

ggsave(
  filename = "model_performance_dashboard.png",
  plot = dashboard,
  width = 14,
  height = 8,
  dpi = 300
)

############################################################
# End
############################################################
