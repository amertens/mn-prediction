
rm(list = ls())

# ---- Packages ----
library(tidyverse)
library(yardstick)   # weighted performance metrics
library(survey)      # survey-weighted prevalence
library(scales)
library(here)



full_res=readRDS(here("results/compiled_predictions.RDS"))
head(full_res)
table(full_res$Y_bin)

full_res %>% filter((Country=="Ghana"  & population=="children")) %>% distinct(outcome )

#full_res <- full_res %>% filter(!(Country=="Ghana" & outcome=="gw_wRBP" & population=="children"))

head(full_res)
full_res <- full_res %>% rename(
  true_conc=Y,
  pred_conc=yhat_full,
  true_def=Y_bin,
  pred_prob=yhat_full_bin
)


full_res %>% distinct(Country, outcome, population)

national_prev <- full_res %>% group_by(Country, outcome, population) %>%
  summarise(
    total_weight = sum(weight, na.rm = TRUE),
    n_cases=sum(true_def, na.rm = TRUE),
    n_pred_cases=round(sum(pred_prob, na.rm = TRUE),0),
    true_prev = sum(weight * true_def, na.rm = TRUE) / total_weight *100,
    pred_prev = sum(weight * pred_prob, na.rm = TRUE) / total_weight *100
  )  #%>%
  # mutate(
  #   who_category = cut(
  #     pred_prev * 100,
  #     breaks = c(-Inf, 5, 20, 40, Inf),
  #     labels = c(
  #       "Not a public health problem",
  #       "Mild",
  #       "Moderate",
  #       "Severe"
  #     ),
  #     right = FALSE
  #   )
  # )

national_prev


#uncertainty estimates
boot_prev_ci <- function(
    data,
    pred_prob,
    weight,
    B = 1000,
    conf = 0.95,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  pred_prob <- dplyr::pull(data, {{ pred_prob }})
  weight    <- dplyr::pull(data, {{ weight }})

  n <- length(pred_prob)

  boot_vals <- replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    sum(weight[idx] * pred_prob[idx], na.rm = TRUE) /
      sum(weight[idx], na.rm = TRUE)
  })

  alpha <- (1 - conf) / 2

  tibble::tibble(
    prev_lwr = quantile(boot_vals, alpha, na.rm = TRUE) *100,
    prev_upr = quantile(boot_vals, 1 - alpha, na.rm = TRUE) *100
  )
}

library(dplyr)

CIs <- full_res %>%
  group_by(Country, outcome, population) %>%  summarise(
    boot_prev_ci(
      data = cur_data(),
      pred_prob = pred_prob,
      weight = weight,
      B = 1000,
      seed = 123
    ),
    .groups = "drop"
  )


national_prev <- national_prev %>% left_join(CIs, by=c("Country", "outcome", "population"))
head(national_prev)

head(national_prev)

unique(national_prev$outcome)


national_prev <- national_prev %>%
  mutate(
    outcome_f = case_when(
      outcome == "gw_LogFerAdj"        ~ "Iron deficiency (Ferritin, BRINDA)",
      outcome == "gw_cFerrAdj"         ~ "Iron deficiency (Ferritin, BRINDA)",
      outcome == "gw_wFerrAdjThurn"    ~ "Iron deficiency (Ferritin, Thurnham)",
      outcome == "fer"                 ~ "Iron deficiency (Ferritin, unadjusted)",

      outcome == "gw_cRBPAdjThurn"     ~ "Vitamin A deficiency (RBP, Thurnham)",
      outcome == "gw_wRBPAdjThurn"     ~ "Vitamin A deficiency (RBP, Thurnham)",
      outcome == "gw_cRBPAdjBrinda"    ~ "Vitamin A deficiency (RBP, BRINDA)",
      outcome == "gw_wRBPAdjBR1"       ~ "Vitamin A deficiency (RBP, BRINDA)",
      outcome == "gw_wRBP"             ~ "Vitamin A deficiency (RBP, unadjusted)",
      outcome == "rbp"                 ~ "Vitamin A deficiency (RBP, unadjusted)",

      outcome == "gw_r_crpagp_sf2"     ~ "Vitamin A deficiency (RBP, CRP/AGP adjusted)",

      outcome == "gw_wFolate"          ~ "Folate deficiency",
      outcome == "gw_wB12"             ~ "Vitamin B12 deficiency",
      outcome == "zn_gdl"              ~ "Zinc deficiency",
      TRUE                             ~ outcome
    ),
    outcome_simple = case_when(
      outcome == "gw_LogFerAdj"        ~ "Iron deficiency",
      outcome == "gw_cFerrAdj"         ~ "Iron deficiency",
      outcome == "gw_wFerrAdjThurn"    ~ "Iron deficiency",
      outcome == "fer"                 ~ "Iron deficiency",
      outcome == "gw_r_crpagp_sf2"                 ~ "Iron deficiency",

      outcome == "gw_cRBPAdjThurn"     ~ "Vitamin A deficiency",
      outcome == "gw_wRBPAdjThurn"     ~ "Vitamin A deficiency",
      outcome == "gw_cRBPAdjBrinda"    ~ "Vitamin A deficiency",
      outcome == "gw_wRBPAdjBR1"       ~ "Vitamin A deficiency",
      outcome == "gw_wRBP"             ~ "Vitamin A deficiency",
      outcome == "rbp"                 ~ "Vitamin A deficiency",

      outcome == "gw_r_crpagp_sf2"     ~ "Vitamin A deficiency",

      outcome == "gw_wFolate"          ~ "Folate deficiency",
      outcome == "gw_wB12"             ~ "Vitamin B12 deficiency",
      outcome == "zn_gdl"              ~ "Zinc deficiency",
      TRUE                             ~ outcome
    )
  )


plot_df <- national_prev %>%
  mutate(
    true_prev = true_prev / 100,   # convert to proportions if currently %
    pred_prev = pred_prev / 100,
    prev_lwr  = prev_lwr  / 100,
    prev_upr  = prev_upr  / 100
  ) %>%
  select(
    Country, outcome_simple, population,
    true_prev, pred_prev, prev_lwr, prev_upr
  ) %>%
  pivot_longer(
    cols = c(true_prev, pred_prev),
    names_to = "type",
    values_to = "prevalence"
  ) %>%
  mutate(
    type = recode(
      type,
      true_prev = "Observed (survey)",
      pred_prev = "Predicted (model)"
    )
  )

plot_df$prev_lwr[plot_df$type=="Observed (survey)"] <- NA
plot_df$prev_upr[plot_df$type=="Observed (survey)"] <- NA




# tab:blue : #1f77b4
#   tab:orange : #ff7f0e
#   tab:green : #2ca02c
#   tab:red : #d62728
#   tab:purple : #9467bd
#   tab:brown : #8c564b
#   tab:pink : #e377c2
#   tab:gray : #7f7f7f
#   tab:olive : #bcbd22
#   tab:cyan : #17becf

plot_df$population <- str_to_title(plot_df$population)
plot_df$type=factor(plot_df$type, levels=c( "Predicted (model)","Observed (survey)"))
plot_df$Xlab=paste0(plot_df$Country," ", plot_df$outcome_simple)
plot_df$Xlab=factor(plot_df$Xlab, levels=rev(unique(plot_df$Xlab)))
ggplot(plot_df, aes(x = paste0(Country," ", outcome_simple), y = prevalence, color=type, group = type)) +
  geom_point(
    size = 2,
    position=position_dodge(0.5)
  ) +
  geom_errorbar(
    #data = subset(plot_df, type == "Predicted (model)"),
    aes(ymin = prev_lwr, ymax = prev_upr),
    width = 0.15,
    linewidth = 2,
    position=position_dodge(0.5)
  ) +
  facet_wrap(
     ~ population,
    scales = "free"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("#1f77b4", "#ff7f0e"),
    name = NULL
  ) +
  labs(
    title = "Observed and model-predicted national micronutrient deficiency prevalence",
    subtitle = "Points show prevalence; bars show bootstrapped 95% confidence intervals for model predictions",
    x = NULL,
    y = "Prevalence"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

