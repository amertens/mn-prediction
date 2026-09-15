


rm(list = ls())

# ---- Packages ----
library(tidyverse)
library(yardstick)   # weighted performance metrics
library(survey)      # survey-weighted prevalence
library(scales)
library(here)
library(pROC)

source(paste0(here::here(),"/src/0-functions.R"))


full_res=readRDS(here("results/compiled_predictions.RDS"))
head(full_res)


metrics = full_res %>% distinct(Country, outcome, population, n, prevalence,roc_auc, roc_ci_low, roc_ci_high,pr_auc, pr_baseline,  pr_gain, pr_norm_auc, pr_ci_low, pr_ci_high,     brier, brier_ref,
                               brier_skill, brier_ci_low, brier_ci_high, brier_skill_ci_low, brier_skill_ci_high) %>% arrange(roc_auc)

summary(metrics$roc_auc)
summary(metrics$pr_auc)
summary(metrics$pr_gain)
summary(metrics$brier_skill )

metrics[metrics$roc_auc>0.6,] %>% arrange(-pr_gain)

res=full_res %>% distinct(Country, outcome, population, null_model, sl_model, null_model_bin, sl_model_bin )
res$R2 <- (1-res$sl_model/res$null_model) *100
res$BSS <- (1-res$sl_model_bin/res$null_model_bin) *100

summary(res$BSS)
summary(res$R2)
#
# head(res)

full_res=readRDS(here("results/compiled_predictions.RDS"))
head(full_res)

full_res %>% group_by(Country, outcome, population) %>%
  summarise(
    brier = mean((yhat_full_bin - Y_bin)^2),
    brier_ref = mean(Y_bin) * (1 - mean(Y_bin))
    ) %>%
  mutate(brier_skill = (1 - (brier / brier_ref))*100  )




