
rm(list = ls())

# ---- Packages ----
library(tidyverse)
library(yardstick)   # weighted performance metrics
library(survey)      # survey-weighted prevalence
library(haven)
library(here)

source(paste0(here::here(),"/src/0-functions.R"))

full_res_ghana=full_res_gambia=full_res_SL=full_res_malawi=NULL

df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))

ghana_df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))  %>%
  select(dataid, Admin1, Admin2, gw_PSUStrat_weight) %>% rename(weight=gw_PSUStrat_weight)
gambia_df <- readRDS(here("data/IPD/Gambia/Gambia_merged_dataset.rds"))  %>% select(dataid, Admin1, Admin2, gw_svy_weight) %>% rename(weight=gw_svy_weight)
SL_df <- readRDS(here("data/IPD/Sierra Leone/SierraLeone_merged_dataset.rds")) %>% select(dataid, Admin1, Admin2, gw_svy_weight) %>% rename(weight=gw_svy_weight)
malawi_df <- readRDS(here("data/Malawi/Malawi_merged_dataset.rds"))  %>% select(dataid, Admin1, Admin2, svy_weight) %>% rename(weight=svy_weight)


#-------------------------------------------------------------------------------
# Ghana
#-------------------------------------------------------------------------------

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_child_vitA_V2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_child_vitA_v2.rds"))


  res_df = left_join(res$res,
                     res_bin$res %>%
                       rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                     by=c("dataid", "clusterid", "population"))
  res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
  res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
  res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
  res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
  res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
  auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
 res_df= bind_cols(res_df, auc_pr_brier_res)
  full_res_ghana=bind_rows(full_res_ghana, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_women_vitA_v2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_women_vitA_v2.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_ghana=bind_rows(full_res_ghana, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_women_b12_v2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_women_b12_v2.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_ghana=bind_rows(full_res_ghana, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_women_folate_v2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_women_folate_v2.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_ghana=bind_rows(full_res_ghana, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_child_iron_v2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_child_iron_v2.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_ghana=bind_rows(full_res_ghana, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Ghana_SL_mom_iron_v2.rds"))
res=readRDS(file = here("results/models/res_GW_Ghana_SL_mom_iron_v2.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_ghana=bind_rows(full_res_ghana, res_df)

full_res_ghana <- left_join(full_res_ghana, ghana_df, by="dataid")


#-------------------------------------------------------------------------------
# Gambia
#-------------------------------------------------------------------------------

res_bin=readRDS(file = here("results/models/res_bin_GW_Gambia_SL_child_iron.rds"))
res=readRDS(file = here("results/models/res_GW_Gambia_SL_child_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_gambia=bind_rows(full_res_gambia, res_df)


res_bin=readRDS(file = here("results/models/res_bin_GW_Gambia_SL_mom_iron.rds"))
res=readRDS(file = here("results/models/res_GW_Gambia_SL_mom_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_gambia=bind_rows(full_res_gambia, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Gambia_SL_child_vitA.rds"))
res=readRDS(file = here("results/models/res_GW_Gambia_SL_child_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_gambia=bind_rows(full_res_gambia, res_df)

res_bin=readRDS(file = here("results/models/res_bin_GW_Gambia_SL_women_vitA.rds"))
res=readRDS(file = here("results/models/res_GW_Gambia_SL_women_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_gambia=bind_rows(full_res_gambia, res_df)

full_res_gambia <- left_join(full_res_gambia, gambia_df, by="dataid")
table(full_res_gambia$Y_bin)

#-------------------------------------------------------------------------------
# Sierra Leone
#-------------------------------------------------------------------------------

#vitA
res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_child_vitA_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_child_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)

res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_women_vitA_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_women_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)

#iron
res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_child_iron_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_child_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)


res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_iron_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)


#iodine
res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_iodine_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)

# #b12 too sparse
# res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_b12_bin.rds"))
# res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_b12.rds"))
# res_df = left_join(res$res,
#                    res_bin$res %>%
#                      rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
#                    by=c("dataid", "clusterid", "population"))
# full_res_SL=bind_rows(full_res_SL, res_df)


#folate
res_bin=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_folate_bin.rds"))
res=readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_mom_folate.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_SL=bind_rows(full_res_SL, res_df)
table(full_res_SL$Y_bin)

full_res_SL <- left_join(full_res_SL, SL_df, by="dataid")


#-------------------------------------------------------------------------------
# Malawi
#-------------------------------------------------------------------------------

res_bin=readRDS(file = here("results/models/res_Malawi_SL_child_vitA_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_child_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)

res_bin=readRDS(file = here("results/models/res_Malawi_SL_women_vitA_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_women_vitA.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)


res_bin=readRDS(file = here("results/models/res_Malawi_SL_child_iron_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_child_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)


res_bin=readRDS(file = here("results/models/res_Malawi_SL_mom_iron_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_mom_iron.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)

res_bin=readRDS(file = here("results/models/res_Malawi_SL_child_zinc_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_child_zinc.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)

res_bin=readRDS(file = here("results/models/res_Malawi_SL_mom_zinc_bin.rds"))
res=readRDS(file = here("results/models/res_Malawi_SL_mom_zinc.rds"))
res_df = left_join(res$res,
                   res_bin$res %>%
                     rename(outcome_bin=outcome ,Y_bin=Y, yhat_full_bin=yhat_full) %>% mutate(Y_bin=as.numeric(unclass(Y_bin))),
                   by=c("dataid", "clusterid", "population"))
res_df$null_model=as.numeric(res$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$null_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[1,]$MSE)
res_df$sl_model_bin=as.numeric(res_bin$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$MSE)
res_df$sl_model_se=as.numeric(res$cv_risk_w_sl_revere[nrow(res$cv_risk_w_sl_revere),]$se)
auc_pr_brier_res = auc_pr_brier_from_sl3(res_bin, pred_source = c("auto"),ci = TRUE,quiet = TRUE)
res_df= bind_cols(res_df, auc_pr_brier_res)
full_res_malawi=bind_rows(full_res_malawi, res_df)

table(full_res_malawi$Y_bin)


full_res_malawi <- left_join(full_res_malawi, malawi_df, by="dataid")



#-------------------------------------------------------------------------------
# compile and save results
#-------------------------------------------------------------------------------


full_res=bind_rows(full_res_ghana %>% mutate(Country="Ghana"),
                   full_res_gambia %>% mutate(Country="Gambia"),
                   full_res_SL %>% mutate(Country="Sierra Leone"),
                   full_res_malawi %>% mutate(Country="Malawi"))


#save file
saveRDS(full_res, file = here("results/compiled_predictions.RDS"))

