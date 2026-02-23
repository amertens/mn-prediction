

rm(list=ls())
library(dplyr)
library(tidyverse)
library(haven)
library(here)
library(purrr)
library(labelled)
library(sl3)
library(origami)
library(tlverse)
library(caret)
library(data.table)
library(ck37r)
library(rdhs)
library(maps)
library(SuperLearner)
library(washb)
options(future.globals.maxSize = 2 * 1024^3)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#read clean and merged data
df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))



dhs_vars <- colnames(df)[grepl("dhs_", colnames(df))]
mics_vars <- colnames(df)[grepl("mics_", colnames(df))]
ihme_vars <- colnames(df)[grepl("ihme_", colnames(df))]
lsms_vars <- colnames(df)[grepl("lsms_", colnames(df))]
map_vars <- colnames(df)[grepl("MAP_", colnames(df))]
wfp_vars <- colnames(df)[grepl("wfp_", colnames(df))]
flunet_vars <- colnames(df)[grepl("flunet_", colnames(df))]
Xvars = c("month",dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars)

df_VitA <- df %>% select(cRBPAdjBrinda, cnum, !!(Xvars)) %>% as.data.frame() %>% filter(!is.na(cRBPAdjBrinda))

res=readRDS(here("results/models/res_GW_Ghana_SL.rds"))

#evaluate performance
sl_performance_res = res$cv_risk_w_sl_revere

#RMSE
sl_performance_res$RMSE = sqrt(sl_performance_res$MSE)
summary(df$cRBPAdjBrinda)
sd(df$cRBPAdjBrinda, na.rm=T)

#RMSE = SD is equal to a null model
Ysd=sd(df$cRBPAdjBrinda, na.rm=T)
sl_RMSE=sl_performance_res$RMSE[length(sl_performance_res$RMSE)]
sl_RMSE/Ysd

# As a rule of thumb:
#
#   RMSE = 0.9 × std dev: Slight improvement over baseline
# RMSE = 0.7 × std dev: Decent model
# RMSE = 0.5 × std dev: Good model
# RMSE = 0.3 × std dev: Very good model

# #Brier skill score - for binomical outcomes
# #*BSS=1→ perfect; 0→ no better than null; <0→ worse than null.*
# prev=mean(df$mod_sev_anemia)
# prev
# null_brier = prev*(1-prev)
# BSS = 1-(0.1907480 /null_brier)
# BSS


# #get outcome predictions
yhat_full <- res$sl_fit$predict_fold(res$task,"validation")
plotdf=data.frame(Y=df_VitA$cRBPAdjBrinda, predY=yhat_full)
ggplot(plotdf, aes(y=predY, x=Y)) + geom_point()




# to do: have missing variables together


# res_vim=try(calc_importance(res$sl_fit))
# res_vim

res_vim_diff=try(calc_importance(res$sl_fit,  eval.fun = loss_squared_error, importance.metric="difference"))
res_vim_diff

res_vim_diff$varimp


demographics=c("sex", "age","age_single", "age_5yrgrp", "rel_to_head", "usual_resident", "slept_last_night","religion", "ethnicity", "dv_selected")
education =  c( "edu_level", "edu_years")
household = c("residence","residence_type","hh_size", "num_under5","floor_material", "wall_material", "roof_material",
              "water_source", "toilet_type", "share_toilet", "cooking_fuel")
wealth = c("wealth_quint", "wealth_score", "has_electricity",
           "has_radio", "has_tv", "has_fridge", "has_bicycle", "has_motorcycle",
           "has_car","has_mobile", "has_computer")

geography = c("ALT_GPS", "region","altitude", "trmm_mean_10km","trmm_mean_25km", "trmm_mean_50km")

pregnancy_parity = c("pregnant_now", "months_preg", "months_last_bir",
                     "total_living_ch", "living_ch_plus")

food_price=c("cassava", "cowpeas", "eggplants", "eggs", "fish_mackerel_fresh",
             "gari", "maize", "maize_yellow", "meat_chicken", "meat_chicken_local",
             "millet", "onions", "peppers_dried", "peppers_fresh", "plantains_apem",
             "plantains_apentu", "rice_imported", "rice_local", "sorghum",
             "soybeans", "tomatoes_local", "tomatoes_navrongo", "yam", "yam_puna",
             "cowpeas_white")


ruralness = c("Global_Travel_Speed_Friction_Surface", "Global_Pv_Parasite_Rate",
              "Global_Travel_Time_to_Cities", "Global_Motorized_Friction_Surface",
              "Global_Motorized_Travel_Time_to_Healthcare", "Global_Walking_Only_Friction_Surface",
              "Global_Walking_Only_Travel_Time_To_Healthcare",
              "pop_mean_10km", "pop_mean_25km",  "pop_mean_50km",
              "nearest_market_distance_km", "nearest_market_id")

malaria_burden =c("Africa_HbC_Allele_Frequency","Global_Duffy_Negativity_Phenotype_Frequency", "Global_G6PDd_Allele_Frequency",
                  "Global_Sickle_Haemoglobin_HbS_Allele_Frequency", "Global_Pf_Reproductive_Number",
                  "Global_Pf_Incidence_Count", "Global_Pf_Incidence_Rate", "Global_Pf_Mortality_Count",
                  "Global_Pf_Mortality_Rate", "Global_Pf_Parasite_Rate", "Global_Pv_Incidence_Count",
                  "Global_Pv_Incidence_Rate")

malaria_treatment =  c( "has_mosquito_net", "used_bednet","Africa_Insecticide_Treated_Net_Access",
                        "Africa_Insecticide_Treated_Net_Use_Rate", "Africa_IRS_Coverage",
                        "Global_Antimalarial_Effective_Treatment")
missing_vars= res$Xvars[grepl("missing_",res$Xvars)]

res_vim_diff_groups=try(calc_importance(res$sl_fit,  eval.fun = loss_squared_error, importance.metric="difference",
                                        covariate.groups=list("demographics" = demographics[demographics %in% res$Xvars],
                                                              "education" =education[education %in% res$Xvars],
                                                              "household" = household[household %in% res$Xvars],
                                                              "wealth" = wealth[wealth %in% res$Xvars],
                                                              "geography" = geography[geography %in% res$Xvars],
                                                              "pregnancy_parity" = pregnancy_parity[pregnancy_parity %in% res$Xvars],
                                                              "food_price" = food_price[food_price %in% res$Xvars],
                                                              "ruralness" = ruralness[ruralness %in% res$Xvars],
                                                              "malaria_burden" = malaria_burden[malaria_burden %in% res$Xvars],
                                                              "malaria_treatment" = malaria_treatment[malaria_treatment %in% res$Xvars],
                                                              "missingness" =missing_vars
                                        )))
res_vim_diff_groups




res_vim_groups=try(calc_importance(res$sl_fit,
                                   covariate.groups=list("demographics" = demographics[demographics %in% res$Xvars],
                                                         "education" =education[education %in% res$Xvars],
                                                         "household" = household[household %in% res$Xvars],
                                                         "wealth" = wealth[wealth %in% res$Xvars],
                                                         "geography" = geography[geography %in% res$Xvars],
                                                         "pregnancy_parity" = pregnancy_parity[pregnancy_parity %in% res$Xvars],
                                                         "food_price" = food_price[food_price %in% res$Xvars],
                                                         "ruralness" = ruralness[ruralness %in% res$Xvars],
                                                         "malaria_burden" = malaria_burden[malaria_burden %in% res$Xvars],
                                                         "malaria_treatment" = malaria_treatment[malaria_treatment %in% res$Xvars],
                                                         "missingness" =missing_vars
                                   )))
res_vim_groups



#to do: group all missingness variables,
# group all food pricing variables
# all asset variables
#combine wealth variables (wealth score, quintile, and assets)

#"months last bir" = interval since last birth

#use DHS_variable_recode to relable variables

#---------------------------------------------------------------------------------
#rename_vars_for_plotting
#---------------------------------------------------------------------------------

extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec))

res_vim_diff_plot_df <- res_vim_diff$varimp
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p=ggplot(res_vim_diff_plot_df %>%
           arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "Top 20 Variables by effect on MSE",
       x = "Covariate",
       y = "MSE_difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#save models
saveRDS(res, here("results/DHS_Ghana_SL_model.RDS"))

#save variable importance results
saveRDS(res_vim_diff, here("results//models/DHS_Ghana_SL_vim_diff.RDS"))
saveRDS(res_vim, here("results//models/DHS_Ghana_SL_vim.RDS"))
saveRDS(res_vim_diff_groups, here("results//models/DHS_Ghana_SL_grouped_vim_diff.RDS"))
saveRDS(res_vim_groups, here("results//models/DHS_Ghana_SL_grouped_vim.RDS"))


#save plots
ggsave(p, file=here("figures/DHS_Ghana_SL_vim_diff.png"),device = "png", width = 10,height = 8,units = "in",dpi = 300)

