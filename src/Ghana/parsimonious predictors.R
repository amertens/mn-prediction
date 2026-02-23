


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
library(future)
library(washb)
library(recipes)

plan(multicore, workers = availableCores()/2)
options(future.globals.maxSize = 5 * 1024^3)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#read clean and merged data
df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))
#dput(colnames(df))


df$gw_wVAD <- as.numeric(df$gw_wVAD)



#run parameters
slmodel=slfull_bin
CVsetting=FALSE
Vfolds=5
id="gw_cnum"


full_model = readRDS( here("results/models/res_full_bin_GW_Ghana_SL_child_vitA_full.rds"))
full_model$cv_risk_w_sl_revere
vim_res = readRDS( here("results/vim_child_vitA_def_proxy_individ.rds"))


Xvars=vim_res$varimp20$covariate
df_child_VitA <- df %>% select(gw_cVAD, gw_cnum, !!(Xvars)) %>% as.data.frame() %>% filter(!is.na(gw_cVAD))

res_top20=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))

full_model$cv_risk_w_sl_revere
res_top20$cv_risk_w_sl_revere


full_model_iron = readRDS( here("results/models/res_full_bin_GW_Ghana_SL_child_iron_full.rds"))
full_model_iron$cv_risk_w_sl_revere
vim_res_iron = readRDS( here("results/vim_child_iron_def_proxy_individ.rds"))
Xvars_iron=vim_res_iron$varimp20$covariate
Xvars_iron=Xvars_iron[!(grepl("missing_",Xvars_iron))]

df_child_iron <- df %>% select(gw_cIDAdjBrinda, gw_cnum, !!(Xvars_iron)) %>% as.data.frame() %>% filter(!is.na(gw_cIDAdjBrinda))
res_top20_iron=try(DHS_SL(d=df_child_iron, outcome="gw_cIDAdjBrinda", Xvars=Xvars_iron, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))

full_model_iron$cv_risk_w_sl_revere
res_top20_iron$cv_risk_w_sl_revere


res_vim_diff_child_iron = readRDS(file= here("results/vim_child_iron_def_proxy_individ.rds"))
res_vim_diff_women_vitA=readRDS(file= here("results/vim_women_vitA_def_proxy_individ.rds"))
res_vim_diff_women_iron=readRDS(file= here("results/vim_women_iron_def_proxy_individ.rds"))
res_vim_diff_women_folate=readRDS(file= here("results/vim_women_folate_def_proxy_individ.rds"))
res_vim_diff_women_b12=readRDS(file= here("results/vim_women_b12_def_proxy_individ.rds"))

