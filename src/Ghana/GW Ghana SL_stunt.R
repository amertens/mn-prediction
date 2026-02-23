

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


#vitamin A data
summary(df$gw_waz06)

gw_vars <- colnames(df)[grepl("gw_", colnames(df))]

dhs_vars <- colnames(df)[grepl("dhs2014_", colnames(df))| grepl("dhs2016_", colnames(df))| grepl("dhs2017_", colnames(df))]
mics_vars <- colnames(df)[grepl("mics_", colnames(df))]
ihme_vars <- colnames(df)[grepl("ihme_", colnames(df))]
lsms_vars <- colnames(df)[grepl("lsms_", colnames(df))]
map_vars <- colnames(df)[grepl("MAP_", colnames(df))]
wfp_vars <- c("nearest_market_id", colnames(df)[grepl("wfp_", colnames(df))])
flunet_vars <- colnames(df)[grepl("flunet_", colnames(df))]
gee_vars <- colnames(df)[grepl("gee_", colnames(df))]

#drop vars related to the outome
gw_vars <-gw_vars[!grepl("haz", gw_vars)]
gw_vars <-gw_vars[!grepl("waz", gw_vars)]
gw_vars <-gw_vars[!grepl("whz", gw_vars)]
gw_vars <-gw_vars[!grepl("Stunt", gw_vars)]
gw_vars <-gw_vars[!grepl("Wast", gw_vars)]
gw_vars <-gw_vars[!grepl("Over", gw_vars)]
gw_vars <-gw_vars[!grepl("Under", gw_vars)]
gw_vars <-gw_vars[!grepl("eight", gw_vars)]

Xvars_full = c("Admin1","gw_month",gw_vars,dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
Xvars_gw = c("Admin1","gw_month",gw_vars)
Xvars = c("Admin1","gw_month",dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)


length(gw_vars)+1
length(dhs_vars)
length(mics_vars)
length(lsms_vars)
length(ihme_vars)
length(map_vars)
length(wfp_vars)
length(flunet_vars)
length(gee_vars)



df_child_waz <- df %>% select(gw_waz06, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_waz06))

#run parameters
slmodel=slmod2
CVsetting=FALSE
Vfolds=2
id="gw_cnum"

res_full=try(DHS_SL(d=df_child_waz, outcome="gw_waz06", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_full, file = here("results/models/res_GW_Ghana_SL_child_waz_full.rds"))

res_gw=try(DHS_SL(d=df_child_waz, outcome="gw_waz06", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_gw, file = here("results/models/res_GW_Ghana_SL_child_waz_gwPreds.rds"))

res_child=try(DHS_SL(d=df_child_waz, outcome="gw_waz06", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_child, file = here("results/models/res_GW_Ghana_SL_child_waz.rds"))



sd(df_child_waz$gw_waz06, na.rm = TRUE)
sqrt(res_full$cv_risk_w_sl_revere$MSE[6])
sqrt(res_gw$cv_risk_w_sl_revere$MSE[6])
sqrt(res_child$cv_risk_w_sl_revere$MSE[6])



