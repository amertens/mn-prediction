

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
df <- readRDS(here("data/IPD/Ghana/Ghana_GMS_admin2_prevalence_merged_dataset.rds"))
#dput(colnames(df))



summary(df$gw_cIDAdjBrinda)
summary(df$gw_wIDAdjBrinda)
summary(df$gw_cVAD)
summary(df$gw_wFolateDef)
summary(df$gw_wB12Def)
summary(df$gw_wVADAdjThurn)





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
gw_vars = gw_vars[!(grepl("RBP", gw_vars))]
gw_vars = gw_vars[!(grepl("rbp", gw_vars))]

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



df_child_VitA <- df %>% select(gw_cVAD, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cVAD))
df_mom_VitA <- df %>% select(gw_wVADAdjThurn, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wVADAdjThurn))

df_child_iron <- df %>% select(gw_cIDAdjBrinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(df$gw_childid),!is.na(gw_cIDAdjBrinda))
df_mom_iron <- df %>% select(gw_wIDAdjBrinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(is.na(df$gw_childid), !is.na(gw_wIDAdjBrinda))

df_mom_b12 <- df %>% select(gw_wB12Def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wB12Def))
df_mom_folate <- df %>% select(gw_wFolateDef, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wFolateDef))




#run parameters
slmodel=slmod2
CVsetting=FALSE
Vfolds=5
id="gw_cnum"

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_vitA_full.rds"))

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_vitA_gwPreds.rds"))

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_vitA.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVADAdjThurn", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_vitA_full.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVADAdjThurn", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_vitA_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVADAdjThurn", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_vitA.rds"))


res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12Def", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_b12_full.rds"))

res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12Def", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_b12_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12Def", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_b12.rds"))


res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolateDef", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_folate_full.rds"))

res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolateDef", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_folate_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolateDef", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_women_folate.rds"))



res=try(DHS_SL(d=df_child_iron, outcome="gw_cIDAdjBrinda", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_iron_full.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="gw_cIDAdjBrinda", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_iron_gwPreds.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="gw_cIDAdjBrinda", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_child_iron.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_wIDAdjBrinda", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_mom_iron_full.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_wIDAdjBrinda", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_mom_iron_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_wIDAdjBrinda", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_admin2_GW_Ghana_SL_mom_iron.rds"))
