

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

#make unique id
df$dataid <- paste0("ghana",1:nrow(df))



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



df_child_VitA <- df %>% select(dataid, gw_cRBPAdjBrinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cRBPAdjBrinda))
df_mom_VitA <- df %>% select(dataid, gw_wRBP, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wRBP))

df_child_iron <- df %>% select(dataid, gw_r_crpagp_sf2, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(df$gw_childid),!is.na(gw_r_crpagp_sf2))
df_mom_iron <- df %>% select(dataid, gw_r_crpagp_sf2, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(is.na(df$gw_childid), !is.na(gw_r_crpagp_sf2))

df_mom_b12 <- df %>% select(dataid, gw_wB12, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wB12))
df_mom_folate <- df %>% select(dataid, gw_wFolate, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wFolate))


#run parameters
slmodel=slfull
CVsetting=FALSE
Vfolds=5
id="gw_cnum"

d=df_child_VitA
outcome="gw_cRBPAdjBrinda"
Xvars=Xvars
id=id
folds=Vfolds
CV=CVsetting
sl=slmodel

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cRBPAdjBrinda", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_vitA_full.rds"))

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cRBPAdjBrinda", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_vitA_gwPreds.rds"))

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cRBPAdjBrinda", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_vitA.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wRBP", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_vitA_full.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wRBP", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_vitA_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wRBP", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_vitA.rds"))


res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_b12_full.rds"))

res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_b12_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_b12.rds"))


res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolate", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_folate_full.rds"))

res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolate", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_folate_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolate", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_women_folate.rds"))



res=try(DHS_SL(d=df_child_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_iron_full.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_iron_gwPreds.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_child_iron.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars_full, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_mom_iron_full.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars_gw, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_mom_iron_gwPreds.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_r_crpagp_sf2", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_GW_Ghana_SL_mom_iron.rds"))
