

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

#make unique id
df$dataid <- paste0("ghana",1:nrow(df))

#vitamin A data
summary(df$gw_wVAD)
summary(df$gw_cVAD)

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
gw_vars = gw_vars[!(grepl("wID", gw_vars))]
gw_vars = gw_vars[!(grepl("cID", gw_vars))]
gw_vars = gw_vars[!(grepl("VAD", gw_vars))]
gw_vars = gw_vars[!(grepl("RBP", gw_vars))]
gw_vars = gw_vars[!(grepl("rbp", gw_vars))]
gw_vars = gw_vars[!(grepl("Ferr", gw_vars))]
gw_vars = gw_vars[!(grepl("TFR", gw_vars))]
gw_vars = gw_vars[!(grepl("crp", gw_vars))]
gw_vars = gw_vars[!(grepl("agp", gw_vars))]
gw_vars = gw_vars[!(grepl("cRBP", gw_vars))]
gw_vars = gw_vars[!(grepl("cVAD", gw_vars))]
gw_vars = gw_vars[!(grepl("VAI", gw_vars))]
gw_vars = gw_vars[!(grepl("RBP", gw_vars))]
gw_vars = gw_vars[!(grepl("Folate", gw_vars))]
gw_vars = gw_vars[!(grepl("B12", gw_vars))]
gw_vars = gw_vars[!(grepl("AGP", gw_vars))]
gw_vars = gw_vars[!(grepl("CRP", gw_vars))]
gw_vars = gw_vars[!(grepl("RDT", gw_vars))]
gw_vars = gw_vars[!(grepl("fer", gw_vars))]
gw_vars = gw_vars[!(grepl("bis", gw_vars))]
gw_vars = gw_vars[!(grepl("nflam", gw_vars))]
gw_vars = gw_vars[!(grepl("nemia", gw_vars))]
gw_vars = gw_vars[!(grepl("globin", gw_vars))]
gw_vars = gw_vars[!(grepl("cHb", gw_vars))]
gw_vars = gw_vars[!(grepl("wHb", gw_vars))]
gw_vars = gw_vars[!(grepl("gchb", gw_vars))]
gw_vars = gw_vars[!(grepl("hbc", gw_vars))]
gw_vars = gw_vars[!(grepl("Anemia", gw_vars))]
gw_vars = gw_vars[!(grepl("wm_wmst", gw_vars))]
gw_vars = gw_vars[gw_vars!="gw_cn" & gw_vars!="gw_hhid"& gw_vars!="gw_childid"]


gsub("gw_", "", gw_vars) #remove gw_ prefix for consistency with other variables

Xvars_full = c("dataid","Admin1","gw_month",gw_vars,dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
Xvars_gw = c("dataid","Admin1","gw_month",gw_vars)
Xvars = c("dataid","Admin1","gw_month",dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)


length(gw_vars)+1
length(dhs_vars)
length(mics_vars)
length(lsms_vars)
length(ihme_vars)
length(map_vars)
length(wfp_vars)
length(flunet_vars)
length(gee_vars)



df$gw_wVAD <- as.numeric(df$gw_wVAD)

df_child_VitA <- df %>% select(gw_cVAD, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cVAD))
df_mom_VitA <- df %>% select(gw_wVAD, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wVAD))

df_child_iron <- df %>% select(gw_cIDAdjBrinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(df$gw_childid),!is.na(gw_cIDAdjBrinda))
df_mom_iron <- df %>% select(gw_wIDAdjBrinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(is.na(df$gw_childid), !is.na(gw_wIDAdjBrinda))

df_mom_b12 <- df %>% select(gw_wB12Def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wB12Def))
df_mom_folate <- df %>% select(gw_wFolateDef, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wFolateDef))


#run parameters
# slmodel=sl_simple_bin
# slmodel=slfull_bin
slmodel=slmod_bin
CVsetting=FALSE
Vfolds=2
id="gw_cnum"

res1=res2=res3=res4=res5=res6=res7=res8=res9=res10=res11=res12=res13=res14=res15=res16=res17=res18=NULL

res3=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res3, file = here("results/models/res_bin_GW_Ghana_SL_child_vitA_V2.rds"))

res6=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVAD", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res6, file = here("results/models/res_bin_GW_Ghana_SL_women_vitA_V2.rds"))

res9=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12Def", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res9, file = here("results/models/res_bin_GW_Ghana_SL_women_b12_V2.rds"))

res12=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolateDef", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res12, file = here("results/models/res_bin_GW_Ghana_SL_women_folate_V2.rds"))

res15=try(DHS_SL(d=df_child_iron, outcome="gw_cIDAdjBrinda", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res15, file = here("results/models/res_bin_GW_Ghana_SL_child_iron_V2.rds"))

res18=try(DHS_SL(d=df_mom_iron, outcome="gw_wIDAdjBrinda", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res18, file = here("results/models/res_bin_GW_Ghana_SL_mom_iron_V2.rds"))

