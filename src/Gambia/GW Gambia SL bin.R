

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
df <- readRDS(here("data/IPD/Gambia/Gambia_merged_dataset.rds"))
df <- sf::st_drop_geometry(df)

#make unique id
df$dataid <- paste0("gambia",1:nrow(df))

#dput(colnames(df))

colnames(df)

#vitamin A data
table(!is.na(df$gw_cVAD_Thurn), df$gw_child_flag)
table(!is.na(df$gw_wVAD_Thurn), df$gw_child_flag)

table((df$gw_cVAD_Thurn), df$gw_child_flag)
table((df$gw_wVAD_Thurn), df$gw_child_flag)


#iron data
table(!is.na(df$gw_wIDA_Brinda), df$gw_child_flag)
table(!is.na(df$gw_cIDA_Brinda), df$gw_child_flag)



colnames(df)[grepl("Brinda",colnames(df))]
colnames(df)[grepl("VAD",colnames(df))]
colnames(df)[grepl("Thurn",colnames(df))]

#iron data
table(is.na(df$gw_LogFerAdj), df$gw_child_flag)


gw_vars <- colnames(df)[grepl("gw_", colnames(df))]

dhs_vars <- colnames(df)[grepl("dhs2019_", colnames(df))]
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

Xvars_full = c("dataid","Admin1","gw_month", gw_vars, dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
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



df_child_VitA <- df %>% filter(gw_child_flag==1) %>% select(gw_cVAD_Thurn, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cVAD_Thurn))
df_mom_VitA <- df %>% filter(gw_child_flag==0) %>% select(gw_wVAD_Thurn, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wVAD_Thurn))

df_child_iron <- df %>% filter(gw_child_flag==1) %>% select(gw_cIDA_Brinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cIDA_Brinda))
df_mom_iron <- df %>% filter(gw_child_flag==0) %>% select(gw_wIDA_Brinda, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wIDA_Brinda))


#run parameters
#slmodel=slfull_bin
slmodel=slmod2_bin
CVsetting=FALSE
Vfolds=5
#Vfolds=2
id="gw_cnum"



res=try(DHS_SL(d=df_child_iron, outcome="gw_cIDA_Brinda", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_bin_GW_Gambia_SL_child_iron.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="gw_wIDA_Brinda", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_bin_GW_Gambia_SL_mom_iron.rds"))

res=try(DHS_SL(d=df_child_VitA, outcome="gw_cVAD_Thurn", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_bin_GW_Gambia_SL_child_vitA.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVAD_Thurn", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_bin_GW_Gambia_SL_women_vitA.rds"))


