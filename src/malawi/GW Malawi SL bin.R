

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
df <- readRDS(here("data/Malawi/Malawi_merged_dataset.rds"))
df <- sf::st_drop_geometry(df)

#make unique id
df$dataid <- paste0("malawi",1:nrow(df))

colnames(df)



# Iron status
table(!is.na(df$iron_def), df$population)

# Vitamin A status
table(!is.na(df$vitA_def), df$population)

# Zinc status
table(!is.na(df$zinc_def), df$population)



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

Xvars_full = c("dataid","Admin1","Admin2",gw_vars,dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
Xvars_gw = c("dataid","Admin1","Admin2",gw_vars)
Xvars = c("dataid","Admin1","Admin2",dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)


length(gw_vars)+1
length(dhs_vars)
length(mics_vars)
length(lsms_vars)
length(ihme_vars)
length(map_vars)
length(wfp_vars)
length(flunet_vars)
length(gee_vars)



df_child_VitA <- df %>% filter(population %in% c("preschool children", "school-age children")) %>%
  select(vitA_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(vitA_def))
df_mom_VitA <- df %>% filter(population %in% c("women")) %>%
  select(vitA_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(vitA_def))

df_child_iron <- df %>% filter(population %in% c("preschool children", "school-age children")) %>%
  select(iron_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(iron_def))
df_mom_iron <- df %>% filter(population %in% c("women")) %>%
  select(iron_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(iron_def))

df_child_zinc <- df %>% filter(population %in% c("preschool children", "school-age children")) %>%
  select(zinc_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(zinc_def))
df_mom_zinc <- df %>% filter(population %in% c("women")) %>%
  select(zinc_def, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(zinc_def))


#run parameters
#slmodel=sl_simple
slmodel=slfull_bin
#slmodel=slmod2_bin
CVsetting=FALSE
#Vfolds=5
Vfolds=4
id="gw_cnum"


res=try(DHS_SL(d=df_child_VitA, outcome="vitA_def", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_vitA_bin.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="vitA_def", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_women_vitA_bin.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="iron_def", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_iron_bin.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="iron_def", population="women",Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_mom_iron_bin.rds"))

res=try(DHS_SL(d=df_child_zinc, outcome="zinc_def", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_zinc_bin.rds"))

res=try(DHS_SL(d=df_mom_zinc, outcome="zinc_def", population="women",Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_mom_zinc_bin.rds"))



