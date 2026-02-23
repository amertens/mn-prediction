

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
# Use: fer (serum ferritin, inflammation-adjusted)
table(!is.na(df$fer), df$population)

# Vitamin A status
# Use: rbp (retinol-binding protein, inflammation-adjusted)
table(!is.na(df$rbp), df$population)

# Zinc status
# Use: zn_gdl (serum zinc, µg/dL)
table(!is.na(df$zn_gdl), df$population)



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
  select(rbp, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(rbp))
df_mom_VitA <- df %>% filter(population %in% c("women")) %>%
  select(rbp, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(rbp))

df_child_iron <- df %>% filter(population %in% c("preschool children", "school-age children")) %>%
  select(fer, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(fer))
df_mom_iron <- df %>% filter(population %in% c("women")) %>%
  select(fer, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(fer))

df_child_zinc <- df %>% filter(population %in% c("preschool children", "school-age children")) %>%
  select(zn_gdl, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(zn_gdl))
df_mom_zinc <- df %>% filter(population %in% c("women")) %>%
  select(zn_gdl, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(zn_gdl))


#run parameters
#slmodel=sl_simple
#slmodel=slfull
slmodel=slmod
CVsetting=FALSE
#Vfolds=5
Vfolds=2
id="gw_cnum"


res=try(DHS_SL(d=df_child_VitA, outcome="rbp", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_vitA.rds"))

res=try(DHS_SL(d=df_mom_VitA, outcome="rbp", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_women_vitA.rds"))

res=try(DHS_SL(d=df_child_iron, outcome="fer", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_iron.rds"))

res=try(DHS_SL(d=df_mom_iron, outcome="fer", population="women",Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_mom_iron.rds"))

res=try(DHS_SL(d=df_child_zinc, outcome="zn_gdl", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_child_zinc.rds"))

res=try(DHS_SL(d=df_mom_zinc, outcome="zn_gdl", population="women",Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res, file = here("results/models/res_Malawi_SL_mom_zinc.rds"))



