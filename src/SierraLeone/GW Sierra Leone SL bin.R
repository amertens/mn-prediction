

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
df <- readRDS(here("data/IPD/Sierra Leone/SierraLeone_merged_dataset.rds"))
df <- sf::st_drop_geometry(df)

colnames(df)

#make unique id
df$dataid <- paste0("SL",1:nrow(df))


#Deficiencies

# | Micronutrient          | Population | Primary variable    |
#   | ---------------------- | ---------- | ------------------- |
#   | Iron deficiency        | Women      | `gw_wFeDefAdjBR1`   |
#   | Iron deficiency        | Children   | `gw_cFeDefAdj`      |
#   | Vitamin A deficiency   | Women      | `gw_wVitADefAdjBR1` |
#   | Vitamin A deficiency   | Children   | `gw_cVitADef`       |
#   | Iodine deficiency      | All        | `gw_UIAdeq`         |
#   | Folate deficiency      | Women      | `gw_wFolDef`        |
#   | Vitamin B12 deficiency | Women      | `gw_wB12DefWHO`     |

attributes(df$gw_cFeDefAdj)
attributes(df$gw_wFeDefAdjBR1)

table(df$gw_wFeDefAdjBR1 )
df$gw_wFeDefAdjBR1 <- ifelse(df$gw_wFeDefAdjBR1==1,1,0)
table(df$gw_wFeDefAdjBR1 )

df$gw_cFeDefAdj <- ifelse(df$gw_cFeDefAdj==1,1,0)
df$gw_wVitADefAdjBR1 <- ifelse(df$gw_wVitADefAdjBR1==1,1,0)
df$gw_cVitADef <- ifelse(df$gw_cVitADef==1,1,0)
df$gw_UIAdeq <- ifelse(df$gw_UIAdeq==1,1,0)
df$gw_wFolDef <- ifelse(df$gw_wFolDef==1,1,0)
df$gw_wB12DefWHO <- ifelse(df$gw_wB12DefWHO==1,1,0)




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

Xvars_full = c("dataid","Admin1",gw_vars,dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
Xvars_gw = c("dataid","Admin1",gw_vars)
Xvars = c("dataid","Admin1",dhs_vars, mics_vars, ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)


length(gw_vars)+1
length(dhs_vars)
length(mics_vars)
length(lsms_vars)
length(ihme_vars)
length(map_vars)
length(wfp_vars)
length(flunet_vars)
length(gee_vars)



df_child_VitA <- df %>% filter(gw_child_flag==1) %>% select(gw_cVitADef, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cVitADef))
df_mom_VitA <- df %>% filter(gw_child_flag==0) %>% select(gw_wVitADefAdjBR1, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wVitADefAdjBR1))

df_child_iron <- df %>% filter(gw_child_flag==1) %>%
  select(gw_cFeDefAdj, gw_cnum, !!(Xvars_full)) %>%
  as.data.frame() %>% filter(!is.na(gw_cFeDefAdj))
df_mom_iron <- df %>% filter(gw_child_flag==0) %>%
  select(gw_wFeDefAdjBR1, gw_cnum, !!(Xvars_full)) %>%
  as.data.frame() %>% filter(!is.na(gw_wFeDefAdjBR1))

df_mom_iodine <- df %>% filter(gw_child_flag==0) %>% select(gw_UIAdeq, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_UIAdeq))
df_mom_b12<- df %>% filter(gw_child_flag==0) %>% select(gw_wB12DefWHO, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wB12DefWHO))
df_mom_folate <- df %>% filter(gw_child_flag==0) %>% select(gw_wFolDef, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wFolDef))




#Xvars=Xvars[950:1000]

#run parameters
#slmodel=sl_simple_bin
#slmodel=slfull_bin
slmodel=slmod2_bin

CVsetting=FALSE
#Vfolds=2 #temp to debug error
Vfolds=5
id="gw_cnum"

d=df_child_VitA
outcome="gw_cFerrAdj"
population="children"
folds=Vfolds
CV=CVsetting
sl=slmodel
prescreen=T

res1=try(DHS_SL(d=df_child_VitA, outcome="gw_cVitADef", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=sl_simple_bin))
saveRDS(res1, file = here("results/models/res_GW_Sierra_Leone_SL_child_vitA_bin.rds"))

res2=try(DHS_SL(d=df_mom_VitA, outcome="gw_wVitADefAdjBR1", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=sl_simple_bin))
saveRDS(res2, file = here("results/models/res_GW_Sierra_Leone_SL_women_vitA_bin.rds"))

res3=try(DHS_SL(d=df_child_iron, outcome="gw_cFeDefAdj", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res3, file = here("results/models/res_GW_Sierra_Leone_SL_child_iron_bin.rds"))

res4=try(DHS_SL(d=df_mom_iron, outcome="gw_wFeDefAdjBR1", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res4, file = here("results/models/res_GW_Sierra_Leone_SL_mom_iron_bin.rds"))

res5=try(DHS_SL(d=df_mom_iodine, outcome="gw_UIAdeq", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res5, file = here("results/models/res_GW_Sierra_Leone_SL_mom_iodine_bin.rds"))

res6=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12DefWHO", population="women", Xvars=Xvars, id=id, folds=2, CV=CVsetting, sl=sl_simple_bin))
# #fails due to sparsity
# table(df_mom_b12$gw_wB12DefWHO)
# set.seed(1232435)
# res6=try(DHS_SL(d=df_mom_b12, outcome="gw_wB12DefWHO", population="women", Xvars=Xvars, id=id, folds=2, CV=CVsetting, sl=slmodel))
saveRDS(res6, file = here("results/models/res_GW_Sierra_Leone_SL_mom_b12_bin.rds"))

res7=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolDef", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res7, file = here("results/models/res_GW_Sierra_Leone_SL_mom_folate_bin.rds"))



