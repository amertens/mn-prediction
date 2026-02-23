

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

#iron data
table(is.na(df$gw_wFeDefAdjBR1), df$gw_child_flag)
table(is.na(df$gw_cFeDefAdj), df$gw_child_flag)

#vitamin A data
table(is.na(df$gw_wVitADefAdjBR1), df$gw_child_flag)
table(is.na(df$gw_cVitADef), df$gw_child_flag)

#iodine
table(is.na(df$gw_UIAdeq), df$gw_child_flag)


#continious values
# | Micronutrient | Population | Continuous variable (preferred)   | Scale         |
#   | ------------- | ---------- | --------------------------------- | ------------- |
#   | Iron          | Women      | `gw_LNwFerAdjBR1` / `gw_wFerrAdj` | Log preferred |
#   | Iron          | Children   | `gw_cFerrAdj`                     | Log preferred |
#   | Vitamin A     | Women      | `gw_wRBPAdjBR1`                   | Log optional  |
#   | Vitamin A     | Children   | `gw_cRBPAdj`                      | Log optional  |
#   | Iodine        | All        | `gw_UIAliqALog`                   | Log           |
#   | Folate        | Women      | `gw_wFolate`                      | Raw or log    |
#   | Vitamin B12   | Women      | `gw_B12`                          | Raw or log    |
#   | Anemia        | All        | `gw_wHb`, `gw_cHb`                | Raw           |
#


#iron data
table(!is.na(df$gw_wFerrAdjThurn), df$gw_child_flag)
table(!is.na(df$gw_cFerrAdj), df$gw_child_flag)

#vitamin A data
table(is.na(df$gw_wRBPAdjBR1), df$gw_child_flag)
table(is.na(df$gw_cRBPAdj), df$gw_child_flag)

#iodine
table(is.na(df$gw_UIAliqALog), df$gw_child_flag)
#b12
table(is.na(df$gw_B12), df$gw_child_flag)
#folate
table(is.na(df$gw_wFolate), df$gw_child_flag)



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



df_child_VitA <- df %>% filter(gw_child_flag==1) %>% select(gw_cFerrAdj, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_cFerrAdj))
df_mom_VitA <- df %>% filter(gw_child_flag==0) %>% select(gw_wRBPAdjBR1, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wRBPAdjBR1))

df_child_iron <- df %>% filter(gw_child_flag==1) %>%
  select(gw_cFerrAdj, gw_cnum, !!(Xvars_full)) %>%
  mutate(gw_cFerrAdj=log(gw_cFerrAdj)) %>%
  as.data.frame() %>% filter(!is.na(gw_cFerrAdj))
df_mom_iron <- df %>% filter(gw_child_flag==0) %>%
  select(gw_wFerrAdjThurn, gw_cnum, !!(Xvars_full)) %>%
  mutate(gw_wFerrAdjThurn=log(gw_wFerrAdjThurn)) %>%
  as.data.frame() %>% filter(!is.na(gw_wFerrAdjThurn))

df_mom_iodine <- df %>% filter(gw_child_flag==0) %>% select(gw_UIAliqALog, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_UIAliqALog))
df_mom_b12<- df %>% filter(gw_child_flag==0) %>% select(gw_B12, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_B12))
df_mom_folate <- df %>% filter(gw_child_flag==0) %>% select(gw_wFolate, gw_cnum, !!(Xvars_full)) %>% as.data.frame() %>% filter(!is.na(gw_wFolate))




#Xvars=Xvars[950:1000]

#run parameters
#slmodel=sl_simple
slmodel=slmod
CVsetting=FALSE
Vfolds=2 #temp to debug error
#Vfolds=5
id="gw_cnum"

d=df_child_VitA
outcome="gw_cFerrAdj"
population="children"
folds=Vfolds
CV=CVsetting
sl=slmodel
prescreen=T

res1=try(DHS_SL(d=df_child_VitA, outcome="gw_cFerrAdj", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res1, file = here("results/models/res_GW_Sierra_Leone_SL_child_vitA.rds"))

res2=try(DHS_SL(d=df_mom_VitA, outcome="gw_wRBPAdjBR1", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res2, file = here("results/models/res_GW_Sierra_Leone_SL_women_vitA.rds"))

res3=try(DHS_SL(d=df_child_iron, outcome="gw_cFerrAdj", population="children", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res3, file = here("results/models/res_GW_Sierra_Leone_SL_child_iron.rds"))

res4=try(DHS_SL(d=df_mom_iron, outcome="gw_wFerrAdjThurn", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res4, file = here("results/models/res_GW_Sierra_Leone_SL_mom_iron.rds"))

res5=try(DHS_SL(d=df_mom_iodine, outcome="gw_UIAliqALog", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res5, file = here("results/models/res_GW_Sierra_Leone_SL_mom_iodine.rds"))

res6=try(DHS_SL(d=df_mom_b12, outcome="gw_B12", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res6, file = here("results/models/res_GW_Sierra_Leone_SL_mom_b12.rds"))

res7=try(DHS_SL(d=df_mom_folate, outcome="gw_wFolate", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res7, file = here("results/models/res_GW_Sierra_Leone_SL_mom_folate.rds"))
