

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
df_gambia <- readRDS(here("data/IPD/Gambia/Gambia_merged_dataset.rds"))  %>% filter(gw_child_flag==0) %>% rename(rbp=gw_wRBP) %>% mutate(country="Gambia") %>% filter(!is.na(rbp))
df_ghana <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))  %>% filter(!is.na(gw_wRBP)) %>% rename(rbp=gw_wRBP) %>% mutate(country="Ghana")%>% filter(!is.na(rbp))
df_malawi <- readRDS(here("data/Malawi/Malawi_merged_dataset.rds"))  %>% filter(population %in% c("women")) %>% mutate(country="Malawi")%>% filter(!is.na(rbp))
df_sierraleone <- readRDS(here("data/IPD/Sierra Leone/SierraLeone_merged_dataset.rds"))   %>% filter(gw_child_flag==0)  %>% rename(rbp=gw_wRBP) %>% mutate(country="Sierra Leone")%>% filter(!is.na(rbp))

df <- bind_rows(df_gambia, df_ghana, df_malawi, df_sierraleone)


gw_vars <- colnames(df)[grepl("gw_", colnames(df))]

gee_vars <- colnames(df)[grepl("gee_", colnames(df))]

gee_vars_gambia <- colnames(df_gambia)[grepl("gee_", colnames(df_gambia))]
gee_vars_ghana <- colnames(df_ghana)[grepl("gee_", colnames(df_ghana))]
gee_vars_malawi <- colnames(df_malawi)[grepl("gee_", colnames(df_malawi))]
gee_vars_sierraleone <- colnames(df_sierraleone)[grepl("gee_", colnames(df_sierraleone))]

gee_vars=gee_vars[gee_vars %in% gee_vars_gambia &
          gee_vars %in% gee_vars_ghana &
          gee_vars %in% gee_vars_malawi &
          gee_vars %in% gee_vars_sierraleone]
length(gee_vars)

df_gambia <- df_gambia %>% select(rbp, dataid, country,  gw_cnum, !!(gee_vars)) %>% as.data.frame()
df_ghana <- df_ghana %>% select(rbp, dataid, country,  gw_cnum, !!(gee_vars)) %>% as.data.frame()
df_malawi <- df_malawi %>% select(rbp, dataid, country,  gw_cnum, !!(gee_vars)) %>% as.data.frame()
df_sierraleone <- df_sierraleone %>% select(rbp, dataid, country,  gw_cnum, !!(gee_vars)) %>% as.data.frame()
df <- df %>% select(rbp, dataid, country,  gw_cnum, !!(gee_vars)) %>% as.data.frame()

Xvars = gee_vars


slmodel=slmod
CVsetting=FALSE
Vfolds=4
id="gw_cnum"

d=df
outcome="rbp"
population="women"
Xvars=Xvars
id="country"
folds=Vfolds
CV=CVsetting
sl=slmodel


res_ghana=res_gambia= res_malawi= res_sierraleone= resfull=NULL

resfull=try(DHS_SL(d=df, outcome="rbp", population="women", Xvars=Xvars, id="country", folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(resfull, file = here("results/models/res_GW_transportability_SL_women_vitA_gee_only.rds"))


res_ghana=try(DHS_SL(d=df_ghana, outcome="rbp", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_ghana, file = here("results/models/res_GW_Ghana_SL_women_vitA_gee_only.rds"))

res_gambia=try(DHS_SL(d=df_gambia, outcome="rbp", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_ghana, file = here("results/models/res_GW_Gambia_SL_women_vitA_gee_only.rds"))

res_malawi=try(DHS_SL(d=df_malawi, outcome="rbp", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_malawi, file = here("results/models/res_GW_Malawi_SL_women_vitA_gee_only.rds"))

res_sierraleone=try(DHS_SL(d=df_sierraleone, outcome="rbp", population="women", Xvars=Xvars, id=id, folds=Vfolds, CV=CVsetting, sl=slmodel))
saveRDS(res_sierraleone, file = here("results/models/res_GW_SierraLeone_SL_women_vitA_gee_only.rds"))






