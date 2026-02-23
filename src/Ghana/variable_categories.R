

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


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#read clean and merged data
df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))
#dput(colnames(df))


gw_vars <- colnames(df)[grepl("gw_", colnames(df))]
dhs_vars <- colnames(df)[grepl("dhs_", colnames(df))]
mics_vars <- colnames(df)[grepl("mics_", colnames(df))]
ihme_vars <- colnames(df)[grepl("ihme_", colnames(df))]
lsms_vars <- colnames(df)[grepl("lsms_", colnames(df))]
map_vars <- colnames(df)[grepl("MAP_", colnames(df))]
wfp_vars <- colnames(df)[grepl("wfp_", colnames(df))]
flunet_vars <- colnames(df)[grepl("flunet_", colnames(df))]
gee_vars <- colnames(df)[grepl("gee_", colnames(df))]


variable_sources=list(
  gw_vars=gw_vars,
  dhs_vars=dhs_vars,
  mics_vars=mics_vars,
  ihme_vars=ihme_vars,
  lsms_vars=lsms_vars,
  map_vars=map_vars,
  wfp_vars=wfp_vars,
  flunet_vars=flunet_vars,
  gee_vars=gee_vars
)

metadata=list(variable_sources=variable_sources)

saveRDS(metadata, file=here("metadata/variable_categories.rds"))


