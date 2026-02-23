

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
library(washb)
options(future.globals.maxSize = 4 * 1024^3)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#read clean and merged data
dfull <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))

#vitamin A data
summary(dfull$r_crpagp_sf2)

mn_vars <- colnames(dfull)[grepl("gw_", colnames(dfull))]
dhs_vars <- colnames(dfull)[grepl("dhs_", colnames(dfull))]
mics_vars <- colnames(dfull)[grepl("mics_", colnames(dfull))]
ihme_vars <- colnames(dfull)[grepl("ihme_", colnames(dfull))]
lsms_vars <- colnames(dfull)[grepl("lsms_", colnames(dfull))]
map_vars <- colnames(dfull)[grepl("MAP_", colnames(dfull))]
wfp_vars <- colnames(dfull)[grepl("wfp_", colnames(dfull))]
flunet_vars <- colnames(dfull)[grepl("flunet_", colnames(dfull))]



Xvars = c("gw_month", mn_vars)


dfull <- dfull %>% select(gw_cRBPAdjBrinda, gw_r_crpagp_sf2, gw_cnum, !!(Xvars)) %>% as.data.frame() %>% filter(!is.na(gw_r_crpagp_sf2))



library(DataExplorer)

#temp
#df <- dfull[,1:200]
df <- dfull

# Comprehensive overview
create_report(df)
