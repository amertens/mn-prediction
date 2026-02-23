


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
library(future)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

res = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_full.rds"))

plotdf=data.frame(Y=res$Y, yhat_full=res$yhat_full)
ggplot(aes(y=Y, x=yhat_full), data=plotdf) + geom_point()


res = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_gwPreds.rds"))


plotdf=data.frame(Y=res$Y, yhat_full=res$yhat_full)
ggplot(aes(y=Y, x=yhat_full), data=plotdf) + geom_point()



res = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA.rds"))

plotdf=data.frame(Y=res$Y, yhat_full=res$yhat_full)
ggplot(aes(y=Y, x=yhat_full), data=plotdf) + geom_point()
