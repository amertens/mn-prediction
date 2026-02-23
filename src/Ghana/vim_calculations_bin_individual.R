


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
options(future.globals.maxSize = 5 * 1024^3)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

res_child_vitA <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA.rds"))
res_women_vitA <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA.rds"))
res_women_b12 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12.rds"))
res_women_folate <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate.rds"))
res_child_iron <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron.rds"))
res_mom_iron <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron.rds"))


start_time <- Sys.time()
res_vim_diff_child_vitA=try(calc_mn_importance(res_child_vitA$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
end_time <- Sys.time()
end_time - start_time

saveRDS(res_vim_diff_child_vitA, file= here("results/vim_child_vitA_def_proxy_individ.rds"))

res_vim_diff_child_iron=try(calc_mn_importance(res_child_iron$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
saveRDS(res_vim_diff_child_iron, file= here("results/vim_child_iron_def_proxy_individ.rds"))

res_vim_diff_women_vitA=try(calc_mn_importance(res_women_vitA$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
saveRDS(res_vim_diff_women_vitA, file= here("results/vim_women_vitA_def_proxy_individ.rds"))

res_vim_diff_women_iron=try(calc_mn_importance(res_women_iron$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
saveRDS(res_vim_diff_women_iron, file= here("results/vim_women_iron_def_proxy_individ.rds"))

res_vim_diff_women_folate=try(calc_mn_importance(res_women_folate$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
saveRDS(res_vim_diff_women_folate, file= here("results/vim_women_folate_def_proxy_individ.rds"))

res_vim_diff_women_b12=try(calc_mn_importance(res_women_b12$sl_fit,  eval_fun = loss_loglik_binomial, importance.metric="difference", covariate_groups=NULL))
saveRDS(res_vim_diff_women_b12, file= here("results/vim_women_b12_def_proxy_individ.rds"))

