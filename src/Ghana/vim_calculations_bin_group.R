


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

metadata= readRDS(here("metadata/variable_categories.rds"))
res = readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate_full.rds"))
dhs_indicators = read.csv( here("data/DHS/clean/dhs_indicators_metadata.csv"))
dhs_var_groups = unique(dhs_indicators$Level1)


sl_fit=res$sl_fit
eval_fun = loss_loglik_binomial
importance.metric="ratio"
n_vars=20





get_vim <- function(res, loss_fun=loss_loglik_binomial){

  X = res$task$nodes$covariates

  #get DHS variables
  dhs_var_list = list()
  for(i in 1:length(dhs_var_groups)){
    dhs_var_group=dhs_var_groups[i]
    vars = dhs_indicators$IndicatorId[dhs_indicators$Level1==dhs_var_group]

    temp = NULL

    for(j in vars){

      temp =c(temp, X[grepl(j, X)])

    }
    dhs_var_list[[i]]=temp
    names(dhs_var_list)[i] = paste0("DHS: ", dhs_var_group)
  }

  covariate_groups=c(dhs_var_list,

                     list(#"DHS" = metadata$dhs_vars[metadata$dhs_vars %in% X],
                       "MICS WASH" = metadata$mics_vars[metadata$mics_vars %in% X],
                       "IHME Morbidity" = metadata$ihme_vars[metadata$ihme_vars %in% X],
                       "GEE Climate" = metadata$gee_vars[metadata$gee_vars %in% X],
                       "LSMS SES" = metadata$lsms_vars[metadata$lsms_vars %in% X],
                       "Malaria" = metadata$map_vars[metadata$map_vars %in% X],
                       "Food prices" = metadata$wfp_vars[metadata$wfp_vars %in% X],
                       "FluNet" = metadata$flunet_vars[metadata$flunet_vars %in% X],
                       "Missingness"= X[grepl("missing_", X)]))


  set.seed(12345)

  res=try(calc_mn_importance(res$sl_fit,  eval_fun = loss_fun, importance.metric="difference",
                             covariate_groups=covariate_groups))
  return(res)

}



#res1 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA_full.rds"))
res3 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA.rds"))
#res4 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA_full.rds"))
res6 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA.rds"))
#res7 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12_full.rds"))
res9 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12.rds"))
#res10 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate_full.rds"))
res12 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate.rds"))
#res13 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron_full.rds"))
res15 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron.rds"))
#res16 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron_full.rds"))
res18 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron.rds"))


# child_vitA_full_vim =get_vim(res=res1)
# saveRDS(child_vitA_full_vim, file = here("results/child_vitA_def_varimp_group.rds"))

child_vitA_vim =get_vim(res=res3)
saveRDS(child_vitA_vim, file = here("results/child_vitA_def_proxy_varimp_group.rds"))


# child_iron_full_vim =get_vim(res=res13)
# saveRDS(child_iron_full_vim, file = here("results/child_iron_def_varimp_group.rds"))

child_iron_vim =get_vim(res=res15)
saveRDS(child_iron_vim , file = here("results/child_iron_def_proxy_varimp_group.rds"))


# women_vitA_full_vim =get_vim(res=res4)
# saveRDS(women_vitA_full_vim, file = here("results/women_vitA_def_varimp_group.rds"))

women_vitA_vim =get_vim(res=res6)
saveRDS(women_vitA_vim, file = here("results/women_vitA_def_proxy_varimp_group.rds"))

# women_b12_full_vim =get_vim(res=res7)
# saveRDS(women_b12_full_vim, file = here("results/women_b12_def_varimp_group.rds"))

women_b12_vim =get_vim(res=res9)
saveRDS(women_b12_vim, file = here("results/women_b12_def_proxy_varimp_group.rds"))

# women_folate_full_vim =get_vim(res = res10)
# saveRDS(women_folate_full_vim, file = here("results/women_folate_def_varimp_group.rds"))

women_folate_vim =get_vim(res = res12)
saveRDS(women_folate_vim, file = here("results/women_folate_def_proxy_varimp_group.rds"))


# women_iron_full_vim =get_vim(res=res16)
# saveRDS(women_b12_full_vim, file = here("results/women_iron_def_varimp_group.rds"))

women_iron_vim =get_vim(res=res18)
saveRDS(women_iron_vim, file = here("results/women_iron_def_proxy_varimp_group.rds"))


