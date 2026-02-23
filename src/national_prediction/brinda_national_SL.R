

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


#to do:
#-add vmnis as outcomes
#merge in metadata file


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

data <- read_dta("C:/Users/andre/OneDrive/Documents/mn-proxies/data/national/combined_dataset.dta")
head(data)
colnames(data)[1] <- "country"

brinda <- readRDS(here("data","clean_brinda_data.rds"))

head(brinda)

brinda_df <- rbindlist(brinda, fill=TRUE)
head(brinda_df)

#get mean from mean (95% CI)
brinda_df$prev_ferr <- as.numeric(str_split(brinda_df$prev_ferr, " \\(", simplify = T)[,1])
brinda_df$prev_stfr <- as.numeric(str_split(brinda_df$prev_stfr, " \\(", simplify = T)[,1])
brinda_df$prev_rol <- as.numeric(str_split(brinda_df$prev_rol, " \\(", simplify = T)[,1])
brinda_df$prev_rbp <- as.numeric(str_split(brinda_df$prev_rbp, " \\(", simplify = T)[,1])
brinda_df$prev_zinc <- as.numeric(str_split(brinda_df$prev_zinc, " \\(", simplify = T)[,1])


head(brinda_df)
d <- left_join(brinda_df,data,by=c("country","year"))
head(d[,1:30])


library(SuperLearner)
Xvars_full =c("year","age", colnames(d)[9:ncol(d)])
saveRDS(Xvars_full, file=here("metadata","BRINDA_national_predictor_variables.RDS"))
# outcome="prev_ferr"
# Xvars_full=Xvars_full
# id="country"
# folds=2
# CV=FALSE
# sl=slfull

run_national_sl <- function(d, outcome, Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull){

  df <- d %>% select(!!(outcome), all_of(Xvars_full), !!(id)) %>%
    mutate(id=factor(!!(id))) %>%
    as.data.frame()
  df <- df[!is.na(df[[outcome]]), ]
  df <- as.data.frame(df)
  cols_to_drop = colnames(df)[nzv(df)]
  df <- df[,!(colnames(df) %in% cols_to_drop)]
  head(df)

  Xvars_full = Xvars_full[Xvars_full %in% colnames(df)]


  res=try(DHS_SL(d=df, outcome=outcome, Xvars=Xvars_full, id=id, folds=folds, CV=CV, sl=sl))
  return(res)
}




res_prev_ferr=try(run_national_sl(d=d, outcome="prev_ferr", Xvars_full=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_ferr

res_prev_stfr=try(run_national_sl(d=d, outcome="prev_stfr", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_stfr

# res_prev_rol=try(run_national_sl(d=d, outcome="prev_rol", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=sl_simple))
# res_prev_rol

res_prev_rbp=try(run_national_sl(d=d,outcome="prev_rbp", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_rbp

res_prev_zinc=try(run_national_sl(d=d, outcome="prev_zinc", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_zinc

save(list=ls(pattern="res_"), file=here("results","BRINDA_national_SL_results.RData"))


#function to calculate model performances
perf_res= bind_rows(
  evaluate_sl_performance(res=res_prev_ferr, outcome="prev_ferr"),
  evaluate_sl_performance(res=res_prev_stfr, outcome="prev_stfr"),
  #evaluate_sl_performance(res=res_prev_rol, outcome="prev_rol"),
  evaluate_sl_performance(res=res_prev_rbp, outcome="prev_rbp"),
  evaluate_sl_performance(res=res_prev_zinc, outcome="prev_zinc"))
saveRDS(perf_res, file=here("results","BRINDA_national_SL_performance.RDS"))


