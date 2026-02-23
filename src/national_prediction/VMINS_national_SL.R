

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

vmnis_full <- haven::read_dta(here("data","VMNIS","VMNISIndicator_long_format.dta"))
colnames(vmnis_full)
vmnis <- vmnis_full %>% rename(country=Country, year=Beginyear) %>%
  select(country, year, Indicator, Representativeness, Population, Mean, Geomean, Median, starts_with("Prevalence"), Dataadjustedfor, Surveymethodology, Methodofanalysis, Season, Agefrom, Ageto)
head(vmnis)


head(data[,1:30])
d <- left_join(vmnis,data,by=c("country","year"))
head(d[,1:30])

colnames(d)

library(SuperLearner)
Xvars_full =c("year", "Agefrom", "Ageto", "Dataadjustedfor", "Season", colnames(d)[32:ncol(d)])

table(d$Indicator)


# d=d_b12
# outcome="Mean"
# Xvars=Xvars_full
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


table(d$Indicator)

load(here("results","VMINS_national_SL_results.RData"))

#-------------------------------------------------------------------------------
# Zinc
#-------------------------------------------------------------------------------

d_zinc <- d %>% filter(Indicator=="Zinc (plasma or serum)")
head(d_zinc[,1:30])

res_mean_zinc=try(run_national_sl(d=d_zinc, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_mean_zinc

res_prev_zinc=try(run_national_sl(d=d_zinc, outcome="Prevalenceofdeficiency", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_zinc


#-------------------------------------------------------------------------------
# B12
#-------------------------------------------------------------------------------

d_b12 <- d %>% filter(Indicator=="Vitamin B12")
head(d_b12)
summary(d_b12$Mean)
summary(d_b12$Geomean )
summary(d_b12$Median)
summary(d_b12$Prevalenceofdeficiency)

#NOTE! SLfull fails- debug
res_mean_b12=try(run_national_sl(d=d_b12, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=sl_simple))
res_mean_b12

res_prev_b12=try(run_national_sl(d=d_b12, outcome="Prevalenceofdeficiency", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=sl_simple))
res_prev_b12

save(list=ls(pattern="res_"), file=here("results","VMINS_national_SL_results.RData"))


#-------------------------------------------------------------------------------
#  Retinol (plasma or serum)
#-------------------------------------------------------------------------------


d_retinol <- d %>% filter(Indicator=="Retinol (plasma or serum)")
head(d_retinol)
summary(d_retinol$Mean)
summary(d_retinol$Geomean )
summary(d_retinol$Median)
summary(d_retinol$Prevalenceofdeficiency)

res_mean_retinol=try(run_national_sl(d=d_retinol, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_mean_retinol

res_prev_retinol=try(run_national_sl(d=d_retinol, outcome="Prevalenceofdeficiency", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_retinol

save(list=ls(pattern="res_"), file=here("results","VMINS_national_SL_results.RData"))


#-------------------------------------------------------------------------------
#  25-Hydroxyvitamin D
#-------------------------------------------------------------------------------

d_vitD <- d %>% filter(Indicator=="25-Hydroxyvitamin D")
head(d_vitD)
summary(d_vitD$Mean)
summary(d_vitD$Geomean )
summary(d_vitD$Median)
summary(d_vitD$Prevalenceofdeficiency)

#NOTE! SLfull fails- debug
res_mean_vitD=try(run_national_sl(d=d_vitD, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=sl_simple))
res_mean_vitD

res_prev_vitD=try(run_national_sl(d=d_vitD, outcome="Prevalenceofdeficiency", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=sl_simple))
res_prev_vitD

save(list=ls(pattern="res_"), file=here("results","VMINS_national_SL_results.RData"))


#-------------------------------------------------------------------------------
#  Folate (plasma or serum)
#-------------------------------------------------------------------------------

d_folate <- d %>% filter(Indicator=="Folate (plasma or serum)")
head(d_folate)
summary(d_folate$Mean)
summary(d_folate$Geomean )
summary(d_folate$Median)
summary(d_folate$Prevalenceofdeficiency)

res_mean_folate=try(run_national_sl(d=d_folate, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_mean_folate

res_prev_folate=try(run_national_sl(d=d_folate, outcome="Prevalenceofdeficiency", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_prev_folate

save(list=ls(pattern="res_"), file=here("results","VMINS_national_SL_results.RData"))


#-------------------------------------------------------------------------------
#  Haemoglobin
#-------------------------------------------------------------------------------

d_haemoglobin <- d %>% filter(Indicator=="Haemoglobin")
head(d_haemoglobin)
summary(d_haemoglobin$Mean)
summary(d_haemoglobin$Geomean )
summary(d_haemoglobin$Median)
summary(d_haemoglobin$Prevalenceofdeficiency)

res_mean_haemoglobin=try(run_national_sl(d=d_haemoglobin, outcome="Mean", Xvars=Xvars_full, id="country", folds=2, CV=FALSE, sl=slfull))
res_mean_haemoglobin


save(list=ls(pattern="res_"), file=here("results","VMINS_national_SL_results.RData"))



#-------------------------------------------------------------------------------
#  Calc model performances
#-------------------------------------------------------------------------------


#function to calculate model performances
perf_res= bind_rows(
  evaluate_sl_performance(res=res_mean_zinc, outcome="mean_zinc"),
  evaluate_sl_performance(res=res_prev_zinc, outcome="prev_zinc"),
  evaluate_sl_performance(res=res_mean_b12, outcome="mean_b12"),
  evaluate_sl_performance(res=res_prev_b12, outcome="prev_b12"),
  evaluate_sl_performance(res=res_mean_retinol, outcome="mean_retinol"),
  evaluate_sl_performance(res=res_prev_retinol, outcome="prev_retinol"),
  evaluate_sl_performance(res=res_mean_vitD, outcome="mean_vitD"),
  evaluate_sl_performance(res=res_prev_vitD, outcome="prev_vitD"),
  evaluate_sl_performance(res=res_mean_folate, outcome="mean_folate"),
  evaluate_sl_performance(res=res_prev_folate, outcome="prev_folate")
  )


saveRDS(perf_res, file=here("results","VMNIS_national_SL_performance.RDS"))

