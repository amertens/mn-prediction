
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
library(sf)
library(terra)
library(caret)
library(readxl)
library(survey)
library(srvyr)     # tidy interface to 'survey'
library(stringr)
library(purrr)
library(viridis)
library(scales)


adm2_prev = readRDS(here("data", "IPD", "Ghana", "Ghana_GMS_admin2_prevalence.rds"))
head(adm2_prev)



df <- readRDS(here("data/IPD/Ghana/Ghana_merged_dataset.rds"))
df$Admin2
df$flag <- 1

#subset to admin 2-level
df_admin_2 <- df %>% distinct(Admin2, .keep_all = TRUE)

dim(df_admin_2)
length(unique(df_admin_2$Admin2))

dim(df_admin_2)
dim(adm2_prev)
d <- left_join(adm2_prev, df_admin_2, by = "Admin2")

saveRDS(d, here("data/IPD/Ghana/Ghana_GMS_admin2_prevalence_merged_dataset.rds"))
