

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
library(terra)
library(caret)
source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

#load data
gps <- read.csv(here("data", "IPD", "Gambia", "MNS_Coordinates.csv")) %>%
  filter(!is.na(Latitude)) %>%
  mutate(
    latitude = as.numeric(Latitude),
    longitude = as.numeric(Longitude)
  ) %>% select(`MICS.Cluster_no`, latitude, longitude,  LGA, District_code, Ward_code) %>%
  rename(MICS_Cluster_Number=`MICS.Cluster_no`)
head(gps)

#save cleaned gps data
write.csv(gps, file=here("data", "IPD", "Gambia", "Gambia_GMS_GPS_cleaned.csv"))


df_women <- read_dta(here("data", "IPD", "Gambia", "Gambia_Woman data_Davis-Berkeley.dta")) %>% mutate(child_flag=0) %>%
  rename(MICS_Cluster_Number=wm1) %>%
  #drop empty rows
  filter(!is.na(cnum)) %>% rename(svy_weight=wmweight)
df_children <- read_dta(here("data", "IPD", "Gambia", "Gambia_Child data_Davis-Berkeley.dta")) %>% mutate(child_flag=1) %>%
  #drop empty rows
  filter(!is.na(cnum)) %>% rename(svy_weight=chweight)

gambia_variables <- makeVlist(df_women)
gambia_variables[grep("weight", gambia_variables$label),]

gambia_variables <- makeVlist(df_children)
gambia_variables[grep("weight", gambia_variables$label),]

table(is.na(df_women$cnum))
table(is.na(df_women$MICS_Cluster_Number))
table(is.na(df_children$MICS_Cluster_Number))
df_children$MICS_Household_Number

head(df_women)
head(df_children)

summary(df_women$cRBPAdjBrinda)
summary(df_children$cRBPAdjBrinda)

colnames(df_women)




#NOTE! Do this better

#merge women and children
df <- bind_rows(df_women, df_children)
table(is.na(df$cnum))
table(is.na(df$MICS_Cluster_Number))

#fill in missing mics number
mics_id <- df %>%
  select(cnum, MICS_Cluster_Number) %>%
  group_by(cnum) %>%
  summarise( MICS_Cluster_Number=first(MICS_Cluster_Number[!is.na(MICS_Cluster_Number)]))
df <-  df %>% subset(., select = -MICS_Cluster_Number)
df <- left_join(df, mics_id, by = "cnum")

df_cnum=unique(df$MICS_Cluster_Number)
gps_cnum=unique(gps$MICS_Cluster_Number)
df_cnum[!(df_cnum%in%gps_cnum)]
gps_cnum[!(gps_cnum%in%df_cnum)]

head(gps)

df <- left_join(df, gps, by = "MICS_Cluster_Number")
summary(df$latitude)
table(is.na(df$latitude))
table(df$child_flag, is.na(df$latitude))

#temp drop those missing locations, but debug later
df <- df %>% filter(!is.na(latitude) & !is.na(longitude))

#check and drop near zero variance columns
nzv_vars <- colnames(df[, nzv(df)])
dput(nzv_vars)


#drop unneeded variables
df <- df %>% subset(., select = -c(wm7, wm9, nkaq2c, nkaq2x, nkaq2y, nkaq3d,
                                   nkaq3x, nkaq3y, ib3, ff2, ff4, ff5, ff6, wdd1,
                                   pa4, pa7, pa7a, pa8, pa8b, npw2, rf1, rs1, dupes_bi1,
                                   dupes_wID, Note, merge_wweights, WM6Y, WM9, WM17,
                                   TA1, TA6, TA8A, TA8B, TA8C, TA8D, TA8X, TA8NR,
                                   TA9, TA10, TA12A, TA12B, TA12C, TA12X, TA12NR,
                                   TA13, TA14, windex10u, wEverBF, wSmokes, wDrinks,
                                   wMICSresult, merge_mics_women, HH5Y, HH10, HH12, HH46,
                                   HH14, HH17, HH33, HH39, HH44, hhaux, hhfin, HC1A,
                                   HC1B, HC7A, HC10F, HC10G, HC12, HC18G, EU3, EU6,
                                   EU7, IR2B, IR2C, IR2X, IR2Z, WS2, WS3, WS10A,
                                   WS10B, WS10D, WS10E, WS10F, WS10X, WS10Z, WS10NR,
                                   HW7B, HW7C, HW7NR, WQ5Y, WQ8, WQ31, WQ11, WQ12,
                                   WQ13, WQ15A, WQ15B, WQ15C, WQ15D, WQ15E, WQ15F,
                                   WQ15X, WQ15Z, WQ15NR, WQ18, WQ19, WQ21, WQ24Y,
                                   WQ29, wHH56, test, hReligion, hNationality, machine_D,
                                   tech_2, stove_97, light_13, water_11, merge_mics_hh,
                                   fer_flag, stfr_flag, rbp_flag, agp_flag, PLATE__, POS_NEG,
                                   Notes, wAnemPrev_Other, wFortOil_PrevBlind, wFortOil_RedMort,
                                   wFortOil_PrevVitDef, wFortOil_Other, wFortOil_DK, wEnrichi,
                                   wFortFlour, wFF_GivesBlood, wFF_Anemia, wFF_Mort, wFF_ImpHealth,
                                   wFF_School, wFF_Work, wFF_Other, wFF_DK, wWorkVigorousMinsPerWeek,
                                   wWalkWork, wMalariaYN,
                                   wweightconstant, StratumNatl, wBMISev, ch3m, ch3y,
                                   ch3t, ch3h, ch3mm, ch13, ch13h, ch13m, ch14, ch9,
                                   ch10, an3, an5, an7, bs6y, ref1, ref2, blood,
                                   C6, dupes_ch0b, Pos, C9, dupes_cID, cinsurance, cLRI,
                                   cEverBF, cConsInfantfomula, cYogurt, cDair2, cMea1,
                                   cMea2, merge_mics_child, hh5y, merge_gmns_hh, cBirthDateGMNS,
                                   cVASCardRecall, cAnemiaSev, cMalariaYN, cOverweightYN,
                                   cSAM_MUAC, `_logcrpcoeffRBP`, `_logagpcoeffRBP`, cweightconstant,
                                   triskin, subskin, `_zts`, `_zss`, `_fwfl`, `_flen`, `_fwei`,
                                   `_fbmi`, `_fhc`, `_fac`, `_fts`, `_fss`, merge_igrowup, cWasteSev,
                                   cUnderSev,  LGA, District_code,
                                   Ward_code))


head(df)



df$year = year(df$year)
df$month = month(as.numeric(df$month))



# #convert numeric variables
# df$B12_pmol_L<- as.numeric(df$B12_pmol_L)
# df$Folate_nmol_L <- as.numeric(df$Folate_nmol_L)
#
# #rename variables
# df <- df %>% rename( hhid=cachnc, hh_label=cachln, childid=caclnr, child_label=cacln, agem=caam1, child_weight=cbacwb, diar14d=cbbcdo, fever14d=cbbcfv)
# df$date

colnames(df)[!(colnames(df) %in% c("longitude","latitude"))] <- paste0("gw_", colnames(df)[!(colnames(df) %in% c("longitude","latitude"))])


saveRDS(df, file=here("data", "IPD", "Gambia", "Gambia_GMS_cleaned.rds"))


