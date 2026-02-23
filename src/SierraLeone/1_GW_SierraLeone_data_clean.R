

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


gps <- readxl::read_xlsx(here("data", "IPD", "Sierra Leone", "Sierra Leone_Cluster GPS.xlsx")) %>%
  rename(cnum=CLUSTER_NU , latitude=Latitude, longitude=Longitude)
head(gps)

#save cleaned gps data
write.csv(gps, file=here("data", "IPD", "Sierra Leone", "Sierra Leone_GMS_GPS_cleaned.csv"))


df_women <- read_dta(here("data", "IPD", "Sierra Leone", "Sierra Leone_Woman data_Davis-Berkeley.dta")) %>% rename(svy_weight=wStatWtComb)
df_children <- read_dta(here("data", "IPD", "Sierra Leone", "Sierra Leone_Child data_Davis-Berkeley.dta")) %>% rename(svy_weight= cStatWt)

SL_variables <- makeVlist(df_women)
SL_variables[grep("weight", SL_variables$label),]

SL_variables <- makeVlist(df_children)
SL_variables[grep("weight", SL_variables$label),]
SL_variables[SL_variables$name=="cFeDefAdj",]
attributes(df_children$cFeDefAdj)

summary(df_women$cRBPAdjBrinda)
summary(df_children$cRBPAdjBrinda)

# #get ID and date dateset from women's dataset
# summary(df_women$intymd)
# table(is.na(df_women$intymd))
# unique(df_children$cachnc)
# unique(df_women$b_hhn)
# unique(df_children$cachln)
# unique(df_women$b_hhln)


unique(df_women$hClusterNumb)
unique(df_children$hClusterNumb)
gps$cnum

df <- bind_rows(df_women %>% mutate(child_flag=0),
                df_children %>% mutate(child_flag=1)) %>%
  rename(cnum=hClusterNumb) #double check this
head(df)
df$indivID

df <- left_join(df, gps, by = "cnum")
summary(df$latitude)
head(df)



#check and drop near zero variance columns
colnames(df[, nzv(df)])
dput(colnames(df[, nzv(df)]))


#drop unneeded variables
df <- df %>% subset(., select = -c(wFinalResult, wOtherResult, wStartPermission, wReadLevelOth,
                                   wWorkOth, wSmokeCigs, wSmokeCigsNumb, wPregNow, wMeatOrgan,
                                   wInsects, wFoodsOth, wIodIntell, wIodOth,
                                   wVitAOth, wPregLab, wUrineCup,
                                   hNumbNPWData_sum, hDOIMonth, hDOIYear, hFinalResult,
                                   hOtherResult, hStartNow, hReligionOth, hLanguageOth,
                                   hFloorOth, hRoofOth, hExtWallsOth, hCookFuelOth, hHasLandline,
                                   hHasCarTruck, hHasCanoe, hHasMotorboat, hHasWheelbarrow,
                                   hHasSprayer, hHasERiceCutter, hCattleYN, hHorsesYN, hHorsesNumb,
                                   hRabbitsYN, hRabbitsNumb, hPigsYN, hRodentsYN, hBirdsYN,
                                   hBeesYN, hBeesNumb, hAnimOth, hDrinkSourceOth, hWaterFilter,
                                   hWaterSolar, hWaterTreatOth, hWashSourceOth, hSanitOth,
                                   hHandWashLiqSoap, hHandWashAsh, hAllowSalt, hGiveSalt,
                                   hBuyOilQuantUnit, hWheatTypeOth, hBreadTypeOth, HRoorGrp2,
                                   HCookFuelGrp2, wFoodGrp1Grains, wFoodGrp5OrgMeat, UICatPreg,
                                   UICGrpPreg, cFinalResult, cOtherResult, cStartPermission, cFatherAlive,
                                   cFastBreathOth, cBFEver, cBFYest, cVegOrange, cMeatOrgan,
                                   cMeat, cInsects, cRUTF, cLabConsent, hHRostLine01,
                                   hHRostLine02, hHRostLine03, hHRostLine04, hHRostLine05,
                                   hHRostLine06, hHRostLine07, hHRostLine08, hHRostLine09,
                                   hHRostLine10, hHRostLine11, hHRostLine12, hHRostLine13,
                                   hHRostLine14, hHRostLine15, hHRostLine16, hHRostLine17,
                                   hHRostLine18, hHRostLine19, hHRostLine20, hHRostLine21,
                                   hHRostLine22, hHRostLine23, hHRostLine24, hHRostLine25,
                                   hHRostLine26, hHRostLine27, hHRostLine28, hHRostLine29,
                                   hHRostLine30, hHRostSex22, hHRostSex24, hHRostSex25,
                                   hHRostSex26, hHRostSex27, hHRostSex28, hHRostSex29, hHRostSex30,
                                   hHRostDOBMo23, hHRostDOBMo24, hHRostDOBMo25, hHRostDOBMo26,
                                   hHRostDOBMo27, hHRostDOBMo28, hHRostDOBMo29, hHRostDOBMo30,
                                   hHRostDOBYr24, hHRostDOBYr25, hHRostDOBYr26, hHRostDOBYr27,
                                   hHRostDOBYr28, hHRostDOBYr29, hHRostDOBYr30, hHRostAge24,
                                   hHRostAge25, hHRostAge26, hHRostAge27, hHRostAge28, hHRostAge29,
                                   hHRostAge30, hHRostWomanYN01, hHRostWomanYN02, hHRostWomanYN03,
                                   hHRostWomanYN04, hHRostWomanYN05, hHRostWomanYN06, hHRostWomanYN07,
                                   hHRostWomanYN08, hHRostWomanYN09, hHRostWomanYN10, hHRostWomanYN11,
                                   hHRostWomanYN12, hHRostWomanYN13, hHRostWomanYN14, hHRostWomanYN15,
                                   hHRostWomanYN16, hHRostWomanYN17, hHRostWomanYN18, hHRostWomanYN19,
                                   hHRostWomanYN20, hHRostWomanYN21, hHRostWomanYN22, hHRostWomanYN23,
                                   hHRostWomanYN24, hHRostWomanYN25, hHRostWomanYN26, hHRostWomanYN27,
                                   hHRostWomanYN28, hHRostWomanYN29, hHRostWomanYN30, hHRostPregYN01,
                                   hHRostPregYN02, hHRostPregYN04, hHRostPregYN05, hHRostPregYN06,
                                   hHRostPregYN07, hHRostPregYN08, hHRostPregYN09, hHRostPregYN10,
                                   hHRostPregYN11, hHRostPregYN12, hHRostPregYN13, hHRostPregYN14,
                                   hHRostPregYN15, hHRostPregYN16, hHRostPregYN17, hHRostPregYN18,
                                   hHRostPregYN19, hHRostPregYN20, hHRostPregYN21, hHRostPregYN22,
                                   hHRostPregYN23, hHRostPregYN24, hHRostPregYN25, hHRostPregYN26,
                                   hHRostPregYN27, hHRostPregYN28, hHRostPregYN29, hHRostPregYN30,
                                   hHRostChildYN01, hHRostChildYN02, hHRostChildYN03, hHRostChildYN04,
                                   hHRostChildYN05, hHRostChildYN06, hHRostChildYN07, hHRostChildYN08,
                                   hHRostChildYN09, hHRostChildYN10, hHRostChildYN11, hHRostChildYN12,
                                   hHRostChildYN13, hHRostChildYN14, hHRostChildYN15, hHRostChildYN16,
                                   hHRostChildYN17, hHRostChildYN18, hHRostChildYN19, hHRostChildYN20,
                                   hHRostChildYN21, hHRostChildYN22, hHRostChildYN23, hHRostChildYN24,
                                   hHRostChildYN25, hHRostChildYN26, hHRostChildYN27, hHRostChildYN28,
                                   hHRostChildYN29, hHRostChildYN30, hHRostCaretakerID01,
                                   hHRostCaretakerID20, hHRostCaretakerID21, hHRostCaretakerID22,
                                   hHRostCaretakerID23, hHRostCaretakerID24, hHRostCaretakerID25,
                                   hHRostCaretakerID26, hHRostCaretakerID27, hHRostCaretakerID28,
                                   hHRostCaretakerID29, hHRostCaretakerID30, hNumChildren_sum,
                                   hNumbChildrenData_sum, PrimaryLast, cAnemiaYN, LBW, BFExclusive,
                                   cAgeGrpLT6, cFeDefUNAdj, cFeRichFood ))


#Looks like no time information in Sierra Leone
# head(df)
#
# df$date= parse_date_time(as.character(df$start),
#                            orders = c("a b d H:M:S z Y", "a b d H:M:S Y", "Y-m-d H:M:S"),
#                            tz = "UTC")
#
# df$year = year(df$date)
# df$month = month(df$date)
# df$week = week(df$date)
#
# df  <- df %>% subset(., select = -c(start, end, deviceid, tn, intymd, amon, ayr, birthymd))



colnames(df)[!(colnames(df) %in% c("longitude","latitude"))] <- paste0("gw_", colnames(df)[!(colnames(df) %in% c("longitude","latitude"))])


saveRDS(df, file=here("data", "IPD", "Sierra Leone", "SierraLeone_GMS_cleaned.rds"))


