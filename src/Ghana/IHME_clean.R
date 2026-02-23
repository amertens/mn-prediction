
# Load necessary libraries
rm(list=ls())
library(tidyverse)
library(readr)
library(fuzzyjoin)
library(here)

# # Load raster datasets
# oncho <- read_csv(here("data/IHME/Africa and Yemen Onchocerciasis Prevalence Geospatial Estimates 2000-2018/Data [CSV]/IHME_AFRICA_ONCHO_2000_2018_PREV_ADMIN_2_Y2020M08D18.CSV"))
# hiv <- read_csv(here("data/IHME/Africa HIV Prevalence Geospatial Estimates 2000-2017/Data [CSV]/IHME_AFRICA_HIV_2000_2017_HIV_PREV_ADMIN_2_Y2019M03D15.CSV"))
# gps <- read.csv(here("data", "IPD", "Ghana", "Ghana_GMS_GPS_cleaned.csv"))
#
# # Filter for Ghana 2017
# oncho_ghana <- oncho %>% filter(ADM0_NAME == "Ghana", year == 2017)
# hiv_ghana <- hiv %>% filter(ADM0_NAME == "Ghana", year == 2017)



ihme <- read.csv(here("data/IHME/IHME_Ghana_data.csv")) %>% filter(adm0_name=="Ghana" & year==2017) %>%
  filter(adm2_name!="", source!="africa hiv prevalence geospatial estimates 2000-2017",
         source!="lmic dbm geospatial estimates 2000-2017", !is.na(mean)) #%>%
#select(adm2_name, age_group_name, age_group_id, sex, measure, metric, mean)
head(ihme)

ihme$measure[ihme$source=="global malaria geospatial estimates 2000-2019" & ihme$measure=="Incidence"] = "malaria incidence"
ihme$measure[ihme$source=="global malaria geospatial estimates 2000-2019" & ihme$measure=="Prevalence"] = "malaria prevalence"

ihme$measure[ihme$source=="lmic under-5 diarrhea geospatial estimates 2000-2019" & ihme$measure=="Incidence"] = "u5 diarrhea incidence"
ihme$measure[ihme$source=="lmic under-5 diarrhea geospatial estimates 2000-2019" & ihme$measure=="Mortality"] = "u5 diarrhea mortality"
ihme$measure[ihme$source=="lmic under-5 diarrhea geospatial estimates 2000-2019" & ihme$measure=="Prevalence"] = "u5 diarrhea prevalence"

ihme$measure[ihme$source=="africa and yemen onchocerciasis prevalence geospatial estimates 2000-2018" & ihme$measure=="Prevalence"] = "onchocerciasis prevalence"
ihme$measure[ihme$source=="global lymphatic filariasis prevalence geospatial estimates 2000-2018" & ihme$measure=="Prevalence"] = "LF prevalence"


#ihme <- ihme %>% filter(measure %in% c("Mortality","Prevalence","Incidence"))


table(ihme$sex)

table(ihme$measure, ihme$sex)
table(ihme$source, ihme$measure)

temp = ihme %>% dplyr::summarise(n = dplyr::n(), .by = c(adm2_name, age_group_name, age_group_id, sex, measure, metric)) |>
  dplyr::filter(n > 1L)
temp

temp = ihme %>% filter( adm2_name=="Adansi North", age_group_name=="All Ages",sex=="Both",
                        measure=="Prevalence",metric=="Rate")
temp
# Create combined column names from age_group_name, sex, and measure
df_wide_tidyr <- ihme %>%
  # Create a unique identifier for each combination
  unite("variable", age_group_name, age_group_id, sex, measure, metric, sep = "_", remove = FALSE) %>%
  # Convert to wide format
  pivot_wider(
    id_cols = adm2_name,
    names_from = variable,
    values_from = mean,
    values_fill = NA  # Fill missing combinations with NA
  ) %>% janitor::clean_names(.)
colnames(df_wide_tidyr) <- paste0("ihme_",colnames(df_wide_tidyr))
colnames(df_wide_tidyr) <- gsub("_x","_",colnames(df_wide_tidyr))
head(df_wide_tidyr)



# Set file



# Save result
write_csv(df_wide_tidyr, here("data/IHME/ghana_2017_merged_IHME_data.csv"))


