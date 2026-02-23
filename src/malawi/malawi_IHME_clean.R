
# Load necessary libraries
rm(list=ls())
library(tidyverse)
library(readr)
library(fuzzyjoin)
library(here)
library(haven)

# # Load raster datasets
ihme <- read_dta(here("data/IHME/IHME_Malawi_data.dta")) %>% filter(adm0_name=="Malawi" & year==2015) %>%
  filter(adm2_name!="", source!="africa hiv prevalence geospatial estimates 2000-2017",
         source!="lmic dbm geospatial estimates 2000-2017", !is.na(mean)) %>%
  mutate(mean=as.numeric(mean)) %>% as.data.frame()
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
table(ihme$source[ihme$measure=="Prevalence"])

temp = ihme %>% dplyr::summarise(n = dplyr::n(), .by = c(adm2_name, age_group_name, age_group_id, sex, measure, metric)) |>
  dplyr::filter(n > 1L)
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
df_wide_tidyr$ihme_15_49_years_24_males_number_of_uncircumcised_men_counts


# Set file



# Save result
write_csv(df_wide_tidyr, here("data/IHME/Malawi_2015_merged_IHME_data.csv"))


