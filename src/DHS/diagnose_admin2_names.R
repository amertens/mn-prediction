# Diagnostic: Compare DHS admin-2 names vs survey admin-2 names
# Run this in RStudio to see unmatched names
library(here)

cat("\n========== GHANA ==========\n")
sp <- readRDS(here("data/DHS/clean/Ghana_2014_dhs_admin2_wide.rds"))
cu <- readRDS(here("data/DHS/clean/Ghana_2014_dhs_custom_admin2_wide.rds"))
survey <- readRDS(here("data/IPD/Ghana/Ghana_GMS_cleaned.rds"))

dhs_names <- sort(unique(c(sp$Admin2, cu$Admin2)))
survey_names <- sort(unique(survey$Admin2))

cat(sprintf("DHS admin2: %d unique names\n", length(dhs_names)))
cat(sprintf("Survey admin2: %d unique names\n", length(survey_names)))

# Exact case-insensitive matches
dhs_clean <- tolower(trimws(dhs_names))
survey_clean <- tolower(trimws(survey_names))
exact_match <- sum(dhs_clean %in% survey_clean)
cat(sprintf("Exact case-insensitive matches: %d\n\n", exact_match))

unmatched_dhs <- dhs_names[!dhs_clean %in% survey_clean]
unmatched_survey <- survey_names[!survey_clean %in% dhs_clean]

cat(sprintf("DHS names NOT in survey (%d):\n", length(unmatched_dhs)))
print(sort(unmatched_dhs))

cat(sprintf("\nSurvey names NOT in DHS (%d):\n", length(unmatched_survey)))
print(sort(unmatched_survey))


cat("\n\n========== MALAWI ==========\n")
sp2 <- readRDS(here("data/DHS/clean/Malawi_2015_dhs_admin2_wide.rds"))
cu2 <- readRDS(here("data/DHS/clean/Malawi_2015_dhs_custom_admin2_wide.rds"))
survey2 <- readRDS(here("data/IPD/Malawi/clean_malawi_mn_data.RDS"))

dhs_names2 <- sort(unique(c(sp2$Admin2, cu2$Admin2)))
survey_names2 <- sort(unique(survey2$Admin2))

cat(sprintf("DHS admin2: %d unique names\n", length(dhs_names2)))
cat(sprintf("Survey admin2: %d unique names\n", length(survey_names2)))

dhs_clean2 <- tolower(trimws(dhs_names2))
survey_clean2 <- tolower(trimws(survey_names2))
exact_match2 <- sum(dhs_clean2 %in% survey_clean2)
cat(sprintf("Exact case-insensitive matches: %d\n\n", exact_match2))

unmatched_dhs2 <- dhs_names2[!dhs_clean2 %in% survey_clean2]
unmatched_survey2 <- survey_names2[!survey_clean2 %in% dhs_clean2]

cat(sprintf("DHS names NOT in survey (%d):\n", length(unmatched_dhs2)))
print(sort(unmatched_dhs2))

cat(sprintf("\nSurvey names NOT in DHS (%d):\n", length(unmatched_survey2)))
print(sort(unmatched_survey2))
