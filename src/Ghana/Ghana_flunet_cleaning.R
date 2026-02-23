# Ghana 2017 FluNet Data Analysis
# Subset and summarize weekly flu incidence data

# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(lubridate)
library(here)

# Read the FluNet data
flunet_data <- read_csv(here("data/FluNet/VIW_FNT.csv"))


# Subset data for Ghana 2017
ghana_2017 <- flunet_data %>%
  filter(COUNTRY_AREA_TERRITORY == "Ghana" & ISO_YEAR == 2017) %>%
  arrange(ISO_WEEK)

# Check the data
cat("Number of records for Ghana 2017:", nrow(ghana_2017), "\n")
cat("Week range:", min(ghana_2017$ISO_WEEK, na.rm = TRUE), "to", max(ghana_2017$ISO_WEEK, na.rm = TRUE), "\n")
cat("Unique locations:", unique(ghana_2017$COUNTRY_AREA_TERRITORY), "\n")

# Convert ISO week start date to Date format
ghana_2017$ISO_WEEKSTARTDATE <- as.Date(ghana_2017$ISO_WEEKSTARTDATE)

# Create summary statistics for weekly flu incidence (national level)
weekly_summary <- ghana_2017 %>%
  select(ISO_WEEK, ISO_WEEKSTARTDATE, SPEC_PROCESSED_NB, SPEC_RECEIVED_NB,
         INF_A, INF_B, INF_ALL, INF_NEGATIVE, ILI_ACTIVITY) %>%
  mutate(
    # Calculate positivity rates
    inf_a_rate = ifelse(SPEC_PROCESSED_NB > 0, INF_A / SPEC_PROCESSED_NB * 100, 0),
    inf_b_rate = ifelse(SPEC_PROCESSED_NB > 0, INF_B / SPEC_PROCESSED_NB * 100, 0),
    inf_all_rate = ifelse(SPEC_PROCESSED_NB > 0, INF_ALL / SPEC_PROCESSED_NB * 100, 0),
    # Handle missing INF_ALL values by summing INF_A and INF_B
    inf_all_calculated = ifelse(is.na(INF_ALL),
                                ifelse(is.na(INF_A), 0, INF_A) + ifelse(is.na(INF_B), 0, INF_B),
                                INF_ALL),
    inf_all_rate_calc = ifelse(SPEC_PROCESSED_NB > 0, inf_all_calculated / SPEC_PROCESSED_NB * 100, 0)
  ) %>%
  arrange(ISO_WEEK)

# Display the weekly summary
print("Weekly Flu Incidence Summary for Ghana 2017:")
print(weekly_summary)

# Overall summary statistics
overall_summary <- ghana_2017 %>%
  summarise(
    total_weeks = n(),
    total_specimens_received = sum(SPEC_RECEIVED_NB, na.rm = TRUE),
    total_specimens_processed = sum(SPEC_PROCESSED_NB, na.rm = TRUE),
    total_inf_a = sum(INF_A, na.rm = TRUE),
    total_inf_b = sum(INF_B, na.rm = TRUE),
    total_inf_all = sum(INF_ALL, na.rm = TRUE),
    total_inf_negative = sum(INF_NEGATIVE, na.rm = TRUE),
    avg_weekly_specimens = mean(SPEC_PROCESSED_NB, na.rm = TRUE),
    avg_weekly_inf_a = mean(INF_A, na.rm = TRUE),
    avg_weekly_inf_b = mean(INF_B, na.rm = TRUE),
    avg_weekly_inf_all = mean(INF_ALL, na.rm = TRUE),
    overall_positivity_rate = ifelse(total_specimens_processed > 0,
                                     total_inf_all / total_specimens_processed * 100, 0),
    inf_a_positivity = ifelse(total_specimens_processed > 0,
                              total_inf_a / total_specimens_processed * 100, 0),
    inf_b_positivity = ifelse(total_specimens_processed > 0,
                              total_inf_b / total_specimens_processed * 100, 0)
  )

print("\nOverall Summary Statistics for Ghana 2017:")
print(overall_summary)

# Create visualizations
# 1. Weekly flu incidence over time
p1 <- ggplot(weekly_summary, aes(x = ISO_WEEK)) +
  geom_line(aes(y = INF_A, color = "Influenza A"), size = 1) +
  geom_line(aes(y = INF_B, color = "Influenza B"), size = 1) +
  geom_line(aes(y = inf_all_calculated, color = "Total Influenza"), size = 1) +
  labs(title = "Weekly Flu Incidence - Ghana 2017",
       x = "ISO Week",
       y = "Number of Positive Cases",
       color = "Flu Type") +
  theme_minimal() +
  scale_color_manual(values = c("Influenza A" = "red",
                                "Influenza B" = "blue",
                                "Total Influenza" = "purple"))

# 2. Weekly positivity rates
p2 <- ggplot(weekly_summary, aes(x = ISO_WEEK)) +
  geom_line(aes(y = inf_a_rate, color = "Influenza A"), size = 1) +
  geom_line(aes(y = inf_b_rate, color = "Influenza B"), size = 1) +
  geom_line(aes(y = inf_all_rate_calc, color = "Total Influenza"), size = 1) +
  labs(title = "Weekly Flu Positivity Rates - Ghana 2017",
       x = "ISO Week",
       y = "Positivity Rate (%)",
       color = "Flu Type") +
  theme_minimal() +
  scale_color_manual(values = c("Influenza A" = "red",
                                "Influenza B" = "blue",
                                "Total Influenza" = "purple"))

# 3. Testing volume over time
p3 <- ggplot(weekly_summary, aes(x = ISO_WEEK)) +
  geom_line(aes(y = SPEC_PROCESSED_NB), color = "darkgreen", size = 1) +
  labs(title = "Weekly Testing Volume - Ghana 2017",
       x = "ISO Week",
       y = "Number of Specimens Processed") +
  theme_minimal()

# Display plots
print(p1)
print(p2)
print(p3)

# Note about subnational data
cat("\n=== SUBNATIONAL DATA AVAILABILITY ===\n")
cat("Based on the data structure, subnational flu surveillance data for Ghana 2017 is NOT available.\n")
cat("All records show 'Ghana' as the location, indicating national-level reporting only.\n")
cat("This is common for many countries in the FluNet system where only national surveillance is conducted.\n")


#clean time data
weekly_summary$year <- year(weekly_summary$ISO_WEEKSTARTDATE)


# Export results
colnames(weekly_summary)
weekly_summary <- weekly_summary %>% select(ISO_WEEK, inf_a_rate, inf_b_rate, inf_all_calculated) %>% rename(week=ISO_WEEK)
colnames(weekly_summary) <- paste0("flunet_", colnames(weekly_summary))
write.csv(weekly_summary, here("data/FluNet/Ghana/ghana_2017_weekly_flu_summary.csv"), row.names = FALSE)
write.csv(overall_summary,  here("data/FluNet/Ghana/ghana_2017_overall_flu_summary.csv"), row.names = FALSE)

