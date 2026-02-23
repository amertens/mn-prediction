

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



d <- readRDS(here("data", "IPD", "Gambia", "Gambia_GMS_cleaned.rds"))


gw_vars <- colnames(d)

#-------------------------------------------------------------------------------
# GEE data
#-------------------------------------------------------------------------------


gee <- read.csv(here("data/gee/gambia2018_buffers_01.08.2026.csv")) %>% select(SR,MICS_Cluster_Number.x, trmm_Jan_10km:grassland_50km) #%>%
#mutate(EACode=as.character(EACode))
colnames(gee)= paste0("gee_",colnames(gee))
head(gee)


gee_vars <- colnames(gee)[-c(1,2)]

table(d$gw_cnum)
table(gee$gee_cnum.x)

table(d$gw_MICS_Cluster_Number)
table(gee$gee_MICS_Cluster_Number.x)

unique(d$gw_EACode)
unique(gee$gee_EACode)

#d <- left_join(d, gee, by = c("longitude"="gee_longnum", "latitude"="gee_latnum"))
d <- left_join(d, gee, by = c("gw_MICS_Cluster_Number"="gee_MICS_Cluster_Number.x"))
table(is.na(gee$gee_grassland_50km))
table(is.na(d$gee_grassland_50km))



#-------------------------------------------------------------------------------
#  Admin 1 and 2 membership
#-------------------------------------------------------------------------------

summary(d$latitude)
summary(d$longitude)
d$lat= as.numeric(d$latitude)
d$lon= as.numeric(d$longitude)

poly.adm <- geodata::gadm(country="GM", level=2, path=tempdir())
poly.adm <- sf::st_as_sf(poly.adm) %>% select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
d_sf <- st_as_sf(d, coords = c("longitude","latitude"), crs = 4326)
poly.adm <- st_transform(poly.adm, crs = 4326)


#update admin 1 naming to old convention to match DHS and MICS
poly.adm <- poly.adm %>%
  mutate(
    Admin1_old = dplyr::case_when(
      Admin1 == "Banjul" ~ "Banjul",
      Admin1 == "Western" ~ "Brikama",
      Admin1 == "Lower River" ~ "Mansakonko",
      Admin1 == "North Bank" ~ "Kerewan",
      Admin1 == "Upper River" ~ "Basse",

      # Central River North → Kuntaur
      Admin1 == "Maccarthy Island" & Admin2 %in% c(
        "Central Baddibu",
        "Fulladu East",
        "Fulladu West",
        "Lower Saloum",
        "Niamina Dankunku",
        "Niamina East",
        "Niamina West",
        "Niani"
      ) ~ "Kuntaur",

      # Central River South → Janjanbureh
      Admin1 == "Maccarthy Island" & Admin2 %in% c(
        "Janjanbureh",
        "Kantora",
        "Nianija",
        "Sami",
        "Upper Saloum"
      ) ~ "Janjanbureh",

      TRUE ~ NA_character_
    )
  )


table(poly.adm$Admin1, (poly.adm$Admin1_old))
table(poly.adm$Admin1, is.na(poly.adm$Admin1_old))


#merge with survey
df <- (st_join(d_sf, poly.adm, join = st_within))


table(df$Admin1)
table(df$Admin1_old)
table(df$Admin2)


#get unique admin-2's for future merges
Admin1=unique(poly.adm$Admin1)
Admin2=unique(poly.adm$Admin2)




#-------------------------------------------------------------------------------
# Food price
#-------------------------------------------------------------------------------

wfp <- read.csv(here("data/food_price/wfp_food_prices_gmb.csv"))
head(wfp)
wfp <- wfp[-1,]


#subset to keep: market name, lat, long. Also only keep one row from long format
markets_subset <- wfp %>% #drop top row, variable description
  group_by(market) %>%
  slice(1) %>%
  ungroup() %>%
  select(market, latitude, longitude)

# Convert to sf
#DHS_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
markets_sf <- st_as_sf(markets_subset, coords = c("longitude", "latitude"), crs = 4326)
st_geometry(df)

# Project to meters for accurate distance, drop z coordinate first
df <- st_zm(df, drop = TRUE, what = "ZM")
DHS_proj <- st_transform(df, 3857)  # Web Mercator
markets_proj <- st_transform(markets_sf, 3857)

# Distance matrix: rows = original points, columns = markets
dist_matrix <- st_distance(DHS_proj, markets_proj)

# Get the minimum distance and index of nearest market
nearest_dist <- apply(dist_matrix, 1, min)
nearest_market_index <- apply(dist_matrix, 1, which.min)

# Add the results to your original data
df$nearest_market_distance_km <- nearest_dist / 1000  # convert to km
df$nearest_market_id <- markets_subset$market[nearest_market_index]


#merge in food pricing data

#get the average food prices the year of sampling (should this be the year before?)
table(wfp$priceflag)
table(wfp$pricetype )
wfp <- wfp %>%
  mutate(year=year(wfp$date),
         usdprice=as.numeric(usdprice),
         priceflag=factor(priceflag, levels=c("actual","actual,aggregate","aggregate")),
         pricetype = factor(pricetype, levels=c("Retail","Wholesale")),
         unit=factor(unit)) %>% filter(year==2017)

levels(wfp$unit) <- c("KG", setdiff(levels(wfp$unit), "KG"))

ave_price <- wfp %>%
  group_by(market, commodity , currency,priceflag, pricetype, unit) %>%
  summarise(price = mean(usdprice, na.rm=T),
            sd_price=sd(usdprice, na.rm=T)) %>%
  group_by(market, commodity , currency) %>%
  arrange(market, commodity , currency, priceflag, pricetype, unit) %>% slice(1) %>%
  ungroup()
head(ave_price)


#get the price df to merge
price_df <- ave_price %>% select(market, commodity, price ) %>%
  # transform commodity to wide
  pivot_wider(names_from = commodity, values_from = price) %>%
  janitor::clean_names() %>%
  rename(nearest_market_id=market)

head(price_df)
colnames(price_df) <- paste0("wfp_", colnames(price_df))
wfp_vars <- colnames(price_df)



head(df)
df <- left_join(df, price_df, by = c("nearest_market_id" = "wfp_nearest_market_id"))
table(is.na(price_df$wfp_cassava ))
table(is.na(df$wfp_cassava ))


#-------------------------------------------------------------------------------
# Malaria Atlas
#-------------------------------------------------------------------------------

#Note! Need to check the different datasets and make sure I'm getting the right year

rasters <- c("Malaria__202206_Global_Pf_Incidence_Count.tif",
             "Malaria__202206_Global_Pf_Incidence_Rate.tif",
             "Malaria__202206_Global_Pf_Mortality_Count.tif",
             "Malaria__202206_Global_Pf_Mortality_Rate.tif",
             "Malaria__202206_Global_Pf_Parasite_Rate.tif",
             "Malaria__202206_Global_Pv_Incidence_Count.tif",
             "Malaria__202206_Global_Pv_Incidence_Rate.tif",
             "Malaria__202206_Global_Pv_Parasite_Rate.tif",
             "Malaria__202406_Global_Pf_Incidence_Count.tif",
             "Malaria__202406_Global_Pf_Incidence_Rate.tif",
             "Malaria__202406_Global_Pf_Mortality_Count.tif",
             "Malaria__202406_Global_Pf_Mortality_Rate.tif",
             "Malaria__202406_Global_Pf_Parasite_Rate.tif",
             "Malaria__202406_Global_Pv_Incidence_Count.tif",
             "Malaria__202406_Global_Pv_Incidence_Rate.tif",
             "Malaria__202406_Global_Pv_Parasite_Rate.tif",
             "Interventions__202106_Africa_Insecticide_Treated_Net_Access.tif",
             "Interventions__202106_Africa_Insecticide_Treated_Net_Use.tif",
             "Interventions__202106_Africa_Insecticide_Treated_Net_Use_Rate.tif",
             "Interventions__202106_Africa_IRS_Coverage.tif",
             "Interventions__202106_Global_Antimalarial_Effective_Treatment.tif",
             "Interventions__202406_Africa_Insecticide_Treated_Net_Access.tif",
             "Interventions__202406_Africa_Insecticide_Treated_Net_Use.tif",
             "Interventions__202406_Africa_Insecticide_Treated_Net_Use_Rate.tif",
             "Interventions__202406_Global_Antimalarial_Effective_Treatment.tif",
             "Malaria__202202_Global_Pf_Reproductive_Number.tif")


pts = data.frame(lon=df$lon, lat=df$lat)


for(i in rasters){
  rast <- rast(here(paste0("data/Malaria Atlas/Gambia/",i)))
  rast_extract <- terra::extract(rast, pts, method="bilinear")
  col_name <- i
  col_name <- gsub(".tif", "", i)
  col_name <- gsub("Malaria__", "MAP_", col_name)
  col_name <- gsub("Interventions__", "MAP_", col_name)
  col_name <- gsub("Global_", "", col_name)
  col_name <- gsub("Africa_", "", col_name)
  df$temp <- rast_extract[, 2]
  colnames(df)[ncol(df)] <- col_name
}

map_vars <- colnames(df)[grepl("MAP_", colnames(df))]

length(colnames(df))
length(unique(colnames(df)))

#-------------------------------------------------------------------------------
# IHME Data
#-------------------------------------------------------------------------------

ihme <- read.csv(here("data/IHME/Gambia_2018_merged_IHME_data.csv"))
head(ihme)

# Load required libraries
library(dplyr)
library(stringdist)
library(stringr)

# Extract unique names for analysis
df_admin2 <- unique(df$Admin2)
ihme_admin2 <- unique(ihme$ihme_adm2_name)

# Function to clean admin names for better matching
clean_admin_names <- function(names) {
  names %>%
    # Remove common administrative suffixes
    gsub("\\s+(Municipal|Metropolis|Metropolitan|Metro)$", "", .) %>%
    # Standardize separators
    gsub("[/-]", " ", .) %>%
    # Remove extra spaces and normalize
    str_squish() %>%
    # Convert to title case
    str_to_title() %>%
    # Handle common abbreviations
    gsub("\\bKma\\b", "Kumasi", .) %>%
    # Remove parentheses and contents
    gsub("\\s*\\([^)]*\\)", "", .) %>%
    # Trim whitespace
    str_trim()
}

# Clean the names
df_clean <- data.frame(
  original = df_admin2,
  cleaned = clean_admin_names(df_admin2),
  stringsAsFactors = FALSE
)

ihme_clean <- data.frame(
  original = ihme_admin2,
  cleaned = clean_admin_names(ihme_admin2),
  stringsAsFactors = FALSE
)

# Function to perform comprehensive fuzzy matching
fuzzy_match_admin2 <- function(df_names, ihme_names, threshold = 0.7) {

  cat("=== PERFORMING FUZZY MATCHING ===\n")

  # Calculate similarity matrix using multiple methods
  methods <- c("jw", "cosine", "jaccard")

  results_list <- list()

  for(method in methods) {
    # Calculate distance matrix
    dist_matrix <- stringdistmatrix(
      df_names$cleaned,
      ihme_names$cleaned,
      method = method
    )

    # Convert to similarity
    sim_matrix <- 1 - dist_matrix

    # Find best matches
    matches <- data.frame(
      df_admin2 = character(),
      ihme_admin2 = character(),
      similarity = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )

    for(i in 1:nrow(df_names)) {
      best_match_idx <- which.max(sim_matrix[i, ])
      best_similarity <- sim_matrix[i, best_match_idx]

      matches <- rbind(matches, data.frame(
        df_admin2 = df_names$original[i],
        ihme_admin2 = ihme_names$original[best_match_idx],
        df_cleaned = df_names$cleaned[i],
        ihme_cleaned = ihme_names$cleaned[best_match_idx],
        similarity = best_similarity,
        method = method,
        stringsAsFactors = FALSE
      ))
    }

    results_list[[method]] <- matches
  }

  return(results_list)
}

# Perform fuzzy matching with multiple methods
matching_results <- fuzzy_match_admin2(df_clean, ihme_clean, threshold = 0.9)

# Combine results and pick best match for each district
combined_matches <- bind_rows(matching_results, .id = "method") %>%
  group_by(df_admin2) %>%
  slice_max(similarity, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    match_quality = case_when(
      similarity >= 0.95 ~ "Excellent",
      similarity >= 0.85 ~ "Very Good",
      similarity >= 0.75 ~ "Good",
      similarity >= 0.65 ~ "Fair",
      TRUE ~ "Poor"
    )
  ) %>%
  arrange(desc(similarity))



# Create lookup table for merging
lookup_table <- combined_matches %>%
  select(
    Admin2 = df_admin2,
    ihme_adm2_name = ihme_admin2,
    similarity,
    match_quality,
    matching_method = method
  )

# Function to merge datasets using the lookup
merge_with_ihme <- function(main_df, ihme_df, lookup_table, min_similarity = 0.70) {

  cat("\n=== MERGING DATASETS ===\n")

  # Filter lookup for acceptable matches
  good_matches <- lookup_table %>%
    filter(similarity >= min_similarity)

  cat("Using", nrow(good_matches), "matches with similarity >=", min_similarity, "\n")

  # Perform the merge
  merged_data <- main_df %>%
    left_join(good_matches, by = "Admin2") %>%
    left_join(ihme_df, by = "ihme_adm2_name")

  # Summary statistics
  cat("Original df rows:", nrow(main_df), "\n")
  cat("Merged rows:", nrow(merged_data), "\n")
  cat("Successfully matched:",
      sum(!is.na(merged_data$ihme_adm2_name)), "districts\n")
  cat("Unmatched:",
      sum(is.na(merged_data$ihme_adm2_name)), "districts\n")

  return(merged_data)
}

ihme_vars <- colnames(ihme)

# merge datasets
df <- merge_with_ihme(df, ihme, lookup_table, min_similarity = 0.90)
table(is.na(df$ihme_2_to_10_years_2_to_10_both_malaria_prevalence_rate ))



#-------------------------------------------------------------------------------
# MICS Data
#-------------------------------------------------------------------------------

mics <- read.csv(here("data/MICS/mics_Gambia_2018_region_summary.csv"))
table(mics$mics_region)
table(df$Admin1)
table(df$Admin1_old)

mics_vars <- colnames(mics)

df <- left_join(df, mics, by = c("Admin1" = "mics_region"))
summary(df$mics_hc4)



#-------------------------------------------------------------------------------
# Conflict Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Food security Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Cadre Harmonise Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# DHS Admin 1 indicators
#-------------------------------------------------------------------------------


dhs2019 <- readRDS(here("data/DHS/clean/Gambia_2019_dhs_aggregation.rds"))
colnames(dhs2019) <- paste0("dhs2019_",colnames(dhs2019))

table(dhs2019$dhs2019_DHSREGEN)
table(df$Admin1)
table(df$Admin1_old)

dhs_vars <- c(colnames(dhs2019))

df <- left_join(as.data.frame(df), dhs2019, by = c("Admin1_old" = "dhs2019_DHSREGEN"))
table(is.na(df$dhs2019_AH_CIGA_M_59C))

#-------------------------------------------------------------------------------
# FluNet Data
#-------------------------------------------------------------------------------

#no flunet in Gambia

#-------------------------------------------------------------------------------
# LSMS
#-------------------------------------------------------------------------------

#TO ADD

# #Individual measures
# lsms <- readRDS(here("data/LSMS/Gambia_LSMS_clean.RDS"))
# head(lsms)
#
# unique(df$Admin1_old)[!(unique(df$Admin1_old) %in% unique(lsms$lsms_admin1))]
# lsms$lsms_admin1[!(unique(lsms$lsms_admin1) %in% unique(df$Admin1_old))]
#
# lsms_vars <- colnames(lsms)
#
# unique(lsms$lsms_admin1)
# unique(df$Admin1)
# df <- left_join(df, lsms, by = c("Admin1_old" = "lsms_admin1"))


#-------------------------------------------------------------------------------
# clean data
#-------------------------------------------------------------------------------


#check for any columns that have many factor levels
# Alternative function that also considers character variables
check_categorical_levels <- function(df, max_levels = 10, include_character = TRUE) {
  # Function to count unique values
  count_unique <- function(x) {
    if (is.factor(x)) {
      return(length(levels(x)))
    } else if (is.character(x)) {
      return(length(unique(x)))
    } else {
      return(NA)
    }
  }

  # Identify categorical columns (factors and optionally characters)
  if (include_character) {
    categorical_cols <- sapply(df, function(x) is.factor(x) || is.character(x))
  } else {
    categorical_cols <- sapply(df, is.factor)
  }

  if (!any(categorical_cols)) {
    message("No categorical variables found in the dataframe.")
    return(list(
      categorical_summary = data.frame(),
      high_level_variables = character(0),
      summary_stats = list()
    ))
  }

  # Extract categorical columns
  categorical_data <- df[, categorical_cols, drop = FALSE]

  # Count levels/unique values for each variable
  level_counts <- sapply(categorical_data, count_unique)
  variable_types <- sapply(categorical_data, class)

  # Create summary dataframe
  categorical_summary <- data.frame(
    Variable = names(level_counts),
    Type = variable_types,
    Num_Levels = level_counts,
    Over_Limit = level_counts > max_levels,
    stringsAsFactors = FALSE
  )

  # Sort by number of levels (descending)
  categorical_summary <- categorical_summary[order(categorical_summary$Num_Levels, decreasing = TRUE), ]
  rownames(categorical_summary) <- NULL

  # Identify variables with more than X levels
  high_level_variables <- categorical_summary$Variable[categorical_summary$Over_Limit]

  # Summary statistics
  summary_stats <- list(
    total_categorical = nrow(categorical_summary),
    variables_over_limit = length(high_level_variables),
    max_levels_found = max(level_counts, na.rm = TRUE),
    mean_levels = round(mean(level_counts, na.rm = TRUE), 2),
    median_levels = median(level_counts, na.rm = TRUE)
  )

  return(list(
    categorical_summary = categorical_summary,
    high_level_variables = high_level_variables,
    summary_stats = summary_stats
  ))
}


very_high_fact <- check_categorical_levels(df, max_levels = 100, include_character = TRUE)
high_fact <- check_categorical_levels(df, max_levels = 10, include_character = TRUE)
print(high_fact$categorical_summary %>% filter(Over_Limit))

#drop unneeded columns
df <- df %>% select(-all_of(very_high_fact$high_level_variables))



#-------------------------------------------------------------------------------
# check data
#-------------------------------------------------------------------------------

#remove any spatial or time (not date) columns
df = st_drop_geometry(df)
df <- df[, !sapply(df, function(x) inherits(x, c("POSIXct", "POSIXt")))]

df <- df %>% subset(., select=-c(dhs2019_geometry))

#make unique id
df$dataid <- paste0("gambia",1:nrow(df))

#-------------------------------------------------------------------------------
# Save data
#-------------------------------------------------------------------------------


saveRDS(df, file=here("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"))

#-------------------------------------------------------------------------------
# Save metadata
#-------------------------------------------------------------------------------

metadata <- list(
  gw_vars = gw_vars[gw_vars %in% colnames(df)],
  dhs_vars = dhs_vars[dhs_vars %in% colnames(df)],
  mics_vars = mics_vars[mics_vars %in% colnames(df)],
  ihme_vars = ihme_vars[ihme_vars %in% colnames(df)],
  #lsms_vars = lsms_vars[lsms_vars %in% colnames(df)],
  map_vars = map_vars[map_vars %in% colnames(df)],
  wfp_vars = wfp_vars[wfp_vars %in% colnames(df)],
  #flunet_vars = flu_vars[flu_vars %in% colnames(df)],
  gee_vars = gee_vars[gee_vars %in% colnames(df)]
)


saveRDS(metadata, here("metadata/gambia_variable_categories.rds"))

