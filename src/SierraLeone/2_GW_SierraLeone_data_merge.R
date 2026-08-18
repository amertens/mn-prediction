

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



d <- readRDS(here("data", "IPD", "Sierra Leone", "SierraLeone_GMS_cleaned.rds"))


gw_vars <- colnames(d)

#-------------------------------------------------------------------------------
# GEE data
#-------------------------------------------------------------------------------

gee_buffer_path <- {
  cands <- Sys.glob(here("data/GEE/SL2012_buffers_*.csv"))
  if (length(cands) == 0) stop("No SL2012_buffers_*.csv found in data/GEE/")
  cands[which.max(file.info(cands)$mtime)]
}
gee <- read.csv(gee_buffer_path) %>% select(SR,cnum.x, trmm_Jan_10km:grassland_50km) %>% rename(cnum=cnum.x)
head(gee)
colnames(gee)= paste0("gee_",colnames(gee))
head(gee)



gee_vars <- colnames(gee)[-c(1,2)]

d$gw_cnum
gee$gee_cnum
d <- left_join(d, gee, by = c("gw_cnum"="gee_cnum"))
table(is.na(gee$gee_grassland_50km))
table(is.na(d$gee_trmm_Jan_10km))

# F-1: Derive seasonality summary features from monthly GEE buffer columns.
source(here("src/GEE/seasonality_features.R"))
d <- add_seasonality_features(d, verbose = TRUE)
gee_vars <- unique(c(gee_vars,
  grep("_(annual_mean|seasonal_sd|seasonal_cv|seasonal_range|peak_month)$",
       colnames(d), value = TRUE)))



#-------------------------------------------------------------------------------
#  Admin 1 and 2 membership
#-------------------------------------------------------------------------------

summary(d$latitude)
summary(d$longitude)
d$lat= as.numeric(d$latitude)
d$lon= as.numeric(d$longitude)

source(here("R", "data_prep.R"))
poly.adm <- load_gadm_cached("SLE", level = 2)
poly.adm <- sf::st_as_sf(poly.adm) %>% select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
d_sf <- st_as_sf(d, coords = c("longitude","latitude"), crs = 4326)
poly.adm <- st_transform(poly.adm, crs = 4326)
df <- (st_join(d_sf, poly.adm, join = st_within))

table(df$Admin1)
table(df$Admin2)

#get unique admin-2's for future merges
Admin1=unique(poly.adm$Admin1)
Admin2=unique(poly.adm$Admin2)


#-------------------------------------------------------------------------------
# GEE Admin-2 zonal-mean predictors
# Built by scripts/build_gee_admin2.R from country & global GEE rasters.
# Provides one column per raster band (prefix gee_a2_), joined on Admin2.
#-------------------------------------------------------------------------------

gee_a2_path <- here("data/GEE/SierraLeone_2013_admin2_gee.csv")
if (file.exists(gee_a2_path)) {
  gee_a2 <- read.csv(gee_a2_path, check.names = FALSE)
  gee_a2$Admin2 <- trimws(gee_a2$Admin2)
  if ("Admin1" %in% colnames(gee_a2)) gee_a2$Admin1 <- trimws(gee_a2$Admin1)
  # Join on (Admin1, Admin2) when both present — disambiguates duplicate
  # Admin2 names across regions (e.g. Malawi TAs sharing names across districts).
  gee_a2_keys <- intersect(c("Admin1", "Admin2"), colnames(gee_a2))
  gee_a2_vars <- setdiff(colnames(gee_a2), gee_a2_keys)
  df <- df %>% dplyr::left_join(gee_a2, by = gee_a2_keys)
  gee_vars <- unique(c(gee_vars, gee_a2_vars))
  cat(sprintf("  GEE admin2 merge: %d gee_a2_ columns added\n",
              length(gee_a2_vars)))
} else {
  warning("Admin2 GEE CSV not found — run scripts/build_gee_admin2.R SierraLeone")
}


#-------------------------------------------------------------------------------
# SoilGrids/ISRIC Admin-2 predictors
# Built by scripts/build_soilgrids_admin2.R. pH, organic carbon, nitrogen,
# clay, sand, silt, CEC at 0-5cm depth, zonal-mean per Admin-2.
#-------------------------------------------------------------------------------

soil_path <- here("data/SoilGrids/SierraLeone_soilgrids_admin2.csv")
if (file.exists(soil_path)) {
  soil_df <- read.csv(soil_path, check.names = FALSE)
  soil_df$Admin2 <- trimws(soil_df$Admin2)
  if ("Admin1" %in% colnames(soil_df)) soil_df$Admin1 <- trimws(soil_df$Admin1)
  soil_keys <- intersect(c("Admin1", "Admin2"), colnames(soil_df))
  soil_vars <- setdiff(colnames(soil_df), soil_keys)
  df <- df %>% dplyr::left_join(soil_df, by = soil_keys)
  cat(sprintf("  SoilGrids merge: %d soil_ columns added\n", length(soil_vars)))
} else {
  soil_vars <- character(0)
  warning("SoilGrids CSV not found — run scripts/build_soilgrids_admin2.R SierraLeone")
}

#-------------------------------------------------------------------------------
# ESPEN helminth / schistosomiasis Admin-2 predictors (iron pathway) — GUARDED
# Built by scripts/build_espen_admin2.R (STH/SCH prevalence -> iron deficiency).
#-------------------------------------------------------------------------------
espen_path <- here("data/ESPEN/SierraLeone_espen_admin2.csv")
if (file.exists(espen_path)) {
  espen_df <- read.csv(espen_path, check.names = FALSE)
  espen_df$Admin2 <- trimws(espen_df$Admin2)
  if ("Admin1" %in% colnames(espen_df)) espen_df$Admin1 <- trimws(espen_df$Admin1)
  espen_keys <- intersect(c("Admin1", "Admin2"), colnames(espen_df))
  espen_vars <- setdiff(colnames(espen_df), espen_keys)
  df <- df %>% dplyr::left_join(espen_df, by = espen_keys)
  cat(sprintf("  ESPEN merge: %d espen_ columns added\n", length(espen_vars)))
} else {
  espen_vars <- character(0)
  message("ESPEN CSV not found — run scripts/build_espen_admin2.R SierraLeone (optional domain)")
}

#-------------------------------------------------------------------------------
# MapSPAM crop-composition Admin-2 predictors (dietary availability) — GUARDED
# Built by scripts/build_mapspam_admin2.R (crop mix -> micronutrient availability).
#-------------------------------------------------------------------------------
spam_path <- here("data/MapSPAM/SierraLeone_spam_admin2.csv")
if (file.exists(spam_path)) {
  spam_df <- read.csv(spam_path, check.names = FALSE)
  spam_df$Admin2 <- trimws(spam_df$Admin2)
  if ("Admin1" %in% colnames(spam_df)) spam_df$Admin1 <- trimws(spam_df$Admin1)
  spam_keys <- intersect(c("Admin1", "Admin2"), colnames(spam_df))
  spam_vars <- setdiff(colnames(spam_df), spam_keys)
  df <- df %>% dplyr::left_join(spam_df, by = spam_keys)
  cat(sprintf("  MapSPAM merge: %d spam_ columns added\n", length(spam_vars)))
} else {
  spam_vars <- character(0)
  message("MapSPAM CSV not found — run scripts/build_mapspam_admin2.R SierraLeone (optional domain)")
}

#-------------------------------------------------------------------------------
# National causal-driver Admin-2 predictors (broadcast national value) — GUARDED
#   vas_ Vitamin A supplementation      (scripts/build_vas_national.R)
#   fao_ FAOSTAT food/nutrient supply   (scripts/build_faostat_supply.R)
#   fpn_ Food Prices for Nutrition      (scripts/build_fpn_affordability.R)
#-------------------------------------------------------------------------------
for (.nd in list(list(v="vas", path=here("data/VAS/SierraLeone_vas_admin2.csv")),
                 list(v="fao", path=here("data/FAOSTAT/SierraLeone_fao_admin2.csv")),
                 list(v="fpn", path=here("data/FPN/SierraLeone_fpn_admin2.csv")))) {
  assign(paste0(.nd$v, "_vars"), character(0))
  if (file.exists(.nd$path)) {
    .ndf <- read.csv(.nd$path, check.names = FALSE)
    .ndf$Admin2 <- trimws(.ndf$Admin2)
    if ("Admin1" %in% colnames(.ndf)) .ndf$Admin1 <- trimws(.ndf$Admin1)
    .keys <- intersect(c("Admin1", "Admin2"), colnames(.ndf))
    assign(paste0(.nd$v, "_vars"), setdiff(colnames(.ndf), .keys))
    df <- df %>% dplyr::left_join(.ndf, by = .keys)
    cat(sprintf("  %s merge: %d %s_ columns added\n", toupper(.nd$v),
                length(setdiff(colnames(.ndf), .keys)), .nd$v))
  } else message(sprintf("%s CSV not found — run its builder (optional national domain)", toupper(.nd$v)))
}


#-------------------------------------------------------------------------------
# Food price
#-------------------------------------------------------------------------------

wfp <- read.csv(here("data/food_price/wfp_food_prices_sle.csv"))
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
         unit=factor(unit)) %>% filter(year==2011)

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

WFP_MAX_DIST_KM <- 100
far <- !is.na(df$nearest_market_distance_km) &
        df$nearest_market_distance_km > WFP_MAX_DIST_KM
if (any(far)) {
  cat(sprintf("  WFP distance cutoff: %d/%d rows beyond %d km — wfp_* set to NA\n",
              sum(far), nrow(df), WFP_MAX_DIST_KM))
  wfp_data_cols <- intersect(wfp_vars, colnames(df))
  for (cc_wfp in wfp_data_cols) df[[cc_wfp]][far] <- NA
}


#
# #subset for methods example
# # d$gw_cRBPAdjBrinda
# #
# # gw_cVAD
#
# library(washb)
#
# # gee_built_surface
# # gee_ndvi_2_NDVI_Mar_50km
# # gee_temp_1_Feb_10km
# # gee_smod_code
# d_clust <- d %>% group_by(gw_cluster_id) %>% summarize(y, .keep_all=TRUE) %>% ungroup()
#
# gee_pre_screened_vars =washb_prescreen(Y=d$gw_cVAD, Ws=d%>%select(all_of(gee_vars, wfp_vars)), family = "gaussian", pval = 0.05, print = TRUE)




#-------------------------------------------------------------------------------
# Malaria Atlas
#-------------------------------------------------------------------------------

#Note! Need to check the different datasets and make sure I'm getting the right year

# Blood-disorder rasters confound iron biomarker interpretation:
#   HbS/HbC: chronic hemolysis raises some Brinda inflammation markers
#   G6PDd:   hemolytic anemia confounder
#   Duffy:   Plasmodium vivax susceptibility (sanity covariate)
# Auto-download any missing files via the malariaAtlas package (idempotent).
ensure_map_blood_disorders <- function(map_dir, iso) {
  ids <- c(
    "Blood_Disorders__201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency",
    "Blood_Disorders__201201_Africa_HbC_Allele_Frequency",
    "Blood_Disorders__201201_Global_G6PDd_Allele_Frequency",
    "Blood_Disorders__201201_Global_Duffy_Negativity_Phenotype_Frequency"
  )
  files <- paste0(ids, ".tif")
  missing <- files[!file.exists(file.path(map_dir, files))]
  if (length(missing) == 0) return(invisible(files))
  if (!requireNamespace("malariaAtlas", quietly = TRUE))
    stop("Install malariaAtlas: install.packages('malariaAtlas')")
  dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)
  shp <- malariaAtlas::getShp(ISO = iso, admin_level = "admin0")
  for (mf in missing) {
    id <- sub("\\.tif$", "", mf)
    r  <- malariaAtlas::getRaster(dataset_id = id, shp = shp)
    if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
    terra::writeRaster(r, file.path(map_dir, mf), overwrite = TRUE)
  }
  invisible(files)
}
blood_disorder_files <- ensure_map_blood_disorders(
  here("data/Malaria Atlas/SierraLeone"), iso = "SLE"
)

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
             "Malaria__202202_Global_Pf_Reproductive_Number.tif",
             blood_disorder_files)


pts = data.frame(lon=df$lon, lat=df$lat)


for(i in rasters){
  rast <- rast(here(paste0("data/Malaria Atlas/SierraLeone/",i)))
  rast_extract <- terra::extract(rast, pts, method="bilinear")
  col_name <- i
  col_name <- gsub(".tif", "", i)
  col_name <- gsub("Malaria__", "MAP_", col_name)
  col_name <- gsub("Interventions__", "MAP_", col_name)
  col_name <- gsub("Blood_Disorders__", "MAP_", col_name)
  col_name <- gsub("Global_", "", col_name)
  col_name <- gsub("Africa_", "", col_name)
  # col_name <- gsub("202106_", "", col_name)
  # col_name <- gsub("202206_", "", col_name)
  # col_name <- gsub("202406_", "", col_name)
  df$temp <- rast_extract[, 2]
  colnames(df)[ncol(df)] <- col_name
}

map_vars <- colnames(df)[grepl("MAP_", colnames(df))]

length(colnames(df))
length(unique(colnames(df)))

table(is.na(df$MAP_202106_Insecticide_Treated_Net_Use))

#-------------------------------------------------------------------------------
# IHME Data
#-------------------------------------------------------------------------------

ihme <- read.csv(here("data/IHME/SierraLeone_2012_merged_IHME_data.csv"))
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

# Admin1-level IHME predictors (built by SL_IHME_clean.R via build_ihme_admin1).
source(here("src/IHME/build_ihme_admin1.R"))
df <- merge_ihme_admin1(
  df,
  admin1_csv  = here("data/IHME/SierraLeone_2012_merged_IHME_admin1_data.csv"),
  admin1_col  = "Admin1"
)
ihme_vars <- unique(c(ihme_vars, grep("^ihme_adm1_", colnames(df), value = TRUE)))


#-------------------------------------------------------------------------------
# Food security (Cadre Harmonisé + HFID / FEWS NET IPC)
#-------------------------------------------------------------------------------

source(here("R", "food_security.R"))
source(here("R", "config.R"))
cc_fsec  <- get_country_configs()[["SierraLeone"]]
hfid_path <- here("data", "HFID", "hfid_hv1.csv")
ch_path   <- here("data", "CadreHarmonise", "cadre_harmonise_caf_ipc_dec25.xlsx")
if (file.exists(hfid_path) || file.exists(ch_path)) {
  df <- as.data.frame(df)
  df <- merge_food_security(df, cc_fsec,
                             hfid_path = if (file.exists(hfid_path)) hfid_path else NULL,
                             ch_path   = if (file.exists(ch_path))   ch_path   else NULL)
  fsec_vars <- grep("^fsec_", colnames(df), value = TRUE)
} else {
  fsec_vars <- character(0)
  warning("Neither HFID nor Cadre Harmonisé file found; skipping fsec_ merge.")
}


#-------------------------------------------------------------------------------
# MICS Data
#-------------------------------------------------------------------------------

#see if I can aggregate at the admin-2 level- only 4 units at admin1
#need to clean the MICS 4 data instead of MICS 6

# mics <- read.csv(here("data/MICS/mics_ghana_2017_region_summary.csv"))
# table(mics$region)
# #mics$region[mics$region=="Brong Ahafo"] <- "Ahafo"
# table(df$Admin1)
#
# mics_vars <- colnames(mics)
#
# df <- left_join(df, mics, by = c("Admin1_old" = "mics_region"))
# summary(df$mics_hc4)



#-------------------------------------------------------------------------------
# Conflict Data
#-------------------------------------------------------------------------------

#()


#-------------------------------------------------------------------------------
# Food security Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Cadre Harmonise Data
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DHS Admin-1 indicators (from surveyPrev direct estimates, pivoted to wide)
#-------------------------------------------------------------------------------

source(here("R", "data_prep.R"))
dhs2013_adm1 <- load_dhs_admin1(
  dhs_dir = here("data", "DHS", "clean"),
  country = "Sierra Leone",
  year    = 2013
)
if (!is.null(dhs2013_adm1)) {
  dhs_vars <- names(dhs2013_adm1)
  df <- left_join(as.data.frame(df), dhs2013_adm1, by = c("Admin1" = "dhs2013_DHSREGEN"))
  cat(sprintf("  DHS admin-1 merge complete: %d dhs2013_ columns added\n",
              sum(grepl("^dhs2013_", colnames(df)))))
} else {
  dhs_vars <- c()
  warning("Admin-1 DHS file not found for Sierra Leone 2013")
}

#-------------------------------------------------------------------------------
# DHS Admin-2 indicators (from surveyPrev FH BYM2 smoothed estimates)
# Run src/DHS/DHS_admin2_aggregation.R first to generate the wide file
#-------------------------------------------------------------------------------

source(here("R", "data_prep.R"))
source(here("R", "food_security.R"))

dhs2013_adm2 <- load_dhs_admin2(
  dhs_dir = here("data", "DHS", "clean"),
  country = "Sierra Leone",
  year    = 2013
)
if (!is.null(dhs2013_adm2)) {
  target_a2 <- unique(df$Admin2)
  source_a2 <- dhs2013_adm2$Admin2
  name_lookup <- build_name_lookup(target_a2, source_a2)
  dhs2013_adm2$Admin2 <- name_lookup[dhs2013_adm2$Admin2]
  n_matched <- sum(!is.na(dhs2013_adm2$Admin2))
  cat(sprintf("  DHS 2013 Admin2 name matching: %d/%d matched\n", n_matched, length(source_a2)))
  dhs2013_adm2 <- dhs2013_adm2[!is.na(dhs2013_adm2$Admin2), ]
  df <- left_join(df, dhs2013_adm2, by = "Admin2")
  cat(sprintf("  DHS 2013 admin-2 merge complete\n"))
} else {
  warning("Admin-2 DHS file not found for Sierra Leone 2013")
}

#-------------------------------------------------------------------------------
# FluNet Data
#-------------------------------------------------------------------------------

# flu <- read.csv(here("data/FluNet/Ghana/ghana_2017_weekly_flu_summary.csv"))
# head(flu)
#
# flu_vars <- colnames(flu)
#
# df <- left_join(df, flu, by = c("gw_week" = "flunet_week"))
# table(is.na(df$flunet_inf_b_rate))

#-------------------------------------------------------------------------------
# LSMS
#-------------------------------------------------------------------------------

# #Individual measures
# lsms <- readRDS(here("data/LSMS/Ghana_LSMS_clean.RDS"))
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
# df <- df %>% subset(., select = -c(dhs2014_REG_ID, dhs2014_REGCODE, dhs2014_REGNAME, dhs2014_REGNOTES, dhs2014_geometry,
#                                    dhs2014_REG_ID, dhs2014_REGNAME,
#                                    dhs2016_REG_ID, dhs2016_REGCODE, dhs2016_REGNAME, dhs2016_REGNOTES, dhs2016_geometry,
#                                    dhs2016_REG_ID, dhs2016_REGNAME,
#                                    dhs2017_REG_ID, dhs2017_REGCODE, dhs2017_REGNAME, dhs2017_REGNOTES, dhs2017_geometry,
#                                    dhs2017_REG_ID, dhs2017_REGNAME,
#                                    ihme_adm2_name , Admin2_old ,
#                                    gw_District, `gw_District Name`, gw_EANAME, gw_EACode))





#-------------------------------------------------------------------------------
# check data
#-------------------------------------------------------------------------------

#remove any spatial or time (not date) columns
df = st_drop_geometry(df)
df <- df[, !sapply(df, function(x) inherits(x, c("POSIXct", "POSIXt")))]

if ("dhs2013_geometry" %in% colnames(df)) {
  df <- df %>% subset(., select=-c(dhs2013_geometry))
}

#make unique id
df$dataid <- paste0("SL",1:nrow(df))

#-------------------------------------------------------------------------------
# Save data
#-------------------------------------------------------------------------------


saveRDS(df, file=here("data", "IPD", "Sierra Leone", "SierraLeone_merged_dataset.rds"))

#-------------------------------------------------------------------------------
# Save metadata
#-------------------------------------------------------------------------------

metadata <- list(
  gw_vars = gw_vars[gw_vars %in% colnames(df)],
  dhs_vars = dhs_vars[dhs_vars %in% colnames(df)],
  #mics_vars = mics_vars[mics_vars %in% colnames(df)],
  ihme_vars = ihme_vars[ihme_vars %in% colnames(df)],
  #lsms_vars = lsms_vars[lsms_vars %in% colnames(df)],
  map_vars = map_vars[map_vars %in% colnames(df)],
  wfp_vars = wfp_vars[wfp_vars %in% colnames(df)],
  #flunet_vars = flu_vars[flu_vars %in% colnames(df)],
  gee_vars = gee_vars[gee_vars %in% colnames(df)],
  fsec_vars = fsec_vars[fsec_vars %in% colnames(df)],
  soil_vars = soil_vars[soil_vars %in% colnames(df)],
  espen_vars = espen_vars[espen_vars %in% colnames(df)],
  spam_vars = spam_vars[spam_vars %in% colnames(df)],
  vas_vars = vas_vars[vas_vars %in% colnames(df)],
  fao_vars = fao_vars[fao_vars %in% colnames(df)],
  fpn_vars = fpn_vars[fpn_vars %in% colnames(df)]
)


saveRDS(metadata, here("metadata/SL_variable_categories.rds"))

