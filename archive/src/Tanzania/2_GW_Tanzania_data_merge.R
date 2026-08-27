# =============================================================================
# src/Tanzania/2_GW_Tanzania_data_merge.R
#
# SCAFFOLD — join the Tanzania (TDHS 2009-10) vitamin A outcome to predictor
# domains and produce data/IPD/Tanzania 2010/Tanzania_merged_dataset.rds
# (the file consumed by R/config.R -> Tanzania$data_path).
#
# Ported from src/SierraLeone/2_GW_SierraLeone_data_merge.R. KEY DIFFERENCE:
# every external-domain join is guarded with file.exists(). Tanzania's proxy
# extractions come online incrementally (see README_TANZANIA_TODO.md), so this
# script runs to completion with whatever exists — at minimum the GADM admin-2
# assignment — instead of hard-stopping on the first missing input. As each
# extraction lands, its block activates automatically.
#
# PREREQUISITE: run src/Tanzania/1_GW_Tanzania_data_clean.R first.
#
# Domains and the scripts that must produce their inputs (none run yet):
#   gee_       data/GEE/TZ2010_buffers_*.csv                 (EE cluster buffers)
#   gee_a2_    data/GEE/Tanzania_2010_admin2_gee.csv         scripts/build_gee_admin2.R
#   soil_      data/SoilGrids/Tanzania_soilgrids_admin2.csv  scripts/build_soilgrids_admin2.R
#   MAP_       data/Malaria Atlas/Tanzania/*.tif             (malariaAtlas auto-download)
#   ihme_      data/IHME/tanzania_2010_merged_IHME_data.csv  src/IHME (cf. SL_IHME_clean.R)
#   dhs2010_   data/DHS/clean (surveyPrev admin-1/2)         src/DHS/DHS_admin2_aggregation.R
#   wfp_/fsec_ optional (several countries skip these)
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(haven)
  library(here)
  library(sf)
  library(terra)
})

# Domain toggles — flip off to skip a heavy step during iteration.
# RUN_MAP=false skips the Malaria Atlas step (it needs the global malaria rasters
# placed under data/Malaria Atlas/Tanzania/; blood-disorder rasters auto-download).
RUN_MAP <- tolower(Sys.getenv("RUN_MAP", "true")) %in% c("1", "true", "yes")

ISO        <- "TZA"
SURVEY_YEAR <- 2010L
ID_PREFIX  <- "TZ"
tz_dir     <- here("data", "IPD", "Tanzania 2010")

d <- readRDS(file.path(tz_dir, "Tanzania_GMS_cleaned.rds"))
gw_vars <- colnames(d)

# Initialise per-domain var lists (stay empty if a domain is skipped).
gee_vars <- character(0); soil_vars <- character(0); map_vars <- character(0)
ihme_vars <- character(0); wfp_vars <- character(0); fsec_vars <- character(0)
dhs_vars <- character(0)

#-------------------------------------------------------------------------------
# GEE cluster-buffer predictors (join on cnum) — GUARDED
#-------------------------------------------------------------------------------
gee_buffer_cands <- Sys.glob(here("data/GEE/TZ2010_buffers_*.csv"))
if (length(gee_buffer_cands) > 0) {
  gee_buffer_path <- gee_buffer_cands[which.max(file.info(gee_buffer_cands)$mtime)]
  gee <- read.csv(gee_buffer_path)
  # Expect a cnum column plus monthly/buffer bands (adapt the select range to TZ).
  colnames(gee) <- paste0("gee_", colnames(gee))
  cnum_col <- grep("cnum", colnames(gee), value = TRUE)[1]
  d <- left_join(d, gee, by = setNames(cnum_col, "gw_cnum"))
  gee_vars <- setdiff(colnames(gee), cnum_col)
  # Seasonality summaries from monthly buffer columns (if present).
  if (file.exists(here("src/GEE/seasonality_features.R"))) {
    source(here("src/GEE/seasonality_features.R"))
    d <- add_seasonality_features(d, verbose = TRUE)
    gee_vars <- unique(c(gee_vars,
      grep("_(annual_mean|seasonal_sd|seasonal_cv|seasonal_range|peak_month)$",
           colnames(d), value = TRUE)))
  }
  cat(sprintf("  GEE buffers: %d gee_ columns added\n", length(gee_vars)))
} else {
  warning("No data/GEE/TZ2010_buffers_*.csv — skipping GEE buffer merge.")
}

#-------------------------------------------------------------------------------
# Admin-1 / Admin-2 membership via GADM (ALWAYS runs — the spatial join)
#-------------------------------------------------------------------------------
d$lat <- as.numeric(d$latitude)
d$lon <- as.numeric(d$longitude)
n_nogps <- sum(is.na(d$lat) | is.na(d$lon))
if (n_nogps > 0) cat(sprintf("  %d rows have no GPS — dropped before spatial join\n", n_nogps))
d <- d[!is.na(d$lat) & !is.na(d$lon), ]

source(here("R", "data_prep.R"))
poly.adm <- load_gadm_cached(ISO, level = 2)   # downloads gadm41_TZA_2 if uncached
poly.adm <- sf::st_as_sf(poly.adm) %>%
  select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
poly.adm <- st_transform(poly.adm, crs = 4326)
d_sf <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
df <- st_join(d_sf, poly.adm, join = st_within)
cat(sprintf("  GADM admin-2 join: %d rows, %d distinct Admin2, %d unmatched\n",
            nrow(df), length(unique(na.omit(df$Admin2))), sum(is.na(df$Admin2))))

#-------------------------------------------------------------------------------
# GEE Admin-2 zonal means (join on Admin1+Admin2) — GUARDED
#-------------------------------------------------------------------------------
gee_a2_path <- here("data/GEE/Tanzania_2010_admin2_gee.csv")
if (file.exists(gee_a2_path)) {
  gee_a2 <- read.csv(gee_a2_path, check.names = FALSE)
  gee_a2$Admin2 <- trimws(gee_a2$Admin2)
  if ("Admin1" %in% colnames(gee_a2)) gee_a2$Admin1 <- trimws(gee_a2$Admin1)
  gee_a2_keys <- intersect(c("Admin1", "Admin2"), colnames(gee_a2))
  gee_a2_vars <- setdiff(colnames(gee_a2), gee_a2_keys)
  df <- df %>% left_join(gee_a2, by = gee_a2_keys)
  gee_vars <- unique(c(gee_vars, gee_a2_vars))
  cat(sprintf("  GEE admin2 merge: %d gee_a2_ columns added\n", length(gee_a2_vars)))
} else {
  warning("Admin2 GEE CSV not found — run scripts/build_gee_admin2.R Tanzania")
}

#-------------------------------------------------------------------------------
# SoilGrids Admin-2 predictors (join on Admin1+Admin2) — GUARDED
#-------------------------------------------------------------------------------
soil_path <- here("data/SoilGrids/Tanzania_soilgrids_admin2.csv")
if (file.exists(soil_path)) {
  soil_df <- read.csv(soil_path, check.names = FALSE)
  soil_df$Admin2 <- trimws(soil_df$Admin2)
  if ("Admin1" %in% colnames(soil_df)) soil_df$Admin1 <- trimws(soil_df$Admin1)
  soil_keys <- intersect(c("Admin1", "Admin2"), colnames(soil_df))
  soil_vars <- setdiff(colnames(soil_df), soil_keys)
  df <- df %>% left_join(soil_df, by = soil_keys)
  cat(sprintf("  SoilGrids merge: %d soil_ columns added\n", length(soil_vars)))
} else {
  warning("SoilGrids CSV not found — run scripts/build_soilgrids_admin2.R Tanzania")
}

#-------------------------------------------------------------------------------
# ESPEN helminth / schistosomiasis Admin-2 predictors (iron pathway) — GUARDED
# Built by scripts/build_espen_admin2.R (STH/SCH prevalence -> iron deficiency).
#-------------------------------------------------------------------------------
espen_path <- here("data/ESPEN/Tanzania_espen_admin2.csv")
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
  message("ESPEN CSV not found — run scripts/build_espen_admin2.R Tanzania (optional domain)")
}

#-------------------------------------------------------------------------------
# MapSPAM crop-composition Admin-2 predictors (dietary availability) — GUARDED
# Built by scripts/build_mapspam_admin2.R (crop mix -> micronutrient availability).
#-------------------------------------------------------------------------------
spam_path <- here("data/MapSPAM/Tanzania_spam_admin2.csv")
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
  message("MapSPAM CSV not found — run scripts/build_mapspam_admin2.R Tanzania (optional domain)")
}

#-------------------------------------------------------------------------------
# National causal-driver Admin-2 predictors (broadcast national value) — GUARDED
#   vas_ Vitamin A supplementation      (scripts/build_vas_national.R)
#   fao_ FAOSTAT food/nutrient supply   (scripts/build_faostat_supply.R)
#   fpn_ Food Prices for Nutrition      (scripts/build_fpn_affordability.R)
#-------------------------------------------------------------------------------
for (.nd in list(list(v="vas", path=here("data/VAS/Tanzania_vas_admin2.csv")),
                 list(v="fao", path=here("data/FAOSTAT/Tanzania_fao_admin2.csv")),
                 list(v="fpn", path=here("data/FPN/Tanzania_fpn_admin2.csv")))) {
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
# WFP food prices (optional) — GUARDED. Nearest-market join + survey-year mean.
#-------------------------------------------------------------------------------
wfp_path <- here("data/food_price/wfp_food_prices_tza.csv")
if (file.exists(wfp_path)) {
  wfp <- read.csv(wfp_path)
  wfp <- wfp[-1, ]   # drop the variable-description row
  markets_subset <- wfp %>% group_by(market) %>% slice(1) %>% ungroup() %>%
    select(market, latitude, longitude)
  markets_sf <- st_as_sf(markets_subset, coords = c("longitude", "latitude"), crs = 4326)
  df <- st_zm(df, drop = TRUE, what = "ZM")
  dist_matrix <- st_distance(st_transform(df, 3857), st_transform(markets_sf, 3857))
  nearest_idx <- apply(dist_matrix, 1, which.min)
  df$nearest_market_distance_km <- apply(dist_matrix, 1, min) / 1000
  df$nearest_market_id <- markets_subset$market[nearest_idx]

  wfp <- wfp %>% mutate(year = year(as.Date(date)), usdprice = as.numeric(usdprice)) %>%
    filter(year == SURVEY_YEAR)
  ave_price <- wfp %>% group_by(market, commodity) %>%
    summarise(price = mean(usdprice, na.rm = TRUE), .groups = "drop")
  price_df <- ave_price %>% pivot_wider(names_from = commodity, values_from = price) %>%
    janitor::clean_names() %>% rename(nearest_market_id = market)
  colnames(price_df) <- paste0("wfp_", colnames(price_df))
  wfp_vars <- setdiff(colnames(price_df), "wfp_nearest_market_id")
  df <- left_join(df, price_df, by = c("nearest_market_id" = "wfp_nearest_market_id"))

  WFP_MAX_DIST_KM <- 100
  far <- !is.na(df$nearest_market_distance_km) & df$nearest_market_distance_km > WFP_MAX_DIST_KM
  if (any(far)) for (cw in intersect(wfp_vars, colnames(df))) df[[cw]][far] <- NA
  cat(sprintf("  WFP merge: %d wfp_ columns added\n", length(wfp_vars)))
} else {
  warning("WFP price CSV not found (optional) — skipping wfp_ merge.")
}

#-------------------------------------------------------------------------------
# Malaria Atlas (extract rasters at cluster points) — GUARDED by RUN_MAP
#-------------------------------------------------------------------------------
if (RUN_MAP) tryCatch({
  map_dir <- here("data/Malaria Atlas/Tanzania")

  ensure_map_blood_disorders <- function(map_dir, iso) {
    ids <- c(
      "Blood_Disorders__201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency",
      "Blood_Disorders__201201_Africa_HbC_Allele_Frequency",
      "Blood_Disorders__201201_Global_G6PDd_Allele_Frequency",
      "Blood_Disorders__201201_Global_Duffy_Negativity_Phenotype_Frequency")
    files <- paste0(ids, ".tif")
    missing <- files[!file.exists(file.path(map_dir, files))]
    if (length(missing) == 0) return(invisible(files))
    if (!requireNamespace("malariaAtlas", quietly = TRUE))
      stop("Install malariaAtlas: install.packages('malariaAtlas')")
    dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)
    shp <- malariaAtlas::getShp(ISO = iso, admin_level = "admin0")
    for (mf in missing) {
      r <- malariaAtlas::getRaster(dataset_id = sub("\\.tif$", "", mf), shp = shp)
      if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
      terra::writeRaster(r, file.path(map_dir, mf), overwrite = TRUE)
    }
    invisible(files)
  }
  blood_disorder_files <- ensure_map_blood_disorders(map_dir, iso = ISO)

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

  pts <- data.frame(lon = df$lon, lat = df$lat)
  for (i in rasters) {
    rpath <- file.path(map_dir, i)
    if (!file.exists(rpath)) next
    rast_extract <- terra::extract(terra::rast(rpath), pts, method = "bilinear")
    col_name <- gsub(".tif", "", i)
    col_name <- gsub("Malaria__|Interventions__|Blood_Disorders__", "MAP_", col_name)
    col_name <- gsub("Global_|Africa_", "", col_name)
    df$temp <- rast_extract[, 2]
    colnames(df)[ncol(df)] <- col_name
  }
  map_vars <- grep("^MAP_", colnames(df), value = TRUE)
  cat(sprintf("  Malaria Atlas: %d MAP_ columns added\n", length(map_vars)))
}, error = function(e) warning("MAP extraction failed/skipped: ", conditionMessage(e)))

#-------------------------------------------------------------------------------
# IHME admin-2 + admin-1 (fuzzy admin-name match) — GUARDED
#-------------------------------------------------------------------------------
ihme_path <- here("data/IHME/tanzania_2010_merged_IHME_data.csv")
if (file.exists(ihme_path)) {
  library(stringdist); library(stringr)
  ihme <- read.csv(ihme_path)

  clean_admin_names <- function(names) {
    names %>%
      gsub("\\s+(Municipal|Metropolis|Metropolitan|Metro)$", "", .) %>%
      gsub("[/-]", " ", .) %>% str_squish() %>% str_to_title() %>%
      gsub("\\s*\\([^)]*\\)", "", .) %>% str_trim()
  }
  df_a2   <- unique(df$Admin2);            df_a2   <- df_a2[!is.na(df_a2)]
  ihme_a2 <- unique(ihme$ihme_adm2_name);  ihme_a2 <- ihme_a2[!is.na(ihme_a2)]
  df_clean   <- data.frame(original = df_a2,   cleaned = clean_admin_names(df_a2),
                           stringsAsFactors = FALSE)
  ihme_clean <- data.frame(original = ihme_a2, cleaned = clean_admin_names(ihme_a2),
                           stringsAsFactors = FALSE)

  best_per_method <- lapply(c("jw", "cosine", "jaccard"), function(method) {
    sim <- 1 - stringdistmatrix(df_clean$cleaned, ihme_clean$cleaned, method = method)
    sim[is.na(sim)] <- -Inf                          # unmatchable pairs
    idx  <- max.col(sim, ties.method = "first")      # always a clean integer vector
    best <- sim[cbind(seq_len(nrow(sim)), idx)]
    data.frame(df_admin2 = df_clean$original,
               ihme_admin2 = ihme_clean$original[idx],
               similarity = ifelse(is.finite(best), best, NA_real_),
               method = method, stringsAsFactors = FALSE)
  })
  lookup_table <- bind_rows(best_per_method) %>%
    group_by(df_admin2) %>% slice_max(similarity, n = 1, with_ties = FALSE) %>%
    ungroup() %>% transmute(Admin2 = df_admin2, ihme_adm2_name = ihme_admin2,
                            similarity, matching_method = method)

  good <- lookup_table %>% filter(similarity >= 0.90)
  cat(sprintf("  IHME admin-2 fuzzy match: %d/%d Admin2 matched >=0.90\n",
              nrow(good), length(unique(df$Admin2))))
  df <- df %>% left_join(good, by = "Admin2") %>% left_join(ihme, by = "ihme_adm2_name")
  ihme_vars <- colnames(ihme)

  ihme_a1_path <- here("data/IHME/tanzania_2010_merged_IHME_admin1_data.csv")
  if (file.exists(ihme_a1_path) && file.exists(here("src/IHME/build_ihme_admin1.R"))) {
    source(here("src/IHME/build_ihme_admin1.R"))
    df <- merge_ihme_admin1(df, admin1_csv = ihme_a1_path, admin1_col = "Admin1")
    ihme_vars <- unique(c(ihme_vars, grep("^ihme_adm1_", colnames(df), value = TRUE)))
  }
} else {
  warning("IHME CSV not found — run the IHME clean step for Tanzania 2010.")
}

#-------------------------------------------------------------------------------
# DHS admin-1 / admin-2 indicators (surveyPrev direct / smoothed) — GUARDED
#-------------------------------------------------------------------------------
source(here("R", "data_prep.R"))
source(here("R", "food_security.R"))  # build_name_lookup() for the admin-2 name crosswalk below
dhs_adm1 <- tryCatch(load_dhs_admin1(dhs_dir = here("data", "DHS", "clean"),
                                     country = "Tanzania", year = SURVEY_YEAR),
                     error = function(e) NULL)
if (!is.null(dhs_adm1)) {
  df <- left_join(as.data.frame(df), dhs_adm1,
                  by = c("Admin1" = paste0("dhs", SURVEY_YEAR, "_DHSREGEN")))
  dhs_vars <- c(dhs_vars, grep(paste0("^dhs", SURVEY_YEAR, "_"), colnames(df), value = TRUE))
  cat(sprintf("  DHS admin-1 merge: %d dhs%d_ columns\n",
              sum(grepl(paste0("^dhs", SURVEY_YEAR, "_"), colnames(df))), SURVEY_YEAR))
} else warning("DHS admin-1 file not found for Tanzania 2010.")

dhs_adm2 <- tryCatch(load_dhs_admin2(dhs_dir = here("data", "DHS", "clean"),
                                     country = "Tanzania", year = SURVEY_YEAR),
                     error = function(e) NULL)
if (!is.null(dhs_adm2)) {
  name_lookup <- build_name_lookup(unique(df$Admin2), dhs_adm2$Admin2)
  dhs_adm2$Admin2 <- name_lookup[dhs_adm2$Admin2]
  dhs_adm2 <- dhs_adm2[!is.na(dhs_adm2$Admin2), ]
  df <- left_join(df, dhs_adm2, by = "Admin2")
  dhs_vars <- unique(c(dhs_vars, grep(paste0("^dhs", SURVEY_YEAR, "_"), colnames(df), value = TRUE)))
  cat("  DHS admin-2 merge complete\n")
} else warning("DHS admin-2 file not found for Tanzania 2010.")

#-------------------------------------------------------------------------------
# Food security (optional) — GUARDED
#-------------------------------------------------------------------------------
if (file.exists(here("R", "food_security.R"))) {
  source(here("R", "food_security.R"))
  source(here("R", "config.R"))
  cc_fsec   <- get_country_configs()[["Tanzania"]]
  hfid_path <- here("data", "HFID", "hfid_hv1.csv")
  ch_path   <- here("data", "CadreHarmonise", "cadre_harmonise_caf_ipc_dec25.xlsx")
  if (file.exists(hfid_path) || file.exists(ch_path)) {
    df <- as.data.frame(df)
    df <- merge_food_security(df, cc_fsec,
                              hfid_path = if (file.exists(hfid_path)) hfid_path else NULL,
                              ch_path   = if (file.exists(ch_path))   ch_path   else NULL)
    fsec_vars <- grep("^fsec_", colnames(df), value = TRUE)
  } else warning("No HFID / Cadre Harmonisé file — skipping fsec_ merge.")
}

#-------------------------------------------------------------------------------
# LSMS (optional, admin-1) — GUARDED. Built by src/Tanzania/TZ_LSMS_clean.R
# from the Tanzania National Panel Survey; joined on Admin1 (region).
#-------------------------------------------------------------------------------
lsms_vars <- character(0)
lsms_path <- here("data", "LSMS", "Tanzania_LSMS_clean.RDS")
if (file.exists(lsms_path)) {
  lsms <- readRDS(lsms_path)
  # Join on Admin1; light-touch (exact) — refine to fuzzy if region spellings
  # differ from GADM (cf. Ghana's Admin1_old handling).
  df <- left_join(as.data.frame(df), lsms, by = "Admin1")
  lsms_vars <- grep("^lsms_", colnames(df), value = TRUE)
  cat(sprintf("  LSMS merge: %d lsms_ columns added\n", length(lsms_vars)))
} else {
  warning("Tanzania_LSMS_clean.RDS not found (optional) — run src/Tanzania/TZ_LSMS_clean.R after downloading the NPS.")
}

#-------------------------------------------------------------------------------
# Clean: drop very-high-cardinality categoricals, geometry, time columns
#-------------------------------------------------------------------------------
check_categorical_levels <- function(df, max_levels = 100) {
  cat_cols <- sapply(df, function(x) is.factor(x) || is.character(x))
  if (!any(cat_cols)) return(character(0))
  lv <- sapply(df[, cat_cols, drop = FALSE],
               function(x) if (is.factor(x)) nlevels(x) else length(unique(x)))
  names(lv)[lv > max_levels]
}
high_card <- check_categorical_levels(df, max_levels = 100)
# Keep the id/key columns even if high-cardinality.
high_card <- setdiff(high_card, c("Admin1", "Admin2", "gw_cnum", "nearest_market_id"))
if (length(high_card)) df <- df %>% select(-all_of(high_card))

if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
# Drop any residual geometry list-columns even when df is no longer classed "sf"
# (a plain data.frame can still carry an sfc column after joins).
geom_cols <- vapply(df, function(x) inherits(x, "sfc"), logical(1))
if (any(geom_cols)) df <- df[, !geom_cols, drop = FALSE]
df <- df[, !sapply(df, function(x) inherits(x, c("POSIXct", "POSIXt")))]

df$dataid <- paste0(ID_PREFIX, seq_len(nrow(df)))

#-------------------------------------------------------------------------------
# Save merged dataset + metadata
#-------------------------------------------------------------------------------
out <- file.path(tz_dir, "Tanzania_merged_dataset.rds")
saveRDS(df, file = out)
cat(sprintf("\nSaved %s  (%d rows, %d cols)\n", out, nrow(df), ncol(df)))

metadata <- list(
  gw_vars   = gw_vars[gw_vars %in% colnames(df)],
  dhs_vars  = dhs_vars[dhs_vars %in% colnames(df)],
  ihme_vars = ihme_vars[ihme_vars %in% colnames(df)],
  map_vars  = map_vars[map_vars %in% colnames(df)],
  wfp_vars  = wfp_vars[wfp_vars %in% colnames(df)],
  gee_vars  = gee_vars[gee_vars %in% colnames(df)],
  fsec_vars = fsec_vars[fsec_vars %in% colnames(df)],
  soil_vars = soil_vars[soil_vars %in% colnames(df)],
  lsms_vars = lsms_vars[lsms_vars %in% colnames(df)],
  espen_vars = espen_vars[espen_vars %in% colnames(df)],
  spam_vars = spam_vars[spam_vars %in% colnames(df)],
  vas_vars = vas_vars[vas_vars %in% colnames(df)],
  fao_vars = fao_vars[fao_vars %in% colnames(df)],
  fpn_vars = fpn_vars[fpn_vars %in% colnames(df)]
)
dir.create(here("metadata"), showWarnings = FALSE)
saveRDS(metadata, here("metadata", "TZ_variable_categories.rds"))
cat("Saved metadata/TZ_variable_categories.rds\n")
