# Utility functions to support DHS data cleaning, merging, and analysis
# for multiple countries using a consistent interface.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(janitor)
  library(terra)
  library(here)
})

source(paste0(here::here(), "/src/0-functions.R"))
source(paste0(here::here(), "/src/0-SL-setup.R"))
source(paste0(here::here(), "/src/DHS/DHS_functions.R"))
source(paste0(here::here(), "/src/DHS/DHS_variable_recode.R"))

#' Read and clean a DHS country file
#'
#' @param dhs_path Path to the DHS RDS file.
#'
clean_country_dhs <- function(dhs_path) {
  if (!file.exists(dhs_path)) {
    stop(sprintf("DHS file not found at %s", dhs_path))
  }

  dhs_object <- readRDS(dhs_path)
  dhs_data <- if (is.list(dhs_object)) dhs_object[[1]] else dhs_object

  df <- clean_DHS(dhs_data)
  df <- df %>% select(!ends_with(".y"))
  colnames(df) <- gsub("\\.x", "", colnames(df))
  df
}

#' Attach compiled GPS coordinates from a shared DHS GPS file
#'
#' @param df Cleaned DHS dataframe containing `cluster`
#' @param gps_path Path to the compiled GPS RDS file
#' @param country Country name to filter
#' @param year Survey year used to subset the GPS data
attach_compiled_gps <- function(df, gps_path, country, year) {
  if (!file.exists(gps_path)) {
    warning(sprintf("GPS file not found at %s; returning data without GPS merge", gps_path))
    return(df)
  }

  geo <- readRDS(gps_path) %>% dplyr::filter(country == !!country, year == !!year)
  geo_df <- as.data.frame(geo) %>%
    select(DHSCLUST, LATNUM, LONGNUM) %>%
    rename(cluster = DHSCLUST, latitude = LATNUM, longitude = LONGNUM)

  left_join(df, geo_df, by = "cluster")
}

#' Merge cluster-level GEE predictors
attach_gee_features <- function(df, gee_path) {
  if (!file.exists(gee_path)) {
    warning(sprintf("GEE file not found at %s; returning data without GEE merge", gee_path))
    return(df)
  }

  gee_df <- read.csv(gee_path) %>%
    rename(cluster = DHSCLUST)

  left_join(df, gee_df, by = "cluster")
}

#' Summarise price data into a wide format and merge
attach_food_prices <- function(df, price_path, market_id_col = "nearest_market_id") {
  if (!file.exists(price_path)) {
    warning(sprintf("Food price file not found at %s; returning data without price merge", price_path))
    return(df)
  }

  wfp <- read.csv(price_path)

  price_df <- wfp %>%
    dplyr::mutate(
      usdprice = suppressWarnings(as.numeric(usdprice)),
      priceflag = factor(priceflag, levels = c("actual", "actual,aggregate", "aggregate")),
      pricetype = factor(pricetype, levels = c("Retail", "Wholesale")),
      unit = factor(unit),
      year = lubridate::year(date)
    ) %>%
    group_by(market, commodity, currency, priceflag, pricetype, unit) %>%
    summarise(price = mean(usdprice, na.rm = TRUE), .groups = "drop") %>%
    group_by(market, commodity, currency) %>%
    arrange(market, commodity, currency, priceflag, pricetype, unit) %>%
    slice(1) %>%
    ungroup() %>%
    select(market, commodity, price) %>%
    pivot_wider(names_from = commodity, values_from = price) %>%
    clean_names() %>%
    rename(!!market_id_col := market)

  left_join(df, price_df, by = market_id_col)
}

#' Extract raster-based predictors from Malaria Atlas layers
#'
#' @param df Dataframe containing `latitude` and `longitude` columns
#' @param raster_dir Directory containing country-specific rasters
#' @param country_prefix Prefix used by the raster filenames (e.g., "GHA")
attach_malaria_atlas <- function(df, raster_dir, country_prefix) {
  if (is.null(raster_dir) || is.null(country_prefix)) {
    return(df)
  }

  if (!dir.exists(raster_dir)) {
    warning(sprintf("Raster directory %s not found; returning data without malaria atlas features", raster_dir))
    return(df)
  }

  raster_files <- list.files(raster_dir, pattern = paste0("^", country_prefix, ".*\\.tif$"), full.names = TRUE)
  if (length(raster_files) == 0) {
    warning(sprintf("No raster files with prefix %s found in %s", country_prefix, raster_dir))
    return(df)
  }

  coords <- df %>% filter(!is.na(longitude), !is.na(latitude))
  if (nrow(coords) == 0) {
    warning("No coordinates available to extract malaria atlas features")
    return(df)
  }

  pts <- terra::vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")

  extracted <- lapply(raster_files, function(path) {
    layer <- rast(path)
    values <- terra::extract(layer, pts)[, 2]
    tibble(value = values, layer = tools::file_path_sans_ext(basename(path)))
  })

  extracted_df <- bind_cols(coords[ , c("cluster")],
                            setNames(
                              as_tibble(do.call(cbind, lapply(extracted, `[[`, "value"))),
                              sapply(extracted, function(x) unique(x$layer))
                            ))

  df %>% left_join(extracted_df, by = "cluster")
}

#' Full workflow wrapper
run_country_workflow <- function(country, iso2, iso3, year,
                                 dhs_path,
                                 gee_path = NULL,
                                 price_path = NULL,
                                 gps_path = here("data/DHS_compiled_gps_data.rds"),
                                 raster_dir = here("data/Malaria Atlas"),
                                 malaria_prefix = NULL,
                                 output_dir = here("data/DHS/clean")) {
  df <- clean_country_dhs(dhs_path)
  df <- attach_compiled_gps(df, gps_path, country, year)

  if (!is.null(gee_path)) {
    df <- attach_gee_features(df, gee_path)
  }

  if (!is.null(price_path)) {
    df <- attach_food_prices(df, price_path)
  }

  if (!is.null(malaria_prefix)) {
    df <- attach_malaria_atlas(df, raster_dir, malaria_prefix)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  outfile <- file.path(output_dir, sprintf("dhs_%s_%s_clean.RDS", gsub(" ", "", country), year))
  saveRDS(df, outfile)
  message(sprintf("Saved cleaned data to %s", outfile))

  invisible(df)
}
