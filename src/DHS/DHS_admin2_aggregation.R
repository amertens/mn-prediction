# =============================================================================
# src/DHS/DHS_admin2_aggregation.R
#
# Download DHS microdata and compute admin-2 indicator prevalences using
# surveyPrev + SUMMER. Produces both direct and Fay-Herriot (BYM2) smoothed
# estimates for all available surveyPrev built-in indicators.
#
# Outputs per country-year:
#   data/DHS/clean/{Country}_{Year}_dhs_admin2_direct.rds    (long format)
#   data/DHS/clean/{Country}_{Year}_dhs_admin2_fh_bym2.rds   (long format)
#   data/DHS/clean/{Country}_{Year}_dhs_admin2_wide.rds      (wide, for merging)
#   data/DHS/clean/{Country}_{Year}_dhs_admin1_direct.rds    (admin-1 reference)
#   data/DHS/clean/{Country}_{Year}_cluster_admin_info.rds   (spatial metadata)
#
# Requirements:
#   - rdhs configured: rdhs::set_rdhs_config(email, project, config_path)
#   - INLA installed (for Fay-Herriot smoothed estimates)
#   - Internet access for DHS data download and GADM boundaries
#
# Usage:
#   source("src/DHS/DHS_admin2_aggregation.R")
#   Or run interactively — set COUNTRIES_TO_RUN to subset if needed.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(sf)
  library(geodata)
  library(surveyPrev)
  library(SUMMER)
  library(rdhs)
})

has_inla <- requireNamespace("INLA", quietly = TRUE)
if (has_inla) library(INLA) else message("INLA not installed — Fay-Herriot smoothed estimates will be skipped.")

# =============================================================================
# 0. Configuration
# =============================================================================

# Country registry: country name, ISO3 code, DHS survey years to process
COUNTRY_REGISTRY <- list(
  list(country = "Gambia",       iso3 = "GMB", years = 2019),
  list(country = "Ghana",        iso3 = "GHA", years = c(2014, 2017)),
  list(country = "Sierra Leone", iso3 = "SLE", years = 2013),
  list(country = "Malawi",       iso3 = "MWI", years = 2015)
)

# Set this to a subset if you only want to re-run specific countries
# e.g., COUNTRIES_TO_RUN <- c("Gambia")
COUNTRIES_TO_RUN <- NULL  # NULL = all countries

# Output directory
OUT_DIR <- here("data", "DHS", "clean")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Built-in surveyPrev indicators
data(surveyPrevIndicators)
ALL_INDICATORS <- surveyPrevIndicators$ID

cat(sprintf("=== DHS Admin-2 Aggregation ===\n"))
cat(sprintf("  Indicators to attempt: %d\n", length(ALL_INDICATORS)))
cat(sprintf("  INLA available: %s\n", has_inla))

# =============================================================================
# 1. Helper functions
# =============================================================================

#' Download DHS data and compute indicator values for one country-year
#'
#' @param country Country name (as used by rdhs/surveyPrev)
#' @param year DHS survey year
#' @param indicators Character vector of surveyPrev indicator IDs
#' @return list(indicator_data, indicator_log)
download_and_process_indicators <- function(country, year, indicators) {

  indicator_data <- list()
  indicator_log <- tibble(indicator = character(), ok = logical(), message = character())

  for (ind_id in indicators) {
    cat(sprintf("    [%s] ", ind_id))

    res <- tryCatch({
      dhsData <- getDHSdata(country = country, indicator = ind_id, year = year)
      dat <- getDHSindicator(dhsData, indicator = ind_id)
      n_valid <- sum(!is.na(dat$value))
      cat(sprintf("OK (n = %d)\n", n_valid))
      list(ok = TRUE, data = dat, msg = "")
    }, error = function(e) {
      cat(sprintf("SKIP (%s)\n", e$message))
      list(ok = FALSE, data = NULL, msg = e$message)
    })

    indicator_log <- bind_rows(indicator_log, tibble(
      indicator = ind_id, ok = res$ok, message = res$msg
    ))
    if (res$ok) indicator_data[[ind_id]] <- res$data
  }

  list(indicator_data = indicator_data, indicator_log = indicator_log)
}


#' Get GPS + GADM boundaries and build cluster/admin info
#'
#' @param country Country name
#' @param iso3 ISO3 country code
#' @param year DHS survey year
#' @return list(geo, poly.adm1, poly.adm2, cluster.info, admin.info1, admin.info2)
build_spatial_info <- function(country, iso3, year) {

  geo <- getDHSgeo(country = country, year = year)
  cat(sprintf("    GPS clusters: %d\n", nrow(geo)))

  poly.adm1 <- geodata::gadm(country = iso3, level = 1, path = tempdir())
  poly.adm1 <- sf::st_as_sf(poly.adm1)

  poly.adm2 <- geodata::gadm(country = iso3, level = 2, path = tempdir())
  poly.adm2 <- sf::st_as_sf(poly.adm2)

  cat(sprintf("    Admin1: %d | Admin2: %d polygons\n", nrow(poly.adm1), nrow(poly.adm2)))

  cluster.info <- clusterInfo(
    geo       = geo,
    poly.adm1 = poly.adm1,
    poly.adm2 = poly.adm2,
    by.adm1   = "NAME_1",
    by.adm2   = "NAME_2"
  )

  n_wrong <- length(cluster.info$wrong.points)
  if (n_wrong > 0) {
    cat(sprintf("    WARNING: %d clusters with invalid GPS excluded\n", n_wrong))
  }

  admin.info1 <- adminInfo(poly.adm = poly.adm1, admin = 1, by.adm = "NAME_1")
  admin.info2 <- adminInfo(poly.adm = poly.adm2, admin = 2,
                           by.adm = "NAME_2", by.adm.upper = "NAME_1")

  list(
    geo          = geo,
    poly.adm1    = poly.adm1,
    poly.adm2    = poly.adm2,
    cluster.info = cluster.info,
    admin.info1  = admin.info1,
    admin.info2  = admin.info2
  )
}


#' Compute direct estimates at admin-2 (and admin-1 for reference)
#'
#' @param indicator_data Named list of processed indicator data frames
#' @param cluster.info From build_spatial_info()
#' @return list(direct_admin2, direct_admin1) — each a list of directEST results
compute_direct_estimates <- function(indicator_data, cluster.info) {

  direct_admin2 <- list()
  direct_admin1 <- list()

  for (ind_id in names(indicator_data)) {
    dat <- indicator_data[[ind_id]]

    # Admin-2 direct
    tryCatch({
      est2 <- directEST(data = dat, cluster.info = cluster.info, admin = 2, aggregation = FALSE)
      n_est <- sum(!is.na(est2$res.admin2$direct.est))
      cat(sprintf("    [%s] Admin2 direct: %d areas\n", ind_id, n_est))
      direct_admin2[[ind_id]] <- est2
    }, error = function(e) {
      cat(sprintf("    [%s] Admin2 direct FAILED: %s\n", ind_id, e$message))
    })

    # Admin-1 direct (for reference)
    tryCatch({
      est1 <- directEST(data = dat, cluster.info = cluster.info, admin = 1)
      direct_admin1[[ind_id]] <- est1
    }, error = function(e) NULL)
  }

  list(direct_admin2 = direct_admin2, direct_admin1 = direct_admin1)
}


#' Compute Fay-Herriot BYM2 smoothed estimates at admin-2
#'
#' @param indicator_data Named list of processed indicator data frames
#' @param direct_admin2 Named list of directEST results (to detect zero-variance areas)
#' @param cluster.info From build_spatial_info()
#' @param admin.info2 From build_spatial_info()
#' @return Named list of fhModel results
compute_fh_estimates <- function(indicator_data, direct_admin2, cluster.info, admin.info2) {

  if (!has_inla) {
    cat("    Skipping Fay-Herriot (INLA not installed)\n")
    return(list())
  }

  fh_results <- list()

  for (ind_id in names(direct_admin2)) {
    dat <- indicator_data[[ind_id]]

    tryCatch({
      # Exclude areas with near-zero design variance (causes FH to fail)
      dir_est <- direct_admin2[[ind_id]]$res.admin2
      bad_admin2 <- dir_est$admin2.name.full[
        !is.na(dir_est$direct.var) & dir_est$direct.var < 1e-30
      ]
      bad_clusters <- cluster.info$data$cluster[
        cluster.info$data$admin2.name.full %in% bad_admin2
      ]

      dat_clean <- if (length(bad_clusters) > 0) {
        subset(dat, !cluster %in% bad_clusters)
      } else {
        dat
      }

      est <- fhModel(
        data         = dat_clean,
        cluster.info = cluster.info,
        admin.info   = admin.info2,
        admin        = 2,
        model        = "bym2",
        aggregation  = FALSE
      )

      n_est <- sum(!is.na(est$res.admin2$mean))
      cat(sprintf("    [%s] FH BYM2: %d areas (removed %d bad clusters)\n",
                  ind_id, n_est, length(bad_clusters)))
      fh_results[[ind_id]] <- est

    }, error = function(e) {
      cat(sprintf("    [%s] FH FAILED: %s\n", ind_id, e$message))
    })
  }

  fh_results
}


#' Compile results into long data frames
#'
#' @param results Named list of directEST or fhModel results
#' @param admin_level 1 or 2
#' @param value_col Column name for the estimate ("direct.est" or "mean")
#' @return tibble in long format
compile_long <- function(results, admin_level = 2, value_col = "direct.est") {
  rows <- list()
  for (ind_id in names(results)) {
    res_name <- paste0("res.admin", admin_level)
    if (!is.null(results[[ind_id]][[res_name]])) {
      df <- results[[ind_id]][[res_name]]
      df$indicator <- ind_id
      rows[[ind_id]] <- df
    }
  }
  if (length(rows) == 0) return(tibble())
  bind_rows(rows)
}


#' Pivot long results to wide format for merging into the pipeline
#'
#' Uses Fay-Herriot smoothed mean where available, falls back to direct estimate.
#' Column naming: {prefix}{indicator}_adm2
#'
#' @param direct_long Long-format direct estimates
#' @param fh_long Long-format FH estimates (can be empty tibble)
#' @param prefix Column name prefix (e.g., "dhs2019_")
#' @return data.frame with Admin2 + prefixed indicator columns
build_wide_output <- function(direct_long, fh_long, prefix) {

  # Start with direct estimates
  if (nrow(direct_long) == 0) return(tibble())

  direct_wide <- direct_long %>%
    select(admin2.name.full, indicator, direct.est) %>%
    rename(est = direct.est)

  # Override with FH smoothed where available
  if (nrow(fh_long) > 0) {
    fh_wide <- fh_long %>%
      select(admin2.name.full, indicator, mean) %>%
      rename(est = mean)

    # Use FH where available, direct as fallback
    direct_wide <- direct_wide %>%
      left_join(fh_wide, by = c("admin2.name.full", "indicator"), suffix = c("_direct", "_fh")) %>%
      mutate(est = coalesce(est_fh, est_direct)) %>%
      select(admin2.name.full, indicator, est)
  }

  # Pivot to wide
  wide <- direct_wide %>%
    pivot_wider(
      names_from   = indicator,
      values_from  = est,
      names_prefix = prefix
    )

  # Add _adm2 suffix to all indicator columns (not admin2.name.full)
  wide <- wide %>%
    rename_with(~ paste0(.x, "_adm2"), -admin2.name.full)

  # Extract Admin2 name from "Admin1_Admin2" format
  wide$Admin2 <- sub("^[^_]+_", "", wide$admin2.name.full)
  wide$admin2.name.full <- NULL

  # Move Admin2 to first column
  wide <- wide %>% select(Admin2, everything())

  wide
}


# =============================================================================
# 2. Main loop over countries
# =============================================================================

for (entry in COUNTRY_REGISTRY) {

  if (!is.null(COUNTRIES_TO_RUN) && !(entry$country %in% COUNTRIES_TO_RUN)) {
    next
  }

  for (yr in entry$years) {

    cat(sprintf("\n%s\n", strrep("=", 70)))
    cat(sprintf("  %s %d\n", entry$country, yr))
    cat(sprintf("%s\n\n", strrep("=", 70)))

    prefix <- paste0("dhs", yr, "_")

    # --- Step 1: Download and process indicators ---
    cat("  [1/5] Downloading DHS microdata and processing indicators...\n")
    ind_res <- tryCatch(
      download_and_process_indicators(entry$country, yr, ALL_INDICATORS),
      error = function(e) {
        cat(sprintf("  FATAL: Could not download data for %s %d: %s\n",
                    entry$country, yr, e$message))
        NULL
      }
    )
    if (is.null(ind_res) || length(ind_res$indicator_data) == 0) {
      cat(sprintf("  Skipping %s %d — no indicators processed\n", entry$country, yr))
      next
    }

    cat(sprintf("  Processed: %d / %d indicators\n",
                length(ind_res$indicator_data), length(ALL_INDICATORS)))

    # --- Step 2: Build spatial info ---
    cat("  [2/5] Building spatial info (GPS + GADM boundaries)...\n")
    spatial <- tryCatch(
      build_spatial_info(entry$country, entry$iso3, yr),
      error = function(e) {
        cat(sprintf("  FATAL: Spatial info failed for %s %d: %s\n",
                    entry$country, yr, e$message))
        NULL
      }
    )
    if (is.null(spatial)) next

    # --- Step 3: Direct estimates ---
    cat("  [3/5] Computing direct estimates (admin-2 + admin-1)...\n")
    direct_res <- compute_direct_estimates(
      ind_res$indicator_data, spatial$cluster.info
    )

    # --- Step 4: Fay-Herriot smoothed estimates ---
    cat("  [4/5] Computing Fay-Herriot BYM2 smoothed estimates...\n")
    fh_res <- compute_fh_estimates(
      ind_res$indicator_data,
      direct_res$direct_admin2,
      spatial$cluster.info,
      spatial$admin.info2
    )

    # --- Step 5: Compile and save ---
    cat("  [5/5] Compiling and saving results...\n")

    tryCatch({
      # Long format tables
      direct_long <- compile_long(direct_res$direct_admin2, admin_level = 2, value_col = "direct.est")
      fh_long     <- compile_long(fh_res, admin_level = 2, value_col = "mean")
      admin1_long <- compile_long(direct_res$direct_admin1, admin_level = 1, value_col = "direct.est")

      # Wide format for pipeline merging (FH preferred, direct fallback)
      wide <- build_wide_output(direct_long, fh_long, prefix)

      # Country-year file prefix
      file_prefix <- paste0(entry$country, "_", yr)

      # Save direct (long)
      if (nrow(direct_long) > 0) {
        f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin2_direct.rds"))
        saveRDS(direct_long, f)
        cat(sprintf("    Saved: %s (%d rows, %d indicators)\n",
                    basename(f), nrow(direct_long), length(unique(direct_long$indicator))))
      }

      # Save FH smoothed (long)
      if (nrow(fh_long) > 0) {
        f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin2_fh_bym2.rds"))
        saveRDS(fh_long, f)
        cat(sprintf("    Saved: %s (%d rows, %d indicators)\n",
                    basename(f), nrow(fh_long), length(unique(fh_long$indicator))))
      }

      # Save wide (for merging)
      if (nrow(wide) > 0) {
        f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin2_wide.rds"))
        saveRDS(wide, f)
        cat(sprintf("    Saved: %s (%d admin2 areas, %d columns)\n",
                    basename(f), nrow(wide), ncol(wide)))
      }

      # Save admin-1 direct (for reference / existing pipeline)
      if (nrow(admin1_long) > 0) {
        f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin1_direct.rds"))
        saveRDS(admin1_long, f)
        cat(sprintf("    Saved: %s\n", basename(f)))
      }

      # Save spatial metadata
      f <- file.path(OUT_DIR, paste0(file_prefix, "_cluster_admin_info.rds"))
      saveRDS(list(
        cluster.info = spatial$cluster.info,
        admin.info1  = spatial$admin.info1,
        admin.info2  = spatial$admin.info2,
        poly.adm1    = spatial$poly.adm1,
        poly.adm2    = spatial$poly.adm2
      ), f)
      cat(sprintf("    Saved: %s\n", basename(f)))

      # Save processing log
      f <- file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin2_log.rds"))
      saveRDS(ind_res$indicator_log, f)

      cat(sprintf("\n  %s %d complete: %d direct, %d FH indicators\n",
                  entry$country, yr,
                  length(direct_res$direct_admin2),
                  length(fh_res)))

    }, error = function(e) {
      cat(sprintf("  ERROR in compile/save for %s %d: %s\n",
                  entry$country, yr, e$message))
      cat("  Attempting to save raw results as fallback...\n")
      file_prefix <- paste0(entry$country, "_", yr)
      tryCatch({
        saveRDS(direct_res, file.path(OUT_DIR, paste0(file_prefix, "_dhs_admin2_raw_results.rds")))
        cat("    Saved raw results fallback.\n")
      }, error = function(e2) {
        cat(sprintf("    Fallback save also failed: %s\n", e2$message))
      })
    })
  }
}

cat("\n=== All countries complete ===\n")
