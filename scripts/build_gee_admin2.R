# =============================================================================
# scripts/build_gee_admin2.R
#
# One-time / refresh script: builds a per-country wide CSV of Admin-2 zonal-
# mean values for every .tif under data/<Country>_GEE_rasters/ plus shared
# rasters under data/GEE/. Writes:
#
#   data/GEE/<Country>_<year>_admin2_gee.csv
#
# These CSVs are consumed by each country's 2_GW_<Country>_data_merge.R via
# a left-join on Admin2, and the resulting columns (prefixed `gee_a2_`) flow
# into the GEE proxy-predictor domain automatically.
#
# Year-of-survey filtering: rasters whose filename embeds a 4-digit year are
# kept only when that year matches the country's survey year (or the closest
# available, configured below). Time-invariant rasters (no year in filename)
# are always kept.
#
# Run:
#   Rscript scripts/build_gee_admin2.R               # all 4 countries
#   Rscript scripts/build_gee_admin2.R Gambia        # one country
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(readr)
  library(here)
})

source(here("src/GEE/extract_gee_admin2.R"))
source(here("R", "data_prep.R"))  # provides load_gadm_cached()

# ── Country configuration ────────────────────────────────────────────────────
# `target_year` is the year-of-survey raster to prefer when a family has
# multiple year-specific files. If that exact year isn't available, the
# orchestrator falls back to the nearest available year for that family.
countries <- list(
  Gambia = list(
    iso2 = "GM",  folder = "Gambia_GEE_rasters",
    target_year = 2018L, tokens = c("Gambia")
  ),
  Ghana = list(
    iso2 = "GH",  folder = "Ghana_GEE rasters",
    target_year = 2017L, tokens = c("Ghana")
  ),
  SierraLeone = list(
    iso2 = "SLE", folder = "Sierra_Leone_GEE_rasters",
    target_year = 2013L,
    tokens = c("Sierra_Leone", "Sierra Leone", "SierraLeone")
  ),
  Malawi = list(
    iso2 = "MW",  folder = "Malawi_GEE_rasters",
    # Survey fieldwork ran Dec 2015 - Feb 2016 (majority Jan-Feb 2016).
    # Target 2016; nearest-year fallback yields 2015 for families whose
    # rasters stop at 2015 (FLDAS, LST_Night, AerosolOptical, etc.).
    target_year = 2016L, tokens = c("Malawi")
  )
)

# ── Helpers ──────────────────────────────────────────────────────────────────
extract_year <- function(fname) {
  m <- regmatches(fname, regexpr("(19|20)\\d{2}", fname))
  if (length(m) == 0) NA_integer_ else as.integer(m)
}

#' Strip year + country tokens to identify the raster "family" (so we can
#' pick the best year per family).
family_key <- function(fname, tokens) {
  key <- sub("\\.tif$", "", basename(fname), ignore.case = TRUE)
  key <- gsub("(19|20)\\d{2}", "", key)
  for (tok in tokens) {
    key <- gsub(paste0("[_ ]*", gsub(" ", "[_ ]", tok), "[_ ]*"),
                "_", key, ignore.case = TRUE)
  }
  key <- gsub("__+", "_", key)
  key <- gsub("^_|_$", "", key)
  tolower(key)
}

#' For each raster family, pick the single best file for the survey year.
#' Priority order:
#'   1. Year-specific file whose year exactly matches target_year.
#'   2. Year-specific file with the nearest year.
#'   3. No-year file (e.g. SoilIron_Gambia.tif) — used only if no year-
#'      specific exists for the family.
#'
#' This avoids duplicate columns when a family has both a no-year canonical
#' file and year-specific variants (e.g. SoilIron has 3 files: no-year + 2017
#' + 2018; we keep just the year matching the survey).
select_rasters_for_year <- function(files, target_year, tokens) {
  if (length(files) == 0) return(character(0))
  yrs <- vapply(basename(files), extract_year, integer(1))
  fams <- vapply(files, family_key, character(1), tokens = tokens)
  selected <- character(0)
  for (fam in unique(fams)) {
    idx <- which(fams == fam)
    fam_files <- files[idx]
    fam_years <- yrs[idx]
    year_idx <- which(!is.na(fam_years))
    if (length(year_idx) > 0) {
      diffs <- abs(fam_years[year_idx] - target_year)
      pick <- year_idx[which.min(diffs)]
      selected <- c(selected, fam_files[pick])
    } else {
      # All files for this family are no-year — keep them all.
      selected <- c(selected, fam_files)
    }
  }
  selected
}

# ── Main loop ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
sel <- if (length(args) == 0) names(countries) else args

global_tifs <- list.files(here("data/GEE"), pattern = "\\.tif$",
                          full.names = TRUE)

for (cn in sel) {
  cc <- countries[[cn]]
  if (is.null(cc)) { message("unknown country: ", cn); next }
  ctry_folder <- here("data", cc$folder)
  if (!dir.exists(ctry_folder)) {
    message("skip ", cn, ": missing folder ", ctry_folder); next
  }

  message("\n=== ", cn, " (target_year = ", cc$target_year, ") ===")

  poly <- load_gadm_cached(cc$iso2, level = 2)
  poly <- sf::st_as_sf(poly) |>
    dplyr::select(NAME_1, NAME_2) |>
    dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)
  poly$Admin2 <- trimws(poly$Admin2)

  country_rasters <- list.files(ctry_folder, pattern = "\\.tif$",
                                full.names = TRUE)
  country_rasters <- select_rasters_for_year(country_rasters,
                                              cc$target_year, cc$tokens)

  all_rasters <- c(country_rasters, global_tifs)
  message(sprintf("  %d rasters selected (%d country + %d global)",
                  length(all_rasters), length(country_rasters),
                  length(global_tifs)))

  res <- extract_gee_admin2(poly, all_rasters,
                             country_tokens = cc$tokens, verbose = TRUE)

  out_path <- here("data/GEE",
                   sprintf("%s_%d_admin2_gee.csv", cn, cc$target_year))
  readr::write_csv(res, out_path)
  message(sprintf("  ✓ wrote %s  (%d admin2 × %d cols)",
                  out_path, nrow(res), ncol(res) - 1L))
}
