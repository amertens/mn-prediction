# =============================================================================
# scripts/build_soilgrids_admin2.R
#
# Per-country Admin-2 zonal-mean extraction of SoilGrids/ISRIC properties:
#   soil_ph, soil_organic_carbon, soil_nitrogen, soil_clay, soil_sand,
#   soil_silt, soil_cec — at 0–5 cm depth, mean values.
#
# Downloads the SoilGrids global VRT (cloud-optimised GeoTIFF) via
# /vsicurl/, crops to country bbox, reprojects to WGS84, zonal-mean per
# Admin-2 polygon. Output: data/SoilGrids/<Country>_soilgrids_admin2.csv.
#
# Cache is shared across countries (downloads are bbox-cropped per country
# but the global VRTs are queried once per property).
#
# Run:
#   Rscript scripts/build_soilgrids_admin2.R               # all 4 countries
#   Rscript scripts/build_soilgrids_admin2.R Gambia        # one country
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(here)
})

source(here("R", "data_prep.R"))     # load_gadm_cached()
source(here("R", "external_data.R")) # extract_soilgrids()

countries <- list(
  Gambia      = list(iso2 = "GM"),
  Ghana       = list(iso2 = "GH"),
  SierraLeone = list(iso2 = "SLE"),
  Malawi      = list(iso2 = "MW"),
  # Tanzania: use "TZA" so the GADM polygons here match those the merge step
  # pulls via load_gadm_cached("TZA"), making the exact (Admin1, Admin2) soil
  # join in 2_GW_Tanzania_data_merge.R line up without fuzzy matching.
  Tanzania    = list(iso2 = "TZA")
)

cache_dir <- here("data", "external_cache")
out_dir   <- here("data", "SoilGrids")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
sel  <- if (length(args) == 0) names(countries) else args

for (cn in sel) {
  cc <- countries[[cn]]
  if (is.null(cc)) { message("unknown country: ", cn); next }

  message("\n=== ", cn, " ===")
  poly <- load_gadm_cached(cc$iso2, level = 2)
  poly <- sf::st_as_sf(poly) |>
    dplyr::select(NAME_1, NAME_2) |>
    dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2)
  poly$Admin1 <- trimws(as.character(poly$Admin1))
  poly$Admin2 <- trimws(as.character(poly$Admin2))

  res <- tryCatch(
    extract_soilgrids(poly, cache_dir = cache_dir),
    error = function(e) {
      message("  extract_soilgrids failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res) || nrow(res) == 0) {
    message("  no result; skipping ", cn); next
  }

  # extract_soilgrids returns Admin2 + soil_* cols, one row per polygon (same
  # order as admin2_sf). Attach Admin1 column-wise — merging by Admin2 alone
  # would fan out for countries where Admin2 names are not unique without
  # their parent Admin1 (e.g. Malawi has "Lake Malawi" in 8 different districts).
  stopifnot(nrow(res) == nrow(poly))
  res <- cbind(data.frame(Admin1 = poly$Admin1, stringsAsFactors = FALSE),
               res)
  res <- res[, c("Admin1", "Admin2",
                 setdiff(colnames(res), c("Admin1", "Admin2")))]

  out_path <- file.path(out_dir, sprintf("%s_soilgrids_admin2.csv", cn))
  readr::write_csv(res, out_path)
  message(sprintf("  ✓ wrote %s  (%d admin2 × %d cols)",
                  out_path, nrow(res), ncol(res) - 2L))
}
