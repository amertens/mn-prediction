# =============================================================================
# scripts/download_external_predictors.R
#
# Downloads and caches ALL external predictor data for the 4 study countries.
# Run this once to populate data/external_cache/, then the pipeline can
# use cached files without internet access.
#
# Usage:
#   source("scripts/download_external_predictors.R")
#
# Prerequisites:
#   - .Renviron with NASA_EARTHDATA_TOKEN (for nightlights)
#   - Optional: ACLED_EMAIL_ADDRESS + ACLED_ACCESS_KEY, IPC_API_KEY
#   - data/external_cache/hvstat_africa_data_v1.0.csv (manual download from Dryad)
#   - R packages: see install block below
# =============================================================================

library(here)

# Source the extraction functions
source(here("R", "external_data.R"))
source(here("R", "config.R"))

cache_dir <- here("data", "external_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Countries and survey years (from survey reports)
# Keys must match get_country_configs() names exactly
countries_years <- list(
  Gambia      = 2018L,  # GMNS fieldwork: 13 Mar – 4 May 2018
  Ghana       = 2017L,  # GMS fieldwork: 27 Apr – 9 Jun 2017
  SierraLeone = 2013L,  # SLMS fieldwork: 11 Nov – 2 Dec 2013
  Malawi      = 2016L   # MNS 2015-16 (exact dates TBD)
)

all_configs <- get_country_configs()

# =============================================================================
# PART 1: Run full extraction for all 4 countries
# =============================================================================

results <- list()
for (cn in names(countries_years)) {
  yr <- countries_years[[cn]]
  cc <- all_configs[[cn]]

  if (is.null(cc)) {
    cat(sprintf("\n[SKIP] No config found for '%s'\n", cn))
    next
  }

  cat(sprintf("\n\n%s\n", paste(rep("=", 70), collapse = "")))
  cat(sprintf("  %s — survey year %d\n", cn, yr))
  cat(sprintf("%s\n\n", paste(rep("=", 70), collapse = "")))

  res <- tryCatch(
    extract_all_external(cc, yr, cache_dir),
    error = function(e) {
      cat(sprintf("\n[ERROR] %s extraction failed: %s\n", cn, e$message))
      NULL
    }
  )

  if (!is.null(res)) {
    # Save per-country result
    out_file <- file.path(cache_dir,
                          sprintf("%s_external_predictors.rds",
                                  tolower(gsub(" ", "_", cn))))
    saveRDS(res, out_file)

    pred_cols <- setdiff(colnames(res), c("Admin1", "Admin2"))
    n_with_data <- sum(colSums(!is.na(res[, pred_cols, drop = FALSE])) > 0)
    cat(sprintf("\n[SAVED] %s: %d Admin-2 x %d predictors (%d with data) -> %s\n",
                cn, nrow(res), length(pred_cols), n_with_data, basename(out_file)))
    results[[cn]] <- res
  }
}


# =============================================================================
# PART 2: Debug failing sources
# =============================================================================

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  DEBUGGING FAILING SOURCES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Use Ghana as test case (best data coverage)
cc_test <- all_configs[["Ghana"]]
admin2_sf <- sf::st_as_sf(geodata::gadm("GHA", level = 2, path = file.path(cache_dir, "gadm")))
admin2_sf$Admin2 <- admin2_sf$NAME_2
admin2_sf$Admin1 <- admin2_sf$NAME_1


# ── Debug 1: VIIRS Nighttime Lights ──────────────────────────────────────────
cat("\n--- DEBUG: VIIRS Nighttime Lights ---\n")

nasa_token <- Sys.getenv("NASA_EARTHDATA_TOKEN", unset = "")
cat(sprintf("  NASA token length: %d chars\n", nchar(nasa_token)))

if (nchar(nasa_token) == 0) {
  cat("  ISSUE: No NASA_EARTHDATA_TOKEN in environment.\n")
  cat("  FIX: Add to .Renviron or run: Sys.setenv(NASA_EARTHDATA_TOKEN = 'your_token')\n")
} else {
  cat("  Token present. Testing blackmarbler...\n")

  if (!requireNamespace("blackmarbler", quietly = TRUE)) {
    cat("  ISSUE: blackmarbler package not installed.\n")
    cat("  FIX: install.packages('blackmarbler')\n")
  } else {
    cat(sprintf("  blackmarbler version: %s\n",
                as.character(packageVersion("blackmarbler"))))

    # Try downloading a small test
    ntl_test <- tryCatch({
      ntl_dir <- file.path(cache_dir, "nightlights", "GHA")
      dir.create(ntl_dir, showWarnings = FALSE, recursive = TRUE)

      cat("  Attempting bm_raster download for Ghana 2017...\n")
      r <- blackmarbler::bm_raster(
        roi_sf     = admin2_sf,
        product_id = "VNP46A4",
        date       = 2017,
        bearer     = nasa_token
      )
      cat(sprintf("  SUCCESS! Class: %s, nlyr: %d\n", class(r)[1], terra::nlyr(r)))

      # Extract and show sample values
      vals <- exactextractr::exact_extract(r, admin2_sf[1:3, ], fun = "mean")
      cat(sprintf("  Sample values (first 3 Admin2): %s\n",
                  paste(round(vals, 2), collapse = ", ")))
      "OK"
    }, error = function(e) {
      cat(sprintf("  FAILED: %s\n", e$message))
      cat("  Possible causes:\n")
      cat("    - Token expired (check https://urs.earthdata.nasa.gov/profile)\n")
      cat("    - Network/firewall blocking NASA LAADS DAAC\n")
      cat("    - blackmarbler API change\n")
      e$message
    })
  }
}


# ── Debug 2: SoilGrids ──────────────────────────────────────────────────────
cat("\n--- DEBUG: SoilGrids (ISRIC REST API) ---\n")

# Test the WCS endpoint first
cat("  Testing ISRIC WCS 2.0 endpoint...\n")
wcs_test <- tryCatch({
  # Try a small WCS GetCoverage request for pH
  url <- paste0(
    "https://maps.isric.org/mapserv?map=/map/phh2o.map",
    "&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=phh2o_0-5cm_mean",
    "&FORMAT=image/tiff",
    "&SUBSET=long(-3.3,1.2)&SUBSET=lat(4.7,11.2)",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  )
  cat(sprintf("  URL: %s\n", substr(url, 1, 80)))
  resp <- httr::GET(url, httr::timeout(30))
  cat(sprintf("  HTTP status: %d\n", httr::status_code(resp)))
  cat(sprintf("  Content type: %s\n", httr::headers(resp)$`content-type`))
  cat(sprintf("  Content length: %d bytes\n", length(httr::content(resp, "raw"))))

  if (httr::status_code(resp) == 200 &&
      grepl("tiff", httr::headers(resp)$`content-type`, ignore.case = TRUE)) {
    # Save and try to read
    tmp <- tempfile(fileext = ".tif")
    writeBin(httr::content(resp, "raw"), tmp)
    r <- terra::rast(tmp)
    cat(sprintf("  SUCCESS! Raster: %d x %d, CRS: %s\n",
                terra::nrow(r), terra::ncol(r), terra::crs(r, describe = TRUE)$code))
    unlink(tmp)
    "WCS_OK"
  } else {
    cat("  WCS returned non-TIFF response. Trying REST API...\n")
    "WCS_FAIL"
  }
}, error = function(e) {
  cat(sprintf("  WCS FAILED: %s\n", e$message))
  "WCS_ERROR"
})

# Test REST API point query
cat("\n  Testing ISRIC REST API point query...\n")
rest_test <- tryCatch({
  # Query a single point (Accra, Ghana)
  url <- "https://rest.isric.org/soilgrids/v2.0/properties/query?lon=-0.19&lat=5.55&property=phh2o&depth=0-5cm&value=mean"
  cat(sprintf("  URL: %s\n", substr(url, 1, 80)))
  resp <- httr::GET(url, httr::timeout(30))
  cat(sprintf("  HTTP status: %d\n", httr::status_code(resp)))

  if (httr::status_code(resp) == 200) {
    body <- httr::content(resp, "parsed")
    cat(sprintf("  Response type: %s\n", class(body)))
    # Try to extract value
    if (is.list(body) && "properties" %in% names(body)) {
      layers <- body$properties$layers
      if (length(layers) > 0) {
        val <- layers[[1]]$depths[[1]]$values$mean
        cat(sprintf("  pH value at Accra: %s\n", val))
        "REST_OK"
      } else {
        cat("  No layers in response\n")
        "REST_EMPTY"
      }
    } else {
      cat(sprintf("  Unexpected response structure: %s\n",
                  paste(names(body), collapse = ", ")))
      "REST_UNEXPECTED"
    }
  } else {
    body <- tryCatch(httr::content(resp, "text"), error = function(e) "")
    cat(sprintf("  Error body: %s\n", substr(body, 1, 200)))
    "REST_HTTP_ERROR"
  }
}, error = function(e) {
  cat(sprintf("  REST FAILED: %s\n", e$message))
  "REST_ERROR"
})


# ── Debug 3: Global Data Lab ────────────────────────────────────────────────
cat("\n--- DEBUG: Global Data Lab (SHDI) ---\n")

# Check if the user has a pre-downloaded CSV
gdl_paths <- c(
  file.path(cache_dir, "gdl_shdi_full.csv"),
  file.path(cache_dir, "GDL-Sub-national-HDI-data.csv"),
  here("data", "raw", "gdl_shdi_full.csv")
)
gdl_found <- FALSE
for (p in gdl_paths) {
  if (file.exists(p)) {
    cat(sprintf("  Found local GDL file: %s\n", p))
    gdl_found <- TRUE
    break
  }
}

if (!gdl_found) {
  cat("  No local GDL CSV found. Trying web download...\n")

  # Try multiple possible download URLs
  gdl_urls <- c(
    "https://globaldatalab.org/assets/2023/09/SHDI-SGDI-Total%207.0.csv",
    "https://globaldatalab.org/assets/2024/01/SHDI-SGDI-Total%208.0.csv",
    "https://globaldatalab.org/shdi/download_file/"
  )

  for (url in gdl_urls) {
    cat(sprintf("  Trying: %s\n", substr(url, 1, 70)))
    dl_test <- tryCatch({
      resp <- httr::GET(url, httr::timeout(30),
                        httr::user_agent("Mozilla/5.0 (academic research)"))
      cat(sprintf("    HTTP status: %d, Content-Type: %s, Size: %d\n",
                  httr::status_code(resp),
                  httr::headers(resp)$`content-type` %||% "unknown",
                  length(httr::content(resp, "raw"))))

      if (httr::status_code(resp) == 200 &&
          length(httr::content(resp, "raw")) > 1000) {
        # Save it
        gdl_save <- file.path(cache_dir, "gdl_shdi_full.csv")
        writeBin(httr::content(resp, "raw"), gdl_save)
        cat(sprintf("    SAVED to %s\n", gdl_save))

        # Try to parse
        df <- read.csv(gdl_save, nrows = 5)
        cat(sprintf("    Columns: %s\n", paste(head(colnames(df), 8), collapse = ", ")))
        cat(sprintf("    Rows (preview): %d\n", nrow(df)))
        "DOWNLOAD_OK"
      } else {
        "DOWNLOAD_FAIL"
      }
    }, error = function(e) {
      cat(sprintf("    FAILED: %s\n", e$message))
      "DOWNLOAD_ERROR"
    })

    if (dl_test == "DOWNLOAD_OK") break
  }

  if (dl_test != "DOWNLOAD_OK") {
    cat("\n  MANUAL DOWNLOAD NEEDED:\n")
    cat("    1. Go to https://globaldatalab.org/shdi/\n")
    cat("    2. Click 'Download' or 'Table' -> export as CSV\n")
    cat("    3. Save as: data/external_cache/gdl_shdi_full.csv\n")
  }
}


# ── Debug 4: WFP HungerMap ──────────────────────────────────────────────────
cat("\n--- DEBUG: WFP HungerMap LIVE ---\n")

# Test various public endpoints
wfp_endpoints <- c(
  country_data  = "https://static.hungermapdata.org/hunger-map-data/v2/adm0/GHA.json",
  adm1_data     = "https://static.hungermapdata.org/hunger-map-data/v2/adm1/GHA.json",
  global_status = "https://static.hungermapdata.org/hunger-map-data/v2/global.json",
  country_iso   = "https://api.hungermapdata.org/v2/info/country?iso3=GHA",
  country_adm0  = "https://api.hungermapdata.org/v2/adm0/GHA"
)

for (ep_name in names(wfp_endpoints)) {
  url <- wfp_endpoints[[ep_name]]
  cat(sprintf("  [%s] %s\n", ep_name, substr(url, 1, 70)))

  resp <- tryCatch({
    r <- httr::GET(url, httr::timeout(15),
                   httr::user_agent("Mozilla/5.0 (academic research)"))
    cat(sprintf("    Status: %d, Type: %s, Size: %d bytes\n",
                httr::status_code(r),
                httr::headers(r)$`content-type` %||% "unknown",
                length(httr::content(r, "raw"))))

    if (httr::status_code(r) == 200) {
      body <- tryCatch(httr::content(r, "parsed"), error = function(e) NULL)
      if (is.list(body)) {
        cat(sprintf("    Top-level keys: %s\n",
                    paste(head(names(body), 10), collapse = ", ")))

        # Look for FCS/rCSI data
        if (any(grepl("fcs|rcsi|food", tolower(names(body))))) {
          cat("    >>> Found food security indicators!\n")
        }
        if ("body" %in% names(body) && is.list(body$body)) {
          cat(sprintf("    body keys: %s\n",
                      paste(head(names(body$body), 10), collapse = ", ")))
        }
      }
    }
    "OK"
  }, error = function(e) {
    cat(sprintf("    FAILED: %s\n", e$message))
    "FAIL"
  })
}


# =============================================================================
# PART 3: Summary
# =============================================================================

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# List saved files
saved_files <- list.files(cache_dir, pattern = "_external_predictors\\.rds$",
                          full.names = TRUE)
for (f in saved_files) {
  res <- readRDS(f)
  pred_cols <- setdiff(colnames(res), c("Admin1", "Admin2"))
  n_data <- sum(colSums(!is.na(res[, pred_cols, drop = FALSE])) > 0)
  cat(sprintf("  %-45s %4d Admin2 x %3d predictors (%3d with data)\n",
              basename(f), nrow(res), length(pred_cols), n_data))
}

# Source status
cat("\n  Source Status:\n")
sources <- c(
  "CHIRPS"      = "chirps_",
  "WorldPop"    = "worldpop_",
  "Nightlights" = "ntl_",
  "MAP"         = "map2_",
  "SoilGrids"   = "soil_",
  "GDL"         = "gdl_",
  "WFP"         = "wfp_",
  "IPC"         = "ipc_",
  "ACLED"       = "acled_",
  "HarvestStat" = "crop_"
)

# Use Ghana as reference
if ("Ghana" %in% names(results)) {
  ref <- results[["Ghana"]]
  for (src_name in names(sources)) {
    pfx <- sources[[src_name]]
    cols <- grep(paste0("^", pfx), colnames(ref), value = TRUE)
    n_data <- if (length(cols) > 0) sum(colSums(!is.na(ref[, cols, drop = FALSE])) > 0) else 0
    status <- if (n_data > 0) sprintf("OK (%d vars)", n_data) else "FAILED"
    cat(sprintf("    %-15s %s\n", src_name, status))
  }
}

cat("\nDone! External predictors cached in:", cache_dir, "\n")
