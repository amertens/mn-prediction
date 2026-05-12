# =============================================================================
# R/external_data.R
#
# Functions to download and extract predictor variables from external data
# sources for the micronutrient prediction pipeline.
#
# Each extract_* function is self-contained and independently callable.
# All functions cache downloads to avoid re-downloading.
# All functions return a data.frame with an Admin2 column + predictor columns.
# If a source fails, functions warn and return an empty data.frame.
#
# External sources:
#   1. CHIRPS rainfall              (chirps_*)
#   2. WorldPop population density  (worldpop_*)
#   3. VIIRS nighttime lights       (ntl_*)
#   4. Malaria Atlas Project        (map2_*)
#   5. SoilGrids / ISRIC           (soil_*)
#   6. Global Data Lab HDI          (gdl_*)
#   7. WFP HungerMap                (wfp_*)
#   8. IPC/Cadre Harmonisé API      (ipc_*)
#   9. ACLED conflict events        (acled_*)
#  10. HarvestStat Africa crops     (crop_*)
#
# Usage:
#   source("R/external_data.R")
#   cc <- get_country_configs()$Gambia
#   df <- extract_all_external(cc, survey_year = 2018,
#                              cache_dir = "data/external_cache")
# =============================================================================


# =============================================================================
# Helper utilities
# =============================================================================

#' Ensure a directory exists, creating it if necessary
#' @param path Directory path
#' @return path (invisibly)
.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Check whether a package is installed, without loading it
#' @param pkg Package name (character)
#' @return logical
.has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

#' Require a package or stop with a clear message
#' @param pkg Package name
#' @param purpose What the package is needed for
.require_pkg <- function(pkg, purpose = NULL) {
  if (!.has_pkg(pkg)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(purpose)) msg <- paste0(msg, " for ", purpose)
    msg <- paste0(msg, ". Install it with: install.packages('", pkg, "')")
    stop(msg, call. = FALSE)
  }
}

#' Build an empty result data.frame with an Admin2 column
#' @param admin2_names Character vector of Admin-2 names
#' @return data.frame with only the Admin2 column
.empty_result <- function(admin2_names = character(0)) {
  data.frame(Admin2 = admin2_names, stringsAsFactors = FALSE)
}

#' Safe zonal extraction: exact_extract with error handling
#'
#' @param raster A SpatRaster or RasterLayer
#' @param polys sf polygons
#' @param fun Aggregation function name (default "mean")
#' @param col_name Name for the resulting column
#' @return Numeric vector of extracted values (one per polygon)
.safe_extract <- function(raster, polys, fun = "mean", col_name = "value") {
  .require_pkg("exactextractr", "raster zonal extraction")
  tryCatch({
    vals <- exactextractr::exact_extract(raster, polys, fun = fun)
    vals
  }, error = function(e) {
    warning(sprintf("Zonal extraction failed for '%s': %s", col_name, e$message))
    rep(NA_real_, nrow(polys))
  })
}


# =============================================================================
# 1. CHIRPS Rainfall
# =============================================================================

#' Extract CHIRPS rainfall variables per Admin-2 polygon
#'
#' Downloads monthly CHIRPS rainfall data for the 12 months of the survey year,
#' extracts zonal means per Admin-2 polygon, and computes annual summary
#' statistics.
#'
#' @param admin2_sf sf object of Admin-2 polygons (must have "Admin2" column)
#' @param country_code ISO3 country code (e.g. "GMB")
#' @param survey_year Integer survey year
#' @param cache_dir Directory for caching downloaded rasters
#' @return data.frame with Admin2 + chirps_* columns
extract_chirps <- function(admin2_sf, country_code, survey_year, cache_dir) {

  admin2_names <- admin2_sf[["Admin2"]]
  if (is.null(admin2_names)) admin2_names <- admin2_sf[[1]]

  result <- tryCatch({

    .require_pkg("terra", "raster handling")
    .require_pkg("exactextractr", "zonal extraction")

    chirps_dir <- file.path(cache_dir, "chirps", country_code, survey_year)
    .ensure_dir(chirps_dir)

    month_names <- tolower(month.abb)  # jan, feb, ..., dec
    monthly_means <- list()

    for (m in 1:12) {
      month_str <- sprintf("%02d", m)
      cached_file <- file.path(chirps_dir, paste0("chirps_", survey_year, "_", month_str, ".tif"))

      if (!file.exists(cached_file)) {
        # CHIRPS monthly data via FTP/HTTP
        # Format: chirps-v2.0.YYYY.MM.tif.gz
        base_url <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs"
        fname <- sprintf("chirps-v2.0.%d.%s.tif.gz", survey_year, month_str)
        url <- file.path(base_url, fname)
        gz_file <- file.path(chirps_dir, fname)

        cat(sprintf("[CHIRPS] Downloading %s ...\n", fname))
        dl_status <- tryCatch({
          utils::download.file(url, gz_file, mode = "wb", quiet = TRUE)
          0L
        }, error = function(e) {
          warning(sprintf("[CHIRPS] Download failed for %s: %s", fname, e$message))
          1L
        })

        if (dl_status == 0L && file.exists(gz_file)) {
          # Decompress .gz to .tif
          R.utils::gunzip(gz_file, destname = cached_file, remove = TRUE, overwrite = TRUE)
        }
      }

      if (file.exists(cached_file)) {
        r <- terra::rast(cached_file)
        # Crop to country extent for speed
        r_crop <- tryCatch(
          terra::crop(r, terra::ext(terra::vect(admin2_sf)) + 0.5),
          error = function(e) r
        )
        monthly_means[[m]] <- .safe_extract(r_crop, admin2_sf, fun = "mean",
                                            col_name = paste0("chirps_precip_", month_names[m]))
      } else {
        monthly_means[[m]] <- rep(NA_real_, nrow(admin2_sf))
      }
    }

    # Build data.frame with monthly columns
    df <- data.frame(Admin2 = admin2_names, stringsAsFactors = FALSE)
    for (m in 1:12) {
      df[[paste0("chirps_precip_", month_names[m])]] <- monthly_means[[m]]
    }

    # Compute annual summaries
    monthly_mat <- do.call(cbind, monthly_means)
    df$chirps_annual_mean  <- rowMeans(monthly_mat, na.rm = TRUE)
    df$chirps_annual_sd    <- apply(monthly_mat, 1, sd, na.rm = TRUE)
    df$chirps_annual_total <- rowSums(monthly_mat, na.rm = TRUE)
    df$chirps_dry_months   <- rowSums(monthly_mat < 50, na.rm = TRUE)

    # Replace NaN with NA
    df[is.nan(df$chirps_annual_mean), "chirps_annual_mean"] <- NA_real_
    df[is.nan(df$chirps_annual_sd), "chirps_annual_sd"] <- NA_real_

    cat(sprintf("[CHIRPS] Extracted %d variables for %d Admin-2 units\n",
                ncol(df) - 1, nrow(df)))
    df

  }, error = function(e) {
    warning(sprintf("[CHIRPS] Failed: %s", e$message))
    .empty_result(admin2_names)
  })

  result
}


# =============================================================================
# 2. WorldPop Population Density
# =============================================================================

#' Extract WorldPop population density per Admin-2 polygon
#'
#' Downloads constrained population density raster from WorldPop for the
#' given country and survey year, then extracts zonal statistics.
#'
#' @param admin2_sf sf object of Admin-2 polygons
#' @param country_iso3 ISO3 country code (e.g. "GMB")
#' @param survey_year Integer survey year
#' @param cache_dir Directory for caching downloaded rasters
#' @return data.frame with Admin2 + worldpop_* columns
extract_worldpop <- function(admin2_sf, country_iso3, survey_year, cache_dir) {

  admin2_names <- admin2_sf[["Admin2"]]
  if (is.null(admin2_names)) admin2_names <- admin2_sf[[1]]

  result <- tryCatch({

    .require_pkg("terra", "raster handling")
    .require_pkg("exactextractr", "zonal extraction")

    wp_dir <- file.path(cache_dir, "worldpop", country_iso3)
    .ensure_dir(wp_dir)

    iso_lower <- tolower(country_iso3)
    cached_file <- file.path(wp_dir, paste0("worldpop_", country_iso3, "_", survey_year, ".tif"))

    if (!file.exists(cached_file)) {
      # WorldPop constrained population density URL pattern
      # Try the standard WorldPop FTP URL pattern
      url <- sprintf(
        "https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/%s/%s_ppp_%d_constrained.tif",
        toupper(country_iso3), iso_lower, survey_year
      )

      cat(sprintf("[WorldPop] Downloading population density for %s %d ...\n",
                  country_iso3, survey_year))

      dl_ok <- tryCatch({
        utils::download.file(url, cached_file, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)

      # Fallback: try unconstrained if constrained fails
      if (!dl_ok || !file.exists(cached_file)) {
        url2 <- sprintf(
          "https://data.worldpop.org/GIS/Population/Global_2000_2020/%d/%s/%s_ppp_%d.tif",
          survey_year, toupper(country_iso3), iso_lower, survey_year
        )
        cat("[WorldPop] Constrained download failed, trying unconstrained ...\n")
        tryCatch(
          utils::download.file(url2, cached_file, mode = "wb", quiet = TRUE),
          error = function(e) {
            warning(sprintf("[WorldPop] Download failed: %s", e$message))
          }
        )
      }
    }

    if (!file.exists(cached_file)) {
      warning("[WorldPop] No raster available — returning empty result")
      return(.empty_result(admin2_names))
    }

    r <- terra::rast(cached_file)
    r_crop <- tryCatch(
      terra::crop(r, terra::ext(terra::vect(admin2_sf)) + 0.5),
      error = function(e) r
    )

    df <- data.frame(
      Admin2              = admin2_names,
      worldpop_pop_density = .safe_extract(r_crop, admin2_sf, fun = "mean",
                                           col_name = "worldpop_pop_density"),
      worldpop_pop_total   = .safe_extract(r_crop, admin2_sf, fun = "sum",
                                           col_name = "worldpop_pop_total"),
      stringsAsFactors = FALSE
    )

    cat(sprintf("[WorldPop] Extracted 2 variables for %d Admin-2 units\n", nrow(df)))
    df

  }, error = function(e) {
    warning(sprintf("[WorldPop] Failed: %s", e$message))
    .empty_result(admin2_names)
  })

  result
}


# =============================================================================
# 3. VIIRS Nighttime Lights (via blackmarbler)
# =============================================================================

#' Extract VIIRS nighttime lights per Admin-2 polygon
#'
#' Uses the blackmarbler package to download annual VIIRS composites from
#' NASA Earthdata and extract zonal statistics.
#'
#' @param admin2_sf sf object of Admin-2 polygons
#' @param country_iso3 ISO3 country code
#' @param survey_year Integer survey year
#' @param bearer_token NASA Earthdata bearer token (character). If NULL, skips.
#' @param cache_dir Directory for caching downloaded rasters
#' @return data.frame with Admin2 + ntl_* columns
extract_nightlights <- function(admin2_sf, country_iso3, survey_year,
                                bearer_token = NULL, cache_dir) {

  admin2_names <- admin2_sf[["Admin2"]]
  if (is.null(admin2_names)) admin2_names <- admin2_sf[[1]]

  if (is.null(bearer_token) || bearer_token == "") {
    warning("[NTL] No NASA Earthdata bearer_token provided — skipping nightlights extraction")
    return(.empty_result(admin2_names))
  }

  result <- tryCatch({

    .require_pkg("blackmarbler", "VIIRS nighttime lights extraction")
    .require_pkg("terra", "raster handling")
    .require_pkg("exactextractr", "zonal extraction")

    ntl_dir <- file.path(cache_dir, "nightlights", country_iso3)
    .ensure_dir(ntl_dir)

    cached_file <- file.path(ntl_dir, paste0("ntl_", country_iso3, "_", survey_year, ".tif"))

    if (!file.exists(cached_file)) {
      cat(sprintf("[NTL] Downloading VIIRS annual composite for %s %d ...\n",
                  country_iso3, survey_year))

      # VIIRS nighttime lights — try multiple sources
      ntl_year <- max(as.integer(survey_year), 2012L)

      r_ntl <- NULL

      # Method 1: EOG direct download (no auth needed)
      # Colorado School of Mines hosts annual VIIRS composites as cloud-free GeoTIFFs
      cat("[NTL] Trying EOG direct download (no auth needed)...\n")
      r_ntl <- tryCatch({
        # EOG hosts annual composites at ~500m resolution
        # URL pattern for v2.2 annual composites (median masked)
        eog_url <- sprintf(
          "https://eogdata.mines.edu/nighttime_light/annual/v22/%d/VNP46A4/vcmslcfg/",
          ntl_year
        )

        # The files are large (~1GB global). Instead, use the average_masked tiles.
        # For Africa, tiles are typically h17v07, h17v08, h18v07, h18v08
        # But easier: use the vsicurl to read directly and crop
        # Try the global low-res version first
        eog_avg_url <- sprintf(
          "/vsicurl/https://eogdata.mines.edu/nighttime_light/annual/v22/%d/VNP46A4/vcmslcfg/VNP46A4_t%d.average_masked.tif",
          ntl_year, ntl_year
        )

        cat(sprintf("[NTL] Trying global composite for %d...\n", ntl_year))
        r_global <- tryCatch(terra::rast(eog_avg_url), error = function(e) NULL)

        if (!is.null(r_global)) {
          bbox_ext <- terra::ext(
            bbox["xmin"] - 0.5, bbox["xmax"] + 0.5,
            bbox["ymin"] - 0.5, bbox["ymax"] + 0.5
          )
          cropped <- terra::crop(r_global, bbox_ext)
          cat(sprintf("[NTL] Cropped EOG raster: %d x %d pixels\n",
                      terra::nrow(cropped), terra::ncol(cropped)))
          terra::writeRaster(cropped, cached_file, overwrite = TRUE)
          cropped
        } else {
          NULL
        }
      }, error = function(e) {
        cat(sprintf("[NTL] EOG download failed: %s\n", e$message))
        NULL
      })

      # Method 2: Use existing CCNL rasters from GEE directory
      if (is.null(r_ntl)) {
        cat("[NTL] Trying existing CCNL raster from GEE directory...\n")
        # Search for CCNL files in the GEE raster directories
        gee_dirs <- list.dirs(here::here("data"), recursive = FALSE)
        gee_dirs <- gee_dirs[grepl("GEE", gee_dirs, ignore.case = TRUE)]

        for (gd in gee_dirs) {
          ccnl_files <- list.files(gd, pattern = "CCNL.*\\.tif$", full.names = TRUE,
                                   ignore.case = TRUE)
          # Also check if the directory name matches our country
          country_match <- grepl(gsub("_", ".", country_iso3), gd, ignore.case = TRUE) ||
                           grepl(gsub("_", " ", gsub("GEE.*", "", basename(gd))),
                                 admin2_sf$Admin1[1], ignore.case = TRUE)
          if (length(ccnl_files) > 0 && country_match) {
            cat(sprintf("[NTL] Found CCNL raster: %s\n", basename(ccnl_files[1])))
            r_ntl <- tryCatch({
              r <- terra::rast(ccnl_files[1])
              if (terra::nlyr(r) > 1) r <- r[[1]]
              terra::writeRaster(r, cached_file, overwrite = TRUE)
              r
            }, error = function(e) {
              cat(sprintf("[NTL] Cannot read CCNL: %s\n", e$message))
              NULL
            })
            if (!is.null(r_ntl)) break
          }
        }
      }

      # Method 3: Try blackmarbler (needs NASA auth)
      if (is.null(r_ntl) && .has_pkg("blackmarbler") &&
          !is.null(bearer_token) && nchar(bearer_token) > 10) {
        cat("[NTL] Trying NASA Black Marble via blackmarbler...\n")
        r_ntl <- tryCatch({
          country_sf <- sf::st_union(admin2_sf)
          country_sf <- sf::st_sf(geometry = country_sf)

          blackmarbler::bm_raster(
            roi_sf     = country_sf,
            product_id = "VNP46A4",
            date       = ntl_year,
            bearer     = bearer_token,
            h5_dir     = ntl_dir,
            quiet      = FALSE
          )
        }, error = function(e) {
          cat(sprintf("[NTL] blackmarbler failed: %s\n", e$message))
          NULL
        })
      }

      # blackmarbler may return a SpatRaster; save for caching
      if (!is.null(r_ntl) && inherits(r_ntl, "SpatRaster")) {
        terra::writeRaster(r_ntl, cached_file, overwrite = TRUE)
        cat(sprintf("[NTL] Saved raster: %s\n", cached_file))
      } else {
        cat("[NTL] bm_raster returned NULL or non-raster object\n")
        cat(sprintf("[NTL] Object class: %s\n", paste(class(r_ntl), collapse = ", ")))
      }
    } else {
      cat(sprintf("[NTL] Using cached: %s\n", basename(cached_file)))
    }

    if (!file.exists(cached_file)) {
      cat("[NTL] No raster file available — returning empty result\n")
      return(.empty_result(admin2_names))
    }

    r <- terra::rast(cached_file)

    df <- data.frame(
      Admin2     = admin2_names,
      ntl_mean   = .safe_extract(r, admin2_sf, fun = "mean",   col_name = "ntl_mean"),
      ntl_median = .safe_extract(r, admin2_sf, fun = "median", col_name = "ntl_median"),
      ntl_max    = .safe_extract(r, admin2_sf, fun = "max",    col_name = "ntl_max"),
      ntl_sd     = .safe_extract(r, admin2_sf, fun = "stdev",  col_name = "ntl_sd"),
      stringsAsFactors = FALSE
    )

    cat(sprintf("[NTL] Extracted 4 variables for %d Admin-2 units\n", nrow(df)))
    df

  }, error = function(e) {
    warning(sprintf("[NTL] Failed: %s", e$message))
    .empty_result(admin2_names)
  })

  result
}


# =============================================================================
# 4. Malaria Atlas Project
# =============================================================================

#' Extract Malaria Atlas Project raster surfaces per Admin-2 polygon
#'
#' Downloads PfPR, ITN coverage, IRS coverage, ACT coverage, and travel
#' time to cities from the Malaria Atlas Project and extracts zonal means.
#'
#' @param admin2_sf sf object of Admin-2 polygons
#' @param survey_year Integer survey year
#' @param cache_dir Directory for caching downloaded rasters
#' @return data.frame with Admin2 + map2_* columns
extract_malaria_atlas <- function(admin2_sf, survey_year, cache_dir,
                                  country_iso3 = NULL) {

  admin2_names <- admin2_sf[["Admin2"]]
  if (is.null(admin2_names)) admin2_names <- admin2_sf[[1]]

  # Determine country code for per-country caching
  if (is.null(country_iso3)) {
    # Try to infer from the sf object's extent
    country_iso3 <- "UNKNOWN"
  }

  result <- tryCatch({

    .require_pkg("malariaAtlas", "Malaria Atlas Project data")
    .require_pkg("terra", "raster handling")
    .require_pkg("exactextractr", "zonal extraction")

    # Per-country cache directory — rasters are clipped to each country's extent
    map_dir <- file.path(cache_dir, "malaria_atlas", country_iso3)
    .ensure_dir(map_dir)

    # List ALL available rasters
    available <- tryCatch(
      malariaAtlas::listRaster(),
      error = function(e) {
        warning("[MAP] Could not list available rasters: ", e$message)
        NULL
      }
    )

    if (is.null(available) || nrow(available) == 0) {
      cat("[MAP] No rasters available from malariaAtlas\n")
      return(.empty_result(admin2_names))
    }

    cat(sprintf("[MAP] Found %d available rasters — downloading ALL\n", nrow(available)))

    # Skip Anopheles mosquito species maps (50+, very niche for malaria entomology)
    # and Macaca monkey distribution maps (irrelevant)
    # Keep Explorer__ datasets — they contain useful accessibility/housing/blood disorder data
    skip_patterns <- c("Anopheles_", "Macaca_", "Pk_SEAsia",
                        "Dominant_Vector", "Secondary.Dominant_Vector")
    keep_mask <- !Reduce("|", lapply(skip_patterns, function(p) {
      grepl(p, available$dataset_id)
    }))
    available <- available[keep_mask, ]
    cat(sprintf("[MAP] After filtering (no mosquito species/Explorer dupes): %d rasters\n",
                nrow(available)))

    # For datasets with multiple years, pick the one closest to (but not after) survey_year
    # Group by base name (strip year from dataset_id)
    available$base_name <- gsub("__\\d{6}_", "__YYYY_", available$dataset_id)

    # For each base_name group, find best year match
    best_rows <- list()
    for (bn in unique(available$base_name)) {
      group <- available[available$base_name == bn, ]

      # Extract year from dataset_id.
      # Malaria__/Interventions__/Blood__/Accessibility__ use 6-digit dates (e.g., "202206")
      # Explorer__ uses 4-digit years (e.g., "2020")
      group$ds_year <- sapply(group$dataset_id, function(did) {
        m6 <- regmatches(did, regexpr("\\d{6}", did))
        if (length(m6) > 0 && nchar(m6) == 6) return(as.integer(substr(m6, 1, 4)))
        # Try 4-digit year after "__"
        m4 <- regmatches(did, regexpr("(?<=__)\\d{4}", did, perl = TRUE))
        if (length(m4) > 0) return(as.integer(m4))
        9999L
      })
      group$ds_year[is.na(group$ds_year)] <- 9999L

      # Prefer years at or before survey_year; otherwise nearest
      prior <- group[group$ds_year <= survey_year, ]
      if (nrow(prior) > 0) {
        best <- prior[which.max(prior$ds_year), ]
      } else {
        best <- group[which.min(abs(group$ds_year - survey_year)), ]
      }
      best_rows[[bn]] <- best[1, ]
    }
    to_download <- do.call(rbind, best_rows)
    cat(sprintf("[MAP] Deduplicated to %d unique surfaces (best year match for %d)\n",
                nrow(to_download), survey_year))

    df <- data.frame(Admin2 = admin2_names, stringsAsFactors = FALSE)
    n_success <- 0L
    n_fail <- 0L

    for (i in seq_len(nrow(to_download))) {
      ds_id <- to_download$dataset_id[i]

      # Create a clean variable name from dataset_id
      col_name <- paste0("map2_", tolower(gsub("[^A-Za-z0-9]+", "_", ds_id)))
      col_name <- gsub("_+", "_", col_name)
      col_name <- sub("_$", "", col_name)
      # Truncate very long names
      if (nchar(col_name) > 60) col_name <- substr(col_name, 1, 60)

      cached_file <- file.path(map_dir, paste0(gsub("[^A-Za-z0-9_]", "_", ds_id), ".tif"))

      extracted <- tryCatch({
        if (!file.exists(cached_file)) {
          cat(sprintf("  [MAP %d/%d] Downloading %s ...\n", i, nrow(to_download), ds_id))

          # Use dataset_id directly.
          # Some datasets are temporal (have year-specific versions), others are static.
          # For temporal datasets (Malaria__, Interventions__), pass the embedded year.
          # For static datasets (Accessibility__, Blood_Disorders__, Explorer__), don't pass year.
          ds_year_str <- regmatches(ds_id, regexpr("\\d{6}", ds_id))
          has_temporal_year <- length(ds_year_str) > 0 && nchar(ds_year_str) == 6 &&
            any(grepl("^(Malaria|Interventions)", ds_id))

          if (has_temporal_year) {
            ds_yr <- as.integer(substr(ds_year_str, 1, 4))
            r <- malariaAtlas::getRaster(
              dataset_id = ds_id,
              shp = admin2_sf,
              year = ds_yr
            )
          } else {
            r <- malariaAtlas::getRaster(
              dataset_id = ds_id,
              shp = admin2_sf
            )
          }

          # getRaster can return: SpatRaster, Raster*, SpatRasterCollection,
          # or a file path string (when year doesn't match available years).
          if (is.character(r)) {
            # File path returned — try to read it as a raster
            if (file.exists(r)) {
              r <- terra::rast(r)
            } else {
              warning(sprintf("[MAP] getRaster returned path that doesn't exist: %s", r))
              next
            }
          }
          if (inherits(r, "Raster")) r <- terra::rast(r)
          if (inherits(r, "SpatRasterCollection")) {
            r <- r[[1]]
          }
          if (!inherits(r, "SpatRaster")) {
            warning(sprintf("[MAP] getRaster returned unexpected type: %s", class(r)[1]))
            next
          }
          terra::writeRaster(r, cached_file, overwrite = TRUE)
        } else {
          cat(sprintf("  [MAP %d/%d] Cached: %s\n", i, nrow(to_download), ds_id))
        }

        if (file.exists(cached_file)) {
          r <- terra::rast(cached_file)
          if (terra::nlyr(r) > 1) r <- r[[1]]
          vals <- .safe_extract(r, admin2_sf, fun = "mean", col_name = col_name)
          n_success <<- n_success + 1L
          vals
        } else {
          rep(NA_real_, nrow(admin2_sf))
        }

      }, error = function(e) {
        cat(sprintf("  [MAP %d/%d] FAILED %s: %s\n", i, nrow(to_download), ds_id, e$message))
        n_fail <<- n_fail + 1L
        rep(NA_real_, nrow(admin2_sf))
      })

      df[[col_name]] <- extracted
    }

    n_valid <- sum(vapply(df[-1], function(x) !all(is.na(x)), logical(1)))
    cat(sprintf("[MAP] Done: %d/%d downloaded, %d failed, %d with data for %d Admin-2 units\n",
                n_success, nrow(to_download), n_fail, n_valid, nrow(df)))
    df

  }, error = function(e) {
    warning(sprintf("[MAP] Failed: %s", e$message))
    .empty_result(admin2_names)
  })

  result
}


# =============================================================================
# 5. SoilGrids (ISRIC)
# =============================================================================

#' Extract SoilGrids soil properties per Admin-2 polygon
#'
#' Uses the ISRIC SoilGrids REST API (or pre-downloaded tiles) to extract
#' soil properties at 0-5cm depth: pH, organic carbon, nitrogen, clay, sand,
#' silt, and CEC.
#'
#' @param admin2_sf sf object of Admin-2 polygons
#' @param cache_dir Directory for caching downloaded rasters
#' @return data.frame with Admin2 + soil_* columns
extract_soilgrids <- function(admin2_sf, cache_dir) {

  admin2_names <- admin2_sf[["Admin2"]]
  if (is.null(admin2_names)) admin2_names <- admin2_sf[[1]]

  result <- tryCatch({

    soil_dir <- file.path(cache_dir, "soilgrids")
    .ensure_dir(soil_dir)

    cached_rds <- file.path(soil_dir,
                            paste0("soilgrids_", digest::digest(sf::st_bbox(admin2_sf)), ".rds"))

    if (file.exists(cached_rds)) {
      cat("[SoilGrids] Loading from cache\n")
      return(readRDS(cached_rds))
    }

    # Download full SoilGrids rasters from the ISRIC file server (VRT/COG).
    # These are Cloud-Optimized GeoTIFFs that terra can read with windowed access.
    # Much faster than per-point API queries (7 downloads vs 260×7).
    # File server: https://files.isric.org/soilgrids/latest/data/
    .require_pkg("terra", "raster handling")
    .require_pkg("exactextractr", "zonal extraction")

    properties <- list(
      list(name = "phh2o",    col = "soil_ph"),
      list(name = "soc",      col = "soil_organic_carbon"),
      list(name = "nitrogen", col = "soil_nitrogen"),
      list(name = "clay",     col = "soil_clay"),
      list(name = "sand",     col = "soil_sand"),
      list(name = "silt",     col = "soil_silt"),
      list(name = "cec",      col = "soil_cec")
    )

    bbox <- sf::st_bbox(admin2_sf)
    df <- data.frame(Admin2 = admin2_names, stringsAsFactors = FALSE)

    for (prop in properties) {
      cached_tif <- file.path(soil_dir, paste0(prop$col, "_0_5cm.tif"))

      if (!file.exists(cached_tif)) {
        cat(sprintf("[SoilGrids] Downloading %s (0-5cm mean) ...\n", prop$name))

        # The ISRIC file server hosts VRT files that point to COG tiles.
        # We can read them directly with terra using GDAL's /vsicurl/ driver,
        # which does windowed reads (only downloads the region we need).
        vrt_url <- sprintf(
          "/vsicurl/https://files.isric.org/soilgrids/latest/data/%s/%s_0-5cm_mean.vrt",
          prop$name, prop$name
        )

        r <- tryCatch({
          # Read the VRT — GDAL will only fetch tiles covering our bbox
          full_r <- terra::rast(vrt_url)

          # SoilGrids uses Homolosine projection (EPSG:152160 / igh).
          # We must transform our bbox to the raster's CRS before cropping.
          raster_crs <- terra::crs(full_r)
          admin2_reproj <- sf::st_transform(admin2_sf, raster_crs)
          bbox_reproj <- sf::st_bbox(admin2_reproj)

          ext <- terra::ext(
            bbox_reproj["xmin"] - 50000, bbox_reproj["xmax"] + 50000,
            bbox_reproj["ymin"] - 50000, bbox_reproj["ymax"] + 50000
          )
          cropped <- terra::crop(full_r, ext)

          cat(sprintf("  [SoilGrids] Cropped: %d x %d pixels\n",
                      terra::nrow(cropped), terra::ncol(cropped)))

          # Reproject to WGS84 for consistent extraction
          cropped_wgs84 <- terra::project(cropped, "EPSG:4326")

          # Save locally for caching
          terra::writeRaster(cropped_wgs84, cached_tif, overwrite = TRUE)
          cat(sprintf("  [SoilGrids] Cached: %s (%d x %d)\n",
                      basename(cached_tif), terra::nrow(cropped_wgs84), terra::ncol(cropped_wgs84)))
          cropped_wgs84
        }, error = function(e) {
          cat(sprintf("  [SoilGrids] VRT/vsicurl failed for %s: %s\n", prop$name, e$message))
          NULL
        })
      } else {
        cat(sprintf("[SoilGrids] Cached: %s\n", prop$name))
        r <- tryCatch(terra::rast(cached_tif), error = function(e) NULL)
      }

      if (!is.null(r)) {
        df[[prop$col]] <- .safe_extract(r, admin2_sf, fun = "mean", col_name = prop$col)
      } else {
        df[[prop$col]] <- NA_real_
      }
    }

    n_valid <- sum(vapply(df[-1], function(x) !all(is.na(x)), logical(1)))
    cat(sprintf("[SoilGrids] Extracted %d/%d properties for %d Admin-2 units\n",
                n_valid, length(properties), nrow(df)))
    if (n_valid > 0) saveRDS(df, cached_rds)
    df

  }, error = function(e) {
    warning(sprintf("[SoilGrids] Failed: %s", e$message))
    .empty_result(admin2_names)
  })

  result
}


# =============================================================================
# 6. Global Data Lab (Subnational HDI)
# =============================================================================

#' Extract Global Data Lab subnational development indices
#'
#' Downloads the SHDI (Subnational Human Development Index) dataset from
#' Global Data Lab and matches to Admin-1 names using fuzzy matching.
#' Returns Admin-1 level data (replicated to Admin-2 via join).
#'
#' @param country_iso3 ISO3 country code (e.g. "GMB")
#' @param admin1_names Character vector of Admin-1 names to match
#' @param survey_year Integer survey year
#' @return data.frame with Admin1 + gdl_* columns
extract_gdl <- function(country_iso3, admin1_names, survey_year) {

  result <- tryCatch({

    .require_pkg("jsonlite", "GDL data parsing")

    # GDL SHDI dataset — look for local file first, then try web download
    cat(sprintf("[GDL] Loading subnational HDI data for %s ...\n", country_iso3))

    # Search for local file with common naming patterns
    cache_dir <- here::here("data", "external_cache")
    local_candidates <- c(
      file.path(cache_dir, "GDL-Subnational-HDI-data.csv"),
      file.path(cache_dir, "gdl_shdi_full.csv"),
      file.path(cache_dir, "SHDI-SGDI-Total 7.0.csv"),
      file.path(cache_dir, "SHDI-SGDI-Total 8.0.csv")
    )
    # Also search for any CSV with "GDL" or "SHDI" or "HDI" in the name
    if (dir.exists(cache_dir)) {
      all_csvs <- list.files(cache_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
      gdl_csvs <- all_csvs[grepl("gdl|shdi|hdi", basename(all_csvs), ignore.case = TRUE)]
      local_candidates <- unique(c(local_candidates, gdl_csvs))
    }

    local_file <- NULL
    for (f in local_candidates) {
      if (file.exists(f)) {
        local_file <- f
        cat(sprintf("[GDL] Found local file: %s\n", basename(f)))
        break
      }
    }

    gdl_data <- if (!is.null(local_file)) {
      tryCatch(
        utils::read.csv(local_file, stringsAsFactors = FALSE),
        error = function(e) {
          cat(sprintf("[GDL] Failed to read local file: %s\n", e$message))
          NULL
        }
      )
    } else {
      # Try web download as fallback
      tryCatch({
        gdl_url <- "https://globaldatalab.org/assets/2024/09/SHDI-SGDI-Total%207.0.csv"
        tmp <- tempfile(fileext = ".csv")
        utils::download.file(gdl_url, tmp, mode = "wb", quiet = TRUE)
        utils::read.csv(tmp, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(sprintf("[GDL] Web download failed: %s\n", e$message))
        NULL
      })
    }

    if (is.null(gdl_data) || nrow(gdl_data) == 0) {
      warning("[GDL] No data downloaded — returning empty result")
      return(data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE))
    }

    # Filter to country and closest year
    # Column names vary across GDL versions; try common patterns
    iso_col <- intersect(c("iso_code", "ISO_Code", "iso3", "country_code"),
                         colnames(gdl_data))
    year_col <- intersect(c("year", "Year"), colnames(gdl_data))
    region_col <- intersect(c("region", "Region", "sub_nat_name", "GDLCODE",
                              "subnational"), colnames(gdl_data))

    if (length(iso_col) == 0 || length(year_col) == 0 || length(region_col) == 0) {
      warning("[GDL] Cannot identify expected columns in GDL data. ",
              "Available columns: ", paste(colnames(gdl_data), collapse = ", "))
      return(data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE))
    }

    iso_col <- iso_col[1]
    year_col <- year_col[1]
    region_col <- region_col[1]

    # Filter country
    country_data <- gdl_data[gdl_data[[iso_col]] == country_iso3, ]
    if (nrow(country_data) == 0) {
      # Try lowercase
      country_data <- gdl_data[toupper(gdl_data[[iso_col]]) == toupper(country_iso3), ]
    }

    if (nrow(country_data) == 0) {
      warning(sprintf("[GDL] No data found for country '%s'", country_iso3))
      return(data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE))
    }

    # Find closest year
    available_years <- sort(unique(country_data[[year_col]]))
    closest_year <- available_years[which.min(abs(available_years - survey_year))]
    cat(sprintf("[GDL] Using year %d (requested %d, available: %s)\n",
                closest_year, survey_year, paste(available_years, collapse = ", ")))

    year_data <- country_data[country_data[[year_col]] == closest_year, ]

    # Identify index columns
    shdi_col <- intersect(c("shdi", "SHDI", "Sub-national_HDI"), colnames(year_data))
    educ_col <- intersect(c("edindex", "Education_index", "educ"), colnames(year_data))
    health_col <- intersect(c("healthindex", "Health_index", "health"), colnames(year_data))
    income_col <- intersect(c("incindex", "Income_index", "income"), colnames(year_data))

    # Fuzzy match region names to Admin-1 names
    gdl_regions <- year_data[[region_col]]

    match_idx <- vapply(admin1_names, function(a1) {
      # Try exact match first
      exact <- match(tolower(a1), tolower(gdl_regions))
      if (!is.na(exact)) return(exact)

      # Fuzzy match via agrep
      fuzzy <- agrep(a1, gdl_regions, max.distance = 0.3, ignore.case = TRUE)
      if (length(fuzzy) == 1) return(fuzzy)
      if (length(fuzzy) > 1) {
        # Pick the closest (shortest edit distance)
        dists <- utils::adist(tolower(a1), tolower(gdl_regions[fuzzy]))[1, ]
        return(fuzzy[which.min(dists)])
      }

      # Try partial matching
      partial <- grep(a1, gdl_regions, ignore.case = TRUE, fixed = FALSE)
      if (length(partial) >= 1) return(partial[1])

      NA_integer_
    }, integer(1))

    df <- data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE)

    .safe_col <- function(data, cols, idx) {
      if (length(cols) == 0) return(rep(NA_real_, length(idx)))
      vals <- data[[cols[1]]]
      out <- rep(NA_real_, length(idx))
      valid <- !is.na(idx)
      out[valid] <- as.numeric(vals[idx[valid]])
      out
    }

    df$gdl_shdi          <- .safe_col(year_data, shdi_col, match_idx)
    df$gdl_education_idx <- .safe_col(year_data, educ_col, match_idx)
    df$gdl_health_idx    <- .safe_col(year_data, health_col, match_idx)
    df$gdl_income_idx    <- .safe_col(year_data, income_col, match_idx)

    n_matched <- sum(!is.na(match_idx))
    cat(sprintf("[GDL] Matched %d/%d Admin-1 regions\n", n_matched, length(admin1_names)))
    df

  }, error = function(e) {
    warning(sprintf("[GDL] Failed: %s", e$message))
    data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE)
  })

  result
}


# =============================================================================
# 7. WFP HungerMap
# =============================================================================

#' Extract WFP food security indicators
#'
#' Queries WFP HungerMap LIVE API or Data Bridges API for food security
#' indicators: FCS (food consumption score), rCSI, and prevalence of
#' insufficient food consumption.
#'
#' @param country_iso3 ISO3 country code (e.g. "GMB")
#' @param admin1_names Character vector of Admin-1 names
#' @param api_key WFP Data Bridges API key (optional). If NULL, tries public endpoints.
#' @return data.frame with Admin1 + wfp_* columns
extract_wfp_hungermap <- function(country_iso3, admin1_names, api_key = NULL) {

  result <- tryCatch({

    .require_pkg("jsonlite", "JSON parsing")

    df <- data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE)

    # Map ISO3 to ISO numeric code (used by some WFP endpoints)
    # Common codes for our countries
    iso_map <- c(GMB = 270, GHA = 288, SLE = 694, MWI = 454)
    country_num <- iso_map[country_iso3]

    got_data <- FALSE

    # ── Try Data Bridges API if key provided ──────────────────────────────
    if (!is.null(api_key) && nzchar(api_key)) {
      cat("[WFP] Querying Data Bridges API ...\n")

      # Food security endpoint
      # https://api.wfp.org/api/v1/foodsecurity/country/{iso3}
      api_base <- "https://api.wfp.org/api/v1"
      headers <- c(
        "Ocp-Apim-Subscription-Key" = api_key,
        "Accept" = "application/json"
      )

      fs_url <- sprintf("%s/foodsecurity/country/%s", api_base, country_iso3)

      resp <- tryCatch({
        tmp <- tempfile()
        utils::download.file(fs_url, tmp, method = "libcurl", quiet = TRUE,
                             headers = headers)
        jsonlite::fromJSON(readLines(tmp, warn = FALSE))
      }, error = function(e) {
        warning(sprintf("[WFP] Data Bridges API query failed: %s", e$message))
        NULL
      })

      if (!is.null(resp) && length(resp) > 0) {
        # Parse response — structure varies by endpoint
        # Try to extract subnational data
        if ("body" %in% names(resp)) resp <- resp$body
        if (is.data.frame(resp) && nrow(resp) > 0) {
          # Match admin names
          region_col <- intersect(c("admin1Name", "adm1_name", "region"),
                                  colnames(resp))
          if (length(region_col) > 0) {
            match_idx <- vapply(admin1_names, function(a1) {
              exact <- match(tolower(a1), tolower(resp[[region_col[1]]]))
              if (!is.na(exact)) return(exact)
              fuzzy <- agrep(a1, resp[[region_col[1]]],
                             max.distance = 0.3, ignore.case = TRUE)
              if (length(fuzzy) >= 1) return(fuzzy[1])
              NA_integer_
            }, integer(1))

            fcs_col <- intersect(c("fcs", "FCS", "fcsPrevalence"), colnames(resp))
            rcsi_col <- intersect(c("rcsi", "rCSI", "rcsiPrevalence"), colnames(resp))
            insuf_col <- intersect(c("prevalence", "insufficientFoodConsumption",
                                     "fc_prevalence"), colnames(resp))

            .get_vals <- function(cols) {
              if (length(cols) == 0) return(rep(NA_real_, length(match_idx)))
              out <- rep(NA_real_, length(match_idx))
              valid <- !is.na(match_idx)
              out[valid] <- as.numeric(resp[[cols[1]]][match_idx[valid]])
              out
            }

            df$wfp_fcs   <- .get_vals(fcs_col)
            df$wfp_rcsi  <- .get_vals(rcsi_col)
            df$wfp_insufficient_food <- .get_vals(insuf_col)
            got_data <- TRUE
          }
        }
      }
    }

    # ── Try public HungerMap LIVE endpoints ───────────────────────────────
    if (!got_data) {
      cat("[WFP] Trying public HungerMap LIVE endpoints ...\n")

      # HungerMap LIVE public API
      # https://static.hungermapdata.org/api-catalog/
      hm_url <- sprintf(
        "https://api.hungermapdata.org/v2/adm0/%s/adm1.json",
        if (!is.null(country_num)) country_num else country_iso3
      )

      # Alternative endpoint
      hm_url2 <- sprintf(
        "https://static.hungermapdata.org/hunger-map-data/adm0/%s/adm1.json",
        if (!is.null(country_num)) country_num else country_iso3
      )

      hm_data <- tryCatch({
        tmp <- tempfile(fileext = ".json")
        tryCatch(
          utils::download.file(hm_url, tmp, quiet = TRUE),
          error = function(e) utils::download.file(hm_url2, tmp, quiet = TRUE)
        )
        jsonlite::fromJSON(readLines(tmp, warn = FALSE))
      }, error = function(e) {
        warning(sprintf("[WFP] HungerMap public API failed: %s", e$message))
        NULL
      })

      if (!is.null(hm_data)) {
        # Parse HungerMap response
        regions <- NULL
        if (is.list(hm_data) && "body" %in% names(hm_data)) {
          regions <- hm_data$body
        } else if (is.data.frame(hm_data)) {
          regions <- hm_data
        } else if (is.list(hm_data) && length(hm_data) > 0) {
          # Try to coerce list of regions to data.frame
          regions <- tryCatch(
            do.call(rbind, lapply(hm_data, as.data.frame, stringsAsFactors = FALSE)),
            error = function(e) NULL
          )
        }

        if (!is.null(regions) && is.data.frame(regions) && nrow(regions) > 0) {
          region_col <- intersect(c("name", "admin1", "adm1_name", "region"),
                                  colnames(regions))
          if (length(region_col) > 0) {
            match_idx <- vapply(admin1_names, function(a1) {
              exact <- match(tolower(a1), tolower(regions[[region_col[1]]]))
              if (!is.na(exact)) return(exact)
              fuzzy <- agrep(a1, regions[[region_col[1]]],
                             max.distance = 0.3, ignore.case = TRUE)
              if (length(fuzzy) >= 1) return(fuzzy[1])
              NA_integer_
            }, integer(1))

            fcs_col <- intersect(c("fcs", "FCS", "metrics.fcs.people",
                                   "food_consumption_score"), colnames(regions))
            rcsi_col <- intersect(c("rcsi", "rCSI", "metrics.rcsi.people",
                                    "reduced_coping_strategy_index"), colnames(regions))

            .get_vals2 <- function(cols) {
              if (length(cols) == 0) return(rep(NA_real_, length(match_idx)))
              out <- rep(NA_real_, length(match_idx))
              valid <- !is.na(match_idx)
              out[valid] <- as.numeric(regions[[cols[1]]][match_idx[valid]])
              out
            }

            df$wfp_fcs   <- .get_vals2(fcs_col)
            df$wfp_rcsi  <- .get_vals2(rcsi_col)
            got_data <- TRUE
          }
        }
      }
    }

    if (!got_data) {
      warning("[WFP] No food security data could be retrieved for ", country_iso3)
      df$wfp_fcs  <- NA_real_
      df$wfp_rcsi <- NA_real_
      df$wfp_insufficient_food <- NA_real_
    }

    # Ensure all expected columns exist
    for (col in c("wfp_fcs", "wfp_rcsi", "wfp_insufficient_food")) {
      if (!col %in% colnames(df)) df[[col]] <- NA_real_
    }

    n_valid <- sum(!is.na(df$wfp_fcs))
    cat(sprintf("[WFP] Got data for %d/%d Admin-1 regions\n",
                n_valid, length(admin1_names)))
    df

  }, error = function(e) {
    warning(sprintf("[WFP] Failed: %s", e$message))
    df <- data.frame(
      Admin1               = admin1_names,
      wfp_fcs              = NA_real_,
      wfp_rcsi             = NA_real_,
      wfp_insufficient_food = NA_real_,
      stringsAsFactors     = FALSE
    )
    df
  })

  result
}


# =============================================================================
# 8. IPC/CH API — Food security phase classifications via ripc
# =============================================================================

#' Extract IPC/Cadre Harmonisé food security data via the ripc package
#'
#' Uses the IPC API to download area-level food security phase classifications
#' at Admin-2. Returns population counts per phase, phase proportions, and
#' combined phase 3-5 (crisis+emergency+famine) proportion.
#'
#' Requires: ripc package (CRAN) + IPC_API_KEY environment variable.
#' Register at https://docs.api.ipcinfo.org to get an API key.
#'
#' @param country_iso3 ISO3 country code
#' @param admin2_sf sf object with Admin2 polygons (for name matching)
#' @param survey_year Integer survey year (will find nearest analysis)
#' @return data.frame with Admin2 + ipc_* columns
extract_ipc <- function(country_iso3, admin2_sf, survey_year) {

  if (!.has_pkg("ripc")) {
    cat("  [ipc] ripc package not installed. Install with: install.packages('ripc')\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Check for API key
  api_key <- Sys.getenv("IPC_API_KEY", unset = "")
  if (nchar(api_key) == 0) {
    cat("  [ipc] No IPC_API_KEY set. Register at https://docs.api.ipcinfo.org\n")
    cat("  [ipc] Then: Sys.setenv(IPC_API_KEY = 'your_key')\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Map ISO3 to IPC country code (usually same)
  iso3_map <- c("GMB" = "GM", "GHA" = "GH", "SLE" = "SL", "MWI" = "MW")
  ipc_country <- country_iso3  # IPC API usually accepts ISO3

  cat(sprintf("  [ipc] Fetching IPC/CH data for %s near year %d...\n",
              country_iso3, survey_year))

  # Get area-level population data
  ipc_data <- tryCatch({
    pop_data <- ripc::ipc_get_population(
      country = ipc_country,
      year = survey_year,
      type = "A"  # Area level
    )

    if (is.null(pop_data) || length(pop_data) == 0) {
      # Try broader year range
      years_to_try <- c(survey_year, survey_year - 1, survey_year + 1,
                        survey_year - 2, survey_year + 2)
      for (yr in years_to_try) {
        pop_data <- tryCatch(
          ripc::ipc_get_population(country = ipc_country, year = yr, type = "A"),
          error = function(e) NULL
        )
        if (!is.null(pop_data) && length(pop_data) > 0) {
          cat(sprintf("  [ipc] Using year %d data (nearest to survey year %d)\n", yr, survey_year))
          break
        }
      }
    }

    # The ripc package returns a list with $areas element
    if (is.list(pop_data) && "areas" %in% names(pop_data)) {
      pop_data$areas
    } else if (is.data.frame(pop_data)) {
      pop_data
    } else {
      NULL
    }
  }, error = function(e) {
    cat(sprintf("  [ipc] API call failed: %s\n", e$message))
    NULL
  })

  if (is.null(ipc_data) || nrow(ipc_data) == 0) {
    cat("  [ipc] No data returned from IPC API\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  cat(sprintf("  [ipc] Got %d area records\n", nrow(ipc_data)))

  # Extract phase population columns
  phase_cols <- c("phase1", "phase2", "phase3", "phase4", "phase5")
  avail_phases <- intersect(phase_cols, colnames(ipc_data))

  # Compute proportions
  result <- data.frame(
    ipc_area_name = ipc_data$area_name %||% ipc_data$name %||% NA_character_,
    stringsAsFactors = FALSE
  )

  for (pc in avail_phases) {
    result[[paste0("ipc_", pc)]] <- as.numeric(ipc_data[[pc]])
  }

  # Total classified
  result$ipc_total <- rowSums(result[, paste0("ipc_", avail_phases), drop = FALSE], na.rm = TRUE)

  # Proportions
  for (pc in avail_phases) {
    result[[paste0("ipc_prop_", pc)]] <- ifelse(
      result$ipc_total > 0,
      result[[paste0("ipc_", pc)]] / result$ipc_total,
      NA_real_
    )
  }

  # Phase 3-5 combined
  p345 <- intersect(c("phase3", "phase4", "phase5"), avail_phases)
  if (length(p345) > 0) {
    result$ipc_prop_phase35 <- rowSums(
      result[, paste0("ipc_prop_", p345), drop = FALSE], na.rm = TRUE
    )
  }

  # Match IPC area names to Admin2 names
  target_names <- admin2_sf[["Admin2"]]
  source_names <- result$ipc_area_name[!is.na(result$ipc_area_name)]

  if (length(source_names) > 0) {
    # Build name lookup using the same fuzzy matching as food_security.R
    lookup <- stats::setNames(rep(NA_character_, length(source_names)), source_names)
    target_clean <- tolower(trimws(target_names))
    source_clean <- tolower(trimws(source_names))

    # Exact match first
    for (i in seq_along(source_names)) {
      exact <- which(target_clean == source_clean[i])
      if (length(exact) == 1) lookup[source_names[i]] <- target_names[exact]
    }

    # Fuzzy match
    unmatched <- is.na(lookup)
    if (any(unmatched)) {
      um <- source_names[unmatched]
      um_clean <- source_clean[unmatched]
      dist_mat <- utils::adist(um_clean, target_clean, ignore.case = TRUE, partial = FALSE)
      len_mat <- outer(nchar(um_clean), nchar(target_clean), pmax)
      norm_dist <- dist_mat / pmax(len_mat, 1)
      for (i in seq_along(um)) {
        best <- which.min(norm_dist[i, ])
        if (length(best) == 1 && norm_dist[i, best] < 0.3) {
          lookup[um[i]] <- target_names[best]
        }
      }
    }

    result$Admin2 <- lookup[result$ipc_area_name]

    # Keep only matched + predictor columns
    pred_cols <- grep("^ipc_prop_|^ipc_total", colnames(result), value = TRUE)
    result <- result[!is.na(result$Admin2), c("Admin2", pred_cols), drop = FALSE]

    # Aggregate if multiple periods per area
    if (nrow(result) > 0) {
      result <- result %>%
        dplyr::group_by(Admin2) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)),
                         .groups = "drop") %>%
        as.data.frame()
    }

    n_matched <- sum(target_names %in% result$Admin2)
    cat(sprintf("  [ipc] Matched %d/%d Admin-2 areas\n", n_matched, length(target_names)))
  } else {
    result <- .empty_result(target_names)
  }

  result
}


# =============================================================================
# 9. ACLED — Armed Conflict Location & Event Data
# =============================================================================

#' Extract ACLED conflict indicators aggregated to Admin-2
#'
#' Downloads geolocated conflict events for a country and time window,
#' then aggregates to Admin-2 by spatial join or the ACLED admin2 field.
#'
#' Requires: acled.api package (CRAN) + ACLED credentials.
#' Register (free academic) at https://developer.acleddata.com/
#' Set: Sys.setenv(ACLED_EMAIL_ADDRESS = "...", ACLED_ACCESS_KEY = "...")
#'
#' @param admin2_sf sf object with Admin2 polygons
#' @param country_name Country name as used by ACLED (e.g., "Ghana", "Gambia")
#' @param survey_year Integer survey year
#' @param window_years How many years before survey_year to include (default 3)
#' @return data.frame with Admin2 + acled_* columns
extract_acled <- function(admin2_sf, country_name, survey_year,
                          window_years = 3L) {

  if (!.has_pkg("acled.api")) {
    cat("  [acled] acled.api package not installed. Install with: install.packages('acled.api')\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Check credentials
  acled_email <- Sys.getenv("ACLED_EMAIL_ADDRESS", unset = "")
  acled_key   <- Sys.getenv("ACLED_ACCESS_KEY", unset = "")
  if (nchar(acled_email) == 0 || nchar(acled_key) == 0) {
    cat("  [acled] No ACLED credentials set. Register at https://developer.acleddata.com/\n")
    cat("  [acled] Then: Sys.setenv(ACLED_EMAIL_ADDRESS='...', ACLED_ACCESS_KEY='...')\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # ACLED country names
  acled_country_map <- c(
    "Gambia" = "Gambia", "Ghana" = "Ghana",
    "Sierra Leone" = "Sierra Leone", "Malawi" = "Malawi"
  )
  acled_name <- acled_country_map[country_name]
  if (is.na(acled_name)) acled_name <- country_name

  start_date <- sprintf("%d-01-01", survey_year - window_years)
  end_date   <- sprintf("%d-12-31", survey_year)

  cat(sprintf("  [acled] Fetching conflict events for %s (%s to %s)...\n",
              acled_name, start_date, end_date))

  events <- tryCatch({
    acled.api::acled.api(
      email.address = acled_email,
      access.key    = acled_key,
      country       = acled_name,
      start.date    = start_date,
      end.date      = end_date,
      all.variables = TRUE
    )
  }, error = function(e) {
    cat(sprintf("  [acled] API call failed: %s\n", e$message))
    NULL
  })

  if (is.null(events) || nrow(events) == 0) {
    cat("  [acled] No events returned\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  cat(sprintf("  [acled] Downloaded %d events\n", nrow(events)))

  # ACLED has admin1 and admin2 columns — use them directly
  admin2_names_target <- admin2_sf[["Admin2"]]

  # Try matching ACLED admin2 to our Admin2 names
  if ("admin2" %in% colnames(events)) {
    acled_a2 <- unique(events$admin2[!is.na(events$admin2) & events$admin2 != ""])

    # Build fuzzy lookup
    target_clean <- tolower(trimws(admin2_names_target))
    source_clean <- tolower(trimws(acled_a2))

    lookup <- stats::setNames(rep(NA_character_, length(acled_a2)), acled_a2)
    for (i in seq_along(acled_a2)) {
      exact <- which(target_clean == source_clean[i])
      if (length(exact) == 1) lookup[acled_a2[i]] <- admin2_names_target[exact]
    }
    unmatched <- is.na(lookup)
    if (any(unmatched)) {
      um <- acled_a2[unmatched]
      um_clean <- source_clean[unmatched]
      if (length(um) > 0 && length(target_clean) > 0) {
        dist_mat <- utils::adist(um_clean, target_clean, ignore.case = TRUE, partial = FALSE)
        len_mat <- outer(nchar(um_clean), nchar(target_clean), pmax)
        norm_dist <- dist_mat / pmax(len_mat, 1)
        for (i in seq_along(um)) {
          best <- which.min(norm_dist[i, ])
          if (length(best) == 1 && norm_dist[i, best] < 0.3) {
            lookup[um[i]] <- admin2_names_target[best]
          }
        }
      }
    }

    events$Admin2_matched <- lookup[events$admin2]
  } else {
    # Fallback: spatial join using lat/lon
    cat("  [acled] No admin2 column — using spatial join...\n")
    if (all(c("latitude", "longitude") %in% colnames(events))) {
      events_sf <- sf::st_as_sf(events,
                                coords = c("longitude", "latitude"),
                                crs = 4326)
      joined <- sf::st_join(events_sf, admin2_sf[, "Admin2"], left = TRUE)
      events$Admin2_matched <- joined$Admin2
    } else {
      events$Admin2_matched <- NA_character_
    }
  }

  # Aggregate by Admin2
  matched_events <- events[!is.na(events$Admin2_matched), ]

  if (nrow(matched_events) == 0) {
    cat("  [acled] No events matched to Admin-2 areas\n")
    return(.empty_result(admin2_names_target))
  }

  # Ensure fatalities is numeric
  if ("fatalities" %in% colnames(matched_events)) {
    matched_events$fatalities <- as.numeric(matched_events$fatalities)
  }

  agg <- matched_events %>%
    dplyr::group_by(Admin2 = Admin2_matched) %>%
    dplyr::summarise(
      acled_total_events   = dplyr::n(),
      acled_fatalities     = sum(fatalities, na.rm = TRUE),
      acled_battles        = sum(event_type == "Battles", na.rm = TRUE),
      acled_violence_civil = sum(event_type == "Violence against civilians", na.rm = TRUE),
      acled_protests       = sum(event_type == "Protests", na.rm = TRUE),
      acled_riots          = sum(event_type == "Riots", na.rm = TRUE),
      acled_explosions     = sum(event_type == "Explosions/Remote violence", na.rm = TRUE),
      acled_strategic      = sum(event_type == "Strategic developments", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  # Normalize by years in window
  n_years <- window_years + 1  # inclusive
  agg$acled_events_per_year <- agg$acled_total_events / n_years
  agg$acled_fatalities_per_year <- agg$acled_fatalities / n_years

  # Fill missing Admin2 areas with 0 (no conflict = 0 events)
  all_admin2 <- data.frame(Admin2 = admin2_names_target, stringsAsFactors = FALSE)
  result <- merge(all_admin2, agg, by = "Admin2", all.x = TRUE)
  acled_cols <- grep("^acled_", colnames(result), value = TRUE)
  for (col in acled_cols) {
    result[[col]][is.na(result[[col]])] <- 0
  }

  n_with_events <- sum(result$acled_total_events > 0)
  cat(sprintf("  [acled] %d/%d Admin-2 areas have conflict events (%d total events)\n",
              n_with_events, nrow(result), sum(result$acled_total_events)))

  result
}


# =============================================================================
# 10. HarvestStat Africa — Subnational crop statistics
# =============================================================================

#' Extract subnational crop yield data from HarvestStat Africa
#'
#' HarvestStat Africa provides subnational crop production, harvested area,
#' and yields for 33 SSA countries (1980-2022). Data from FEWS NET.
#'
#' Download the CSV from: https://datadryad.org/dataset/doi:10.5061/dryad.vq83bk42w
#' Place in cache_dir as "harveststat_africa.csv"
#'
#' @param admin2_sf sf object with Admin2 polygons (for name matching)
#' @param country_name Country name (e.g., "Ghana", "Gambia")
#' @param survey_year Integer survey year
#' @param cache_dir Cache directory (looks for harveststat_africa.csv here)
#' @return data.frame with Admin2 + crop_* columns
extract_harveststat <- function(admin2_sf, country_name, survey_year, cache_dir) {

  # Look for the CSV in cache_dir or data/raw/
  possible_paths <- c(
    file.path(cache_dir, "hvstat_africa_data_v1.0.csv"),
    file.path(cache_dir, "harveststat_africa.csv"),
    file.path(cache_dir, "HarvestStat_Africa.csv"),
    here::here("data", "raw", "hvstat_africa_data_v1.0.csv"),
    here::here("data", "raw", "harveststat_africa.csv"),
    here::here("data", "raw", "HarvestStat_Africa.csv")
  )

  csv_path <- NULL
  for (p in possible_paths) {
    if (file.exists(p)) { csv_path <- p; break }
  }

  if (is.null(csv_path)) {
    cat("  [harveststat] CSV not found. Download from:\n")
    cat("    https://datadryad.org/dataset/doi:10.5061/dryad.vq83bk42w\n")
    cat("  Place as: ", possible_paths[1], "\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  cat(sprintf("  [harveststat] Loading from %s\n", basename(csv_path)))
  hs <- data.table::fread(csv_path, showProgress = FALSE)

  # Filter to country
  # HarvestStat v1.0 uses "country" column with names like "Ghana", "Malawi"
  # and "country_code" with ISO2 codes like "GH", "MW"
  hs_country_map <- c(
    "Gambia" = "Gambia", "Ghana" = "Ghana",
    "Sierra Leone" = "Sierra Leone", "Malawi" = "Malawi"
  )
  hs_name <- hs_country_map[country_name]
  if (is.na(hs_name)) hs_name <- country_name

  # Find country column
  country_col <- NULL
  for (cc in c("country", "admin0", "Country", "ADMIN0", "admin_0")) {
    if (cc %in% colnames(hs)) { country_col <- cc; break }
  }
  if (is.null(country_col)) {
    cat("  [harveststat] Cannot find country column in CSV\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  hs <- hs[hs[[country_col]] == hs_name, ]
  if (nrow(hs) == 0) {
    cat(sprintf("  [harveststat] No data for %s\n", hs_name))
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Filter to nearest years (survey year or prior)
  # HarvestStat v1.0 has harvest_year and planting_year
  year_col <- NULL
  for (yc in c("harvest_year", "year", "Year", "season_year")) {
    if (yc %in% colnames(hs)) { year_col <- yc; break }
  }
  if (!is.null(year_col)) {
    hs[[year_col]] <- as.integer(hs[[year_col]])
    avail_years <- sort(unique(hs[[year_col]]))
    # Use 3-year average centered on survey year
    target_years <- avail_years[avail_years >= (survey_year - 2) & avail_years <= survey_year]
    if (length(target_years) == 0) {
      # Fallback: nearest available years
      target_years <- avail_years[order(abs(avail_years - survey_year))][1:min(3, length(avail_years))]
    }
    hs <- hs[hs[[year_col]] %in% target_years, ]
    cat(sprintf("  [harveststat] Using years: %s\n", paste(target_years, collapse = ", ")))
  }

  cat(sprintf("  [harveststat] %d records for %s\n", nrow(hs), hs_name))

  # Find admin columns and crop data columns
  # HarvestStat v1.0 uses: admin_1, admin_2, product, yield, production, area
  admin_col <- NULL
  for (ac in c("admin_1", "admin1", "Admin1", "ADMIN1", "region", "Region")) {
    if (ac %in% colnames(hs)) { admin_col <- ac; break }
  }

  admin2_src_col <- NULL
  for (ac2 in c("admin_2", "admin2", "Admin2", "ADMIN2")) {
    if (ac2 %in% colnames(hs)) { admin2_src_col <- ac2; break }
  }

  crop_col <- NULL
  for (crc in c("product", "crop", "Crop", "commodity", "Commodity")) {
    if (crc %in% colnames(hs)) { crop_col <- crc; break }
  }

  yield_col <- NULL
  for (yld in c("yield", "Yield", "yield_mt_ha", "Yield_MT_HA")) {
    if (yld %in% colnames(hs)) { yield_col <- yld; break }
  }

  production_col <- NULL
  for (prd in c("production", "Production", "production_mt", "Production_MT")) {
    if (prd %in% colnames(hs)) { production_col <- prd; break }
  }

  area_col <- NULL
  for (arc in c("area", "harvested_area", "Harvested_Area", "area_ha", "Area_HA")) {
    if (arc %in% colnames(hs)) { area_col <- arc; break }
  }

  if (is.null(admin_col)) {
    cat("  [harveststat] Cannot find admin column\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Aggregate: for each admin area, compute mean yield across major crops
  # Focus on staple crops relevant to nutrition
  staple_crops <- c("maize", "rice", "sorghum", "millet", "wheat",
                    "cassava", "yam", "groundnut", "cowpea", "beans",
                    "Maize", "Rice", "Sorghum", "Millet", "Wheat",
                    "Cassava", "Yam", "Groundnut", "Cowpea", "Beans")

  if (!is.null(crop_col)) {
    # Overall stats (all crops)
    overall <- hs %>%
      dplyr::group_by(admin_area = .data[[admin_col]]) %>%
      dplyr::summarise(
        crop_n_crops     = dplyr::n_distinct(.data[[crop_col]]),
        crop_mean_yield  = if (!is.null(yield_col)) mean(as.numeric(.data[[yield_col]]), na.rm = TRUE) else NA_real_,
        crop_total_production = if (!is.null(production_col)) sum(as.numeric(.data[[production_col]]), na.rm = TRUE) else NA_real_,
        crop_total_area  = if (!is.null(area_col)) sum(as.numeric(.data[[area_col]]), na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      as.data.frame()

    # Staple-specific stats
    hs_staples <- hs[tolower(hs[[crop_col]]) %in% tolower(staple_crops), ]
    if (nrow(hs_staples) > 0) {
      staple_agg <- hs_staples %>%
        dplyr::group_by(admin_area = .data[[admin_col]]) %>%
        dplyr::summarise(
          crop_staple_yield = if (!is.null(yield_col)) mean(as.numeric(.data[[yield_col]]), na.rm = TRUE) else NA_real_,
          crop_staple_n     = dplyr::n_distinct(.data[[crop_col]]),
          .groups = "drop"
        ) %>%
        as.data.frame()
      overall <- merge(overall, staple_agg, by = "admin_area", all.x = TRUE)
    }
  } else {
    overall <- hs %>%
      dplyr::group_by(admin_area = .data[[admin_col]]) %>%
      dplyr::summarise(
        crop_mean_yield = if (!is.null(yield_col)) mean(as.numeric(.data[[yield_col]]), na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      as.data.frame()
  }

  # Match admin names to Admin1 (HarvestStat is typically Admin-1 level)
  # Then broadcast to Admin-2 via the Admin1 column in admin2_sf
  target_a1 <- unique(admin2_sf[["Admin1"]])
  source_a1 <- unique(overall$admin_area)

  # Fuzzy match
  target_clean <- tolower(trimws(target_a1))
  source_clean <- tolower(trimws(source_a1))

  lookup <- stats::setNames(rep(NA_character_, length(source_a1)), source_a1)
  for (i in seq_along(source_a1)) {
    exact <- which(target_clean == source_clean[i])
    if (length(exact) == 1) lookup[source_a1[i]] <- target_a1[exact]
  }
  unmatched <- is.na(lookup)
  if (any(unmatched) && length(target_a1) > 0) {
    um <- source_a1[unmatched]
    um_clean <- source_clean[unmatched]
    dist_mat <- utils::adist(um_clean, target_clean, ignore.case = TRUE, partial = FALSE)
    len_mat <- outer(nchar(um_clean), nchar(target_clean), pmax)
    norm_dist <- dist_mat / pmax(len_mat, 1)
    for (i in seq_along(um)) {
      best <- which.min(norm_dist[i, ])
      if (length(best) == 1 && norm_dist[i, best] < 0.3) {
        lookup[um[i]] <- target_a1[best]
      }
    }
  }

  overall$Admin1 <- lookup[overall$admin_area]
  overall <- overall[!is.na(overall$Admin1), ]
  overall$admin_area <- NULL

  if (nrow(overall) == 0) {
    cat("  [harveststat] No admin areas matched\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Broadcast to Admin-2 via Admin1
  admin2_a1 <- data.frame(
    Admin2 = admin2_sf[["Admin2"]],
    Admin1 = admin2_sf[["Admin1"]],
    stringsAsFactors = FALSE
  )
  result <- merge(admin2_a1, overall, by = "Admin1", all.x = TRUE)
  result$Admin1 <- NULL

  n_matched <- sum(!is.na(result$crop_mean_yield))
  cat(sprintf("  [harveststat] %d/%d Admin-2 areas with crop data\n",
              n_matched, nrow(result)))

  result
}


# =============================================================================
# 11. WFP Staple Food Prices
# =============================================================================

#' Extract WFP staple food commodity prices aggregated to Admin-2
#'
#' Reads WFP food price CSV data, filters to the survey year (± 1 year),
#' aggregates mean USD prices per commodity per Admin-2, and returns a wide
#' data.frame with one column per commodity.
#'
#' @param admin2_sf sf object with Admin-2 polygons (must have Admin2 column)
#' @param country_name Country name as used in config (e.g., "Gambia")
#' @param survey_year Integer survey year
#' @param price_dir Directory containing wfp_food_prices_*.csv files
#' @return data.frame with Admin2 + wfp_price_* columns
extract_wfp_foodprices <- function(admin2_sf, country_name, survey_year,
                                    price_dir = here::here("data", "food_price")) {

  country_codes <- c(
    "Gambia"       = "gmb",
    "Ghana"        = "gha",
    "Sierra Leone" = "sle",
    "Malawi"       = "mwi"
  )

  code <- country_codes[[country_name]]
  if (is.null(code)) {
    cat(sprintf("[WFP prices] No food price CSV mapping for '%s'\n", country_name))
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  csv_path <- file.path(price_dir, paste0("wfp_food_prices_", code, ".csv"))
  if (!file.exists(csv_path)) {
    cat(sprintf("[WFP prices] File not found: %s\n", csv_path))
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  cat(sprintf("[WFP prices] Loading %s\n", basename(csv_path)))
  raw <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

  # Remove HXL metadata row (starts with #)
  raw <- raw[!grepl("^#", raw[[1]]), ]
  cat(sprintf("[WFP prices] %d rows after removing metadata\n", nrow(raw)))

  # Parse year from date column
  raw$year <- as.integer(substr(raw$date, 1, 4))

  # Filter to survey year ± 1 for temporal alignment
  raw <- raw[!is.na(raw$year) &
             raw$year >= (survey_year - 1L) &
             raw$year <= survey_year, ]
  cat(sprintf("[WFP prices] %d rows in %d-%d\n", nrow(raw),
              survey_year - 1L, survey_year))

  if (nrow(raw) == 0) {
    cat("[WFP prices] No data for survey period\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Convert usdprice to numeric
  raw$usdprice <- suppressWarnings(as.numeric(raw$usdprice))
  raw <- raw[!is.na(raw$usdprice) & raw$usdprice > 0, ]

  # Exclude non-food items (fuel, exchange rates, wage labour, etc.)
  if ("category" %in% colnames(raw)) {
    raw <- raw[!grepl("non-food|miscellaneous", raw$category, ignore.case = TRUE), ]
  }

  cat(sprintf("[WFP prices] %d food price observations\n", nrow(raw)))
  if (nrow(raw) == 0) return(.empty_result(admin2_sf[["Admin2"]]))

  # Prioritize actual > aggregate prices, Retail > Wholesale
  if ("priceflag" %in% colnames(raw)) {
    raw <- raw[order(raw$priceflag), ]
  }
  if ("pricetype" %in% colnames(raw)) {
    raw <- raw[order(raw$pricetype), ]
  }

  # Aggregate: mean USD price per admin2 × commodity
  # Use admin2 from WFP CSV; fall back to admin1 if admin2 is missing
  admin_col <- if ("admin2" %in% colnames(raw) &&
                   sum(nchar(trimws(raw$admin2)) > 0) > nrow(raw) / 2) {
    "admin2"
  } else if ("admin1" %in% colnames(raw)) {
    cat("[WFP prices] admin2 sparse — using admin1 for aggregation\n")
    "admin1"
  } else {
    cat("[WFP prices] No admin columns available\n")
    return(.empty_result(admin2_sf[["Admin2"]]))
  }

  # Remove rows with empty admin names
  raw <- raw[nchar(trimws(raw[[admin_col]])) > 0, ]
  if (nrow(raw) == 0) return(.empty_result(admin2_sf[["Admin2"]]))

  # Clean commodity names for column naming (don't use make_clean_names —
  # it deduplicates and would create one unique name per row)
  raw$commodity_clean <- tolower(gsub("[^a-zA-Z0-9]+", "_", raw$commodity))
  raw$commodity_clean <- gsub("_+", "_", raw$commodity_clean)
  raw$commodity_clean <- gsub("^_|_$", "", raw$commodity_clean)
  raw$wfp_admin <- raw[[admin_col]]

  # Average price per admin area × commodity
  agg <- stats::aggregate(
    usdprice ~ wfp_admin + commodity_clean,
    data = raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  # Pivot to wide format
  wide <- tidyr::pivot_wider(
    agg,
    names_from  = commodity_clean,
    values_from = usdprice
  )

  # Add wfp_price_ prefix to commodity columns
  price_cols <- setdiff(colnames(wide), "wfp_admin")
  for (pc in price_cols) {
    colnames(wide)[colnames(wide) == pc] <- paste0("wfp_price_", pc)
  }

  cat(sprintf("[WFP prices] %d admin areas × %d commodities\n",
              nrow(wide), length(price_cols)))

  # ── Match WFP admin names to GADM Admin2 names ──
  gadm_names <- admin2_sf[["Admin2"]]
  gadm_admin1 <- admin2_sf[["Admin1"]]
  wfp_names <- wide$wfp_admin

  # Build name map: GADM name → WFP name
  name_map <- setNames(rep(NA_character_, length(gadm_names)), gadm_names)

  # If WFP uses admin1 aggregation, match via Admin1 instead
  match_target <- if (admin_col == "admin1") gadm_admin1 else gadm_names

  # Exact matches first
  for (i in seq_along(match_target)) {
    nm <- match_target[i]
    if (nm %in% wfp_names) {
      name_map[gadm_names[i]] <- nm
    }
  }

  # Fuzzy matches for unmatched
  unmatched <- gadm_names[is.na(name_map)]
  for (nm in unique(match_target[gadm_names %in% unmatched])) {
    fuzzy <- agrep(nm, wfp_names, max.distance = 0.2, value = TRUE)
    if (length(fuzzy) >= 1) {
      # Pick closest by string length
      best <- fuzzy[which.min(abs(nchar(fuzzy) - nchar(nm)))]
      idx <- which(match_target == nm & is.na(name_map))
      name_map[gadm_names[idx]] <- best
    }
  }

  n_matched <- sum(!is.na(name_map))
  cat(sprintf("[WFP prices] Matched %d/%d GADM Admin-2 units to WFP areas\n",
              n_matched, length(gadm_names)))

  if (n_matched == 0) return(.empty_result(gadm_names))

  # Build result: GADM Admin2 spine with matched WFP prices
  result <- data.frame(Admin2 = gadm_names, stringsAsFactors = FALSE)
  result$.wfp_admin <- name_map[gadm_names]

  wfp_price_cols <- grep("^wfp_price_", colnames(wide), value = TRUE)
  result <- merge(result, wide[, c("wfp_admin", wfp_price_cols), drop = FALSE],
                  by.x = ".wfp_admin", by.y = "wfp_admin",
                  all.x = TRUE, sort = FALSE)
  result$.wfp_admin <- NULL

  result
}


# =============================================================================
# 12. Master function: extract_all_external
# =============================================================================

#' Extract all external predictor variables for a country
#'
#' Loads Admin-2 polygons, calls each extract_* function, and joins all
#' results into a single data.frame. Failures in individual sources are
#' handled gracefully — the other sources still contribute their columns.
#'
#' @param cc Country config list from get_country_configs()
#' @param survey_year Integer survey year
#' @param cache_dir Directory for caching downloads (default: "data/external_cache")
#' @param nasa_token NASA Earthdata bearer token for nightlights (optional)
#' @param wfp_key WFP Data Bridges API key (optional)
#' @return data.frame with Admin2 + all external predictor columns
extract_all_external <- function(cc, survey_year, cache_dir = "data/external_cache",
                                 nasa_token = NULL, wfp_key = NULL) {

  .require_pkg("geodata", "downloading GADM boundaries")
  .require_pkg("sf", "spatial operations")

  cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
  .ensure_dir(cache_dir)

  # Read tokens from environment if not passed explicitly
  if (is.null(nasa_token) || nasa_token == "") {
    nasa_token <- Sys.getenv("NASA_EARTHDATA_TOKEN", unset = "")
    if (nchar(nasa_token) == 0) nasa_token <- NULL
  }
  if (is.null(wfp_key) || wfp_key == "") {
    wfp_key <- Sys.getenv("WFP_API_KEY", unset = "")
    if (nchar(wfp_key) == 0) wfp_key <- NULL
  }

  country_iso3 <- cc$gadm_code
  admin2_col   <- cc$admin2_col  # typically "Admin2"
  admin1_col   <- cc$admin1_col  # typically "Admin1"

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Extracting external predictors: %s (%s), year = %d\n",
              cc$country, country_iso3, survey_year))
  cat(sprintf("Cache directory: %s\n", cache_dir))
  cat(sprintf("========================================\n\n"))

  # ── Load GADM Admin-2 boundaries ─────────────────────────────────────────
  gadm_cache <- file.path(cache_dir, "gadm")
  .ensure_dir(gadm_cache)

  admin2_sf <- tryCatch({
    gadm <- geodata::gadm(country_iso3, level = 2, path = gadm_cache)
    admin2 <- sf::st_as_sf(gadm)

    # Ensure Admin2 name column exists
    # GADM uses NAME_2 for admin level 2 and NAME_1 for admin level 1
    if (!"Admin2" %in% colnames(admin2) && "NAME_2" %in% colnames(admin2)) {
      admin2$Admin2 <- admin2$NAME_2
    }
    if (!"Admin1" %in% colnames(admin2) && "NAME_1" %in% colnames(admin2)) {
      admin2$Admin1 <- admin2$NAME_1
    }

    cat(sprintf("[GADM] Loaded %d Admin-2 polygons for %s\n",
                nrow(admin2), country_iso3))
    admin2
  }, error = function(e) {
    stop(sprintf("Cannot load GADM boundaries for '%s': %s", country_iso3, e$message))
  })

  admin2_names <- admin2_sf[["Admin2"]]
  admin1_names <- unique(admin2_sf[["Admin1"]])

  # ── Call each extraction function with tryCatch ──────────────────────────

  # 1. CHIRPS rainfall
  cat("\n--- [1/11] CHIRPS Rainfall ---\n")
  chirps_df <- tryCatch(
    extract_chirps(admin2_sf, country_iso3, survey_year, cache_dir),
    error = function(e) {
      warning(sprintf("CHIRPS extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 2. WorldPop population density
  cat("\n--- [2/11] WorldPop Population ---\n")
  worldpop_df <- tryCatch(
    extract_worldpop(admin2_sf, country_iso3, survey_year, cache_dir),
    error = function(e) {
      warning(sprintf("WorldPop extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 3. Nighttime lights
  cat("\n--- [3/11] VIIRS Nighttime Lights ---\n")
  ntl_df <- tryCatch(
    extract_nightlights(admin2_sf, country_iso3, survey_year, nasa_token, cache_dir),
    error = function(e) {
      warning(sprintf("Nightlights extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 4. Malaria Atlas Project
  cat("\n--- [4/11] Malaria Atlas Project ---\n")
  map_df <- tryCatch(
    extract_malaria_atlas(admin2_sf, survey_year, cache_dir, country_iso3 = country_iso3),
    error = function(e) {
      warning(sprintf("Malaria Atlas extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 5. SoilGrids
  cat("\n--- [5/11] SoilGrids ---\n")
  soil_df <- tryCatch(
    extract_soilgrids(admin2_sf, cache_dir),
    error = function(e) {
      warning(sprintf("SoilGrids extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 6. Global Data Lab (Admin-1 level)
  cat("\n--- [6/11] Global Data Lab ---\n")
  gdl_df <- tryCatch(
    extract_gdl(country_iso3, admin1_names, survey_year),
    error = function(e) {
      warning(sprintf("GDL extraction failed: %s", e$message))
      data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE)
    }
  )

  # 7. WFP HungerMap (Admin-1 level)
  cat("\n--- [7/11] WFP HungerMap ---\n")
  wfp_df <- tryCatch(
    extract_wfp_hungermap(country_iso3, admin1_names, wfp_key),
    error = function(e) {
      warning(sprintf("WFP extraction failed: %s", e$message))
      data.frame(Admin1 = admin1_names, stringsAsFactors = FALSE)
    }
  )

  # 8. IPC/Cadre Harmonisé API
  cat("\n--- [8/11] IPC/CH Food Security ---\n")
  ipc_df <- tryCatch(
    extract_ipc(country_iso3, admin2_sf, survey_year),
    error = function(e) {
      warning(sprintf("IPC extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 9. ACLED Conflict Data
  cat("\n--- [9/11] ACLED Conflict Events ---\n")
  acled_df <- tryCatch(
    extract_acled(admin2_sf, cc$country, survey_year),
    error = function(e) {
      warning(sprintf("ACLED extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 10. HarvestStat Africa Crop Statistics
  cat("\n--- [10/11] HarvestStat Africa Crop Data ---\n")
  harvest_df <- tryCatch(
    extract_harveststat(admin2_sf, cc$country, survey_year, cache_dir),
    error = function(e) {
      warning(sprintf("HarvestStat extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # 11. WFP Staple Food Prices
  cat("\n--- [11/11] WFP Staple Food Prices ---\n")
  foodprice_df <- tryCatch(
    extract_wfp_foodprices(admin2_sf, cc$country, survey_year),
    error = function(e) {
      warning(sprintf("WFP food prices extraction failed: %s", e$message))
      .empty_result(admin2_names)
    }
  )

  # ── Join all Admin-2 level results ───────────────────────────────────────
  cat("\n--- Joining results ---\n")

  # Start with Admin-2 names as the spine
  combined <- data.frame(
    Admin2 = admin2_names,
    Admin1 = admin2_sf[["Admin1"]],
    stringsAsFactors = FALSE
  )

  # Join Admin-2 level data — deduplicate sources first to prevent cartesian joins
  admin2_sources <- list(chirps_df, worldpop_df, ntl_df, map_df, soil_df,
                         ipc_df, acled_df, harvest_df, foodprice_df)
  for (src_df in admin2_sources) {
    if (is.null(src_df) || !is.data.frame(src_df)) next
    if (!("Admin2" %in% colnames(src_df)) || ncol(src_df) <= 1) next
    if (nrow(src_df) == 0) next

    # Deduplicate by Admin2 — take first row for each (some sources may have dupes)
    src_df <- src_df[!duplicated(src_df$Admin2), , drop = FALSE]
    extra_cols <- setdiff(colnames(src_df), "Admin2")
    combined <- merge(combined, src_df[, c("Admin2", extra_cols), drop = FALSE],
                      by = "Admin2", all.x = TRUE, sort = FALSE)

    # Safety check: combined should never grow beyond original spine
    if (nrow(combined) > length(admin2_names)) {
      cat(sprintf("  WARNING: Join inflated rows from %d to %d — deduplicating\n",
                  length(admin2_names), nrow(combined)))
      combined <- combined[!duplicated(combined$Admin2), , drop = FALSE]
    }
  }

  # Join Admin-1 level data (GDL, WFP) — broadcast to Admin-2
  admin1_sources <- list(gdl_df, wfp_df)
  for (src_df in admin1_sources) {
    if (is.null(src_df) || !is.data.frame(src_df)) next
    if (!("Admin1" %in% colnames(src_df)) || ncol(src_df) <= 1) next
    if (nrow(src_df) == 0) next

    src_df <- src_df[!duplicated(src_df$Admin1), , drop = FALSE]
    extra_cols <- setdiff(colnames(src_df), "Admin1")
    combined <- merge(combined, src_df[, c("Admin1", extra_cols), drop = FALSE],
                      by = "Admin1", all.x = TRUE, sort = FALSE)

    if (nrow(combined) > length(admin2_names)) {
      cat(sprintf("  WARNING: Admin1 join inflated rows — deduplicating\n"))
      combined <- combined[!duplicated(combined$Admin2), , drop = FALSE]
    }
  }

  # ── Report summary ──────────────────────────────────────────────────────
  predictor_cols <- setdiff(colnames(combined), c("Admin1", "Admin2"))
  n_complete <- sum(stats::complete.cases(combined[, predictor_cols, drop = FALSE]))
  n_any_data <- sum(rowSums(!is.na(combined[, predictor_cols, drop = FALSE])) > 0)

  cat(sprintf("\n========================================\n"))
  cat(sprintf("External predictor extraction complete\n"))
  cat(sprintf("  Admin-2 units:  %d\n", nrow(combined)))
  cat(sprintf("  Total predictors: %d\n", length(predictor_cols)))
  cat(sprintf("  Units with all data: %d\n", n_complete))
  cat(sprintf("  Units with any data: %d\n", n_any_data))
  cat(sprintf("  Predictor columns:\n"))

  # Group by prefix
  prefixes <- unique(sub("_.*", "", predictor_cols))
  for (pfx in prefixes) {
    pfx_cols <- predictor_cols[startsWith(predictor_cols, pfx)]
    n_non_na <- sum(colSums(!is.na(combined[, pfx_cols, drop = FALSE])) > 0)
    cat(sprintf("    %s: %d cols (%d with data)\n", pfx, length(pfx_cols), n_non_na))
  }
  cat(sprintf("========================================\n\n"))

  combined
}


# =============================================================================
# Target-compatible wrapper
# =============================================================================

#' Target-compatible wrapper for external predictor extraction
#'
#' Designed to be called from _targets.R as a target function. Extracts the
#' country-specific survey year from params and calls extract_all_external().
#'
#' @param cc Country config list from get_country_configs()
#' @param params Pipeline params from get_pipeline_params()
#' @param cache_dir Directory for caching downloads
#' @param nasa_token NASA Earthdata bearer token (optional)
#' @param wfp_key WFP Data Bridges API key (optional)
#' @return data.frame with Admin2 + all external predictor columns
external_predictors_target <- function(cc, params, cache_dir = "data/external_cache",
                                       nasa_token = NULL, wfp_key = NULL) {

  # Get survey year from config (centralized) or fall back to hardcoded
  survey_year <- cc$survey_year
  if (is.null(survey_year)) {
    survey_years <- c(
      "Gambia"       = 2018L,
      "Ghana"        = 2017L,
      "Sierra Leone" = 2013L,
      "Malawi"       = 2016L
    )
    survey_year <- survey_years[cc$country]
  }
  if (is.na(survey_year) || is.null(survey_year)) {
    warning(sprintf("No survey year defined for country '%s'. Using 2018 as default.",
                    cc$country))
    survey_year <- 2018L
  }

  cat(sprintf("[external_predictors_target] %s (survey year %d, mode = %s)\n",
              cc$country, survey_year, params$mode))

  # Read tokens from environment variables if not provided
  if (is.null(nasa_token)) {
    nasa_token <- Sys.getenv("NASA_EARTHDATA_TOKEN", unset = "")
    if (nchar(nasa_token) == 0) nasa_token <- NULL
  }
  if (is.null(wfp_key)) {
    wfp_key <- Sys.getenv("WFP_API_KEY", unset = "")
    if (nchar(wfp_key) == 0) wfp_key <- NULL
  }

  extract_all_external(
    cc          = cc,
    survey_year = survey_year,
    cache_dir   = cache_dir,
    nasa_token  = nasa_token,
    wfp_key     = wfp_key
  )
}
