# =============================================================================
# src/GEE/extract_gee_admin2.R
#
# Generic Admin-2 zonal-mean extractor for GEE-style raster files.
# Used by scripts/build_gee_admin2.R to pre-compute a wide per-country CSV
# (one row per Admin-2 polygon × one column per raster band) which the
# country merge scripts then left-join on Admin2.
#
# All output columns are prefixed `gee_a2_` so they enter the GEE proxy-
# predictor domain (R/config.R: GEE = list(prefix = "gee_")) automatically.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
})

#' Detect whether a raster's bands are temporal sub-divisions of a year.
#' When TRUE, the extractor reduces multi-band output to a single annual-mean
#' column per polygon (representing the survey-year average).
#'
#' Three signals, OR'd:
#'   1. terra::time() returns multiple dates -> temporal by metadata.
#'   2. Filename contains monthly / daily / 8day / hourly.
#'   3. Band count > 31 (more than 1 month of dailies) and filename does NOT
#'      look like a variable stack (TerraClimate, LandCover, GPW_Demographic,
#'      FLDAS_Annual_Mean, etc.) — heuristic for rasters with unlabelled
#'      daily bands such as Atmosphere_<Year>.
is_temporal_multiband <- function(r, fname) {
  bn <- basename(fname)

  # Signal 1: terra time metadata
  tt <- tryCatch(terra::time(r), error = function(e) NULL)
  if (!is.null(tt) && length(tt) > 1 && !all(is.na(tt))) return(TRUE)

  # Signal 2: filename keyword
  if (grepl("monthly|daily|8day|hourly", bn, ignore.case = TRUE)) return(TRUE)

  # Signal 3: many bands and not a known variable-stack family
  variable_stack_patterns <- c(
    "TerraClimate", "LandCover", "ESA_WorldCereal", "GPW_Demographic",
    "FLDAS.*Annual_Mean", "WAPOR", "Productivity"
  )
  is_var_stack <- any(vapply(variable_stack_patterns,
                              function(p) grepl(p, bn, ignore.case = TRUE),
                              logical(1)))
  if (!is_var_stack && terra::nlyr(r) > 31L) return(TRUE)

  FALSE
}

#' Extract Admin-2 zonal mean from a set of raster files.
#'
#' For each raster:
#' - Single-band  -> one column.
#' - Multi-band where filename suggests temporal subdivision (Monthly/Daily/
#'   8days) -> aggregated to one annual-mean column per polygon (mean across
#'   bands).
#' - Multi-band otherwise (e.g. TerraClimate variable stacks, LandCover class
#'   layers, GPW age bins) -> one column per band; band names used as
#'   suffix where available, else b01..bN.
#'
#' @param admin2_sf      sf object with an `Admin2` column. One row per polygon.
#' @param raster_files   character vector of full paths to .tif files.
#' @param country_tokens character vector of strings to strip from raster
#'   filenames before forming column names (e.g. c("Gambia"), or
#'   c("Sierra_Leone", "Sierra Leone", "SierraLeone")).
#' @param verbose        if TRUE, prints per-raster progress.
#' @return data.frame with `Admin2` + extracted columns, all prefixed
#'   `gee_a2_<cleaned_family>[_<band_suffix>]`.
extract_gee_admin2 <- function(admin2_sf, raster_files,
                                country_tokens = character(0),
                                verbose = TRUE) {

  stopifnot(inherits(admin2_sf, "sf"), "Admin2" %in% names(admin2_sf))

  poly_vect <- terra::vect(admin2_sf)
  # Preserve Admin1 alongside Admin2 so country merge scripts can join on the
  # (Admin1, Admin2) pair — needed for countries (e.g. Malawi) where Admin2
  # names are not unique without their parent Admin1 (TA names like
  # "TA Lundu", "Lake Malawi" appear in multiple districts).
  out <- if ("Admin1" %in% names(admin2_sf)) {
    data.frame(Admin1 = trimws(as.character(admin2_sf$Admin1)),
               Admin2 = trimws(as.character(admin2_sf$Admin2)),
               stringsAsFactors = FALSE)
  } else {
    data.frame(Admin2 = trimws(as.character(admin2_sf$Admin2)),
               stringsAsFactors = FALSE)
  }

  clean_family <- function(fname) {
    family <- sub("\\.tif$", "", basename(fname), ignore.case = TRUE)
    for (tok in country_tokens) {
      family <- gsub(paste0("[_ ]*", gsub(" ", "[_ ]", tok), "[_ ]*"),
                     "_", family, ignore.case = TRUE)
    }
    family <- gsub("__+", "_", family)
    family <- gsub("^_|_$", "", family)
    family <- gsub("[^A-Za-z0-9_]", "_", family)
    family
  }

  for (f in raster_files) {
    if (!file.exists(f)) next
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r)) {
      if (verbose) message("  skip (read error): ", basename(f))
      next
    }

    # Reproject polygons to raster CRS if needed (avoids reprojecting the
    # raster, which is much more expensive).
    pv <- tryCatch(
      if (!terra::same.crs(poly_vect, r))
        terra::project(poly_vect, terra::crs(r)) else poly_vect,
      error = function(e) poly_vect
    )

    vals <- tryCatch(
      terra::extract(r, pv, fun = mean, na.rm = TRUE,
                     weights = TRUE, ID = FALSE),
      error = function(e1) {
        tryCatch(
          terra::extract(r, pv, fun = mean, na.rm = TRUE, ID = FALSE),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(vals) || nrow(vals) != nrow(admin2_sf)) {
      if (verbose) message("  skip (extract error/mismatch): ", basename(f))
      next
    }

    family <- clean_family(f)
    nbands <- ncol(vals)

    temporal <- is_temporal_multiband(r, f)
    if (nbands == 1L) {
      colnames(vals) <- paste0("gee_a2_", family)
    } else if (temporal) {
      # Temporal sub-bands — collapse to annual mean per polygon so the
      # column represents the survey-year average.
      annual_mean <- rowMeans(as.matrix(vals), na.rm = TRUE)
      vals <- data.frame(annual_mean = annual_mean)
      colnames(vals) <- paste0("gee_a2_", family, "_annual_mean")
    } else {
      # Non-temporal multi-band (variable stack) — keep per-band columns.
      bn <- names(r)
      if (length(bn) == nbands && all(nzchar(bn)) && !any(duplicated(bn))) {
        suffix <- gsub("[^A-Za-z0-9]", "_", bn)
      } else {
        suffix <- sprintf("b%02d", seq_len(nbands))
      }
      colnames(vals) <- paste0("gee_a2_", family, "_", suffix)
    }

    # If columns collide with existing, suffix with a counter
    dupes <- intersect(colnames(out), colnames(vals))
    if (length(dupes) > 0) {
      for (d in dupes) {
        colnames(vals)[colnames(vals) == d] <-
          paste0(d, "_dup", sum(grepl(paste0("^", d), colnames(out))))
      }
    }
    out <- cbind(out, vals)
    if (verbose) {
      tag <- if (temporal) sprintf("temporal:%d->1", nbands)
             else if (nbands == 1L) "1"
             else sprintf("stack:%d", nbands)
      message(sprintf("  + %-60s [%s]", basename(f), tag))
    }
  }
  out
}
