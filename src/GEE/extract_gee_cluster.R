# =============================================================================
# src/GEE/extract_gee_cluster.R
#
# Point-level (cluster/EA) zonal-mean extractor for the SAME local GEE raster
# files used by extract_gee_admin2.R — but the zonal geometry is a small BUFFER
# around each survey cluster's GPS centroid instead of the Admin-2 polygon.
# This recovers within-district covariate variation that the Admin-2 aggregation
# throws away (see archive/docs/ghana_geostat_pilot.md).
#
# Output columns are prefixed `gee_cl_` (one row per cluster x raster band).
# Buffer radius: 2 km urban / 5 km rural by default (robust whether the EA GPS
# is exact — these GroundWork surveys document no displacement — or modestly
# displaced). No new GEE export is needed: same rasters, different geometry.
# =============================================================================

suppressPackageStartupMessages({ library(sf); library(terra) })
`%||%` <- function(a, b) if (is.null(a)) b else a
# reuse is_temporal_multiband() from the Admin-2 extractor (load if not present)
if (!exists("is_temporal_multiband", mode = "function")) {
  .src <- if (requireNamespace("here", quietly = TRUE))
            here::here("src", "GEE", "extract_gee_admin2.R") else "src/GEE/extract_gee_admin2.R"
  if (file.exists(.src)) source(.src)
}

#' UTM (or southern-UTM) EPSG for a metric buffer, from data centroid.
utm_epsg <- function(lon, lat) {
  zone <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  (if (mean(lat, na.rm = TRUE) >= 0) 32600 else 32700) + zone
}

#' Extract point-level (buffered-cluster) zonal means from raster files.
#'
#' @param cluster_sf   sf POINT object; must have `id_col`. Optional `Admin1`,
#'   `Admin2` (carried through for joining/validation) and `urban_col`.
#' @param raster_files character vector of .tif paths (same set as Admin-2).
#' @param id_col       cluster id column name (default "cl").
#' @param urban_col    optional urban/rural flag column; if given, buffers are
#'   `buffer_urban`/`buffer_rural`, else a uniform `buffer_m`.
#' @param country_tokens strings to strip from raster filenames (as Admin-2).
#' @return data.frame: id (+ Admin1/Admin2 if present) + `gee_cl_*` columns.
extract_gee_cluster <- function(cluster_sf, raster_files, id_col = "cl",
                                buffer_m = 3000, urban_col = NULL,
                                buffer_urban = 2000, buffer_rural = 5000,
                                country_tokens = character(0), verbose = TRUE) {
  stopifnot(inherits(cluster_sf, "sf"), id_col %in% names(cluster_sf))
  if (any(sf::st_geometry_type(cluster_sf) != "POINT"))
    cluster_sf <- sf::st_centroid(cluster_sf)

  # --- buffer points in a local metric CRS, then back to lon/lat -------------
  xy  <- sf::st_coordinates(sf::st_transform(cluster_sf, 4326))
  epsg <- utm_epsg(xy[, 1], xy[, 2])
  radius <- if (!is.null(urban_col) && urban_col %in% names(cluster_sf)) {
    ifelse(tolower(as.character(cluster_sf[[urban_col]])) %in% c("urban", "u", "1", "yes"),
           buffer_urban, buffer_rural)
  } else rep(buffer_m, nrow(cluster_sf))
  buf <- sf::st_buffer(sf::st_transform(cluster_sf, epsg), dist = radius)
  poly_vect <- terra::vect(buf)
  if (verbose) message(sprintf("Buffered %d clusters (EPSG %d, radius %d-%d m)",
                               nrow(cluster_sf), epsg, min(radius), max(radius)))

  out <- data.frame(.cl = cluster_sf[[id_col]], stringsAsFactors = FALSE)
  names(out)[1] <- id_col
  for (a in intersect(c("Admin1", "Admin2"), names(cluster_sf)))
    out[[a]] <- trimws(as.character(cluster_sf[[a]]))

  clean_family <- function(fname) {
    fam <- sub("\\.tif$", "", basename(fname), ignore.case = TRUE)
    for (tok in country_tokens)
      fam <- gsub(paste0("[_ ]*", gsub(" ", "[_ ]", tok), "[_ ]*"), "_", fam, ignore.case = TRUE)
    fam <- gsub("__+", "_", fam); fam <- gsub("^_|_$", "", fam)
    gsub("[^A-Za-z0-9_]", "_", fam)
  }

  for (f in raster_files) {
    if (!file.exists(f)) next
    r <- tryCatch(terra::rast(f), error = function(e) NULL); if (is.null(r)) next
    pv <- tryCatch(if (!terra::same.crs(poly_vect, r))
                     terra::project(poly_vect, terra::crs(r)) else poly_vect,
                   error = function(e) poly_vect)
    vals <- tryCatch(
      terra::extract(r, pv, fun = mean, na.rm = TRUE, weights = TRUE, ID = FALSE),
      error = function(e1) tryCatch(
        terra::extract(r, pv, fun = mean, na.rm = TRUE, ID = FALSE),
        error = function(e2) NULL))
    if (is.null(vals) || nrow(vals) != nrow(cluster_sf)) {
      if (verbose) message("  skip (extract error/mismatch): ", basename(f)); next
    }
    family <- clean_family(f); nbands <- ncol(vals)
    if (nbands == 1L) {
      colnames(vals) <- paste0("gee_cl_", family)
    } else if (is_temporal_multiband(r, f)) {
      vals <- data.frame(annual_mean = rowMeans(as.matrix(vals), na.rm = TRUE))
      colnames(vals) <- paste0("gee_cl_", family, "_annual_mean")
    } else {
      bn <- names(r)
      colnames(vals) <- if (length(bn) == nbands && all(nzchar(bn)) && !any(duplicated(bn)))
        paste0("gee_cl_", family, "_", bn) else paste0("gee_cl_", family, "_b", seq_len(nbands))
    }
    out <- cbind(out, vals)
    if (verbose) message("  ok: ", basename(f), " (", nbands, " band(s))")
  }
  out
}
