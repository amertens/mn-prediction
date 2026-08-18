# =============================================================================
# src/GEE/extract_gee_ee_api.R
#
# Earth Engine API auto-extraction (via rgee) — replaces the manual EE Code
# Editor export. Pulls the SAME layers as the legacy .tif exports (defined in
# src/GEE/gee_layer_manifest.R) directly from the EE catalog and reduces them
# server-side over (a) admin-2 polygons and (b) buffered survey clusters,
# returning data frames written to the CSV paths the merge step expects.
#
# REQUIRES Earth Engine API access (see src/GEE/README_GEE_API_SETUP.md):
#   - a Google Cloud project with the Earth Engine API enabled
#   - rgee + the Python earthengine-api installed (rgee::ee_install() once)
#   - authentication: rgee::ee_Initialize(user = "you@gmail.com", drive = TRUE)
#
# NOTE ON COLUMN NAMING (important): EE band names differ from the legacy
# terra-extracted names in data/GEE/<Country>_<year>_admin2_gee.csv. So the
# columns this produces will NOT be byte-identical to the existing four
# countries. That is fine for WITHIN-Tanzania modelling. For the POOLED / LOCO
# common-vocabulary set, either (a) re-extract all countries through this same
# script so names match, or (b) add a name crosswalk. The admin-2 builder prints
# a parity report vs Sierra Leone so you can see the overlap.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(here)
})
source(here("src", "GEE", "gee_layer_manifest.R"))

#' Build the reduced ee$Image for one manifest layer at a target year.
#' @param layer a GEE_LAYER_MANIFEST entry.
#' @param year  survey year (used for filter_year / nearest_epoch modes).
#' @return an ee$Image (single- or multi-band), or NULL on failure.
ee_image_for_layer <- function(layer, year) {
  ee <- rgee::ee
  # Select the requested bands. Do NOT swallow a bad band name — let it error so
  # the outer tryCatch skips the layer, rather than silently keeping ALL bands
  # (which is how Atmosphere ballooned to 628 columns).
  sel <- function(x) if (!is.null(layer$bands)) x$select(layer$bands) else x
  # A layer is a collection if kind says so OR its temporal mode requires one
  # (filter_year / nearest_epoch). Some manifest entries set kind to a temporal
  # value (e.g. "nearest_epoch") — this keeps them working either way.
  is_collection <- identical(layer$kind, "ImageCollection") ||
    (!is.null(layer$temporal) && layer$temporal %in% c("filter_year", "nearest_epoch"))

  # Build the band-selected, reduced, scaled image for a specific year.
  build_year_img <- function(yr) {
    if (is_collection) {
      ic <- ee$ImageCollection(layer$asset)
      if (identical(layer$temporal, "filter_year")) {
        ic <- ic$filterDate(sprintf("%d-01-01", yr), sprintf("%d-12-31", yr))
      } else if (identical(layer$temporal, "nearest_epoch")) {
        epoch <- max(layer$avail[1], min(layer$avail[2], round(yr / 5) * 5))
        ic <- ic$filterDate(sprintf("%d-01-01", epoch), sprintf("%d-12-31", epoch + 4L))
      }
      ic <- sel(ic)   # select bands on the collection before reducing
      im <- if (isTRUE(layer$to_bands)) ee$Image(ic$toBands()) else ee$Image(ic$mean())
    } else {
      im <- sel(ee$Image(layer$asset))
    }
    if (!is.null(layer$scale_factor)) im <- im$multiply(layer$scale_factor)
    im
  }

  img <- tryCatch(build_year_img(gee_nearest_year(layer, year)),
                  error = function(e) { message("  [", layer$family, "] build failed: ",
                                                conditionMessage(e)); NULL })
  if (is.null(img)) return(NULL)

  # ── Year-over-year change band(s): current - prior year ──────────────────
  # Recovers the legacy *_yoy features (e.g. LAI8days_yoy was a top-5 predictor).
  # Only for filter_year layers flagged yoy = TRUE (nearest_epoch would snap both
  # years to the same 5-yr epoch -> a useless zero delta).
  if (isTRUE(layer$yoy) && identical(layer$temporal, "filter_year")) {
    yr <- gee_nearest_year(layer, year)
    if (yr - 1L >= layer$avail[1]) {
      prev <- tryCatch(build_year_img(yr - 1L), error = function(e) NULL)
      if (!is.null(prev)) {
        bn <- tryCatch(img$bandNames()$getInfo(), error = function(e) NULL)
        if (!is.null(bn)) {
          yoy <- img$subtract(prev)$rename(paste0(bn, "_yoy"))
          img <- tryCatch(img$addBands(yoy), error = function(e) img)
        }
      }
    }
  }
  img
}

#' ee_extract in batches to stay under EE's 10 MB request-payload limit.
#' Splits features into chunks, extracts each, and rbinds in order. Returns NULL
#' if any chunk errors (the caller then skips that layer).
.ee_extract_batched <- function(img, geom_sf, reducer, scale, batch, verbose,
                                label = "") {
  n <- nrow(geom_sf)
  idx <- split(seq_len(n), ceiling(seq_len(n) / batch))
  parts <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    if (verbose) message(sprintf("    %-20s batch %d/%d ...", label, k, length(idx)))
    part <- tryCatch(
      rgee::ee_extract(x = img, y = geom_sf[idx[[k]], ], fun = reducer,
                       scale = scale, sf = FALSE),
      error = function(e) { message("    batch ", k, "/", length(idx),
                                    " failed: ", conditionMessage(e)); NULL })
    if (is.null(part)) return(NULL)
    parts[[k]] <- part
  }
  # bind_rows (not rbind): batches can return different columns when ee_extract
  # drops a band that is fully masked over one batch's districts. bind_rows
  # aligns by name and fills the gaps with NA, preserving row order.
  dplyr::bind_rows(parts)
}

#' Extract manifest layers over an sf geometry (polygons or buffered points).
#'
#' @param geom_sf     sf with an id column `id_col` (+ optionally Admin1/Admin2).
#' @param year        survey year.
#' @param manifest    subset of GEE_LAYER_MANIFEST to run (default: all).
#' @param id_col      feature id column carried through.
#' @param col_prefix  prefix for value columns (e.g. "gee_a2_" or "gee_cl_").
#' @param simplify_deg Douglas-Peucker tolerance (decimal degrees) applied to the
#'   geometry once, up front, to shrink the request payload. ~0.01 deg (~1 km) is
#'   negligible for admin-2 zonal means. Set 0/NULL to disable.
#' @param batch       features per ee_extract request (keeps each < 10 MB and
#'   under getInfo's ~5000-element cap: batch x max_bands should stay < 5000).
#' @param scale_floor minimum reduction scale (m). A district/buffer MEAN is
#'   unchanged by sub-250 m resolution, but native scales (soil 30 m, cropland
#'   10 m, WSF 10 m) make the server scan ~70x more pixels. Floors every layer's
#'   scale to this, hugely speeding the high-res layers. Set 0 to use native.
#' @param cache_dir if set, each successfully-extracted layer is cached to
#'   `<cache_dir>/<family>.rds`; on re-run cached layers are loaded (not
#'   re-extracted), so a re-run only pulls the missing/failed layers. Failed
#'   layers are NOT cached, so they retry. Delete the dir (or force=TRUE) to
#'   rebuild from scratch.
#' @param force re-extract every layer, ignoring (and refreshing) the cache.
#' @return data.frame: id (+ Admin1/Admin2) + one column per layer-band.
extract_gee_layers <- function(geom_sf, year, manifest = GEE_LAYER_MANIFEST,
                               id_col = "Admin2", col_prefix = "gee_a2_",
                               verbose = TRUE, simplify_deg = 0.01, batch = 40,
                               scale_floor = 250, cache_dir = NULL, force = FALSE) {
  stopifnot(inherits(geom_sf, "sf"), id_col %in% names(geom_sf))
  ee <- rgee::ee

  # Simplify once up front — GADM admin-2 borders are far more detailed than
  # zonal means require, and the raw GeoJSON exceeds EE's 10 MB request limit.
  if (!is.null(simplify_deg) && simplify_deg > 0) {
    geom_sf <- suppressWarnings(
      sf::st_make_valid(sf::st_simplify(geom_sf, dTolerance = simplify_deg,
                                        preserveTopology = TRUE)))
  }

  keep <- intersect(c(id_col, "Admin1", "Admin2"), names(geom_sf))
  out  <- sf::st_drop_geometry(geom_sf)[, keep, drop = FALSE]

  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  for (layer in manifest) {
    cache_file <- if (!is.null(cache_dir))
      file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_]", "_", layer$family), ".rds")) else NULL

    # ── Load from cache (order-based; verified by id column) ────────────────
    if (!force && !is.null(cache_file) && file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (!is.null(cached) && nrow(cached) == nrow(out) && id_col %in% names(cached) &&
          identical(as.character(cached[[id_col]]), as.character(out[[id_col]]))) {
        add <- cached[, setdiff(names(cached), id_col), drop = FALSE]
        out <- cbind(out, add)
        if (verbose) message(sprintf("  cached: %-20s (%d col[s])", layer$family, ncol(add)))
        next
      }
    }

    img <- ee_image_for_layer(layer, year)
    if (is.null(img)) next
    eff_scale <- max(layer$scale, scale_floor)
    vals <- .ee_extract_batched(img, geom_sf, ee$Reducer$mean(),
                                eff_scale, batch, verbose, label = layer$family)
    if (is.null(vals) || nrow(vals) != nrow(geom_sf)) {
      if (verbose) message("  [", layer$family, "] skipped (extract failed)")
      next  # not cached -> retried next run
    }
    # Drop carried-through id/admin cols ee_extract echoes back; keep value cols.
    val_cols <- setdiff(names(vals), c(keep, id_col, "Admin1", "Admin2"))
    if (length(val_cols) == 0) next
    newnames <- paste0(col_prefix, layer$family,
                       if (length(val_cols) > 1) paste0("_", val_cols) else "")
    add <- vals[, val_cols, drop = FALSE]; names(add) <- newnames

    if (!is.null(cache_file))  # cache id + this layer's columns for fast re-runs
      saveRDS(cbind(out[, id_col, drop = FALSE], add), cache_file)

    out <- cbind(out, add)
    if (verbose) message(sprintf("  ok: %-22s (%d band[s], scale %dm)",
                                  layer$family, length(val_cols), eff_scale))
  }
  out
}

#' Buffer survey-cluster points at a set of radii (km) and extract each.
#' Produces a wide df with a `_<r>km` suffix per layer-band per radius — the
#' cluster-buffer product (legacy gee_ cluster columns).
#'
#' @param cluster_sf sf POINTs with `id_col` (cluster id).
#' @param year       survey year.
#' @param radii_km   buffer radii (default 10/25/50 km, matching legacy files).
extract_gee_cluster_buffers <- function(cluster_sf, year, id_col = "cnum",
                                        radii_km = c(10, 25, 50),
                                        manifest = GEE_LAYER_MANIFEST,
                                        verbose = TRUE, cache_dir = NULL,
                                        force = FALSE) {
  stopifnot(inherits(cluster_sf, "sf"), id_col %in% names(cluster_sf))
  xy <- sf::st_coordinates(sf::st_transform(cluster_sf, 4326))
  zone <- floor((mean(xy[, 1], na.rm = TRUE) + 180) / 6) + 1
  epsg <- (if (mean(xy[, 2], na.rm = TRUE) >= 0) 32600 else 32700) + zone
  base <- sf::st_drop_geometry(cluster_sf)[, id_col, drop = FALSE]
  res <- base
  for (rk in radii_km) {
    buf <- sf::st_buffer(sf::st_transform(cluster_sf, epsg), dist = rk * 1000)
    buf <- sf::st_transform(buf, 4326)
    if (verbose) message(sprintf("== buffer %d km ==", rk))
    # Per-radius cache subdir so the same family cached at 10 km isn't reused at 25 km.
    rk_cache <- if (!is.null(cache_dir)) file.path(cache_dir, paste0("buf", rk, "km")) else NULL
    df <- extract_gee_layers(buf, year, manifest = manifest, id_col = id_col,
                             col_prefix = "gee_", verbose = verbose,
                             cache_dir = rk_cache, force = force)
    val_cols <- setdiff(names(df), id_col)
    names(df)[match(val_cols, names(df))] <- paste0(val_cols, "_", rk, "km")
    res <- dplyr::left_join(res, df, by = id_col)
  }
  res
}
