# =============================================================================
# src/GEE/extract_gee_legacy_parity.R
#
# Earth Engine extraction that reproduces the LEGACY admin-2 gee_* column names.
#
# WHY THIS EXISTS
# ---------------
# The four original countries get their admin-2 GEE covariates from local .tif
# exports in data/<Country>_GEE_rasters/, zonal-averaged by
# .append_gee_zonal_cols() in R/admin2_analysis.R. That function derives each
# column name from the raster FILENAME, so the shared cross-country vocabulary
# is a set of names like `gee_ndvi_2013`, `gee_popdensity_2016`,
# `gee_soilzinc_mean_0_20`, `gee_gpw_demographic_2010_annual_mean`.
#
# A new country added through the Earth Engine API (src/GEE/extract_gee_ee_api.R)
# gets EE band names instead (`gee_a2_NDVI`, `gee_a2_SoilZinc_mean_0_20`, ...).
# Those overlap the legacy vocabulary ZERO percent, so the new country silently
# contributes no covariates to any pooled/LOCO analysis -- and, worse, collapses
# the pooled intersection to nothing for every other country too.
#
# This script closes that gap WITHOUT re-exporting rasters: it pulls the same EE
# assets and emits the columns under the legacy names, so a new country joins
# the existing cross-country vocabulary directly.
#
# The spec below is derived from the intersection of the four legacy countries'
# gee_admin2_* targets (149 columns as of 2026-08). Regenerate the reference
# list and check parity with scripts/check_gee_legacy_parity.R.
#
# Usage:
#   source(here::here("src/GEE/extract_gee_legacy_parity.R"))
#   df <- extract_legacy_parity_admin2("TZA", verbose = TRUE,
#           cache_dir = here::here("data/GEE/.cache_parity_TZA"))
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(here)
})
source(here("src", "GEE", "gee_layer_manifest.R"))
source(here("src", "GEE", "extract_gee_ee_api.R"))

# iSDAsoil rasters were exported with both the depth means and their standard
# deviations, so the legacy zonal step produced four bands per soil property.
SOIL_PARITY_BANDS <- c("mean_0_20", "mean_20_50", "stdev_0_20", "stdev_20_50")

SOIL_PARITY_FAMILIES <- c(
  SoilAluminium  = "soilaluminium",  SoilCalcium     = "soilcalcium",
  SoilCEC        = "soilcec",        SoilIron        = "soiliron",
  SoilMagnesium  = "soilmagnesium",  SoilNitrogen    = "soilnitrogen",
  SoilPhosphorus = "soilphosphorus", SoilPotassium   = "soilpotassium",
  SoilSulfur     = "soilsulfur",     SoilTotalCarbon = "soiltotalcarbon",
  SoilZinc       = "soilzinc"
)

# Order matters: .append_gee_zonal_cols() emits mean, sd, min, max, range.
.parity_summary_suffixes <- function()
  c("_annual_mean", "_annual_sd", "_annual_min", "_annual_max", "_annual_range")

#' Spec of every legacy admin-2 column this extractor reproduces.
#'
#' `mode` controls how a layer's extracted bands become columns, mirroring
#' .append_gee_zonal_cols():
#'   "year"      one extraction per year          -> gee_<stem>_<year>
#'   "single"    one extraction                   -> gee_<stem>
#'   "collapse"  many bands -> gee_<stem>_annual_mean + _annual_sd only
#'               (the legacy >24-band path)
#'   "summaries" many bands -> the five gee_<stem>_annual_* summaries
#'               (per-band names are country-specific, so they are dropped)
#'   "bands"     per-band gee_<stem>_<band> columns PLUS the five summaries
#'               (the legacy <=24-band path)
#' Families deliberately EXCLUDED from the parity set.
#'
#' Each was extracted, compared against the raster-derived values for Sierra
#' Leone, and could not be reproduced on the same scale. Shipping a column whose
#' values are on a different scale than the other countries' is worse than
#' omitting it: it manufactures a country effect that the pooled model reads as
#' real signal. See results/sensitivity/gee_legacy_parity_validation_*.csv.
#' Columns deliberately EXCLUDED from the parity set.
#'
#' Each was extracted, compared against the raster-derived values for Sierra
#' Leone, and could not be reproduced on the same scale. Shipping a column whose
#' values sit on a different scale than the other countries' is worse than
#' omitting it: it manufactures a country effect the pooled model reads as real
#' signal. See results/sensitivity/gee_legacy_parity_validation_*.csv.
GEE_PARITY_EXCLUDED <- c(
  ESA_WorldCereal = paste("The EE classification band is masked outside cropland,",
                          "so a zonal mean returns the class code rather than the",
                          "legacy cropland fraction, and the legacy per-tile bands",
                          "are country-specific (only the across-band summaries",
                          "ever intersected)."),
  Accessibility = paste("Oxford/MAP/accessibility_to_healthcare_2019 is flagged",
                        "corrupted server-side by Earth Engine."),
  SoilTotalCarbon_spread = paste("carbon_total's stdev bands do not match the legacy",
                                 "values under the iSDAsoil back-transform that works",
                                 "for the other ten soil properties (~0.25x), so only",
                                 "its two depth means are kept.")
)

#' Soil properties whose stdev bands (and the summaries derived from them) are
#' dropped -- see GEE_PARITY_EXCLUDED.
SOIL_PARITY_MEANS_ONLY <- "SoilTotalCarbon"

gee_legacy_parity_spec <- function() {
  # The legacy .tif exports carry the base band ONLY, and were exported at the
  # raster resolutions visible in data/<Country>_GEE_rasters/ (NDVI 0.0449 deg
  # ~ 5 km, TRMM 0.0898 deg ~ 10 km). Two overrides matter for both:
  #   yoy = FALSE      the manifest's year-over-year delta band is an addition
  #                    for the modern extractor; averaging it into the base band
  #                    halves the value (this is what made TRMM read ~0.5x).
  #   scale_factor     NDVI rasters hold the RAW integer NDVI (mean ~2900), not
  #                    the 1e-4-scaled reflectance.
  # With those set, both reproduce the legacy values at r = 1.000, ratio = 1.00.
  spec <- list(
    list(family = "NDVI", stem = "ndvi", mode = "year", years = 2010:2022,
         override = list(yoy = FALSE, scale_factor = NULL, scale = 5000,
                         avail = c(1981L, 2022L))),
    list(family = "TRMM", stem = "trmm", mode = "year", years = 2010:2019,
         override = list(yoy = FALSE, scale = 10000)),
    list(family = "PopDensity", stem = "popdensity", mode = "year", years = 2010:2023),
    list(family = "CCNL",                    stem = "ccnl_2013",                    mode = "single", year = 2013L),
    list(family = "Elevation",               stem = "elevation_2000",               mode = "single", year = 2000L),
    list(family = "GHSL_SMOD",               stem = "ghsl_smod_2015",               mode = "single", year = 2015L),
    list(family = "GlobalHumanModification", stem = "globalhumanmodification_2016", mode = "single", year = 2016L),
    # WSF2015 is a mask: built pixels hold 255 and everything else is MASKED, so
    # a server-side mean over unmasked pixels is a constant 255. The legacy
    # raster carried explicit zeros, making the zonal mean a built-up density.
    # unmask(0) restores that.
    list(family = "WSF", stem = "wsf_2015", mode = "single", year = 2015L,
         ee_post = "unmask0"),
    list(family = "GPW_Demographic", stem = "gpw_demographic_2010", mode = "collapse", year = 2010L)
  )
  for (fam in names(SOIL_PARITY_FAMILIES)) {
    means_only <- fam %in% SOIL_PARITY_MEANS_ONLY
    spec[[length(spec) + 1L]] <- list(
      family = fam, stem = SOIL_PARITY_FAMILIES[[fam]],
      mode = if (means_only) "bands_only" else "bands",
      year = 2010L,
      bands = if (means_only) grep("^mean_", SOIL_PARITY_BANDS, value = TRUE)
              else SOIL_PARITY_BANDS,
      # iSDAsoil stores ln(x + 1) * 10 as an integer; the legacy exports were
      # back-transformed before zonal averaging. Do the same server-side so the
      # mean is over back-transformed pixels, not the mean of the log values.
      ee_post = "isda_expm1"
    )
  }
  spec
}

#' Every column name gee_legacy_parity_spec() is expected to produce.
gee_legacy_parity_colnames <- function(spec = gee_legacy_parity_spec()) {
  unlist(lapply(spec, function(s) {
    switch(s$mode,
      year       = paste0("gee_", s$stem, "_", s$years),
      single     = paste0("gee_", s$stem),
      collapse   = paste0("gee_", s$stem, c("_annual_mean", "_annual_sd")),
      summaries  = paste0("gee_", s$stem, .parity_summary_suffixes()),
      bands      = c(paste0("gee_", s$stem, "_", s$bands),
                     paste0("gee_", s$stem, .parity_summary_suffixes())),
      bands_only = paste0("gee_", s$stem, "_", s$bands)
    )
  }), use.names = FALSE)
}

#' Append the five across-band summaries exactly as .append_gee_zonal_cols() does.
#'
#' NOTE (faithful reproduction of a questionable legacy step): the summaries are
#' taken across BANDS even when the bands are distinct quantities -- for soil
#' that means averaging depth means together with their standard deviations. The
#' resulting columns are not physically meaningful, but they are part of the
#' shared cross-country vocabulary and are fed to the models for the other four
#' countries, so a new country must reproduce them to be poolable.
.parity_band_summaries <- function(mat) {
  rmin <- suppressWarnings(apply(mat, 1, min, na.rm = TRUE))
  rmax <- suppressWarnings(apply(mat, 1, max, na.rm = TRUE))
  data.frame(
    annual_mean  = rowMeans(mat, na.rm = TRUE),
    annual_sd    = apply(mat, 1, stats::sd, na.rm = TRUE),
    annual_min   = rmin,
    annual_max   = rmax,
    annual_range = rmax - rmin
  )
}

#' Look up one manifest entry by family name, applying any parity overrides.
.parity_layer <- function(family, override = NULL) {
  hit <- Filter(function(l) identical(l$family, family), GEE_LAYER_MANIFEST)
  if (!length(hit)) stop("[parity] family not in GEE_LAYER_MANIFEST: ", family)
  layer <- hit[[1]]
  # NULL in an override means "remove this field" (e.g. drop scale_factor), so
  # assign through a name loop rather than modifyList(), which drops NULLs.
  for (k in names(override)) layer[[k]] <- override[[k]]
  layer
}

#' Server-side band transforms needed to match the legacy raster values.
.parity_ee_post <- function(img, what) {
  if (is.null(what)) return(img)
  switch(what,
    # iSDAsoil: stored as ln(x + 1) * 10.
    isda_expm1 = img$divide(10)$exp()$subtract(1),
    # Mask layers whose legacy raster carried explicit zeros outside the mask.
    unmask0    = img$unmask(0),
    stop("[parity] unknown ee_post: ", what))
}

#' Extract legacy-named admin-2 GEE covariates for one country.
#'
#' @param gadm_code ISO3 code (e.g. "TZA"); admin-2 polygons come from the
#'   shared GADM cache, matching what extract_gee_admin2() uses.
#' @param polys optional sf with Admin1/Admin2 to use instead of the GADM
#'   download (used by the validation run against a legacy country).
#' @param spec extraction spec (default gee_legacy_parity_spec()).
#' @param cache_dir per-(family, year) cache; a re-run only pulls what is
#'   missing or previously failed.
#' @return data.frame with Admin1, Admin2 and the legacy gee_* columns.
extract_legacy_parity_admin2 <- function(gadm_code, polys = NULL,
                                         spec = gee_legacy_parity_spec(),
                                         verbose = TRUE, simplify_deg = 0.01,
                                         batch = 40, scale_floor = 250,
                                         cache_dir = NULL, force = FALSE) {
  ee <- rgee::ee

  if (is.null(polys)) {
    polys <- load_gadm_cached(gadm_code, level = 2) |>
      sf::st_as_sf() |>
      dplyr::transmute(Admin1 = trimws(NAME_1), Admin2 = trimws(NAME_2))
  }
  if (!is.null(simplify_deg) && simplify_deg > 0) {
    polys <- suppressWarnings(
      sf::st_make_valid(sf::st_simplify(polys, dTolerance = simplify_deg,
                                        preserveTopology = TRUE)))
  }

  out <- sf::st_drop_geometry(polys)[, c("Admin1", "Admin2"), drop = FALSE]
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # One extraction unit = one (family, year) pair, so the cache survives a
  # partial run of the multi-year families.
  units <- unlist(lapply(spec, function(s) {
    yrs <- if (identical(s$mode, "year")) s$years else s$year
    lapply(yrs, function(y) c(s[setdiff(names(s), "years")], list(year = y)))
  }), recursive = FALSE)

  for (u in units) {
    tag <- sprintf("%s_%d", gsub("[^A-Za-z0-9_]", "_", u$family), u$year)
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(tag, ".rds")) else NULL

    if (!force && !is.null(cache_file) && file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (!is.null(cached) && nrow(cached) == nrow(out) &&
          identical(as.character(cached$Admin2), as.character(out$Admin2))) {
        out <- cbind(out, cached[, setdiff(names(cached), "Admin2"), drop = FALSE])
        if (verbose) message(sprintf("  cached: %-28s (%d col[s])", tag, ncol(cached) - 1L))
        next
      }
    }

    layer <- .parity_layer(u$family, u$override)
    if (!is.null(u$bands)) layer$bands <- u$bands
    # ee_image_for_layer() still applies the manifest's availability clamp, so a
    # requested year outside the asset's coverage falls back to the nearest one.
    img <- ee_image_for_layer(layer, u$year)
    if (is.null(img)) { if (verbose) message("  [", tag, "] image build failed"); next }
    img <- tryCatch(.parity_ee_post(img, u$ee_post),
                    error = function(e) { message("  [", tag, "] ee_post failed: ",
                                                  conditionMessage(e)); NULL })
    if (is.null(img)) next

    vals <- .ee_extract_batched(img, polys, ee$Reducer$mean(),
                                max(layer$scale, scale_floor), batch, verbose,
                                label = tag)
    if (is.null(vals) || nrow(vals) != nrow(polys)) {
      if (verbose) message("  [", tag, "] skipped (extract failed)")
      next
    }
    val_cols <- setdiff(names(vals), c("Admin1", "Admin2"))
    if (!length(val_cols)) { if (verbose) message("  [", tag, "] no value columns"); next }
    mat <- as.matrix(vals[, val_cols, drop = FALSE])
    storage.mode(mat) <- "double"

    # Single-value modes must use the BASE band only. ee_image_for_layer() can
    # append a `<band>_yoy` change band, and averaging it into the base value is
    # what previously made TRMM read at half the legacy scale.
    base_col <- val_cols[!grepl("_yoy$", val_cols)][1]
    if (is.na(base_col)) base_col <- val_cols[1]

    # EE returns one column per requested band, sometimes reordered or suffixed.
    match_bands <- function(bands) vapply(bands, function(b) {
      hit <- grep(paste0("(^|_)", b, "$"), val_cols)
      if (length(hit)) val_cols[hit[1]] else NA_character_
    }, character(1))

    add <- switch(u$mode,
      year      = stats::setNames(data.frame(mat[, base_col]),
                                  paste0("gee_", u$stem, "_", u$year)),
      single    = stats::setNames(data.frame(mat[, base_col]),
                                  paste0("gee_", u$stem)),
      collapse  = stats::setNames(
                    .parity_band_summaries(mat)[, c("annual_mean", "annual_sd")],
                    paste0("gee_", u$stem, c("_annual_mean", "_annual_sd"))),
      summaries = stats::setNames(.parity_band_summaries(mat),
                                  paste0("gee_", u$stem, .parity_summary_suffixes())),
      bands     = {
        band_cols <- match_bands(u$bands)
        if (anyNA(band_cols)) {
          if (verbose) message("  [", tag, "] missing band(s): ",
                               paste(u$bands[is.na(band_cols)], collapse = ", "))
          NULL
        } else {
          bm <- mat[, band_cols, drop = FALSE]
          cbind(
            stats::setNames(as.data.frame(bm), paste0("gee_", u$stem, "_", u$bands)),
            stats::setNames(.parity_band_summaries(bm),
                            paste0("gee_", u$stem, .parity_summary_suffixes()))
          )
        }
      },
      bands_only = {
        band_cols <- match_bands(u$bands)
        if (anyNA(band_cols)) {
          if (verbose) message("  [", tag, "] missing band(s): ",
                               paste(u$bands[is.na(band_cols)], collapse = ", "))
          NULL
        } else {
          stats::setNames(as.data.frame(mat[, band_cols, drop = FALSE]),
                          paste0("gee_", u$stem, "_", u$bands))
        }
      }
    )
    if (is.null(add)) next

    # min()/max() on all-NA rows return +/-Inf; normalise those to NA.
    add[] <- lapply(add, function(x) { x[!is.finite(x)] <- NA_real_; x })

    if (!is.null(cache_file))
      saveRDS(cbind(out[, "Admin2", drop = FALSE], add), cache_file)
    out <- cbind(out, add)
    if (verbose) message(sprintf("  ok: %-28s (%d col[s])", tag, ncol(add)))
  }
  out
}
