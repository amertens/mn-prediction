# =============================================================================
# sandbox_parsimony/R/merge_tanzania_chunks.R
#
# Earth Engine caps a direct download at 50 MB, so FLDAS (143 MB) and the
# 12-month LST night stack (90 MB) had to be fetched in band chunks
# (fetch_tanzania_rasters3.py). This stitches them back into single files whose
# NAME and BAND SET match the other countries' exports.
#
# WHY THE BAND NAMES ARE SET EXPLICITLY. Band names do not survive
# terra::writeRaster to a plain GeoTIFF, and .append_gee_zonal_cols() derives
# one column per band FROM THE BAND NAME. Left unset, Tanzania's FLDAS columns
# come out as gee_fldas_..._chunk_0_1 and match nothing in the cross-country
# name intersection -- the same silent-vocabulary-mismatch failure that made
# Tanzania contribute nothing in the first place.
#
# Idempotent: if the chunks are already gone it repairs the merged files in
# place, so it is safe to re-run.
# =============================================================================
suppressPackageStartupMessages(library(terra))
`%||%` <- function(a, b) if (is.null(a)) b else a

DIR <- "data/Tanzania_GEE_rasters"

# Exact band list from the source collection, not a guess:
#   ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001").first().bandNames()
FLDAS_BANDS <- c("Evap_tavg", "LWdown_f_tavg", "Lwnet_tavg", "Psurf_f_tavg",
                 "Qair_f_tavg", "Qg_tavg", "Qh_tavg", "Qle_tavg", "Qs_tavg",
                 "Qsb_tavg", "RadT_tavg", "Rainf_f_tavg", "SnowCover_inst",
                 "SnowDepth_inst", "Snowf_tavg", "SoilMoi00_10cm_tavg",
                 "SoilMoi10_40cm_tavg", "SoilMoi40_100cm_tavg",
                 "SoilMoi100_200cm_tavg", "SoilTemp00_10cm_tavg",
                 "SoilTemp10_40cm_tavg", "SoilTemp40_100cm_tavg",
                 "SoilTemp100_200cm_tavg", "SWdown_f_tavg", "SWE_inst",
                 "Swnet_tavg", "Tair_f_tavg", "Wind_f_tavg")

GROUPS <- list(
  list(out = "FLDAS_2010_Annual_Mean_Tanzania.tif",
       pat = "^_chunk_FLDAS_2010_Annual_Mean_Tanzania_[0-9]+\\.tif$",
       bands = FLDAS_BANDS),
  list(out = "LST_Night_2010_Monthly_Tanzania.tif",
       pat = "^_chunk_LST_Night_2010_Monthly_Tanzania_[0-9]+\\.tif$",
       bands = sprintf("2010_%02d_Mean", 1:12))
)

for (g in GROUPS) {
  dest <- file.path(DIR, g$out)
  fs <- list.files(DIR, pattern = g$pat, full.names = TRUE)
  if (length(fs)) {
    # numeric order, so chunk 10 does not sort before chunk 2
    idx <- as.integer(sub(".*_([0-9]+)\\.tif$", "\\1", basename(fs)))
    r <- do.call(c, lapply(fs[order(idx)], terra::rast))
  } else if (file.exists(dest)) {
    r <- terra::rast(dest)          # repair a previously merged file in place
  } else {
    message("nothing to do for ", g$out); next
  }
  if (terra::nlyr(r) == length(g$bands)) {
    names(r) <- g$bands
  } else {
    warning(sprintf("%s: %d bands but %d names; band names NOT set",
                    g$out, terra::nlyr(r), length(g$bands)))
  }
  # terra refuses to write onto the file it is reading, so an in-place repair
  # goes via a temporary file.
  tmp <- file.path(DIR, paste0("_tmp_", basename(dest)))
  terra::writeRaster(r, tmp, overwrite = TRUE)
  rm(r); gc(verbose = FALSE)
  file.rename(tmp, dest)
  r <- terra::rast(dest)
  if (length(fs)) file.remove(fs)
  message(sprintf("%-42s %2d bands: %s", g$out, terra::nlyr(r),
                  paste(utils::head(names(r), 4), collapse = ", ")))
}

# Multi-band files written by the single-shot fetcher also lost their band
# names to writeRaster. Restore the ones that CAN line up with the other
# countries' vocabulary.
#
# ESA WorldCereal is deliberately absent from this list: its band names carry
# the AEZ id of the tile (Ghana's read "32121_TC-MAIZE-MAIN_..."), so they are
# country-specific by construction and will never match across countries no
# matter what is done here. Its columns stay Tanzania-only.
SINGLE_SHOT <- list(
  list(f = "LandCoverType_Tanzania_2010.tif",           # Ghana: 2016_01_01_LC_Type1..5
       bands = sprintf("2010_01_01_LC_Type%d", 1:5)),
  list(f = "Productivity_Tanzania_2010.tif",            # Ghana: 01_01_Gpp / 01_01_Npp
       bands = c("01_01_Gpp", "01_01_Npp"))
)
for (g in SINGLE_SHOT) {
  p <- file.path(DIR, g$f)
  if (!file.exists(p)) next
  r <- terra::rast(p)
  if (terra::nlyr(r) != length(g$bands)) {
    warning(sprintf("%s: %d bands but %d names; skipped", g$f,
                    terra::nlyr(r), length(g$bands)))
    next
  }
  if (identical(names(r), g$bands)) next        # already correct; nothing to do
  names(r) <- g$bands
  tmp <- file.path(DIR, paste0("_tmp_", g$f))
  terra::writeRaster(r, tmp, overwrite = TRUE)
  rm(r); gc(verbose = FALSE)
  file.rename(tmp, p)
  message(sprintf("renamed bands in %-38s -> %s", g$f,
                  paste(utils::head(g$bands, 3), collapse = ", ")))
}

# Pass 1 of the fetcher wrote two files whose names do not follow the other
# countries' convention, and one duplicates a pass-2 file exactly. Remove them
# so the same measurement cannot enter the covariate table twice under two
# different derived column names.
STALE <- c("LST_Night_Annual_Mean_Tanzania_2010.tif",  # dup of LST_Night_2010_Annual_Mean_*
           "LST_Day_Monthly_Tanzania_2010.tif")        # single band, named "Monthly"
for (f in STALE) {
  p <- file.path(DIR, f)
  if (file.exists(p)) { file.remove(p); message("removed stale ", f) }
}

# Sanity value per file, so an all-NA or nodata-contaminated layer cannot pass
# unnoticed -- an earlier pass shipped a GHSPOP layer whose mean was -10.96.
message("\nfinal contents:")
for (f in sort(list.files(DIR, pattern = "\\.tif$", full.names = TRUE))) {
  r <- tryCatch(terra::rast(f), error = function(e) NULL)
  if (is.null(r)) { message(sprintf("  %-44s UNREADABLE", basename(f))); next }
  v <- tryCatch(terra::global(r[[1]], "mean", na.rm = TRUE)[1, 1],
                error = function(e) NA_real_)
  message(sprintf("  %-44s %3d band(s)  b1=%-22s mean = %s", basename(f),
                  terra::nlyr(r), utils::head(names(r), 1),
                  if (is.finite(v)) sprintf("%.3f", v) else "ALL NA"))
}
