# =============================================================================
# src/Tanzania/TZ_MAP_download.R
#
# Download Malaria Atlas Project rasters clipped to Tanzania (ISO=TZA), year
# 2010, into data/Malaria Atlas/Tanzania/. Adapted from
# src/SierraLeone/SL_MAP_download.R. The merge (2_GW_Tanzania_data_merge.R,
# RUN_MAP=true) then extracts these at the Tanzania cluster points -> MAP_ cols.
#
#   Rscript src/Tanzania/TZ_MAP_download.R
# =============================================================================
rm(list = ls())
suppressPackageStartupMessages({ library(here); library(terra); library(malariaAtlas) })

out_dir <- here("data", "Malaria Atlas", "Tanzania")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# The malariaAtlas catalogue has grown well beyond the layer set the pipeline
# was built on (it now also serves ~100 "Explorer__" vector-occurrence layers
# and "Accessibility__" surfaces that no other country in this project has).
# The merge (2_GW_Tanzania_data_merge.R) only maps the three canonical prefixes
# below to MAP_ columns, so restrict the download to exactly that set. This keeps
# Tanzania byte-consistent with the other four countries (each 43 rasters) and
# makes re-runs reproducible regardless of future catalogue additions.
CANONICAL_PREFIXES <- "^(Malaria|Interventions|Blood_Disorders)__"
map_rasters <- listRaster()
map_rasters <- map_rasters[grepl(CANONICAL_PREFIXES, map_rasters$dataset_id), , drop = FALSE]
cat(sprintf("Canonical MAP layers to ensure present: %d\n", nrow(map_rasters)))
tza_shp <- getShp(ISO = "TZA", admin_level = "admin0")

ok <- 0L
for (i in seq_len(nrow(map_rasters))) {
  did <- map_rasters$dataset_id[i]
  dst <- file.path(out_dir, paste0(did, ".tif"))
  if (file.exists(dst)) { ok <- ok + 1L; next }        # resume: skip already-downloaded
  r <- try(getRaster(dataset_id = did, shp = tza_shp, year = 2010), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) {
    if (!inherits(r, "SpatRaster")) r <- try(terra::rast(r), silent = TRUE)
    if (inherits(r, "SpatRaster")) {
      try(terra::writeRaster(r, dst, overwrite = TRUE), silent = TRUE)
      if (file.exists(dst)) { ok <- ok + 1L; cat("  ok:", did, "\n") }
    }
  }
}
cat(sprintf("\nDownloaded/present %d of %d MAP rasters -> %s\n", ok, nrow(map_rasters), out_dir))
