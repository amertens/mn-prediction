# =============================================================================
# scripts/extract_civ_gee_admin2.R
#
# Re-extract Côte d'Ivoire Admin-2 GEE features using the *same* extraction
# logic the targets pipeline uses for training countries — i.e.
# .append_gee_zonal_cols() in R/admin2_analysis.R.
#
# This produces a gee_admin2 data.frame with column names AND value scales
# identical to gee_admin2_<country> targets in _targets_full/. The earlier
# scripts/civ_external_validation.R bypassed this helper and computed its own
# zonal means, producing features 4×–290× off scale from training and breaking
# the spatial+soil GAM.
#
# The pipeline's own oos_gee_civ target uses extract_gee_for_country() in
# R/oos_prediction.R, which collapses multi-band rasters to band 1 — so it
# produces gee_soilaluminium but NOT gee_soilaluminium_stdev_0_20 (one of the
# eight spatial+soil predictors). That's why none of the existing OOS paths
# can run the spatial+soil winning model on CIV. This script closes the gap.
#
# Inputs:
#   - data/Cote_dIvoire_GEE_rasters/ (83 .tif files, already on disk)
#   - data/gadm/ (GADM polygons; downloaded on first run)
#
# Output:
#   - results/transportability/civ_gee_admin2.rds
#       Schema: Admin1, Admin2, gee_*  (matches training-country gee_admin2_*)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(terra)
  library(exactextractr)
  library(geodata)
  library(dplyr)
})

# Source the pipeline's extraction helpers
source(here::here("R", "admin2_analysis.R"))

# Construct a CIV "country config" matching the cc shape that extract_gee_admin2
# expects (uses cc$gadm_code, cc$raster_dir, cc$country).
cc_civ <- list(
  country    = "Cote dIvoire",
  gadm_code  = "CIV",
  raster_dir = here::here("data", "Cote_dIvoire_GEE_rasters")
)

stopifnot(dir.exists(cc_civ$raster_dir))

cat("[civ-extract] Starting CIV admin-2 GEE extraction\n")
cat(sprintf("[civ-extract] gadm_code  = %s\n", cc_civ$gadm_code))
cat(sprintf("[civ-extract] raster_dir = %s\n", cc_civ$raster_dir))
cat(sprintf("[civ-extract] n .tif     = %d\n",
            length(list.files(cc_civ$raster_dir, pattern = "\\.tif$"))))

t0 <- proc.time()
civ_gee_admin2 <- extract_gee_admin2(cc_civ)
elapsed <- (proc.time() - t0)[["elapsed"]]

cat(sprintf("[civ-extract] Done in %.1f min — %d polygons × %d cols\n",
            elapsed / 60, nrow(civ_gee_admin2), ncol(civ_gee_admin2)))

out_dir <- here::here("results", "transportability")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "civ_gee_admin2.rds")
saveRDS(civ_gee_admin2, out_file)
cat(sprintf("[civ-extract] Saved: %s\n", out_file))

# Quick sanity check — confirm all 8 spatial+soil default predictors are present.
spatial_plus_soil_vars <- c(
  "gee_soilaluminium_annual_min",
  "gee_soilaluminium_stdev_20_50",
  "gee_soilcalcium_stdev_0_20",
  "gee_soilaluminium_stdev_0_20",
  "gee_soilmagnesium_stdev_0_20",
  "gee_soilcalcium_stdev_20_50",
  "gee_soiltotalcarbon_mean_0_20",
  "gee_soilzinc_mean_20_50"
)
present <- spatial_plus_soil_vars %in% colnames(civ_gee_admin2)
cat("\n[civ-extract] Spatial+soil predictor coverage:\n")
for (i in seq_along(spatial_plus_soil_vars)) {
  cat(sprintf("  %s %s\n",
              if (present[i]) "OK " else "!! ",
              spatial_plus_soil_vars[i]))
}
cat(sprintf("[civ-extract] %d / %d spatial+soil predictors present\n",
            sum(present), length(present)))
