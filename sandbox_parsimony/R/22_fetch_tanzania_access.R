# =============================================================================
# sandbox_parsimony/R/22_fetch_tanzania_access.R
#
# Tanzania has no data/Tanzania_GEE_rasters/, so it lacks the single covariate
# with the strongest prior claim on dietary diversity: TRAVEL TIME TO CITIES
# (Weiss et al. 2018, distributed by the Malaria Atlas Project). The other four
# countries already have it on disk as Accessibility_<Country>_2019.tif.
#
# This fetches the MAP accessibility surface for Tanzania's bounding box only
# (not the global raster) and aggregates it to Admin-2 polygons.
#
# SOURCE:  Malaria Atlas Project, via the malariaAtlas R package
#          "Accessibility__201501_Global_Travel_Time_to_Cities" and
#          "Accessibility__202001_Global_Motorized_Travel_Time_to_Healthcare"
# LICENCE: CC BY 4.0
# SIZE:    ~10-40 MB per surface for the Tanzania bounding box
#
# Writes sandbox_parsimony/out/tanzania_access_admin2.csv, which
# 21_extend_covariates.R picks up automatically if present. Safe to skip: the
# rest of the sandbox runs without it.
# =============================================================================
suppressPackageStartupMessages({library(terra); library(sf); library(dplyr)})

OUT <- "sandbox_parsimony/out/tanzania_access_admin2.csv"
CACHE <- "sandbox_parsimony/out/_tz_access_cache"
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

ac <- tryCatch(readRDS("_targets_full/objects/area_covariates_tanzania"),
               error = function(e) NULL)
if (is.null(ac) || is.null(ac$polygons)) {
  message("No Tanzania polygons available; nothing to do.")
} else {
  pg <- ac$polygons
  bb <- sf::st_bbox(sf::st_transform(pg, 4326))
  message(sprintf("Tanzania bbox: [%.2f, %.2f] x [%.2f, %.2f]",
                  bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))

  SURFACES <- c(
    tt_city   = "Accessibility__201501_Global_Travel_Time_to_Cities",
    tt_health = "Accessibility__202001_Global_Motorized_Travel_Time_to_Healthcare")

  get_surface <- function(key, dataset_id) {
    f <- file.path(CACHE, paste0(key, ".tif"))
    if (file.exists(f)) {
      message("  cached: ", key)
      return(tryCatch(terra::rast(f), error = function(e) NULL))
    }
    if (!requireNamespace("malariaAtlas", quietly = TRUE)) {
      message("  malariaAtlas not installed; skipping ", key); return(NULL)
    }
    message("  downloading ", dataset_id, " for the Tanzania bounding box ...")
    r <- tryCatch(
      malariaAtlas::getRaster(dataset_id = dataset_id,
                              shp = sf::as_Spatial(sf::st_as_sfc(bb))),
      error = function(e) { message("    getRaster failed: ", conditionMessage(e)); NULL })
    if (is.null(r)) return(NULL)
    r <- tryCatch(terra::rast(r), error = function(e) r)
    tryCatch({ terra::writeRaster(r, f, overwrite = TRUE)
               message(sprintf("    saved %s (%.1f MB)", basename(f),
                               file.size(f) / 1e6)) },
             error = function(e) message("    could not cache: ", conditionMessage(e)))
    r
  }

  res <- data.frame(Admin2 = as.character(pg$Admin2), stringsAsFactors = FALSE)
  got <- character(0)
  for (k in names(SURFACES)) {
    r <- get_surface(k, SURFACES[[k]])
    if (is.null(r)) next
    v <- tryCatch({
      pv <- terra::vect(sf::st_transform(pg, terra::crs(r)))
      # travel time is heavily right-skewed and a district mean is dominated by
      # its remotest pixel, so take the median as well as the mean
      terra::extract(r, pv, fun = function(x, ...)
        c(mean(x, na.rm = TRUE), stats::median(x, na.rm = TRUE)))
    }, error = function(e) { message("    extract failed: ", conditionMessage(e)); NULL })
    if (is.null(v)) next
    vv <- as.data.frame(v)
    res[[paste0("gee_", k, "_mean")]]   <- suppressWarnings(as.numeric(vv[[2]]))
    res[[paste0("gee_", k, "_median")]] <- suppressWarnings(as.numeric(vv[[min(3, ncol(vv))]]))
    got <- c(got, k)
  }

  if (!length(got)) {
    message("\nNo surfaces retrieved. Tanzania keeps its reduced covariate set;\n",
            "21_extend_covariates.R will run without them.")
  } else {
    res <- res[!duplicated(res$Admin2), , drop = FALSE]
    write.csv(res, OUT, row.names = FALSE)
    message(sprintf("\nSaved %s (%d districts, surfaces: %s)",
                    OUT, nrow(res), paste(got, collapse = ", ")))
    print(summary(res[, -1, drop = FALSE]))
  }
}
