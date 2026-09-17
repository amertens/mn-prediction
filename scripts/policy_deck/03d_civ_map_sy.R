# =============================================================================
# scripts/policy_deck/03d_civ_map_sy.R   [CV-01]
#
# Malaria Atlas layers for Cote d'Ivoire at the reference year, the same nine
# products and releases as scripts/protocol_v2/52_map_survey_year.R, cached as
#   data/external_cache/malaria_atlas_sy/CoteDIvoire/<product>_<year>.tif
# Every WCS request takes ~2.5 minutes; run in the background. Script 03b
# aggregates whatever is cached. CIV has no biomarker survey, so the reference
# year is CIV_REF_YEAR (default 2016, the MICS 2016 year; see 03b).
# =============================================================================
suppressPackageStartupMessages({library(sf); library(terra); library(malariaAtlas)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
CACHE <- "data/external_cache/malaria_atlas_sy/CoteDIvoire"; dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
YR <- as.integer(Sys.getenv("CIV_REF_YEAR", "2016"))
PRODUCTS <- c(
  map_sy_pf_parasite_rate       = "Malaria__202508_Global_Pf_Parasite_Rate",
  map_sy_pf_incidence_rate      = "Malaria__202508_Global_Pf_Incidence_Rate",
  map_sy_pf_mortality_rate      = "Malaria__202508_Global_Pf_Mortality_Rate",
  map_sy_pf_reproductive_number = "Malaria__202202_Global_Pf_Reproductive_Number",
  map_sy_itn_access             = "Interventions__202508_Africa_Insecticide_Treated_Net_Access",
  map_sy_itn_use                = "Interventions__202508_Africa_Insecticide_Treated_Net_Use",
  map_sy_itn_use_rate           = "Interventions__202508_Africa_Insecticide_Treated_Net_Use_Rate",
  map_sy_irs_coverage           = "Interventions__202508_Africa_IRS_Coverage",
  map_sy_effective_treatment    = "Interventions__202508_Global_Antimalarial_Effective_Treatment")
P <- sf::st_read("data/external_cache/gee_geoms/civ_admin2.gpkg", quiet = TRUE)
bbox <- st_as_sf(st_as_sfc(st_bbox(st_buffer(st_union(P), 0.05))))
avail <- tryCatch(listRaster(printed = FALSE), error = function(e) NULL)
yr_range <- function(id) { r <- avail[avail$dataset_id == id, ]; if (!nrow(r)) c(NA, NA) else c(r$min_raster_year[1], r$max_raster_year[1]) }
use_year <- function(id, yr) { rg <- yr_range(id); if (all(is.finite(rg))) min(max(yr, rg[1]), rg[2]) else yr }
first_layer <- function(r) { x <- if (inherits(r, "SpatRasterCollection")) r[1] else if (inherits(r, "SpatRaster")) r else tryCatch(rast(r), error = function(e) NULL)
  if (is.null(x)) return(NULL); x <- x[[1]]; if (!hasValues(x)) NULL else x }
for (col in names(PRODUCTS)) { id <- PRODUCTS[[col]]; uy <- use_year(id, YR); f <- file.path(CACHE, sprintf("%s_%d.tif", id, uy))
  if (file.exists(f)) { cat(sprintf("  %-32s %d cached\n", col, uy)); next }
  t0 <- Sys.time()
  r <- tryCatch(getRaster(dataset_id = id, shp = as(bbox, "Spatial"), year = uy), error = function(e) { cat(sprintf("  %-32s FAILED: %s\n", col, conditionMessage(e))); NULL })
  x <- if (is.null(r)) NULL else first_layer(r)
  if (is.null(x)) { cat(sprintf("  %-32s %d: no cell values returned\n", col, uy)); flush.console(); next }
  writeRaster(x, f, overwrite = TRUE)
  cat(sprintf("  %-32s %d written in %.0f s (mean %.4f)\n", col, uy, as.numeric(Sys.time() - t0, units = "secs"), global(x, "mean", na.rm = TRUE)[1, 1])); flush.console() }
# static blood-disorder surfaces (2012 release), the three the shared set carries
STATIC <- c(map_blooddisorders201201africahbcallelefrequency = "Blood_Disorders__201201_Africa_HbC_Allele_Frequency",
            map_blooddisorders201201globalg6pddallelefrequency = "Blood_Disorders__201201_Global_G6PDd_Allele_Frequency",
            map_blooddisorders201201globalsicklehaemoglobinhbsallelefrequency = "Blood_Disorders__201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency")
for (col in names(STATIC)) { id <- STATIC[[col]]; f <- file.path(CACHE, paste0(id, ".tif"))
  if (file.exists(f)) { cat(sprintf("  %-32s cached\n", col)); next }
  t0 <- Sys.time()
  r <- tryCatch(getRaster(dataset_id = id, shp = as(bbox, "Spatial")), error = function(e) { cat(sprintf("  %-32s FAILED: %s\n", col, conditionMessage(e))); NULL })
  x <- if (is.null(r)) NULL else first_layer(r)
  if (is.null(x)) { cat(sprintf("  %-32s no cell values returned\n", col)); next }
  writeRaster(x, f, overwrite = TRUE); cat(sprintf("  %-32s written in %.0f s\n", col, as.numeric(Sys.time() - t0, units = "secs"))) }
cat("DONE\n")
