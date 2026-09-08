# =============================================================================
# scripts/protocol_v2/52_map_survey_year.R   [MP-01]
#
# MALARIA ATLAS LAYERS AT THE SURVEY YEAR, ONE RELEASE PER PRODUCT
#
# The existing map_* columns were extracted from rasters requested at the
# year embedded in the dataset id, which is the RELEASE date (202206, 202406,
# 202508), 4-12 years after the surveys, and three releases of the same
# product sit in the vocabulary as separate columns (AU-01 finding 4). Here
# each of the nine time-varying products is fetched from its latest release
# for each country at that country's survey year (clamped to the product's
# year range), cached as one GeoTIFF per country, product and year, and
# aggregated to the spine polygons by the same area-weighted mean the
# original block used. The static blood-disorder surfaces (2012 release) are
# unchanged.
#
# Every WCS request takes about 2.5 minutes whatever its extent, and a
# single request over the four-country bounding box comes back without cell
# values, so the 36 country x product requests are made per country and the
# four countries run in parallel:
#
#   MAP_COUNTRY=Gambia Rscript scripts/protocol_v2/52_map_survey_year.R   (download only, x4)
#   Rscript scripts/protocol_v2/52_map_survey_year.R                      (all countries: fills any gap, then aggregates)
# -> data/external_cache/malaria_atlas_sy/<country>/<product>_<year>.tif
#    data/covariates/harmonized/predictors_admin2_map_sy.csv (+ _metadata)
# Script 08 then drops the release-year map_malaria* / map_interventions*
# columns and appends this block.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf); library(terra); library(malariaAtlas)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"; CACHE <- "data/external_cache/malaria_atlas_sy"
source("R/survey_years.R"); SURVEY_YEAR <- survey_years()   # single source: metadata/survey_years.csv (Gambia 2018, Ghana 2017, Malawi 2016, Sierra Leone 2013)
LC <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
SHARD <- Sys.getenv("MAP_COUNTRY", ""); SHARD <- if (nzchar(SHARD)) trimws(strsplit(SHARD, ",")[[1]]) else names(SURVEY_YEAR)
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
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
polys <- function(cn) { b <- st_as_sf(BND[[LC[[cn]]]]); b <- st_make_valid(st_transform(b, 4326))
  st_sf(country = cn, Admin1 = as.character(b$Admin1), Admin2 = as.character(b$Admin2), geometry = st_geometry(b)) }
P <- lapply(names(SURVEY_YEAR), polys); names(P) <- names(SURVEY_YEAR)
avail <- tryCatch(listRaster(printed = FALSE), error = function(e) NULL)
yr_range <- function(id) { r <- avail[avail$dataset_id == id, ]; if (!nrow(r)) c(NA, NA) else c(r$min_raster_year[1], r$max_raster_year[1]) }
use_year <- function(id, yr) { rg <- yr_range(id); if (all(is.finite(rg))) min(max(yr, rg[1]), rg[2]) else yr }
layer_file <- function(cn, id, yr) file.path(CACHE, cn, sprintf("%s_%d.tif", id, yr))
# NB: terra::rast(<SpatRaster>) makes an EMPTY template (no values): only convert raster-package objects
first_layer <- function(r) { x <- if (inherits(r, "SpatRasterCollection")) r[1] else if (inherits(r, "SpatRaster")) r else tryCatch(rast(r), error = function(e) NULL)
  if (is.null(x)) return(NULL); x <- x[[1]]; if (!hasValues(x)) NULL else x }

# ── downloads (one request per country x product) ────────────────────────────
for (cn in SHARD) { yr <- SURVEY_YEAR[[cn]]; dir.create(file.path(CACHE, cn), showWarnings = FALSE, recursive = TRUE)
  bbox <- st_as_sf(st_as_sfc(st_bbox(st_buffer(st_union(P[[cn]]), 0.05))))
  for (col in names(PRODUCTS)) { id <- PRODUCTS[[col]]; uy <- use_year(id, yr); f <- layer_file(cn, id, uy)
    if (file.exists(f)) { cat(sprintf("  %-12s %-32s %d cached\n", cn, col, uy)); flush.console(); next }
    t0 <- Sys.time()
    r <- tryCatch(getRaster(dataset_id = id, shp = as(bbox, "Spatial"), year = uy), error = function(e) { cat(sprintf("  %-12s %-32s FAILED: %s\n", cn, col, conditionMessage(e))); NULL })
    x <- if (is.null(r)) NULL else first_layer(r)
    if (is.null(x)) { cat(sprintf("  %-12s %-32s %d: no cell values returned\n", cn, col, uy)); flush.console(); next }
    writeRaster(x, f, overwrite = TRUE)
    cat(sprintf("  %-12s %-32s %d written in %.0f s (mean %.4f)\n", cn, col, uy, as.numeric(Sys.time() - t0, units = "secs"), global(x, "mean", na.rm = TRUE)[1, 1])); flush.console() } }
if (nzchar(Sys.getenv("MAP_COUNTRY", ""))) { cat("\nshard", paste(SHARD, collapse = ","), "DONE\n"); quit(save = "no") }

# ── aggregation (all countries) ──────────────────────────────────────────────
rows <- list(); meta <- list()
for (cn in names(SURVEY_YEAR)) { yr <- SURVEY_YEAR[[cn]]; out <- st_drop_geometry(P[[cn]])
  for (col in names(PRODUCTS)) { id <- PRODUCTS[[col]]; uy <- use_year(id, yr); f <- layer_file(cn, id, uy)
    if (!file.exists(f)) { cat(sprintf("  %-12s %-32s no layer for %d\n", cn, col, uy)); next }
    r <- rast(f); v <- tryCatch({ vp <- vect(st_transform(P[[cn]], crs(r))); terra::extract(r[[1]], vp, fun = mean, na.rm = TRUE, ID = FALSE)[, 1] }, error = function(e) NULL)
    if (is.null(v) || length(v) != nrow(out)) { cat(sprintf("  %-12s %-32s extraction failed\n", cn, col)); next }
    out[[col]] <- as.numeric(v)
    meta[[length(meta) + 1L]] <- data.frame(country = cn, column = col, dataset_id = id, year_requested = yr, year_used = uy, finite_share = round(mean(is.finite(v)), 3), stringsAsFactors = FALSE) }
  cat(sprintf("%-12s survey %d: %d columns\n", cn, yr, ncol(out) - 3)); rows[[cn]] <- out }
OUT <- bind_rows(rows); MD <- bind_rows(meta)
write.csv(OUT, file.path(HDIR, "predictors_admin2_map_sy.csv"), row.names = FALSE); write.csv(MD, file.path(HDIR, "predictors_admin2_map_sy_metadata.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %d rows x %d map_sy columns\n", nrow(OUT), ncol(OUT) - 3))
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
pairs <- c(map_sy_pf_parasite_rate = "map_malaria202406globalpfparasiterate", map_sy_pf_incidence_rate = "map_malaria202406globalpfincidencerate",
           map_sy_itn_use = "map_interventions202406africainsecticidetreatednetuse", map_sy_irs_coverage = "map_interventions202106africairscoverage")
pairs <- pairs[pairs %in% names(S)]
j <- inner_join(OUT, S[, c("country", "Admin1", "Admin2", unname(pairs))], by = c("country", "Admin1", "Admin2"))
cat("\nSpearman between the survey-year layer and the release-year column, by country:\n")
for (nm in names(pairs)) if (nm %in% names(j)) cat(sprintf("  %-26s %s\n", nm, paste(sapply(names(SURVEY_YEAR), function(cn) { k <- j$country == cn; sprintf("%s %.2f", substr(cn, 1, 6), suppressWarnings(cor(j[[nm]][k], j[[pairs[[nm]]]][k], method = "spearman", use = "complete.obs"))) }), collapse = " | ")))
cat("\nDONE\n")
