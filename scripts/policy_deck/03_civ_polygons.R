# =============================================================================
# scripts/policy_deck/03_civ_polygons.R   [CV-01, 2026-09-17]
#
# The Cote d'Ivoire Admin-2 polygons every CIV extraction joins on, keyed
# exactly as results/transportability/civ_canonical_admin2_full.rds (GADM 4.1
# level 2: 33 regions, NAME_1 -> Admin1, NAME_2 -> Admin2), written once as
#   data/external_cache/gee_geoms/civ_admin2.gpkg          (full geometry)
#   data/external_cache/gee_geoms/civ_admin2_simplified.geojson   (for Earth Engine)
# so that scripts 03b (raster domains) and 03c (climatology) use one spine.
# =============================================================================
suppressPackageStartupMessages({library(sf); library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/data_prep.R")   # load_gadm_cached()
CIV <- readRDS("results/transportability/civ_canonical_admin2_full.rds")
g <- sf::st_as_sf(load_gadm_cached("CIV", level = 2)); g <- sf::st_make_valid(sf::st_transform(g, 4326))
g <- g |> transmute(country = "CoteDIvoire", Admin1 = as.character(NAME_1), Admin2 = as.character(NAME_2))
miss <- setdiff(paste(CIV$Admin1, CIV$Admin2), paste(g$Admin1, g$Admin2))
extra <- setdiff(paste(g$Admin1, g$Admin2), paste(CIV$Admin1, CIV$Admin2))
cat(sprintf("GADM CIV level 2: %d polygons; canonical rds: %d rows; unmatched rds %d, unmatched gadm %d\n", nrow(g), nrow(CIV), length(miss), length(extra)))
if (length(miss)) { print(miss); print(extra) }
stopifnot(length(miss) == 0)
g <- g[match(paste(CIV$Admin1, CIV$Admin2), paste(g$Admin1, g$Admin2)), ]
dir.create("data/external_cache/gee_geoms", showWarnings = FALSE, recursive = TRUE)
sf::st_write(g, "data/external_cache/gee_geoms/civ_admin2.gpkg", delete_dsn = TRUE, quiet = TRUE)
gs <- sf::st_simplify(g, dTolerance = 0.01); sf::st_write(gs, "data/external_cache/gee_geoms/civ_admin2_simplified.geojson", delete_dsn = TRUE, quiet = TRUE)
cat("written civ_admin2.gpkg and civ_admin2_simplified.geojson\n")
