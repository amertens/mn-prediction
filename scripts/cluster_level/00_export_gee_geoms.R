# Export simplified Admin-2 polygons and cluster points as GeoJSON for the Earth Engine reducers.
suppressPackageStartupMessages({library(sf); library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
cmap <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
polys <- bind_rows(lapply(names(cmap), function(k) {
  b <- BND[[k]]; if (is.null(b)) return(NULL); b <- st_as_sf(b)
  st_sf(country = cmap[[k]], Admin1 = as.character(b$Admin1), Admin2 = as.character(b$Admin2), geometry = st_geometry(b))
}))
polys <- st_make_valid(st_transform(polys, 4326))
# simplify in a metric CRS (500 m tolerance), then back to lon/lat; keep the original where simplification empties a polygon
simp <- st_transform(st_simplify(st_transform(polys, 3857), dTolerance = 500, preserveTopology = TRUE), 4326)
simp <- st_make_valid(simp); bad <- st_is_empty(simp) | !st_is_valid(simp)
if (any(bad)) st_geometry(simp)[bad] <- st_geometry(polys)[bad]
dir.create("data/external_cache/gee_geoms", showWarnings = FALSE, recursive = TRUE)
f1 <- "data/external_cache/gee_geoms/admin2_simplified.geojson"; unlink(f1); st_write(simp, f1, driver = "GeoJSON", quiet = TRUE)
cat("admin2 polygons:", nrow(simp), "| file MB:", round(file.size(f1) / 1e6, 2), "\n"); print(table(simp$country))
TC <- read.csv("results/tables/cluster_level/targets_cluster.csv") |> distinct(country, cluster, lat, lon) |> filter(is.finite(lat), is.finite(lon))
RURAL_KM <- as.numeric(Sys.getenv("CL_RURAL_KM", "5")); URBAN_KM <- as.numeric(Sys.getenv("CL_URBAN_KM", "2")); OT <- Sys.getenv("CL_OUT_TAG", "")
PCU <- read.csv("data/covariates/cluster/predictors_cluster.csv", check.names = FALSE)[, c("country", "cluster", "urban")]
TC$cluster <- as.character(TC$cluster); PCU$cluster <- as.character(PCU$cluster); TC <- left_join(TC, PCU, by = c("country", "cluster"))
TC$radius_km <- ifelse(TC$urban %in% c(1, TRUE), URBAN_KM, RURAL_KM)
pts <- st_as_sf(TC, coords = c("lon", "lat"), crs = 4326)
f2 <- paste0("data/external_cache/gee_geoms/clusters", OT, ".geojson"); unlink(f2); st_write(pts, f2, driver = "GeoJSON", quiet = TRUE)
cat("clusters:", nrow(pts), "\n"); print(table(pts$country, pts$radius_km))
cat("DONE\n")
