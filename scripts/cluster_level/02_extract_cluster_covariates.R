# =============================================================================
# scripts/cluster_level/02_extract_cluster_covariates.R   [CL-02]
#
# COVARIATES AT THE CLUSTER, FROM RASTERS ALREADY ON DISK
#
# Buffer means around each survey cluster's GPS point, following the DHS
# convention (2 km for urban clusters, 5 km for rural; CL_URBAN_KM /
# CL_RURAL_KM override) with the urban flag read from GHSL SMOD at the point
# (code >= 21: urban centre, dense or semi-dense cluster). No Earth Engine
# call: every layer comes from data/<Country>_GEE_rasters/, the external
# caches (CHIRPS, WorldPop, VIIRS night lights, Malaria Atlas), the IHME
# rasters and the GLW4 livestock surface. The recipe is the one in
# docs/findings/CLUSTER_LEVEL_DESIGN_2026-09.md: ONE spatial summary per layer
# (SD only for elevation), and for dynamic layers a climatology mean,
# seasonal amplitude and peak month, plus fieldwork-window value, anomaly and
# preceding-three-month sum (columns ending _fw, _fw_anom, _prev3), which need
# the survey's dates (FW-01) and are adjusters, not transportable predictors.
#
# Extraction uses area-weighted means (terra::extract exact = TRUE), so a
# 5 km buffer on a 10 km rainfall pixel still gets a value.
#
#   Rscript scripts/cluster_level/02_extract_cluster_covariates.R
# -> data/covariates/cluster/predictors_cluster.csv
# -> data/covariates/cluster/predictors_cluster_metadata.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(terra); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/cluster_level"; CDIR <- "data/covariates/cluster"; dir.create(CDIR, showWarnings = FALSE, recursive = TRUE)
RURAL_KM <- as.numeric(Sys.getenv("CL_RURAL_KM", "5")); URBAN_KM <- as.numeric(Sys.getenv("CL_URBAN_KM", "2"))
OT <- Sys.getenv("CL_OUT_TAG", "")   # suffix for the output files, e.g. "_r10" for the 10 km rural sensitivity
CFG <- list(
  Gambia      = list(dir = "data/Gambia_GEE_rasters",       tag = "Gambia",       iso = "GMB", year = 2018, utm = 32628),
  Ghana       = list(dir = "data/Ghana_GEE rasters",        tag = "Ghana",        iso = "GHA", year = 2017, utm = 32630),
  Malawi      = list(dir = "data/Malawi_GEE_rasters",       tag = "Malawi",       iso = "MWI", year = 2015, utm = 32736),
  # Sierra Leone's rasters are named both "Sierra_Leone" and "Sierra Leone";
  # the tag is a regex so both are found (the first run lost 44 columns to this)
  SierraLeone = list(dir = "data/Sierra_Leone_GEE_rasters", tag = "Sierra[ _]Leone", iso = "SLE", year = 2013, utm = 32629))
TC <- read.csv(file.path(OUTDIR, "targets_cluster.csv"), stringsAsFactors = FALSE)
CL <- TC |> filter(is.finite(lat), is.finite(lon)) |> group_by(country, cluster) |> summarise(lat = lat[1], lon = lon[1], month_med = month_med[1], .groups = "drop")
san <- function(x) { x <- gsub("[^A-Za-z0-9]+", "_", tolower(x)); gsub("^_|_$", "", gsub("_+", "_", x)) }
find_file <- function(dir, pattern) { f <- list.files(dir, pattern = pattern, full.names = TRUE); if (!length(f)) NA_character_ else f[1] }
pick_year <- function(dir, stem, tag, year) { fs <- list.files(dir, pattern = paste0("^", stem, "_", tag, "_[0-9]{4}[.]tif$"), full.names = TRUE)
  if (!length(fs)) return(NA_character_); yrs <- as.integer(sub(".*_([0-9]{4})[.]tif$", "\\1", fs)); fs[which.min(abs(yrs - year))] }
# Timing lesson from the first run: cropping a 30 m country raster and
# extracting with exact = TRUE took ~150 s for eight buffers, and the 10 m
# settlement footprint never finished. So: no crop for country-sized rasters
# (terra reads each polygon's window on its own); exact area weights only for
# coarse rasters (>= ~500 m), where a 5 km buffer covers few cells and the
# weights matter; cell-centre inclusion for fine rasters, where a buffer holds
# thousands of cells and the weights do not.
ex <- function(f, v, bands = NULL, fun = "mean") {
  if (is.na(f) || !file.exists(f)) return(NULL)
  t0 <- Sys.time()
  r <- tryCatch(rast(f), error = function(e) NULL); if (is.null(r)) return(NULL)
  if (!is.null(bands)) { k <- grep(bands, names(r)); if (!length(k)) return(NULL); r <- r[[k]] }
  fine <- res(r)[1] < 0.005
  if (!fine) r <- tryCatch(crop(r, ext(v) + 0.5), error = function(e) r)
  m <- tryCatch({ if (fun == "mean") terra::extract(r, v, fun = mean, na.rm = TRUE, exact = !fine, ID = FALSE)
                  else terra::extract(r, v, fun = fun, na.rm = TRUE, ID = FALSE) }, error = function(e) NULL)
  cat(sprintf("      %-44s %3d band(s) %6.1f s%s\n", substr(basename(f), 1, 44), nlyr(r), as.numeric(difftime(Sys.time(), t0, units = "secs")), if (is.null(m)) "  FAILED" else ""))
  if (is.null(m)) return(NULL); m <- as.matrix(m); colnames(m) <- names(r); m }
mstats <- function(M, prefix, months) {   # M: n x 12, January..December
  out <- data.frame(row.names = seq_len(nrow(M)))
  out[[paste0(prefix, "_mean")]] <- rowMeans(M, na.rm = TRUE)
  out[[paste0(prefix, "_amp")]]  <- apply(M, 1, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE) - min(z, na.rm = TRUE))
  out[[paste0(prefix, "_peak")]] <- apply(M, 1, function(z) if (all(is.na(z))) NA_real_ else which.max(z))
  fw <- vapply(seq_len(nrow(M)), function(i) if (is.na(months[i])) NA_real_ else M[i, months[i]], 0)
  out[[paste0(prefix, "_fw")]] <- fw; out[[paste0(prefix, "_fw_anom")]] <- fw - out[[paste0(prefix, "_mean")]]
  out[[paste0(prefix, "_prev3")]] <- vapply(seq_len(nrow(M)), function(i) if (is.na(months[i])) NA_real_ else sum(M[i, ((months[i] - 1 - 1:3) %% 12) + 1], na.rm = TRUE), 0)
  out }
ihme_file <- function(folder, pattern, year) { fs <- list.files(folder, pattern = pattern, full.names = TRUE); if (!length(fs)) return(NA_character_)
  b <- basename(fs); y <- suppressWarnings(as.integer(sub(".*_([0-9]{4})_Y[0-9]{4}M.*", "\\1", b)))
  if (all(is.na(y))) y <- vapply(regmatches(b, gregexpr("(19|20)[0-9]{2}", b)), function(z) { z <- as.integer(z); z <- z[z <= 2025]; if (length(z)) z[length(z)] else NA_integer_ }, 0L)
  if (all(is.na(y))) return(fs[1]); fs[which.min(abs(y - year))] }

META <- list(); ALL <- list()
add <- function(X, m, prefix, domain, role, source) {
  if (is.null(m)) { cat("    (no data)", prefix, "\n"); return(X) }
  # band names must be the same in every country: drop a leading year
  # ("2018_crops-coverfraction") and reduce WorldCereal tile/date names
  # ("32121_TC-MAIZE-MAIN_ACTIVECROPLAND_20210327_20211119_classification")
  # to the crop ("maize_main"), averaging tiles that share a crop
  cn0 <- colnames(m); if (!is.null(cn0)) { cn0 <- sub("^20[0-9]{2}_", "", cn0); tc <- grepl("_TC-", cn0); cn0[tc] <- tolower(gsub("-", "_", sub(".*_TC-([A-Z-]+)_.*", "\\1", cn0[tc])))
    if (any(duplicated(cn0))) { m <- sapply(unique(cn0), function(k) rowMeans(m[, cn0 == k, drop = FALSE], na.rm = TRUE)); m <- matrix(m, ncol = length(unique(cn0))); cn0 <- unique(cn0) }
    colnames(m) <- cn0 }
  nm <- if (ncol(m) == 1 && (is.null(colnames(m)) || colnames(m)[1] %in% c("", prefix))) prefix else paste0(prefix, "_", san(colnames(m)))
  for (j in seq_len(ncol(m))) X[[nm[j]]] <- as.numeric(m[, j])
  META[[length(META) + 1L]] <<- data.frame(column = nm, domain = domain, role = role, source = source, stringsAsFactors = FALSE); X }
addf <- function(X, df, domain, role, source) { for (nm in names(df)) X[[nm]] <- df[[nm]]
  META[[length(META) + 1L]] <<- data.frame(column = names(df), domain = domain, role = ifelse(grepl("_fw|_prev3", names(df)), "fieldwork", role), source = source, stringsAsFactors = FALSE); X }

for (cn in names(CFG)) { cf <- CFG[[cn]]; cl <- CL[CL$country == cn, ]; if (!nrow(cl)) next
  t0 <- Sys.time(); d <- cf$dir; tg <- cf$tag; Y <- cf$year
  pts <- st_as_sf(cl, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  smod <- find_file(d, paste0("^GHSL_SMOD_.*", tg, "[.]tif$"))
  sm <- if (!is.na(smod)) as.numeric(terra::extract(rast(smod), vect(pts), ID = FALSE)[, 1]) else rep(NA_real_, nrow(cl))
  urban <- as.integer(is.finite(sm) & sm >= 21); rad <- ifelse(urban == 1, URBAN_KM, RURAL_KM) * 1000
  buf <- st_transform(st_buffer(st_transform(pts, cf$utm), rad), 4326); v <- vect(buf)
  X <- data.frame(country = cn, cluster = cl$cluster, lat = cl$lat, lon = cl$lon, month_med = cl$month_med, smod = sm, urban = urban, buffer_km = rad / 1000, stringsAsFactors = FALSE)
  cat(sprintf("== %s: %d clusters, %d urban (2 km), %d rural (%g km)\n", cn, nrow(cl), sum(urban == 1), sum(urban == 0), RURAL_KM))
  # ── static ───────────────────────────────────────────────────────────────
  fe <- find_file(d, paste0("^Elevation_", tg)); X <- add(X, ex(fe, v), "elev", "Terrain", "static", "SRTM (GEE)")
  es <- ex(fe, v, fun = "sd"); if (!is.null(es)) { X$elev_sd <- as.numeric(es[, 1]); META[[length(META) + 1L]] <- data.frame(column = "elev_sd", domain = "Terrain", role = "static", source = "SRTM (GEE)") }
  X <- add(X, ex(find_file(d, paste0("^Accessibility_", tg)), v), "access", "Built environment", "static", "MAP accessibility (GEE)")
  X <- add(X, ex(find_file(d, paste0("^CCNL_.*", tg)), v), "ntl_ccnl", "Ruralness, population density, built environment", "slow", "CCNL (GEE)")
  X <- add(X, ex(find_file(d, paste0("^GHSBUILTS_", tg)), v), "ghs", "Ruralness, population density, built environment", "slow", "GHSL (GEE)")
  X <- add(X, ex(find_file(d, paste0("^GHSPOP_", tg)), v), "ghs", "Ruralness, population density, built environment", "slow", "GHSL (GEE)")
  X <- add(X, ex(find_file(d, paste0("^GlobalHumanModification_.*", tg)), v), "ghm", "Built environment", "slow", "gHM (GEE)")
  X <- add(X, ex(find_file(d, paste0("^WSF_", tg)), v), "wsf", "Ruralness, population density, built environment", "slow", "WSF (GEE)")
  X <- add(X, ex(pick_year(d, "LandCoverLayers", tg, Y), v, "coverfraction"), "lcover", "Land cover", "static", "Copernicus land cover (GEE)")
  X <- add(X, ex(find_file(d, paste0("^ESA_WorldCereal_.*", tg)), v), "cereal", "Agricultural production, land use", "static", "ESA WorldCereal (GEE)")
  # Sierra Leone's space-named soil files crash GDAL part-way through a read
  # (silent exit, no R error); single-band LZW copies written by
  # scratchpad/clean_sl_soil.R are used wherever a *_clean.tif exists
  for (el in c("Aluminium", "Calcium", "CEC", "Iron", "Magnesium", "Nitrogen", "Phosphorus", "Potassium", "Sulfur", "TotalCarbon", "Zinc")) {
    f <- find_file(d, paste0("^Soil", el, "_", tg, "_clean[.]tif$")); if (is.na(f)) f <- find_file(d, paste0("^Soil", el, "_", tg, "[.]tif$"))
    X <- add(X, ex(f, v, "^mean_0_20$"), paste0("soil_", tolower(el)), "Soil characteristics", "static", "iSDA soil (GEE)") }
  X <- add(X, ex("data/GEE/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif", v), "glw4_livestock", "Agricultural production, land use", "static", "GLW4")
  # ── dynamic: survey year, climatology, fieldwork window ───────────────────
  X <- add(X, ex(pick_year(d, "TerraClimate", tg, Y), v), "tclim", "Climate and weather", "dynamic", "TerraClimate (GEE)")
  # night LST: the monthly file closest to the survey year that carries all 12 months
  lfs <- list.files(d, pattern = paste0("^LST_Night_[0-9]{4}_Monthly_", tg, "[.]tif$"), full.names = TRUE)
  if (length(lfs)) { lyr <- as.integer(sub(".*LST_Night_([0-9]{4})_.*", "\\1", lfs)); nb <- vapply(lfs, function(f) length(grep("_Mean$", names(rast(f)))), 0L)
    ok <- which(nb == 12); if (length(ok)) { f <- lfs[ok][which.min(abs(lyr[ok] - Y))]; M <- ex(f, v, "_Mean$"); if (!is.null(M) && ncol(M) == 12) X <- addf(X, mstats(M, "lst_night", cl$month_med), "Climate and weather", "dynamic", "MODIS LST (GEE)") } }
  fl <- find_file(d, paste0("^FLDAS_", Y, "_Monthly_", tg)); if (is.na(fl)) fl <- find_file(d, paste0("^FLDAS_[0-9]{4}_Monthly_", tg))
  for (vv in c("Rainf_f_tavg", "Tair_f_tavg", "SoilMoi00_10cm_tavg", "Evap_tavg")) { M <- ex(fl, v, paste0("_", vv, "$"))
    if (!is.null(M) && ncol(M) == 12) X <- addf(X, mstats(M, paste0("fldas_", san(sub("_tavg$", "", vv))), cl$month_med), "Climate and weather", "dynamic", "FLDAS (GEE)") }
  ch <- list.files(file.path("data/external_cache/chirps", cf$iso), pattern = "^chirps_[0-9]{4}_[0-9]{2}[.]tif$", full.names = TRUE, recursive = TRUE)
  if (length(ch) == 12) { M <- do.call(cbind, lapply(sort(ch), function(f) ex(f, v))); if (!is.null(M) && ncol(M) == 12) X <- addf(X, mstats(M, "chirps", cl$month_med), "Climate and weather", "dynamic", "CHIRPS") }
  X <- add(X, ex(pick_year(d, "NDVI", tg, Y), v), "ndvi_y", "Ecosystem productivity/greenness", "dynamic", "MODIS NDVI (GEE)")
  ndf <- list.files(d, pattern = paste0("^NDVI_", tg, "_[0-9]{4}[.]tif$"), full.names = TRUE)
  if (length(ndf) >= 3) { A <- do.call(cbind, lapply(ndf, function(f) ex(f, v))); if (!is.null(A)) {
    X <- addf(X, data.frame(ndvi_clim_mean = rowMeans(A, na.rm = TRUE), ndvi_clim_sd = apply(A, 1, stats::sd, na.rm = TRUE)), "Ecosystem productivity/greenness", "static", "MODIS NDVI (GEE)")
    if ("ndvi_y_ndvi" %in% names(X)) X <- addf(X, data.frame(ndvi_y_anom = X$ndvi_y_ndvi - X$ndvi_clim_mean), "Ecosystem productivity/greenness", "dynamic", "MODIS NDVI (GEE)") } }
  X <- add(X, ex(pick_year(d, "DailyEVI", tg, Y), v), "evi", "Ecosystem productivity/greenness", "dynamic", "MODIS EVI (GEE)")
  X <- add(X, ex(pick_year(d, "LAI8days", tg, Y), v), "lai", "Ecosystem productivity/greenness", "dynamic", "MODIS LAI (GEE)")
  X <- add(X, ex(pick_year(d, "Productivity", tg, Y), v, "_(Gpp|Npp)$"), "prod", "Ecosystem productivity/greenness", "dynamic", "MODIS GPP/NPP (GEE)")
  M <- ex(pick_year(d, "WAPOR", tg, Y), v, "NPP"); if (!is.null(M) && ncol(M) %in% c(12, 36)) { if (ncol(M) == 36) M <- sapply(1:12, function(k) rowMeans(M[, (3 * (k - 1) + 1):(3 * k), drop = FALSE], na.rm = TRUE)); X <- addf(X, mstats(M, "wapor_npp", cl$month_med), "Ecosystem productivity/greenness", "dynamic", "WaPOR (GEE)") }
  X <- add(X, ex(pick_year(d, "TRMM", tg, Y), v), "trmm_y", "Climate and weather", "dynamic", "TRMM (GEE)")
  trf <- list.files(d, pattern = paste0("^TRMM_", tg, "_[0-9]{4}[.]tif$"), full.names = TRUE)
  if (length(trf) >= 3) { A <- do.call(cbind, lapply(trf, function(f) ex(f, v))); if (!is.null(A)) X <- addf(X, data.frame(trmm_clim_mean = rowMeans(A, na.rm = TRUE), trmm_clim_cv = apply(A, 1, stats::sd, na.rm = TRUE) / pmax(rowMeans(A, na.rm = TRUE), 1e-6)), "Climate and weather", "static", "TRMM (GEE)") }
  X <- add(X, ex(pick_year(d, "AerosolOptical", tg, Y), v), "aod", "Climate and weather", "dynamic", "MODIS AOD (GEE)")
  X <- add(X, ex(pick_year(d, "PopDensity", tg, Y), v), "popdens", "Ruralness, population density, built environment", "slow", "GPW (GEE)")
  m <- ex(find_file(file.path("data/external_cache/worldpop", cf$iso), "^worldpop_"), v); if (!is.null(m)) { colnames(m) <- ""; X <- add(X, m, "worldpop", "Ruralness, population density, built environment", "slow", "WorldPop") }   # one shared column; the band name is the country file name, which gave worldpop_<file> per country (all-NA elsewhere) in the 2026-09-04 extraction
  X <- add(X, ex(find_file(file.path("data/external_cache/nightlights", cf$iso), "^ntl_"), v), "ntl_viirs", "Ruralness, population density, built environment", "slow", "VIIRS")
  # ── modelled surfaces and malaria ────────────────────────────────────────
  for (f in list.files(file.path("data/external_cache/malaria_atlas", cf$iso), pattern = "^(Malaria|Interventions|Blood_Disorders)__.*[.]tif$", full.names = TRUE)) {
    m <- ex(f, v); if (!is.null(m)) { colnames(m) <- ""; X <- add(X, m, paste0("map_", san(sub("[.]tif$", "", basename(f)))), if (grepl("^Interventions", basename(f))) "Malaria incidence and treatment" else "Malaria", "slow", "Malaria Atlas") } }
  for (k in c("Stunting", "Wasting", "Underweight")) { f <- ihme_file(sprintf("data/IHME/CGF/GeoTIFF/%s Prevalence [GeoTIFF]", k), "_MEAN_", Y); m <- ex(f, v); if (!is.null(m)) { colnames(m) <- ""; X <- add(X, m, paste0("ihme_", tolower(k)), "Nutrition status (MODELLED SURFACE)", "slow", "IHME LBD") } }
  f <- ihme_file("data/IHME/Anemia/GeoTIFF/3 - All Anemia [GeoTIFF]/Mean", "[.](tif|TIF)$", Y); m <- ex(f, v); if (!is.null(m)) { colnames(m) <- ""; X <- add(X, m, "ihme_anemia_all", "Nutrition status (MODELLED SURFACE)", "slow", "IHME LBD") }
  ALL[[cn]] <- X
  cat(sprintf("   %s done: %d columns in %.1f min\n", cn, ncol(X) - 8, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
P <- bind_rows(ALL); write.csv(P, file.path(CDIR, paste0("predictors_cluster", OT, ".csv")), row.names = FALSE)
MD <- bind_rows(META) |> distinct(column, .keep_all = TRUE) |> filter(column %in% names(P))
MD$completeness <- round(vapply(MD$column, function(cc) mean(is.finite(P[[cc]])), 0), 3)
MD$countries <- vapply(MD$column, function(cc) paste(names(which(tapply(is.finite(P[[cc]]), P$country, mean) > 0.5)), collapse = "|"), "")
MD$n_countries <- lengths(strsplit(MD$countries, "[|]"))
write.csv(MD, file.path(CDIR, paste0("predictors_cluster", OT, "_metadata.csv")), row.names = FALSE)
cat("\n===== CL-02: cluster covariates =====\n")
cat(sprintf("clusters %d | columns %d | by domain:\n", nrow(P), nrow(MD))); print(table(MD$domain, MD$role))
cat(sprintf("columns present in all four countries: %d of %d\n", sum(MD$n_countries == 4), nrow(MD)))
cat("\nDONE\n")
