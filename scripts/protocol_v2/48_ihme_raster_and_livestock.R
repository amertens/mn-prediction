# =============================================================================
# scripts/protocol_v2/48_ihme_raster_and_livestock.R   [IH-02, LV-01]
#
# 1. IHME MODELLED SURFACES BY ZONAL MEAN, NOT BY NAME.
#    The extra-source block (script 08) joins IHME's Admin-2 rollups to the
#    spine by district NAME. Ghana's 2019 district splits leave 26% of its
#    surveyed districts without a value, and Malawi can only be joined at the
#    district (this project's Admin-1) and broadcast to its Traditional
#    Authorities. The same products exist as 5 km GeoTIFF surfaces on disk
#    (data/IHME/*/GeoTIFF), so here every indicator that has a surface is
#    rebuilt as a WorldPop-weighted zonal mean over the spine polygons, at the
#    nearest available year to each survey, under the SAME column names the
#    tabular block uses (ihme_stuntingprevalence, ihme_allanemia, ...). Script
#    08 swaps these in where present and keeps the tabular value only for the
#    six indicators that have no surface (education, ORS, circumcision).
#    The two versions are compared where both exist (Spearman by country), and
#    the raster scale is put on the tabular scale automatically.
#
# 2. LIVESTOCK DENSITY (Gridded Livestock of the World 4, 2015, 5 arc-minutes,
#    dasymetric; Harvard Dataverse, CC-BY 4.0; data/external_cache/glw4/).
#    Animal-source foods are the main dietary route for iron, zinc, B12 and
#    preformed vitamin A; the vocabulary so far carries only DHS livestock
#    OWNERSHIP, which is missing for Sierra Leone. Columns: head per km2 for
#    cattle, sheep, goats, pigs, chickens; tropical livestock units per km2
#    (0.7 / 0.1 / 0.1 / 0.2 / 0.01); TLU per person (WorldPop); ruminant
#    share of TLU.
#
# Both blocks are also extracted at the survey clusters (2 km urban / 5 km
# rural buffers) for the cluster track.
#
#   Rscript scripts/protocol_v2/48_ihme_raster_and_livestock.R
# -> data/covariates/harmonized/predictors_admin2_ihme_raster.csv (+ _metadata)
#    data/covariates/harmonized/predictors_admin2_livestock.csv
#    data/covariates/cluster/predictors_cluster_ihme_raster.csv
#    data/covariates/cluster/predictors_cluster_livestock.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf); library(terra); library(exactextractr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"; CDIR <- "data/covariates/cluster"
source("R/survey_years.R"); SURVEY_YEAR <- survey_years()   # single source: metadata/survey_years.csv (Gambia 2018, Ghana 2017, Malawi 2016, Sierra Leone 2013)
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
LC  <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
polys <- function(cn) { b <- st_as_sf(BND[[LC[[cn]]]]); b <- st_make_valid(st_transform(b, 4326))
  st_sf(country = cn, Admin1 = as.character(b$Admin1), Admin2 = as.character(b$Admin2), geometry = st_geometry(b)) }
pop_of <- function(cn) rast(list.files(file.path("data/external_cache/worldpop", ISO[[cn]]), pattern = "[.]tif$", full.names = TRUE)[1])
TC <- read.csv("results/tables/cluster_level/targets_cluster.csv") |> distinct(country, cluster, lat, lon) |> filter(is.finite(lat), is.finite(lon))
# buffer radius: the GHSL SMOD urban flag the cluster extraction (script 02) uses, 2 km urban / 5 km rural by default
RURAL_KM <- as.numeric(Sys.getenv("CL_RURAL_KM", "5")); URBAN_KM <- as.numeric(Sys.getenv("CL_URBAN_KM", "2")); OT <- Sys.getenv("CL_OUT_TAG", "")
PCU <- read.csv(file.path(CDIR, "predictors_cluster.csv"), check.names = FALSE)[, c("country", "cluster", "urban")]
TC$cluster <- as.character(TC$cluster); PCU$cluster <- as.character(PCU$cluster); TC <- left_join(TC, PCU, by = c("country", "cluster"))
TC$radius_km <- ifelse(TC$urban %in% c(1, TRUE), URBAN_KM, RURAL_KM)
buffers <- function(cn) { t <- TC[TC$country == cn, ]; p <- st_as_sf(t, coords = c("lon", "lat"), crs = 4326)
  st_geometry(p) <- st_buffer(st_geometry(p), t$radius_km * 1000); p }

# ── IHME surfaces: column name -> where the mean surface lives ───────────────
IH <- list(
  ihme_stuntingprevalence    = list(dir = "data/IHME/CGF/GeoTIFF", pat = "STUNTING_PREV_MEAN_%d_", yrs = 2000:2017),
  ihme_wastingprevalence     = list(dir = "data/IHME/CGF/GeoTIFF", pat = "WASTING_PREV_MEAN_%d_", yrs = 2000:2017),
  ihme_underweightprevalence = list(dir = "data/IHME/CGF/GeoTIFF", pat = "UNDERWEIGHT_PREV_MEAN_%d_", yrs = 2000:2017),
  ihme_overweightprevalence  = list(dir = "data/IHME/DBM/GeoTIFF", pat = "OVERWEIGHT_PREV_MEAN_%d_", yrs = 2000:2017),
  ihme_allanemia      = list(dir = "data/IHME/Anemia/GeoTIFF", pat = "ALL_ANEMIA_PREV_PERCENT_MEAN_%d_", yrs = 2000:2019),
  ihme_mildanemia     = list(dir = "data/IHME/Anemia/GeoTIFF", pat = "MILD_ANEMIA_PREV_MEAN_%d_", yrs = 2000:2019),
  ihme_moderateanemia = list(dir = "data/IHME/Anemia/GeoTIFF", pat = "MOD(ERATE)?_ANEMIA_PREV_MEAN_%d_", yrs = 2000:2019),
  ihme_severeanemia   = list(dir = "data/IHME/Anemia/GeoTIFF", pat = "SEV(ERE)?_ANEMIA_PREV_MEAN_%d_", yrs = 2000:2019),
  ihme_ebfprevalence  = list(dir = "data/IHME/EBF/GeoTIFF", pat = "EBF_2000_2019_PREV_PERCENT_MEAN_%d_", yrs = 2000:2019),
  ihme_hivprevalence  = list(dir = "data/IHME/HIV/HIV Prevalence [GeoTIFF]/Both sexes 15-49", pat = "HIV_PREVALENCE_MEAN_15_49_BOTH_%d_", yrs = 2000:2018),
  ihme_simp      = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_S_IMP_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_simpother = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_S_IMP_OTHER_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_sod       = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_S_OD_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_spiped    = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_S_PIPED_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_sunimp    = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_S_UNIMP_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_wimp      = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_W_IMP_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_wimpother = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_W_IMP_OTHER_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_wpiped    = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_W_PIPED_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_wsurface  = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_W_SURFACE_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_wunimp    = list(dir = "data/IHME/WASH access/GeoTIFF", pat = "_W_UNIMP_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_incidence  = list(dir = "data/IHME/u5 diarrhea/GeoTIFF", pat = "DIARRHEA_2000_2017_INC_RT_MEAN_%d_", yrs = 2000:2017),
  ihme_deaths     = list(dir = "data/IHME/u5 diarrhea/GeoTIFF", pat = "DIARRHEA_2000_2017_MORT_RT_MEAN_%d_", yrs = 2000:2017),
  ihme_prevalence = list(dir = "data/IHME/u5 diarrhea/GeoTIFF", pat = "DIARRHEA_2000_2017_PREV_RT_MEAN_%d_", yrs = 2000:2017),
  ihme_mcvcoverage = list(file = "data/IHME/mcv1/Data [GeoTIFF]/IHME_LMIC_MCV1_2000_2019_MEAN_Y2020M12D16.TIF", band_year0 = 1999, yrs = 2000:2019)
)
find_file <- function(spec, yr) {
  if (!is.null(spec$file)) return(spec$file)
  f <- list.files(spec$dir, pattern = sprintf(spec$pat, yr), recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[grepl("[.]tif$", f, ignore.case = TRUE) & !grepl("LOWER|UPPER", basename(f))]
  if (length(f)) f[1] else NA_character_ }
load_surface <- function(spec, yr, bbox) {
  f <- find_file(spec, yr); if (is.na(f)) return(NULL)
  r <- rast(f); if (!is.null(spec$band_year0)) { k <- yr - spec$band_year0; if (k < 1 || k > nlyr(r)) return(NULL); r <- r[[k]] }
  r <- crop(r, ext(bbox) + 0.25); names(r) <- "v"; r }
# WorldPop resampled onto a surface's grid, cached per grid, so a 100 m population raster is aggregated once per resolution
WCACHE <- new.env()
weights_for <- function(r, pop, cn) {
  # full-precision key: IHME products share a nominal 5 km grid but not an origin (offsets of 1e-5 degrees), and
  # exactextractr requires the weight grid to match the value grid exactly, tighter than terra's compareGeom tolerance
  key <- paste(cn, paste(res(r), collapse = "x"), paste(as.vector(ext(r)), collapse = "_"), paste(dim(r)[1:2], collapse = "x"))
  if (!is.null(WCACHE[[key]])) return(WCACHE[[key]])
  pc <- crop(pop, ext(r), snap = "out"); fact <- max(1L, round(res(r)[1] / res(pc)[1]))
  pa <- if (fact > 1) aggregate(pc, fact = fact, fun = "sum", na.rm = TRUE) else pc
  w <- resample(pa, r, method = "near"); w <- subst(w, NA, 0); names(w) <- "w"; WCACHE[[key]] <- w; w }
zonal_w <- function(r, w, P) {
  aligned <- isTRUE(all.equal(as.vector(ext(w)), as.vector(ext(r)), tolerance = 1e-12)) && all(dim(w)[1:2] == dim(r)[1:2])
  if (!aligned) w <- subst(resample(w, r, method = "near"), NA, 0)
  v <- exact_extract(r, P, "weighted_mean", weights = w, progress = FALSE)
  m <- exact_extract(r, P, "mean", progress = FALSE); v[!is.finite(v)] <- m[!is.finite(v)]; v }

EX <- read.csv(file.path(HDIR, "predictors_admin2_extrasrc.csv"), check.names = FALSE)
A2 <- list(); CL <- list(); META <- list()
cat("===== IH-02: IHME surfaces, zonal means =====\n")
for (cn in names(SURVEY_YEAR)) {
  P <- polys(cn); B <- buffers(cn); pop <- pop_of(cn)
  oa <- st_drop_geometry(P); oc <- st_drop_geometry(B)[, c("country", "cluster")]
  for (key in names(IH)) { spec <- IH[[key]]; yr <- spec$yrs[which.min(abs(spec$yrs - SURVEY_YEAR[[cn]]))]
    r <- load_surface(spec, yr, P); if (is.null(r)) { cat(sprintf("  %-12s %-28s no surface for %d\n", cn, key, yr)); next }
    w <- weights_for(r, pop, cn)
    oa[[key]] <- zonal_w(r, w, P); oc[[key]] <- zonal_w(r, w, B)
    META[[length(META) + 1L]] <- data.frame(country = cn, column = key, year = yr, file = basename(sources(r)[1]), stringsAsFactors = FALSE) }
  cat(sprintf("  %-12s %d surfaces extracted (%d districts, %d clusters)\n", cn, ncol(oa) - 3, nrow(oa), nrow(oc)))
  A2[[cn]] <- oa; CL[[cn]] <- oc }
IA <- bind_rows(A2); IC <- bind_rows(CL)
# ── put the raster on the tabular scale, and compare the two where both exist ─
cat("\n-- raster vs tabular (name-joined) IHME, same column, matched districts --\n")
cmp <- list()
for (key in setdiff(names(IA), c("country", "Admin1", "Admin2"))) {
  if (!key %in% names(EX)) next
  j <- inner_join(IA[, c("country", "Admin1", "Admin2", key)], EX[, c("country", "Admin1", "Admin2", key)], by = c("country", "Admin1", "Admin2"), suffix = c(".r", ".t"))
  ok <- is.finite(j[[paste0(key, ".r")]]) & is.finite(j[[paste0(key, ".t")]]); if (sum(ok) < 10) next
  ratio <- median(j[[paste0(key, ".t")]][ok]) / median(j[[paste0(key, ".r")]][ok])
  scale <- if (is.finite(ratio) && ratio > 0.005 && ratio < 0.02) 0.01 else if (is.finite(ratio) && ratio > 50 && ratio < 200) 100 else 1
  if (scale != 1) { IA[[key]] <- IA[[key]] * scale; IC[[key]] <- IC[[key]] * scale; j[[paste0(key, ".r")]] <- j[[paste0(key, ".r")]] * scale }
  by_cn <- j[ok, ] |> group_by(country) |> summarise(n = n(), rho = suppressWarnings(cor(.data[[paste0(key, ".r")]], .data[[paste0(key, ".t")]], method = "spearman")), .groups = "drop")
  cmp[[key]] <- data.frame(column = key, scale_applied = scale, median_ratio = round(ratio, 3), n = sum(ok), rho_all = round(suppressWarnings(cor(j[[paste0(key, ".r")]][ok], j[[paste0(key, ".t")]][ok], method = "spearman")), 3),
                           rho_by_country = paste(sprintf("%s %.2f (n=%d)", by_cn$country, by_cn$rho, by_cn$n), collapse = "; "), stringsAsFactors = FALSE) }
CMP <- bind_rows(cmp); print(CMP, row.names = FALSE)
write.csv(IA, file.path(HDIR, "predictors_admin2_ihme_raster.csv"), row.names = FALSE)
write.csv(IC, file.path(CDIR, paste0("predictors_cluster_ihme_raster", OT, ".csv")), row.names = FALSE)
MD <- bind_rows(META) |> left_join(CMP[, c("column", "scale_applied", "rho_all")], by = "column")
write.csv(MD, file.path(HDIR, "predictors_admin2_ihme_raster_metadata.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %d Admin-2 rows x %d IHME columns; %d cluster rows\n", nrow(IA), ncol(IA) - 3, nrow(IC)))
for (cn in names(SURVEY_YEAR)) { m <- as.matrix(IA[IA$country == cn, setdiff(names(IA), c("country", "Admin1", "Admin2"))]); cat(sprintf("  %-12s finite share %.3f\n", cn, mean(is.finite(m)))) }

# ── LV-01: livestock density ─────────────────────────────────────────────────
# GLW_YEAR=2015: Harvard Dataverse dasymetric HEAD COUNTS per 5-arc-min cell (data/external_cache/glw4/)
# GLW_YEAR=2020 (default): FAO catalog GLW4-2020 D-DA rasters, already head per km2 (data/external_cache/glw4_2020/)
# 2020 is the block of record since 2026-09-07; district rankings agree with 2015 at Spearman 0.99-1.00
GLW_YEAR <- Sys.getenv("GLW_YEAR", "2020"); GLW_TAG <- Sys.getenv("GLW_TAG", "")
cat(sprintf("\n===== LV-01: Gridded Livestock of the World 4 (%s) =====\n", GLW_YEAR))
TLU <- c(cattle = 0.7, sheep = 0.1, goats = 0.1, pigs = 0.2, chickens = 0.01)
glw_raster <- function(sp) if (GLW_YEAR == "2020") rast(sprintf("data/external_cache/glw4_2020/glw4_2020_%s_density.tif", sp)) else rast(sprintf("data/external_cache/glw4/glw4_%s_2015_Da.tif", sp))
LA <- list(); LCL <- list()
for (cn in names(SURVEY_YEAR)) {
  P <- polys(cn); B <- buffers(cn); pop <- pop_of(cn)
  oa <- st_drop_geometry(P); oc <- st_drop_geometry(B)[, c("country", "cluster")]
  pop_a2 <- exact_extract(pop, P, "sum", progress = FALSE); pop_cl <- exact_extract(pop, B, "sum", progress = FALSE)
  tlu_a2 <- 0; tlu_cl <- 0; heads_tlu_a2 <- 0; heads_tlu_cl <- 0; rum_a2 <- 0; rum_cl <- 0
  for (sp in names(TLU)) {
    r <- crop(glw_raster(sp), ext(P) + 0.5)
    if (GLW_YEAR == "2020") { dens <- r; r <- dens * cellSize(dens, unit = "km") } else dens <- r / cellSize(r, unit = "km")   # r = heads per cell, dens = heads per km2
    names(dens) <- "d"
    da <- exact_extract(dens, P, "mean", progress = FALSE); dc <- exact_extract(dens, B, "mean", progress = FALSE)
    ha <- exact_extract(r, P, "sum", progress = FALSE); hc <- exact_extract(r, B, "sum", progress = FALSE)
    oa[[paste0("glw_", sp, "_km2")]] <- da; oc[[paste0("glw_", sp, "_km2")]] <- dc
    tlu_a2 <- tlu_a2 + TLU[[sp]] * da; tlu_cl <- tlu_cl + TLU[[sp]] * dc
    heads_tlu_a2 <- heads_tlu_a2 + TLU[[sp]] * ha; heads_tlu_cl <- heads_tlu_cl + TLU[[sp]] * hc
    if (sp %in% c("cattle", "sheep", "goats")) { rum_a2 <- rum_a2 + TLU[[sp]] * da; rum_cl <- rum_cl + TLU[[sp]] * dc } }
  oa$glw_tlu_km2 <- tlu_a2; oc$glw_tlu_km2 <- tlu_cl
  oa$glw_tlu_per_capita <- ifelse(pop_a2 > 0, heads_tlu_a2 / pop_a2, NA_real_); oc$glw_tlu_per_capita <- ifelse(pop_cl > 0, heads_tlu_cl / pop_cl, NA_real_)
  oa$glw_ruminant_share <- ifelse(tlu_a2 > 0, rum_a2 / tlu_a2, NA_real_); oc$glw_ruminant_share <- ifelse(tlu_cl > 0, rum_cl / tlu_cl, NA_real_)
  cat(sprintf("  %-12s TLU/km2 median %.1f (IQR %.1f-%.1f) | TLU per person median %.2f | ruminant share median %.2f\n", cn, median(oa$glw_tlu_km2, na.rm = TRUE), quantile(oa$glw_tlu_km2, 0.25, na.rm = TRUE), quantile(oa$glw_tlu_km2, 0.75, na.rm = TRUE), median(oa$glw_tlu_per_capita, na.rm = TRUE), median(oa$glw_ruminant_share, na.rm = TRUE)))
  LA[[cn]] <- oa; LCL[[cn]] <- oc }
LVA <- bind_rows(LA); LVC <- bind_rows(LCL)
write.csv(LVA, file.path(HDIR, paste0("predictors_admin2_livestock", GLW_TAG, ".csv")), row.names = FALSE)
write.csv(LVC, file.path(CDIR, paste0("predictors_cluster_livestock", GLW_TAG, OT, ".csv")), row.names = FALSE)
cat(sprintf("written: %d Admin-2 rows x %d livestock columns; %d cluster rows\n", nrow(LVA), ncol(LVA) - 3, nrow(LVC)))
cat("\nDONE\n")
