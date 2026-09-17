# =============================================================================
# scripts/policy_deck/03b_civ_extra_domains.R   [CV-01, 2026-09-17]
#
# THE OPEN-RASTER DOMAINS FOR COTE D'IVOIRE, UNDER THE SHARED COLUMN NAMES
#
# results/transportability/civ_canonical_admin2_full.rds holds 132 columns:
# the GEE climate, soil and greenness layers. The five-domain transport
# candidate (DA-03: climate, soil, anaemia surfaces, agriculture, infection)
# and the climate-normals block need more. Everything here is an open raster
# or a public grid, extracted for the 33 CIV Admin-2 polygons (script 03) by
# the same code paths the training countries use:
#
#   IHME surfaces      script 48's list, WorldPop-weighted zonal means at the
#                      nearest year to CIV_REF_YEAR         -> ihme_*
#   livestock (GLW4)   script 48's block, 2020 densities      -> glw_*
#   Malaria Atlas      script 03d's cache (survey-year + static) -> map_sy_*, map_blooddisorders*
#   MapSPAM            scripts/build_mapspam_admin2.R CoteDIvoire -> spam_*
#   climate normals    script 03c's GEE climatology + build_climate_normals_block.R formulas -> clim_*
#
# CIV has no biomarker survey, so CIV_REF_YEAR (default 2016, the MICS 2016
# year) stands in for the survey year wherever a layer is time-varying.
# Blocks whose inputs are missing are skipped and listed, never invented.
#
#   Rscript scripts/policy_deck/03b_civ_extra_domains.R
# -> results/transportability/civ_canonical_admin2_full_v2.rds   (132 + new columns)
#    results/tables/policy_deck/civ_extra_domains_manifest.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf); library(terra); library(exactextractr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
REF <- as.integer(Sys.getenv("CIV_REF_YEAR", "2016"))
CIV <- readRDS("results/transportability/civ_canonical_admin2_full.rds")
P <- sf::st_read("data/external_cache/gee_geoms/civ_admin2.gpkg", quiet = TRUE)
stopifnot(identical(paste(P$Admin1, P$Admin2), paste(CIV$Admin1, CIV$Admin2)))
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
pop <- rast("data/external_cache/worldpop/CIV/worldpop_CIV_2017.tif")
OUT <- data.frame(Admin1 = CIV$Admin1, Admin2 = CIV$Admin2, stringsAsFactors = FALSE)
manifest <- list()
note <- function(block, column, status, detail = "") manifest[[length(manifest) + 1L]] <<- data.frame(block = block, column = column, status = status, detail = detail, stringsAsFactors = FALSE)

# ── IHME surfaces (the list of script 48) ────────────────────────────────────
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
  ihme_u5_diarrhoea_prev = list(dir = "data/IHME/u5 diarrhea/GeoTIFF", pat = "DIARRHEA_2000_2017_PREV_RT_MEAN_%d_", yrs = 2000:2017),
  ihme_oralrehydrationsolution = list(dir = "data/IHME/oral rehydration/GeoTIFF", pat = "_ORS_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_orsorrhf                = list(dir = "data/IHME/oral rehydration/GeoTIFF", pat = "_ORT_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_recommendedhomefluids   = list(dir = "data/IHME/oral rehydration/GeoTIFF", pat = "_RHF_PERCENT_MEAN_%d_", yrs = 2000:2017),
  ihme_meanyearsofattainment = list(file = "data/IHME/education/GeoTIFF/IHME_LMIC_EDU_2000_2017_MEAN_15_49_FEMALE_MEAN_Y2019M12D24.TIF", band_year0 = 1999, yrs = 2000:2017),
  ihme_edu_0y_share          = list(file = "data/IHME/education/GeoTIFF/IHME_LMIC_EDU_2000_2017_ZEROPROP_15_49_FEMALE_MEAN_Y2019M12D24.TIF", band_year0 = 1999, yrs = 2000:2017),
  ihme_edu_6_11y_share       = list(file = "data/IHME/education/GeoTIFF/IHME_LMIC_EDU_2000_2017_PRIMARYPROP_15_49_FEMALE_MEAN_Y2019M12D24.TIF", band_year0 = 1999, yrs = 2000:2017),
  ihme_edu_12plus_share      = list(file = "data/IHME/education/GeoTIFF/IHME_LMIC_EDU_2000_2017_SECONDARYPROP_1549_FEMALE_MEAN_Y2019M12D24.TIF", band_year0 = 1999, yrs = 2000:2017),
  ihme_mcvcoverage = list(file = "data/IHME/mcv1/Data [GeoTIFF]/IHME_LMIC_MCV1_2000_2019_MEAN_Y2020M12D16.TIF", band_year0 = 1999, yrs = 2000:2019)
)
find_file <- function(spec, yr) {
  if (!is.null(spec$file)) return(if (file.exists(spec$file)) spec$file else NA_character_)
  f <- list.files(spec$dir, pattern = sprintf(spec$pat, yr), recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[grepl("[.]tif$", f, ignore.case = TRUE) & !grepl("LOWER|UPPER", basename(f))]
  if (length(f)) f[1] else NA_character_ }
load_surface <- function(spec, yr, bbox) {
  f <- find_file(spec, yr); if (is.na(f)) return(NULL)
  r <- rast(f); if (!is.null(spec$band_year0)) { k <- yr - spec$band_year0; if (k < 1 || k > nlyr(r)) return(NULL); r <- r[[k]] }
  r <- crop(r, ext(bbox) + 0.25); names(r) <- "v"; r }
WCACHE <- new.env()
weights_for <- function(r) {
  key <- paste(paste(res(r), collapse = "x"), paste(as.vector(ext(r)), collapse = "_"), paste(dim(r)[1:2], collapse = "x"))
  if (!is.null(WCACHE[[key]])) return(WCACHE[[key]])
  pc <- crop(pop, ext(r), snap = "out"); fact <- max(1L, round(res(r)[1] / res(pc)[1]))
  pa <- if (fact > 1) aggregate(pc, fact = fact, fun = "sum", na.rm = TRUE) else pc
  w <- resample(pa, r, method = "near"); w <- subst(w, NA, 0); names(w) <- "w"; WCACHE[[key]] <- w; w }
zonal_w <- function(r, w) {
  aligned <- isTRUE(all.equal(as.vector(ext(w)), as.vector(ext(r)), tolerance = 1e-12)) && all(dim(w)[1:2] == dim(r)[1:2])
  if (!aligned) w <- subst(resample(w, r, method = "near"), NA, 0)
  v <- exact_extract(r, P, "weighted_mean", weights = w, progress = FALSE)
  m <- exact_extract(r, P, "mean", progress = FALSE); v[!is.finite(v)] <- m[!is.finite(v)]; v }
# the scale each IHME column carries in the shared set (script 48 put the raster on the tabular scale)
IHMD <- tryCatch(read.csv("data/covariates/harmonized/predictors_admin2_ihme_raster_metadata.csv"), error = function(e) NULL)
scale_of <- function(key) { if (is.null(IHMD) || !"scale_applied" %in% names(IHMD)) return(1); s <- IHMD$scale_applied[IHMD$column == key]; s <- s[is.finite(s)]; if (length(s)) s[1] else 1 }
cat("===== IHME surfaces =====\n")
for (key in names(IH)) {
  if (!key %in% MD$column) { note("IHME", key, "skipped", "not in the shared set"); next }
  spec <- IH[[key]]; yr <- spec$yrs[which.min(abs(spec$yrs - REF))]
  r <- load_surface(spec, yr, P); if (is.null(r)) { note("IHME", key, "missing", sprintf("no surface for %d", yr)); cat(sprintf("  %-30s no surface for %d\n", key, yr)); next }
  v <- zonal_w(r, weights_for(r)) * scale_of(key); OUT[[key]] <- v
  note("IHME", key, "ok", sprintf("year %d, scale %g, finite %.2f", yr, scale_of(key), mean(is.finite(v)))) }
cat(sprintf("  %d IHME columns\n", sum(grepl("^ihme_", names(OUT)))))

# ── livestock (GLW4 2020, the block of record) ──────────────────────────────
cat("===== GLW4 livestock =====\n")
TLU <- c(cattle = 0.7, sheep = 0.1, goats = 0.1, pigs = 0.2, chickens = 0.01)
pop_a2 <- exact_extract(pop, P, "sum", progress = FALSE); tlu <- 0; heads_tlu <- 0; rum <- 0
for (sp in names(TLU)) {
  f <- sprintf("data/external_cache/glw4_2020/glw4_2020_%s_density.tif", sp)
  if (!file.exists(f)) { note("GLW4", paste0("glw_", sp, "_km2"), "missing", f); next }
  dens <- crop(rast(f), ext(P) + 0.5); r <- dens * cellSize(dens, unit = "km"); names(dens) <- "d"
  da <- exact_extract(dens, P, "mean", progress = FALSE); ha <- exact_extract(r, P, "sum", progress = FALSE)
  OUT[[paste0("glw_", sp, "_km2")]] <- da; note("GLW4", paste0("glw_", sp, "_km2"), "ok")
  tlu <- tlu + TLU[[sp]] * da; heads_tlu <- heads_tlu + TLU[[sp]] * ha; if (sp %in% c("cattle", "sheep", "goats")) rum <- rum + TLU[[sp]] * da }
OUT$glw_tlu_km2 <- tlu; OUT$glw_tlu_per_capita <- ifelse(pop_a2 > 0, heads_tlu / pop_a2, NA_real_); OUT$glw_ruminant_share <- ifelse(tlu > 0, rum / tlu, NA_real_)
for (k in c("glw_tlu_km2", "glw_tlu_per_capita", "glw_ruminant_share")) note("GLW4", k, "ok")

# ── Malaria Atlas (script 03d's cache) ───────────────────────────────────────
cat("===== Malaria Atlas =====\n")
MCACHE <- "data/external_cache/malaria_atlas_sy/CoteDIvoire"
mp <- MD$column[grepl("^map_", MD$column)]
for (col in mp) {
  fs <- list.files(MCACHE, pattern = "[.]tif$", full.names = TRUE)
  f <- if (grepl("^map_sy_", col)) {
    stem <- switch(col, map_sy_pf_parasite_rate = "Pf_Parasite_Rate", map_sy_pf_incidence_rate = "Pf_Incidence_Rate", map_sy_pf_mortality_rate = "Pf_Mortality_Rate",
                   map_sy_pf_reproductive_number = "Pf_Reproductive_Number", map_sy_itn_access = "Net_Access", map_sy_itn_use = "Net_Use_", map_sy_itn_use_rate = "Net_Use_Rate",
                   map_sy_irs_coverage = "IRS_Coverage", map_sy_effective_treatment = "Effective_Treatment", NA)
    if (col == "map_sy_itn_use") fs[grepl("Net_Use_[0-9]", basename(fs))] else fs[grepl(stem, basename(fs), fixed = TRUE)]
  } else fs[tolower(gsub("[^A-Za-z0-9]", "", sub("[.]tif$", "", basename(fs)))) == sub("^map_", "", col)]
  if (!length(f)) { note("MAP", col, "missing", "not in the CIV cache (run 03d)"); next }
  r <- rast(f[1]); v <- tryCatch(terra::extract(r[[1]], vect(st_transform(P, crs(r))), fun = mean, na.rm = TRUE, ID = FALSE)[, 1], error = function(e) NULL)
  if (is.null(v)) { note("MAP", col, "failed"); next }
  OUT[[col]] <- as.numeric(v); note("MAP", col, "ok", basename(f[1])) }
cat(sprintf("  %d Malaria Atlas columns\n", sum(grepl("^map_", names(OUT)))))

# ── MapSPAM ──────────────────────────────────────────────────────────────────
cat("===== MapSPAM =====\n")
spf <- "data/MapSPAM/CoteDIvoire_spam_admin2.csv"
if (file.exists(spf)) { sp <- read.csv(spf, check.names = FALSE); j <- match(paste(OUT$Admin1, OUT$Admin2), paste(sp$Admin1, sp$Admin2))
  for (col in intersect(MD$column[grepl("^spam_", MD$column)], names(sp))) { OUT[[col]] <- sp[[col]][j]; note("MapSPAM", col, "ok") }
  cat(sprintf("  %d columns, %d districts matched\n", sum(grepl("^spam_", names(OUT))), sum(is.finite(j))))
} else note("MapSPAM", "spam_*", "missing", "run scripts/build_mapspam_admin2.R CoteDIvoire")

# ── climate normals (script 03c's GEE climatology, build_climate_normals_block.R formulas) ─
cat("===== climate normals =====\n")
cf <- "data/external_cache/gee_geoms/civ/gee_climatology_admin2.csv"
if (file.exists(cf)) {
  C <- read.csv(cf, check.names = FALSE); j <- match(paste(OUT$Admin1, OUT$Admin2), paste(C$Admin1, C$Admin2)); C <- C[j, ]
  aw <- function(v) C[[paste0(v, "_aw")]]; mon <- function(p) sapply(sprintf("%s_m%02d", p, 1:12), aw)
  prm <- mon("pr"); txm <- mon("tmax"); ldm <- mon("lstd"); lnm <- mon("lstn")
  top3 <- t(apply(prm, 1, function(x) { s <- sort(x, decreasing = TRUE); c(sum(s[1:3]), sum(x)) }))
  CL <- data.frame(clim_pr_ann_mean = aw("pr_ann_mean"), clim_pr_ann_sd = aw("pr_ann_sd"), clim_pr_cv = aw("pr_ann_sd") / aw("pr_ann_mean"),
    clim_pr_season_cv = apply(prm, 1, sd) / rowMeans(prm), clim_pr_top3_share = top3[, 1] / top3[, 2],
    clim_pr_sy_anom_z = (aw("pr_sy") - aw("pr_ann_mean")) / aw("pr_ann_sd"), clim_pr_win_anom_z = (aw("pr_win") - aw("pr_ann_mean")) / aw("pr_ann_sd"),
    clim_tmax_ann = rowMeans(txm), clim_tmax_range = apply(txm, 1, max) - apply(txm, 1, min), clim_tmax_sy_anom = aw("tmax_sy") - rowMeans(txm),
    clim_tmin_ann = aw("tmin_ann"), clim_pet_ann = aw("pet_ann"), clim_def_ann = aw("def_ann"), clim_aet_ann = aw("aet_ann"),
    clim_soilm_ann = aw("soilm_ann"), clim_vpd_ann = aw("vpd_ann"), clim_srad_ann = aw("srad_ann"),
    clim_lstd_ann = rowMeans(ldm), clim_lstn_ann = rowMeans(lnm), clim_lst_diurnal = rowMeans(ldm) - rowMeans(lnm))
  for (col in intersect(MD$column, names(CL))) { OUT[[col]] <- CL[[col]]; note("climate normals", col, "ok") }
  cat(sprintf("  %d clim_ columns\n", sum(grepl("^clim_", names(OUT)))))
} else note("climate normals", "clim_*", "missing", "run scripts/policy_deck/03c_civ_climatology_gee.py")

# ── assemble ─────────────────────────────────────────────────────────────────
new <- setdiff(names(OUT), c("Admin1", "Admin2", names(CIV)))
V2 <- cbind(CIV, OUT[, new, drop = FALSE])
saveRDS(V2, "results/transportability/civ_canonical_admin2_full_v2.rds")
MF <- bind_rows(manifest); write.csv(MF, "results/tables/policy_deck/civ_extra_domains_manifest.csv", row.names = FALSE)
cat(sprintf("\nwritten civ_canonical_admin2_full_v2.rds: %d rows x %d columns (%d new)\n", nrow(V2), ncol(V2), length(new)))
print(table(MF$block, MF$status))
dom <- MD$domain[match(new, MD$column)]; print(sort(table(dom), decreasing = TRUE))
cat("DONE\n")
