# =============================================================================
# scripts/protocol_v2/55_engineer_climate_soil.R   [FE-01]
#
# ENGINEERED CLIMATE, TERRAIN AND SOIL BLOCKS FROM THE SCRIPT-54 EXTRACTS
#
# Turns the raw zonal tables into six add-on blocks for the harness (script
# 39), each in an area-weighted (aw) and a population-weighted (pw) version:
#
#   climatology   replaces the 30 calendar slices and the one-year TerraClimate
#                 columns: annual means, Walsh-Lawler seasonality, first-harmonic
#                 amplitude and phase (sin, cos), wet / dry month counts,
#                 inter-annual CV, survey-year and fieldwork-window anomalies,
#                 diurnal range, hottest month, aridity, LST day / night
#                 climatology and range.
#   terrain       relief from MERIT Hydro and Geomorpho90m: elevation, HAND,
#                 log upstream area, slope, TRI, TPI, roughness, VRM, CTI.
#   soil          depth-weighted 0-50 cm properties in natural units, pH-
#                 conditioned zinc / iron / phosphorus availability, cation
#                 ratios, texture and structure indices, a base-saturation proxy
#                 and a fertility score; no within-district dispersion columns.
#
#   Rscript scripts/protocol_v2/55_engineer_climate_soil.R
# -> data/covariates/harmonized/predictors_admin2_fe_{climatology,terrain,soil}_{aw,pw}.csv (+ _metadata)
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"; KEY <- c("country", "Admin1", "Admin2")
rd <- function(f) read.csv(file.path(HDIR, f), check.names = FALSE, stringsAsFactors = FALSE)
pick <- function(D, suffix) { cols <- grep(paste0("_", suffix, "$"), names(D), value = TRUE); out <- D[, c(KEY, cols)]; names(out) <- c(KEY, sub(paste0("_", suffix, "$"), "", cols)); out }
sig <- function(x, centre, slope = 2) 1 / (1 + exp(slope * (x - centre)))      # 1 below `centre`, 0 above
zc  <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s }
write_block <- function(D, name, domain, source, note) {
  cols <- setdiff(names(D), KEY)
  write.csv(D, file.path(HDIR, sprintf("predictors_admin2_fe_%s.csv", name)), row.names = FALSE)
  md <- data.frame(column = cols, domain = domain, source = source, note = note, completeness = round(vapply(cols, function(v) mean(is.finite(D[[v]])), 0), 3), stringsAsFactors = FALSE)
  write.csv(md, file.path(HDIR, sprintf("predictors_admin2_fe_%s_metadata.csv", name)), row.names = FALSE)
  cat(sprintf("  %-18s %d rows x %2d columns | completeness min %.2f | %s\n", name, nrow(D), length(cols), min(md$completeness), paste(head(cols, 6), collapse = ", ")))
}

# ── climatology ──────────────────────────────────────────────────────────────
engineer_clim <- function(C) {
  P <- as.matrix(C[, sprintf("pr_m%02d", 1:12)]); TX <- as.matrix(C[, sprintf("tmax_m%02d", 1:12)])
  LD <- as.matrix(C[, sprintf("lstd_m%02d", 1:12)]); LN <- as.matrix(C[, sprintf("lstn_m%02d", 1:12)])
  ang <- 2 * pi * (1:12) / 12
  pr_ann <- rowSums(P); pr_bar <- pr_ann / 12
  si <- rowSums(abs(P - pr_bar)) / pmax(pr_ann, 1)                                   # Walsh-Lawler seasonality index
  a1 <- P %*% cos(ang) * 2 / 12; b1 <- P %*% sin(ang) * 2 / 12                        # first harmonic
  amp <- sqrt(a1^2 + b1^2) / pmax(pr_bar, 1); ph <- atan2(b1, a1)
  tmax_ann <- rowMeans(TX); tmax_hot <- apply(TX, 1, max); tmax_rng <- tmax_hot - apply(TX, 1, min)
  data.frame(C[, KEY],
    fe_pr_ann = pr_ann, fe_pr_si = as.numeric(si), fe_pr_harm_amp = as.numeric(amp), fe_pr_phase_sin = as.numeric(sin(ph)), fe_pr_phase_cos = as.numeric(cos(ph)),
    fe_pr_wet_months = rowSums(P > 100), fe_pr_dry_months = rowSums(P < 50),
    fe_pr_iacv = C$pr_ann_sd / pmax(C$pr_ann_mean, 1),
    fe_pr_anom_sy = (C$pr_sy - C$pr_ann_mean) / pmax(C$pr_ann_sd, 1), fe_pr_anom_win = (C$pr_win - C$pr_ann_mean) / pmax(C$pr_ann_sd, 1),
    fe_tmax_ann = tmax_ann, fe_tmin_ann = C$tmin_ann, fe_dtr = tmax_ann - C$tmin_ann, fe_tmax_hot = tmax_hot, fe_tmax_range = tmax_rng,
    fe_tmax_anom_sy = C$tmax_sy - tmax_ann,
    fe_pet_ann = C$pet_ann, fe_aridity = pr_ann / pmax(12 * C$pet_ann, 1), fe_def_ann = C$def_ann, fe_aet_ann = C$aet_ann,
    fe_soilm_ann = C$soilm_ann, fe_vpd_ann = C$vpd_ann, fe_srad_ann = C$srad_ann,
    fe_lst_day_ann = rowMeans(LD), fe_lst_night_ann = rowMeans(LN), fe_lst_day_range = apply(LD, 1, max) - apply(LD, 1, min),
    fe_lst_night_range = apply(LN, 1, max) - apply(LN, 1, min), fe_lst_diurnal = rowMeans(LD) - rowMeans(LN), fe_lst_day_hot = apply(LD, 1, max),
    check.names = FALSE)
}

# ── terrain ──────────────────────────────────────────────────────────────────
engineer_terrain <- function(T) data.frame(T[, KEY],
  fe_elev = T$elv, fe_hand = T$hnd, fe_log_upa = log1p(T$upa), fe_slope = T$slope, fe_tri = T$tri, fe_tpi = T$tpi,
  fe_roughness = T$roughness, fe_vrm = T$vrm, fe_cti = T$cti, fe_elev_sd90 = T$elev_stdev, check.names = FALSE)

# ── soil ─────────────────────────────────────────────────────────────────────
engineer_soil <- function(S) {
  d50 <- function(v) 0.4 * S[[paste0(v, "_0_20")]] + 0.6 * S[[paste0(v, "_20_50")]]    # depth-weighted 0-50 cm
  ph <- d50("ph"); clay <- d50("clay"); sand <- d50("sand"); silt <- d50("silt"); oc <- d50("oc"); cec <- d50("cec")
  zn <- d50("zn"); fe <- d50("fe"); ca <- d50("ca"); mg <- d50("mg"); k <- d50("k"); p <- d50("p"); s <- d50("s"); al <- d50("al"); n <- d50("ntot"); bd <- d50("bd")
  # availability: extractable zinc and iron fall away above neutral pH; phosphorus is fixed by Al / Fe below 5.5 and by Ca above 7.5
  zn_av <- log1p(zn) + log(sig(ph, 7.0)); fe_av <- log1p(fe) + log(sig(ph, 6.8)); p_av <- log1p(p) + log(sig(ph, 7.5) * (1 - sig(ph, 5.5)))
  # base-saturation proxy: exchangeable Ca + Mg + K (cmol(+)/kg from ppm) over CEC
  bases <- ca / 200 + mg / 120 + k / 390; bs <- pmin(bases / pmax(cec, 0.1), 1.5)
  fert <- (zc(log1p(oc)) + zc(log1p(cec)) + zc(log1p(n)) + zc(bs)) / 4
  data.frame(S[, KEY],
    fe_ph = ph, fe_clay = clay, fe_sand = sand, fe_silt = silt, fe_fines = clay + silt, fe_bd = bd, fe_oc = oc, fe_ntot = n, fe_cec = cec,
    fe_zn = zn, fe_fe = fe, fe_ca = ca, fe_mg = mg, fe_k = k, fe_p = p, fe_s = s, fe_al = al,
    fe_zn_avail = zn_av, fe_fe_avail = fe_av, fe_p_avail = p_av,
    fe_ca_mg = log((ca + 1) / (mg + 1)), fe_k_mg = log((k + 1) / (mg + 1)), fe_oc_clay = log1p(oc) - log1p(clay), fe_cn = log1p(oc) - log1p(n),
    fe_base_sat = bs, fe_fertility = fert, fe_ph_topsub = S$ph_0_20 - S$ph_20_50, fe_oc_topsub = log1p(S$oc_0_20) - log1p(S$oc_20_50),
    check.names = FALSE)
}

cat("[FE-01] engineered blocks\n")
C <- rd("gee_climatology_admin2.csv"); T <- rd("gee_terrain_admin2.csv"); S <- rd("gee_isda_admin2.csv")
for (suf in c("aw", "pw")) {
  write_block(engineer_clim(pick(C, suf)), paste0("climatology_", suf), "Climatology (engineered)", "TerraClimate 1991-2020, MODIS MOD11A2 2003-2020 (Earth Engine)", if (suf == "aw") "area-weighted zonal mean" else "WorldPop-weighted zonal mean at the survey year")
  write_block(engineer_terrain(pick(T, suf)), paste0("terrain_", suf), "Terrain relief (engineered)", "MERIT Hydro v1.0.1, Geomorpho90m (Earth Engine)", if (suf == "aw") "area-weighted zonal mean" else "WorldPop-weighted zonal mean at the survey year")
  write_block(engineer_soil(pick(S, suf)), paste0("soil_", suf), "Soil bioavailability (engineered)", "iSDAsoil Africa v1 (Earth Engine), back-transformed to natural units", if (suf == "aw") "area-weighted zonal mean" else "WorldPop-weighted zonal mean at the survey year")
}
cat("\nsanity (aw, medians over districts):\n")
Sa <- pick(S, "aw"); for (v in c("ph_0_20", "clay_0_20", "oc_0_20", "cec_0_20", "zn_0_20", "fe_0_20", "ca_0_20", "mg_0_20", "p_0_20", "bd_0_20")) cat(sprintf("  %-10s %.2f\n", v, median(Sa[[v]], na.rm = TRUE)))
Ca <- pick(C, "aw"); cat(sprintf("  pr_ann %.0f mm | tmax_ann %.1f C | lst_night %.1f C | pr_sy anomaly z %.2f (mean over districts)\n", median(rowSums(Ca[, sprintf("pr_m%02d", 1:12)])), median(rowMeans(Ca[, sprintf("tmax_m%02d", 1:12)])), median(rowMeans(Ca[, sprintf("lstn_m%02d", 1:12)])), mean((Ca$pr_sy - Ca$pr_ann_mean) / Ca$pr_ann_sd, na.rm = TRUE)))
cat("DONE\n")
