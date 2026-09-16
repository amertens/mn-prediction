# =============================================================================
# scripts/covariates/build_climate_normals_block.R   [CN-01, 2026-09-15]
#
# LONG-RUN CLIMATE NORMALS AND THE SURVEY-YEAR ANOMALY, FOR THE SHARED SET
#
# Retiring the TRMM year series (precip_y2010-2019) left the shared set with
# survey-year climate only (tclim_*_t0): no long-run rainfall, no inter-annual
# variability, no seasonality. Long-run climate is what transports (DA-01,
# DA-03), and the regional climate + soil arm fell 0.448 -> 0.418 when the
# series went (RR-11). This block restores it the way the NDVI switch did -
# "normal + anomaly" - from the 30-year climatology script 54 already
# extracted (data/covariates/harmonized/gee_climatology_admin2.csv: TerraClimate
# 1991-2020 monthly, MODIS LST 2003-2020; area-weighted _aw columns):
#
#   clim_pr_ann_mean       mean annual precipitation (mm), 1991-2020
#   clim_pr_ann_sd         s.d. of the 30 annual totals (inter-annual variability)
#   clim_pr_cv             sd / mean
#   clim_pr_season_cv      coefficient of variation of the 12 monthly normals
#   clim_pr_top3_share     share of annual rain in the wettest three months
#   clim_pr_sy_anom_z      (survey-year total - mean) / sd
#   clim_pr_win_anom_z     (12 months to the fieldwork median month - mean) / sd
#   clim_tmax_ann          mean of the 12 monthly tmax normals (C)
#   clim_tmax_range        warmest minus coolest monthly tmax normal
#   clim_tmax_sy_anom      survey-year mean tmax minus the normal
#   clim_tmin_ann, clim_pet_ann, clim_def_ann, clim_aet_ann, clim_soilm_ann,
#   clim_vpd_ann, clim_srad_ann   annual normals (TerraClimate units)
#   clim_lstd_ann, clim_lstn_ann  mean monthly day / night LST (C)
#   clim_lst_diurnal       lstd - lstn
#
#   Rscript -e "source('scripts/covariates/build_climate_normals_block.R')"
# -> data/covariates/harmonized/predictors_admin2_clim.csv (+ _metadata.csv)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"
SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
spine <- SH[, c("country", "Admin1", "Admin2")]
C <- read.csv(file.path(HDIR, "gee_climatology_admin2.csv"), check.names = FALSE, stringsAsFactors = FALSE)
aw <- function(v) C[[paste0(v, "_aw")]]
mon <- function(p) sapply(sprintf("%s_m%02d", p, 1:12), aw)   # 554 x 12
prm <- mon("pr"); txm <- mon("tmax"); ldm <- mon("lstd"); lnm <- mon("lstn")
top3 <- t(apply(prm, 1, function(x) { s <- sort(x, decreasing = TRUE); c(sum(s[1:3]), sum(x)) }))
OUT <- data.frame(country = C$country, Admin1 = C$Admin1, Admin2 = C$Admin2,
  clim_pr_ann_mean = aw("pr_ann_mean"), clim_pr_ann_sd = aw("pr_ann_sd"), clim_pr_cv = aw("pr_ann_sd") / aw("pr_ann_mean"),
  clim_pr_season_cv = apply(prm, 1, sd) / rowMeans(prm), clim_pr_top3_share = top3[, 1] / top3[, 2],
  clim_pr_sy_anom_z = (aw("pr_sy") - aw("pr_ann_mean")) / aw("pr_ann_sd"), clim_pr_win_anom_z = (aw("pr_win") - aw("pr_ann_mean")) / aw("pr_ann_sd"),
  clim_tmax_ann = rowMeans(txm), clim_tmax_range = apply(txm, 1, max) - apply(txm, 1, min), clim_tmax_sy_anom = aw("tmax_sy") - rowMeans(txm),
  clim_tmin_ann = aw("tmin_ann"), clim_pet_ann = aw("pet_ann"), clim_def_ann = aw("def_ann"), clim_aet_ann = aw("aet_ann"),
  clim_soilm_ann = aw("soilm_ann"), clim_vpd_ann = aw("vpd_ann"), clim_srad_ann = aw("srad_ann"),
  clim_lstd_ann = rowMeans(ldm), clim_lstn_ann = rowMeans(lnm), clim_lst_diurnal = rowMeans(ldm) - rowMeans(lnm), stringsAsFactors = FALSE)
OUT <- left_join(spine, OUT, by = c("country", "Admin1", "Admin2"))
stopifnot(nrow(OUT) == nrow(spine))
cols <- grep("^clim_", names(OUT), value = TRUE)
for (v in cols) OUT[[v]][!is.finite(OUT[[v]])] <- NA
write.csv(OUT, file.path(HDIR, "predictors_admin2_clim.csv"), row.names = FALSE)

desc <- c(clim_pr_ann_mean = "Mean annual precipitation 1991-2020 (mm)", clim_pr_ann_sd = "SD of the 30 annual precipitation totals (mm)", clim_pr_cv = "Inter-annual coefficient of variation of annual precipitation",
          clim_pr_season_cv = "Coefficient of variation of the 12 monthly precipitation normals (seasonality)", clim_pr_top3_share = "Share of mean annual precipitation falling in the three wettest months",
          clim_pr_sy_anom_z = "Survey-year precipitation anomaly in SD units of the 1991-2020 annual totals", clim_pr_win_anom_z = "Precipitation over the 12 months to the fieldwork median month, anomaly in SD units",
          clim_tmax_ann = "Mean of the monthly maximum-temperature normals (C)", clim_tmax_range = "Warmest minus coolest monthly tmax normal (C)", clim_tmax_sy_anom = "Survey-year mean tmax minus its normal (C)",
          clim_tmin_ann = "Annual mean of monthly minimum temperature (C)", clim_pet_ann = "Annual mean potential evapotranspiration (mm/month)", clim_def_ann = "Annual mean climatic water deficit (mm/month)",
          clim_aet_ann = "Annual mean actual evapotranspiration (mm/month)", clim_soilm_ann = "Annual mean soil moisture (mm)", clim_vpd_ann = "Annual mean vapour-pressure deficit (kPa)", clim_srad_ann = "Annual mean downward shortwave radiation (W/m2)",
          clim_lstd_ann = "Mean monthly daytime land-surface temperature 2003-2020 (C)", clim_lstn_ann = "Mean monthly night-time land-surface temperature 2003-2020 (C)", clim_lst_diurnal = "Day minus night LST (C)")
md <- data.frame(column = cols, source = "GEE", domain = "Climate and weather", subnational = TRUE,
                 assumption = paste(unname(desc[cols]), "TerraClimate 1991-2020 monthly and MODIS LST 2003-2020 climatologies, area-weighted over the GADM Admin-2 polygon (scripts/protocol_v2/54_extract_climatology_terrain_soil.py); normals are time-invariant, the anomalies are at the survey year (CN-01)."), stringsAsFactors = FALSE)
md$n_countries <- vapply(md$column, function(v) sum(tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
md$countries <- vapply(md$column, function(v) { s <- tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))); paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
md$completeness <- round(vapply(md$column, function(v) mean(is.finite(OUT[[v]])), 0), 3)
write.csv(md, file.path(HDIR, "predictors_admin2_clim_metadata.csv"), row.names = FALSE)
cat(sprintf("climate normals block: %d columns, completeness %s\n", length(cols), paste(range(md$completeness), collapse = "-")))
print(OUT |> group_by(country) |> summarise(across(c(clim_pr_ann_mean, clim_pr_cv, clim_pr_top3_share, clim_pr_sy_anom_z, clim_tmax_ann, clim_lst_diurnal), ~ round(mean(.x, na.rm = TRUE), 2))), width = 200)
cat("DONE\n")
