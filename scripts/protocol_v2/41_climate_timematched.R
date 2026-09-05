# =============================================================================
# scripts/protocol_v2/41_climate_timematched.R   [SC-01]
#
# CLIMATE MATCHED TO THE FIELDWORK WINDOW, AND SEASONALITY THAT TRANSPORTS
#
# Night-time land-surface temperature carries women's iron, and every
# climate column in the vocabulary is an annual summary. Ferritin and retinol
# respond to infection, which follows rainfall and temperature with a lag of
# weeks to months, so the climate of the months around the blood draw is a
# different quantity from the annual mean. Two kinds of column are built:
#
#   TIME-MATCHED (need the survey's fieldwork dates; nuisance adjusters for
#   the training targets, or predictors only when a future survey's timing
#   is known):
#     tmc_lst_night_win        night LST, mean over the fieldwork months
#     tmc_lst_night_anom       the same minus the survey-year annual mean
#     tmc_lst_night_lag1/2     mean over the window shifted back 1 / 2 months
#     tmc_precip_win           CHIRPS rainfall, mean over the fieldwork months
#     tmc_precip_anom          minus the survey-year monthly mean
#     tmc_precip_prev3         rainfall in the 3 months before the window
#                              (transmission lag)
#     tmc_wet_share            share of window months wetter than the year's
#                              monthly mean
#     tmc_fieldwork_month      the calendar month itself
#   SEASONALITY (fieldwork-independent, transportable):
#     sea_precip_amp           max - min monthly rainfall of the survey year
#     sea_precip_cv            CV of monthly rainfall
#     sea_precip_peak_month    month of peak rainfall
#     sea_lst_night_peak_month month of peak night LST
#     sea_precip_conc          rainfall concentration: share in the wettest 3
#                              months
#
# Sources: lst_night_m01..m12_t0 already in the shared table (survey year);
# CHIRPS monthly rasters cached under data/external_cache/chirps/<ISO>/ for
# the survey year, zonal means over the project's Admin-2 polygons. Malawi's
# fieldwork spans Dec 2015 - Feb 2016 and only 2016 is cached, so December
# uses 2016's December (flagged in the metadata).
#
#   Rscript scripts/protocol_v2/41_climate_timematched.R
# -> data/covariates/harmonized/predictors_admin2_climate_tm.csv (+ _metadata.csv)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"; HDIR <- "data/covariates/harmonized"
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
LC  <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi", SierraLeone = "sierraleone")
S  <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
FW <- read.csv(file.path(OUTDIR, "fieldwork_windows_admin2.csv"), stringsAsFactors = FALSE)
FW$date_first <- as.Date(FW$date_first); FW$date_last <- as.Date(FW$date_last)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
lst_cols <- sprintf("lst_night_m%02d_t0", 1:12); stopifnot(all(lst_cols %in% names(S)))
mon_of <- function(d) as.integer(format(d, "%m"))
win_months <- function(d1, d2) { a <- mon_of(d1); b <- mon_of(d2); n <- as.integer(round(as.numeric(d2 - d1) / 30.4)) + 1L; ((a - 1 + seq_len(max(1L, n)) - 1L) %% 12) + 1L }
shift <- function(ms, k) ((ms - 1 - k) %% 12) + 1

rows <- list()
for (cn in names(ISO)) {
  b <- BND[[LC[[cn]]]]; keys <- data.frame(country = cn, Admin1 = as.character(sf::st_drop_geometry(b)$Admin1), Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), stringsAsFactors = FALSE)
  # CHIRPS monthly zonal means for the survey year
  # layout written by R/external_data.R: chirps/<ISO>/<year>/chirps_<year>_<MM>.tif
  dir <- file.path("data/external_cache/chirps", ISO[[cn]]); fs <- list.files(dir, pattern = "^chirps_[0-9]{4}_[0-9]{2}\\.tif$", full.names = TRUE, recursive = TRUE)
  yrs <- unique(substr(basename(fs), 8, 11)); cat(cn, "CHIRPS years on disk:", paste(yrs, collapse = " "), "\n")
  P <- matrix(NA_real_, nrow = nrow(keys), ncol = 12)
  if (length(fs)) { yr <- yrs[1]
    v <- terra::vect(sf::st_transform(b, 4326))
    for (mth in 1:12) { f <- file.path(dir, yr, sprintf("chirps_%s_%02d.tif", yr, mth)); if (!file.exists(f)) next
      r <- terra::rast(f); r[r < 0] <- NA
      P[, mth] <- tryCatch(terra::extract(r, v, fun = mean, na.rm = TRUE, ID = FALSE)[, 1], error = function(e) rep(NA_real_, nrow(keys))) } }
  sx <- S[S$country == cn, c("Admin1", "Admin2", lst_cols)]; j <- match(paste(keys$Admin1, keys$Admin2), paste(sx$Admin1, sx$Admin2))
  L <- as.matrix(sx[j, lst_cols]); L[!is.finite(L)] <- NA
  fw <- FW[FW$country == cn, ]; k <- match(paste(keys$Admin1, keys$Admin2), paste(fw$Admin1, fw$Admin2))
  out <- keys
  for (i in seq_len(nrow(keys))) {
    ms <- if (is.na(k[i])) integer(0) else win_months(fw$date_first[k[i]], fw$date_last[k[i]])
    lm_ <- L[i, ]; pm_ <- P[i, ]
    g <- function(x, idx) if (!length(idx) || all(is.na(x[idx]))) NA_real_ else mean(x[idx], na.rm = TRUE)
    out$tmc_lst_night_win[i]  <- g(lm_, ms); out$tmc_lst_night_anom[i] <- g(lm_, ms) - mean(lm_, na.rm = TRUE)
    out$tmc_lst_night_lag1[i] <- g(lm_, shift(ms, 1)); out$tmc_lst_night_lag2[i] <- g(lm_, shift(ms, 2))
    out$tmc_precip_win[i] <- g(pm_, ms); out$tmc_precip_anom[i] <- g(pm_, ms) - mean(pm_, na.rm = TRUE)
    prev3 <- if (length(ms)) unique(shift(rep(ms[1], 3), 1:3)) else integer(0); out$tmc_precip_prev3[i] <- if (length(prev3)) sum(pm_[prev3], na.rm = TRUE) else NA_real_
    out$tmc_wet_share[i] <- if (length(ms) && any(is.finite(pm_))) mean(pm_[ms] > mean(pm_, na.rm = TRUE), na.rm = TRUE) else NA_real_
    out$tmc_fieldwork_month[i] <- if (is.na(k[i])) NA_real_ else fw$month_med[k[i]]
    out$sea_precip_amp[i] <- if (any(is.finite(pm_))) max(pm_, na.rm = TRUE) - min(pm_, na.rm = TRUE) else NA_real_
    out$sea_precip_cv[i] <- if (any(is.finite(pm_)) && mean(pm_, na.rm = TRUE) > 0) stats::sd(pm_, na.rm = TRUE) / mean(pm_, na.rm = TRUE) else NA_real_
    out$sea_precip_peak_month[i] <- if (any(is.finite(pm_))) which.max(pm_) else NA_real_
    out$sea_precip_conc[i] <- if (any(is.finite(pm_)) && sum(pm_, na.rm = TRUE) > 0) sum(sort(pm_, decreasing = TRUE)[1:3], na.rm = TRUE) / sum(pm_, na.rm = TRUE) else NA_real_
    out$sea_lst_night_peak_month[i] <- if (any(is.finite(lm_))) which.max(lm_) else NA_real_
  }
  rows[[cn]] <- out
  cat(sprintf("%-12s units %3d | dated %3d | CHIRPS finite %3d | LST finite %3d\n", cn, nrow(out), sum(!is.na(k)), sum(is.finite(out$sea_precip_amp)), sum(is.finite(out$tmc_lst_night_win))))
}
CT <- bind_rows(rows); newcols <- setdiff(names(CT), c("country", "Admin1", "Admin2"))
CT[newcols] <- lapply(CT[newcols], function(z) round(z, 4))
write.csv(CT, file.path(HDIR, "predictors_admin2_climate_tm.csv"), row.names = FALSE)
MDT <- data.frame(column = newcols, domain = ifelse(grepl("^sea_", newcols), "Seasonality (transportable)", "TM climate"),
  source = "lst_night_m01-12_t0 (shared) + CHIRPS monthly (external_cache) x FW-01 fieldwork windows",
  completeness = round(vapply(newcols, function(cl) mean(is.finite(CT[[cl]])), 0), 3), subnational = TRUE,
  note = ifelse(grepl("^tmc_", newcols), "needs fieldwork dates; Malawi December uses the 2016 raster", "fieldwork-independent"), stringsAsFactors = FALSE)
write.csv(MDT, file.path(HDIR, "predictors_admin2_climate_tm_metadata.csv"), row.names = FALSE)
cat("\n===== SC-01: time-matched climate add-on =====\n"); print(MDT[, c("column", "domain", "completeness")], row.names = FALSE)
print(as.data.frame(CT |> group_by(country) |> summarise(across(c(tmc_lst_night_anom, tmc_precip_anom, tmc_precip_prev3, tmc_wet_share, sea_precip_amp, sea_precip_conc), ~ round(stats::median(.x, na.rm = TRUE), 3)), .groups = "drop")), row.names = FALSE)
cat("\nDONE\n")
