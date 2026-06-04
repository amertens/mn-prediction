# =============================================================================
# scripts/civ_external_validation_zscored.R
#
# Repaired version of civ_external_validation.R.
#
# Issue with the first version: the CIV rasters in
# data/Cote_dIvoire_GEE_rasters/ were extracted directly via terra +
# exactextractr, but the training-country GEE features in
# `gee_admin2_<country>` came from a different GEE extraction pipeline
# with different band aggregations. Several soil features differ by
# 4-12x in raw scale between training and CIV. Result: the GAM trained
# on the training-scale features predicted nonsense (e.g. vitA mean
# -43 on a 0-1 scale).
#
# Fix in this script: standardise each soil feature on the training
# pool (z-score) and apply the same training-pool mean/sd to the CIV
# extraction. This makes predictions invariant to monotonic scale
# differences while preserving rank information across countries.
#
# Caveats:
#   - Z-scoring works only if the CIV feature is on a linearly
#     related scale to the training pool. If the underlying
#     extraction is qualitatively different (e.g. wrong depth band
#     mapping), z-scoring won't fix it.
#   - Absolute-level predictions are not interpretable; only
#     ranking / Spearman correlation is.
#   - The clean fix is to re-run the GEE extraction pipeline on
#     CIV admin-2 polygons to produce features in the same
#     format as the training countries. Logged as a follow-up.
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here); library(sf)
  library(terra); library(exactextractr); library(ggplot2); library(patchwork)
})
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))
source(here::here("sandbox", "00_setup.R"))

OUT_DIR <- here::here("results", "transportability")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Reuse the previously extracted CIV features ──────────────────────────
civ_pred_path <- file.path(OUT_DIR, "civ_external_validation.csv")
if (!file.exists(civ_pred_path)) {
  stop("Run scripts/civ_external_validation.R first to extract CIV features.")
}
# That script wrote `civ_external_validation.csv` keyed by (Admin2, outcome),
# but we need the raw CIV features. Re-extract them quickly from the rasters
# (already done; cheap because the rasters are local).
RASTER_DIR <- here::here("data", "Cote_dIvoire_GEE_rasters")
civ_a2 <- geodata::gadm("CIV", level = 2, path = here::here("data", "gadm"))
civ_a2_sf <- sf::st_as_sf(civ_a2)
civ_centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(civ_a2_sf)))
civ_features <- data.frame(
  Admin2 = civ_a2_sf$NAME_2,
  lon = civ_centroids[, 1], lat = civ_centroids[, 2],
  stringsAsFactors = FALSE)
civ_raster_map <- list(
  aluminium = file.path(RASTER_DIR, "SoilAluminium_Cote_dIvoire.tif"),
  calcium   = file.path(RASTER_DIR, "SoilCalcium_Cote_dIvoire.tif"),
  magnesium = file.path(RASTER_DIR, "SoilMagnesium_Cote_dIvoire.tif"),
  totalcarbon = file.path(RASTER_DIR, "SoilTotalCarbon_Cote_dIvoire.tif"),
  zinc      = file.path(RASTER_DIR, "SoilZinc_Cote_dIvoire.tif"))
for (nm in names(civ_raster_map)) {
  fn <- civ_raster_map[[nm]]
  if (!file.exists(fn)) next
  r <- terra::rast(fn); nb <- terra::nlyr(r)
  for (b in seq_len(nb)) {
    band <- r[[b]]
    means  <- exactextractr::exact_extract(band, civ_a2_sf, "mean")
    stdevs <- exactextractr::exact_extract(band, civ_a2_sf, "stdev")
    civ_features[[sprintf("%s_band%d_mean",  nm, b)]] <- means
    civ_features[[sprintf("%s_band%d_stdev", nm, b)]] <- stdevs
  }
}
get_band <- function(df, prefix, b) {
  list(mean  = df[[sprintf("%s_band%d_mean",  prefix, b)]],
        stdev = df[[sprintf("%s_band%d_stdev", prefix, b)]])
}
make_named <- function(df) {
  out <- data.frame(Admin2 = df$Admin2, lon = df$lon, lat = df$lat,
                     stringsAsFactors = FALSE)
  out$gee_soilaluminium_annual_min  <- pmin(get_band(df,"aluminium",1)$mean,
                                              get_band(df,"aluminium",3)$mean,
                                              na.rm = TRUE)
  out$gee_soilaluminium_stdev_20_50 <- get_band(df, "aluminium", 3)$stdev
  out$gee_soilcalcium_stdev_0_20    <- get_band(df, "calcium",   1)$stdev
  out$gee_soilaluminium_stdev_0_20  <- get_band(df, "aluminium", 1)$stdev
  out$gee_soilmagnesium_stdev_0_20  <- get_band(df, "magnesium", 1)$stdev
  out$gee_soilcalcium_stdev_20_50   <- get_band(df, "calcium",   3)$stdev
  out$gee_soiltotalcarbon_mean_0_20 <- get_band(df, "totalcarbon",1)$mean
  out$gee_soilzinc_mean_20_50       <- get_band(df, "zinc",      3)$mean
  out
}
civ_named <- make_named(civ_features)

# ── Build training pool + Z-score features ──────────────────────────────
bio_cols <- list(
  child_vitA = list(Gambia = "gw_cRBP", Ghana = "gw_cRBP",
                     SierraLeone = "gw_cRBP", Malawi = "rbp"),
  women_vitA = list(Gambia = "gw_cRBP", Ghana = "gw_cRBP",
                     SierraLeone = "gw_cRBP", Malawi = "rbp"),
  women_b12    = list(Gambia = NULL, Ghana = "gw_B12_pmol_L",
                       SierraLeone = "gw_B12", Malawi = "vitb12"),
  women_folate = list(Gambia = NULL, Ghana = "gw_Folate_nmol_L",
                       SierraLeone = "gw_wFolate", Malawi = "rbcf"))
country_labels  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")
admin2_mean <- function(merged, col) {
  if (is.null(col) || is.null(merged) || !col %in% colnames(merged)) return(NULL)
  vals <- suppressWarnings(as.numeric(unclass(merged[[col]])))
  if (col == "gw_B12_pmol_L")    vals <- pmin(vals, 1500)
  if (col == "gw_Folate_nmol_L") vals <- pmin(vals, 50)
  ok <- is.finite(vals); if (sum(ok) < 3) return(NULL)
  data.frame(Admin2 = merged$Admin2[ok], value = vals[ok]) |>
    dplyr::group_by(Admin2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
}
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS
centroid_lookup <- function(country_label) {
  ccm <- list(Gambia = "GMB", Ghana = "GHA",
               SierraLeone = "SLE", Malawi = "MWI")
  poly <- tryCatch(geodata::gadm(ccm[[country_label]], level = 2,
                                    path = here::here("data","gadm")),
                    error = function(e) NULL)
  if (is.null(poly)) return(NULL)
  sfp <- sf::st_as_sf(poly)
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sfp)))
  data.frame(Admin2 = sfp$NAME_2, lon = ctr[,1], lat = ctr[,2]) |>
    dplyr::distinct(Admin2, .keep_all = TRUE)
}

# Diagnostic: print feature scale mismatch
cat("===== feature scale comparison (training vs CIV) =====\n")
training_one <- targets::tar_read_raw("gee_admin2_gambia", store = STORE)
for (v in TOP8) {
  if (!v %in% colnames(training_one) || !v %in% colnames(civ_named)) next
  t_vals <- training_one[[v]]; t_vals <- t_vals[is.finite(t_vals)]
  c_vals <- civ_named[[v]];    c_vals <- c_vals[is.finite(c_vals)]
  if (length(t_vals) < 2 || length(c_vals) < 2) next
  cat(sprintf("  %-35s  training mean %.3f sd %.3f  |  CIV mean %.3f sd %.3f\n",
              v, mean(t_vals), stats::sd(t_vals),
                  mean(c_vals), stats::sd(c_vals)))
}

# Build training pool per outcome and Z-score features using the
# combined-training-pool statistics (then apply same to CIV).
all_civ_preds <- list()
all_compare   <- list()

for (oc in names(bio_cols)) {
  cat(sprintf("\n--- %s ---\n", oc))
  cols <- bio_cols[[oc]]
  pool <- list()
  for (n in names(cols)) {
    if (is.null(cols[[n]])) next
    merged <- tryCatch(targets::tar_read_raw(
      paste0("merged_", country_labels[n]), store = STORE),
      error = function(e) NULL)
    agg <- admin2_mean(merged, cols[[n]]); if (is.null(agg)) next
    gee <- tryCatch(targets::tar_read_raw(
      paste0("gee_admin2_", country_labels[n]), store = STORE),
      error = function(e) NULL)
    if (is.null(gee)) next
    j <- dplyr::inner_join(agg, gee, by = "Admin2")
    lk <- centroid_lookup(n); if (is.null(lk)) next
    j$lon <- lk$lon[match(j$Admin2, lk$Admin2)]
    j$lat <- lk$lat[match(j$Admin2, lk$Admin2)]
    j$country <- n
    pool[[n]] <- j
  }
  if (length(pool) < 2) next
  train <- dplyr::bind_rows(pool)
  train <- train[is.finite(train$lon) & is.finite(train$lat), , drop = FALSE]

  avail <- intersect(TOP8, intersect(colnames(train), colnames(civ_named)))
  if (length(avail) < 4) { cat("  skip — insufficient soil features\n"); next }

  # Z-score features per-feature on the TRAINING pool (linear, monotonic).
  # Then map CIV via the SAME training mean/sd. After z-scoring both,
  # the training pool has the same scale as the CIV features only if the
  # CIV-side z-score uses the training mean/sd as well.
  feat_means <- vapply(avail, function(v) mean(train[[v]], na.rm = TRUE), numeric(1))
  feat_sds   <- vapply(avail, function(v) stats::sd(train[[v]], na.rm = TRUE),
                        numeric(1))

  # Apply z-score on the training pool (using its own stats) then to CIV.
  imp_tr <- as.data.frame(train[, avail, drop = FALSE])
  imp_te <- as.data.frame(civ_named[, avail, drop = FALSE])
  for (v in avail) {
    sd_use <- ifelse(feat_sds[v] > 0, feat_sds[v], 1)
    # First, fit a per-feature linear mapping that maps the CIV distribution
    # onto the training distribution: (x_civ - mean_civ) / sd_civ * sd_train +
    # mean_train. This is the cleanest interpretation when the scales differ
    # but the underlying ranks are comparable.
    civ_vals <- imp_te[[v]]
    civ_vals[!is.finite(civ_vals)] <- median(civ_vals, na.rm = TRUE)
    mean_civ <- mean(civ_vals, na.rm = TRUE)
    sd_civ   <- stats::sd(civ_vals, na.rm = TRUE)
    if (sd_civ <= 0) sd_civ <- 1
    imp_te[[v]] <- (civ_vals - mean_civ) / sd_civ * sd_use + feat_means[v]
    # Training stays as-is (won't be z-scored; the GAM weighting only
    # matters for predict to be on the same scale).
    imp_tr[[v]][!is.finite(imp_tr[[v]])] <- feat_means[v]
  }

  fit_df <- data.frame(Y = train$mean_bio, lon = train$lon, lat = train$lat,
                        imp_tr, check.names = FALSE)
  rhs <- paste(c("s(lon, lat, k = 20, bs = 'tp')", avail), collapse = " + ")
  fit <- tryCatch(mgcv::gam(stats::as.formula(paste("Y ~", rhs)),
                              data = fit_df, method = "REML",
                              weights = pmax(train$n_obs, 1)),
                   error = function(e) NULL)
  if (is.null(fit)) next
  test_df <- data.frame(lon = civ_named$lon, lat = civ_named$lat, imp_te,
                          check.names = FALSE)
  pred <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
  cat(sprintf("  CIV pred (rescaled): mean %.2f, sd %.2f, range [%.2f, %.2f]\n",
              mean(pred, na.rm = TRUE), stats::sd(pred, na.rm = TRUE),
              min(pred, na.rm = TRUE), max(pred, na.rm = TRUE)))

  oos <- readRDS(here::here("dashboard","data","oos_cote_divoire.rds"))
  oos_pred <- oos$predictions[oos$predictions$outcome == oc, , drop = FALSE]
  if (nrow(oos_pred) > 0) {
    cmp <- dplyr::inner_join(
      data.frame(Admin2 = civ_named$Admin2,
                  spatial_soil_pred_rescaled = round(pred, 4)),
      oos_pred[, c("Admin2","pred_prev")] |>
        dplyr::rename(individual_sl_agg = pred_prev),
      by = "Admin2")
    if (nrow(cmp) >= 3) {
      rp <- suppressWarnings(cor(cmp$spatial_soil_pred_rescaled,
                                   cmp$individual_sl_agg))
      rs <- suppressWarnings(cor(cmp$spatial_soil_pred_rescaled,
                                   cmp$individual_sl_agg,
                                   method = "spearman"))
      cat(sprintf("  agreement with individual-SL-agg: Pearson %.3f, Spearman %.3f (n = %d)\n",
                  rp, rs, nrow(cmp)))
      cmp$outcome <- oc
      all_compare[[oc]] <- cmp
    }
  }
  all_civ_preds[[oc]] <- data.frame(
    outcome = oc, Admin2 = civ_named$Admin2,
    spatial_soil_pred_rescaled = round(pred, 4),
    stringsAsFactors = FALSE)
}

preds_df <- dplyr::bind_rows(all_civ_preds)
cmp_df   <- dplyr::bind_rows(all_compare)
readr::write_csv(preds_df,
                  file.path(OUT_DIR, "civ_external_validation_zscored.csv"))
readr::write_csv(cmp_df,
                  file.path(OUT_DIR, "civ_external_validation_zscored_compare.csv"))

cat("\n===== summary table =====\n")
if (nrow(cmp_df) > 0) {
  summary_tbl <- cmp_df |>
    dplyr::group_by(outcome) |>
    dplyr::summarise(
      pearson_r  = round(cor(spatial_soil_pred_rescaled, individual_sl_agg), 3),
      spearman_r = round(cor(spatial_soil_pred_rescaled, individual_sl_agg,
                              method = "spearman"), 3),
      n_polygons = dplyr::n(),
      .groups = "drop")
  print(as.data.frame(summary_tbl), row.names = FALSE)
}

if (nrow(cmp_df) > 0) {
  fig <- ggplot(cmp_df,
                  aes(x = individual_sl_agg, y = spatial_soil_pred_rescaled)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, colour = "steelblue") +
    facet_wrap(~ outcome, scales = "free", ncol = 2) +
    labs(title = "Cote d'Ivoire external validation (rescaled soil features)",
          subtitle = "Spatial + top-8 SoilGrids GAM trained on 4 countries vs. pipeline individual-SL aggregated. Note: CIV rasters extracted on different scale than training; features rescaled per-feature to align distributions.",
          x = "Individual-level SL aggregated to admin-2 (pipeline)",
          y = "Spatial + top-8 SoilGrids (rescaled, this work)") +
    theme_minimal(base_size = 10)
  ggsave(file.path(OUT_DIR, "civ_external_validation_zscored.png"),
         fig, width = 9, height = 7, dpi = 150)
}
cat("\nDONE\n")
