# =============================================================================
# scripts/civ_external_validation.R
#
# External validation of the spatial + top-8 SoilGrids recipe on
# Côte d'Ivoire — a fifth West African country not in the training set.
#
# Procedure:
#   1. Load the 33 CIV admin-2 polygons from GADM and compute centroids.
#   2. Extract the eight SoilGrids covariates at admin-2 resolution from
#      the rasters in data/Cote_dIvoire_GEE_rasters/ using exactextractr.
#   3. For each of the 4 continuous biomarkers (child_vitA, women_vitA,
#      women_b12, women_folate), train the spatial+top-8-soil GAM on the
#      4-country pool (Gambia, Ghana, Sierra Leone, Malawi) and predict
#      at CIV admin-2 polygons.
#   4. Compare predicted prevalences to (a) the existing
#      oos_cote_divoire.rds pipeline output (which used individual-level
#      SL aggregated to admin-2), and (b) any VMNIS national-level
#      estimates available for CIV.
#
# Output:
#   results/transportability/civ_external_validation.csv
#   results/transportability/civ_external_validation.png
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

RASTER_DIR <- here::here("data", "Cote_dIvoire_GEE_rasters")

# ── 1. CIV admin-2 polygons ─────────────────────────────────────────────
cat("===== STEP 1: CIV admin-2 polygons + centroids =====\n")
civ_a2 <- geodata::gadm("CIV", level = 2, path = here::here("data", "gadm"))
civ_a2_sf <- sf::st_as_sf(civ_a2)
cat(sprintf("CIV admin-2: %d polygons\n", nrow(civ_a2_sf)))
civ_centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(civ_a2_sf)))
civ_centroids_df <- data.frame(
  Admin1 = civ_a2_sf$NAME_1, Admin2 = civ_a2_sf$NAME_2,
  lon = civ_centroids[, 1], lat = civ_centroids[, 2])

# ── 2. Extract SoilGrids covariates at CIV admin-2 ──────────────────────
cat("\n===== STEP 2: SoilGrids extraction at CIV admin-2 =====\n")
# Map each top-8 covariate (in our trained model) to its CIV raster file.
# The raster files contain SoilGrids data at multiple depths in band layers.
# The top-8 from our trained model:
#   gee_soilaluminium_annual_min
#   gee_soilaluminium_stdev_20_50
#   gee_soilcalcium_stdev_0_20
#   gee_soilaluminium_stdev_0_20
#   gee_soilmagnesium_stdev_0_20
#   gee_soilcalcium_stdev_20_50
#   gee_soiltotalcarbon_mean_0_20
#   gee_soilzinc_mean_20_50
# CIV rasters have one element per: Aluminium, Calcium, Magnesium,
# TotalCarbon, Zinc. We'll compute mean and stdev at admin-2 from
# the available bands and produce the named features.

civ_raster_map <- list(
  aluminium = file.path(RASTER_DIR, "SoilAluminium_Cote_dIvoire.tif"),
  calcium   = file.path(RASTER_DIR, "SoilCalcium_Cote_dIvoire.tif"),
  magnesium = file.path(RASTER_DIR, "SoilMagnesium_Cote_dIvoire.tif"),
  totalcarbon = file.path(RASTER_DIR, "SoilTotalCarbon_Cote_dIvoire.tif"),
  zinc      = file.path(RASTER_DIR, "SoilZinc_Cote_dIvoire.tif")
)

# Extract mean and stdev per band per polygon.
civ_features <- civ_centroids_df
for (nm in names(civ_raster_map)) {
  fn <- civ_raster_map[[nm]]
  if (!file.exists(fn)) next
  r <- terra::rast(fn)
  nb <- terra::nlyr(r)
  cat(sprintf("  %s (%d bands): extracting mean + stdev per band\n", nm, nb))
  for (b in seq_len(nb)) {
    band <- r[[b]]
    means  <- exactextractr::exact_extract(band, civ_a2_sf, "mean")
    stdevs <- exactextractr::exact_extract(band, civ_a2_sf, "stdev")
    bandname <- sprintf("%s_band%d", nm, b)
    civ_features[[paste0(bandname, "_mean")]]  <- means
    civ_features[[paste0(bandname, "_stdev")]] <- stdevs
  }
}
cat(sprintf("\nCIV features assembled: %d polygons x %d cols\n",
            nrow(civ_features), ncol(civ_features)))

# Construct named CIV covariate columns matching trained-model names.
# SoilGrids typically has 4-6 depth bands per file (0-5, 5-15, 15-30,
# 30-60 cm or similar). We approximate "0_20" ~ band 1 mean of 0-5/5-15,
# and "20_50" ~ band 3 mean of 15-30/30-60. To stay simple and reproducible,
# we use band 1 as the "0_20" proxy and band 3 as the "20_50" proxy.
# `annual_min` is approximated by min across years/bands.
get_band <- function(df, prefix, b) {
  col_m <- sprintf("%s_band%d_mean",  prefix, b)
  col_s <- sprintf("%s_band%d_stdev", prefix, b)
  list(mean = df[[col_m]], stdev = df[[col_s]])
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
cat("\nCIV final feature columns:\n")
print(colnames(civ_named))
cat("\nfeature ranges in CIV:\n")
for (v in setdiff(colnames(civ_named), c("Admin2","lon","lat"))) {
  cat(sprintf("  %-35s: [%.3f, %.3f]  (mean %.3f)\n",
              v, min(civ_named[[v]], na.rm=TRUE),
                  max(civ_named[[v]], na.rm=TRUE),
                  mean(civ_named[[v]], na.rm=TRUE)))
}

# ── 3. For each outcome, train on 4-country pool and predict CIV ────────
cat("\n===== STEP 3: train spatial+soil on 4 countries, predict CIV =====\n")
# Country biomarker columns (from earlier audit) - continuous mean per admin2.
bio_cols <- list(
  child_vitA = list(Gambia = "gw_cRBP", Ghana = "gw_cRBP",
                     SierraLeone = "gw_cRBP", Malawi = "rbp"),
  women_vitA = list(Gambia = "gw_cRBP", Ghana = "gw_cRBP",
                     SierraLeone = "gw_cRBP", Malawi = "rbp"),
  women_b12    = list(Gambia = NULL, Ghana = "gw_B12_pmol_L",
                       SierraLeone = "gw_B12", Malawi = "vitb12"),
  women_folate = list(Gambia = NULL, Ghana = "gw_Folate_nmol_L",
                       SierraLeone = "gw_wFolate", Malawi = "rbcf")
)
country_labels  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")

admin2_mean <- function(merged, col) {
  if (is.null(col) || is.null(merged) || !col %in% colnames(merged)) return(NULL)
  vals <- suppressWarnings(as.numeric(unclass(merged[[col]])))
  # Apply data-quality caps from sandbox/20
  if (col == "gw_B12_pmol_L")    vals <- pmin(vals, 1500)
  if (col == "gw_Folate_nmol_L") vals <- pmin(vals, 50)
  ok <- is.finite(vals)
  if (sum(ok) < 3) return(NULL)
  data.frame(Admin2 = merged$Admin2[ok], value = vals[ok]) |>
    dplyr::group_by(Admin2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
}
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS
cc_list <- get_country_configs()

# Per-country centroid lookup (fresh, dedup)
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
    agg <- admin2_mean(merged, cols[[n]])
    if (is.null(agg)) next
    gee <- tryCatch(targets::tar_read_raw(
      paste0("gee_admin2_", country_labels[n]), store = STORE),
      error = function(e) NULL)
    if (is.null(gee)) next
    j <- dplyr::inner_join(agg, gee, by = "Admin2")
    lk <- centroid_lookup(n)
    if (is.null(lk)) next
    j$lon <- lk$lon[match(j$Admin2, lk$Admin2)]
    j$lat <- lk$lat[match(j$Admin2, lk$Admin2)]
    j$country <- n
    pool[[n]] <- j
  }
  if (length(pool) < 2) { cat("  skip — too few countries\n"); next }
  train <- dplyr::bind_rows(pool)
  train <- train[is.finite(train$lon) & is.finite(train$lat), , drop = FALSE]
  cat(sprintf("  training pool: %d polygons across %d countries\n",
              nrow(train), length(pool)))

  # Build GAM: spatial spline + linear top-8 soil.
  avail <- intersect(TOP8, colnames(train))
  if (length(avail) < 4) { cat("  skip — insufficient soil features\n"); next }
  imp_tr <- as.data.frame(train[, avail, drop = FALSE])
  med <- vapply(imp_tr, median, numeric(1), na.rm = TRUE)
  for (v in avail) imp_tr[[v]][!is.finite(imp_tr[[v]])] <- med[v]
  fit_df <- data.frame(Y = train$mean_bio, lon = train$lon, lat = train$lat,
                        imp_tr, check.names = FALSE)
  rhs <- paste(c("s(lon, lat, k = 20, bs = 'tp')", avail), collapse = " + ")
  fit <- tryCatch(mgcv::gam(stats::as.formula(paste("Y ~", rhs)),
                              data = fit_df, method = "REML",
                              weights = pmax(train$n_obs, 1)),
                   error = function(e) {cat("    gam error\n"); NULL})
  if (is.null(fit)) next

  # Predict at CIV admin-2.
  imp_te <- civ_named[, avail, drop = FALSE]
  for (v in avail) imp_te[[v]][!is.finite(imp_te[[v]])] <- med[v]
  test_df <- data.frame(lon = civ_named$lon, lat = civ_named$lat, imp_te,
                          check.names = FALSE)
  pred <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
  cat(sprintf("  CIV prediction summary: mean %.2f, sd %.2f, range [%.2f, %.2f]\n",
              mean(pred, na.rm = TRUE), stats::sd(pred, na.rm = TRUE),
              min(pred, na.rm = TRUE), max(pred, na.rm = TRUE)))
  all_civ_preds[[oc]] <- data.frame(
    outcome = oc, Admin1 = civ_named$Admin2 |>
      (\(x) civ_centroids_df$Admin1[match(x, civ_centroids_df$Admin2)])(),
    Admin2 = civ_named$Admin2,
    spatial_soil_pred = round(pred, 4),
    stringsAsFactors = FALSE)

  # Compare to existing oos_cote_divoire.rds (individual-level SL agg)
  oos <- readRDS(here::here("dashboard","data","oos_cote_divoire.rds"))
  if (is.list(oos) && !is.null(oos$predictions)) {
    oos_pred <- oos$predictions[oos$predictions$outcome == oc, , drop = FALSE]
    if (nrow(oos_pred) > 0) {
      cmp <- dplyr::inner_join(
        all_civ_preds[[oc]][, c("Admin2","spatial_soil_pred")],
        oos_pred[, c("Admin2","pred_prev")] |>
          dplyr::rename(individual_sl_agg = pred_prev),
        by = "Admin2")
      if (nrow(cmp) >= 3) {
        r <- suppressWarnings(cor(cmp$spatial_soil_pred,
                                    cmp$individual_sl_agg))
        cat(sprintf("  agreement with individual-SL-agg (oos_cote_divoire.rds):\n"))
        cat(sprintf("    Pearson r = %.3f (across %d admin-2 polygons)\n",
                    r, nrow(cmp)))
        cmp$outcome <- oc
        all_compare[[oc]] <- cmp
      }
    }
  }
}

# ── Save outputs ─────────────────────────────────────────────────────────
preds_df <- dplyr::bind_rows(all_civ_preds)
cmp_df   <- dplyr::bind_rows(all_compare)
readr::write_csv(preds_df,
                  file.path(OUT_DIR, "civ_external_validation.csv"))
readr::write_csv(cmp_df,
                  file.path(OUT_DIR, "civ_external_validation_compare.csv"))
cat(sprintf("\n===== summary =====\nspatial+soil predictions: %d rows\nagreement comparison: %d rows\n",
            nrow(preds_df), nrow(cmp_df)))

# ── Figure ───────────────────────────────────────────────────────────────
if (nrow(cmp_df) > 0) {
  fig <- ggplot(cmp_df, aes(x = individual_sl_agg, y = spatial_soil_pred)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, colour = "steelblue") +
    facet_wrap(~ outcome, scales = "free", ncol = 2) +
    labs(
      title = "Cote d'Ivoire external validation",
      subtitle = "Spatial + top-8 SoilGrids GAM trained on 4 countries vs. pipeline individual-SL aggregated",
      x = "Individual-level SL aggregated to admin-2 (existing pipeline)",
      y = "Spatial + top-8 SoilGrids GAM (this work)") +
    theme_minimal(base_size = 11)
  ggsave(file.path(OUT_DIR, "civ_external_validation.png"),
         fig, width = 9, height = 7, dpi = 150)
  cat(sprintf("figure: %s\n",
              file.path(OUT_DIR, "civ_external_validation.png")))
}

cat("\nDONE\n")
