# =============================================================================
# scripts/civ_external_validation_pipeline.R
#
# CIV external validation using PIPELINE-EXTRACTED admin-2 features.
#
# Differs from scripts/civ_external_validation.R only in step 2: instead of
# computing zonal means on its own (which produced features 4×–290× off the
# training scale), this script loads the features extracted by
# scripts/extract_civ_gee_admin2.R (which uses the same .append_gee_zonal_cols
# logic the targets pipeline uses for every other country). Steps 1, 3, 4 are
# unchanged.
#
# Output:
#   results/transportability/civ_external_validation_pipeline.csv
#   results/transportability/civ_external_validation_pipeline_compare.csv
#   results/transportability/civ_external_validation_pipeline.png
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here); library(sf)
  library(ggplot2); library(patchwork)
})
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))
source(here::here("sandbox", "00_setup.R"))

OUT_DIR <- here::here("results", "transportability")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. CIV admin-2 polygons (for centroids) ─────────────────────────────
cat("===== STEP 1: CIV admin-2 polygons + centroids =====\n")
civ_a2 <- geodata::gadm("CIV", level = 2, path = here::here("data", "gadm"))
civ_a2_sf <- sf::st_as_sf(civ_a2)
civ_centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(civ_a2_sf)))
civ_centroids_df <- data.frame(
  Admin1 = civ_a2_sf$NAME_1, Admin2 = civ_a2_sf$NAME_2,
  lon = civ_centroids[, 1], lat = civ_centroids[, 2],
  stringsAsFactors = FALSE)
cat(sprintf("CIV admin-2: %d polygons\n", nrow(civ_centroids_df)))

# ── 2. Load pipeline-extracted CIV features ─────────────────────────────
cat("\n===== STEP 2: load pipeline-extracted CIV gee_admin2 =====\n")
civ_feat_file <- file.path(OUT_DIR, "civ_gee_admin2.rds")
if (!file.exists(civ_feat_file)) {
  stop("Missing ", civ_feat_file,
       " — run scripts/extract_civ_gee_admin2.R first")
}
civ_gee <- readRDS(civ_feat_file)
cat(sprintf("CIV features: %d polygons x %d cols (incl. Admin1/Admin2)\n",
            nrow(civ_gee), ncol(civ_gee)))

# Join centroids (lon/lat) to the pipeline features on Admin2
civ_named <- dplyr::left_join(
  civ_gee, civ_centroids_df[, c("Admin2", "lon", "lat")], by = "Admin2"
)
cat(sprintf("Joined centroids: %d rows have lon/lat (of %d)\n",
            sum(is.finite(civ_named$lon)), nrow(civ_named)))

# Verify spatial+soil features present and report ranges (sanity-check scale)
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS
cat("\nSpatial+soil predictor ranges in CIV (pipeline extraction):\n")
for (v in TOP8) {
  if (!v %in% colnames(civ_named)) {
    cat(sprintf("  !! MISSING: %s\n", v)); next
  }
  vals <- civ_named[[v]]; vals <- vals[is.finite(vals)]
  if (!length(vals)) { cat(sprintf("  !! all-NA: %s\n", v)); next }
  cat(sprintf("  %-35s: [%.3f, %.3f]  (mean %.3f, sd %.3f)\n",
              v, min(vals), max(vals), mean(vals), stats::sd(vals)))
}

# ── 3. Train spatial+soil GAM on 4 countries, predict CIV ───────────────
cat("\n===== STEP 3: train spatial+soil, predict CIV =====\n")
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
  if (col == "gw_B12_pmol_L")    vals <- pmin(vals, 1500)
  if (col == "gw_Folate_nmol_L") vals <- pmin(vals, 50)
  ok <- is.finite(vals)
  if (sum(ok) < 3) return(NULL)
  data.frame(Admin2 = merged$Admin2[ok], value = vals[ok]) |>
    dplyr::group_by(Admin2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
}

centroid_lookup <- function(country_label) {
  ccm <- list(Gambia = "GMB", Ghana = "GHA",
               SierraLeone = "SLE", Malawi = "MWI")
  poly <- tryCatch(geodata::gadm(ccm[[country_label]], level = 2,
                                    path = here::here("data","gadm")),
                    error = function(e) NULL)
  if (is.null(poly)) return(NULL)
  sfp <- sf::st_as_sf(poly)
  ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sfp)))
  data.frame(Admin2 = sfp$NAME_2, lon = ctr[,1], lat = ctr[,2],
             stringsAsFactors = FALSE) |>
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

  avail <- intersect(TOP8, colnames(train))
  avail <- intersect(avail, colnames(civ_named))  # must be in BOTH
  if (length(avail) < 4) { cat("  skip — insufficient soil features\n"); next }
  cat(sprintf("  using %d/%d top-8 soil predictors available in both train+CIV\n",
              length(avail), length(TOP8)))

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

  imp_te <- civ_named[, avail, drop = FALSE]
  for (v in avail) imp_te[[v]][!is.finite(imp_te[[v]])] <- med[v]
  test_df <- data.frame(lon = civ_named$lon, lat = civ_named$lat, imp_te,
                          check.names = FALSE)
  pred <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
  cat(sprintf("  CIV prediction summary: mean %.2f, sd %.2f, range [%.2f, %.2f]\n",
              mean(pred, na.rm = TRUE), stats::sd(pred, na.rm = TRUE),
              min(pred, na.rm = TRUE), max(pred, na.rm = TRUE)))
  all_civ_preds[[oc]] <- data.frame(
    outcome = oc,
    Admin1 = civ_named$Admin1,
    Admin2 = civ_named$Admin2,
    spatial_soil_pred = round(pred, 4),
    stringsAsFactors = FALSE)

  # Compare to existing oos_cote_divoire.rds (individual-level SL agg)
  oos_path <- here::here("dashboard","data","oos_cote_divoire.rds")
  if (file.exists(oos_path)) {
    oos <- readRDS(oos_path)
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
}

# ── Save outputs ─────────────────────────────────────────────────────────
preds_df <- dplyr::bind_rows(all_civ_preds)
cmp_df   <- dplyr::bind_rows(all_compare)
readr::write_csv(preds_df,
                  file.path(OUT_DIR, "civ_external_validation_pipeline.csv"))
readr::write_csv(cmp_df,
                  file.path(OUT_DIR, "civ_external_validation_pipeline_compare.csv"))
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
      title = "Cote d'Ivoire external validation (pipeline-extracted features)",
      subtitle = "Spatial + top-8 SoilGrids GAM trained on 4 countries vs. pipeline individual-SL aggregated",
      x = "Individual-level SL aggregated to admin-2 (existing pipeline)",
      y = "Spatial + top-8 SoilGrids GAM (this work)") +
    theme_minimal(base_size = 11)
  ggsave(file.path(OUT_DIR, "civ_external_validation_pipeline.png"),
         fig, width = 9, height = 7, dpi = 150)
  cat(sprintf("figure: %s\n",
              file.path(OUT_DIR, "civ_external_validation_pipeline.png")))
}

cat("\nDONE\n")
