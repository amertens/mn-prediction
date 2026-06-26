# =============================================================================
# sandbox/04_spatial_gam_coords.R
#
# Hypothesis H3: spatial autocorrelation is enough. Fit a thin-plate spline
# on (longitude, latitude) of admin-2 polygon centroids, conditioning on a
# few key covariates (or none). If micronutrient deficiency is mainly a
# regional/agroecological gradient, this should capture most of the
# transportable signal *without any feature engineering on the 152 GEE
# rasters*.
#
# Rationale: many widely-cited disease-mapping methods (BYM2, kriging, INLA-
# SPDE) are essentially spatial smoothers. Within-country admin-2 ranking
# may be driven by smooth geographic gradients (cereal-belt latitude,
# coastal vs interior, urban gravity) — let the spline learn that directly
# instead of inferring it from a wide proxy panel.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

# Pull admin-2 centroids from GADM (cached).
ensure_centroids <- function() {
  cache <- file.path(SANDBOX_CACHE, "admin2_centroids.rds")
  if (file.exists(cache)) return(readRDS(cache))
  source(here::here("R", "config.R"))
  cc_list <- get_country_configs()
  out <- list()
  for (cc in cc_list) {
    poly <- tryCatch(
      geodata::gadm(cc$gadm_code, level = 2,
                     path = here::here("data", "gadm")),
      error = function(e) NULL)
    if (is.null(poly)) next
    sf_poly <- sf::st_as_sf(poly)
    suppressWarnings({
      ctr <- sf::st_coordinates(sf::st_centroid(sf_poly))
    })
    out[[cc$country]] <- data.frame(
      country = cc$country, Admin2 = sf_poly$NAME_2,
      lon = ctr[, 1], lat = ctr[, 2], stringsAsFactors = FALSE)
  }
  centroids <- dplyr::bind_rows(out)
  # Normalize country naming to match svy_admin2_list display names
  # ("Gambia", "Ghana", "SierraLeone", "Malawi"). Robust to either
  # capitalisation in cc$country.
  norm_lookup <- function(x) {
    canonical <- c(gambia = "Gambia", ghana = "Ghana",
                   sierraleone = "SierraLeone", "sierra leone" = "SierraLeone",
                   malawi = "Malawi")
    out <- canonical[gsub(" ", "", tolower(x))]
    out[is.na(out)] <- x[is.na(out)]   # fallback to original
    out
  }
  centroids$country <- norm_lookup(centroids$country)
  saveRDS(centroids, cache)
  centroids
}

# Augment pooled_data with lon/lat from centroids.
add_coords <- function(pd, centroids) {
  dplyr::left_join(pd, centroids, by = c("country", "Admin2"))
}

fit_predict_gam_spatial <- function(train, test, vars = NULL,
                                       k_spline = 30, extra_covars = character()) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  # train/test must already have lon, lat (added via add_coords upstream).
  if (any(is.na(train$lon)) || any(is.na(train$lat))) return(NULL)
  Y  <- train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(train)), 1)
  if (length(extra_covars) > 0) {
    imp <- impute_median(train, test, extra_covars)
    fit_df <- data.frame(Y = Y, lon = train$lon, lat = train$lat, imp$train)
    test_df <- data.frame(lon = test$lon, lat = test$lat, imp$test)
    rhs <- paste(c("te(lon, lat, k = c(k_spline, k_spline))",
                    extra_covars), collapse = " + ")
    rhs <- gsub("k_spline", as.character(k_spline), rhs)
  } else {
    fit_df  <- data.frame(Y = Y, lon = train$lon, lat = train$lat)
    test_df <- data.frame(    lon = test$lon, lat = test$lat)
    rhs <- sprintf("s(lon, lat, k = %d, bs = 'tp')", k_spline)
  }
  form <- stats::as.formula(paste("Y ~", rhs))
  fit <- tryCatch(mgcv::gam(form, data = fit_df, weights = wt, method = "REML"),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  p <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
  pmin(pmax(p, 0), 1)
}

centroids <- ensure_centroids()
cat("centroids loaded:", nrow(centroids), "polygons\n")

pooled_all <- load_all_pooled()
# Augment with coords.
pooled_all <- lapply(pooled_all, function(p) {
  p$pooled_data <- add_coords(p$pooled_data, centroids); p
})

res <- list()
configs <- list(
  list(k_spline = 15, extra_covars = character(),                              name = "gam_coords_k15"),
  list(k_spline = 30, extra_covars = character(),                              name = "gam_coords_k30"),
  list(k_spline = 50, extra_covars = character(),                              name = "gam_coords_k50"),
  list(k_spline = 30,
        extra_covars = c("gee_accessibility_2019","gee_elevation_2000"),         name = "gam_coords_k30_2cov")
)
for (cfg in configs) {
  cat(sprintf("\n[H3] %s ...\n", cfg$name))
  r <- loco_evaluate_all(pooled_all,
    predict_fn  = function(train, test, vars)
      fit_predict_gam_spatial(train, test, vars,
                                k_spline = cfg$k_spline,
                                extra_covars = cfg$extra_covars),
    method_name = cfg$name)
  res[[cfg$name]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "04_spatial_gam.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "04_spatial_gam.csv")))
