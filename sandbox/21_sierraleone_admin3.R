# =============================================================================
# sandbox/21_sierraleone_admin3.R
#
# Build a Sierra Leone admin-3 (chiefdom) parallel pipeline alongside the
# existing admin-2 pipeline. NEVER deletes the admin-2 results — admin-3
# is an ADDITIONAL comparison track.
#
# Steps:
#   1. Load DHS cluster coordinates for Sierra Leone (434 clusters with
#      lat/long).
#   2. Spatially join clusters to GADM admin-3 polygons (153 chiefdoms).
#   3. For each biomarker measurement in merged_sierraleone, look up its
#      cluster's admin-3 polygon.
#   4. Aggregate biomarker outcomes (binary svy_prev, continuous mean) per
#      admin-3 polygon.
#   5. Look up each admin-3 polygon's parent admin-2 and inherit the
#      admin-2 GEE features. (Avoids re-running heavy GEE raster extraction;
#      isolates the "does finer SPATIAL resolution help?" question from
#      "does finer GEE help?".)
#   6. Build a pooled LOCO frame where SL uses admin-3, other countries
#      use admin-2. Run LOCO and compare to admin-2-only baseline.
#
# Outputs:
#   results/transportability/sierraleone_admin3_svy.rds
#   results/transportability/sierraleone_admin3_gee.rds
#   sandbox/results/21_sl_admin3_vs_admin2.csv
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here); library(sf)
})
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))
source(here::here("sandbox", "00_setup.R"))

OUT_DIR <- here::here("results", "transportability")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load SL DHS cluster coordinates ───────────────────────────────────
cat("===== STEP 1: load SL cluster coordinates =====\n")
sl_cluster_info <- readRDS(here::here("data", "DHS", "clean",
                                        "Sierra Leone_2013_cluster_admin_info.rds"))
ci <- sl_cluster_info$cluster.info$data
cat(sprintf("DHS clusters in SL: %d\n", nrow(ci)))
cat(sprintf("  with lat/long: %d\n",
            sum(!is.na(ci$LONGNUM) & !is.na(ci$LATNUM))))

# ── 2. Load GADM admin-3 polygons + spatial-join clusters ────────────────
cat("\n===== STEP 2: spatially join clusters to admin-3 =====\n")
sl_a3 <- geodata::gadm("SLE", level = 3, path = here::here("data", "gadm"))
sl_a3_sf <- sf::st_as_sf(sl_a3)
cat(sprintf("GADM admin-3: %d chiefdoms\n", nrow(sl_a3_sf)))

# Also load admin-2 to map admin-3 → admin-2
sl_a2 <- geodata::gadm("SLE", level = 2, path = here::here("data", "gadm"))
sl_a2_sf <- sf::st_as_sf(sl_a2)
cat(sprintf("GADM admin-2: %d districts\n", nrow(sl_a2_sf)))

# Build cluster sf points
ci_sf <- sf::st_as_sf(ci, coords = c("LONGNUM", "LATNUM"),
                       crs = sf::st_crs(sl_a3_sf), na.fail = FALSE)
# Spatial join clusters → admin-3 polygons
cluster_a3 <- suppressMessages(suppressWarnings(
  sf::st_join(ci_sf, sl_a3_sf[, c("NAME_2", "NAME_3")],
              join = sf::st_within)))
cluster_a3_df <- sf::st_drop_geometry(cluster_a3)
cat(sprintf("clusters assigned to admin-3: %d / %d\n",
            sum(!is.na(cluster_a3_df$NAME_3)), nrow(cluster_a3_df)))
cat("  sample of admin-3 names found:\n")
print(head(sort(unique(cluster_a3_df$NAME_3)), 10))
n_a3_with_clusters <- length(unique(cluster_a3_df$NAME_3[!is.na(cluster_a3_df$NAME_3)]))
cat(sprintf("\n%d admin-3 polygons have at least one cluster\n",
            n_a3_with_clusters))

# ── 3. Aggregate biomarker measurements by admin-3 ───────────────────────
cat("\n===== STEP 3: aggregate biomarker measurements per admin-3 =====\n")
merged_sle <- targets::tar_read_raw("merged_sierraleone", store = STORE)
cat(sprintf("merged_sle: %d rows, unique clusters %d\n",
            nrow(merged_sle), length(unique(merged_sle$gw_wClusterNumb))))

# Attach admin-3 label to each row via cluster number.
clust_to_a3 <- cluster_a3_df[, c("cluster", "NAME_3", "NAME_2")]
merged_sle$cluster <- as.numeric(unclass(merged_sle$gw_wClusterNumb))
merged_sle_a3 <- merged_sle |>
  dplyr::left_join(clust_to_a3, by = "cluster")
n_no_match <- sum(is.na(merged_sle_a3$NAME_3))
cat(sprintf("rows without admin-3 assignment: %d / %d\n",
            n_no_match, nrow(merged_sle_a3)))

# ── 4. Build admin-3 svy_prev frames per outcome ─────────────────────────
# Use the harmonized binary deficiency columns where available.
# Sierra Leone outcome columns (from earlier audit):
outcome_cols_sle <- list(
  child_vitA   = "gw_cRBPAdj",      # not directly binary; we use the binary derived from threshold
  women_vitA   = "gw_wRBPAdj",
  child_iron   = "gw_wFerritin",
  women_iron   = "gw_wFerritin",
  women_b12    = "gw_B12",
  women_folate = "gw_wFolate"
)
# But the BINARY harmonised outcomes live in svy_admin2 directly. Easier:
# use the SAME definition as the existing svy_admin2_sierraleone_<outcome>
# targets by reading them and inferring the threshold per outcome.

# Alternative simple approach: compute admin-3 MEAN of the continuous
# biomarker (matches sandbox/19 framing — works for vitA / B12 / folate
# transportability and is unit-agnostic via Spearman).

continuous_biomarker_cols <- list(
  child_vitA   = "gw_cRBP",
  women_vitA   = "gw_cRBP",
  women_b12    = "gw_B12",
  women_folate = "gw_wFolate"
  # iron: ferritin vs sTFR vs hemoglobin across countries is the
  # apples-to-oranges problem — skip.
)

admin3_mean <- function(df, col, id_col = "NAME_3") {
  if (!col %in% colnames(df)) return(NULL)
  vals <- suppressWarnings(as.numeric(unclass(df[[col]])))
  ok <- is.finite(vals) & !is.na(df[[id_col]])
  if (sum(ok) < 3) return(NULL)
  d <- data.frame(NAME_3 = df[[id_col]][ok],
                   NAME_2 = df[["NAME_2"]][ok],
                   value  = vals[ok])
  d |> dplyr::group_by(NAME_3, NAME_2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
}

sl_admin3_biomarkers <- list()
for (oc in names(continuous_biomarker_cols)) {
  col <- continuous_biomarker_cols[[oc]]
  agg <- admin3_mean(merged_sle_a3, col)
  if (!is.null(agg)) {
    cat(sprintf("  %s (%s): %d admin-3 polygons, mean %.2f, sd %.2f\n",
                oc, col, nrow(agg),
                mean(agg$mean_bio, na.rm = TRUE),
                stats::sd(agg$mean_bio, na.rm = TRUE)))
    sl_admin3_biomarkers[[oc]] <- agg
  }
}

# ── 5. Inherit admin-2 GEE features at admin-3 + add centroids ────────────
cat("\n===== STEP 5: build admin-3 features (admin-2 GEE replicated + admin-3 centroids) =====\n")
gee_a2_sle <- targets::tar_read_raw("gee_admin2_sierraleone", store = STORE)
# Map each admin-3 polygon to its admin-2 NAME_2 (already present from spatial join).
# Compute admin-3 centroid coordinates from GADM.
sl_a3_centroids <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sl_a3_sf)))
sl_a3_centroids_df <- data.frame(
  NAME_2 = sl_a3_sf$NAME_2,
  NAME_3 = sl_a3_sf$NAME_3,
  lon = sl_a3_centroids[, 1],
  lat = sl_a3_centroids[, 2])
cat(sprintf("admin-3 polygons with centroids: %d\n", nrow(sl_a3_centroids_df)))

# Join admin-2 GEE features by NAME_2 (Admin2 in gee_a2_sle).
# gee_a2_sle has "Admin2" as the polygon name. Map "NAME_2" → "Admin2".
sl_a3_features <- sl_a3_centroids_df |>
  dplyr::left_join(gee_a2_sle, by = c("NAME_2" = "Admin2"))
cat(sprintf("admin-3 features assembled: %d rows × %d cols\n",
            nrow(sl_a3_features), ncol(sl_a3_features)))

# ── 6. Build pooled LOCO frame: SL at admin-3, other countries at admin-2 ─
cat("\n===== STEP 6: pooled LOCO frame (SL admin-3 + others admin-2) =====\n")
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS

# Load other-country biomarker data (admin-2 means) for each outcome.
other_country_admin2_bio <- function(outcome_tag) {
  bio_cols <- list(
    Gambia = if (outcome_tag %in% c("child_vitA","women_vitA")) "gw_cRBP" else NULL,
    Ghana  = list(child_vitA = "gw_cRBP", women_vitA = "gw_cRBP",
                   women_b12 = "gw_B12_pmol_L",
                   women_folate = "gw_Folate_nmol_L")[[outcome_tag]],
    Malawi = list(child_vitA = "rbp", women_vitA = "rbp",
                   women_b12 = "vitb12",
                   women_folate = "rbcf")[[outcome_tag]]
  )
  out <- list()
  cc_map <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi")
  for (n in names(bio_cols)) {
    col <- bio_cols[[n]]
    if (is.null(col) || length(col) == 0) next
    merged <- tryCatch(targets::tar_read_raw(paste0("merged_", cc_map[n]),
                                               store = STORE),
                        error = function(e) NULL)
    if (is.null(merged) || !col %in% colnames(merged)) next
    vals <- suppressWarnings(as.numeric(unclass(merged[[col]])))
    # Apply data-quality caps from sandbox/20.
    if (n == "Ghana" && col == "gw_B12_pmol_L") vals <- pmin(vals, 1500)
    if (n == "Ghana" && col == "gw_Folate_nmol_L") vals <- pmin(vals, 50)
    ok <- is.finite(vals)
    df <- data.frame(Admin2 = merged$Admin2[ok], value = vals[ok])
    agg <- df |> dplyr::group_by(Admin2) |>
      dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                        n_obs    = dplyr::n(), .groups = "drop")
    gee <- tryCatch(targets::tar_read_raw(paste0("gee_admin2_", cc_map[n]),
                                            store = STORE),
                     error = function(e) NULL)
    if (is.null(gee)) next
    j <- dplyr::inner_join(agg, gee, by = "Admin2")
    j$country <- n
    j$polygon_id <- j$Admin2
    out[[n]] <- j
  }
  out
}

# Run LOCO at admin-2 (sandbox/19 baseline) and admin-3 (this script's
# question). For admin-3, swap SL's admin-2 frame for the admin-3 frame.

run_loco_mixed <- function(outcome_tag, sl_resolution = c("admin2","admin3")) {
  sl_resolution <- match.arg(sl_resolution)
  cc_list <- get_country_configs()
  others <- other_country_admin2_bio(outcome_tag)
  # SL: admin-2 baseline OR admin-3 expansion.
  if (sl_resolution == "admin2") {
    bio_col <- continuous_biomarker_cols[[outcome_tag]]
    merged_full <- merged_sle  # for SL specifically
    vals <- suppressWarnings(as.numeric(unclass(merged_full[[bio_col]])))
    ok <- is.finite(vals)
    sl_agg <- data.frame(Admin2 = merged_full$Admin2[ok], value = vals[ok]) |>
      dplyr::group_by(Admin2) |>
      dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                        n_obs    = dplyr::n(), .groups = "drop")
    sl_j <- dplyr::inner_join(sl_agg, gee_a2_sle, by = "Admin2")
    sl_j$country <- "SierraLeone"
    sl_j$polygon_id <- sl_j$Admin2
  } else {
    bio <- sl_admin3_biomarkers[[outcome_tag]]
    if (is.null(bio)) return(NULL)
    sl_j <- bio |>
      dplyr::left_join(sl_a3_features, by = c("NAME_3", "NAME_2"))
    sl_j$country <- "SierraLeone"
    sl_j$polygon_id <- sl_j$NAME_3
  }
  # Add lon/lat to ALL frames before binding.
  # OTHER countries -> admin-2 GADM centroids via add_admin2_centroids().
  # SL admin-2 case -> add_admin2_centroids on its own frame.
  # SL admin-3 case -> lon/lat already present from sl_a3_centroids_df.
  svy_list_other <- setNames(
    lapply(c("gambia","ghana","malawi"), function(c)
      tryCatch(targets::tar_read_raw(paste0("svy_admin2_", c, "_child_vitA"),
                                       store = STORE),
                error = function(e) NULL)),
    c("Gambia","Ghana","Malawi"))
  # Fresh per-country centroid lookup (avoids stale cache duplicates).
  centroid_lookup_for <- function(country_label) {
    ccm <- list(Gambia = "GMB", Ghana = "GHA",
                 SierraLeone = "SLE", Malawi = "MWI")
    code <- ccm[[country_label]]
    if (is.null(code)) return(NULL)
    poly <- tryCatch(
      geodata::gadm(code, level = 2, path = here::here("data","gadm")),
      error = function(e) NULL)
    if (is.null(poly)) return(NULL)
    sfp <- sf::st_as_sf(poly)
    ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sfp)))
    data.frame(Admin2 = sfp$NAME_2, lon = ctr[, 1], lat = ctr[, 2],
                stringsAsFactors = FALSE) |>
      dplyr::distinct(Admin2, .keep_all = TRUE)
  }
  ensure_coords <- function(df, country_label) {
    if (nrow(df) == 0) return(df)
    if (all(c("lon","lat") %in% colnames(df)) &&
        !all(is.na(df$lon))) return(df)
    lookup <- centroid_lookup_for(country_label)
    if (is.null(lookup)) {
      df$lon <- NA_real_; df$lat <- NA_real_; return(df)
    }
    key_col <- if ("polygon_id" %in% colnames(df)) "polygon_id" else "Admin2"
    idx <- match(df[[key_col]], lookup$Admin2)
    df$lon <- lookup$lon[idx]
    df$lat <- lookup$lat[idx]
    df
  }
  others_with_coords <- Map(function(df, n) ensure_coords(df, n),
                              others, names(others))
  if (sl_resolution == "admin2") {
    sl_j <- ensure_coords(sl_j, "SierraLeone")
  }
  shared_meta <- c("polygon_id","country","mean_bio","n_obs","lon","lat")
  all_gee <- Reduce(intersect, c(
    lapply(others_with_coords, function(d) grep("^gee_", colnames(d), value = TRUE)),
    list(grep("^gee_", colnames(sl_j), value = TRUE))
  ))
  thin_frame <- function(d) {
    if (!"polygon_id" %in% colnames(d)) {
      d$polygon_id <- d$Admin2 %||% NA_character_
    }
    for (m in shared_meta) if (!m %in% colnames(d)) d[[m]] <- NA
    d[, c(shared_meta, intersect(all_gee, colnames(d))), drop = FALSE]
  }
  others_thin <- lapply(others_with_coords, thin_frame)
  sl_thin     <- thin_frame(sl_j)
  all_data <- dplyr::bind_rows(c(others_thin, list(SierraLeone = sl_thin)))
  ok <- !is.na(all_data$lon) & !is.na(all_data$lat) &
        is.finite(all_data$mean_bio)
  all_data <- all_data[ok, , drop = FALSE]

  # Spatial+top-8-soil LOCO over country.
  countries <- unique(all_data$country)
  results <- list()
  for (held in countries) {
    train <- all_data[all_data$country != held, , drop = FALSE]
    test  <- all_data[all_data$country == held, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    avail <- intersect(TOP8, colnames(train))
    Y <- train$mean_bio
    wt <- pmax(train$n_obs, 1)
    imp_tr <- as.data.frame(train[, avail, drop = FALSE])
    imp_te <- as.data.frame(test [, avail, drop = FALSE])
    med <- vapply(imp_tr, median, numeric(1), na.rm = TRUE)
    for (v in avail) {
      imp_tr[[v]][!is.finite(imp_tr[[v]])] <- med[v]
      imp_te[[v]][!is.finite(imp_te[[v]])] <- med[v]
    }
    fit_df <- data.frame(Y = Y, lon = train$lon, lat = train$lat,
                          imp_tr, check.names = FALSE)
    test_df <- data.frame(   lon = test$lon, lat = test$lat,
                              imp_te, check.names = FALSE)
    rhs <- paste(c("s(lon, lat, k = 20, bs = 'tp')", avail), collapse = " + ")
    fit <- tryCatch(mgcv::gam(stats::as.formula(paste("Y ~", rhs)),
                                data = fit_df, method = "REML", weights = wt),
                     error = function(e) NULL)
    if (is.null(fit)) next
    p <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
    ok2 <- is.finite(p) & is.finite(test$mean_bio)
    if (sum(ok2) < 3) next
    rs <- suppressWarnings(cor(test$mean_bio[ok2], p[ok2], method = "spearman"))
    rp <- suppressWarnings(cor(test$mean_bio[ok2], p[ok2]))
    results[[held]] <- data.frame(
      outcome = outcome_tag, sl_resolution = sl_resolution,
      held_out = held, n_test = sum(ok2),
      pearson_r = round(rp, 3), spearman_r = round(rs, 3))
  }
  dplyr::bind_rows(results)
}

cat("\n===== STEP 6: LOCO admin-2 vs admin-3 =====\n")
all_results <- list()
for (oc in names(continuous_biomarker_cols)) {
  cat(sprintf("\n--- %s ---\n", oc))
  r2 <- run_loco_mixed(oc, sl_resolution = "admin2")
  r3 <- run_loco_mixed(oc, sl_resolution = "admin3")
  if (!is.null(r2) && nrow(r2) > 0) cat("admin2: "); print(r2, row.names=FALSE)
  if (!is.null(r3) && nrow(r3) > 0) cat("admin3: "); print(r3, row.names=FALSE)
  all_results[[oc]] <- dplyr::bind_rows(r2, r3)
}
final <- dplyr::bind_rows(all_results)
readr::write_csv(final, file.path(SANDBOX_RESULTS,
                                    "21_sl_admin3_vs_admin2.csv"))

# Save the admin-3 data products for future reuse (NEVER overwriting admin-2).
saveRDS(sl_admin3_biomarkers, file.path(OUT_DIR, "sierraleone_admin3_svy.rds"))
saveRDS(sl_a3_features,       file.path(OUT_DIR, "sierraleone_admin3_gee.rds"))
saveRDS(cluster_a3_df,        file.path(OUT_DIR, "sierraleone_cluster_to_admin3.rds"))

cat("\n===== summary: SL admin-3 vs admin-2 LOCO =====\n")
sl_only <- final[final$held_out == "SierraLeone", ]
print(sl_only, row.names = FALSE)

cat(sprintf("\nwritten:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(SANDBOX_RESULTS, "21_sl_admin3_vs_admin2.csv"),
            file.path(OUT_DIR, "sierraleone_admin3_svy.rds"),
            file.path(OUT_DIR, "sierraleone_admin3_gee.rds"),
            file.path(OUT_DIR, "sierraleone_cluster_to_admin3.rds")))
