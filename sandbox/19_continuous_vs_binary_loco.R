# =============================================================================
# sandbox/19_continuous_vs_binary_loco.R
#
# Tests whether predicting the CONTINUOUS biomarker concentration (admin-2
# mean) transports better than predicting the binary prevalence above WHO
# threshold, for the outcomes that failed binary LOCO:
#
#   women_vitA  — RBP concentration (all 4 countries)
#   women_folate — folate concentration (3 countries: Ghana, SL, Malawi)
#   women_b12   — B12 concentration (3 countries: Ghana, SL, Malawi)
#
# Evaluation: Spearman rank correlation between predicted and observed
# admin-2 mean. Rank correlation is unit-agnostic, so countries that
# measure on different scales (Ghana B12 in pmol/L, Malawi B12 in
# unclear units) still produce comparable rankings.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))
source(here::here("R", "area_level_comparison.R"))
source(here::here("R", "config.R"))

countries_lower <- c("gambia","ghana","sierraleone","malawi")
country_labels  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")

# Per-country biomarker columns identified from sandbox/18 audit.
# NULL = country lacks this biomarker.
BIOMARKER_COLS <- list(
  women_vitA   = list(Gambia = "gw_cRBP",         Ghana = "gw_cRBP",
                       SierraLeone = "gw_cRBP",    Malawi = "rbp"),
  child_vitA   = list(Gambia = "gw_cRBP",         Ghana = "gw_cRBP",
                       SierraLeone = "gw_cRBP",    Malawi = "rbp"),
  women_b12    = list(Gambia = NULL,              Ghana = "gw_B12_pmol_L",
                       SierraLeone = "gw_B12",     Malawi = "vitb12"),
  women_folate = list(Gambia = NULL,              Ghana = "gw_Folate_nmol_L",
                       SierraLeone = "gw_wFolate", Malawi = "rbcf")
)

# Filter merged_<country> to women / children at risk per outcome. The DHS
# data may have multiple respondents per household; we want per-admin-2
# means weighted by sample weight if available, else unweighted.

compute_admin2_mean_biomarker <- function(merged, biomarker_col, sex_filter = NULL) {
  if (is.null(merged) || is.null(biomarker_col) || !biomarker_col %in% colnames(merged))
    return(NULL)
  if (!"Admin2" %in% colnames(merged)) return(NULL)
  vals <- suppressWarnings(as.numeric(unclass(merged[[biomarker_col]])))
  ok <- is.finite(vals)
  df <- data.frame(Admin2 = merged$Admin2[ok], value = vals[ok],
                    stringsAsFactors = FALSE)
  agg <- df |> dplyr::group_by(Admin2) |>
    dplyr::summarise(mean_bio = mean(value, na.rm = TRUE),
                      n_obs    = dplyr::n(), .groups = "drop")
  agg
}

build_pooled_continuous <- function(outcome_tag) {
  bm <- BIOMARKER_COLS[[outcome_tag]]
  per_country <- list()
  for (n in names(bm)) {
    ctry <- country_labels[[n]]
    merged <- tryCatch(targets::tar_read_raw(paste0("merged_", ctry),
                                               store = STORE),
                        error = function(e) NULL)
    gee    <- tryCatch(targets::tar_read_raw(paste0("gee_admin2_", ctry),
                                               store = STORE),
                        error = function(e) NULL)
    if (is.null(merged) || is.null(gee)) next
    bio <- compute_admin2_mean_biomarker(merged, bm[[n]])
    if (is.null(bio) || nrow(bio) < 3) next
    df <- dplyr::inner_join(bio, gee, by = "Admin2")
    df$country <- n
    cat(sprintf("  %s (%s): %d admin-2 polygons, biomarker mean = %.2f (sd %.2f)\n",
                n, bm[[n]], nrow(df), mean(df$mean_bio, na.rm = TRUE),
                stats::sd(df$mean_bio, na.rm = TRUE)))
    per_country[[n]] <- df
  }
  if (length(per_country) < 2) return(NULL)
  # Pool with common GEE vars only.
  gee_vars <- Reduce(intersect, lapply(per_country, function(x)
    grep("^gee_", colnames(x), value = TRUE)))
  pooled <- dplyr::bind_rows(lapply(per_country, function(x)
    x[, c("country", "Admin2", "mean_bio", "n_obs", gee_vars)]))
  list(pooled_data = pooled, gee_vars = gee_vars,
        country_names = names(per_country))
}

# Per-country centroids.
cc_list <- get_country_configs()

# Use the winning spatial+top-8-soil recipe to predict the continuous outcome.
TOP8 <- .SPATIAL_PLUS_SOIL_DEFAULT_VARS

predict_continuous_spatial_plus_soil <- function(train, test, vars) {
  # Use mgcv::gam with thin-plate spline on (lon, lat) + linear soil.
  if (!all(c("lon","lat") %in% colnames(train))) return(NULL)
  available_soil <- intersect(TOP8, colnames(train))
  # Impute missing soil vars with column means on training set.
  imp_tr <- as.data.frame(train[, available_soil, drop = FALSE])
  imp_te <- as.data.frame(test [, available_soil, drop = FALSE])
  med <- vapply(imp_tr, median, numeric(1), na.rm = TRUE)
  for (v in available_soil) {
    imp_tr[[v]][!is.finite(imp_tr[[v]])] <- med[v]
    imp_te[[v]][!is.finite(imp_te[[v]])] <- med[v]
  }
  wt <- pmax(train$n_obs, 1)
  Y  <- train$mean_bio
  fit_df  <- data.frame(Y = Y, lon = train$lon, lat = train$lat, imp_tr,
                          check.names = FALSE)
  test_df <- data.frame(     lon = test$lon,  lat = test$lat,  imp_te,
                          check.names = FALSE)
  rhs <- paste(c("s(lon, lat, k = 20, bs = 'tp')", available_soil),
                collapse = " + ")
  form <- stats::as.formula(paste("Y ~", rhs))
  fit <- tryCatch(mgcv::gam(form, data = fit_df, weights = wt,
                              method = "REML"),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  as.numeric(mgcv::predict.gam(fit, newdata = test_df))
}

# ── Main: run LOCO on continuous biomarker for each outcome ─────────────
results <- list()
for (oc in names(BIOMARKER_COLS)) {
  cat(sprintf("\n===== outcome: %s (continuous biomarker mean) =====\n", oc))
  p <- build_pooled_continuous(oc)
  if (is.null(p) || length(p$country_names) < 2) {
    cat("  skip — too few countries with biomarker data\n"); next
  }
  # Filter to common svy_admin2_list for centroid join.
  svy_list <- setNames(lapply(country_labels[p$country_names], function(c)
    tryCatch(targets::tar_read_raw(paste0("svy_admin2_", c, "_", oc),
                                     store = STORE),
              error = function(e) NULL)),
    p$country_names)
  # Drop nulls
  ok <- !vapply(svy_list, is.null, logical(1))
  svy_list <- svy_list[ok]
  p$pooled_data <- add_admin2_centroids(p$pooled_data, cc_list, svy_list)
  # LOCO over countries.
  for (held_out in p$country_names) {
    train <- p$pooled_data[p$pooled_data$country != held_out, , drop = FALSE]
    test  <- p$pooled_data[p$pooled_data$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    preds <- predict_continuous_spatial_plus_soil(train, test, p$gee_vars)
    if (is.null(preds) || length(preds) != nrow(test)) next
    obs <- test$mean_bio
    ok <- is.finite(obs) & is.finite(preds)
    obs <- obs[ok]; preds <- preds[ok]
    if (length(obs) < 3) next
    r_p <- suppressWarnings(cor(obs, preds))
    r_s <- suppressWarnings(cor(obs, preds, method = "spearman"))
    cat(sprintf("  held_out %s (n_test=%d): Pearson %+.3f, Spearman %+.3f\n",
                held_out, length(obs), r_p, r_s))
    results[[paste(oc, held_out)]] <- data.frame(
      outcome = oc, held_out = held_out, target = "continuous_mean",
      n_test = length(obs),
      pearson_r = round(r_p, 3),
      spearman_r = round(r_s, 3),
      stringsAsFactors = FALSE)
  }
}

out <- dplyr::bind_rows(results)
cat("\n===== continuous biomarker LOCO summary =====\n")
print(as.data.frame(out), row.names = FALSE)

# Compare to binary LOCO from the existing transportability_area_loco_metrics.
binary_path <- here::here("results", "tables",
                            "transportability_area_loco_metrics.csv")
if (file.exists(binary_path)) {
  binary <- read.csv(binary_path)
  cat("\n===== binary prevalence LOCO (from area-transport pipeline) =====\n")
  binary$held_out_norm <- gsub(" ", "", binary$held_out)
  binary_sub <- binary[binary$outcome %in% names(BIOMARKER_COLS) &
                          binary$held_out_norm %in% out$held_out, ]
  print(binary_sub[, c("outcome", "held_out", "pearson_r", "spearman_r")],
        row.names = FALSE)

  # Side-by-side comparison
  cmp <- merge(
    out[, c("outcome","held_out","pearson_r","spearman_r")],
    binary_sub[, c("outcome","held_out_norm","pearson_r","spearman_r")] |>
      setNames(c("outcome","held_out","binary_pearson","binary_spearman")),
    by = c("outcome","held_out"), all = TRUE)
  cat("\n===== continuous vs binary LOCO Spearman r =====\n")
  print(as.data.frame(cmp[, c("outcome","held_out",
                                "spearman_r", "binary_spearman")]),
        row.names = FALSE)
}

readr::write_csv(out, file.path(SANDBOX_RESULTS,
                                  "19_continuous_biomarker_loco.csv"))
cat(sprintf("\nwritten: %s\n",
            file.path(SANDBOX_RESULTS, "19_continuous_biomarker_loco.csv")))
