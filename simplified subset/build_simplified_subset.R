# =============================================================================
# build_simplified_subset.R
# -----------------------------------------------------------------------------
# Builds the self-contained "simplified subset": ALL micronutrient-deficiency
# outcomes for the four study countries, aggregated to three geographic levels
# (cluster / Admin-2 / Admin-1), in WIDE format — one row per area, with a
# {prevalence, n, n_deficient, design_variance} block PER outcome plus a shared
# set of 16 proxy predictors. The wide multi-outcome layout lets a collaborator
# run per-outcome models AND multi-outcome / sequential models (e.g. predict
# vitamin A using predicted iron as a covariate), exploiting the covariance
# among deficiencies.
#
# Reads the project's _targets_full outcome-data objects (real pipeline output;
# no values fabricated). Does NOT modify the production pipeline or dashboard.
# Run:
#   "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla \
#       "simplified subset/build_simplified_subset.R"
# All outputs are written atomically (temp file + rename) for OneDrive safety.
# =============================================================================

suppressWarnings(suppressMessages({ library(dplyr); library(tidyr) }))

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
if (!dir.exists(file.path(ROOT, "_targets_full"))) ROOT <- getwd()
STORE   <- file.path(ROOT, "_targets_full", "objects")
SUB     <- file.path(ROOT, "simplified subset")
OUT_DIR <- file.path(SUB, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "R", "config.R"))

clip01 <- function(x) pmin(pmax(x, 0), 1)
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
write_csv_atomic <- function(df, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(df, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) file.remove(path)
  file.rename(tmp, path); invisible(path)
}

# ── Countries (config key / object suffix / display label) ───────────────────
COUNTRIES <- tibble::tribble(
  ~key,          ~low,          ~label,
  "Gambia",      "gambia",      "Gambia",
  "Ghana",       "ghana",       "Ghana",
  "SierraLeone", "sierraleone", "Sierra Leone",
  "Malawi",      "malawi",      "Malawi"
)

# ── Outcome display order + human-readable labels ────────────────────────────
# All binary micronutrient-deficiency outcomes configured in the pipeline.
# Availability varies by country (NA columns where an outcome was not measured).
OUTCOME_ORDER  <- c("women_iron", "women_vitA", "women_b12", "women_folate",
                    "women_zinc", "child_iron", "child_vitA", "child_zinc")
OUTCOME_LABELS <- c(
  women_iron   = "Iron deficiency, women 15-49 y (ferritin-based)",
  women_vitA   = "Vitamin A deficiency, women (RBP-based)",
  women_b12    = "Vitamin B12 deficiency, women (serum B12)",
  women_folate = "Folate deficiency, women (serum folate)",
  women_zinc   = "Zinc deficiency, women (serum zinc)",
  child_iron   = "Iron deficiency, children 6-59 m (ferritin-based)",
  child_vitA   = "Vitamin A deficiency, children (RBP-based)",
  child_zinc   = "Zinc deficiency, children (serum zinc)"
)
DEFF_AREA <- 1.5   # DHS-typical design effect, applied at Admin-2 / Admin-1

# ── The 16 shared proxy predictors (importance-ranked, non-redundant) ────────
# Cross-country-available geospatial proxies (cluster-GPS buffer extractions),
# selected from the pipeline's permutation importance for iron in women and kept
# for all outcomes as a common proxy set. See data_dictionary.csv / README.md.
pred_meta <- tibble::tribble(
  ~raw,                                          ~clean,                       ~description, ~units, ~source,
  "gee_trmm_50km_annual_mean",  "rain_annual_mean_50km",   "Mean annual precipitation, 50 km buffer around the cluster", "mm/month (TRMM units)", "TRMM precipitation (GEE cluster-buffer extraction)",
  "gee_trmm_Jul_25km",          "rain_july_25km",          "July precipitation (wet-season peak in much of the region), 25 km buffer", "mm/month", "TRMM precipitation (GEE cluster-buffer extraction)",
  "gee_trmm_10km_seasonal_sd",  "rain_seasonality_sd_10km","Within-year standard deviation of monthly precipitation (rainfall seasonality), 10 km buffer", "mm/month", "TRMM precipitation (GEE cluster-buffer extraction)",
  "gee_temp_May_25km",          "temp_may_25km",           "May land-surface temperature, 25 km buffer", "deg C", "MODIS land-surface temperature (GEE cluster-buffer extraction)",
  "gee_temp_Oct_10km",          "temp_oct_10km",           "October land-surface temperature, 10 km buffer", "deg C", "MODIS land-surface temperature (GEE cluster-buffer extraction)",
  "gee_temp_10km_seasonal_sd",  "temp_seasonality_sd_10km","Within-year standard deviation of monthly land-surface temperature (temperature seasonality), 10 km buffer", "deg C", "MODIS land-surface temperature (GEE cluster-buffer extraction)",
  "gee_temp_25km_peak_month",   "temp_peak_month_25km",    "Calendar month (1-12) of maximum land-surface temperature, 25 km buffer", "month index 1-12", "MODIS land-surface temperature (GEE cluster-buffer extraction)",
  "gee_ndvi_Apr_50km",          "ndvi_april_50km",         "April NDVI (vegetation greenness), 50 km buffer", "NDVI index (-1 to 1)", "NDVI (GEE cluster-buffer extraction)",
  "gee_ndvi_Nov_25km",          "ndvi_november_25km",      "November NDVI (vegetation greenness), 25 km buffer", "NDVI index (-1 to 1)", "NDVI (GEE cluster-buffer extraction)",
  "gee_ndvi_Dec_25km",          "ndvi_december_25km",      "December NDVI (vegetation greenness), 25 km buffer", "NDVI index (-1 to 1)", "NDVI (GEE cluster-buffer extraction)",
  "gee_globalhumanmodification_2016", "human_modification_index", "Global Human Modification index: cumulative human pressure on the landscape (built-up, agriculture, transport, etc.), 2016", "0-1 (0 = wild, 1 = fully modified)", "CSP Global Human Modification (GEE)",
  "gee_grassland_50km",         "grassland_fraction_50km", "Grassland land-cover fraction, 50 km buffer", "fraction 0-1", "GPW grasslands (GEE cluster-buffer extraction)",
  "gee_built_surface_nres_10km","built_surface_nonres_10km","Non-residential built-up surface (GHSL), 10 km buffer", "m^2 built per cell", "GHSL Built-Up surface (GEE cluster-buffer extraction)",
  "gee_fpp_1km_10km",           "pop_density_cropland_10km","Population density on cropland (1 km source), 10 km buffer", "persons/km^2", "GHSL population x cropland mask (GEE cluster-buffer extraction)",
  "MAP_202106_Antimalarial_Effective_Treatment", "malaria_effective_treatment", "Proportion of malaria cases receiving effective antimalarial treatment (Malaria Atlas Project, 2021.06 release)", "proportion 0-1", "Malaria Atlas Project (cluster-GPS raster extraction)",
  "MAP_202202_Pf_Reproductive_Number", "malaria_pf_repro_number", "P. falciparum reproductive number under control (malaria transmission intensity; MAP 2022.02 release)", "reproductive number", "Malaria Atlas Project (cluster-GPS raster extraction)"
)
PRED_RAW   <- pred_meta$raw
PRED_CLEAN <- pred_meta$clean

# ── 1. Pooled individual-level long frame (country x outcome) ────────────────
recs <- list()
for (i in seq_len(nrow(COUNTRIES))) {
  cr  <- COUNTRIES[i, ]
  cfg <- get_country_configs()[[cr$key]]
  clid <- cfg$cluster_id; wcol <- cfg$weight_col
  for (ocn in names(cfg$outcomes)) {
    oc <- cfg$outcomes[[ocn]]; yb <- oc$binary
    p  <- file.path(STORE, paste0("outcome_data_", cr$low, "_", ocn))
    if (!file.exists(p)) next
    d  <- tryCatch(readRDS(p)$data, error = function(e) NULL)
    if (is.null(d) || is.null(yb) || !yb %in% colnames(d)) next
    df <- data.frame(
      country    = cr$label,
      admin1     = trimws(as.character(d$Admin1)),
      admin2     = trimws(as.character(d$Admin2)),
      cluster_id = as.character(d[[clid]]),
      weight     = if (wcol %in% colnames(d)) to_num(d[[wcol]]) else 1,
      outcome    = ocn,
      y          = to_num(d[[yb]]),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(PRED_RAW)) {
      v <- PRED_RAW[j]
      df[[PRED_CLEAN[j]]] <- if (v %in% colnames(d)) to_num(d[[v]]) else NA_real_
    }
    df <- df[is.finite(df$y), , drop = FALSE]
    if (nrow(df) > 0) recs[[paste0(cr$low, "_", ocn)]] <- df
  }
}
indiv <- dplyr::bind_rows(recs)
present_oc <- intersect(OUTCOME_ORDER, unique(indiv$outcome))
cat(sprintf("Pooled %d individual records; outcomes present: %s\n",
            nrow(indiv), paste(present_oc, collapse = ", ")))

# ── 2. Aggregator: wide table at a given level ───────────────────────────────
wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) mean(x, na.rm = TRUE) else sum(x[ok] * w[ok]) / sum(w[ok])
}

build_level <- function(df, keys, deff, with_nclusters) {
  # Predictors (cluster-constant): dedup to cluster, then mean over the area's
  # clusters. Outcome-agnostic — one predictor value per area, shared by all
  # outcomes.
  clu <- df %>%
    group_by(across(all_of(c("country", "admin1", "admin2", "cluster_id")))) %>%
    summarise(across(all_of(PRED_CLEAN), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  preds <- clu %>%
    group_by(across(all_of(keys))) %>%
    summarise(n_clusters = dplyr::n_distinct(cluster_id),
              across(all_of(PRED_CLEAN), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  # Per-outcome aggregates at this level.
  oc_agg <- df %>%
    group_by(across(all_of(c(keys, "outcome")))) %>%
    summarise(prev    = wmean(y, weight),
              n       = sum(is.finite(y)),
              ndef    = sum(y == 1, na.rm = TRUE),
              .groups = "drop")
  eff_n <- pmax(oc_agg$n / deff, 1); pp <- clip01(oc_agg$prev)
  oc_agg$var <- pp * (1 - pp) / eff_n
  ocw <- oc_agg %>%
    pivot_wider(id_cols = all_of(keys), names_from = outcome,
                values_from = c(prev, n, ndef, var),
                names_glue = "{.value}_{outcome}")

  out <- dplyr::left_join(preds, ocw, by = keys)
  for (v in PRED_CLEAN) out[[v]][!is.finite(out[[v]])] <- NA_real_

  # Column order: ids -> [n_clusters] -> predictors -> per-outcome blocks.
  oc_cols <- unlist(lapply(present_oc, function(o)
    paste0(c("prev_", "n_", "ndef_", "var_"), o)))
  oc_cols <- intersect(oc_cols, colnames(out))
  extra   <- if (with_nclusters) "n_clusters" else character(0)
  out[, c(keys, extra, PRED_CLEAN, oc_cols)]
}

cluster_df <- build_level(indiv, c("country", "admin1", "admin2", "cluster_id"),
                          deff = 1,         with_nclusters = FALSE)
admin2_df  <- build_level(indiv, c("country", "admin1", "admin2"),
                          deff = DEFF_AREA, with_nclusters = TRUE)
admin1_df  <- build_level(indiv, c("country", "admin1"),
                          deff = DEFF_AREA, with_nclusters = TRUE)

p1 <- file.path(OUT_DIR, "mn_cluster.csv")
p2 <- file.path(OUT_DIR, "mn_admin2.csv")
p3 <- file.path(OUT_DIR, "mn_admin1.csv")
write_csv_atomic(cluster_df, p1)
write_csv_atomic(admin2_df,  p2)
write_csv_atomic(admin1_df,  p3)
cat(sprintf("\nWrote:\n  %s  (%d x %d)\n  %s  (%d x %d)\n  %s  (%d x %d)\n",
            basename(p1), nrow(cluster_df), ncol(cluster_df),
            basename(p2), nrow(admin2_df),  ncol(admin2_df),
            basename(p3), nrow(admin1_df),  ncol(admin1_df)))
cat("Rows per country (admin2):\n"); print(table(admin2_df$country, dnn = NULL))
cat("Outcome availability (admin2, # areas with a non-NA prevalence):\n")
for (o in present_oc)
  cat(sprintf("  %-13s %d\n", o, sum(!is.na(admin2_df[[paste0("prev_", o)]]))))

# ── 3. Data dictionary ───────────────────────────────────────────────────────
id_rows <- tibble::tribble(
  ~variable,    ~description, ~type, ~units, ~source_derivation, ~role,
  "country",    "Study country", "character", "", "Survey metadata", "id",
  "admin1",     "First-level administrative unit (region/province)", "character", "", "Survey cluster geolocation -> GADM Admin-1", "id",
  "admin2",     "Second-level administrative unit (district)", "character", "", "Survey cluster geolocation -> GADM Admin-2", "id",
  "cluster_id", "Survey cluster (primary sampling unit) identifier; cluster file only", "character", "", "Source biomarker survey (gw_cnum)", "id",
  "n_clusters", "Number of distinct survey clusters contributing (Admin-2 / Admin-1 files only)", "integer", "count", "Distinct cluster_id within the area", "id"
)
# One block of 4 rows per outcome (prev / n / ndef / var).
oc_rows <- do.call(rbind, lapply(present_oc, function(o) {
  lab <- OUTCOME_LABELS[[o]]
  data.frame(
    variable = paste0(c("prev_", "n_", "ndef_", "var_"), o),
    description = c(
      sprintf("Survey-weighted prevalence of %s (the PRIMARY outcome to model for this nutrient; proportion 0-1)", lab),
      sprintf("Number of surveyed individuals with a non-missing outcome for %s (denominator)", lab),
      sprintf("Number of those classified deficient for %s (numerator; unweighted prevalence = ndef/n)", lab),
      sprintf("Design-based sampling variance of prev_%s: p(1-p)/effective_n, effective_n = n/deff (deff=1.5 Admin-2/Admin-1, 1 at cluster). Use as the known sampling-error variance, or 1/var as a precision weight.", o)
    ),
    type = c("numeric", "integer", "integer", "numeric"),
    units = c("proportion 0-1", "count", "count", "variance of a proportion"),
    source_derivation = c("Survey-weighted mean of the binary flag",
                          "Count of survey respondents", "Sum of the binary flag",
                          "Derived"),
    role = c("outcome", "outcome_denominator", "outcome_denominator", "weight"),
    stringsAsFactors = FALSE
  )
}))
pred_rows <- data.frame(
  variable = pred_meta$clean, description = pred_meta$description,
  type = "numeric", units = pred_meta$units,
  source_derivation = paste0(pred_meta$source,
    "; aggregated to this level as the mean over the area's survey clusters"),
  role = "predictor", stringsAsFactors = FALSE
)
dict <- dplyr::bind_rows(id_rows, oc_rows, pred_rows)
pd <- file.path(SUB, "data_dictionary.csv")
write_csv_atomic(dict, pd)
cat(sprintf("\nWrote data dictionary: %s (%d variables)\n", basename(pd), nrow(dict)))
cat("\nBuild complete.\n")
