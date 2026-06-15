# =============================================================================
# global.R — loaded once at app startup, shared across all sessions
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(leaflet)
  library(plotly)
  library(reactable)
  library(htmltools)
})

# ── Data paths ─────────────────────────────────────────────────────────────
DATA_DIR <- "data"

# ── Load all dashboard data once ──────────────────────────────────────────
admin2_pred <- readRDS(file.path(DATA_DIR, "admin2_predictions.rds"))
admin2_pop  <- readRDS(file.path(DATA_DIR, "admin2_population.rds"))

# Full-coverage area-level SAE predictions (every district, surveyed +
# unsurveyed) — the alternate "Area-level SAE" map layer. NULL if not built.
admin2_area_pred <- {
  .p <- file.path(DATA_DIR, "admin2_area_predictions.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}
admin2_bnds <- readRDS(file.path(DATA_DIR, "admin2_boundaries.rds"))
admin1_bnds <- readRDS(file.path(DATA_DIR, "admin1_boundaries.rds"))
natl_est    <- readRDS(file.path(DATA_DIR, "national_estimates.rds"))
cv_perf     <- readRDS(file.path(DATA_DIR, "cv_performance.rds"))
meta        <- readRDS(file.path(DATA_DIR, "metadata.rds"))

# GBD comparison data (placeholder by default; actual GBD data if available)
gbd_obj <- readRDS(file.path(DATA_DIR, "gbd_estimates.rds"))
gbd_data <- list(data = gbd_obj$data)
gbd_meta <- gbd_obj$meta

# Côte d'Ivoire OOS predictions (may be empty if pipeline not yet rebuilt)
oos_civ_path <- file.path(DATA_DIR, "oos_cote_divoire.rds")
oos_civ <- if (file.exists(oos_civ_path)) {
  readRDS(oos_civ_path)
} else {
  list(predictions = data.frame(), boundaries = NULL)
}

# Transportability (leave-one-country-out) Admin-2 predictions, used by the
# map explorer's "transportability error" difference layer. NULL if the
# data-prep step has not produced it yet (layer then shows as "no data").
loco_path <- file.path(DATA_DIR, "transportability_loco.rds")
loco_pred <- if (file.exists(loco_path)) readRDS(loco_path) else NULL

# Transport-calibration (predicted vs true prevalence under LOCO, multi-level
# and multi-approach). NULL components if the data-prep step has not produced
# the source tables yet — the module then renders an empty state.
transport_cal_path <- file.path(DATA_DIR, "transport_calibration.rds")
transport_cal <- if (file.exists(transport_cal_path)) {
  readRDS(transport_cal_path)
} else {
  list(nat_indiv = NULL, nat_area = NULL, adm2_area = NULL, build_time = NULL)
}

# Cluster-resolution comparison (admin-2 vs survey-cluster). NULL components if
# the data-prep step has not produced the cluster CSVs yet — the module then
# renders an empty state rather than erroring.
cluster_res_path <- file.path(DATA_DIR, "cluster_resolution.rds")
cluster_res <- if (file.exists(cluster_res_path)) {
  readRDS(cluster_res_path)
} else {
  list(comparison = NULL, loco = NULL, build_time = NULL)
}

# Importance / SHAP / domain ablation
importance_path <- file.path(DATA_DIR, "importance.rds")
importance_data <- if (file.exists(importance_path)) {
  readRDS(importance_path)
} else {
  list(shap = NULL, varimp = NULL, ablation = NULL)
}

# Model diagnostics: ROC / PR / calibration curves + binary/continuous metrics
# and the Platt-recalibrated binary metrics. NULL components if the data-prep
# step has not produced the source CSVs yet — modules render an empty state.
model_diag_path <- file.path(DATA_DIR, "model_diagnostics.rds")
model_diag <- if (file.exists(model_diag_path)) {
  readRDS(model_diag_path)
} else {
  list(binary = NULL, continuous = NULL, calibrated = NULL,
       roc = NULL, pr = NULL, calibration = NULL)
}

# Methods comparison (corrected vs production): the parallel _targets_corrected
# pipeline's P1-P8 head-to-head bundle. NULL if the corrected pipeline has not
# been run yet — the module then renders empty states.
methods_comp_path <- file.path(DATA_DIR, "methods_comparison.rds")
methods_comp <- if (file.exists(methods_comp_path)) {
  readRDS(methods_comp_path)
} else {
  list(cv_compare = NULL, calibration = NULL, admin2_error = NULL,
       decision = NULL, trust = NULL, area_pp = NULL, interval_summary = NULL)
}

# Method benchmarks: SuperLearner vs small-area-estimation methods, individual-
# vs area-level SL, per-district accuracy, prescreened SL LOCO, and enriched
# area-transport LOCO metrics. NULL components -> empty state in the module.
benchmarks_path <- file.path(DATA_DIR, "benchmarks.rds")
benchmarks_data <- if (file.exists(benchmarks_path)) {
  readRDS(benchmarks_path)
} else {
  list(benchmarks = NULL, area_comparison = NULL, admin2_error = NULL,
       sl_prescreened = NULL, area_transport = NULL)
}

# Useful operator
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Build version stamp from data file modification time
data_build_time <- format(file.info(file.path(DATA_DIR, "metadata.rds"))$mtime,
                            "%Y-%m-%d %H:%M")

# ── Source helpers and modules ────────────────────────────────────────────
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

# ── Convenience lookups ───────────────────────────────────────────────────
country_choices <- setNames(names(meta$countries), meta$countries)
outcome_choices <- setNames(names(meta$outcome_labels), meta$outcome_labels)

# Color palette for WHO public health classes (color-blind safe)
who_colors <- c(
  "Low"      = "#2c7bb6",  # blue
  "Mild"     = "#abd9e9",  # light blue
  "Moderate" = "#fdae61",  # orange
  "Severe"   = "#d7191c",  # red
  "No data"  = "#cccccc"   # gray
)
