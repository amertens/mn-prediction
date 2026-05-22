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

# Importance / SHAP / domain ablation
importance_path <- file.path(DATA_DIR, "importance.rds")
importance_data <- if (file.exists(importance_path)) {
  readRDS(importance_path)
} else {
  list(shap = NULL, varimp = NULL, ablation = NULL)
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
