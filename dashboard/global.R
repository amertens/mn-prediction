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

# Full-coverage Fay-Herriot / empirical-Bayes layer (every district) WITH
# design-based 95% intervals (tight for surveyed, wide for unsurveyed). NULL if
# not built.
admin2_fh_pred <- {
  .p <- file.path(DATA_DIR, "admin2_fh_predictions.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}

# Full-coverage SL -> BYM2 layer (every district) WITH spatially-smoothed 95%
# credible intervals. Tighter, better-calibrated intervals than the Fay-Herriot
# layer (see sandbox/sae_sl_hybrid_prototype.R). NULL if not built.
admin2_bym2_pred <- {
  .p <- file.path(DATA_DIR, "admin2_bym2_predictions.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}
admin2_bnds <- readRDS(file.path(DATA_DIR, "admin2_boundaries.rds"))
admin1_bnds <- readRDS(file.path(DATA_DIR, "admin1_boundaries.rds"))

# Admin-3 (chiefdom) layer — currently Sierra Leone only (153 chiefdoms vs 14
# districts), built by data-raw/_build_sl_admin3.R via a Fay-Herriot/empirical-
# Bayes area model on chiefdom GEE covariates. NULL components if not built; the
# map explorer then simply does not offer the Admin-3 level.
admin3_pred <- {
  .p <- file.path(DATA_DIR, "admin3_predictions.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}
admin3_bnds <- {
  .p <- file.path(DATA_DIR, "admin3_boundaries.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}
# Countries that have an Admin-3 layer available (lowercased boundary-list keys).
admin3_countries <- if (!is.null(admin3_bnds)) names(admin3_bnds) else character(0)
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

# Per-outcome biomarker / data-quality caveats (surfaced in the map, district
# profiles and Methods tab). Keyed by outcome tag; NULL = no specific caveat.
biomarker_caveats <- list(
  women_vitA   = paste("Vitamin A in women is measured by retinol-binding protein",
                       "(RBP), a weaker proxy for vitamin A status in women that is",
                       "affected by inflammation. Treat women's vitamin A estimates",
                       "with extra caution."),
  child_vitA   = paste("Vitamin A is measured by RBP, which is influenced by",
                       "inflammation; inflammation-adjustment is applied but some",
                       "residual bias is possible."),
  women_b12    = paste("B12 is measured by serum B12; holo-transcobalamin (active",
                       "B12) is a more specific marker that is still being explored.",
                       "Interpret B12 estimates as indicative."),
  women_folate = paste("Folate deficiency is often low-prevalence, which is harder",
                       "to estimate precisely; small district differences may not",
                       "be reliable."),
  child_zinc   = paste("Zinc biomarkers and thresholds are less standardized, so",
                       "zinc estimates carry more uncertainty than iron or vitamin A."),
  women_zinc   = paste("Zinc biomarkers and thresholds are less standardized, so",
                       "zinc estimates carry more uncertainty than iron or vitamin A.")
)
GENERAL_CAVEAT <- paste(
  "Surveys span different years (2013–2018) and laboratories, so biomarker",
  "values are not perfectly comparable across countries. Predictions transported",
  "to a country with no survey are the least reliable and may be no better than",
  "chance for some outcomes — see the Transportability and Decision value tabs."
)

# ── Shared popover content (used by both the internal and public apps) ──────
about_content <- div(
  h5("About this dashboard", style = "margin-top: 0;"),
  p("This dashboard visualizes sub-national micronutrient deficiency ",
    "predictions for The Gambia, Ghana, Sierra Leone, and Malawi, with ",
    "out-of-sample predictions for Côte d'Ivoire. It is designed for ",
    "ministries of health, funders, and researchers planning interventions ",
    "or surveys."),
  h6("Methodology"),
  p("Predictions come from SuperLearner ensemble machine learning models ",
    "trained on individual-level biomarker data linked to ~1,000 routinely ",
    "available proxy indicators (satellite imagery, climate reanalysis, ",
    "modeled disease burden, food security and pricing, and household ",
    "surveys). Uncertainty is shown as an approximate 95% range around each ",
    "prediction; in our own out-of-sample checks these ranges covered the ",
    "true value about 90% of the time."),
  h6("Citation"),
  p(em("Mertens et al. (in prep.). Sub-national micronutrient deficiency prediction ",
        "from proxy data: a multi-country machine learning framework.")),
  h6("Data sources"),
  p("Biomarker surveys: GMNS (Gambia 2018), GMS (Ghana 2017), SLMS (Sierra ",
    "Leone 2013), MNS (Malawi 2015–16). Proxy data: CHIRPS, WorldPop, MAP, ",
    "SoilGrids, GDL, WFP HungerMap, WFP food prices, IPC/CH, ACLED, ",
    "HarvestStat, DHS Admin-2 indicators, IHME GBD."),
  h6("Limitations"),
  tags$ul(
    tags$li("Single survey timepoint per country; predictions assume the ",
            "proxy–outcome relationship is stable over time."),
    tags$li("Cross-country generalizability is limited (LOCO AUC 0.50–0.73)."),
    tags$li("Population denominators carry their own modelling uncertainty."),
    tags$li("GBD comparison data is currently placeholder; an RA task is ",
            "open to source actual GBD Results Tool exports.")
  ),
  p(em(sprintf("Data build: %s", data_build_time)),
    style = "color: #888; font-size: 0.85em;")
)

glossary_content <- div(
  h5("Glossary", style = "margin-top: 0;"),
  tags$dl(
    tags$dt("AUC (ROC-AUC)"),
    tags$dd("Area under the receiver-operating characteristic curve. ",
            "Measures discrimination: 0.5 = random, 1.0 = perfect. ",
            "Above 0.7 is considered fair, above 0.8 good."),
    tags$dt("Brier score"),
    tags$dd("Mean squared error between predicted probability and observed ",
            "outcome. Lower is better; depends on outcome prevalence."),
    tags$dt("Prediction interval (conformal)"),
    tags$dd("An approximate 95% range around a prediction, built from the ",
            "model's past errors rather than an assumed bell curve. In our ",
            "checks it covered the true value about 90% of the time, so read ",
            "it as indicative, not exact."),
    tags$dt("Percentage points (pp)"),
    tags$dd("The plain gap between two percentages. Moving from 20% to 23% ",
            "is a 3 percentage-point (3 pp) increase, not a 3% increase."),
    tags$dt("District and region (Admin-2 / Admin-1)"),
    tags$dd("\"District\" is the second administrative level (Admin-2); ",
            "\"region\" is the first (Admin-1). Predictions are made for districts."),
    tags$dt("Cluster-blocked cross-validation"),
    tags$dd("Cross-validation where individuals from the same survey ",
            "cluster are always assigned to the same fold. Prevents ",
            "optimistic estimates from spatial autocorrelation."),
    tags$dt("LOCO cross-validation"),
    tags$dd("Leave-one-country-out: train on three countries, test on the ",
            "fourth. Estimates how well the model would perform if applied ",
            "to a country with no biomarker data."),
    tags$dt("Proxy-only model"),
    tags$dd("A model that excludes all variables from the original ",
            "biomarker survey as predictors, ensuring it can be applied ",
            "to areas where no survey exists."),
    tags$dt("WHO public health classification"),
    tags$dd("Standard severity tiers based on prevalence: Low / Mild / ",
            "Moderate / Severe. Thresholds are well-established for vitamin A ",
            "and iron deficiency; less standardized for B12 and zinc."),
    tags$dt("Hidden burden"),
    tags$dd("Population living in districts where local prevalence is clearly ",
            "different from the country average. It measures the cost of ",
            "relying on national-level estimates for sub-national targeting.")
  )
)
