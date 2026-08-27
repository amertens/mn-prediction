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
# layer (see archive/sandbox/sae_sl_hybrid_prototype.R). NULL if not built.
admin2_bym2_pred <- {
  .p <- file.path(DATA_DIR, "admin2_bym2_predictions.rds")
  if (file.exists(.p)) readRDS(.p) else NULL
}

# Full-coverage AREA-LEVEL RECIPE layer (docs/AREA_LEVEL_RECIPE_SPEC.md): area
# prevalence under a binomial loss with a leakage-free near-null pre-filter,
# trained on surveyed districts and applied to every polygon (universal/gee
# features). The recommended primary district estimator. NULL if not built; build
# with dashboard/data-raw/build_recipe_predictions.R. (Intervals: TODO — the SAE
# layer is for interval width only; point level from national-anchor calibration.)
admin2_recipe_pred <- {
  .p <- file.path(DATA_DIR, "admin2_recipe_predictions.rds")
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

# ── Which model produces the numbers, defined once ─────────────────────────
# Every tab reads DEFAULT_PRED_MODEL, so the same district cannot show one
# figure on the map and a different one on its profile page.
#
# The default is the district prevalence model (`recipe`).
# docs/AREA_LEVEL_RECIPE_SPEC.md settles this: modelling each district's
# prevalence directly beats modelling individual people and averaging up
# (district-ranking correlation 0.06 -> 0.5-0.8 on the outcomes that carry
# signal). The same spec is explicit that the small-area models are there for
# uncertainty ranges, not for the level -- the covariates rank districts well
# but do not pin down how high prevalence is, and shrinking the level toward
# the covariate mean makes absolute error worse. So Fay-Herriot and BYM2 stay
# available, and stay off by default.
DEFAULT_PRED_MODEL <- "recipe"

PRED_MODEL_INFO <- list(
  recipe = list(
    label = "Recommended (all districts)",
    blurb = paste("Estimates each district's prevalence directly from satellite,",
                  "climate and geospatial indicators. Best at telling you which",
                  "districts are worst off. Read the percentage as a planning",
                  "figure rather than a measurement. No uncertainty range yet.")),
  bym2 = list(
    label = "Uncertainty ranges, spatially smoothed",
    blurb = paste("Adds a 95% range to every district by borrowing information",
                  "from its neighbours. Narrow where a survey exists, wider",
                  "where none does. Use it to see how firm a number is.")),
  fh = list(
    label = "Uncertainty ranges, blended with the survey",
    blurb = paste("Blends each district's own survey result, where there is one,",
                  "with the indicator-based estimate, and gives every district a",
                  "95% range. Ranges run wider than the spatial model's.")),
  area = list(
    label = "Alternative fit, for comparison",
    blurb = paste("A different machine-learning fit over the same indicators.",
                  "Kept so you can check whether a district's standing depends",
                  "on the method. Point estimates only.")),
  sl = list(
    label = "Person-level fit (surveyed districts only)",
    blurb = paste("Models individual people and averages up to the district.",
                  "Covers only districts that were surveyed, so most of the map",
                  "is blank. Kept as a sensitivity check."))
)

# Which layers actually carry an interval — the map's confidence views and the
# district profile's error bars need this, and two layers have none.
PRED_MODEL_HAS_CI <- c(recipe = FALSE, bym2 = TRUE, fh = TRUE,
                       area = FALSE, sl = TRUE)

#' The prediction table behind a model key, falling back to the person-level
#' table if that layer was not built.
pred_model_data <- function(key) {
  d <- switch(key %||% DEFAULT_PRED_MODEL,
              recipe = admin2_recipe_pred, bym2 = admin2_bym2_pred,
              fh     = admin2_fh_pred,     area = admin2_area_pred,
              sl     = admin2_pred,        NULL)
  if (is.null(d)) admin2_pred else d
}

#' Named choices for a model selector, in the order given, dropping any layer
#' that was not built so the app never offers an empty option.
pred_model_choices <- function(keys = names(PRED_MODEL_INFO)) {
  keys <- Filter(function(k) {
    d <- switch(k, recipe = admin2_recipe_pred, bym2 = admin2_bym2_pred,
                fh = admin2_fh_pred, area = admin2_area_pred, sl = admin2_pred)
    !is.null(d) && nrow(d) > 0
  }, keys)
  setNames(keys, vapply(keys, function(k) PRED_MODEL_INFO[[k]]$label, character(1)))
}

pred_model_has_ci <- function(key) {
  isTRUE(unname(PRED_MODEL_HAS_CI[key %||% DEFAULT_PRED_MODEL]))
}

#' Does this country x outcome actually vary between districts?
#'
#' For some cells the model finds no district-level signal and every district
#' comes back with the same number -- the near-null pre-filter drops every
#' candidate predictor and the fit reduces to an intercept. Measured on the
#' 2026-08-27 build, 4 of 24 cells are exactly constant (Ghana women's vitamin A
#' is 1.53% in all 260 districts; Sierra Leone child vitamin A 11.93% in all 14)
#' and 2 more vary by under a percentage point.
#'
#' That is a legitimate result -- it says the proxies carry nothing for this
#' nutrient here -- but a choropleth drawn over it implies variation that is not
#' there, and a "worst districts" table becomes an arbitrary ordering of
#' identical values. Callers should say so instead of drawing the map.
#'
#' @param prev numeric vector of district estimates
#' @param min_spread minimum max-min, as a proportion, to count as informative
#' @return TRUE if the estimates carry usable district-level variation
has_district_signal <- function(prev, min_spread = 0.01) {
  v <- prev[is.finite(prev)]
  length(unique(round(v, 6))) > 2 && diff(range(v)) >= min_spread
}

#' Give a district table a 95% range when its own model does not carry one.
#'
#' The recommended model produces a level but no interval yet. Rather than
#' present its numbers as though they were exact, borrow the spatial model's
#' interval WIDTH for the same district and centre it on the recommended
#' estimate. This is the split docs/AREA_LEVEL_RECIPE_SPEC.md prescribes: the
#' covariates rank districts, the small-area layer supplies the uncertainty.
#'
#' Districts with no interval get NA, and callers must cope with that.
#'
#' Fay-Herriot is the source rather than BYM2. Both are sane as rebuilt on
#' 2026-08-27 -- FH sits 4x or further from the national survey figure on 1 of
#' 24 country x outcome cells, BYM2 on 2 -- but FH is the marginally cleaner of
#' the two and its intervals are narrower (mean width 14.8 pp against 22.9),
#' so it is the less generous choice and the one that fails more quietly.
#'
#' @param df data.frame with Admin2 and pred_prev
#' @param country_label display label ("Ghana"), matching the layer's `country`
#' @param oc outcome tag
#' @return `df` with prev_lo / prev_hi added
attach_prevalence_range <- function(df, country_label, oc) {
  df$prev_lo <- NA_real_; df$prev_hi <- NA_real_
  src <- admin2_fh_pred %||% admin2_bym2_pred
  if (is.null(src)) return(df)
  s <- src[src$country == country_label & src$outcome == oc, , drop = FALSE]
  if (!nrow(s) || !all(c("ci_lo", "ci_hi", "pred_prev") %in% names(s))) return(df)
  s <- s[!duplicated(s$Admin2), ]
  i <- match(df$Admin2, s$Admin2)
  lo_gap <- s$pred_prev[i] - s$ci_lo[i]
  hi_gap <- s$ci_hi[i]  - s$pred_prev[i]
  df$prev_lo <- pmax(0, df$pred_prev - lo_gap)
  df$prev_hi <- pmin(1, df$pred_prev + hi_gap)
  df
}

#' The "how these are built" note, shared by the map and the district profiles
#' so the two cannot drift apart.
pred_model_help <- function(keys = names(PRED_MODEL_INFO)) {
  keys <- intersect(keys, names(PRED_MODEL_INFO))
  tags$details(
    style = "font-size:0.82em; color:#555; margin:-2px 0 6px;",
    tags$summary("What do these options mean?"),
    tags$ul(
      style = "padding-left:1.1em; margin-bottom:0;",
      lapply(keys, function(k) tags$li(
        strong(PRED_MODEL_INFO[[k]]$label, ": "), PRED_MODEL_INFO[[k]]$blurb))
    )
  )
}

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

# ── Site banner ────────────────────────────────────────────────────────────
# This used to read "not for citation or external distribution" on every page.
# The app is served without a login, so that was either untrue or an instruction
# nobody could follow — and a caveat people learn to scroll past is worse than
# none, because it takes the real warnings down with it. What replaces it says
# what the estimates will and will not carry, which is the thing a reader
# actually needs before using a number.
# Held as plain strings, not just as HTML, because the PDF/HTML review reports
# in dashboard/report/ print the same wording. One source, so the app and the
# printed brief cannot end up saying different things about what these numbers
# will carry.
SITE_SCOPE_HEADLINE <- "Preliminary working estimates, not official statistics."
SITE_SCOPE_BODY <- paste(
  "Reliable enough to rank districts by need, but not precise enough to quote",
  "as a measured prevalence.")
SITE_SCOPE_POINTER <- "Decision value shows where they are solid enough to act on."

site_banner <- div(
  class = "alert alert-info",
  style = paste("margin:0 0 10px;border-radius:0;text-align:center;",
                "font-size:0.88em;padding:6px 12px;"),
  bsicons::bs_icon("info-circle"), " ",
  strong(SITE_SCOPE_HEADLINE), " ", SITE_SCOPE_BODY,
  tags$span(style = "color:#4a6b7c;", " ", SITE_SCOPE_POINTER)
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
