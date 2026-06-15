# =============================================================================
# Micronutrient Burden Dashboard
# =============================================================================
# Interactive dashboard for exploring sub-national micronutrient deficiency
# predictions across The Gambia, Ghana, Sierra Leone, and Malawi, with
# out-of-sample predictions for Côte d'Ivoire.
#
# To run locally:
#   setwd("dashboard"); shiny::runApp()
# To deploy:
#   Rscript dashboard/deploy.R
# =============================================================================

source("global.R")

# ── About / citation popover content ──────────────────────────────────────
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

# ── Glossary popover content ──────────────────────────────────────────────
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

# ── UI ─────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "Micronutrient Burden",
  theme = bs_theme(
    version = 5,
    bootswatch = "cosmo",
    primary = "#2c7bb6",
    base_font = font_google("Source Sans Pro")
  ),
  # Only map-heavy tabs should fill the viewport. Tabs with stacked cards +
  # methods notes need to size to their content — otherwise each card is
  # cropped to a fraction of viewport height and shows a scrollbar.
  fillable = c("Map explorer", "Côte d'Ivoire (OOS)"),
  underline = TRUE,

  nav_panel(
    title = "Map explorer",
    icon  = bsicons::bs_icon("map"),
    mod_map_explorer_ui("map")
  ),

  nav_panel(
    title = "District profiles",
    icon  = bsicons::bs_icon("file-earmark-medical"),
    mod_district_profile_ui("district")
  ),

  nav_panel(
    title = "National burden",
    icon  = bsicons::bs_icon("graph-up"),
    mod_national_burden_ui("burden")
  ),

  nav_panel(
    title = "Decision value",
    icon  = bsicons::bs_icon("bullseye"),
    mod_decision_value_ui("decision")
  ),

  nav_panel(
    title = "Scenarios",
    icon  = bsicons::bs_icon("sliders"),
    mod_scenarios_ui("scenarios")
  ),

  nav_panel(
    title = "Importance",
    icon  = bsicons::bs_icon("bar-chart-line"),
    mod_importance_ui("importance")
  ),

  nav_panel(
    title = "Model diagnostics",
    icon  = bsicons::bs_icon("clipboard-data"),
    mod_diagnostics_ui("diagnostics")
  ),

  nav_panel(
    title = "Benchmarks",
    icon  = bsicons::bs_icon("trophy"),
    mod_benchmarks_ui("benchmarks")
  ),

  nav_panel(
    title = "Methods (corrected)",
    icon  = bsicons::bs_icon("shield-check"),
    mod_methods_comparison_ui("methods_corrected")
  ),

  nav_panel(
    title = "Resolution",
    icon  = bsicons::bs_icon("bullseye"),
    mod_resolution_ui("resolution")
  ),

  nav_panel(
    title = "Transportability",
    icon  = bsicons::bs_icon("globe-americas"),
    mod_transport_calibration_ui("transport")
  ),

  nav_panel(
    title = "Côte d'Ivoire (OOS)",
    icon  = bsicons::bs_icon("compass"),
    mod_oos_civ_ui("oos_civ")
  ),

  nav_panel(
    title = "GBD comparison",
    icon  = bsicons::bs_icon("globe2"),
    mod_gbd_compare_ui("gbd")
  ),

  nav_panel(
    title = "Methods",
    icon  = bsicons::bs_icon("info-circle"),
    mod_methods_ui("methods")
  ),

  nav_spacer(),

  nav_item(
    popover(
      tags$button(class = "btn btn-link",
                  bsicons::bs_icon("book"), " Glossary"),
      glossary_content,
      title = NULL,
      placement = "bottom",
      options = list(html = TRUE, container = "body")
    )
  ),

  nav_item(
    popover(
      tags$button(class = "btn btn-link",
                  bsicons::bs_icon("info-circle"), " About"),
      about_content,
      title = NULL,
      placement = "bottom",
      options = list(html = TRUE, container = "body")
    )
  ),

  footer = div(
    style = "text-align: center; color: #888; font-size: 0.8em; padding: 8px;",
    sprintf("Micronutrient Burden Dashboard | Data built %s | ",
            data_build_time),
    tags$a(href = "https://github.com/", target = "_blank", "Source")
  )
)

# ── Server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  mod_map_explorer_server("map")
  mod_district_profile_server("district")
  mod_national_burden_server("burden")
  mod_decision_value_server("decision")
  mod_scenarios_server("scenarios")
  mod_importance_server("importance")
  mod_diagnostics_server("diagnostics")
  mod_benchmarks_server("benchmarks")
  mod_methods_comparison_server("methods_corrected")
  mod_resolution_server("resolution")
  mod_transport_calibration_server("transport")
  mod_oos_civ_server("oos_civ")
  mod_gbd_compare_server("gbd")
  mod_methods_server("methods")
}

shinyApp(ui, server)
