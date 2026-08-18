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

# About/Glossary popover content is defined in global.R (shared with app_public.R).

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

  # Shown on every tab: in-development banner + a faint fixed watermark.
  header = tagList(
    tags$style(HTML(
      ".indev-watermark{position:fixed;bottom:10px;right:12px;z-index:1030;
       font-weight:700;font-size:12px;letter-spacing:1px;color:rgba(200,0,0,0.33);
       border:2px solid rgba(200,0,0,0.28);border-radius:5px;padding:2px 8px;
       transform:rotate(-4deg);pointer-events:none;background:rgba(255,255,255,0.45);}")),
    div(class = "alert alert-warning",
        style = "margin:0 0 10px;border-radius:0;text-align:center;font-size:0.88em;padding:6px 10px;",
        bsicons::bs_icon("cone-striped"), " ",
        strong("In development"),
        " — preliminary results for internal review; not for citation or external distribution."),
    div(class = "indev-watermark", "IN DEVELOPMENT")
  ),

  nav_panel(
    title = "Start here",
    icon  = bsicons::bs_icon("signpost-2"),
    mod_start_here_ui("start")
  ),

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

  # GBD comparison hidden for now (placeholder data). Restore by uncommenting
  # this nav_panel and the mod_gbd_compare_server() call in the server.
  # nav_panel(
  #   title = "GBD comparison",
  #   icon  = bsicons::bs_icon("globe2"),
  #   mod_gbd_compare_ui("gbd")
  # ),

  nav_panel(
    title = "Methods",
    icon  = bsicons::bs_icon("info-circle"),
    mod_methods_ui("methods")
  ),

  nav_panel(
    title = "Methods comparison",
    icon  = bsicons::bs_icon("clipboard2-check"),
    mod_methods_comparison_ui("methods_comp")
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
  mod_start_here_server("start")
  mod_map_explorer_server("map")
  mod_district_profile_server("district")
  mod_national_burden_server("burden")
  mod_decision_value_server("decision")
  mod_scenarios_server("scenarios")
  mod_importance_server("importance")
  mod_diagnostics_server("diagnostics")
  mod_benchmarks_server("benchmarks")
  mod_resolution_server("resolution")
  mod_transport_calibration_server("transport")
  mod_oos_civ_server("oos_civ")
  # mod_gbd_compare_server("gbd")   # GBD comparison hidden for now
  mod_methods_server("methods")
  mod_methods_comparison_server("methods_comp")
}

shinyApp(ui, server)
