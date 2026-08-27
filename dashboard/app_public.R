# =============================================================================
# app_public.R — lean, policymaker-facing build of the Micronutrient dashboard
# =============================================================================
# Shares all modules and data with the full app (app.R) via global.R, but shows
# only the policymaker-facing tabs (no analyst/technical or low-confidence tabs).
# Fewer simultaneous heavy renders -> more stable for external users, and this is
# the shareable artifact. Deploy with: Rscript dashboard/deploy_public.R
# =============================================================================

source("global.R")

ui <- page_navbar(
  title = "Micronutrient Burden",
  id = "main_nav",
  theme = bs_theme(
    version = 5,
    bootswatch = "cosmo",
    primary = "#2c7bb6",
    base_font = font_google("Source Sans Pro")
  ),
  fillable = c("Map explorer"),
  underline = TRUE,

  # Same banner as app.R — defined in global.R so the two cannot drift.
  header = site_banner,

  nav_panel(title = "Start here", icon = bsicons::bs_icon("signpost-2"),
            mod_start_here_ui("start")),
  nav_panel(title = "Map explorer", icon = bsicons::bs_icon("map"),
            mod_map_explorer_ui("map")),
  nav_panel(title = "District profiles", icon = bsicons::bs_icon("file-earmark-medical"),
            mod_district_profile_ui("district")),
  nav_panel(title = "National burden", icon = bsicons::bs_icon("graph-up"),
            mod_national_burden_ui("burden")),
  nav_panel(title = "Decision value", icon = bsicons::bs_icon("bullseye"),
            mod_decision_value_ui("decision")),
  nav_panel(title = "Scenarios", icon = bsicons::bs_icon("sliders"),
            mod_scenarios_ui("scenarios")),
  nav_panel(title = "What drives the estimate",
            icon = bsicons::bs_icon("bar-chart-line"),
            mod_importance_ui("importance")),
  nav_panel(title = "Methods", icon = bsicons::bs_icon("info-circle"),
            mod_methods_ui("methods")),

  nav_spacer(),
  nav_item(popover(tags$button(class = "btn btn-link",
                               bsicons::bs_icon("book"), " Glossary"),
                   glossary_content, title = NULL, placement = "bottom",
                   options = list(html = TRUE, container = "body"))),
  nav_item(popover(tags$button(class = "btn btn-link",
                               bsicons::bs_icon("info-circle"), " About"),
                   about_content, title = NULL, placement = "bottom",
                   options = list(html = TRUE, container = "body"))),

  footer = div(
    style = "text-align: center; color: #888; font-size: 0.8em; padding: 8px;",
    sprintf("Micronutrient Burden Dashboard | Data built %s", data_build_time)
  )
)

server <- function(input, output, session) {
  go_to <- function(tab, country = NULL, outcome = NULL) {
    if (!is.null(country)) {
      updateSelectInput(session, "map-country",         selected = country)
      updateSelectInput(session, "district-country",    selected = country)
      updateSelectInput(session, "decision-ta_country", selected = country)
    }
    if (!is.null(outcome)) {
      updateSelectInput(session, "map-outcome",         selected = outcome)
      updateSelectInput(session, "decision-ta_outcome", selected = outcome)
    }
    nav_select("main_nav", tab)
  }

  mod_start_here_server("start", go_to = go_to)
  mod_map_explorer_server("map")
  mod_district_profile_server("district")
  mod_national_burden_server("burden")
  mod_decision_value_server("decision")
  mod_scenarios_server("scenarios")
  mod_importance_server("importance")
  mod_methods_server("methods")
}

shinyApp(ui, server)
