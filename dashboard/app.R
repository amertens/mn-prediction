# =============================================================================
# Micronutrient Burden Dashboard
# =============================================================================
# District rankings of micronutrient deficiency from public data, scored
# against four national biomarker surveys (protocol v2). Rebuilt 2026-09-13.
#
# To run locally:
#   Rscript dashboard/data-raw/05_build_protocol_v2_bundles.R   # data, once
#   setwd("dashboard"); shiny::runApp()
# To deploy:
#   Rscript dashboard/deploy.R
# =============================================================================

source("global.R")

ui <- page_navbar(
  title = "Micronutrient Burden",
  id = "main_nav",
  theme = bs_theme(version = 5, bootswatch = "cosmo", primary = "#0F7B8A",
                   base_font = font_google("Source Sans Pro")),
  fillable = c("Map explorer", "Cote d'Ivoire"),
  underline = TRUE,
  header = site_banner,

  # Navigation follows the three questions the work answers, in the order a
  # reader asks them, with the two tabs a programme person opens first kept at
  # the top level.
  nav_panel(title = "Start here", icon = bsicons::bs_icon("signpost-2"), mod_start_here_ui("start")),

  nav_menu(
    title = "Where is deficiency?", icon = bsicons::bs_icon("map"),
    nav_panel(title = "Map explorer", icon = bsicons::bs_icon("map"), mod_map_explorer_ui("map")),
    nav_panel(title = "District profiles", icon = bsicons::bs_icon("file-earmark-medical"), mod_district_ui("district")),
    nav_panel(title = "Cote d'Ivoire", icon = bsicons::bs_icon("compass"), mod_civ_ui("civ"))
  ),

  nav_menu(
    title = "What drives it?", icon = bsicons::bs_icon("bar-chart-line"),
    nav_panel(title = "What drives the estimate", icon = bsicons::bs_icon("bar-chart-line"), mod_importance_ui("importance")),
    nav_panel(title = "What tracks which nutrient", icon = bsicons::bs_icon("clipboard2-pulse"), value = "nutrient",
              mod_nutrient_signal_ui("nutrient")),
    nav_panel(title = "Predictor catalogue", icon = bsicons::bs_icon("journal-text"), value = "catalogue",
              mod_catalogue_ui("catalogue"))
  ),

  nav_menu(
    title = "Can we trust it?", icon = bsicons::bs_icon("shield-check"),
    nav_panel(title = "How well it works", icon = bsicons::bs_icon("speedometer2"), value = "trust", mod_trust_ui("trust")),
    nav_panel(title = "What the ranking buys", icon = bsicons::bs_icon("bullseye"), value = "targeting", mod_targeting_ui("targeting")),
    nav_panel(title = "Methods", icon = bsicons::bs_icon("info-circle"), mod_methods_ui("methods"))
  ),

  nav_panel(title = "Plan a survey", icon = bsicons::bs_icon("cash-coin"), value = "planning", mod_survey_design_ui("planning")),

  nav_spacer(),
  nav_item(popover(tags$button(class = "btn btn-link", bsicons::bs_icon("book"), " Glossary"),
                   glossary_content, title = NULL, placement = "bottom", options = list(html = TRUE, container = "body"))),
  nav_item(popover(tags$button(class = "btn btn-link", bsicons::bs_icon("info-circle"), " About"),
                   about_content, title = NULL, placement = "bottom", options = list(html = TRUE, container = "body"))),

  footer = div(style = "text-align: center; color: #888; font-size: 0.8em; padding: 8px;",
               sprintf("Micronutrient Burden Dashboard | %s | data built %s", PROTOCOL_LABEL, data_build_time))
)

server <- function(input, output, session) {
  go_to <- function(tab, country = NULL, outcome = NULL) {
    if (!is.null(country)) {
      updateSelectInput(session, "map-country", selected = country)
      updateSelectInput(session, "district-country", selected = country)
    }
    if (!is.null(outcome)) {
      updateSelectInput(session, "map-outcome", selected = outcome)
      updateSelectInput(session, "targeting-outcome", selected = outcome)
    }
    nav_select("main_nav", tab)
  }
  mod_start_here_server("start", go_to = go_to)
  mod_map_explorer_server("map")
  mod_district_server("district")
  mod_civ_server("civ")
  mod_importance_server("importance")
  mod_nutrient_signal_server("nutrient")
  mod_catalogue_server("catalogue")
  mod_trust_server("trust")
  mod_targeting_server("targeting")
  mod_methods_server("methods")
  mod_survey_design_server("planning")
}

shinyApp(ui, server)
