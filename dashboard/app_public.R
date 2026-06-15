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
  theme = bs_theme(
    version = 5,
    bootswatch = "cosmo",
    primary = "#2c7bb6",
    base_font = font_google("Source Sans Pro")
  ),
  fillable = c("Map explorer"),
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
  nav_panel(title = "Importance", icon = bsicons::bs_icon("bar-chart-line"),
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
  mod_start_here_server("start")
  mod_map_explorer_server("map")
  mod_district_profile_server("district")
  mod_national_burden_server("burden")
  mod_decision_value_server("decision")
  mod_scenarios_server("scenarios")
  mod_importance_server("importance")
  mod_methods_server("methods")
}

shinyApp(ui, server)
