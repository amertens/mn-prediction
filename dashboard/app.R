# =============================================================================
# Micronutrient Burden Dashboard
# =============================================================================
# District rankings of micronutrient deficiency from public data, scored
# against four national biomarker surveys (protocol v2). Rebuilt 2026-09-13;
# policymaker revamp 2026-09-27: headline-tier fit, external validation
# (XV-01/02), stability ensembles and WHO exceedance on every district (UE-01),
# survey planner with SP-01 validation, searchable importance with ranges,
# What-more-data-buys, Malawi selenium/iodine, country briefs.
#
# To run locally:
#   Rscript dashboard/data-raw/05_build_protocol_v2_bundles.R   # data, once
#   Rscript dashboard/data-raw/06_build_uncertainty_ensembles.R # stability, once
#   Rscript dashboard/data-raw/07_build_country_briefs.R        # briefs, once
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
    nav_panel(title = "The data behind it", icon = bsicons::bs_icon("journal-text"),
              navset_card_tab(
                nav_panel(title = "What tracks which nutrient", icon = bsicons::bs_icon("clipboard2-pulse"), mod_nutrient_signal_ui("nutrient")),
                nav_panel(title = "Predictor catalogue", icon = bsicons::bs_icon("journal-text"), mod_catalogue_ui("catalogue"))))
  ),

  nav_menu(
    title = "Can we trust it?", icon = bsicons::bs_icon("shield-check"),
    nav_panel(title = "How well it works", icon = bsicons::bs_icon("speedometer2"), mod_trust_ui("trust")),
    nav_panel(title = "Tested in six more countries", icon = bsicons::bs_icon("globe-americas"), mod_external_ui("external")),
    nav_panel(title = "What the ranking buys", icon = bsicons::bs_icon("bullseye"), mod_targeting_ui("targeting")),
    nav_panel(title = "What more data buys", icon = bsicons::bs_icon("graph-up-arrow"), mod_roadmap_ui("roadmap")),
    nav_panel(title = "Methods", icon = bsicons::bs_icon("info-circle"), mod_methods_ui("methods"))
  ),

  nav_panel(title = "Plan a survey", icon = bsicons::bs_icon("cash-coin"), mod_survey_design_ui("planning")),

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

  # ── URL state ──────────────────────────────────────────────────────────────
  # The address bar tracks the open tab and its main selectors, so any view can
  # be shared by copying the link. A pasted link is applied twice, one flush
  # apart, because choosing a country repopulates its outcome choices.
  q0 <- isolate(parseQueryString(session$clientData$url_search))
  apply_query <- function(q) {
    if (!is.null(q$tab)) try(nav_select("main_nav", q$tab), silent = TRUE)
    pairs <- c(country = "map-country", outcome = "map-outcome", layer = "map-layer", level = "map-admin_level",
               dcountry = "district-country", district = "district-district", doutcome = "district-outcome",
               civ = "civ-outcome",
               pcountry = "planning-country", poutcome = "planning-outcome", preset = "planning-preset")
    for (nm in names(pairs)) if (!is.null(q[[nm]])) updateSelectInput(session, pairs[[nm]], selected = q[[nm]])
    if (!is.null(q$k)) updateSliderInput(session, "planning-k", value = suppressWarnings(as.integer(q$k)))
  }
  url_done <- reactiveVal(!length(q0))
  url_tries <- 0L
  observe({
    if (isTRUE(url_done())) return()
    apply_query(q0)
    url_tries <<- url_tries + 1L
    # a country update repopulates its outcome choices on the NEXT client
    # round-trip, which resets the outcome we just set; re-apply on a timer
    # until the wiring has settled, then hand the address bar to the writer
    if (url_tries < 3L) invalidateLater(1400) else url_done(TRUE)
  })
  observe({
    req(isTRUE(url_done()))
    tab <- input$main_nav
    enc <- function(nm, id, default = NULL) {
      v <- input[[id]]
      if (is.null(v) || length(v) != 1 || !nzchar(as.character(v))) return(NULL)
      if (!is.null(default) && identical(as.character(v), default)) return(NULL)
      sprintf("%s=%s", nm, utils::URLencode(as.character(v), reserved = TRUE))
    }
    parts <- if (!is.null(tab) && tab != "Start here") c(sprintf("tab=%s", utils::URLencode(tab, reserved = TRUE))) else character(0)
    parts <- c(parts, switch(tab %||% "",
      "Map explorer" = c(enc("country", "map-country"), enc("outcome", "map-outcome"),
                         enc("layer", "map-layer", "priority"), enc("level", "map-admin_level", "admin2")),
      "District profiles" = c(enc("dcountry", "district-country"), enc("district", "district-district"), enc("doutcome", "district-outcome")),
      "Cote d'Ivoire" = enc("civ", "civ-outcome"),
      "Plan a survey" = c(enc("pcountry", "planning-country"), enc("poutcome", "planning-outcome"),
                          enc("preset", "planning-preset", "shrink"), enc("k", "planning-k")),
      NULL))
    updateQueryString(if (length(parts)) paste0("?", paste(parts, collapse = "&")) else "?", mode = "replace", session = session)
  })
  mod_start_here_server("start", go_to = go_to)
  mod_map_explorer_server("map")
  mod_district_server("district")
  mod_civ_server("civ")
  mod_importance_server("importance")
  mod_nutrient_signal_server("nutrient")
  mod_catalogue_server("catalogue")
  mod_trust_server("trust")
  mod_external_server("external")
  mod_targeting_server("targeting")
  mod_roadmap_server("roadmap", go_to = go_to)
  mod_methods_server("methods")
  mod_survey_design_server("planning")
}

shinyApp(ui, server)
