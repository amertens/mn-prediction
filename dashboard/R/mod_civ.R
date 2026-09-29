# =============================================================================
# Module: Cote d'Ivoire (concise version)
# =============================================================================

mod_civ_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330, title = "Cote d'Ivoire",
      selectInput(ns("outcome"), "Outcome", choices = NULL),
      radioButtons(ns("layer"), "What to show",
                   choices = c("Priority score" = "priority",
                               "Rank range when re-estimated" = "rank_width",
                               "How often in the worst third" = "p_worst3rd"),
                   selected = "priority"),
      hr(),
      uiOutput(ns("summary")),
      hr(),
      p(style = "font-size:0.85em; color:#555;",
        "Ranked from climate and soil data by a model built on the four surveyed countries. With a surveyed country",
        " left out, that model scored ", strong(fmt_num(Q$cs)), ". No national biomarker survey exists here to check it.")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Districts ranked from public data alone"),
           card_body(padding = 0, leafletOutput(ns("map"), height = "600px")),
           card_footer(tags$small(class = "text-muted", "100 = ranked worst of 33 districts."))),
      card(card_header("The ranked list"),
           card_body(reactableOutput(ns("table")),
                     p(style = "font-size:0.85em; color:#666;", "The range shows where the rank fell in 90% of 400 re-estimates.")))
    )
  )
}

mod_civ_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    observe({
      req(CIV)
      ocs <- intersect(names(outcome_short), unique(CIV$ranking$outcome))
      updateSelectInput(session, "outcome", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
    })

    dat <- reactive({
      req(CIV, input$outcome)
      r <- CIV$ranking[CIV$ranking$outcome == input$outcome, ]
      b <- CIV$boundaries
      i <- match(.key(b$Admin1, b$Admin2), .key(r$Admin1, r$Admin2))
      b$priority <- r$priority[i]; b$rank_worst <- r$rank_worst[i]; b$index <- r$index[i]
      b$rank_width <- NA_real_; b$p_worst3rd <- NA_real_; b$rank_lo <- NA_real_; b$rank_hi <- NA_real_
      u <- NULL
      if (!is.null(CIV$uncertainty_all)) {
        u <- CIV$uncertainty_all[CIV$uncertainty_all$outcome == input$outcome & CIV$uncertainty_all$domain_set == "cs", ]
        if (!nrow(u)) u <- NULL
      }
      if (is.null(u) && !is.null(CIV$uncertainty) && input$outcome == "child_iron") u <- CIV$uncertainty
      if (!is.null(u)) {
        j <- match(.key(b$Admin1, b$Admin2), .key(u$Admin1, u$Admin2))
        b$rank_width <- u$rank_width[j]; b$p_worst3rd <- u$p_worst3rd[j]; b$rank_lo <- u$rank_lo[j]; b$rank_hi <- u$rank_hi[j]
      }
      b
    })

    output$summary <- renderUI({
      d <- sf::st_drop_geometry(dat()); req(nrow(d) > 0)
      worst <- d[order(d$rank_worst), ][1:min(5, nrow(d)), ]
      tagList(
        h5(outcome_short[[input$outcome]], style = "margin-top:0;"),
        p(strong("Ranked worst: "), paste(worst$Admin2, collapse = ", ")),
        tags$ul(style = "font-size:0.84em; padding-left:1.1em;",
          tags$li("A 2007 survey of women's B12 in nine zones agrees with this ranking (0.95)."),
          tags$li("The WFP and MIMI vitamin A map also puts the northern districts worst."))
      )
    })

    output$map <- renderLeaflet({
      d <- dat(); req(nrow(d) > 0)
      layer <- input$layer
      if (layer != "priority" && !any(is.finite(d[[layer]]))) layer <- "priority"
      vals <- d[[layer]]
      pal <- switch(layer,
                    priority = colorNumeric("YlOrRd", domain = c(0, 100), na.color = "#d9d9d9"),
                    rank_width = colorNumeric(c("#252525", "#bdbdbd", "#f7f7f7"), domain = range(vals, na.rm = TRUE), na.color = "#d9d9d9"),
                    p_worst3rd = colorNumeric(c("#f7fbff", "#9ecae1", "#08519c"), domain = c(0, 1), na.color = "#d9d9d9"))
      title <- switch(layer, priority = "Priority score", rank_width = "Places the rank moves", p_worst3rd = "In worst third (share of runs)")
      labels <- sprintf("<strong>%s</strong><br/>%s<br/>Rank %s of %d<br/>%s", d$Admin2, d$Admin1, d$rank_worst, nrow(d),
                        ifelse(is.finite(d$rank_lo), sprintf("Rank range %s to %s", fmt_num(d$rank_lo, 0), fmt_num(d$rank_hi, 0)), "")) |> lapply(HTML)
      leaflet(d) |> addProviderTiles(providers$Esri.WorldGrayCanvas) |>
        addPolygons(fillColor = pal(vals), fillOpacity = 0.78, color = "#666", weight = 0.7,
                    highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE),
                    label = labels, labelOptions = labelOptions(textsize = "13px")) |>
        addLegend(pal = pal, values = vals[is.finite(vals)], title = title, position = "bottomright", opacity = 0.78,
                  labFormat = if (layer == "p_worst3rd") labelFormat(transform = function(x) round(100 * x), suffix = "%") else labelFormat())
    })

    output$table <- renderReactable({
      d <- sf::st_drop_geometry(dat()); req(nrow(d) > 0)
      d <- d[order(d$rank_worst), ]
      t <- data.frame(Rank = d$rank_worst, District = d$Admin2, Region = d$Admin1,
                      `Rank range (90%)` = ifelse(is.finite(d$rank_lo), sprintf("%s to %s", fmt_num(d$rank_lo, 0), fmt_num(d$rank_hi, 0)), "—"),
                      check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 33, pagination = FALSE, height = 420)
    })
  })
}
