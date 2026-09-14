# =============================================================================
# Module: Cote d'Ivoire
# =============================================================================
# The case the project exists for: a country with the full proxy database and
# no biomarker survey, ranked from climate and soil layers by an index fitted
# on the four surveyed countries. It has never seen a Cote d'Ivoire biomarker.

mod_civ_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 340, title = "Cote d'Ivoire",
      selectInput(ns("outcome"), "Outcome", choices = NULL),
      radioButtons(ns("layer"), "What to show",
                   choices = c("Priority score" = "priority",
                               "How firmly placed (children's iron)" = "rank_width",
                               "Chance of the worst third (children's iron)" = "p_worst3rd"),
                   selected = "priority"),
      hr(),
      uiOutput(ns("summary")),
      hr(),
      p(style = "font-size:0.85em; color:#555;",
        "The index was fitted on The Gambia, Ghana, Sierra Leone and Malawi from climate and soil layers",
        " only, then applied here. Its accuracy in a country it has never seen is ",
        strong(fmt_num(Q$cs)), " across districts on the four surveyed countries, against ", fmt_num(Q$null_d),
        " for chance. There is no Cote d'Ivoire ground truth, so that is the only check available.")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Districts ranked from public data alone"),
           card_body(padding = 0, leafletOutput(ns("map"), height = "600px")),
           card_footer(tags$small(class = "text-muted",
                                  "Priority score 100 = ranked worst of 33 districts. The firmness layers come from",
                                  " refitting the index 400 times on resampled training districts."))),
      card(card_header("The list"),
           card_body(reactableOutput(ns("table")),
                     methods_note("Rank 1 is the district the model puts worst. For children's iron, the range is where",
                                  " the rank fell in 90 percent of 400 refits; a narrow range means the position does",
                                  " not depend much on which training districts the model learned from. It does not test",
                                  " whether the model transports to Cote d'Ivoire at all; the surveyed-country accuracy is",
                                  " the external bound on that."),
                     uiOutput(ns("guards"))))
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
      if (!is.null(CIV$uncertainty) && input$outcome == "child_iron") {
        u <- CIV$uncertainty; j <- match(.key(b$Admin1, b$Admin2), .key(u$Admin1, u$Admin2))
        b$rank_width <- u$rank_width[j]; b$p_worst3rd <- u$p_worst3rd[j]; b$rank_lo <- u$rank_lo[j]; b$rank_hi <- u$rank_hi[j]
      }
      b
    })

    output$summary <- renderUI({
      d <- sf::st_drop_geometry(dat()); req(nrow(d) > 0)
      worst <- d[order(d$rank_worst), ][1:min(6, nrow(d)), ]
      tagList(
        h5(outcome_short[[input$outcome]], style = "margin-top:0;"),
        p(strong("Ranked worst: "), paste(worst$Admin2, collapse = ", ")),
        if (input$outcome == "child_iron" && any(is.finite(d$p_worst3rd)))
          p(sprintf("%d of %d districts have at least an 80 percent chance of being in the worst third; %d have at most 20 percent. The median rank moves %s places across refits.",
                    sum(d$p_worst3rd >= 0.8, na.rm = TRUE), nrow(d), sum(d$p_worst3rd <= 0.2, na.rm = TRUE),
                    fmt_num(stats::median(d$rank_width, na.rm = TRUE), 0)), style = "font-size:0.9em;")
        else p("Firmness was computed for children's iron; choose it to see the range of ranks.", style = "font-size:0.85em; color:#777;")
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
      title <- switch(layer, priority = "Priority score", rank_width = "Places the rank can move", p_worst3rd = "Chance of worst third")
      labels <- sprintf("<strong>%s</strong><br/>%s<br/>Rank %s of %d<br/>%s", d$Admin2, d$Admin1, d$rank_worst, nrow(d),
                        ifelse(is.finite(d$rank_lo), sprintf("Rank range %s to %s; chance of worst third %s", fmt_num(d$rank_lo, 0), fmt_num(d$rank_hi, 0), fmt_pct(d$p_worst3rd, 0)), "")) |> lapply(HTML)
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
                      `Priority` = round(d$priority),
                      `Rank range (90%)` = ifelse(is.finite(d$rank_lo), sprintf("%s to %s", fmt_num(d$rank_lo, 0), fmt_num(d$rank_hi, 0)), "—"),
                      `Chance of worst third` = ifelse(is.finite(d$p_worst3rd), fmt_pct(d$p_worst3rd, 0), "—"), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 33, pagination = FALSE, height = 420)
    })

    output$guards <- renderUI({
      g <- CIV$guards; if (is.null(g)) return(NULL)
      tags$details(style = "font-size:0.82em; margin-top:8px;", tags$summary("Two checks run before this map was drawn"),
                   p("The same code with each surveyed country held out reproduces the published transport accuracy, and",
                     " adding Cote d'Ivoire to the column intersection costs a few climate and soil columns without degrading it."),
                   tags$pre(style = "font-size:0.8em; white-space:pre-wrap;", paste(capture.output(print(g, row.names = FALSE)), collapse = "\n")))
    })
  })
}
