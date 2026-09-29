# =============================================================================
# Module: Cote d'Ivoire
# =============================================================================
# A country with the full public-data layers and no national biomarker survey,
# ranked from climate and soil layers by a model built on the four surveyed
# countries. The model has never seen a Cote d'Ivoire biomarker.

mod_civ_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 340, title = "Cote d'Ivoire",
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
        "The model was built from climate and soil data for The Gambia, Ghana, Sierra Leone and Malawi, then applied",
        " here. Tested on those four countries with each one left out in turn, it ranked districts at ",
        strong(fmt_num(Q$cs)), " (a random ranking stays below ", fmt_num(Q$null_d), "). The climate-and-soil",
        " combination was chosen after seeing those results, so this figure may flatter it. Cote d'Ivoire has no",
        " national biomarker survey to check the map against.")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Districts ranked from public data alone"),
           card_body(padding = 0, leafletOutput(ns("map"), height = "600px")),
           card_footer(tags$small(class = "text-muted",
                                  "Priority score 100 = ranked worst of 33 districts. The rank ranges come from re-estimating",
                                  " the model 400 times on different samples of the training districts."))),
      card(card_header("The ranked list"),
           card_body(reactableOutput(ns("table")),
                     methods_note("Rank 1 is the district the model places worst. The range shows where the rank fell in",
                                  " 90% of 400 re-estimates; a narrow range means the position does not depend much on which",
                                  " training districts the model learned from. It does not show whether the model is right for",
                                  " Cote d'Ivoire. For that, see the tests with a country left out and the six-country check",
                                  " under Can we trust it?"),
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
      # rank ranges for every outcome (climate-and-soil version); the older
      # single-outcome file remains the fallback for older bundles
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
      worst <- d[order(d$rank_worst), ][1:min(6, nrow(d)), ]
      tagList(
        h5(outcome_short[[input$outcome]], style = "margin-top:0;"),
        p(strong("Ranked worst: "), paste(worst$Admin2, collapse = ", ")),
        if (any(is.finite(d$p_worst3rd)))
          p(sprintf("%d of %d districts were in the worst third in at least 80%% of re-estimates, and %d in at most 20%%. The median rank moves %s places.",
                    sum(d$p_worst3rd >= 0.8, na.rm = TRUE), nrow(d), sum(d$p_worst3rd <= 0.2, na.rm = TRUE),
                    fmt_num(stats::median(d$rank_width, na.rm = TRUE), 0)), style = "font-size:0.9em;")
        else p("Rank ranges have not been computed for this outcome.", style = "font-size:0.85em; color:#777;"),
        tags$details(style = "font-size:0.84em; margin-top:8px;",
          tags$summary(strong("Other evidence that agrees with this map")),
          tags$ul(style = "padding-left:1.1em; margin-top:4px;",
            tags$li("A 2007 survey measured women's vitamin B12 in nine zones of Cote d'Ivoire. Its ranking of the zones",
                    " agrees with this map's at 0.95. That is one outcome and nine zones, so it is encouraging rather than conclusive."),
            tags$li("The one published external map (WFP and MIMI vitamin A intake, Nature Food 2026) also places the",
                    " northern savanna districts worst."),
            tags$li("Two versions of the model built from different data groups agree with each other at 0.88 to 0.97.",
                    " This shows the map does not hinge on one choice of data; it does not show that the map is correct.")))
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
                        ifelse(is.finite(d$rank_lo), sprintf("Rank range %s to %s; in the worst third in %s of runs", fmt_num(d$rank_lo, 0), fmt_num(d$rank_hi, 0), fmt_pct(d$p_worst3rd, 0)), "")) |> lapply(HTML)
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
                      `In worst third (share of runs)` = ifelse(is.finite(d$p_worst3rd), fmt_pct(d$p_worst3rd, 0), "—"), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 33, pagination = FALSE, height = 420)
    })

    output$guards <- renderUI({
      g <- CIV$guards; if (is.null(g)) return(NULL)
      tags$details(style = "font-size:0.82em; margin-top:8px;", tags$summary("Technical checks run before this map was drawn"),
                   p("Run with each surveyed country left out, the same code reproduces the published accuracy for new",
                     " countries. Adding Cote d'Ivoire meant dropping a few climate and soil layers that it lacks, and this",
                     " did not reduce accuracy."),
                   tags$pre(style = "font-size:0.8em; white-space:pre-wrap;", paste(capture.output(print(g, row.names = FALSE)), collapse = "\n")))
    })
  })
}
