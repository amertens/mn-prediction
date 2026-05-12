# =============================================================================
# Module: Map Explorer
# =============================================================================
# Interactive Admin-2 choropleth map with layer toggles, district drill-down,
# and a side panel showing detailed information for the clicked district.

mod_map_explorer_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 320,
      title = "Map controls",

      selectInput(ns("country"), "Country",
                  choices = country_choices,
                  selected = "ghana"),

      selectInput(ns("outcome"), "Outcome",
                  choices = outcome_choices,
                  selected = "women_iron"),

      radioButtons(ns("layer"), "Display layer",
                   choices = c("Predicted prevalence" = "pred_prev",
                               "Survey-observed prevalence (where available)" = "obs_prev",
                               "Confidence interval width" = "ci_width",
                               "Population at risk (count)" = "pop_at_risk",
                               "WHO public health classification" = "who_class"),
                   selected = "pred_prev"),

      hr(),

      htmlOutput(ns("headline")),

      hr(),

      htmlOutput(ns("district_detail"))
    ),

    layout_columns(
      col_widths = 12,
      card(
        full_screen = TRUE,
        card_header(
          "Subnational micronutrient deficiency map",
          downloadButton(ns("download_map_data"),
                          "Download CSV",
                          class = "btn-sm btn-outline-primary float-end")
        ),
        card_body(
          padding = 0,
          leafletOutput(ns("map"), height = "650px")
        ),
        card_footer(
          textOutput(ns("map_caption"))
        )
      )
    )
  )
}


mod_map_explorer_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive: joined sf for selected country × outcome
    map_data <- reactive({
      req(input$country, input$outcome)
      get_country_admin2(input$country, input$outcome,
                          admin2_bnds, admin2_pred, admin2_pop)
    })

    # National headline
    output$headline <- renderUI({
      df <- map_data()
      req(df)
      natl <- national_aggregate(df)
      ctry_label <- meta$countries[input$country]
      oc_label <- meta$outcome_labels[input$outcome]

      tagList(
        h5("National summary", style = "margin-top: 0;"),
        p(strong(ctry_label), " | ", em(oc_label)),
        p(
          "Estimated prevalence: ",
          tags$span(style = "font-size: 1.2em; color: #d7191c;",
                    strong(fmt_pct(natl$pred_prev_natl)))
        ),
        p(
          "Population affected: ",
          strong(fmt_count(natl$pop_at_risk_natl)),
          " of ", fmt_count(natl$pop_total), " total"
        ),
        if (!is.na(natl$ci_lo_natl)) {
          p("95% conformal CI: ",
            sprintf("[%s, %s]",
                    fmt_pct(natl$ci_lo_natl),
                    fmt_pct(natl$ci_hi_natl)),
            style = "font-size: 0.9em; color: #555;")
        }
      )
    })

    # Map output
    output$map <- renderLeaflet({
      df <- map_data()
      req(df, nrow(df) > 0)

      # Determine which column drives color and the palette
      layer <- input$layer
      vals <- df[[layer]]

      if (layer == "who_class") {
        pal <- colorFactor(palette = unname(who_colors),
                           levels = names(who_colors),
                           na.color = "#cccccc")
        fill_col <- pal(df$who_class)
        legend_title <- "WHO classification"
      } else if (layer == "pop_at_risk") {
        # Log-scale for population counts
        nz <- vals[!is.na(vals) & vals > 0]
        if (length(nz) > 0) {
          dom <- range(nz, na.rm = TRUE)
          pal <- colorNumeric("YlOrRd", domain = log10(pmax(dom, 1)),
                              na.color = "#cccccc")
          fill_col <- ifelse(is.na(vals) | vals == 0, "#cccccc",
                             pal(log10(pmax(vals, 1))))
        } else {
          fill_col <- rep("#cccccc", nrow(df))
          pal <- NULL
        }
        legend_title <- "Population at risk"
      } else {
        # Continuous prevalence-like layer
        nz <- vals[!is.na(vals)]
        if (length(nz) > 0) {
          dom <- range(nz, na.rm = TRUE)
          pal <- colorNumeric("YlOrRd", domain = dom, na.color = "#cccccc")
          fill_col <- pal(vals)
        } else {
          fill_col <- rep("#cccccc", nrow(df))
          pal <- NULL
        }
        legend_title <- switch(layer,
          "pred_prev" = "Predicted prevalence",
          "obs_prev"  = "Observed prevalence",
          "ci_width"  = "CI width (pp)",
          "Layer"
        )
      }

      # Build hover labels
      labels <- sprintf(
        "<strong>%s</strong><br/>%s<br/>Predicted: %s<br/>Population: %s<br/>WHO class: %s",
        df$Admin2,
        ifelse(is.na(df$Admin1), "", df$Admin1),
        fmt_pct(df$pred_prev),
        fmt_count(df$population),
        df$who_class
      ) |> lapply(HTML)

      m <- leaflet(df) |>
        addProviderTiles(providers$CartoDB.Positron) |>
        addPolygons(
          fillColor = fill_col,
          fillOpacity = 0.75,
          color = "white",
          weight = 1,
          highlightOptions = highlightOptions(
            weight = 3, color = "#333", bringToFront = TRUE
          ),
          label = labels,
          labelOptions = labelOptions(
            style = list("font-weight" = "normal", padding = "6px 8px"),
            textsize = "13px",
            direction = "auto"
          ),
          layerId = df$Admin2
        )

      if (!is.null(pal) && layer != "who_class") {
        legend_vals <- if (layer == "pop_at_risk") {
          legend_vals_actual <- vals[!is.na(vals) & vals > 0]
          log10(pmax(legend_vals_actual, 1))
        } else {
          vals[!is.na(vals)]
        }
        m <- m |> addLegend(
          pal = pal, values = legend_vals,
          opacity = 0.75, title = legend_title, position = "bottomright",
          labFormat = if (layer == "pop_at_risk") {
            labelFormat(transform = function(x) round(10^x))
          } else if (layer %in% c("pred_prev", "obs_prev")) {
            labelFormat(suffix = "", transform = function(x) round(x * 100, 1))
          } else {
            labelFormat()
          }
        )
      } else if (layer == "who_class") {
        m <- m |> addLegend(
          colors = unname(who_colors),
          labels = names(who_colors),
          opacity = 0.75, title = legend_title, position = "bottomright"
        )
      }

      m
    })

    # Reactive value for clicked district
    clicked <- reactiveVal(NULL)

    observeEvent(input$map_shape_click, {
      clicked(input$map_shape_click$id)
    })

    output$district_detail <- renderUI({
      district <- clicked()
      df <- map_data()
      if (is.null(district) || is.null(df)) {
        return(p(em("Click a district to see details."),
                 style = "color: #888;"))
      }

      row <- df[df$Admin2 == district, , drop = FALSE]
      if (nrow(row) == 0) return(NULL)

      conf <- confidence_badge(row$ci_width)

      tagList(
        h5(district, style = "margin-top: 0;"),
        if (!is.na(row$Admin1[1])) p(em(row$Admin1[1])),
        p(strong("Predicted prevalence: "),
          fmt_pct(row$pred_prev[1])),
        if (!is.na(row$ci_lo[1])) {
          p("95% CI: [",
            fmt_pct(row$ci_lo[1]), ", ", fmt_pct(row$ci_hi[1]), "]",
            br(),
            tags$span(class = "badge",
                      style = "background: #2c7bb6; color: white; padding: 2px 8px;",
                      conf))
        },
        if (!is.na(row$obs_prev[1])) {
          p(strong("Survey observed: "), fmt_pct(row$obs_prev[1]),
            br(),
            tags$small(sprintf("(n=%s individuals surveyed)",
                                fmt_count(row$n_survey[1]))))
        } else {
          p(em("No survey data for this district"),
            style = "color: #888;")
        },
        p(strong("WHO classification: "), row$who_class[1]),
        p(strong("Population: "), fmt_count(row$population[1])),
        p(strong("Population at risk: "), fmt_count(row$pop_at_risk[1]))
      )
    })

    output$map_caption <- renderText({
      ctry_label <- meta$countries[input$country]
      sy <- meta$survey_years[input$country]
      sprintf("Data: %s, survey year %s. Districts in gray have no model coverage (unsurveyed clusters).",
              ctry_label, sy)
    })

    output$download_map_data <- downloadHandler(
      filename = function() {
        sprintf("map_%s_%s_%s.csv",
                input$country, input$outcome, Sys.Date())
      },
      content = function(file) {
        df <- map_data()
        # Drop the geometry column for CSV export
        df_export <- as.data.frame(df)
        df_export$geometry <- NULL
        write.csv(df_export, file, row.names = FALSE)
      }
    )

  })
}
