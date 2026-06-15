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

      radioButtons(ns("admin_level"), "Geographic level",
                   choices = c("Admin 2 (district)" = "admin2",
                               "Admin 1 (region)"  = "admin1"),
                   selected = "admin2", inline = TRUE),

      radioButtons(ns("pred_model"), "Prediction model",
                   choices = c("Area-level SAE — HAL (all districts)" = "area",
                               "Fay-Herriot SAE (all districts, with intervals)" = "fh",
                               "Individual SuperLearner (surveyed districts)" = "sl"),
                   selected = "area"),
      tags$details(
        style = "font-size:0.82em; color:#555; margin-top:-2px; margin-bottom:6px;",
        tags$summary("How are these predictions built?"),
        tags$ul(
          style = "padding-left:1.1em; margin-bottom:0;",
          tags$li(strong("Area-level SAE (HAL): "),
                  "predicts every district's prevalence from satellite and ",
                  "geospatial indicators using a Highly Adaptive Lasso. Full ",
                  "coverage; point estimates only (no per-district interval)."),
          tags$li(strong("Fay-Herriot SAE: "),
                  "a small-area model that blends each district's survey ",
                  "estimate (where one exists) with the indicator-based ",
                  "prediction. Full coverage with a 95% interval for every ",
                  "district — narrow where a survey exists, wide where it does not."),
          tags$li(strong("Individual SuperLearner: "),
                  "a person-level machine-learning ensemble, averaged up to ",
                  "the surveyed districts. Covers surveyed districts only, ",
                  "with conformal intervals.")
        )
      ),

      radioButtons(ns("layer"), "Display layer",
                   choices = c("Predicted prevalence" = "pred_prev",
                               "Survey-observed prevalence (where available)" = "obs_prev",
                               "Modeled − survey difference (pp)" = "diff_prev",
                               "Transportability error: held-out country (pp)" = "loco_diff",
                               "Confidence interval width" = "ci_width",
                               "Population affected (count)" = "pop_at_risk",
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

    # Which prediction source feeds the map: individual SuperLearner (surveyed
    # districts only) or the full-coverage area-level SAE model. Falls back to
    # SuperLearner if the area bundle is not present.
    pred_df <- reactive({
      m <- input$pred_model %||% "area"
      if (m == "fh" && exists("admin2_fh_pred") && !is.null(admin2_fh_pred)) {
        admin2_fh_pred
      } else if (m == "area" && exists("admin2_area_pred") && !is.null(admin2_area_pred)) {
        admin2_area_pred
      } else {
        admin2_pred
      }
    })

    # Reactive: joined sf for selected country × outcome × admin level
    map_data <- reactive({
      req(input$country, input$outcome, input$admin_level)
      pd <- pred_df()
      if (input$admin_level == "admin1") {
        get_country_admin1(input$country, input$outcome,
                            admin1_bnds, admin2_bnds, pd, admin2_pop)
      } else {
        get_country_admin2(input$country, input$outcome,
                            admin2_bnds, pd, admin2_pop)
      }
    })

    # Convenience: the area-name column for the active level
    area_col <- reactive(
      if (isTRUE(input$admin_level == "admin1")) "Admin1" else "Admin2"
    )

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
          " of ", fmt_count(natl$pop_total), " ",
          if (startsWith(input$outcome, "child_")) "children 6–59 months"
          else "women 15–49 years",
          tags$br(),
          tags$span(
            style = "font-size: 0.82em; color: #777;",
            switch(input$pred_model %||% "area",
              sl = "Counts cover surveyed districts only — switch to an all-district model above.",
              fh = "Counts cover all districts (Fay-Herriot model, with intervals).",
              "Counts cover all districts (area-level HAL model; point estimates, no interval)."))
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
      } else if (layer %in% c("diff_prev", "loco_diff")) {
        # Diverging palette centered at 0: blue = model UNDER-predicts vs survey,
        # red = model OVER-predicts. Domain is symmetric around 0.
        nz <- vals[is.finite(vals)]
        if (length(nz) > 0) {
          M <- max(abs(nz), na.rm = TRUE)
          if (!is.finite(M) || M == 0) M <- 0.01
          pal <- colorNumeric(c("#2166ac", "#f7f7f7", "#b2182b"),
                              domain = c(-M, M), na.color = "#cccccc")
          fill_col <- ifelse(is.finite(vals), pal(vals), "#cccccc")
        } else {
          fill_col <- rep("#cccccc", nrow(df)); pal <- NULL
        }
        legend_title <- if (layer == "diff_prev")
          "Modeled − survey (pp)" else "Transport − survey (pp)"
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

      # Build hover labels. At Admin-1 the subtitle line is hidden (no parent).
      area_name <- df[[area_col()]]
      sub_line <- if (area_col() == "Admin2") {
        ifelse(is.na(df$Admin1), "", df$Admin1)
      } else {
        rep("", nrow(df))
      }
      fmt_pp <- function(x) ifelse(is.na(x), "—", sprintf("%+.1f pp", x * 100))
      # Difference line shown for the difference layers (and harmless otherwise)
      diff_line <- if (layer == "loco_diff") {
        sprintf("Transport modeled: %s<br/>Survey (harmonized): %s<br/><strong>Transport − survey: %s</strong><br/>",
                fmt_pct(df$loco_pred_prev), fmt_pct(df$loco_survey_prev), fmt_pp(df$loco_diff))
      } else if (layer == "diff_prev") {
        sprintf("Modeled: %s<br/>Survey: %s<br/><strong>Modeled − survey: %s</strong><br/>",
                fmt_pct(df$pred_prev), fmt_pct(df$obs_prev), fmt_pp(df$diff_prev))
      } else {
        sprintf("Predicted: %s<br/>", fmt_pct(df$pred_prev))
      }
      labels <- sprintf(
        "<strong>%s</strong><br/>%s<br/>%sPopulation: %s<br/>WHO class: %s",
        area_name, sub_line,
        diff_line,
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
          layerId = df[[area_col()]]
        )

      if (!is.null(pal) && layer != "who_class") {
        legend_vals <- if (layer == "pop_at_risk") {
          legend_vals_actual <- vals[!is.na(vals) & vals > 0]
          log10(pmax(legend_vals_actual, 1))
        } else if (layer %in% c("diff_prev", "loco_diff")) {
          M <- max(abs(vals[is.finite(vals)]), na.rm = TRUE)
          c(-M, M)  # symmetric, centered at 0
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
          } else if (layer %in% c("diff_prev", "loco_diff")) {
            labelFormat(suffix = " pp", transform = function(x) round(x * 100, 1))
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
      area <- clicked()
      df <- map_data()
      level_label <- if (isTRUE(input$admin_level == "admin1"))
                       "region" else "district"
      if (is.null(area) || is.null(df)) {
        return(p(em(sprintf("Click a %s to see details.", level_label)),
                 style = "color: #888;"))
      }

      row <- df[df[[area_col()]] == area, , drop = FALSE]
      if (nrow(row) == 0) return(NULL)

      conf <- confidence_badge(row$ci_width)

      tagList(
        h5(area, style = "margin-top: 0;"),
        if (area_col() == "Admin2" && !is.na(row$Admin1[1]))
          p(em(row$Admin1[1])),
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
          tagList(
            p(strong("Survey observed: "), fmt_pct(row$obs_prev[1]),
              br(),
              tags$small(sprintf("(n=%s individuals surveyed)",
                                  fmt_count(row$n_survey[1])))),
            if (!is.na(row$diff_prev[1]))
              p(strong("Modeled − survey: "),
                sprintf("%+.1f pp", row$diff_prev[1] * 100)),
            if (!is.null(row$loco_diff) && !is.na(row$loco_diff[1]))
              p(strong("Transportability error: "),
                sprintf("%+.1f pp", row$loco_diff[1] * 100),
                br(),
                tags$small("(model trained on other countries only)"))
          )
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
      lvl <- if (isTRUE(input$admin_level == "admin1"))
               "Admin-1 (region)" else "Admin-2 (district)"
      layer_note <- switch(input$layer,
        diff_prev = " Difference layer: modeled (within-country) minus survey-observed prevalence; red = model over-predicts, blue = under-predicts.",
        loco_diff = " Transportability layer: prevalence predicted by the universal model trained ONLY on the other countries, minus the local survey, using a harmonized cross-country outcome (uniform WHO cutoffs: ferritin<12/15 = iron deficiency, RBP<0.70 = VAD). Red = over-predicts, blue = under-predicts. Gray = country/outcome not covered.",
        ""
      )
      sprintf(
        "Data: %s, survey year %s. %s view. %s in gray have no model coverage. Admin-1 values are population-weighted aggregates of Admin-2 predictions.%s",
        ctry_label, sy, lvl,
        if (isTRUE(input$admin_level == "admin1")) "Regions" else "Districts",
        layer_note
      )
    })

    # Clear the click selection when the user switches admin level
    observeEvent(input$admin_level, { clicked(NULL) }, ignoreInit = TRUE)

    output$download_map_data <- downloadHandler(
      filename = function() {
        sprintf("map_%s_%s_%s_%s.csv",
                input$country, input$outcome, input$admin_level, Sys.Date())
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
