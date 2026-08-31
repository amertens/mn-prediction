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

      uiOutput(ns("outcome_caveat")),

      radioButtons(ns("admin_level"), "Level",
                   choices = c("District" = "admin2", "Region" = "admin1"),
                   selected = "admin2", inline = TRUE),

      # Grouped so the four everyday views sit apart from the two that only
      # mean something to an analyst.
      selectInput(ns("layer"), "What to show",
                  choices = list(
                    "Deficiency" = c(
                      "Estimated prevalence"          = "pred_prev",
                      "Survey-measured prevalence"    = "obs_prev",
                      "WHO severity class"            = "who_class"),
                    "People" = c(
                      "Number of people affected"     = "pop_at_risk"),
                    "How firm is the number" = c(
                      "Width of the 95% range"        = "ci_width",
                      "Could the severity class change?" = "misclass"),
                    "Model checks" = c(
                      "Estimate minus survey (pp)"    = "diff_prev",
                      "Error when borrowed from other countries (pp)" = "loco_diff")),
                  selected = "pred_prev"),

      # The model choice is a statistical judgement, not a program one, so it
      # starts closed on a sensible default rather than asking every reader to
      # make it. DEFAULT_PRED_MODEL lives in global.R and is shared with the
      # district profiles, so both tabs open on the same estimator.
      tags$details(
        style = "margin:2px 0 8px;",
        tags$summary(style = "cursor:pointer; font-size:0.86em; color:#2c7bb6;",
                     "Advanced: choose a different model"),
        div(
          style = "padding:6px 0 0 2px;",
          radioButtons(ns("pred_model"), NULL,
                       choices  = pred_model_choices(),
                       selected = DEFAULT_PRED_MODEL),
          pred_model_help()
        )
      ),

      hr(),

      uiOutput(ns("no_signal")),
      uiOutput(ns("layer_unavailable")),
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
          textOutput(ns("map_caption")),
          # The shading a reader takes off this map is a level, and the level
          # comes from the survey, not the covariates. Say so where the map is,
          # not only on Start Here -- and point at the evidence.
          tags$small(
            class = "text-muted",
            "Shading is anchored to each region's design-based survey total, ",
            "so the level is the survey's and the within-region pattern is the ",
            "model's. See ", tags$strong("Resolution & anchoring → Where ",
            "the level comes from"), " for what that anchoring is worth."
          )
        )
      )
    )
  )
}


mod_map_explorer_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Which table feeds the map. Resolved in global.R so every tab that shows a
    # prevalence resolves it the same way.
    pred_df <- reactive(pred_model_data(input$pred_model %||% DEFAULT_PRED_MODEL))

    # The two "how firm is the number" views need a 95% range, and the default
    # model does not carry one yet. Say so plainly and name the model that does,
    # rather than drawing an empty map.
    output$layer_unavailable <- renderUI({
      needs_ci <- (input$layer %||% "") %in% c("ci_width", "misclass")
      if (!needs_ci || pred_model_has_ci(input$pred_model)) return(NULL)
      div(class = "alert alert-warning",
          style = "font-size:0.82em; padding:8px 10px; margin:6px 0 0;",
          bsicons::bs_icon("exclamation-triangle"), " ",
          "This view needs a 95% range, which the recommended model does not ",
          "produce yet. Open ", strong("Advanced"), " above and pick a model ",
          "with uncertainty ranges to see it.")
    })

    # Some country x outcome combinations come back the same in every district:
    # the model found nothing to separate them. The map then looks like a map
    # but is a single flat colour, which reads as "no variation in deficiency"
    # rather than "no information". Say which it is.
    output$no_signal <- renderUI({
      df <- tryCatch(map_data(), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
      if (has_district_signal(df$pred_prev)) return(NULL)
      v <- df$pred_prev[is.finite(df$pred_prev)]
      div(class = "alert alert-secondary",
          style = "font-size:0.82em; padding:8px 10px; margin:6px 0 0;",
          bsicons::bs_icon("dash-circle"), " ",
          strong("No district-level signal for this combination. "),
          sprintf("Every district returns about %s, so the map is one flat colour. ",
                  fmt_pct(stats::median(v))),
          "The routinely available indicators carry nothing that separates ",
          "districts for this nutrient here. Use the national figure and do ",
          "not rank districts on it.")
    })

    # Offer the Admin-3 (chiefdom) level only for countries that have an
    # Admin-3 layer built (currently Sierra Leone). Update the choices when the
    # country changes, preserving the current selection where still valid.
    observeEvent(input$country, {
      base <- c("Admin 2 (district)" = "admin2", "Admin 1 (region)" = "admin1")
      has_a3 <- !is.null(admin3_pred) && (input$country %in% admin3_countries)
      choices <- if (has_a3)
        c(base, "Admin 3 (chiefdom)" = "admin3") else base
      sel <- if ((input$admin_level %||% "admin2") %in% choices)
        input$admin_level else "admin2"
      updateRadioButtons(session, "admin_level", choices = choices,
                         selected = sel, inline = TRUE)
    })

    # Reactive: joined sf for selected country × outcome × admin level
    map_data <- reactive({
      req(input$country, input$outcome, input$admin_level)
      pd <- pred_df()
      if (input$admin_level == "admin3") {
        get_country_admin3(input$country, input$outcome,
                            admin3_bnds, admin3_pred)
      } else if (input$admin_level == "admin1") {
        get_country_admin1(input$country, input$outcome,
                            admin1_bnds, admin2_bnds, pd, admin2_pop)
      } else {
        get_country_admin2(input$country, input$outcome,
                            admin2_bnds, pd, admin2_pop)
      }
    })

    # Convenience: the area-name column for the active level
    area_col <- reactive(
      switch(input$admin_level %||% "admin2",
             admin1 = "Admin1", admin3 = "Admin3", "Admin2")
    )

    # Per-outcome biomarker / data caveat shown under the outcome selector.
    output$outcome_caveat <- renderUI({
      cv <- biomarker_caveats[[input$outcome %||% ""]]
      if (is.null(cv)) return(NULL)
      div(style = paste("font-size:0.78em; color:#8a6d3b; background:#fcf8e3;",
                        "border-left:3px solid #e0c97f; padding:5px 8px;",
                        "border-radius:3px; margin:-2px 0 6px;"),
          bsicons::bs_icon("info-circle"), " ", cv)
    })

    # National headline
    output$headline <- renderUI({
      df <- map_data()
      req(df)
      natl <- national_aggregate(df)
      ctry_label <- meta$countries[input$country]
      oc_label <- meta$outcome_labels[input$outcome]

      is_a3 <- isTRUE(input$admin_level == "admin3")

      tagList(
        h5("National summary", style = "margin-top: 0;"),
        p(strong(ctry_label), " | ", em(oc_label)),
        p(
          "Estimated prevalence: ",
          tags$span(style = "font-size: 1.2em; color: #d7191c;",
                    strong(fmt_pct(natl$pred_prev_natl))),
          if (is_a3) tags$span(style = "font-size:0.82em; color:#777;",
                               " (unweighted mean across chiefdoms)")
        ),
        if (is_a3) {
          p(tags$span(
            style = "font-size: 0.82em; color: #777;",
            "Chiefdom (Admin-3) layer: 153 chiefdoms via a Fay-Herriot / ",
            "empirical-Bayes area model, narrow where surveyed and wider ",
            "elsewhere. Population counts are not available at chiefdom level."))
        } else {
          p(
            "Population affected: ",
            strong(fmt_count(natl$pop_at_risk_natl)),
            " of ", fmt_count(natl$pop_total), " ",
            if (startsWith(input$outcome, "child_")) "children 6–59 months"
            else "women 15–49 years",
            tags$br(),
            tags$span(
              style = "font-size: 0.82em; color: #777;",
              if ((input$pred_model %||% DEFAULT_PRED_MODEL) == "sl")
                "Covers surveyed districts only — pick an all-district model under Advanced."
              else
                sprintf("Covers every district. Model: %s.",
                        tolower(PRED_MODEL_INFO[[input$pred_model %||%
                                                 DEFAULT_PRED_MODEL]]$label)))
          )
        },
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

      if (layer == "misclass") {
        # WHO-class certainty: is the 95% interval wide enough to cross a WHO
        # severity threshold (i.e. could the district's class change)?
        th <- meta$who_thresholds[[input$outcome]]
        cuts <- if (!is.null(th)) sort(unique(as.numeric(th))) else numeric(0)
        spans <- vapply(seq_len(nrow(df)), function(i) {
          lo <- df$ci_lo[i]; hi <- df$ci_hi[i]
          if (is.na(lo) || is.na(hi) || length(cuts) == 0) return(NA_integer_)
          as.integer(any(cuts > lo & cuts < hi))
        }, integer(1))
        df$misclass <- ifelse(is.na(spans), "No interval",
                       ifelse(spans == 1L, "Class uncertain", "Class stable"))
        mc_cols <- c("Class stable" = "#1a9850", "Class uncertain" = "#d7191c",
                     "No interval" = "#cccccc")
        fill_col <- unname(mc_cols[df$misclass])
        pal <- NULL
        legend_title <- "WHO-class certainty"
      } else if (layer == "who_class") {
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
      # Border encodes whether a district has actual survey data (vs modeled).
      surveyed   <- !is.na(df$obs_prev)
      border_col <- ifelse(surveyed, "#1a1a1a", "#9aa0a6")
      border_wt  <- ifelse(surveyed, 1.8, 0.5)
      cert_line <- if (layer == "misclass")
        sprintf("Class certainty: %s<br/>", df$misclass) else rep("", nrow(df))
      labels <- sprintf(
        "<strong>%s</strong><br/>%s<br/>%sPopulation: %s<br/>WHO class: %s<br/>%sSurvey data: %s",
        area_name, sub_line,
        diff_line,
        fmt_count(df$population),
        df$who_class,
        cert_line,
        ifelse(surveyed, "yes", "no (modeled)")
      ) |> lapply(HTML)

      m <- leaflet(df) |>
        addProviderTiles(providers$CartoDB.Positron) |>
        addPolygons(
          fillColor = fill_col,
          fillOpacity = 0.75,
          color = border_col,
          weight = border_wt,
          opacity = 1,
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
        ) |>
        addControl(
          position = "bottomleft",
          html = HTML(paste0(
            "<div style='background:rgba(255,255,255,0.85);padding:4px 8px;",
            "border-radius:4px;font-size:11px;line-height:1.5;'>",
            "<span style='display:inline-block;width:16px;border-top:3px solid #1a1a1a;",
            "vertical-align:middle;'></span> has survey data<br/>",
            "<span style='display:inline-block;width:16px;border-top:2px solid #9aa0a6;",
            "vertical-align:middle;'></span> modeled only</div>")))

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
      } else if (layer == "misclass") {
        m <- m |> addLegend(
          colors = unname(mc_cols),
          labels = names(mc_cols),
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
      level_label <- switch(input$admin_level %||% "admin2",
                            admin1 = "region", admin3 = "chiefdom", "district")
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
        if (!is.na(row$population[1]))
          p(strong("Population: "), fmt_count(row$population[1])),
        if (!is.na(row$pop_at_risk[1]))
          p(strong("Population at risk: "), fmt_count(row$pop_at_risk[1]))
      )
    })

    output$map_caption <- renderText({
      ctry_label <- meta$countries[input$country]
      sy <- meta$survey_years[input$country]
      lvl <- switch(input$admin_level %||% "admin2",
                    admin1 = "Admin-1 (region)",
                    admin3 = "Admin-3 (chiefdom)",
                    "Admin-2 (district)")
      layer_note <- switch(input$layer,
        diff_prev = " Difference layer: modeled (within-country) minus survey-observed prevalence; red = model over-predicts, blue = under-predicts.",
        loco_diff = " Transportability layer: prevalence predicted by the universal model trained ONLY on the other countries, minus the local survey, using a harmonized cross-country outcome (uniform WHO cutoffs: ferritin<12/15 = iron deficiency, RBP<0.70 = VAD). Red = over-predicts, blue = under-predicts. Gray = country/outcome not covered.",
        ""
      )
      area_word <- switch(input$admin_level %||% "admin2",
                          admin1 = "Regions", admin3 = "Chiefdoms", "Districts")
      agg_note <- switch(input$admin_level %||% "admin2",
        admin1 = " Admin-1 values are population-weighted aggregates of Admin-2 predictions.",
        admin3 = " Chiefdom (Admin-3) values come from a Fay-Herriot / empirical-Bayes area model on chiefdom geospatial indicators; no population denominator at this level.",
        "")
      sprintf(
        "Data: %s, survey year %s. %s view. %s in gray have no model coverage.%s%s",
        ctry_label, sy, lvl, area_word, agg_note, layer_note
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
