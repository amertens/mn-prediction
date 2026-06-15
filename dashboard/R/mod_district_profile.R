# =============================================================================
# Module: District Profile Cards
# =============================================================================
# Detailed view for one district at a time, showing all available outcomes,
# population breakdown, and comparison to country average. This is the
# primary deliverable for sub-national health officers.

mod_district_profile_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Select district",

      selectInput(ns("country"), "Country",
                  choices = country_choices,
                  selected = "ghana"),

      radioButtons(ns("dp_model"), "Prediction model",
                   choices = c("Fay-Herriot SAE (all districts, with intervals)" = "fh",
                               "Area-level SAE — HAL (all districts)" = "area",
                               "Individual SuperLearner (surveyed districts)" = "sl"),
                   selected = "fh"),

      uiOutput(ns("district_picker")),

      hr(),

      htmlOutput(ns("district_basics"))
    ),

    layout_columns(
      col_widths = 12,

      card(
        card_header("Predicted micronutrient deficiency profile"),
        card_body(
          plotlyOutput(ns("outcome_plot"), height = "400px"),
          methods_note(
            "Each row shows the model-predicted prevalence (point) for one ",
            "micronutrient outcome in this district, with its 95% interval ",
            "where available (horizontal line). The dashed vertical line marks the ",
            "country-level average for comparison. Outcomes where the district ",
            "estimate is meaningfully above the country average indicate ",
            "elevated local burden that may warrant targeted intervention. ",
            tags$br(), tags$br(),
            "Intervals, where shown, come from split-conformal prediction on ",
            "cross-validated residuals. Where an interval is not available for a ",
            "district, only the point estimate is shown — the Map explorer's ",
            "Fay-Herriot layer provides an interval for every district."
          )
        )
      ),

      card(
        card_header("Outcome detail"),
        card_body(
          reactableOutput(ns("outcome_table")),
          downloadButton(ns("download_district"), "Download district CSV",
                          class = "btn-sm btn-outline-primary",
                          style = "margin-top: 8px;"),
          methods_note(
            "“Predicted” is the model output for this district; “Observed” is ",
            "the survey-weighted prevalence from the original biomarker survey ",
            "(only available where the survey actually sampled this district — ",
            "blank for unsurveyed districts). Population-affected counts use ",
            "WorldPop population estimates for the survey year multiplied by ",
            "predicted prevalence; they should be interpreted as planning ",
            "estimates rather than census-precise counts. ",
            tags$br(), tags$br(),
            "WHO classifications follow the standard public health thresholds ",
            "for severity of micronutrient deficiency at the population level. ",
            "These thresholds are designed for prevalence interpretation, not ",
            "individual diagnosis: a district classified “Severe” has prevalence ",
            "high enough to warrant population-level intervention."
          )
        )
      )
    )
  )
}


mod_district_profile_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Selected prediction source (mirrors the Map explorer's model toggle).
    dp_pred <- reactive({
      m <- input$dp_model %||% "fh"
      if (m == "fh" && exists("admin2_fh_pred") && !is.null(admin2_fh_pred)) {
        admin2_fh_pred
      } else if (m == "area" && exists("admin2_area_pred") && !is.null(admin2_area_pred)) {
        admin2_area_pred
      } else {
        admin2_pred
      }
    })

    # Dynamic district picker
    output$district_picker <- renderUI({
      req(input$country)
      ctry_label <- meta$countries[input$country]
      pd <- dp_pred()
      districts <- sort(unique(pd$Admin2[pd$country == ctry_label]))
      if (length(districts) == 0) {
        return(p(em("No districts available for this country."),
                 style = "color: #888;"))
      }
      selectInput(ns("district"), "District",
                  choices = districts,
                  selected = districts[1])
    })

    district_data <- reactive({
      req(input$country, input$district)
      ctry_label <- meta$countries[input$country]

      pd <- dp_pred()
      df <- pd[pd$country == ctry_label & pd$Admin2 == input$district, ,
               drop = FALSE]
      if (!"Admin1" %in% colnames(df)) df$Admin1 <- NA_character_

      pop_row <- admin2_pop[admin2_pop$country == ctry_label &
                              admin2_pop$Admin2 == input$district, ,
                              drop = FALSE]
      pop <- if (nrow(pop_row) > 0) pop_row$population[1] else NA_real_

      df$population <- pop
      df$pop_affected <- df$pred_prev * pop
      df$Outcome_label <- meta$outcome_labels[df$outcome]

      df
    })

    # National averages (for comparison line in plot)
    natl_avgs <- reactive({
      req(input$country)
      ctry_label <- meta$countries[input$country]

      pd <- dp_pred()
      avgs <- list()
      for (oc_key in unique(pd$outcome[pd$country == ctry_label])) {
        df <- get_country_admin2(input$country, oc_key,
                                  admin2_bnds, pd, admin2_pop)
        if (is.null(df)) next
        natl <- national_aggregate(df)
        avgs[[oc_key]] <- natl$pred_prev_natl
      }
      avgs
    })

    output$district_basics <- renderUI({
      df <- district_data()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      tagList(
        h5(input$district, style = "margin-top: 0.5em;"),
        if (!is.null(df$Admin1) && length(df$Admin1) > 0 &&
            !is.na(df$Admin1[1])) p(em(df$Admin1[1])),
        p(strong("Population: "), fmt_count(df$population[1]))
      )
    })

    output$outcome_plot <- renderPlotly({
      df <- district_data()
      avgs <- natl_avgs()
      req(nrow(df) > 0)

      df$natl_avg <- vapply(df$outcome, function(o) {
        if (is.null(avgs[[o]])) NA_real_ else avgs[[o]]
      }, numeric(1))

      # Numeric y positions (consistent across layers) with outcome labels as
      # tick text. Mixing numeric and categorical y axes — or drawing all-NA
      # interval segments — crashes plotly, so we avoid both.
      df <- df[order(df$pred_prev), , drop = FALSE]
      df$ypos <- seq_len(nrow(df))
      has_ci <- any(is.finite(df$ci_lo) & is.finite(df$ci_hi))
      xmax <- suppressWarnings(max(c(df$pred_prev, df$ci_hi, df$natl_avg),
                                   na.rm = TRUE))
      if (!is.finite(xmax) || xmax <= 0) xmax <- 0.5

      p <- plot_ly()
      avg_df <- df[is.finite(df$natl_avg), , drop = FALSE]
      if (nrow(avg_df) > 0) {
        p <- add_segments(p, data = avg_df,
          x = ~natl_avg, xend = ~natl_avg,
          y = ~ypos - 0.3, yend = ~ypos + 0.3,
          line = list(color = "#888", width = 1, dash = "dash"),
          name = "Country average", hoverinfo = "text",
          text = ~paste("Country avg:", fmt_pct(natl_avg)))
      }
      if (has_ci) {
        ci_df <- df[is.finite(df$ci_lo) & is.finite(df$ci_hi), , drop = FALSE]
        p <- add_segments(p, data = ci_df,
          x = ~ci_lo, xend = ~ci_hi, y = ~ypos, yend = ~ypos,
          line = list(color = "#2c7bb6", width = 4),
          name = "95% interval", hoverinfo = "text",
          text = ~sprintf("Interval: %s – %s", fmt_pct(ci_lo), fmt_pct(ci_hi)))
      }
      p <- add_markers(p, data = df,
        x = ~pred_prev, y = ~ypos,
        marker = list(color = "#d7191c", size = 12,
                      line = list(color = "white", width = 1.5)),
        name = "Predicted", hoverinfo = "text",
        text = ~sprintf("%s<br>%s<br>WHO: %s",
                        Outcome_label, fmt_pct(pred_prev), who_class))

      p |> layout(
        xaxis = list(title = "Prevalence", tickformat = ".0%",
                     range = c(0, xmax * 1.08)),
        yaxis = list(title = "", tickvals = df$ypos,
                     ticktext = as.character(df$Outcome_label)),
        margin = list(l = 200),
        legend = list(orientation = "h", y = -0.15)
      ) |>
        config(displayModeBar = FALSE)
    })

    output$outcome_table <- renderReactable({
      df <- district_data()
      req(nrow(df) > 0)

      out <- data.frame(
        Outcome    = df$Outcome_label,
        Predicted  = df$pred_prev,
        CI_low     = df$ci_lo,
        CI_high    = df$ci_hi,
        Observed   = df$obs_prev,
        WHO_class  = df$who_class,
        Affected   = df$pop_affected,
        stringsAsFactors = FALSE
      )

      reactable(
        out,
        compact = TRUE, striped = TRUE,
        defaultPageSize = 10,
        columns = list(
          Outcome    = colDef(name = "Outcome", minWidth = 200),
          Predicted  = colDef(name = "Predicted",
                              format = colFormat(percent = TRUE, digits = 1)),
          CI_low     = colDef(name = "CI low",
                              format = colFormat(percent = TRUE, digits = 1)),
          CI_high    = colDef(name = "CI high",
                              format = colFormat(percent = TRUE, digits = 1)),
          Observed   = colDef(name = "Survey observed",
                              format = colFormat(percent = TRUE, digits = 1)),
          WHO_class  = colDef(name = "WHO class",
                              style = function(value) {
                                if (is.na(value)) return(NULL)
                                col <- who_colors[value]
                                if (is.na(col)) return(NULL)
                                list(color = "white", background = unname(col),
                                     fontWeight = "bold", textAlign = "center")
                              }),
          Affected   = colDef(name = "Pop. affected",
                              format = colFormat(separators = TRUE, digits = 0))
        )
      )
    })

    output$download_district <- downloadHandler(
      filename = function() {
        sprintf("district_%s_%s_%s.csv",
                input$country,
                gsub("[^A-Za-z0-9]", "_", input$district),
                Sys.Date())
      },
      content = function(file) {
        df <- district_data()
        write.csv(df, file, row.names = FALSE)
      }
    )

  })
}
