# =============================================================================
# Module: What the ranking buys (concise version)
# =============================================================================

mod_targeting_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Burden reached", icon = bsicons::bs_icon("bullseye"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Share of deficient people in the fifth of districts each method ranks worst"),
             plotlyOutput(ns("burden"), height = "340px")),
        card(card_header("By outcome and country"),
             selectInput(ns("outcome"), "Outcome", choices = outcome_choices, selected = "child_vitA"),
             plotlyOutput(ns("burden_cells"), height = "280px"))),
      p(style = "font-size:0.85em; color:#666;",
        sprintf("Targeting the model's worst fifth reaches %s of deficient people: more than regional averages (%s) or no information (%s), far short of perfect knowledge (%s).",
                fmt_pct(Q$cap_index, 0), fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_null, 0), fmt_pct(Q$cap_oracle, 0)))
    ),
    nav_panel(
      title = "How reliable is the worst-fifth signal?", icon = bsicons::bs_icon("patch-check"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("How often a district was placed in the worst fifth, and how often it really was there"),
             plotlyOutput(ns("calibration"), height = "340px")),
        card(card_body(uiOutput(ns("calibration_text"))))))
    ,
    nav_panel(
      title = "WHO severity bands", icon = bsicons::bs_icon("grid-3x3"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Districts placed in the correct WHO vitamin A band"),
             reactableOutput(ns("bands"))),
        card(card_body(p(sprintf("Correct band for %s of districts and within one band for %s. The survey's regional averages: %s and %s.",
                                 fmt_pct(Q$band_exact, 0), fmt_pct(Q$band_w1, 0), fmt_pct(Q$band_exact_jk, 0), fmt_pct(Q$band_w1_jk, 0)))))))
  )
}

mod_targeting_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    lab <- c(oracle_ceiling = "Perfect knowledge", domain_index = "This model", spatial_plus_domain = "Neighbouring districts + this model",
             spatial = "Neighbouring districts' average", region_mean_jk = "Survey's regional averages", null_train_mean = "No information")
    cols <- c("Perfect knowledge" = "#8c8c8c", "This model" = PROXY_COL, "Neighbouring districts + this model" = "#7fb3bd",
              "Neighbouring districts' average" = "#a9c9cf", "Survey's regional averages" = SURVEY_COL, "No information" = "#bdbdbd")

    output$burden <- renderPlotly({
      NT <- EV$targeting_summary; validate(need(!is.null(NT), "Targeting summary not built."))
      d <- NT[NT$estimand == "infill" & NT$arm %in% names(lab), ]; d$Model <- lab[d$arm]
      d <- d[order(d$mean_capture), ]; d$Model <- factor(d$Model, levels = d$Model)
      plot_ly(d, x = ~mean_capture, y = ~Model, type = "bar", orientation = "h", marker = list(color = unname(cols[as.character(d$Model)])),
              text = ~sprintf("%s: %.0f%% of deficient people reached", Model, 100 * mean_capture), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of deficient people reached", tickformat = ".0%", range = c(0, 0.55)), yaxis = list(title = ""),
               shapes = list(list(type = "line", x0 = 0.2, x1 = 0.2, y0 = -0.5, y1 = nrow(d) - 0.5, line = list(dash = "dash", color = "#666"))),
               margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$burden_cells <- renderPlotly({
      TCc <- EV$targeting_cells; validate(need(!is.null(TCc), "Targeting results not built."))
      d <- TCc[TCc$estimand == "infill" & TCc$outcome == input$outcome & TCc$arm %in% c("domain_index", "region_mean_jk", "oracle_ceiling"), ] |>
        group_by(country, arm) |> summarise(v = mean(capture_top20, na.rm = TRUE), .groups = "drop")
      validate(need(nrow(d) > 0, "This outcome has no targeting results."))
      d$Model <- lab[d$arm]
      plot_ly(d, x = ~v, y = ~country, color = ~Model, colors = cols, type = "bar", orientation = "h",
              text = ~sprintf("%s, %s: %.0f%%", country, Model, 100 * v), hoverinfo = "text") |>
        layout(barmode = "group", xaxis = list(title = "Share of deficient people reached", tickformat = ".0%"), yaxis = list(title = ""),
               legend = list(orientation = "h", y = -0.3), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$calibration <- renderPlotly({
      CAL <- EV$worst_fifth_calibration; validate(need(!is.null(CAL), "Table not built."))
      CAL$band <- factor(CAL$band, levels = CAL$band)
      plot_ly(CAL, x = ~band, y = ~share_in_survey_worst_fifth, type = "bar", marker = list(color = PROXY_COL),
              text = ~sprintf("%s of runs: %d districts, %.0f%% in the survey's worst fifth", band, districts, 100 * share_in_survey_worst_fifth), hoverinfo = "text") |>
        layout(xaxis = list(title = "How often the model placed the district in the worst fifth"),
               yaxis = list(title = "Share in the survey's worst fifth", tickformat = ".0%", range = c(0, 0.6)),
               shapes = list(list(type = "line", x0 = -0.5, x1 = nrow(CAL) - 0.5, y0 = 0.2, y1 = 0.2, line = list(dash = "dash", color = "#666"))),
               margin = list(l = 10, r = 10, t = 10, b = 60)) |> config(displayModeBar = FALSE)
    })

    output$calibration_text <- renderUI({
      p(sprintf(paste("Districts placed in the worst fifth in at least 80%% of runs were in the survey's worst fifth %s of the",
                      "time (20%% by chance). The signal points the right way but overstates certainty."), fmt_pct(Q$cal_top, 0)))
    })

    output$bands <- renderReactable({
      RC <- EV$risk_summary; validate(need(!is.null(RC), "Band accuracy not built."))
      lab2 <- c(domain_index = "This model", region_mean_jk = "Survey's regional averages",
                transport_anchored_full = "Country left out: model + national figure")
      d <- RC[RC$scheme == "who_vitA" & RC$arm %in% names(lab2), ]
      t <- data.frame(Method = lab2[d$arm], `Correct band` = fmt_pct(d$exact_admin2, 0), `Within one band` = fmt_pct(d$within1_admin2, 0), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE)
    })
  })
}
