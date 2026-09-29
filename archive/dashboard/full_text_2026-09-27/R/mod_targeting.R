# =============================================================================
# Module: What the ranking buys
# =============================================================================
# The ranking in programme terms: how much of a country's deficiency burden the
# worst-ranked fifth of districts holds, how far to trust the worst-fifth
# signal, and how often a district lands in the right WHO band. All from the
# tests with districts hidden.

mod_targeting_ui <- function(id) {
  ns <- NS(id)
  gap <- if (all(is.finite(c(Q$cap_index, Q$cap_null, Q$cap_oracle))))
    (Q$cap_index - Q$cap_null) / (Q$cap_oracle - Q$cap_null) else NA_real_
  navset_card_tab(
    nav_panel(
      title = "Burden reached", icon = bsicons::bs_icon("bullseye"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Share of a country's deficient people living in the fifth of districts each method ranks worst"),
             plotlyOutput(ns("burden"), height = "340px")),
        card(card_header("By outcome and country"),
             selectInput(ns("outcome"), "Outcome", choices = outcome_choices, selected = "child_vitA"),
             plotlyOutput(ns("burden_cells"), height = "280px"))),
      methods_note(sprintf(paste("A programme directed to the fifth of districts the model ranks worst would reach %s of a",
                                 "country's deficient people. The fifth chosen with the survey's regional averages reaches %s,",
                                 "a random fifth %s, and the true worst fifth %s. The model therefore closes about %s of the gap",
                                 "between no information and perfect knowledge: likely better than untargeted, but far from",
                                 "complete. Deficiency is spread across many districts, which is why even the true worst fifth",
                                 "holds under half of it. In a country left out of the model this gain disappears, so there the",
                                 "ranking should not be presented as a targeting gain."),
                           fmt_pct(Q$cap_index, 0), fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_null, 0), fmt_pct(Q$cap_oracle, 0),
                           fmt_pct(gap, 0)))
    ),
    nav_panel(
      title = "How reliable is the worst-fifth signal?", icon = bsicons::bs_icon("patch-check"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Districts grouped by how often the model placed them in the worst fifth, and how many were in the survey's worst fifth"),
             plotlyOutput(ns("calibration"), height = "340px")),
        card(card_body(uiOutput(ns("calibration_text"))))),
      methods_note("The model was re-estimated 40 times on different random splits of the surveyed districts, each time",
                   " with the district hidden, and we counted how often it landed in the worst fifth. Districts are grouped",
                   " by that share and compared with the survey's own worst fifth.")
    ),
    nav_panel(
      title = "WHO severity bands", icon = bsicons::bs_icon("grid-3x3"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Share of districts placed in the correct WHO vitamin A band, and within one band"),
             reactableOutput(ns("bands"))),
        card(card_body(p(sprintf(paste("Using the WHO public health bands for vitamin A, the model places %s of districts in",
                                       "the correct band and %s within one band. The survey's own regional averages do slightly",
                                       "better (%s and %s), so the model adds most where regional averages are not available.",
                                       "A model built without the country, combined with that country's national prevalence,",
                                       "does nearly as well; that is the situation of a country without district data."),
                                 fmt_pct(Q$band_exact, 0), fmt_pct(Q$band_w1, 0), fmt_pct(Q$band_exact_jk, 0), fmt_pct(Q$band_w1_jk, 0)))))),
      methods_note("Vitamin A bands: under 2%, 2% to 10%, 10% to 20%, and 20% or more. Predictions are for districts",
                   " the model had not seen, averaged over ten runs; for countries left out of the model, that country's own",
                   " national prevalence sets the level. Malawi's districts all fall in the lowest band, where classification",
                   " is trivial, so they are excluded.")
    )
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
              text = ~sprintf("%s<br>%.0f%% of deficient people reached (%d country-outcome pairs)<br>prevalence in the chosen fifth %.0f%%, against %.0f%% nationally",
                              Model, 100 * mean_capture, cells, mean_prev_top20, mean_prev_national), hoverinfo = "text") |>
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
              text = ~sprintf("Placed in the worst fifth in %s of runs: %d districts, %.0f%% in the survey's worst fifth", band, districts, 100 * share_in_survey_worst_fifth), hoverinfo = "text") |>
        layout(xaxis = list(title = "How often the model placed the district in the worst fifth"),
               yaxis = list(title = "Share actually in the survey's worst fifth", tickformat = ".0%", range = c(0, 0.6)),
               shapes = list(list(type = "line", x0 = -0.5, x1 = nrow(CAL) - 0.5, y0 = 0.2, y1 = 0.2, line = list(dash = "dash", color = "#666"))),
               margin = list(l = 10, r = 10, t = 10, b = 60)) |> config(displayModeBar = FALSE)
    })

    output$calibration_text <- renderUI({
      d <- idx_districts
      tagList(
        p(sprintf(paste("Districts the model placed in the worst fifth in at least 80%% of runs were in the survey's worst fifth",
                        "%s of the time. Districts placed there in at most 20%% of runs were in it %s of the time. By chance it",
                        "would be 20%%."),
                  fmt_pct(Q$cal_top, 0), fmt_pct(Q$cal_low, 0))),
        p("So the signal separates districts in the right direction, but it overstates certainty: a district placed in the",
          " worst fifth in 80% of runs does not have an 80% chance of being there. Part of the gap is because the survey's",
          " own worst fifth is uncertain when districts have only one or two clusters."),
        p(sprintf("In this data, %d surveyed districts are placed in the worst fifth in at least 80%% of runs and %d in at most 20%%; the rest are uncertain.",
                  sum(d$p_worst_fifth >= 0.8, na.rm = TRUE), sum(d$p_worst_fifth <= 0.2, na.rm = TRUE)))
      )
    })

    output$bands <- renderReactable({
      RC <- EV$risk_summary; validate(need(!is.null(RC), "Band accuracy not built."))
      lab2 <- c(domain_index = "This model, district hidden", spatial_plus_domain = "Neighbouring districts + this model, district hidden",
                spatial = "Neighbouring districts' average, district hidden", region_mean_jk = "Survey's regional averages, district hidden",
                transport_anchored_full = "Country left out: all data layers + national prevalence",
                transport_anchored_climate_soil = "Country left out: climate and soil + national prevalence")
      d <- RC[RC$scheme == "who_vitA" & RC$arm %in% names(lab2), ]
      t <- data.frame(Method = lab2[d$arm], `Correct band` = fmt_pct(d$exact_admin2, 0), `Within one band` = fmt_pct(d$within1_admin2, 0),
                      `Country-outcome pairs` = d$cells, check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE)
    })
  })
}
