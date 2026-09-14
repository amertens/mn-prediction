# =============================================================================
# Module: What the ranking buys
# =============================================================================
# The ranking translated into programme terms: how much of a country's
# deficiency burden the worst-ranked fifth of districts holds, whether the
# model's confidence means what it says, and how often a district lands in the
# right WHO band. All from the protocol's out-of-fold targeting tests.

mod_targeting_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Burden reached", icon = bsicons::bs_icon("bullseye"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Share of a country's deficient people in the fifth of districts each method picks"),
             plotlyOutput(ns("burden"), height = "340px")),
        card(card_header("By outcome and country"),
             selectInput(ns("outcome"), "Outcome", choices = outcome_choices, selected = "child_vitA"),
             plotlyOutput(ns("burden_cells"), height = "280px"))),
      methods_note(sprintf(paste("Direct effort at the fifth of districts the model ranks worst and you reach %s of the country's deficient",
                                 "people; the survey's own regional averages reach %s, no information %s, and perfect knowledge %s. So the",
                                 "model closes about a fifth of the gap between guessing and knowing. Burden is spread across districts,",
                                 "which is why even the true worst fifth holds under half of it. Under country transport the burden",
                                 "margin disappears: the transported product is an ordering, and should not be restated as a targeting gain."),
                           fmt_pct(Q$cap_index, 0), fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_null, 0), fmt_pct(Q$cap_oracle, 0)))
    ),
    nav_panel(
      title = "Is the model's confidence honest?", icon = bsicons::bs_icon("patch-check"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Districts the model puts in the worst fifth with a given confidence: how many were there in the survey"),
             plotlyOutput(ns("calibration"), height = "340px")),
        card(card_body(uiOutput(ns("calibration_text"))))),
      methods_note("The model was refitted on 40 random splits, each surveyed district hidden in turn, and the chance of the",
                   " worst fifth is the share of refits in which the district landed there. Districts are grouped by that",
                   " chance and compared with the survey's own worst fifth. The survey's worst fifth is itself noisy at one",
                   " cluster per district, so 100 percent agreement is not the ceiling.")
    ),
    nav_panel(
      title = "WHO severity bands", icon = bsicons::bs_icon("grid-3x3"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Share of districts placed in the correct WHO vitamin A band, and within one band"),
             reactableOutput(ns("bands"))),
        card(card_body(p(sprintf(paste("Classified into the WHO public-health bands for vitamin A, the model places a district in the",
                                       "correct band %s of the time and within one band %s. A ranking transported to a new country and",
                                       "anchored to that country's national prevalence does nearly as well, which is the form a",
                                       "programme in an unsurveyed country would receive."), fmt_pct(Q$band_exact, 0), fmt_pct(Q$band_w1, 0)))))),
      methods_note("Vitamin A bands: under 2 percent, 2 to 10, 10 to 20, 20 and above. Out-of-fold predictions averaged over ten",
                   " draws; the transported rows use the held-out country's own national prevalence as the anchor. Malawi's",
                   " districts all fall in the lowest band, where the classification is trivial, and are excluded.")
    )
  )
}

mod_targeting_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    lab <- c(oracle_ceiling = "Perfect knowledge", domain_index = "Proxy index", spatial_plus_domain = "Neighbour smoother + proxies",
             spatial = "Neighbour smoother", region_mean_jk = "Survey's regional averages", null_train_mean = "No information")
    cols <- c("Perfect knowledge" = "#8c8c8c", "Proxy index" = PROXY_COL, "Neighbour smoother + proxies" = "#7fb3bd",
              "Neighbour smoother" = "#a9c9cf", "Survey's regional averages" = SURVEY_COL, "No information" = "#bdbdbd")

    output$burden <- renderPlotly({
      NT <- EV$targeting_summary; validate(need(!is.null(NT), "Targeting summary not built."))
      d <- NT[NT$estimand == "infill" & NT$arm %in% names(lab), ]; d$Model <- lab[d$arm]
      d <- d[order(d$mean_capture), ]; d$Model <- factor(d$Model, levels = d$Model)
      plot_ly(d, x = ~mean_capture, y = ~Model, type = "bar", orientation = "h", marker = list(color = unname(cols[as.character(d$Model)])),
              text = ~sprintf("%s<br>%.0f%% of burden reached (%d cells)<br>prevalence in the picked fifth %.0f%% against %.0f%% nationally",
                              Model, 100 * mean_capture, cells, mean_prev_top20, mean_prev_national), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of deficient people reached", tickformat = ".0%", range = c(0, 0.55)), yaxis = list(title = ""),
               shapes = list(list(type = "line", x0 = 0.2, x1 = 0.2, y0 = -0.5, y1 = nrow(d) - 0.5, line = list(dash = "dash", color = "#666"))),
               margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$burden_cells <- renderPlotly({
      TCc <- EV$targeting_cells; validate(need(!is.null(TCc), "Targeting cells not built."))
      d <- TCc[TCc$estimand == "infill" & TCc$outcome == input$outcome & TCc$arm %in% c("domain_index", "region_mean_jk", "oracle_ceiling"), ] |>
        group_by(country, arm) |> summarise(v = mean(capture_top20, na.rm = TRUE), .groups = "drop")
      validate(need(nrow(d) > 0, "This outcome has no targeting cells."))
      d$Model <- lab[d$arm]
      plot_ly(d, x = ~v, y = ~country, color = ~Model, colors = cols, type = "bar", orientation = "h",
              text = ~sprintf("%s, %s: %.0f%%", country, Model, 100 * v), hoverinfo = "text") |>
        layout(barmode = "group", xaxis = list(title = "Share of deficient people reached", tickformat = ".0%"), yaxis = list(title = ""),
               legend = list(orientation = "h", y = -0.3), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$calibration <- renderPlotly({
      CAL <- EV$worst_fifth_calibration; validate(need(!is.null(CAL), "Calibration table not built."))
      CAL$band <- factor(CAL$band, levels = CAL$band)
      plot_ly(CAL, x = ~band, y = ~share_in_survey_worst_fifth, type = "bar", marker = list(color = PROXY_COL),
              text = ~sprintf("%s confidence: %d districts, %.0f%% in the survey's worst fifth", band, districts, 100 * share_in_survey_worst_fifth), hoverinfo = "text") |>
        layout(xaxis = list(title = "Model's chance that the district is in the worst fifth"),
               yaxis = list(title = "Share actually in the survey's worst fifth", tickformat = ".0%", range = c(0, 0.6)),
               shapes = list(list(type = "line", x0 = -0.5, x1 = nrow(CAL) - 0.5, y0 = 0.2, y1 = 0.2, line = list(dash = "dash", color = "#666"))),
               margin = list(l = 10, r = 10, t = 10, b = 60)) |> config(displayModeBar = FALSE)
    })

    output$calibration_text <- renderUI({
      d <- idx_districts
      tagList(
        p(sprintf(paste("Across all countries and outcomes, districts the model puts at 80 percent or above are in the survey's worst",
                        "fifth %s of the time, against %s for districts at 20 percent or below and a base rate of 20 percent."),
                  fmt_pct(Q$cal_top, 0), fmt_pct(Q$cal_low, 0))),
        p(sprintf("In this build %d surveyed districts have a chance of 80 percent or more and %d a chance of 20 percent or less; the rest are the uncertain middle.",
                  sum(d$p_worst_fifth >= 0.8, na.rm = TRUE), sum(d$p_worst_fifth <= 0.2, na.rm = TRUE))),
        p("This is the form a priority list should take: the firm districts named, the uncertain middle shown as uncertain.")
      )
    })

    output$bands <- renderReactable({
      RC <- EV$risk_summary; validate(need(!is.null(RC), "Band accuracy not built."))
      lab2 <- c(domain_index = "Proxy index, district hidden", spatial_plus_domain = "Neighbour smoother + proxies, district hidden",
                spatial = "Neighbour smoother, district hidden", region_mean_jk = "Survey's regional average, district hidden",
                transport_anchored_full = "Transported ranking + national anchor, full index",
                transport_anchored_climate_soil = "Transported ranking + national anchor, climate and soil")
      d <- RC[RC$scheme == "who_vitA" & RC$arm %in% names(lab2), ]
      t <- data.frame(Method = lab2[d$arm], `Correct band` = fmt_pct(d$exact_admin2, 0), `Within one band` = fmt_pct(d$within1_admin2, 0), Cells = d$cells, check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE)
    })
  })
}
