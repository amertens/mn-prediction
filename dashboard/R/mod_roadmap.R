# =============================================================================
# Module: What more data buys (concise version)
# =============================================================================

mod_roadmap_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    p(class = "lead", "What would make the estimates better:"),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE,
        card_header("Better surveys: the best achievable score and what the model reaches"),
        card_body(plotlyOutput(ns("headroom"), height = "520px"),
          p(style = "font-size:0.85em; color:#666;",
            "Black marks: the best score any model could reach, given noise in the survey's own district figures.",
            " More clusters per district would raise them."))),
      div(
        card(card_header("More surveys"),
          card_body(plotlyOutput(ns("curve"), height = "220px"),
            p(style = "font-size:0.85em; color:#555;",
              sprintf("Accuracy in a country left out: %s with one training survey, %s with three.", f2(Q$lc1), f2(Q$lc_max))))),
        card(card_header("Other regions"),
          card_body(p(style = "font-size:0.9em;",
            sprintf("The ranking also worked in Pakistan and India (%s). A global soil layer works as well as the African one.", f2(Q$xv_off)),
            " ", actionLink(ns("go_external"), "See the test."))))
      )
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header("Related biomarkers"),
        card_body(p(style = "font-size:0.88em;",
          sprintf("Adding the same nutrient's biomarker from the other population group raised accuracy from %s to %s.",
                  f2(Q$xo_base), f2(Q$xo_same))))),
      card(card_header("New outcomes"),
        card_body(p(style = "font-size:0.88em;",
          sprintf("Selenium, measured in Malawi, ranks at %s to %s, better than the survey's regional averages.",
                  f2(Q$se_idx_lo), f2(Q$se_idx_hi))))),
      card(card_header("Tested in advance"),
        card_body(p(style = "font-size:0.88em;",
          "Two versions of the model are registered for testing on the next national survey.")))
    )
  )
}

mod_roadmap_server <- function(id, go_to = NULL) {
  moduleServer(id, function(input, output, session) {
    output$headroom <- renderPlotly({
      HR <- EV$headroom; validate(need(!is.null(HR), "Table not built."))
      h <- HR[is.finite(HR$r_max_emp), ]
      h$oc <- ifelse(h$outcome %in% names(outcome_short), outcome_short[h$outcome], h$outcome)
      h$cell <- paste(h$country, "-", h$oc)
      h <- h[order(h$r_max_emp), ]; h$cell <- factor(h$cell, levels = h$cell)
      plot_ly(h) |>
        add_segments(x = ~r_max_emp_lo, xend = ~r_max_emp_hi, y = ~cell, yend = ~cell,
                     line = list(color = "#d0d0d0", width = 5), hoverinfo = "none", showlegend = FALSE) |>
        add_markers(x = ~r_max_emp, y = ~cell, marker = list(symbol = "line-ns-open", size = 14, color = "#1a1a1a", line = list(width = 2)),
                    text = ~sprintf("%s<br>best achievable score %.2f", cell, r_max_emp), hoverinfo = "text", name = "Best achievable score") |>
        add_markers(x = ~spearman, y = ~cell, marker = list(color = PROXY_COL, size = 9),
                    text = ~sprintf("%s<br>model reaches %.2f", cell, spearman), hoverinfo = "text", name = "What the model reaches") |>
        layout(xaxis = list(title = "Ranking accuracy", zeroline = TRUE),
               yaxis = list(title = "", tickfont = list(size = 10)),
               legend = list(orientation = "h", y = -0.08), margin = list(l = 10, r = 10, t = 10, b = 40)) |>
        config(displayModeBar = FALSE)
    })
    output$curve <- renderPlotly({
      TC <- EV$training_curve; validate(need(!is.null(TC), "Table not built."))
      a <- TC[TC$arm == "domain_index" & TC$target == "level", ]
      d <- a |> group_by(n_train_countries) |> summarise(m = mean(spearman, na.rm = TRUE), .groups = "drop")
      plot_ly(d, x = ~n_train_countries, y = ~m, type = "scatter", mode = "lines+markers",
              line = list(color = PROXY_COL), marker = list(color = PROXY_COL, size = 10),
              text = ~sprintf("%d training surveys: %.2f", n_train_countries, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Number of training surveys", dtick = 1),
               yaxis = list(title = "Accuracy, country left out", rangemode = "tozero"),
               margin = list(l = 10, r = 10, t = 10, b = 35)) |> config(displayModeBar = FALSE)
    })
    if (!is.null(go_to)) observeEvent(input$go_external, go_to("Tested in six more countries"))
  })
}
