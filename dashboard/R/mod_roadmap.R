# =============================================================================
# Module: What more data buys
# =============================================================================
# The optimistic case, with a receipt for every sentence: the ceiling is set by
# the survey's noise rather than the method; each survey added improves every
# other country's map; the recipe now works on two continents; a related
# biomarker from any population helps; new outcomes ride the same blood draw.

mod_roadmap_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_body(
      p(class = "lead", style = "margin-bottom:4px;",
        "The accuracy this dashboard shows is not the method's limit. It is the limit of the data it was allowed to learn from",
        " — and every one of the four levers below has been measured, not assumed."),
      p(style = "color:#555; font-size:0.92em;",
        sprintf(paste("Where the surveys resolve districts well, the model already reads %s of what is there to read",
                      "(%s of %s in-country tests clear 0.5 and average %s). The rest of the gap is survey noise, and",
                      "survey noise is a budget choice."),
                "most", Q$strong_n, Q$infill_cells, f2(Q$strong_mean))))),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE,
        card_header("Lever 1 - better surveys raise the ceiling: what each cell can reach vs what the model reaches"),
        card_body(plotlyOutput(ns("headroom"), height = "520px"),
          methods_note("For each country and outcome: the teal dot is the model's ranking accuracy with each district",
                       " hidden; the black tick is the most any predictor could score against this survey's own district",
                       " values, given how few clusters they rest on (the measurability ceiling, with its 95% range in",
                       " grey). Where the dot sits near the tick, better predictors cannot help - only a survey with",
                       " more clusters per district can move the tick. Cells at or below zero have no measurable",
                       " district signal in the survey itself (Malawi zinc, vitamin A at under 1 percent prevalence)."))),
      div(
        card(card_header("Lever 2 - every survey added improves every other country's map"),
          card_body(plotlyOutput(ns("curve"), height = "240px"),
            p(style = "font-size:0.85em; color:#555;",
              sprintf(paste("Ranking accuracy in a country held out of training, by number of training surveys: %s with",
                            "one, %s with three, about %s per survey added, and the curve has not flattened. A fifth",
                            "survey is the cheapest accuracy on this page."),
                      f2(Q$lc1), f2(Q$lc_max), f2(Q$lc_step))))),
        card(card_header("Lever 3 - the recipe is no longer bound to Africa"),
          card_body(p(style = "font-size:0.9em;",
            sprintf(paste("Scored against WHO-deposited results in six countries the model never saw, the ranking reaches",
                          "%s in Africa and %s in South Asia - and a globally available soil layer matches the Africa-only",
                          "one (%s against %s on the same cells). Any country with public climate and soil data can be"),
                    f2(Q$xv_level), f2(Q$xv_off), f2(Q$xv_sg_level), f2(Q$xv_level)),
            "given a first map. ", actionLink(ns("go_external"), "See the test."))))
      )
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header("Lever 4 - any related biomarker helps"),
        card_body(p(style = "font-size:0.88em;",
          sprintf(paste("Borrowing the same nutrient's biomarker from the other population (children's iron for women's",
                        "iron, and vice versa) lifts a transported ranking from %s to %s on its own, and adding it to the",
                        "index gains about %s. Surveys that assay one population still improve the map for the other."),
                  f2(Q$xo_base), f2(Q$xo_same), f2(Q$xo_added))),
          p(style = "font-size:0.8em; color:#777;", "Cross-outcome borrowing test (XO-01), iron and vitamin A cells, country held out of training."))),
      card(card_header("New outcomes ride the same blood draw"),
        card_body(p(style = "font-size:0.88em;",
          sprintf(paste("Malawi's survey also assayed selenium: configured as a new outcome, the model ranks its districts",
                        "at %s with each district hidden - the strongest in-country signal on this dashboard, from data",
                        "that cost one extra assay on blood already drawn."), f2(Q$se_mean))),
          p(style = "font-size:0.8em; color:#777;", "Selenium follows soil geology, which is exactly what the proxy layers see."))),
      card(card_header("The next country is a prediction, not a retrofit"),
        card_body(p(style = "font-size:0.88em;",
          "Two candidate indices (climate + soil, and a five-domain extension) are pre-registered for the next national",
          " survey, with their falsification conditions written down before any data are seen. If they hold, the case for",
          " scaling is made on someone else's terms; if they fail, that is worth knowing before anyone fields a survey on it."),
          p(style = "font-size:0.8em; color:#777;", "Pre-registration: docs/findings/PREREGISTRATION_NEW_COUNTRIES_2026-09.md.")))
    )
  )
}

mod_roadmap_server <- function(id, go_to = NULL) {
  moduleServer(id, function(input, output, session) {
    output$headroom <- renderPlotly({
      HR <- EV$headroom; validate(need(!is.null(HR), "Headroom table not built."))
      h <- HR[is.finite(HR$r_max_emp), ]
      h$oc <- ifelse(h$outcome %in% names(outcome_short), outcome_short[h$outcome], h$outcome)
      h$cell <- paste(h$country, "-", h$oc)
      h <- h[order(h$r_max_emp), ]; h$cell <- factor(h$cell, levels = h$cell)
      plot_ly(h) |>
        add_segments(x = ~r_max_emp_lo, xend = ~r_max_emp_hi, y = ~cell, yend = ~cell,
                     line = list(color = "#d0d0d0", width = 5), hoverinfo = "none", showlegend = FALSE) |>
        add_markers(x = ~r_max_emp, y = ~cell, marker = list(symbol = "line-ns-open", size = 14, color = "#1a1a1a", line = list(width = 2)),
                    text = ~sprintf("%s<br>measurability ceiling %.2f (%.2f to %.2f)", cell, r_max_emp, r_max_emp_lo, r_max_emp_hi),
                    hoverinfo = "text", name = "What this survey can measure") |>
        add_markers(x = ~spearman, y = ~cell, marker = list(color = PROXY_COL, size = 9),
                    text = ~sprintf("%s<br>model reaches %.2f of a possible %.2f (%d districts)", cell, spearman, r_max_emp, n_areas),
                    hoverinfo = "text", name = "What the model reaches") |>
        layout(xaxis = list(title = "Ranking accuracy", zeroline = TRUE),
               yaxis = list(title = "", tickfont = list(size = 10)),
               legend = list(orientation = "h", y = -0.08), margin = list(l = 10, r = 10, t = 10, b = 40)) |>
        config(displayModeBar = FALSE)
    })
    output$curve <- renderPlotly({
      TC <- EV$training_curve; validate(need(!is.null(TC), "Training curve not built."))
      a <- TC[TC$arm == "domain_index" & TC$target == "level", ]
      d <- a |> group_by(n_train_countries) |>
        summarise(m = mean(spearman, na.rm = TRUE), lo = quantile(spearman, 0.25, na.rm = TRUE),
                  hi = quantile(spearman, 0.75, na.rm = TRUE), .groups = "drop")
      plot_ly(d, x = ~n_train_countries, y = ~m, type = "scatter", mode = "lines+markers",
              line = list(color = PROXY_COL), marker = list(color = PROXY_COL, size = 10),
              error_y = ~list(array = hi - m, arrayminus = m - lo, thickness = 1, color = "#9a9a9a"),
              text = ~sprintf("%d training surveys: %.2f", n_train_countries, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Training surveys", dtick = 1),
               yaxis = list(title = "Accuracy, country held out", rangemode = "tozero"),
               margin = list(l = 10, r = 10, t = 10, b = 35)) |> config(displayModeBar = FALSE)
    })
    if (!is.null(go_to)) observeEvent(input$go_external, go_to("Tested in six more countries"))
  })
}
