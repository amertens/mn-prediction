# =============================================================================
# Module: What more data buys
# =============================================================================
# What would make the estimates better, each point backed by a test on the
# existing data: better surveys raise the best achievable score, each survey
# added has raised accuracy in countries left out, the approach also worked in
# South Asia, a related biomarker helps, and new outcomes can come from blood
# already collected.

mod_roadmap_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_body(
      p(class = "lead", style = "margin-bottom:4px;",
        "The accuracy shown on this dashboard is limited by the data the model learns from. Each of the improvements below",
        " is supported by a test on the existing surveys."),
      p(style = "color:#555; font-size:0.92em;",
        sprintf(paste("In its best %s of %s within-country tests the model reaches %s on average; the other %s average %s.",
                      "The rest of this page shows what limits the weaker results and what would help."),
                Q$strong_n, Q$infill_cells, f2(Q$strong_mean), Q$weak_n, f2(Q$weak_mean))))),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE,
        card_header("Better surveys: the best achievable score for each country and outcome, and what the model reaches"),
        card_body(plotlyOutput(ns("headroom"), height = "520px"),
          methods_note("For each country and outcome, the teal dot is the model's ranking accuracy with each district hidden,",
                       " and the black mark is the best score any model could reach against the survey's district figures,",
                       " with its 95% range in grey. The best achievable score is limited by how few clusters each district",
                       " has. Where the dot is close to the mark, better data layers cannot help much; only a survey with more",
                       " clusters per district can raise the mark. Where the mark is near or below zero (Malawi zinc, and",
                       " vitamin A where under 1% of people are deficient), the survey shows no measurable difference between",
                       " districts."))),
      div(
        card(card_header("More surveys: each survey added has raised accuracy in countries left out of training"),
          card_body(plotlyOutput(ns("curve"), height = "240px"),
            p(style = "font-size:0.85em; color:#555;",
              sprintf(paste("Ranking accuracy in a country left out of training was %s with one training survey and %s with",
                            "three, about %s for each survey added. With only three points we cannot tell yet where this levels off."),
                      f2(Q$lc1), f2(Q$lc_max), f2(Q$lc_step))))),
        card(card_header("Other regions: the approach also worked in South Asia"),
          card_body(p(style = "font-size:0.9em;",
            sprintf(paste("Compared with survey results held by the WHO for six countries not used in training, the ranking",
                          "scored %s (biomarker levels) in four African countries and %s (prevalence) in Pakistan and India.",
                          "A global soil layer worked as well as the Africa-only one (%s against %s on the same tests), so the",
                          "same public data can be assembled for any country."),
                    f2(Q$xv_level), f2(Q$xv_off), f2(Q$xv_sg_level), f2(Q$xv_level)),
            " ", actionLink(ns("go_external"), "See the test."))))
      )
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header("Related biomarkers help"),
        card_body(p(style = "font-size:0.88em;",
          sprintf(paste("Using the same nutrient's biomarker from the other population group (for example, children's iron",
                        "when ranking women's iron) raised accuracy in a country left out of training from %s to %s on its",
                        "own, and adding it to the full model gained about %s. A survey that tests one group can therefore",
                        "improve the map for the other."),
                  f2(Q$xo_base), f2(Q$xo_same), f2(Q$xo_added))),
          p(style = "font-size:0.8em; color:#777;", "Iron and vitamin A tests, with the country left out of training."))),
      card(card_header("New outcomes from the same blood samples"),
        card_body(p(style = "font-size:0.88em;",
          sprintf(paste("Malawi's survey also measured selenium. Treated as a new outcome, it is ranked at %s to %s with each",
                        "district hidden, better than the survey's own regional averages (%s to %s). Selenium in food follows",
                        "soil geology, which the data layers capture."),
                  f2(Q$se_idx_lo), f2(Q$se_idx_hi), f2(Q$se_jk_lo), f2(Q$se_jk_hi))),
          p(style = "font-size:0.8em; color:#777;", "Malawi only, so this cannot yet be checked in another country."))),
      card(card_header("The next test is set in advance"),
        card_body(p(style = "font-size:0.88em;",
          "Two versions of the model (climate and soil only, and a version with five data groups) have been registered",
          " in advance for the next national survey, together with the results that would count against them. This means",
          " the approach will be judged on data it has not seen."),
          p(style = "font-size:0.8em; color:#777;", "The pre-registration document is in the project repository.")))
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
                    text = ~sprintf("%s<br>best achievable score %.2f (%.2f to %.2f)", cell, r_max_emp, r_max_emp_lo, r_max_emp_hi),
                    hoverinfo = "text", name = "Best achievable score") |>
        add_markers(x = ~spearman, y = ~cell, marker = list(color = PROXY_COL, size = 9),
                    text = ~sprintf("%s<br>model reaches %.2f of a possible %.2f (%d districts)", cell, spearman, r_max_emp, n_areas),
                    hoverinfo = "text", name = "What the model reaches") |>
        layout(xaxis = list(title = "Ranking accuracy", zeroline = TRUE),
               yaxis = list(title = "", tickfont = list(size = 10)),
               legend = list(orientation = "h", y = -0.08), margin = list(l = 10, r = 10, t = 10, b = 40)) |>
        config(displayModeBar = FALSE)
    })
    output$curve <- renderPlotly({
      TC <- EV$training_curve; validate(need(!is.null(TC), "Table not built."))
      a <- TC[TC$arm == "domain_index" & TC$target == "level", ]
      d <- a |> group_by(n_train_countries) |>
        summarise(m = mean(spearman, na.rm = TRUE), lo = quantile(spearman, 0.25, na.rm = TRUE),
                  hi = quantile(spearman, 0.75, na.rm = TRUE), .groups = "drop")
      plot_ly(d, x = ~n_train_countries, y = ~m, type = "scatter", mode = "lines+markers",
              line = list(color = PROXY_COL), marker = list(color = PROXY_COL, size = 10),
              error_y = ~list(array = hi - m, arrayminus = m - lo, thickness = 1, color = "#9a9a9a"),
              text = ~sprintf("%d training surveys: %.2f", n_train_countries, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Number of training surveys", dtick = 1),
               yaxis = list(title = "Accuracy, country left out", rangemode = "tozero"),
               margin = list(l = 10, r = 10, t = 10, b = 35)) |> config(displayModeBar = FALSE)
    })
    if (!is.null(go_to)) observeEvent(input$go_external, go_to("Tested in six more countries"))
  })
}
