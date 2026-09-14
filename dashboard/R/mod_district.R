# =============================================================================
# Module: District profiles
# =============================================================================
# One district across every outcome its country measured: rank, how sure,
# planning prevalence, and the survey's own figure, then the predictors behind
# the score for a chosen outcome.

mod_district_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 320, title = "Choose a district",
      selectInput(ns("country"), "Country", choices = country_choices, selected = "ghana"),
      selectizeInput(ns("district"), "District", choices = NULL),
      selectInput(ns("outcome"), "Outcome for the drivers", choices = outcome_choices, selected = "child_vitA"),
      hr(),
      uiOutput(ns("summary"))
    ),
    layout_columns(
      col_widths = c(12, 12),
      card(card_header("Where this district ranks, by outcome"),
           card_body(plotlyOutput(ns("profile"), height = "320px"),
                     reactableOutput(ns("table")),
                     methods_note("Priority score 100 means ranked worst in the country. The chance of being in the worst fifth",
                                  " exists for surveyed districts only. The planning prevalence is the ranking anchored to the",
                                  " national survey; the survey estimate carries a 95 percent range from its own effective sample size."))),
      card(card_header("What drives the score for the chosen outcome"),
           card_body(plotlyOutput(ns("drivers"), height = "380px"),
                     methods_note("The model's score is a weighted sum of the district's predictors, so each bar is that",
                                  " predictor's exact contribution to the district's distance from the country's surveyed mean,",
                                  " on the model's own scale. Red pushes toward more deficiency, blue toward less. These are",
                                  " markers of where deficiency is, not causes.")))
    )
  )
}

mod_district_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$country, {
      d <- idx_districts[idx_districts$country_key == input$country, ]
      d <- d[!duplicated(paste(d$Admin1, d$Admin2)), ]
      d <- d[order(d$Admin1, d$Admin2), ]
      ch <- setNames(paste(d$Admin1, d$Admin2, sep = "|"), sprintf("%s (%s)", d$Admin2, d$Admin1))
      updateSelectizeInput(session, "district", choices = ch, selected = ch[[1]], server = TRUE)
      oc <- outcomes_for(input$country)
      updateSelectInput(session, "outcome", choices = oc, selected = if ((input$outcome %||% "") %in% oc) input$outcome else oc[[1]])
    })

    rows <- reactive({
      req(input$country, input$district)
      k <- strsplit(input$district, "|", fixed = TRUE)[[1]]
      d <- idx_districts[idx_districts$country_key == input$country & idx_districts$Admin1 == k[1] & idx_districts$Admin2 == k[2], ]
      d <- d[order(match(d$outcome, names(meta$outcome_labels))), ]
      d$label <- outcome_short[d$outcome]
      d
    })

    output$summary <- renderUI({
      d <- rows(); req(nrow(d) > 0)
      tagList(
        h5(d$Admin2[1], style = "margin-top:0;"), p(em(d$Admin1[1])),
        p(if (isTRUE(d$surveyed[1])) sprintf("Surveyed: %s respondents in %s clusters.", fmt_count(d$n_resp[1]), d$n_clusters[1])
          else "No survey clusters in this district; every figure is the model's."),
        p(sprintf("In the worst fifth for %d of %d outcomes.", sum(d$rank_worst <= ceiling(d$n_districts / 5)), nrow(d))),
        if (is.finite(d$population[1])) p(sprintf("Children 6 to 59 months: %s; women 15 to 49: %s.",
                                                  fmt_count(g1(d$population[startsWith(d$outcome, "child_")])),
                                                  fmt_count(g1(d$population[startsWith(d$outcome, "women_")]))), style = "font-size:0.9em; color:#555;")
      )
    })

    output$profile <- renderPlotly({
      d <- rows(); req(nrow(d) > 0)
      d$label <- factor(d$label, levels = rev(d$label))
      plot_ly(d) |>
        add_segments(x = 0, xend = 100, y = ~label, yend = ~label, line = list(color = "#e6e6e6", width = 6), showlegend = FALSE, hoverinfo = "none") |>
        add_markers(x = ~priority, y = ~label, marker = list(color = PROXY_COL, size = 13),
                    text = ~sprintf("%s<br>priority %.0f, rank %d of %d<br>planning prevalence %s%s", label, priority, rank_worst, n_districts,
                                    fmt_pct(prev_anchored), ifelse(is.finite(p_worst_fifth), sprintf("<br>chance of worst fifth %s", fmt_pct(p_worst_fifth, 0)), "")),
                    hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Priority score (100 = ranked worst in the country; dotted line = worst fifth)", range = c(-2, 102)),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50),
               shapes = list(list(type = "line", x0 = 80, x1 = 80, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(color = "#b2182b", dash = "dot")))) |> config(displayModeBar = FALSE)
    })

    output$table <- renderReactable({
      d <- rows(); req(nrow(d) > 0)
      t <- data.frame(Outcome = d$label, Rank = sprintf("%d of %d", d$rank_worst, d$n_districts),
                      `Chance of worst fifth` = ifelse(is.finite(d$p_worst_fifth), fmt_pct(d$p_worst_fifth, 0), "—"),
                      `Planning prevalence` = fmt_pct(d$prev_anchored),
                      `Survey estimate` = ifelse(is.finite(d$survey_prev), sprintf("%s (%s to %s)", fmt_pct(d$survey_prev), fmt_pct(d$survey_lo, 0), fmt_pct(d$survey_hi, 0)), "not surveyed"),
                      `WHO class` = d$who_class, check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 8)
    })

    output$drivers <- renderPlotly({
      d <- rows(); req(nrow(d) > 0, input$outcome)
      dec <- decompose_district(input$country, input$outcome, d$Admin1[1], d$Admin2[1])
      validate(need(!is.null(dec) && nrow(dec) > 0, "No fit for this country and outcome."))
      top <- head(dec, 12); top$label <- unique_labels(top$label, top$column); top$label <- factor(top$label, levels = rev(top$label))
      plot_ly(top, x = ~contribution, y = ~label, type = "bar", orientation = "h",
              marker = list(color = ifelse(top$contribution > 0, "#b2182b", "#2166ac")),
              text = ~sprintf("%s<br>%s<br>contribution %+.2f", column, source, contribution), hoverinfo = "text") |>
        layout(xaxis = list(title = sprintf("Contribution to the score, %s (model scale; total %+.2f)", outcome_short[[input$outcome]], attr(dec, "total"))),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })
  })
}
