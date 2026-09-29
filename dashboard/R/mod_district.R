# =============================================================================
# Module: District profiles (concise version)
# =============================================================================

mod_district_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 320, title = "Choose a district",
      selectInput(ns("country"), "Country", choices = country_choices, selected = "ghana"),
      selectizeInput(ns("district"), "District", choices = NULL),
      selectInput(ns("outcome"), "Outcome for the breakdown", choices = outcome_choices, selected = "child_vitA"),
      hr(),
      uiOutput(ns("summary"))
    ),
    layout_columns(
      col_widths = c(12, 12),
      card(card_header("Where this district ranks, by outcome"),
           card_body(plotlyOutput(ns("profile"), height = "320px"),
                     reactableOutput(ns("table")),
                     p(style = "font-size:0.85em; color:#666;", "100 = ranked worst in the country. Ranges are explained in Technical notes."))),
      card(card_header("What moves this district's score"),
           card_body(plotlyOutput(ns("drivers"), height = "380px"),
                     p(style = "font-size:0.85em; color:#666;",
                       "Each bar is one data layer's contribution. Red pushes toward more deficiency, blue toward less.",
                       " They show where deficiency is likely, not what causes it.")))
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
      d$rho_train <- idx_national$rho_train[match(paste(d$country_key, d$outcome), paste(idx_national$country_key, idx_national$outcome))]
      d$level_skill <- level_skill(d$rho_train)$band
      d$rank_lo <- d$rank_hi <- NA_real_
      if (length(UE) && !is.null(UE$cells)) {
        u <- UE$cells[UE$cells$country_key == input$country & UE$cells$Admin1 == k[1] & UE$cells$Admin2 == k[2], ]
        j <- match(d$outcome, u$outcome)
        for (cc in c("rank_lo", "rank_hi")) d[[cc]] <- u[[cc]][j]
      }
      d
    })

    output$summary <- renderUI({
      d <- rows(); req(nrow(d) > 0)
      tagList(
        h5(d$Admin2[1], style = "margin-top:0;"), p(em(d$Admin1[1])),
        p(if (isTRUE(d$surveyed[1])) sprintf("Surveyed: %s people in %s clusters.", fmt_count(d$n_resp[1]), d$n_clusters[1])
          else "Not surveyed: all figures come from the model."),
        p(sprintf("In the worst fifth for %d of %d outcomes.", sum(d$rank_worst <= ceiling(d$n_districts / 5)), nrow(d)))
      )
    })

    output$profile <- renderPlotly({
      d <- rows(); req(nrow(d) > 0)
      d$label <- factor(d$label, levels = rev(d$label))
      plot_ly(d) |>
        add_segments(x = 0, xend = 100, y = ~label, yend = ~label, line = list(color = "#e6e6e6", width = 6), showlegend = FALSE, hoverinfo = "none") |>
        add_markers(x = ~priority, y = ~label, marker = list(color = PROXY_COL, size = 13),
                    text = ~sprintf("%s<br>priority %.0f, rank %d of %d<br>estimated prevalence %s", label, priority, rank_worst, n_districts, fmt_pct(prev_anchored)),
                    hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Priority score (100 = ranked worst; dotted line = worst fifth)", range = c(-2, 102)),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50),
               shapes = list(list(type = "line", x0 = 80, x1 = 80, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(color = "#b2182b", dash = "dot")))) |> config(displayModeBar = FALSE)
    })

    output$table <- renderReactable({
      d <- rows(); req(nrow(d) > 0)
      t <- data.frame(Outcome = d$label, Rank = sprintf("%d of %d", d$rank_worst, d$n_districts),
                      `Rank range` = ifelse(is.finite(d$rank_lo), sprintf("%d to %d", round(d$rank_lo), round(d$rank_hi)), "—"),
                      `Estimated prevalence (checked range)` = ifelse(is.finite(d$prev_cal_lo), sprintf("%s (%s to %s)", fmt_pct(d$prev_anchored), fmt_pct(d$prev_cal_lo), fmt_pct(d$prev_cal_hi)), fmt_pct(d$prev_anchored)),
                      `Percentages` = skill_word[d$level_skill],
                      `Survey estimate` = ifelse(is.finite(d$survey_prev), fmt_pct(d$survey_prev), "not surveyed"),
                      check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 11, rownames = FALSE)
    })

    output$drivers <- renderPlotly({
      d <- rows(); req(nrow(d) > 0, input$outcome)
      dec <- decompose_district(input$country, input$outcome, d$Admin1[1], d$Admin2[1])
      validate(need(!is.null(dec) && nrow(dec) > 0, "No model for this country and outcome."))
      top <- head(dec, 12); top$label <- unique_labels(top$label, top$column); top$label <- factor(top$label, levels = rev(top$label))
      plot_ly(top, x = ~contribution, y = ~label, type = "bar", orientation = "h",
              marker = list(color = ifelse(top$contribution > 0, "#b2182b", "#2166ac")),
              text = ~sprintf("%s<br>%s<br>contribution %+.2f", column, source, contribution), hoverinfo = "text") |>
        layout(xaxis = list(title = sprintf("Contribution to the score, %s", outcome_short[[input$outcome]])),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })
  })
}
