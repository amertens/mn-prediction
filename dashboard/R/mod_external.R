# =============================================================================
# Module: Tested in six more countries (concise version)
# =============================================================================

mod_external_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = c(4, 8),
    div(
      card(card_header("What this test is"),
        card_body(
          p("We compared the climate-and-soil ranking with survey results held by the WHO for Zambia, Ethiopia, Sudan,",
            " Nigeria, Pakistan and India. None of these surveys was used to build the model."),
          p(sprintf("The ranking matched at %s in the four African countries (biomarker levels) and %s in Pakistan and India (prevalence). A random ranking stays below %s.",
                    f2(Q$xv_level), f2(Q$xv_off), f2(Q$xv_level_null))))),
      card(card_header("Limits"),
        card_body(tags$ul(style = "font-size:0.88em; color:#555; padding-left:1.1em;",
          tags$li("Regional results, not districts."),
          tags$li("Rankings only, not prevalence levels."),
          tags$li("Nigeria has only six zones.")))))
    ,
    card(full_screen = TRUE,
      card_header("Each country and outcome, compared with a random ranking"),
      card_body(
        radioButtons(ns("target"), NULL, inline = TRUE,
                     choices = c("Biomarker level (African surveys)" = "level", "Prevalence (all six countries)" = "prev")),
        plotlyOutput(ns("cells"), height = "430px"),
        p(style = "font-size:0.85em; color:#666;", "Dots: one country and outcome. Diamonds: country averages. Dashed line: random ranking.")))
  )
}

mod_external_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$cells <- renderPlotly({
      XC <- EV$xv_cells; validate(need(!is.null(XC), "External check tables not built yet."))
      tg <- input$target %||% "level"
      d <- XC[XC$arm == "domain_index" & XC$target == tg & is.finite(XC$spearman), ]
      if ("thin" %in% names(d)) d <- d[!d$thin, ]
      d <- d[(d$arm_group == "africa" & d$soil == "isda") | (d$arm_group == "offcontinent" & d$soil == "sgrid"), ]
      validate(need(nrow(d) > 0, "No results for this outcome measure."))
      d$oc <- ifelse(d$outcome %in% names(outcome_short), outcome_short[d$outcome], d$outcome)
      mn <- d |> group_by(country) |> summarise(m = mean(spearman), .groups = "drop") |> arrange(m)
      d$country <- factor(d$country, levels = mn$country)
      null95 <- if (tg == "level") Q$xv_level_null else max(Q$xv_off_null, Q$xv_level_null, na.rm = TRUE)
      plot_ly() |>
        add_markers(data = d, x = ~spearman, y = ~country, color = ~oc,
                    marker = list(size = 9, opacity = 0.75),
                    text = ~sprintf("%s, %s<br>match %.2f over %d regions", country, oc, spearman, n_units),
                    hoverinfo = "text") |>
        add_markers(data = mn, x = ~m, y = ~country, marker = list(symbol = "diamond", size = 14, color = "#1a1a1a"),
                    text = ~sprintf("%s average %.2f", country, m), hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Match with the survey's ranking of regions (0 = random, 1 = perfect)", zeroline = TRUE),
               yaxis = list(title = ""),
               shapes = list(list(type = "line", x0 = null95, x1 = null95, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(dash = "dash", color = "#666"))),
               legend = list(orientation = "h", y = -0.18), margin = list(l = 10, r = 10, t = 10, b = 40)) |>
        config(displayModeBar = FALSE)
    })
  })
}
