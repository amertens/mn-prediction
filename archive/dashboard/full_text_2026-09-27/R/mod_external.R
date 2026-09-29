# =============================================================================
# Module: Tested in six more countries
# =============================================================================
# The climate-and-soil ranking compared with sub-national survey results that
# other teams deposited with the WHO, for six countries outside the training
# data on two continents. Labels this project did not collect or process.

mod_external_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = c(4, 8),
    div(
      card(card_header("What this test is"),
        card_body(
          p("All other accuracy figures on this dashboard come from the four surveys the model learned from. In",
            " September 2026 we also compared the climate-and-soil ranking with survey results that other teams had",
            " submitted to the WHO micronutrient database: ",
            strong("Zambia, Ethiopia, Sudan and Nigeria"), " in Africa, and ", strong("Pakistan and India"),
            ". None of these surveys was used to build the model, and we did not process their data."),
          p(sprintf(paste("In the four African countries the ranking matched the survey results at %s for biomarker levels,",
                          "above zero in %s of %s country-outcome tests (a random ranking stays below %s in 95%% of tries).",
                          "In Pakistan and India it reached %s for prevalence, above zero in %s of %s. Inside the four",
                          "training countries, the same climate-and-soil test at the regional level scores %s, so the",
                          "results in new countries are a little lower, as expected, but of the same size."),
                    f2(Q$xv_level), Q$xv_level_pos, Q$xv_level_cells, f2(Q$xv_level_null),
                    f2(Q$xv_off), Q$xv_off_pos, Q$xv_off_cells, f2(Q$cs_a1))),
          p(sprintf(paste("A globally available soil layer worked as well as the Africa-only one (%s against %s on the same",
                          "tests), so the same public data can be assembled for any country."),
                    f2(Q$xv_sg_level), f2(Q$xv_level))))),
      card(card_header("Limits of this check"),
        card_body(tags$ul(style = "font-size:0.88em; color:#555; padding-left:1.1em;",
          tags$li("These surveys report results for regions or provinces, not districts."),
          tags$li("Only rankings are compared. Cut-offs and laboratory methods differ between surveys, so prevalence",
                  " levels are not compared."),
          tags$li("Nigeria's results cover only six zones, too few to test on their own, so Nigeria counts only in the",
                  " combined result."),
          tags$li("Vitamin A did poorly in the African surveys but was the strongest outcome in South Asia. We do not",
                  " yet know why."),
          tags$li("This is an independent check. It does not replace the planned test on a new country's own survey",
                  " data, which needs the full survey records.")))))
    ,
    card(full_screen = TRUE,
      card_header("Each country and outcome, compared with a random ranking"),
      card_body(
        radioButtons(ns("target"), NULL, inline = TRUE,
                     choices = c("Biomarker level (African surveys)" = "level", "Prevalence (all six countries)" = "prev")),
        plotlyOutput(ns("cells"), height = "430px"),
        methods_note("Each dot is one country and outcome: the match (Spearman correlation) between the climate-and-soil",
                     " ranking and that country's survey results across its regions. Diamonds are country averages. The",
                     " dashed line is what a random ranking reaches in 95% of tries, when all outcomes from one country are",
                     " shuffled together. African countries use the African soil layer; Pakistan and India use a global",
                     " soil layer. Source: WHO Vitamin and Mineral Nutrition Information System.")),
      card_footer(tags$small(class = "text-muted",
        "The conference presentations and the manuscript use the same result tables."))
    )
  )
}

mod_external_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$cells <- renderPlotly({
      XC <- EV$xv_cells; validate(need(!is.null(XC), "External check tables not built yet."))
      tg <- input$target %||% "level"
      d <- XC[XC$arm == "domain_index" & XC$target == tg & is.finite(XC$spearman), ]
      if ("thin" %in% names(d)) d <- d[!d$thin, ]
      # African countries on the pre-registered (iSDA) soil block; South Asia exists only on the global block
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
