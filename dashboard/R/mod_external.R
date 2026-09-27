# =============================================================================
# Module: Tested in six more countries
# =============================================================================
# XV-01/02 (22 Sep 2026): the transported climate + soil ranking scored against
# WHO-deposited sub-national survey results in six countries outside the
# training panel, on two continents — labels this project never collected,
# harmonised or weighted. The strongest honest evidence the approach travels.

mod_external_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = c(4, 8),
    div(
      card(card_header("The one test we could not rig"),
        card_body(
          p("Every accuracy number elsewhere on this dashboard is computed inside the four training surveys.",
            " In September 2026 the same climate-and-soil ranking was scored against sub-national deficiency",
            " results that other teams deposited with the WHO micronutrients database: ",
            strong("Zambia, Ethiopia, Sudan and Nigeria"), " in Africa, and ", strong("Pakistan and India"),
            " off the continent. None was used in training; none was touched by our processing."),
          p(sprintf(paste("Across the African deposits the ranking matches at %s on the biomarker level, positive in",
                          "%s of %s country-outcome tests (chance stays under %s). Off the continent it reaches %s,",
                          "positive in %s of %s. The internal held-out figure at the same regional grain is lower,",
                          "not higher, so these countries did not get a diluted version of the claim."),
                    f2(Q$xv_level), Q$xv_level_pos, Q$xv_level_cells, f2(Q$xv_level_null),
                    f2(Q$xv_off), Q$xv_off_pos, Q$xv_off_cells)),
          p(strong("Why it matters for a new country: "),
            "the recipe needs no Africa-only data. A global soil layer performs as well as the African one",
            sprintf(" (%s against %s on the same cells), so the same public recipe can be built anywhere.",
                    f2(Q$xv_sg_level), f2(Q$xv_level))))),
      card(card_header("Read the fine print"),
        card_body(tags$ul(style = "font-size:0.88em; color:#555; padding-left:1.1em;",
          tags$li("These deposits resolve regions or provinces (first administrative level), not districts."),
          tags$li("Rankings only: deposited cut-offs and assays differ, so levels are never compared."),
          tags$li("Nigeria contributes six zones — too few to test alone; it counts only in the pool."),
          tags$li("Vitamin A fails in the African deposits and is the strongest outcome in South Asia;",
                  " recorded, not yet explained."),
          tags$li("This is an independent check consistent with the pre-registered new-country test,",
                  " not a substitute for it: that test needs a survey's own microdata.")))))
    ,
    card(full_screen = TRUE,
      card_header("Every held-out country and outcome, against what chance would give"),
      card_body(
        radioButtons(ns("target"), NULL, inline = TRUE,
                     choices = c("Biomarker level (Africa deposits)" = "level", "Prevalence (all six)" = "prev")),
        plotlyOutput(ns("cells"), height = "430px"),
        methods_note("Each dot is one country and outcome: the rank agreement (Spearman) between the transported",
                     " climate-and-soil index and that country's deposited sub-national results, over its regions.",
                     " Diamonds are country means. The dashed line is the 95th percentile of a country-block",
                     " permutation null — what 'no information' reaches. African countries are scored with the",
                     " African soil layer, Pakistan and India with the global one (which ties it on the same cells).",
                     " Sources: WHO VMNIS deposits; details in the project's XV-01/02 note.")),
      card_footer(tags$small(class = "text-muted",
        "The Micronutrient Forum decks and the manuscript quote these same tables (results/tables/external_validation)."))
    )
  )
}

mod_external_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$cells <- renderPlotly({
      XC <- EV$xv_cells; validate(need(!is.null(XC), "External-validation tables not built yet."))
      tg <- input$target %||% "level"
      d <- XC[XC$arm == "domain_index" & XC$target == tg & is.finite(XC$spearman), ]
      if ("thin" %in% names(d)) d <- d[!d$thin, ]
      # African countries on the pre-registered (iSDA) soil block; South Asia exists only on the global block
      d <- d[(d$arm_group == "africa" & d$soil == "isda") | (d$arm_group == "offcontinent" & d$soil == "sgrid"), ]
      validate(need(nrow(d) > 0, "No cells for this target."))
      d$oc <- ifelse(d$outcome %in% names(outcome_short), outcome_short[d$outcome], d$outcome)
      mn <- d |> group_by(country) |> summarise(m = mean(spearman), .groups = "drop") |> arrange(m)
      d$country <- factor(d$country, levels = mn$country)
      null95 <- if (tg == "level") Q$xv_level_null else max(Q$xv_off_null, Q$xv_level_null, na.rm = TRUE)
      plot_ly() |>
        add_markers(data = d, x = ~spearman, y = ~country, color = ~oc,
                    marker = list(size = 9, opacity = 0.75),
                    text = ~sprintf("%s, %s<br>rank agreement %.2f over %d units", country, oc, spearman, n_units),
                    hoverinfo = "text") |>
        add_markers(data = mn, x = ~m, y = ~country, marker = list(symbol = "diamond", size = 14, color = "#1a1a1a"),
                    text = ~sprintf("%s mean %.2f", country, m), hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Rank agreement with the deposited survey (0 = chance, 1 = perfect)", zeroline = TRUE),
               yaxis = list(title = ""),
               shapes = list(list(type = "line", x0 = null95, x1 = null95, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(dash = "dash", color = "#666"))),
               legend = list(orientation = "h", y = -0.18), margin = list(l = 10, r = 10, t = 10, b = 40)) |>
        config(displayModeBar = FALSE)
    })
  })
}
