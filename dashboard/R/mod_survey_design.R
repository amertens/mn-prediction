# =============================================================================
# Module: Plan a survey
# =============================================================================
# The anchor-and-rank design (AR-01): what a small national biomarker sample
# plus the transported ranking gives, against spending the same sample on a
# district or regional survey. A design result on these four surveys, not a
# pilot.

mod_survey_design_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330, title = "Design choices",
      selectInput(ns("rank_from"), "Ranking learned from", choices = NULL),
      selectInput(ns("set"), "Predictors in the ranking", choices = NULL),
      sliderInput(ns("fraction"), "Share of a full survey's sample", min = 0.05, max = 1, value = 0.05, step = 0.05),
      hr(),
      uiOutput(ns("at_fraction"))
    ),
    layout_columns(
      col_widths = c(12, 12),
      card(card_header("District error against sample size, four designs"),
           card_body(plotlyOutput(ns("mae"), height = "380px"),
                     methods_note(sprintf(paste("Median district error in percentage points over 22 country-outcome combinations. The national anchor",
                                                "plus ranking measures only a national prevalence from the given share of a full survey's respondents and",
                                                "orders districts by the transported ranking; the district and regional surveys spend the same sample on",
                                                "direct estimates. At 5 percent of the sample the anchored ranking gives %s points, a regional survey %s and",
                                                "a district survey %s; the district survey needs about %s of the full sample to match."),
                                          fmt_num(Q$ar_a1, 1), fmt_num(Q$ar_c, 1), fmt_num(Q$ar_b, 1), fmt_pct(Q$ar_b_match, 0))))),
      card(card_header("Burden captured by the worst-ranked fifth, same designs"),
           card_body(plotlyOutput(ns("capture"), height = "320px"),
                     methods_note("The ranking wins on level error only. A district survey of any size captures more of the burden,",
                                  " because burden sits in populous districts that surveys measure precisely. The design is a level",
                                  " calibration for a ranking, and it says where a survey should spend its clusters; it does not",
                                  " replace the survey.")))
    )
  )
}

mod_survey_design_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    AR <- EV$design_summary
    design_lab <- c(A1_anchor_rank = "National anchor + ranking", A2_region_anchor_rank = "Regional anchors + ranking",
                    B_district_survey = "District survey", C_regional_survey = "Regional survey")
    design_col <- c("National anchor + ranking" = PROXY_COL, "Regional anchors + ranking" = "#7fb3bd",
                    "District survey" = SURVEY_COL, "Regional survey" = "#e0a878")
    observe({
      req(AR)
      rf <- unique(AR$rank_from); st <- unique(AR$set)
      updateSelectInput(session, "rank_from", choices = setNames(rf, c(prev = "Prevalence", level = "Biomarker level")[rf]), selected = if ("prev" %in% rf) "prev" else rf[1])
      updateSelectInput(session, "set", choices = setNames(st, c(climate_soil = "Climate and soil (pre-registered)", full = "Full index")[st] %||% st),
                        selected = if ("climate_soil" %in% st) "climate_soil" else st[1])
      fr <- sort(unique(AR$fraction))
      updateSliderInput(session, "fraction", min = min(fr), max = max(fr), value = min(fr), step = min(diff(fr)))
    })
    dat <- reactive({
      req(AR, input$rank_from, input$set)
      d <- AR[AR$rank_from == input$rank_from & AR$set == input$set & AR$design %in% names(design_lab), ]
      d$Design <- design_lab[d$design]; d
    })
    output$mae <- renderPlotly({
      d <- dat(); validate(need(nrow(d) > 0, "No design rows for this choice."))
      plot_ly(d, x = ~fraction, y = ~mae, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%% of the sample: %.1f points", Design, 100 * fraction, mae), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the full survey sample", tickformat = ".0%"), yaxis = list(title = "Median district error (percentage points)", rangemode = "tozero"),
               shapes = list(list(type = "line", x0 = input$fraction, x1 = input$fraction, y0 = 0, y1 = max(d$mae, na.rm = TRUE), line = list(dash = "dot", color = "#888"))),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$capture <- renderPlotly({
      d <- dat(); validate(need(nrow(d) > 0 && "capture" %in% names(d), "No capture rows."))
      d <- d[is.finite(d$capture), ]
      plot_ly(d, x = ~fraction, y = ~capture, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%%: %.0f%% of burden in the picked fifth", Design, 100 * fraction, 100 * capture), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the full survey sample", tickformat = ".0%"), yaxis = list(title = "Burden in the worst-ranked fifth", tickformat = ".0%"),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$at_fraction <- renderUI({
      d <- dat(); req(nrow(d) > 0)
      fr <- d$fraction[which.min(abs(d$fraction - input$fraction))]
      s <- d[d$fraction == fr, ]
      tagList(
        h6(sprintf("At %s of a full survey's sample", fmt_pct(fr, 0))),
        tags$table(class = "table table-sm", style = "font-size:0.88em;",
                   tags$thead(tags$tr(tags$th("Design"), tags$th("Error, points"), tags$th("Ranking accuracy"), tags$th("Burden reached"))),
                   tags$tbody(lapply(seq_len(nrow(s)), function(i) tags$tr(tags$td(s$Design[i]), tags$td(fmt_num(s$mae[i], 1)), tags$td(fmt_num(s$spearman[i])), tags$td(fmt_pct(s$capture[i], 0)))))),
        p(style = "font-size:0.82em; color:#666;", "Fifty to seventy respondents is about 5 percent of these surveys' samples.")
      )
    })
  })
}
