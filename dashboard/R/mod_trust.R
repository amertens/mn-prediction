# =============================================================================
# Module: How well it works (concise version)
# =============================================================================

mod_trust_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Three tests", icon = bsicons::bs_icon("check2-square"),
      layout_columns(col_widths = c(3, 9),
        div(radioButtons(ns("estimand"), "Test", choices = c("A district hidden inside a surveyed country" = "infill",
                                                              "A whole region hidden" = "region",
                                                              "A whole country left out" = "country"), selected = "infill"),
            radioButtons(ns("target"), "Outcome measured as", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
            uiOutput(ns("test_text"))),
        plotlyOutput(ns("forest"), height = "420px")),
      p(style = "font-size:0.85em; color:#666;", "0 = random order, 1 = the survey's order. Dashed line: what a random ranking reaches.")
    ),
    nav_panel(
      title = "Best achievable score", icon = bsicons::bs_icon("arrow-bar-up"),
      layout_columns(col_widths = c(3, 9),
        div(radioButtons(ns("ceil_target"), "Outcome measured as", choices = c("Prevalence" = "prev", "Biomarker level" = "level"), selected = "prev"),
            uiOutput(ns("ceiling_text"))),
        plotlyOutput(ns("ceiling"), height = "380px"))
    ),
    nav_panel(
      title = "Each survey added", icon = bsicons::bs_icon("graph-up-arrow"),
      layout_columns(col_widths = c(3, 9),
        div(uiOutput(ns("curve_text"))),
        plotlyOutput(ns("curve"), height = "380px"))
    ),
    nav_panel(
      title = "Compared with a geostatistical model", icon = bsicons::bs_icon("geo"),
      layout_columns(col_widths = c(5, 7),
        div(uiOutput(ns("geostat_text"))),
        reactableOutput(ns("geostat")))
    )
  )
}

mod_trust_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$test_text <- renderUI({
      switch(input$estimand,
        infill = p(style = "font-size:0.88em; color:#555;",
                   sprintf("Inside a surveyed country the model scores %s, above the survey's regional averages (%s) and close to averaging neighbouring districts.",
                           fmt_num(Q$infill), fmt_num(Q$infill_jk))),
        region = p(style = "font-size:0.88em; color:#555;",
                   sprintf("With a whole region hidden, the model scores %s.", fmt_num(Q$region))),
        country = p(style = "font-size:0.88em; color:#555;",
                    sprintf("With a whole country left out, the model scores %s (%s with climate and soil only); a random ranking stays below %s.",
                            fmt_num(Q$tr), fmt_num(Q$cs), fmt_num(Q$null_d))))
    })

    output$forest <- renderPlotly({
      BC <- EV$benchmarks_cells; validate(need(!is.null(BC), "Test results not built."))
      d <- BC[BC$estimand == input$estimand & BC$target == input$target & BC$arm %in% names(arm_label), ]
      if (input$estimand == "country" && !is.null(EV$nested_domains)) {
        cs <- EV$nested_domains[EV$nested_domains$arm == "fixed_cs" & EV$nested_domains$target == input$target, ]
        if (nrow(cs)) d <- bind_rows(d, data.frame(arm = "climate_soil", spearman = cs$spearman))
      }
      validate(need(nrow(d) > 0, "No results for this test."))
      s <- cell_ci(d, "spearman", "arm")
      lab <- c(arm_label, climate_soil = "This model, climate and soil only")
      s$Model <- lab[s$arm]; s <- s[order(-s$est), ]
      null <- if (input$estimand == "country") Q$null_d else 0
      forest_plotly(s, "Model", null = null)
    })

    output$ceiling_text <- renderUI({
      p(style = "font-size:0.88em; color:#555;",
        sprintf(paste("Noise in the survey's own district figures limits any model. The best achievable score averages %s for",
                      "prevalence and %s for biomarker levels; the model reaches %s and %s."),
                fmt_num(Q$ceiling_prev), fmt_num(Q$ceiling_level), fmt_num(Q$infill_prev), fmt_num(Q$infill)))
    })

    output$ceiling <- renderPlotly({
      VC <- EV$ceiling; validate(need(!is.null(VC), "Table not built."))
      v <- VC[VC$rung == "admin2" & VC$target == input$ceil_target, ]
      d <- v |> group_by(country) |> summarise(`Best achievable score` = mean(ceiling_vc, na.rm = TRUE),
                                               `This model` = mean(achieved_spearman, na.rm = TRUE), .groups = "drop") |>
        pivot_longer(-country, names_to = "quantity", values_to = "v")
      d <- d[is.finite(d$v), ]
      cols <- c("Best achievable score" = PROXY_COL, "This model" = SURVEY_COL)
      plot_ly(d, x = ~v, y = ~country, color = ~quantity, colors = cols, type = "scatter", mode = "markers", marker = list(size = 13),
              text = ~sprintf("%s<br>%s: %.2f", country, quantity, v), hoverinfo = "text") |>
        layout(xaxis = list(title = "Ranking accuracy (average over outcomes)", range = c(0, 1)), yaxis = list(title = ""),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$curve_text <- renderUI({
      p(style = "font-size:0.88em; color:#555;",
        sprintf("Accuracy in a country left out rose from %s with one training survey to %s with three.", fmt_num(Q$lc1), fmt_num(Q$lc_max)))
    })

    output$curve <- renderPlotly({
      TC <- EV$training_curve; validate(need(!is.null(TC), "Table not built."))
      a <- TC[TC$arm == "domain_index" & TC$target == "level", c("n_train_countries", "spearman")]; a$set <- "All data layers"
      b <- EV$training_curve_cs
      if (!is.null(b)) { b <- b[b$set == "climate_soil" & b$target == "level", c("n_train_countries", "spearman")]; b$set <- "Climate and soil only"; a <- rbind(a, b) }
      d <- a |> group_by(set, n_train_countries) |> summarise(m = mean(spearman, na.rm = TRUE), lo = quantile(spearman, 0.25, na.rm = TRUE),
                                                              hi = quantile(spearman, 0.75, na.rm = TRUE), .groups = "drop")
      plot_ly(d, x = ~n_train_countries, y = ~m, color = ~set, colors = c("All data layers" = PROXY_COL, "Climate and soil only" = "#8c510a"),
              type = "scatter", mode = "lines+markers", error_y = ~list(array = hi - m, arrayminus = m - lo, thickness = 1),
              text = ~sprintf("%s<br>%d training surveys: %.2f", set, n_train_countries, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Number of training surveys", dtick = 1), yaxis = list(title = "Ranking accuracy in the country left out", rangemode = "tozero"),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$geostat_text <- renderUI({
      p(style = "font-size:0.9em;",
        sprintf(paste("A geostatistical model estimates percentages slightly better (%s against %s points of error) but ranks",
                      "districts worse (%s against %s), and it cannot be used for a country without a survey."),
                fmt_num(Q$mbg_err, 1), fmt_num(Q$mbg_err_index, 1), fmt_num(Q$mbg_rank), fmt_num(Q$mbg_rank_index)))
    })

    output$geostat <- renderReactable({
      MB <- EV$geostat_cells; validate(need(!is.null(MB), "Comparison table not built."))
      lab <- c(domain_index = "This model", mbg = "Geostatistical, our data layers", mbg_dhs = "Geostatistical, DHS list of layers",
               mbg_gp = "Geostatistical, location only", spatial = "Neighbouring districts' average", region_mean_jk = "Survey's regional averages")
      m <- MB[MB$arm %in% names(lab) & MB$estimand %in% c("infill", "region"), ] |>
        group_by(estimand, target, arm) |> summarise(rho = mean(rho_agg, na.rm = TRUE), err = mean(wmae_agg, na.rm = TRUE), n = n(), .groups = "drop")
      t <- m |> mutate(Model = lab[arm], Test = c(infill = "District hidden", region = "Region hidden")[estimand],
                       Target = c(level = "Biomarker level", prev = "Prevalence")[target]) |>
        transmute(Model, Test, `Measured as` = Target, `Ranking accuracy` = round(rho, 2),
                  `Error (points)` = ifelse(target == "prev", round(err, 1), NA)) |>
        arrange(Test, `Measured as`, desc(`Ranking accuracy`))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, groupBy = c("Test", "Measured as"))
    })
  })
}
