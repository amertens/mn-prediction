# =============================================================================
# Module: How well it works
# =============================================================================
# The protocol's three tests, the ceiling the survey's own noise sets, what
# each added survey buys, and the geostatistical comparator. Every figure is
# drawn from the committed protocol tables.

mod_trust_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Three tests", icon = bsicons::bs_icon("check2-square"),
      layout_columns(col_widths = c(3, 9),
        div(radioButtons(ns("estimand"), "Test", choices = c("A district hidden inside a surveyed country" = "infill",
                                                              "A whole region hidden" = "region",
                                                              "A whole country hidden" = "country"), selected = "infill"),
            radioButtons(ns("target"), "Target", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
            uiOutput(ns("test_text"))),
        plotlyOutput(ns("forest"), height = "420px")),
      methods_note("Ranking accuracy is the Spearman correlation between the model's order of districts and the survey's,",
                   " mean over country-outcome combinations with a 95 percent interval across them; random splits are",
                   " repeated ten times. Each comparator saw exactly what the model saw: the survey's regional average is",
                   " computed without the scored district, and the smoothers run on the same folds. The dashed line is",
                   " the 95th percentile of a permutation null. Differences under 0.03 are ties.")
    ),
    nav_panel(
      title = "The ceiling", icon = bsicons::bs_icon("arrow-bar-up"),
      layout_columns(col_widths = c(3, 9),
        div(radioButtons(ns("ceil_target"), "Target", choices = c("Prevalence" = "prev", "Biomarker level" = "level"), selected = "prev"),
            uiOutput(ns("ceiling_text"))),
        plotlyOutput(ns("ceiling"), height = "380px")),
      methods_note("Most districts hold one or two survey clusters, so the survey's own district value is an estimate",
                   " with noise of its own. A variance model with respondents in clusters in districts in regions gives",
                   " the correlation a perfect predictor of the true district value could reach. The split-half version",
                   " counts the cluster as geography and is too generous; the honest ceiling removes it.")
    ),
    nav_panel(
      title = "Each survey added", icon = bsicons::bs_icon("graph-up-arrow"),
      layout_columns(col_widths = c(3, 9),
        div(uiOutput(ns("curve_text"))),
        plotlyOutput(ns("curve"), height = "380px")),
      methods_note("Every subset of training countries scored on each held-out country, biomarker level. Bars are the",
                   " interquartile range over fits. The climate-and-soil index was chosen on these same four countries,",
                   " so it enters the next country as a pre-registered prediction rather than a result.")
    ),
    nav_panel(
      title = "A better number, a worse ranking", icon = bsicons::bs_icon("geo"),
      layout_columns(col_widths = c(5, 7),
        div(uiOutput(ns("geostat_text"))),
        reactableOutput(ns("geostat"))),
      methods_note("The DHS Program's model-based geostatistics, a spatial field plus covariates fitted at survey cluster",
                   " locations, scored on the same district folds as the index, once with our covariates and once with",
                   " the DHS Program's own ten-layer list. It needs survey clusters inside the country, so it has no",
                   " country-transport test.")
    ),
    nav_panel(
      title = "What else was tried", icon = bsicons::bs_icon("list-check"),
      uiOutput(ns("tried"))
    )
  )
}

mod_trust_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$test_text <- renderUI({
      switch(input$estimand,
        infill = p(style = "font-size:0.88em; color:#555;",
                   sprintf("Inside a surveyed country the proxy index (%s), the neighbour smoother with proxies and the smoother alone are a three-way tie, all above the survey's own regional averages (%s). Where a survey exists, geography does most of the work.",
                           fmt_num(Q$infill), fmt_num(Q$infill_jk))),
        region = p(style = "font-size:0.88em; color:#555;",
                   sprintf("With a whole region hidden the smoother must extrapolate past the edge of its data. The index scores %s on the level and stays close to the smoother; the survey's own regional average cannot be computed for a region it did not sample.", fmt_num(Q$region))),
        country = p(style = "font-size:0.88em; color:#555;",
                    sprintf("With a whole country hidden only the proxy methods produce anything: %s with everything, %s from climate and soil alone, against %s for chance. The ranking crosses borders; the prevalence level does not, so a new country needs one national survey number to anchor it.",
                            fmt_num(Q$tr), fmt_num(Q$cs), fmt_num(Q$null_d))))
    })

    output$forest <- renderPlotly({
      BC <- EV$benchmarks_cells; validate(need(!is.null(BC), "Benchmark cells not built."))
      d <- BC[BC$estimand == input$estimand & BC$target == input$target & BC$arm %in% names(arm_label), ]
      if (input$estimand == "country" && !is.null(EV$nested_domains)) {
        cs <- EV$nested_domains[EV$nested_domains$arm == "fixed_cs" & EV$nested_domains$target == input$target, ]
        if (nrow(cs)) d <- bind_rows(d, data.frame(arm = "climate_soil", spearman = cs$spearman))
      }
      validate(need(nrow(d) > 0, "No cells for this test."))
      s <- cell_ci(d, "spearman", "arm")
      lab <- c(arm_label, climate_soil = "Climate and soil index")
      s$Model <- lab[s$arm]; s <- s[order(-s$est), ]
      null <- if (input$estimand == "country") Q$null_d else 0
      forest_plotly(s, "Model", null = null)
    })

    output$ceiling_text <- renderUI({
      p(style = "font-size:0.88em; color:#555;",
        sprintf(paste("On prevalence the honest ceiling averages %s and the model reaches %s; on the biomarker level %s and %s.",
                      "The model is at its ceiling in %s of %s cells. More clusters per district in the next survey would raise",
                      "the ceiling itself and close as much of the remaining third as better predictors would."),
                fmt_num(Q$ceiling_prev), fmt_num(Q$infill_prev), fmt_num(Q$ceiling_level), fmt_num(Q$infill), Q$vc_at, Q$vc_n))
    })

    output$ceiling <- renderPlotly({
      VC <- EV$ceiling; validate(need(!is.null(VC), "Ceiling table not built."))
      v <- VC[VC$rung == "admin2" & VC$target == input$ceil_target, ]
      d <- v |> group_by(country) |> summarise(`Split-half ceiling` = mean(ceiling_within_vc, na.rm = TRUE),
                                               `Honest ceiling` = mean(ceiling_vc, na.rm = TRUE),
                                               `Achieved, district hidden` = mean(achieved_spearman, na.rm = TRUE), .groups = "drop") |>
        pivot_longer(-country, names_to = "quantity", values_to = "v")
      d <- d[is.finite(d$v), ]
      cols <- c("Split-half ceiling" = "#bdbdbd", "Honest ceiling" = PROXY_COL, "Achieved, district hidden" = SURVEY_COL)
      plot_ly(d, x = ~v, y = ~country, color = ~quantity, colors = cols, type = "scatter", mode = "markers", marker = list(size = 13),
              text = ~sprintf("%s<br>%s: %.2f", country, quantity, v), hoverinfo = "text") |>
        layout(xaxis = list(title = "Ranking accuracy (mean over outcomes)", range = c(0, 1)), yaxis = list(title = ""),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$curve_text <- renderUI({
      p(style = "font-size:0.88em; color:#555;",
        sprintf(paste("Each training country added buys about %s of ranking accuracy in a country never seen (%s with one training country,",
                      "%s with three) and the curve has not flattened. A survey a country runs for itself improves every other",
                      "country's map. The climate-and-soil index starts higher and is what a fifth country will be scored on."),
                fmt_num(Q$lc_step), fmt_num(Q$lc1), fmt_num(Q$lc_max)))
    })

    output$curve <- renderPlotly({
      TC <- EV$training_curve; validate(need(!is.null(TC), "Training curve not built."))
      a <- TC[TC$arm == "domain_index" & TC$target == "level", c("n_train_countries", "spearman")]; a$set <- "Full index"
      b <- EV$training_curve_cs
      if (!is.null(b)) { b <- b[b$set == "climate_soil" & b$target == "level", c("n_train_countries", "spearman")]; b$set <- "Climate and soil index"; a <- rbind(a, b) }
      d <- a |> group_by(set, n_train_countries) |> summarise(m = mean(spearman, na.rm = TRUE), lo = quantile(spearman, 0.25, na.rm = TRUE),
                                                              hi = quantile(spearman, 0.75, na.rm = TRUE), .groups = "drop")
      plot_ly(d, x = ~n_train_countries, y = ~m, color = ~set, colors = c("Full index" = PROXY_COL, "Climate and soil index" = "#8c510a"),
              type = "scatter", mode = "lines+markers", error_y = ~list(array = hi - m, arrayminus = m - lo, thickness = 1),
              text = ~sprintf("%s<br>%d training countries: %.2f (IQR %.2f to %.2f)", set, n_train_countries, m, lo, hi), hoverinfo = "text") |>
        layout(xaxis = list(title = "Countries used for training", dtick = 1), yaxis = list(title = "Ranking accuracy in the held-out country", rangemode = "tozero"),
               shapes = list(list(type = "rect", x0 = 0.5, x1 = 3.5, y0 = 0, y1 = Q$null_d, fillcolor = "#e0e0e0", line = list(width = 0), layer = "below")),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$geostat_text <- renderUI({
      tagList(
        p(style = "font-size:0.9em;",
          sprintf(paste("The geostatistical model ranks districts worse than the index (%s against %s inside a surveyed country)",
                        "and gets the absolute number better (%s against %s percentage points of error), because it estimates a",
                        "level and shrinks toward it."), fmt_num(Q$mbg_rank), fmt_num(Q$mbg_rank_index), fmt_num(Q$mbg_err, 1), fmt_num(Q$mbg_err_index, 1))),
        p(style = "font-size:0.9em;", "So if a country needs a prevalence figure to publish for a surveyed district, that is the better tool,",
          " and it comes with an interval. If it needs to know which districts to look at first, or anything about a",
          " country with no survey, the ranking is."))
    })

    output$geostat <- renderReactable({
      MB <- EV$geostat_cells; validate(need(!is.null(MB), "Comparator table not built."))
      lab <- c(domain_index = "Proxy index", mbg = "Geostatistical, our covariates", mbg_dhs = "Geostatistical, DHS covariate list",
               mbg_gp = "Geostatistical, spatial field only", spatial = "Neighbour smoother", region_mean_jk = "Survey's regional average")
      m <- MB[MB$arm %in% names(lab) & MB$estimand %in% c("infill", "region"), ] |>
        group_by(estimand, target, arm) |> summarise(rho = mean(rho_agg, na.rm = TRUE), err = mean(wmae_agg, na.rm = TRUE), n = n(), .groups = "drop")
      t <- m |> mutate(Model = lab[arm], Test = c(infill = "District hidden", region = "Region hidden")[estimand],
                       Target = c(level = "level", prev = "prevalence")[target]) |>
        transmute(Model, Test, Target, `Ranking accuracy` = round(rho, 2), `Error, points` = ifelse(target == "prev", round(err, 1), NA), Cells = n) |>
        arrange(Test, Target, desc(`Ranking accuracy`))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, groupBy = c("Test", "Target"))
    })

    output$tried <- renderUI({
      WS <- EV$weight_sources
      rows <- if (!is.null(WS)) {
        lab <- c(domain_index = "The index (pooled rank correlation weights)", ridge_min = "Ridge regression on the same components",
                 domain_enet = "Elastic net on the components", lasso_min = "Lasso on the components",
                 index_decor = "Decorrelated weights", index_soft1 = "Soft-thresholded weights", sparse20 = "Twenty-predictor equal-weight composite")
        w <- WS[WS$arm %in% names(lab) & WS$target == "level", ] |> select(arm, estimand, mean_spearman) |>
          pivot_wider(names_from = estimand, values_from = mean_spearman)
        w$Method <- lab[w$arm]; w <- w[order(match(w$arm, names(lab))), ]
        tags$table(class = "table table-sm", style = "font-size:0.9em;",
                   tags$thead(tags$tr(tags$th("Weighting of the same components"), tags$th("District hidden"), tags$th("Region hidden"), tags$th("Country hidden"))),
                   tags$tbody(lapply(seq_len(nrow(w)), function(i) tags$tr(tags$td(w$Method[i]), tags$td(fmt_num(w$infill[i])), tags$td(fmt_num(w$region[i])), tags$td(fmt_num(w$country[i]))))))
      } else NULL
      tagList(
        p(class = "lead", "Nothing more complicated does better, and the reasons are sample size and survey noise."),
        tags$ul(
          tags$li(sprintf("Machine-learning ensembles (SuperLearner, twelve to sixteen learners) tie the untuned index at best and trail it under squared-error tuning. With 14 to 87 districts per country, every method that learns its own settings loses. Person-level prediction of who is deficient reaches an AUC of about %s, a coin toss, so it is not offered here.", fmt_num(Q$il_auc))),
          tags$li("Ridge, lasso and elastic-net weightings of the same domain components lose inside a country; the index is the infinite-penalty limit of that family. Two variants gain about 0.03 across borders and are pre-registered for the fifth country rather than adopted."),
          tags$li("Fitting at the survey cluster instead of the district does not beat the district fit on any test."),
          tags$li("Rebuilt on 16 September 2026 on 570 layers with the DHS aggregates held out: nothing inside a country changed, the district ranking for a country with no survey rose from 0.26 to 0.31, and restoring the soil layers that a review had dropped brought the climate-and-soil result back to 0.38 across districts and 0.45 across regions."),
          tags$li("Adding data blocks (livestock, water and coast distance, helminths, re-extracted IHME surfaces, fieldwork-month prices and temperatures) moves the index by less than 0.01 each. More kinds of data no longer help; the two remotely sensed domains carry the model.")
        ),
        rows,
        methods_note("An earlier version of this work rested on one random split and a baseline that had seen the answer;",
                     " a withdrawn regional headline rested on three countries. An audit found each, and everything on this",
                     " dashboard comes from the corrected protocol.")
      )
    })
  })
}
