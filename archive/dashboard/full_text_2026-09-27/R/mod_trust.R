# =============================================================================
# Module: How well it works
# =============================================================================
# The three tests, the best score the survey's own noise allows, what each
# added survey buys, the comparison with a geostatistical model, and what else
# was tried. Every figure is drawn from the committed result tables.

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
      methods_note("Ranking accuracy is the correlation (Spearman) between each method's order of districts and the",
                   " survey's, averaged over country-outcome pairs, with a 95% interval. Random splits are repeated ten",
                   " times. Every method gets the same information as the model: for example, the survey's regional",
                   " average leaves out the district being predicted. The dashed line is what a random ranking reaches",
                   " in 95% of tries. Differences below 0.03 are ties.")
    ),
    nav_panel(
      title = "Best achievable score", icon = bsicons::bs_icon("arrow-bar-up"),
      layout_columns(col_widths = c(3, 9),
        div(radioButtons(ns("ceil_target"), "Outcome measured as", choices = c("Prevalence" = "prev", "Biomarker level" = "level"), selected = "prev"),
            uiOutput(ns("ceiling_text"))),
        plotlyOutput(ns("ceiling"), height = "380px")),
      methods_note("Most districts have one or two survey clusters, so the survey's figure for a district is itself",
                   " uncertain. A statistical model of how results vary between people, clusters, districts and regions",
                   " gives the best score any model could reach against those figures. The 'split-half' estimate counts",
                   " differences between clusters as differences between places and is too high; the lower estimate",
                   " removes them.")
    ),
    nav_panel(
      title = "Each survey added", icon = bsicons::bs_icon("graph-up-arrow"),
      layout_columns(col_widths = c(3, 9),
        div(uiOutput(ns("curve_text"))),
        plotlyOutput(ns("curve"), height = "380px")),
      methods_note("Every combination of training countries, scored on each country left out (biomarker levels). Bars",
                   " show the middle half of the results. The climate-and-soil version was chosen using these same four",
                   " countries, so its result on a fifth country is a prediction that has not yet been tested.")
    ),
    nav_panel(
      title = "Compared with a geostatistical model", icon = bsicons::bs_icon("geo"),
      layout_columns(col_widths = c(5, 7),
        div(uiOutput(ns("geostat_text"))),
        reactableOutput(ns("geostat"))),
      methods_note("A geostatistical model combines a smooth map of location with data layers and is fitted to the",
                   " locations of survey clusters. It was scored on the same districts as this model, once with our data",
                   " layers and once with the DHS Program's own list of ten layers. It needs survey clusters inside the",
                   " country, so it cannot be used for a country without a survey.")
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
                   sprintf(paste("Inside a surveyed country, this model (%s), averaging neighbouring districts, and the two",
                                 "combined score about the same, and all do better than the survey's own regional averages (%s).",
                                 "Where a survey exists, location explains much of the pattern."),
                           fmt_num(Q$infill), fmt_num(Q$infill_jk))),
        region = p(style = "font-size:0.88em; color:#555;",
                   sprintf(paste("With a whole region hidden, averaging neighbours has to reach beyond the edge of the data.",
                                 "The model scores %s on biomarker levels, close to neighbour averaging. The survey's regional",
                                 "average cannot be used for a region the survey did not sample."), fmt_num(Q$region))),
        country = p(style = "font-size:0.88em; color:#555;",
                    sprintf(paste("With a whole country left out, only methods that use public data can produce a ranking:",
                                  "%s with all data layers and %s with climate and soil only, compared with %s for a random",
                                  "ranking. The ranking carries over to a new country; prevalence levels do not, so a new",
                                  "country needs at least a national survey figure."),
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
        sprintf(paste("On prevalence, the best achievable score averages %s and the model reaches %s; on biomarker levels,",
                      "%s and %s. The best achievable score is limited by noise in the survey's own district figures, so",
                      "more survey clusters per district would raise it."),
                fmt_num(Q$ceiling_prev), fmt_num(Q$infill_prev), fmt_num(Q$ceiling_level), fmt_num(Q$infill)))
    })

    output$ceiling <- renderPlotly({
      VC <- EV$ceiling; validate(need(!is.null(VC), "Table not built."))
      v <- VC[VC$rung == "admin2" & VC$target == input$ceil_target, ]
      d <- v |> group_by(country) |> summarise(`Split-half estimate (too high)` = mean(ceiling_within_vc, na.rm = TRUE),
                                               `Best achievable score` = mean(ceiling_vc, na.rm = TRUE),
                                               `This model, district hidden` = mean(achieved_spearman, na.rm = TRUE), .groups = "drop") |>
        pivot_longer(-country, names_to = "quantity", values_to = "v")
      d <- d[is.finite(d$v), ]
      cols <- c("Split-half estimate (too high)" = "#bdbdbd", "Best achievable score" = PROXY_COL, "This model, district hidden" = SURVEY_COL)
      plot_ly(d, x = ~v, y = ~country, color = ~quantity, colors = cols, type = "scatter", mode = "markers", marker = list(size = 13),
              text = ~sprintf("%s<br>%s: %.2f", country, quantity, v), hoverinfo = "text") |>
        layout(xaxis = list(title = "Ranking accuracy (average over outcomes)", range = c(0, 1)), yaxis = list(title = ""),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$curve_text <- renderUI({
      p(style = "font-size:0.88em; color:#555;",
        sprintf(paste("Accuracy in a country left out of training rose from %s with one training survey to %s with three,",
                      "about %s for each survey added. With only three points we cannot tell yet where this levels off.",
                      "The climate-and-soil version starts higher and is the version that will be tested on the next country."),
                fmt_num(Q$lc1), fmt_num(Q$lc_max), fmt_num(Q$lc_step)))
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
              text = ~sprintf("%s<br>%d training surveys: %.2f (middle half %.2f to %.2f)", set, n_train_countries, m, lo, hi), hoverinfo = "text") |>
        layout(xaxis = list(title = "Number of training surveys", dtick = 1), yaxis = list(title = "Ranking accuracy in the country left out", rangemode = "tozero"),
               shapes = list(list(type = "rect", x0 = 0.5, x1 = 3.5, y0 = 0, y1 = Q$null_d, fillcolor = "#e0e0e0", line = list(width = 0), layer = "below")),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$geostat_text <- renderUI({
      tagList(
        p(style = "font-size:0.9em;",
          sprintf(paste("A geostatistical model, the kind the DHS Program uses, ranks districts less accurately than this model",
                        "(%s against %s inside a surveyed country) but estimates the percentage more accurately (%s against %s",
                        "percentage points of error). It estimates a level directly and pulls uncertain districts toward it."),
                  fmt_num(Q$mbg_rank), fmt_num(Q$mbg_rank_index), fmt_num(Q$mbg_err, 1), fmt_num(Q$mbg_err_index, 1))),
        p(style = "font-size:0.9em;", "For a published prevalence figure in a surveyed district, the geostatistical model is the",
          " better tool, and it comes with an interval. For deciding which districts to look at first, or for a country",
          " without a survey, the ranking is more useful."))
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
                  `Error (percentage points)` = ifelse(target == "prev", round(err, 1), NA), `Country-outcome pairs` = n) |>
        arrange(Test, `Measured as`, desc(`Ranking accuracy`))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, groupBy = c("Test", "Measured as"))
    })

    output$tried <- renderUI({
      WS <- EV$weight_sources
      rows <- if (!is.null(WS)) {
        lab <- c(domain_index = "This model (weights from rank correlations)", ridge_min = "Ridge regression",
                 domain_enet = "Elastic net", lasso_min = "Lasso",
                 index_decor = "Weights adjusted for overlap between groups", index_soft1 = "Small weights set to zero",
                 sparse20 = "Equal-weight score from 20 layers")
        w <- WS[WS$arm %in% names(lab) & WS$target == "level", ] |> select(arm, estimand, mean_spearman) |>
          pivot_wider(names_from = estimand, values_from = mean_spearman)
        w$Method <- lab[w$arm]; w <- w[order(match(w$arm, names(lab))), ]
        tags$table(class = "table table-sm", style = "font-size:0.9em;",
                   tags$thead(tags$tr(tags$th("Ways of weighting the same data groups"), tags$th("District hidden"), tags$th("Region hidden"), tags$th("Country left out"))),
                   tags$tbody(lapply(seq_len(nrow(w)), function(i) tags$tr(tags$td(w$Method[i]), tags$td(fmt_num(w$infill[i])), tags$td(fmt_num(w$region[i])), tags$td(fmt_num(w$country[i]))))))
      } else NULL
      tagList(
        p(class = "lead", "More complex methods did not do better. The likely reasons are the small number of districts per country and the noise in the survey's district figures."),
        tags$ul(
          tags$li(sprintf(paste("Machine-learning ensembles (SuperLearner, with 12 to 16 learners) at best matched this model, and did",
                                "worse when tuned to minimise squared error. With 14 to 87 districts per country, methods that tune",
                                "their own settings did not do better. Predicting which individual people are deficient reached an",
                                "AUC of %s, close to a coin toss, so it is not offered here."), fmt_num(Q$il_auc))),
          tags$li("Other ways of weighting the same data groups (ridge regression, lasso and elastic net) did slightly worse within a",
                  " country. Two variants did about 0.03 better in new countries; they will be tested on the next country rather",
                  " than adopted now."),
          tags$li("Fitting the model to individual survey clusters instead of districts did not do better in any test."),
          tags$li("Adding further data (livestock, distance to water and coast, intestinal worms, updated health maps, and food",
                  " prices and temperatures at the time of fieldwork) changed accuracy by less than 0.01 each. Climate and soil",
                  " carry most of the signal.")
        ),
        rows,
        methods_note("An earlier version of this work relied on a single random split and on a comparison method that had",
                     " seen the answer, and an earlier regional result was withdrawn because a spelling mismatch had dropped",
                     " one country. An audit found these problems. Everything on this dashboard comes from the corrected analysis.")
      )
    })
  })
}
