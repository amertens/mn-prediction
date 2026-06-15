# =============================================================================
# Module: National Burden Estimator
# =============================================================================
# Country-level summary table with population-weighted prevalence and
# absolute counts of population at risk. Includes a "hidden burden" indicator
# showing how much of the national burden is concentrated in specific
# districts that diverge from the national average.

mod_national_burden_ui <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = 12,

    # ── Lead headline: Hidden burden indicator ──
    div(
      class = "alert alert-primary",
      style = "margin-bottom: 1em; border-left: 5px solid #2c7bb6;",
      h4(bsicons::bs_icon("bullseye"),
         " Hidden burden",
         style = "margin-top: 0; color: #2c7bb6;"),
      p(strong("How much burden is hidden by national averages?"),
        " For each outcome, this section counts the people living in districts ",
        "where local prevalence is clearly different from the country average, ",
        "where a single national number misrepresents what they actually face."),
      htmlOutput(ns("headline_hidden"))
    ),

    card(
      card_header("Hidden burden: districts diverging from the national average"),
      card_body(
        p("Districts above the threshold are where the national estimate ",
          tags$em("undercounts"), " local burden; districts below are where ",
          "the national estimate may ", tags$em("overcount"), " local burden. ",
          "When the “% population mischaracterized” is large, district-level ",
          "estimates substantially improve resource allocation."),
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("country"), "Country",
                      choices = country_choices,
                      selected = "ghana"),
          radioButtons(ns("pop_year"), "Population year",
                       choices = c("Survey year (matched to data)" = "survey",
                                   "2023 projection (current burden)" = "2023"),
                       selected = "survey",
                       inline = FALSE)
        ),
        sliderInput(ns("threshold"), "Deviation threshold (percentage points)",
                    min = 1, max = 20, value = 5, step = 1),
        reactableOutput(ns("hidden_table")),
        downloadButton(ns("download_hidden"), "Download CSV",
                        class = "btn-sm btn-outline-primary",
                        style = "margin-top: 8px;"),
        methods_note(
          "For each country × outcome, we compute the population-weighted national prevalence and ",
          "compare each district's predicted prevalence to that average. Districts whose ",
          "predicted prevalence differs from the national average by more than the chosen ",
          "threshold are flagged. We then sum the populations of flagged districts to ",
          "estimate the “hidden burden population”: people whose local risk is ",
          "poorly represented by a single national number. ",
          tags$br(), tags$br(),
          "Population denominators are subgroup-specific: counts for child outcomes use ",
          "estimated children 6–59 months populations; counts for women's outcomes use ",
          "women 15–49 years. Subgroup shares are taken from UN World Population Prospects ",
          "2022 country-level estimates and applied uniformly to districts within each country. ",
          tags$br(), tags$br(),
          "This is a useful indicator for funders and policymakers because it quantifies ",
          "the cost of relying on national averages for sub-national targeting decisions. ",
          "When the hidden burden population is large, district-level estimates substantially ",
          "improve resource allocation; when it is small, national estimates are adequate."
        )
      )
    ),

    card(
      card_header("National prevalence and burden estimates"),
      card_body(
        p("Population-weighted prevalence and estimated counts of affected individuals, ",
          "aggregated from district-level model predictions. Population denominators ",
          "are subgroup-specific (children 6–59 months for child outcomes; ",
          "women 15–49 years for women's outcomes)."),
        reactableOutput(ns("burden_table")),
        downloadButton(ns("download_burden"), "Download CSV",
                        class = "btn-sm btn-outline-primary",
                        style = "margin-top: 8px;"),
        methods_note(
          "National prevalence is computed as a subgroup-population-weighted average of ",
          "district-level predicted prevalences. The conformal 95% prediction intervals are ",
          "propagated by aggregating district-level CIs to Admin-1 regions, then combining ",
          "across regions with weighted variance assuming inter-region independence — a ",
          "conservative assumption that may slightly overstate national-level uncertainty. ",
          tags$br(), tags$br(),
          "“Population affected” is the product of national prevalence and the relevant ",
          "subgroup population (children 6–59 months or women 15–49 years), summed across ",
          "districts. The subgroup population is computed by applying country-specific ",
          "demographic shares from UN World Population Prospects 2022 to the WorldPop ",
          "total Admin-2 population. ",
          tags$br(), tags$br(),
          "Both prevalence and population denominators carry their own uncertainty (conformal ",
          "interval for prevalence; WorldPop modeled estimates for population), so affected ",
          "counts should be interpreted as planning estimates rather than precise counts."
        )
      )
    ),

    card(
      card_header("Comparison: pipeline national estimates vs. survey-weighted observed"),
      card_body(
        p("Pipeline-predicted national prevalence compared to the survey-weighted observed ",
          "prevalence (using the official survey weights from the original biomarker survey). ",
          "Discrepancies indicate either model bias or differences between the predicted and ",
          "surveyed populations."),
        reactableOutput(ns("compare_table")),
        downloadButton(ns("download_compare"), "Download CSV",
                        class = "btn-sm btn-outline-primary",
                        style = "margin-top: 8px;"),
        methods_note(
          "The “Survey observed” column uses the original survey weights and design ",
          "(stratification, clustering) to produce a design-based national prevalence ",
          "estimate, computed via the ", tags$code("srvyr"), " package. This is the ",
          "best direct measure we have, though district estimates from small ",
          "samples carry their own sampling error. ",
          tags$br(), tags$br(),
          "The “Model predicted” column is the population-weighted aggregation of ",
          "district-level model predictions. The “Difference” column is the model ",
          "estimate minus the survey estimate, expressed in percentage points. ",
          tags$br(), tags$br(),
          "Systematic positive differences across countries indicate the model ",
          "tends to overpredict; differences scattered around zero indicate ",
          "well-calibrated national-level performance even where district-level ",
          "errors exist (errors that average out when aggregated)."
        )
      )
    )
  )
}


mod_national_burden_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Headline summary at top of the tab
    output$headline_hidden <- renderUI({
      req(input$country, input$threshold)
      ctry_label <- meta$countries[input$country]
      thr <- input$threshold / 100

      # Compute total mischaracterized population across all outcomes
      total_pop_misch <- 0
      total_pop <- 0
      n_outcomes_assessed <- 0
      for (oc_key in names(meta$outcome_labels)) {
        df <- get_country_admin2(input$country, oc_key,
                                  admin2_bnds, admin2_pred, admin2_pop,
                                  pop_year = input$pop_year %||% "survey")
        if (is.null(df) || all(is.na(df$pred_prev))) next
        n_outcomes_assessed <- n_outcomes_assessed + 1
        df_clean <- as.data.frame(df)
        df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                              !is.na(df_clean$population), , drop = FALSE]
        natl <- national_aggregate(df_clean)
        natl_prev <- natl$pred_prev_natl
        if (is.na(natl_prev)) next

        diverging <- abs(df_clean$pred_prev - natl_prev) > thr
        total_pop_misch <- total_pop_misch +
          sum(df_clean$population[diverging], na.rm = TRUE)
        total_pop <- total_pop + sum(df_clean$population, na.rm = TRUE)
      }

      avg_pop_misch_per_outcome <- if (n_outcomes_assessed > 0)
        total_pop_misch / n_outcomes_assessed else 0
      pct_misch <- if (total_pop > 0) total_pop_misch / total_pop else 0

      tagList(
        p(strong(ctry_label), " — averaged across ", n_outcomes_assessed,
          " outcomes at the ", input$threshold, "pp threshold:"),
        p(
          tags$span(style = "font-size: 1.4em; color: #2c7bb6; font-weight: bold;",
                    fmt_count(avg_pop_misch_per_outcome)),
          " individuals per outcome live in districts whose burden is ",
          "materially mischaracterized by national averages — ",
          tags$strong(sprintf("%.0f%% of population", pct_misch * 100)),
          " on average across the assessed outcomes."
        )
      )
    })

    # Compute burden table for selected country
    burden_data <- reactive({
      req(input$country)
      ctry_label <- meta$countries[input$country]

      rows <- list()
      for (oc_key in names(meta$outcome_labels)) {
        df <- get_country_admin2(input$country, oc_key,
                                  admin2_bnds, admin2_pred, admin2_pop,
                                  pop_year = input$pop_year %||% "survey")
        if (is.null(df) || all(is.na(df$pred_prev))) next

        natl <- national_aggregate(df)

        rows[[oc_key]] <- data.frame(
          Outcome      = meta$outcome_labels[oc_key],
          Population   = natl$pop_total,
          Prevalence   = natl$pred_prev_natl,
          CI_low       = natl$ci_lo_natl,
          CI_high      = natl$ci_hi_natl,
          Affected     = natl$pop_at_risk_natl,
          stringsAsFactors = FALSE
        )
      }
      bind_rows(rows)
    })

    output$burden_table <- renderReactable({
      df <- burden_data()
      req(nrow(df) > 0)

      reactable(
        df,
        compact = TRUE,
        striped = TRUE,
        defaultPageSize = 10,
        columns = list(
          Outcome    = colDef(name = "Outcome", minWidth = 200),
          Population = colDef(name = "Total population",
                              format = colFormat(separators = TRUE, digits = 0)),
          Prevalence = colDef(name = "Prevalence",
                              format = colFormat(percent = TRUE, digits = 1),
                              style = function(value) {
                                if (is.na(value)) return(NULL)
                                color <- if (value > 0.20) "#d7191c"
                                  else if (value > 0.10) "#fdae61" else "#2c7bb6"
                                list(color = color, fontWeight = "bold")
                              }),
          CI_low     = colDef(name = "95% CI low",
                              format = colFormat(percent = TRUE, digits = 1)),
          CI_high    = colDef(name = "95% CI high",
                              format = colFormat(percent = TRUE, digits = 1)),
          Affected   = colDef(name = "Population affected",
                              format = colFormat(separators = TRUE, digits = 0))
        )
      )
    })

    # Hidden burden table
    hidden_data <- reactive({
      req(input$country, input$threshold)
      ctry_label <- meta$countries[input$country]
      thr <- input$threshold / 100

      rows <- list()
      for (oc_key in names(meta$outcome_labels)) {
        df <- get_country_admin2(input$country, oc_key,
                                  admin2_bnds, admin2_pred, admin2_pop,
                                  pop_year = input$pop_year %||% "survey")
        if (is.null(df) || all(is.na(df$pred_prev))) next

        natl <- national_aggregate(df)
        natl_prev <- natl$pred_prev_natl
        if (is.na(natl_prev)) next

        df_clean <- as.data.frame(df)
        df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                              !is.na(df_clean$population), , drop = FALSE]

        above_idx <- df_clean$pred_prev > natl_prev + thr
        below_idx <- df_clean$pred_prev < natl_prev - thr

        n_above <- sum(above_idx, na.rm = TRUE)
        n_below <- sum(below_idx, na.rm = TRUE)
        pop_above <- sum(df_clean$population[above_idx], na.rm = TRUE)
        pop_below <- sum(df_clean$population[below_idx], na.rm = TRUE)
        total_pop <- sum(df_clean$population, na.rm = TRUE)

        rows[[oc_key]] <- data.frame(
          Outcome     = meta$outcome_labels[oc_key],
          National_prev = natl_prev,
          N_above     = n_above,
          Pop_above   = pop_above,
          N_below     = n_below,
          Pop_below   = pop_below,
          Pct_hidden  = (pop_above + pop_below) / total_pop,
          stringsAsFactors = FALSE
        )
      }
      bind_rows(rows)
    })

    output$hidden_table <- renderReactable({
      df <- hidden_data()
      req(nrow(df) > 0)

      reactable(
        df,
        compact = TRUE,
        striped = TRUE,
        defaultPageSize = 10,
        columns = list(
          Outcome       = colDef(name = "Outcome", minWidth = 200),
          National_prev = colDef(name = "National prev.",
                                  format = colFormat(percent = TRUE, digits = 1)),
          N_above       = colDef(name = "Districts above"),
          Pop_above     = colDef(name = "Pop. above",
                                  format = colFormat(separators = TRUE, digits = 0)),
          N_below       = colDef(name = "Districts below"),
          Pop_below     = colDef(name = "Pop. below",
                                  format = colFormat(separators = TRUE, digits = 0)),
          Pct_hidden    = colDef(name = "% pop. mischaracterized",
                                  format = colFormat(percent = TRUE, digits = 1),
                                  style = function(value) {
                                    if (is.na(value)) return(NULL)
                                    color <- if (value > 0.4) "#d7191c"
                                      else if (value > 0.2) "#fdae61" else "#2c7bb6"
                                    list(color = color, fontWeight = "bold")
                                  })
        )
      )
    })

    # Comparison table
    output$compare_table <- renderReactable({
      req(input$country)
      ctry_label <- meta$countries[input$country]

      df <- natl_est[natl_est$country == ctry_label, , drop = FALSE]

      # Standardize columns
      obs_col <- intersect(c("obs_prev", "observed_prev"), colnames(df))[1]
      pred_col <- intersect(c("pred_prev", "predicted_prev"), colnames(df))[1]
      ci_lo_col <- intersect(c("ci_lo", "lower"), colnames(df))[1]
      ci_hi_col <- intersect(c("ci_hi", "upper"), colnames(df))[1]

      out <- data.frame(
        Outcome   = meta$outcome_labels[df$outcome],
        N         = df$n,
        Observed  = if (!is.na(obs_col)) df[[obs_col]] else NA_real_,
        Predicted = if (!is.na(pred_col)) df[[pred_col]] else NA_real_,
        CI_low    = if (!is.na(ci_lo_col)) df[[ci_lo_col]] else NA_real_,
        CI_high   = if (!is.na(ci_hi_col)) df[[ci_hi_col]] else NA_real_,
        stringsAsFactors = FALSE
      )
      out$Difference_pp <- (out$Predicted - out$Observed) * 100

      reactable(
        out,
        compact = TRUE,
        striped = TRUE,
        defaultPageSize = 10,
        columns = list(
          Outcome       = colDef(name = "Outcome", minWidth = 200),
          N             = colDef(name = "Survey n"),
          Observed      = colDef(name = "Survey observed",
                                  format = colFormat(percent = TRUE, digits = 1)),
          Predicted     = colDef(name = "Model predicted",
                                  format = colFormat(percent = TRUE, digits = 1)),
          CI_low        = colDef(name = "Pred. CI low",
                                  format = colFormat(percent = TRUE, digits = 1)),
          CI_high       = colDef(name = "Pred. CI high",
                                  format = colFormat(percent = TRUE, digits = 1)),
          Difference_pp = colDef(name = "Difference (pp)",
                                  format = colFormat(digits = 1),
                                  style = function(value) {
                                    if (is.na(value)) return(NULL)
                                    color <- if (abs(value) > 5) "#d7191c"
                                      else if (abs(value) > 2) "#fdae61" else "#2c7bb6"
                                    list(color = color)
                                  })
        )
      )
    })

    # ── Download handlers ─────────────────────────────────────────────────
    output$download_burden <- downloadHandler(
      filename = function() {
        sprintf("burden_%s_%s.csv", input$country, Sys.Date())
      },
      content = function(file) {
        write.csv(burden_data(), file, row.names = FALSE)
      }
    )
    output$download_hidden <- downloadHandler(
      filename = function() {
        sprintf("hidden_burden_%s_thr%dpp_%s.csv",
                input$country, input$threshold, Sys.Date())
      },
      content = function(file) {
        write.csv(hidden_data(), file, row.names = FALSE)
      }
    )
    output$download_compare <- downloadHandler(
      filename = function() {
        sprintf("national_compare_%s_%s.csv", input$country, Sys.Date())
      },
      content = function(file) {
        ctry_label <- meta$countries[input$country]
        write.csv(natl_est[natl_est$country == ctry_label, , drop = FALSE],
                   file, row.names = FALSE)
      }
    )
  })
}
