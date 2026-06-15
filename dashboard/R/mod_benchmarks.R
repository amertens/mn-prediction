# =============================================================================
# Module: Method benchmarks
# =============================================================================
# Exposes the area-level model-comparison tables the pipeline computes but that
# were not previously displayed anywhere in the dashboard:
#   (1) SAE method comparison    — benchmarks_all: SuperLearner vs country-mean,
#       OLS GLM, elastic-net, Fay-Herriot, BYM2 under leave-one-country-out.
#   (2) Individual vs area SL     — area_comparison_all: within-country accuracy
#       of individual-level SL (aggregated to Admin-2) vs area-level SL.
#   (3) Admin-2 accuracy          — admin2_error_all: per country × outcome
#       district-level prediction accuracy.
#   (4) Prescreened SL (LOCO)     — sl_prescreened: optimized primary SL.
#   (5) Area-transport (LOCO)     — transportability_area_loco_metrics: enriched
#       parsimonious area-transport model metrics per held-out country.
# Reads the optional `benchmarks_data` object built by
# 01_prepare_dashboard_data.R; each panel renders an empty state when missing.

.bench_empty <- function(msg) {
  plotly::plotly_empty(type = "scatter", mode = "markers") |>
    plotly::layout(annotations = list(text = msg, showarrow = FALSE,
                                      font = list(size = 14)))
}

.bench_pretty_outcome <- function(oc) {
  lab <- meta$outcome_labels[oc]
  ifelse(is.na(lab), oc, lab)
}

mod_benchmarks_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "SAE method comparison",
      icon = bsicons::bs_icon("trophy"),
      div(
        p("Out-of-sample accuracy of the area-level SuperLearner against ",
          "standard small-area-estimation (SAE) methods under ",
          "leave-one-country-out (LOCO) transportability: a country-mean ",
          "baseline, OLS GLM, elastic-net penalised GLM, Fay-Herriot, and BYM2 ",
          "(spatial). Lower mean-absolute-error (MAE) is better."),
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("bench_outcome"), "Outcome", choices = NULL),
          radioButtons(ns("bench_metric"), "Metric",
                       choices = c("MAE (pp)" = "mae_pp",
                                   "RMSE (pp)" = "rmse_pp",
                                   "Pearson r" = "pearson_r"),
                       selected = "mae_pp", inline = TRUE)
        ),
        plotlyOutput(ns("bench_plot"), height = "440px"),
        reactableOutput(ns("bench_table")),
        methods_note(
          "Each bar is the mean across held-out countries for one method; the ",
          "table below gives the full per-held-out-country breakdown. Metrics ",
          "are in percentage points (pp) of prevalence except Pearson r. ",
          "Methods are evaluated on identical proxy-only Admin-2 covariates."
        )
      )
    ),

    nav_panel(
      title = "Individual vs area SL",
      icon = bsicons::bs_icon("layers"),
      div(
        p("Within-country accuracy of two modelling strategies at Admin-2: the ",
          "individual-level SuperLearner aggregated up to district prevalence, ",
          "versus the area-level SuperLearner that predicts district prevalence ",
          "directly from aggregated covariates (under MSE and logit/NLL loss)."),
        reactableOutput(ns("area_comp_table")),
        methods_note(
          "Computed on each country's surveyed Admin-2 districts. ",
          tags$b("MAE/RMSE"), " are in percentage points; ", tags$b("r"),
          " is the Pearson correlation between predicted and survey-observed ",
          "district prevalence; ", tags$b("bias"), " is mean signed error."
        )
      )
    ),

    nav_panel(
      title = "Admin-2 accuracy",
      icon = bsicons::bs_icon("geo"),
      div(
        p("Per country × outcome accuracy of the within-country Admin-2 ",
          "predictions against survey-observed district prevalence."),
        reactableOutput(ns("admin2_err_table")),
        methods_note(
          "Summarises how closely the model reproduces survey-observed Admin-2 ",
          "prevalence in the surveyed districts. ", tags$b("MAE/RMSE"),
          " in percentage points, ", tags$b("r"), " Pearson correlation, ",
          tags$b("bias"), " mean signed error (positive = over-prediction)."
        )
      )
    ),

    nav_panel(
      title = "Prescreened SL (LOCO)",
      icon = bsicons::bs_icon("funnel"),
      div(
        p("The optimized primary SuperLearner workflow: five-stage feature ",
          "prescreening + spatial coordinates + a fast learner library, ",
          "evaluated under leave-one-country-out transportability."),
        reactableOutput(ns("sl_pre_table")),
        methods_note(
          "One row per held-out country. This is the fast main-analysis SL ",
          "(~30 s per LOCO fold) used in place of the historical full SL. ",
          tags$b("MAE/RMSE"), " in percentage points; ", tags$b("calib. slope"),
          " near 1 indicates well-scaled transported predictions."
        )
      )
    ),

    nav_panel(
      title = "Area-transport (LOCO)",
      icon = bsicons::bs_icon("box-arrow-up-right"),
      div(
        p("The enriched, parsimonious within-country-centered elastic-net ",
          "area-transport model (harmonized GEE + IHME + Malaria-Atlas + ",
          "food-security proxies), evaluated per held-out country."),
        reactableOutput(ns("area_trans_table")),
        methods_note(
          "Metrics from the ", tags$code("area_transport_summary"), " target. ",
          tags$b("n selected"), " is the number of predictors retained after ",
          "within-country centering and elastic-net selection; ",
          tags$b("MAE/RMSE"), " in percentage points; ", tags$b("calib. slope"),
          " summarises transport calibration."
        )
      )
    )
  )
}


mod_benchmarks_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── SAE method comparison ────────────────────────────────────────────
    bench_df <- reactive({
      b <- benchmarks_data$benchmarks
      if (is.null(b) || nrow(b) == 0) return(NULL)
      b
    })

    observe({
      b <- bench_df()
      if (is.null(b)) return()
      ocs <- intersect(names(meta$outcome_labels), unique(b$outcome))
      choices <- setNames(ocs, meta$outcome_labels[ocs])
      updateSelectInput(session, "bench_outcome", choices = choices,
                        selected = choices[1])
    })

    output$bench_plot <- renderPlotly({
      b <- bench_df()
      if (is.null(b))
        return(.bench_empty("Benchmark data not yet built. Re-run the pipeline / data-prep."))
      req(input$bench_outcome)
      metric <- input$bench_metric %||% "mae_pp"
      d <- b[b$outcome == input$bench_outcome, , drop = FALSE]
      if (nrow(d) == 0)
        return(.bench_empty("No benchmark rows for this outcome."))
      # Mean across held-out countries per method.
      agg <- aggregate(d[[metric]], by = list(method = d$method),
                       FUN = function(x) mean(x, na.rm = TRUE))
      colnames(agg)[2] <- "value"
      asc <- metric != "pearson_r"   # lower is better except correlation
      agg <- agg[order(agg$value, decreasing = !asc), ]
      agg$method <- factor(agg$method, levels = agg$method)
      ylab <- c(mae_pp = "Mean MAE (pp)", rmse_pp = "Mean RMSE (pp)",
                pearson_r = "Mean Pearson r")[metric]
      plot_ly(agg, x = ~method, y = ~value, type = "bar",
              marker = list(color = "#2c7bb6"),
              hoverinfo = "text",
              text = ~sprintf("%s<br>%s: %.2f", method, ylab, value)) |>
        layout(xaxis = list(title = "Method"),
               yaxis = list(title = ylab))
    })

    output$bench_table <- renderReactable({
      b <- bench_df()
      validate(need(!is.null(b), ""))
      req(input$bench_outcome)
      d <- b[b$outcome == input$bench_outcome, , drop = FALSE]
      keep <- intersect(c("method", "model_type", "held_out", "n_train",
                          "n_test", "obs_prev", "pred_prev", "mae_pp",
                          "rmse_pp", "pearson_r", "mean_bias_pp", "calib_slope"),
                        colnames(d))
      reactable(d[, keep, drop = FALSE], compact = TRUE, striped = TRUE,
                searchable = TRUE, defaultPageSize = 12,
                defaultColDef = colDef(format = colFormat(digits = 3)))
    })

    # ── Individual vs area SL (within-country) ────────────────────────────
    output$area_comp_table <- renderReactable({
      d <- benchmarks_data$area_comparison
      validate(need(!is.null(d) && nrow(d) > 0,
                    "Area comparison not yet built. Re-run the pipeline / data-prep."))
      d$Outcome <- .bench_pretty_outcome(d$outcome)
      keep <- intersect(c("country", "Outcome", "approach", "n_admin2",
                          "mae_pp", "rmse_pp", "pearson_r", "mean_bias_pp"),
                        colnames(d))
      dd <- d[, keep, drop = FALSE]
      nm <- c(country = "Country", approach = "Approach", n_admin2 = "N districts",
              mae_pp = "MAE (pp)", rmse_pp = "RMSE (pp)", pearson_r = "r",
              mean_bias_pp = "Bias (pp)")
      colnames(dd) <- ifelse(colnames(dd) %in% names(nm), nm[colnames(dd)], colnames(dd))
      reactable(dd, compact = TRUE, striped = TRUE, searchable = TRUE,
                defaultPageSize = 15, groupBy = "Outcome",
                defaultColDef = colDef(format = colFormat(digits = 3)),
                columns = list(
                  Country = colDef(format = NULL),
                  Approach = colDef(format = NULL),
                  `N districts` = colDef(format = colFormat(digits = 0))
                ))
    })

    # ── Admin-2 accuracy ─────────────────────────────────────────────────
    output$admin2_err_table <- renderReactable({
      d <- benchmarks_data$admin2_error
      validate(need(!is.null(d) && nrow(d) > 0,
                    "Admin-2 error table not yet built. Re-run the pipeline / data-prep."))
      d$Outcome <- .bench_pretty_outcome(d$outcome)
      keep <- intersect(c("country", "Outcome", "n_admin2", "mae_pp", "rmse_pp",
                          "pearson_r", "mean_bias_pp"), colnames(d))
      dd <- d[, keep, drop = FALSE]
      nm <- c(country = "Country", n_admin2 = "N districts", mae_pp = "MAE (pp)",
              rmse_pp = "RMSE (pp)", pearson_r = "r", mean_bias_pp = "Bias (pp)")
      colnames(dd) <- ifelse(colnames(dd) %in% names(nm), nm[colnames(dd)], colnames(dd))
      reactable(dd, compact = TRUE, striped = TRUE, searchable = TRUE,
                defaultPageSize = 15,
                defaultColDef = colDef(format = colFormat(digits = 3)),
                columns = list(
                  Country = colDef(format = NULL),
                  Outcome = colDef(format = NULL, minWidth = 180),
                  `N districts` = colDef(format = colFormat(digits = 0))
                ))
    })

    # ── Prescreened SL (LOCO) ────────────────────────────────────────────
    output$sl_pre_table <- renderReactable({
      d <- benchmarks_data$sl_prescreened
      validate(need(!is.null(d) && nrow(d) > 0,
                    "Prescreened-SL LOCO table not yet built. Re-run the pipeline / data-prep."))
      d$Outcome <- .bench_pretty_outcome(d$outcome)
      keep <- intersect(c("Outcome", "held_out", "n_train", "n_test", "n_vars",
                          "obs_prev", "pred_prev", "mae_pp", "rmse_pp",
                          "pearson_r", "spearman_r", "calib_slope"),
                        colnames(d))
      dd <- d[, keep, drop = FALSE]
      nm <- c(held_out = "Held-out", n_train = "N train", n_test = "N test",
              n_vars = "N vars", obs_prev = "Obs prev", pred_prev = "Pred prev",
              mae_pp = "MAE (pp)", rmse_pp = "RMSE (pp)", pearson_r = "r",
              spearman_r = "Spearman", calib_slope = "Calib. slope")
      colnames(dd) <- ifelse(colnames(dd) %in% names(nm), nm[colnames(dd)], colnames(dd))
      reactable(dd, compact = TRUE, striped = TRUE, searchable = TRUE,
                defaultPageSize = 15,
                defaultColDef = colDef(format = colFormat(digits = 3)),
                columns = list(
                  Outcome = colDef(format = NULL, minWidth = 180),
                  `Held-out` = colDef(format = NULL)
                ))
    })

    # ── Area-transport (LOCO) ────────────────────────────────────────────
    output$area_trans_table <- renderReactable({
      d <- benchmarks_data$area_transport
      validate(need(!is.null(d) && nrow(d) > 0,
                    "Area-transport LOCO metrics not yet built. Re-run the pipeline / data-prep."))
      d$Outcome <- .bench_pretty_outcome(d$outcome)
      keep <- intersect(c("Outcome", "held_out", "n_train", "n_test", "n_pred",
                          "n_selected", "pearson_r", "spearman_r", "rmse_pp",
                          "mae_pp", "nat_bias_pp", "calib_slope"),
                        colnames(d))
      dd <- d[, keep, drop = FALSE]
      nm <- c(held_out = "Held-out", n_train = "N train", n_test = "N test",
              n_pred = "N pred", n_selected = "N selected", pearson_r = "r",
              spearman_r = "Spearman", rmse_pp = "RMSE (pp)", mae_pp = "MAE (pp)",
              nat_bias_pp = "Nat. bias (pp)", calib_slope = "Calib. slope")
      colnames(dd) <- ifelse(colnames(dd) %in% names(nm), nm[colnames(dd)], colnames(dd))
      reactable(dd, compact = TRUE, striped = TRUE, searchable = TRUE,
                defaultPageSize = 15,
                defaultColDef = colDef(format = colFormat(digits = 3)),
                columns = list(
                  Outcome = colDef(format = NULL, minWidth = 180),
                  `Held-out` = colDef(format = NULL)
                ))
    })

  })
}
