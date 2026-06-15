# =============================================================================
# Module: Model diagnostics
# =============================================================================
# Surfaces the full per-model discrimination & calibration diagnostics produced
# by the pipeline's `diagnostics_all` / `diagnostics_calibrated_all` targets:
#   (1) ROC curves              — sensitivity vs 1-specificity, by outcome
#   (2) Precision–Recall curves — more informative than ROC for rare outcomes
#   (3) Calibration curves      — binned predicted vs observed prevalence
#   (4) Binary metrics table    — AUC, PR-AUC, Brier skill, ECE/MCE, calibration
#   (5) Continuous metrics      — RMSE, MAE, R^2, r
#   (6) Recalibration (Platt)   — Brier-skill / slope before vs after Platt
# Reads the optional `model_diag` object built by 01_prepare_dashboard_data.R;
# each panel renders an empty state when its source is missing.

.diag_empty <- function(msg) {
  plotly::plotly_empty(type = "scatter", mode = "markers") |>
    plotly::layout(annotations = list(text = msg, showarrow = FALSE,
                                      font = list(size = 14)))
}

# Country selector used by the curve panels (one curve per outcome within a
# selected country, mirroring the importance module's country picker).
.diag_country_choices <- function() country_choices

mod_diagnostics_ui <- function(id) {
  ns <- NS(id)

  tagList(
    p(em("For analysts: model diagnostics (ROC, calibration, and similar). The policy tabs above don't require these."),
      style = "color:#6c757d; margin:0 0 0.6em;"),
    navset_card_tab(

    nav_panel(
      title = "Discrimination (ROC)",
      icon = bsicons::bs_icon("activity"),
      div(
        layout_columns(
          col_widths = 6,
          selectInput(ns("roc_country"), "Country",
                      choices = .diag_country_choices(), selected = "ghana")
        ),
        plotlyOutput(ns("roc_plot"), height = "500px"),
        methods_note(
          "Receiver-operating-characteristic curves from cluster-blocked ",
          "cross-validated out-of-fold predictions. Each line is one binary ",
          "deficiency model for the selected country; the dashed diagonal is a ",
          "non-informative classifier (AUC 0.5). Curves bowing toward the ",
          "top-left corner indicate better discrimination. ROC-AUC values are ",
          "listed in the “Binary metrics” sub-tab."
        )
      )
    ),

    nav_panel(
      title = "Precision–Recall",
      icon = bsicons::bs_icon("graph-down"),
      div(
        layout_columns(
          col_widths = 6,
          selectInput(ns("pr_country"), "Country",
                      choices = .diag_country_choices(), selected = "ghana")
        ),
        plotlyOutput(ns("pr_plot"), height = "500px"),
        methods_note(
          "Precision–recall curves for the same out-of-fold predictions. For ",
          "rare outcomes (low prevalence) PR curves are more informative than ",
          "ROC because they ignore the large pool of true negatives. The ",
          "horizontal reference for each model is its outcome prevalence (the ",
          "precision of a random classifier)."
        )
      )
    ),

    nav_panel(
      title = "Calibration",
      icon = bsicons::bs_icon("rulers"),
      div(
        layout_columns(
          col_widths = 6,
          selectInput(ns("cal_country"), "Country",
                      choices = .diag_country_choices(), selected = "ghana")
        ),
        plotlyOutput(ns("cal_plot"), height = "500px"),
        methods_note(
          "Reliability (calibration) curves: predicted probabilities are ",
          "binned and the mean predicted probability (x) is plotted against ",
          "the observed event rate (y) in each bin. Points on the dashed 45° ",
          "line are perfectly calibrated; points below the line indicate ",
          "over-prediction, above indicate under-prediction. See the ",
          "“Recalibration (Platt)” sub-tab for the effect of post-hoc ",
          "logistic recalibration."
        )
      )
    ),

    nav_panel(
      title = "Binary metrics",
      icon = bsicons::bs_icon("table"),
      div(
        reactableOutput(ns("binary_table")),
        methods_note(
          "Per country × outcome binary-model diagnostics on cross-validated ",
          "out-of-fold predictions. ", tags$b("Brier skill"), " compares the ",
          "model's Brier score against an always-predict-prevalence baseline ",
          "(positive = better than baseline). ", tags$b("ECE/MCE"), " are the ",
          "expected and maximum calibration errors. ", tags$b("Calib. slope"),
          " near 1 and ", tags$b("intercept"), " near 0 indicate good ",
          "calibration."
        )
      )
    ),

    nav_panel(
      title = "Continuous metrics",
      icon = bsicons::bs_icon("table"),
      div(
        reactableOutput(ns("continuous_table")),
        methods_note(
          "Per country × outcome continuous-model diagnostics (predicting the ",
          "underlying biomarker value rather than the binary deficiency flag). ",
          tags$b("R²"), " is the cross-validated coefficient of determination; ",
          tags$b("r"), " the Pearson correlation between predicted and observed ",
          "values; ", tags$b("RMSE/MAE"), " the root-mean-square and mean ",
          "absolute error in the biomarker units."
        )
      )
    ),

    nav_panel(
      title = "Recalibration (Platt)",
      icon = bsicons::bs_icon("sliders2"),
      div(
        reactableOutput(ns("recal_table")),
        methods_note(
          "Effect of post-hoc ", tags$b("Platt (logistic) recalibration"),
          " on the binary models. The pipeline re-passes the out-of-fold ",
          "predicted probabilities through a fitted logistic recalibrator and ",
          "recomputes the diagnostics. This repairs models with negative Brier ",
          "skill (mis-scaled probabilities) without changing their ranking, so ",
          "ROC-AUC and PR-AUC are unchanged while Brier skill and the ",
          "calibration slope improve. Columns show the metric before → after ",
          "recalibration."
        )
      )
    )
  )
  )
}


mod_diagnostics_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    pretty_outcome <- function(oc) {
      lab <- meta$outcome_labels[oc]
      ifelse(is.na(lab), oc, lab)
    }

    # ── ROC curves ───────────────────────────────────────────────────────
    output$roc_plot <- renderPlotly({
      roc <- model_diag$roc
      ctry <- meta$countries[input$roc_country]
      if (is.null(roc) || nrow(roc) == 0)
        return(.diag_empty("ROC curve data not yet built. Re-run the pipeline / data-prep."))
      d <- roc[roc$country == ctry, , drop = FALSE]
      if (nrow(d) == 0)
        return(.diag_empty(sprintf("No ROC data for %s.", ctry)))
      d$Outcome <- pretty_outcome(d$outcome)
      plot_ly() |>
        add_trace(data = d, x = ~fpr, y = ~tpr, color = ~Outcome,
                  type = "scatter", mode = "lines",
                  hoverinfo = "text",
                  text = ~sprintf("%s<br>FPR: %.2f<br>TPR: %.2f",
                                  Outcome, fpr, tpr)) |>
        add_segments(x = 0, xend = 1, y = 0, yend = 1,
                     line = list(dash = "dot", color = "#999"),
                     showlegend = FALSE, hoverinfo = "skip", inherit = FALSE) |>
        layout(xaxis = list(title = "False positive rate (1 − specificity)",
                            range = c(0, 1)),
               yaxis = list(title = "True positive rate (sensitivity)",
                            range = c(0, 1)),
               legend = list(title = list(text = "Outcome")))
    })

    # ── Precision–Recall curves ──────────────────────────────────────────
    output$pr_plot <- renderPlotly({
      pr <- model_diag$pr
      ctry <- meta$countries[input$pr_country]
      if (is.null(pr) || nrow(pr) == 0)
        return(.diag_empty("Precision–recall data not yet built. Re-run the pipeline / data-prep."))
      d <- pr[pr$country == ctry, , drop = FALSE]
      if (nrow(d) == 0)
        return(.diag_empty(sprintf("No PR data for %s.", ctry)))
      d$Outcome <- pretty_outcome(d$outcome)
      plot_ly(d, x = ~recall, y = ~precision, color = ~Outcome,
              type = "scatter", mode = "lines",
              hoverinfo = "text",
              text = ~sprintf("%s<br>recall: %.2f<br>precision: %.2f",
                              Outcome, recall, precision)) |>
        layout(xaxis = list(title = "Recall", range = c(0, 1)),
               yaxis = list(title = "Precision", range = c(0, 1)),
               legend = list(title = list(text = "Outcome")))
    })

    # ── Calibration (reliability) curves ─────────────────────────────────
    output$cal_plot <- renderPlotly({
      cal <- model_diag$calibration
      ctry <- meta$countries[input$cal_country]
      if (is.null(cal) || nrow(cal) == 0)
        return(.diag_empty("Calibration data not yet built. Re-run the pipeline / data-prep."))
      d <- cal[cal$country == ctry, , drop = FALSE]
      if (nrow(d) == 0)
        return(.diag_empty(sprintf("No calibration data for %s.", ctry)))
      d$Outcome <- pretty_outcome(d$outcome)
      plot_ly() |>
        add_trace(data = d, x = ~mean_pred, y = ~mean_obs, color = ~Outcome,
                  type = "scatter", mode = "lines+markers",
                  hoverinfo = "text",
                  text = ~sprintf("%s<br>bin: %s<br>pred: %.2f<br>obs: %.2f<br>n: %d",
                                  Outcome, bin, mean_pred, mean_obs, n)) |>
        add_segments(x = 0, xend = 1, y = 0, yend = 1,
                     line = list(dash = "dot", color = "#999"),
                     showlegend = FALSE, hoverinfo = "skip", inherit = FALSE) |>
        layout(xaxis = list(title = "Mean predicted probability"),
               yaxis = list(title = "Observed event rate"),
               legend = list(title = list(text = "Outcome")))
    })

    # ── Binary metrics table ─────────────────────────────────────────────
    output$binary_table <- renderReactable({
      df <- model_diag$binary
      validate(need(!is.null(df) && nrow(df) > 0,
                    "Binary diagnostics not yet built. Re-run the pipeline / data-prep."))
      df$Outcome <- pretty_outcome(df$outcome)
      keep <- c("country", "Outcome", "n", "prevalence", "roc_auc", "pr_auc",
                "brier", "brier_skill", "calib_slope", "calib_int", "ece", "mce")
      keep <- intersect(keep, colnames(df))
      d <- df[, keep, drop = FALSE]
      nm <- c(country = "Country", n = "N", prevalence = "Prevalence",
              roc_auc = "ROC-AUC", pr_auc = "PR-AUC", brier = "Brier",
              brier_skill = "Brier skill", calib_slope = "Calib. slope",
              calib_int = "Calib. intercept", ece = "ECE", mce = "MCE")
      colnames(d) <- ifelse(colnames(d) %in% names(nm), nm[colnames(d)], colnames(d))
      reactable(
        d, compact = TRUE, striped = TRUE, searchable = TRUE,
        defaultPageSize = 15,
        defaultColDef = colDef(format = colFormat(digits = 3)),
        columns = list(
          Country = colDef(format = NULL),
          Outcome = colDef(format = NULL, minWidth = 180),
          N = colDef(format = colFormat(separators = TRUE, digits = 0)),
          `ROC-AUC` = colDef(format = colFormat(digits = 2),
            style = function(value) {
              if (is.na(value)) return(NULL)
              color <- if (value >= 0.80) "#1a9850"
                else if (value >= 0.70) "#91cf60"
                else if (value >= 0.60) "#fdae61" else "#d7191c"
              list(color = color, fontWeight = "bold")
            }),
          `Brier skill` = colDef(format = colFormat(digits = 3),
            style = function(value) {
              if (is.na(value)) return(NULL)
              list(color = if (value < 0) "#d7191c" else "#1a9850")
            })
        )
      )
    })

    # ── Continuous metrics table ─────────────────────────────────────────
    output$continuous_table <- renderReactable({
      df <- model_diag$continuous
      validate(need(!is.null(df) && nrow(df) > 0,
                    "Continuous diagnostics not yet built. Re-run the pipeline / data-prep."))
      df$Outcome <- pretty_outcome(df$outcome)
      keep <- intersect(c("country", "Outcome", "n", "rmse", "mae", "r2", "r"),
                        colnames(df))
      d <- df[, keep, drop = FALSE]
      nm <- c(country = "Country", n = "N", rmse = "RMSE", mae = "MAE",
              r2 = "R²", r = "r")
      colnames(d) <- ifelse(colnames(d) %in% names(nm), nm[colnames(d)], colnames(d))
      reactable(
        d, compact = TRUE, striped = TRUE, searchable = TRUE,
        defaultPageSize = 15,
        defaultColDef = colDef(format = colFormat(digits = 3)),
        columns = list(
          Country = colDef(format = NULL),
          Outcome = colDef(format = NULL, minWidth = 180),
          N = colDef(format = colFormat(separators = TRUE, digits = 0))
        )
      )
    })

    # ── Recalibration (Platt) before vs after ────────────────────────────
    output$recal_table <- renderReactable({
      raw <- model_diag$binary
      cal <- model_diag$calibrated
      validate(need(!is.null(cal) && nrow(cal) > 0,
                    "Calibrated diagnostics not yet built. Re-run the pipeline / data-prep."))
      cal$Outcome <- pretty_outcome(cal$outcome)
      d <- data.frame(
        Country = cal$country, Outcome = cal$Outcome,
        bs_after = cal$brier_skill, slope_after = cal$calib_slope,
        stringsAsFactors = FALSE
      )
      if (!is.null(raw) && nrow(raw) > 0) {
        key_raw <- paste(raw$country, raw$outcome)
        key_cal <- paste(cal$country, cal$outcome)
        idx <- match(key_cal, key_raw)
        d$bs_before    <- raw$brier_skill[idx]
        d$slope_before <- raw$calib_slope[idx]
      } else {
        d$bs_before <- NA_real_; d$slope_before <- NA_real_
      }
      d <- d[, c("Country", "Outcome", "bs_before", "bs_after",
                 "slope_before", "slope_after")]
      reactable(
        d, compact = TRUE, striped = TRUE, searchable = TRUE,
        defaultPageSize = 15,
        defaultColDef = colDef(format = colFormat(digits = 3)),
        columns = list(
          Country = colDef(format = NULL),
          Outcome = colDef(format = NULL, minWidth = 180),
          bs_before    = colDef(name = "Brier skill (raw)"),
          bs_after     = colDef(name = "Brier skill (calibrated)",
            style = function(value) {
              if (is.na(value)) return(NULL)
              list(color = if (value < 0) "#d7191c" else "#1a9850")
            }),
          slope_before = colDef(name = "Calib. slope (raw)"),
          slope_after  = colDef(name = "Calib. slope (calibrated)")
        )
      )
    })

  })
}
