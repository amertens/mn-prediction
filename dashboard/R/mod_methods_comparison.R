# =============================================================================
# Module: Methods comparison (corrected vs production)
# =============================================================================
# Surfaces the parallel "_targets_corrected" pipeline (P1-P8 fixes) head-to-head
# against the production methods, reading the bundle written by
# R_corrected/comparison.R::build_methods_comparison ->
# dashboard/data/methods_comparison.rds. Each panel renders an empty state when
# the bundle (or a component) is absent so the tab is safe before the corrected
# pipeline has run.

.mc_empty <- function(msg) {
  div(class = "text-muted", style = "padding:1em;", em(msg))
}

mod_methods_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    p(em("For analysts: corrected-vs-production methods comparison. The policy tabs above don't require these."),
      style = "color:#6c757d; margin:0 0 0.6em;"),
    navset_card_tab(
    nav_panel(
      title = "CV honesty (P1)",
      icon = bsicons::bs_icon("shield-check"),
      div(
        p("Out-of-fold discrimination under the leakage-free corrected fit ",
          "(all preprocessing fit inside each fold) versus the optimistic ",
          "regime (preprocessing on the full data + random CV) and the ",
          "production number, on identical data and learner library."),
        plotlyOutput(ns("cv_plot"), height = "330px"),
        reactableOutput(ns("cv_tbl")),
        methods_note(
          "Honest cluster- and region/spatial-block CV should sit ", tags$b("below"),
          " the optimistic and production estimates; the gap is the optimism ",
          "removed by fixing P1. Spatial-block is the most deployment-honest."
        )
      )
    ),
    nav_panel(
      title = "Calibration (P2)",
      icon = bsicons::bs_icon("rulers"),
      div(
        p("In-sample Platt recalibration (the production anti-pattern) vs ",
          "out-of-fold Platt. A slope of exactly 1.00 in-sample is an illusion; ",
          "the honest out-of-fold slope is the real calibration."),
        reactableOutput(ns("cal_tbl")),
        methods_note(
          "Calibration slope near 1 and intercept near 0 = well calibrated. ",
          "In-sample calibration is fit and scored on the same predictions, so ",
          "it is optimistic by construction."
        )
      )
    ),
    nav_panel(
      title = "Decision value (P3)",
      icon = bsicons::bs_icon("bullseye"),
      div(reactableOutput(ns("dec_tbl")),
          methods_note(
            "Targeting accuracy / lift computed directly from corrected ",
            "admin-2 predictions as first-class pipeline targets."))
    ),
    nav_panel(
      title = "Intervals (P4)",
      icon = bsicons::bs_icon("distribute-vertical"),
      div(reactableOutput(ns("int_tbl")),
          methods_note(
            "Split-conformal intervals use a held-out calibration split; ",
            "empirical held-out coverage is reported. Design-based district CIs ",
            "come from the survey directly (no Admin-1 -> Admin-2 broadcast)."))
    ),
    nav_panel(
      title = "Error (P5)",
      icon = bsicons::bs_icon("calculator"),
      div(reactableOutput(ns("err_tbl")),
          methods_note(
            "Naive RMSE treats the survey estimate as exact truth; the ",
            "sampling-adjusted RMSE removes the survey's own sampling variance ",
            "from the apparent error."))
    ),
    nav_panel(
      title = "Trust flags (P6)",
      icon = bsicons::bs_icon("shield-exclamation"),
      div(reactableOutput(ns("trust_tbl")),
          methods_note(
            "Per-district reliability from covariate out-of-support distance, ",
            "survey sample size and interval width. 'Do not rely' districts are ",
            "outside the surveyed covariate support."))
    ),
    nav_panel(
      title = "Area model (P7)",
      icon = bsicons::bs_icon("diagram-3"),
      div(reactableOutput(ns("area_tbl")),
          methods_note(
            "Design-aware partial-pooling area estimator (Fay-Herriot, or ",
            "empirical-Bayes shrinkage fallback) as a first-class estimator, ",
            "beside the raw direct survey estimate and the corrected ML estimate."))
    )
  )
  )
}

mod_methods_comparison_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    mc <- if (exists("methods_comp")) methods_comp else NULL

    num_cols <- function(df) names(df)[vapply(df, is.numeric, logical(1))]
    rt <- function(df) reactable(df, compact = TRUE, striped = TRUE,
                                 searchable = TRUE, defaultPageSize = 12,
                                 defaultColDef = colDef(
                                   format = colFormat(digits = 3)))

    output$cv_plot <- renderPlotly({
      d <- mc$cv_compare
      if (is.null(d) || nrow(d) == 0)
        return(plotly_empty() |> layout(
          annotations = list(text = "Run the corrected pipeline first.",
                             showarrow = FALSE)))
      d$lab <- paste0(d$scheme)
      plot_ly(d, x = ~auc, y = ~reorder(lab, auc), type = "bar",
              orientation = "h",
              color = ~ifelse(honest, "honest", "optimistic/production"),
              colors = c(honest = "#1a9850",
                         `optimistic/production` = "#d7191c"),
              hoverinfo = "text",
              text = ~sprintf("%s<br>AUC %.3f", scheme, auc)) |>
        layout(xaxis = list(title = "OOF ROC-AUC"),
               yaxis = list(title = ""), margin = list(l = 10),
               legend = list(orientation = "h"))
    })
    output$cv_tbl <- renderReactable({
      d <- mc$cv_compare; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$cal_tbl <- renderReactable({
      d <- mc$calibration; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$dec_tbl <- renderReactable({
      d <- mc$decision; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$int_tbl <- renderReactable({
      d <- mc$interval_summary; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$err_tbl <- renderReactable({
      d <- mc$admin2_error; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$trust_tbl <- renderReactable({
      d <- mc$trust; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
    output$area_tbl <- renderReactable({
      d <- mc$area_pp; validate(need(!is.null(d), "No corrected results yet."))
      rt(d)
    })
  })
}
