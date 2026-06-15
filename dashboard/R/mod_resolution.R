# =============================================================================
# Module: Spatial resolution comparison (Admin-2 vs survey cluster)
# =============================================================================
# Surfaces the cluster-resolution sensitivity analysis:
#   (1) Within-country comparison — for each country × outcome, cross-validated
#       predictive skill (Pearson r) and error (MAE, pp) when prevalence is
#       modeled at administrative-2 resolution vs. the finer survey-cluster GPS
#       buffers. A scatter against the identity line shows where finer
#       resolution helps, hurts, or makes no difference.
#   (2) Pooled leave-one-country-out transportability at each resolution.
#
# Both views read the optional `cluster_res` object (built by
# 01_prepare_dashboard_data.R). When a piece is missing they render an
# explanatory empty state instead of erroring — so the panel ships before the
# underlying analyses finish recomputing.

# Pretty country label: the comparison CSV stores config keys ("SierraLeone");
# space out the camel-case for display without depending on the metadata map.
.res_pretty_country <- function(x) {
  gsub("([a-z])([A-Z])", "\\1 \\2", x)
}

.res_empty <- function(msg) {
  plotly::plotly_empty(type = "scatter", mode = "markers") |>
    plotly::layout(annotations = list(text = msg, showarrow = FALSE,
                                       font = list(size = 14)))
}

mod_resolution_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "Admin-2 vs cluster",
      icon = bsicons::bs_icon("bullseye"),
      div(
        p("Each survey was modeled at two spatial resolutions from the same ",
          "proxy predictors: administrative-2 areas (district / chiefdom / TA) ",
          "and the finer survey-cluster GPS buffers (2 km). Points above the ",
          "diagonal favor cluster resolution; points below favor Admin-2. The ",
          "spread shows that finer resolution is a sample-size / stability ",
          "trade-off, not a uniform accuracy gain."),
        layout_columns(
          col_widths = c(5, 7),
          radioButtons(ns("metric"), "Metric",
                       choices = c("Predictive skill (Pearson r)" = "r",
                                   "Error (MAE, percentage points)" = "mae"),
                       selected = "r"),
          div()
        ),
        plotlyOutput(ns("scatter"), height = "460px"),
        reactable::reactableOutput(ns("table")),
        methods_note(
          "Engine: a compact area-level SuperLearner (glmnet / ranger / ",
          "xgboost) on logit-transformed prevalence, leave-one-out CV for ",
          "n ≤ 50 areas else 5-fold, run identically at both resolutions. ",
          "Predictors are proxy-only (remote-sensing / geospatial ", tags$code("gee_*"),
          " covariates); survey variables are excluded. ",
          tags$br(), tags$br(),
          "Admin-2 covariates are zonal means over the district polygon; ",
          "cluster covariates are means over a 2 km buffer around each survey ",
          "cluster's GPS point. ", tags$strong("r"), " is the Pearson ",
          "correlation between cross-validated predicted and survey-observed ",
          "prevalence; ", tags$strong("MAE"), " is the mean absolute error in ",
          "percentage points. Δ columns are cluster − admin-2 (positive ",
          "Δr = cluster more skillful; positive ΔMAE = cluster less accurate). ",
          "Blank cluster cells indicate too few labeled clusters to fit a stable ",
          "model (e.g. rare outcomes)."
        )
      )
    ),

    nav_panel(
      title = "Pooled transportability (LOCO)",
      icon = bsicons::bs_icon("compass"),
      div(
        p("Leave-one-country-out: train on three countries, predict the held-out ",
          "fourth, at each resolution and across area-level transport methods ",
          "(baseline / penalized / spatial+soil / CORAL). Bars are the mean ",
          "out-of-country Pearson r across outcomes and held-out countries."),
        plotlyOutput(ns("loco_plot"), height = "440px"),
        reactable::reactableOutput(ns("loco_table")),
        methods_note(
          "Pooled LOCO uses the same area-level comparators at both ",
          "resolutions. A higher bar means the model transports better to an ",
          "unseen country at that resolution. This view populates once ",
          tags$code("results/cluster_vs_admin2_LOCO.csv"), " has been produced."
        )
      )
    )
  )
}


mod_resolution_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cmp <- cluster_res$comparison
    loco <- cluster_res$loco

    # ── Within-country scatter ───────────────────────────────────────────
    output$scatter <- renderPlotly({
      if (is.null(cmp) || nrow(cmp) == 0) {
        return(.res_empty(paste(
          "Cluster-resolution comparison not built yet.",
          "Run the analysis, then re-run 01_prepare_dashboard_data.R.")))
      }
      metric <- input$metric %||% "r"
      ax <- if (metric == "r") "r_admin2"   else "mae_admin2"
      ay <- if (metric == "r") "r_cluster"  else "mae_cluster"

      d <- cmp[is.finite(cmp[[ax]]) & is.finite(cmp[[ay]]), , drop = FALSE]
      if (nrow(d) == 0) return(.res_empty("No paired values for this metric."))

      d$country_lab <- .res_pretty_country(d$country)
      d$oc_lab <- unname(meta$outcome_short[d$outcome])
      d$oc_lab[is.na(d$oc_lab)] <- d$outcome[is.na(d$oc_lab)]

      rng <- range(c(d[[ax]], d[[ay]]), na.rm = TRUE)
      pad <- diff(rng) * 0.08
      lim <- c(rng[1] - pad, rng[2] + pad)
      axis_title <- if (metric == "r") "Pearson r" else "MAE (pp)"

      plot_ly(d) |>
        add_markers(
          x = ~get(ax), y = ~get(ay), color = ~country_lab,
          size = I(110), hoverinfo = "text",
          text = ~sprintf("%s — %s<br>Admin-2: %.3f<br>Cluster: %.3f",
                          country_lab, oc_lab, get(ax), get(ay))) |>
        add_segments(x = lim[1], xend = lim[2], y = lim[1], yend = lim[2],
                     line = list(dash = "dot", color = "#999"),
                     showlegend = FALSE, hoverinfo = "skip", inherit = FALSE) |>
        layout(
          xaxis = list(title = paste0("Admin-2 ", axis_title), range = lim),
          yaxis = list(title = paste0("Cluster ", axis_title), range = lim),
          legend = list(title = list(text = "Country")),
          margin = list(l = 70, b = 60))
    })

    # ── Within-country table ─────────────────────────────────────────────
    output$table <- reactable::renderReactable({
      validate(need(!is.null(cmp) && nrow(cmp) > 0,
                    "Comparison table not built yet."))
      tab <- data.frame(
        Country   = .res_pretty_country(cmp$country),
        Outcome   = ifelse(is.na(unname(meta$outcome_short[cmp$outcome])),
                           cmp$outcome, unname(meta$outcome_short[cmp$outcome])),
        n_admin2  = cmp$n_admin2,  r_admin2 = cmp$r_admin2,  mae_admin2 = cmp$mae_admin2,
        n_cluster = cmp$n_cluster, r_cluster = cmp$r_cluster, mae_cluster = cmp$mae_cluster,
        d_r       = cmp$delta_r,   d_mae    = cmp$delta_mae,
        stringsAsFactors = FALSE, check.names = FALSE)

      num1 <- reactable::colFormat(digits = 3)
      num2 <- reactable::colFormat(digits = 2)
      reactable::reactable(
        tab, compact = TRUE, striped = TRUE, defaultPageSize = 16,
        defaultColDef = reactable::colDef(align = "center"),
        columns = list(
          Country    = reactable::colDef(minWidth = 110, align = "left"),
          Outcome    = reactable::colDef(minWidth = 110, align = "left"),
          r_admin2   = reactable::colDef(name = "r (A2)",  format = num1),
          mae_admin2 = reactable::colDef(name = "MAE (A2)", format = num2),
          r_cluster  = reactable::colDef(name = "r (clus)", format = num1),
          mae_cluster= reactable::colDef(name = "MAE (clus)", format = num2),
          d_r = reactable::colDef(
            name = "Δr", format = num1,
            style = function(v) if (!is.na(v) && v > 0)
              list(color = "#1a7d3c") else list(color = "#b3261e")),
          d_mae = reactable::colDef(
            name = "ΔMAE", format = num2,
            style = function(v) if (!is.na(v) && v <= 0)
              list(color = "#1a7d3c") else list(color = "#b3261e"))
        ))
    })

    # ── Pooled LOCO plot ─────────────────────────────────────────────────
    output$loco_plot <- renderPlotly({
      if (is.null(loco) || nrow(loco) == 0) {
        return(.res_empty(paste(
          "Pooled cluster-vs-admin2 LOCO not built yet.",
          "Produce results/cluster_vs_admin2_LOCO.csv, then re-run data prep.")))
      }
      agg <- aggregate(pearson_r ~ method + resolution, data = loco,
                       FUN = function(x) mean(x, na.rm = TRUE))
      plot_ly(agg, x = ~method, y = ~pearson_r, color = ~resolution,
              type = "bar",
              hoverinfo = "text",
              text = ~sprintf("%s · %s<br>mean r: %.3f",
                              method, resolution, pearson_r)) |>
        layout(barmode = "group",
               xaxis = list(title = "Transport method"),
               yaxis = list(title = "Mean out-of-country Pearson r"),
               legend = list(title = list(text = "Resolution")))
    })

    output$loco_table <- reactable::renderReactable({
      validate(need(!is.null(loco) && nrow(loco) > 0,
                    "LOCO table not built yet."))
      reactable::reactable(loco, compact = TRUE, striped = TRUE,
                           defaultPageSize = 12, filterable = TRUE)
    })

  })
}
