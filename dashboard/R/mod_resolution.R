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
#   (3) Anchoring — where a district map's LEVEL comes from. The proxy
#       covariates carry the spatial PATTERN but not the level, so production
#       maps re-scale predictions to the design-based survey total within each
#       Admin-1 region. This view scores that choice against not anchoring, and
#       against anchoring to the (too coarse) national total.
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

# Arm labels for the anchoring panel. The CSV's arm strings are analyst-facing
# ("admin-2 fit + ADMIN-1 benchmark (hard)"); the dashboard names them for what
# they DO to the map, and the order below is the order of the argument: no
# anchor, an anchor too coarse to help, the two that work, then the alternative
# of not modeling districts at all.
.RES_ARM_ORDER <- c("No anchor", "National anchor", "Admin-1 anchor (shrunk)",
                    "Admin-1 anchor (hard)", "Fit at Admin-1")

.res_arm_label <- function(x) {
  x <- trimws(gsub('"', "", x))
  out <- rep(NA_character_, length(x))
  out[grepl("unbenchmarked", x, fixed = TRUE)]        <- "No anchor"
  out[grepl("national benchmark", x, fixed = TRUE)]   <- "National anchor"
  out[grepl("(shrunk)", x, fixed = TRUE)]             <- "Admin-1 anchor (shrunk)"
  out[grepl("(hard)", x, fixed = TRUE)]               <- "Admin-1 anchor (hard)"
  out[grepl("^ADMIN-1 fit", x)]                       <- "Fit at Admin-1"
  ifelse(is.na(out), x, out)
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
    ),

    nav_panel(
      title = "Where the level comes from",
      icon = bsicons::bs_icon("align-bottom"),
      div(
        p(tags$strong("Two different things come out of a district model, and they are not equally trustworthy."),
          " The ", tags$em("pattern"), " — which districts are worse than ",
          "which — is what the proxy covariates carry. The ", tags$em("level"),
          " — how high prevalence actually sits — they carry badly. ",
          "Every district map elsewhere in this dashboard therefore takes its ",
          "level from the survey and only its pattern from the model. This ",
          "panel is the evidence for that choice."),
        p("Each arm below re-scales the same district predictions so they ",
          "average, within each region, to the design-based survey estimate for ",
          "that region — the ", tags$em("anchor"), ". Anchoring to Admin-1 ",
          "more than doubles rank agreement and cuts absolute bias by more than ",
          "half. Anchoring to the national total instead adds nothing to skill ",
          "and makes bias worse, because one number cannot correct a level that ",
          "is wrong region by region — it moves every district by the same ",
          "amount, including the ones that were already right."),
        layout_columns(
          col_widths = c(5, 7),
          radioButtons(ns("anchor_metric"), "Metric",
                       choices = c("Predictive skill (Pearson r)" = "pearson_r",
                                   "Bias (percentage points)" = "bias_pp",
                                   "Error (MAE, percentage points)" = "mae_pp"),
                       selected = "pearson_r"),
          div()
        ),
        plotlyOutput(ns("anchor_plot"), height = "400px"),
        tags$br(),
        reactable::reactableOutput(ns("anchor_table")),
        methods_note(
          "Arms are scored on the same country × outcome cells with the same ",
          "leave-one-region-out folds, so differences are attributable to the ",
          "anchor and not to the fit. ", tags$strong("Unbenchmarked"),
          " is the raw model output. ", tags$strong("Hard"),
          " forces each region's predicted mean onto the survey estimate; ",
          tags$strong("shrunk"), " moves it part of the way, in proportion to ",
          "how precisely that region was surveyed. ", tags$strong("ADMIN-1 fit"),
          " skips districts entirely and fits at the region. ",
          tags$code("r_share"), " is r as a fraction of the noise ceiling ",
          tags$code("r_max"), " — the most any model could reach given how ",
          "noisily the survey measures each district.",
          tags$br(), tags$br(),
          "The practical consequence: a district map for a country that HAS a ",
          "survey can be read as prevalence, because the survey supplies the ",
          "level. A transported map for a country with no survey of its own has ",
          "no anchor, and should be read as a ranking only."
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
    anchor <- anchoring$arms

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

    # ── Anchoring: where the level comes from ─────────────────
    # One bar per arm, one point per country × outcome cell behind it. The
    # points matter as much as the bars: the anchor helps in 20 of 24 cells,
    # and a bar alone would hide that the remaining four exist.
    output$anchor_plot <- renderPlotly({
      if (is.null(anchor) || nrow(anchor) == 0) {
        return(.res_empty(paste(
          "Anchoring comparison not built yet.",
          "Produce results/tables/admin1_arms.csv, then re-run data prep.")))
      }
      metric <- input$anchor_metric %||% "pearson_r"
      d <- anchor
      d$value <- suppressWarnings(as.numeric(d[[metric]]))
      # Bias is signed; what matters for a level claim is how far off it is.
      if (metric == "bias_pp") d$value <- abs(d$value)
      d <- d[is.finite(d$value), , drop = FALSE]
      if (nrow(d) == 0) return(.res_empty("No values for this metric."))

      d$arm_lab <- .res_arm_label(d$arm)
      d$arm_lab <- factor(d$arm_lab, levels = .RES_ARM_ORDER)
      d <- d[!is.na(d$arm_lab), , drop = FALSE]
      d <- d[order(d$arm_lab), , drop = FALSE]

      agg <- stats::aggregate(value ~ arm_lab, data = d,
                              FUN = function(x) mean(x, na.rm = TRUE))
      ttl <- switch(metric,
                    pearson_r = "Pearson r (higher is better)",
                    bias_pp   = "Absolute bias, pp (lower is better)",
                    mae_pp    = "MAE, pp (lower is better)")

      d$oc_lab <- unname(meta$outcome_short[d$outcome])
      d$oc_lab[is.na(d$oc_lab)] <- d$outcome[is.na(d$oc_lab)]

      plot_ly() |>
        add_bars(data = agg, x = ~arm_lab, y = ~value,
                 marker = list(color = "#c9d6e8"),
                 hoverinfo = "text",
                 text = ~sprintf("%s<br>mean: %.3f", arm_lab, value),
                 name = "mean") |>
        add_markers(data = d, x = ~arm_lab, y = ~value,
                    marker = list(color = "#1f4e79", size = 7, opacity = 0.65),
                    hoverinfo = "text",
                    text = ~sprintf("%s · %s<br>%s: %.3f",
                                    .res_pretty_country(country), oc_lab,
                                    ttl, value),
                    name = "cell") |>
        layout(showlegend = FALSE,
               xaxis = list(title = "", tickangle = -18),
               yaxis = list(title = ttl))
    })

    output$anchor_table <- reactable::renderReactable({
      validate(need(!is.null(anchor) && nrow(anchor) > 0,
                    "Anchoring table not built yet."))
      num2 <- reactable::colFormat(digits = 2)
      num3 <- reactable::colFormat(digits = 3)
      d <- anchor
      d$country <- .res_pretty_country(d$country)
      d$arm <- .res_arm_label(d$arm)
      reactable::reactable(
        d, compact = TRUE, striped = TRUE, filterable = TRUE,
        defaultPageSize = 12, groupBy = "country",
        columns = list(
          country  = reactable::colDef(name = "Country", minWidth = 110),
          outcome  = reactable::colDef(name = "Outcome", minWidth = 110),
          arm      = reactable::colDef(name = "Arm", minWidth = 190),
          n_admin2 = reactable::colDef(name = "Districts", aggregate = "max"),
          pearson_r = reactable::colDef(name = "r", format = num3),
          r_max    = reactable::colDef(name = "Ceiling", format = num3),
          r_share  = reactable::colDef(name = "% of ceiling", format = num2),
          mae_pp   = reactable::colDef(name = "MAE (pp)", format = num2),
          bias_pp  = reactable::colDef(name = "Bias (pp)", format = num2)
        ))
    })

  })
}
