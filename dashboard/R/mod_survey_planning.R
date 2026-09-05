# =============================================================================
# Module: How much survey do you need?
# =============================================================================
# Replaces the old "Scenarios" tab.
#
# WHY SCENARIOS WENT
# ------------------
# Scenarios let a user pick a coverage and an effect size and watch cases-averted
# change. The effect sizes were literature defaults, not estimated here, so the
# tab's output was a restatement of the user's own assumptions with this
# project's district map used as decoration. Its "what-if explorer" mode was
# worse: it scaled prevalence by an arbitrary shift, which reads as a forecast.
# Neither answered a question this project has evidence about.
#
# WHAT THIS ANSWERS INSTEAD
# -------------------------
# The budget question a survey planner actually faces: how many regions to
# visit, and how many clusters to field in each. Every number here is measured
# (scripts/accuracy_impact/ws5_anchoring_budget.R, 40 replicates over a 9 x 6
# design grid), not assumed.
#
# THE FINDING THE TAB IS BUILT AROUND
# -----------------------------------
# District-level accuracy is almost flat in survey size: 10.6 pp mean absolute
# error at 5% of a full survey against 9.2 pp at 100%. Twenty times the money
# buys 1.4 pp. Region-level accuracy, by contrast, responds steeply: 6.8 pp at
# 5% falling to under 3 pp past half a survey.
#
# That is not a defect in the design curve, it is the point. Assigning a
# district its region's mean cannot distinguish districts WITHIN a region, no
# matter how precisely that mean is estimated. So survey money buys regional
# precision, and district resolution has to come from somewhere else - which is
# the case for the covariate model, made in the units a budget holder uses.
#
# HONESTY REQUIREMENTS
#   - mae_a1 = 0 at the full survey is definitional (the estimate is compared
#     with itself), not a claim of perfect accuracy. Flagged in the UI.
#   - The regional mean is JACKKNIFED. The un-jackknifed anchoring result
#     elsewhere in this project was withdrawn for that leak.
#   - This is the region-mean estimator only. It is the baseline a survey buys,
#     not the model's performance.

mod_survey_planning_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "What a survey buys",
      icon = bsicons::bs_icon("cash-coin"),
      div(
        p(class = "lead",
          "How accurate are prevalence estimates if you field only part of a ",
          "full biomarker survey? Every point below is measured, not modelled ",
          "from an assumption."),

        layout_columns(
          col_widths = c(6, 6),
          value_box(
            title = "District estimates barely improve with budget",
            value = textOutput(ns("vb_a2")),
            showcase = bsicons::bs_icon("pin-map"),
            theme = "warning",
            p(class = "small mb-0",
              "Mean absolute error at district level, smallest survey vs. full ",
              "survey. A region's mean cannot tell its districts apart.")),
          value_box(
            title = "Regional estimates improve steeply",
            value = textOutput(ns("vb_a1")),
            showcase = bsicons::bs_icon("bullseye"),
            theme = "success",
            p(class = "small mb-0",
              "Mean absolute error at region level, from the smallest design ",
              "to the largest one short of a full census of clusters. This is ",
              "what survey money actually buys."))
        ),

        tags$br(),
        plotlyOutput(ns("curve"), height = "420px"),
        p(class = "text-muted mt-2",
          "Each point is one design on the grid, pooled across countries and ",
          "outcomes (median over 40 replicates). Horizontal axis is the share ",
          "of a full survey's clusters actually fielded. The full-survey point ",
          "is omitted from the region series, where it is a comparison of the ",
          "estimate with itself."),

        div(class = "alert alert-info border mt-3",
          tags$strong("How to read this. "),
          "Both curves describe the same estimator: give every district its ",
          "region's survey mean. The regional curve falls steeply because more ",
          "survey means a better-estimated regional mean. The district curve is ",
          "nearly flat because the limit is not precision, it is ",
          tags$em("resolution"),
          " — every district in a region gets the same number however well ",
          "that number is measured. Buying district resolution needs ",
          "district-varying information, which is the case for the covariate ",
          "model, not for a bigger survey.")
      )
    ),

    nav_panel(
      title = "Design a survey",
      icon = bsicons::bs_icon("sliders"),
      div(
        p("Pick a design and read off what it costs and what it delivers. ",
          "Cost is expressed as a share of a full survey, because that is the ",
          "unit these estimates are calibrated in."),
        layout_columns(
          col_widths = c(6, 6),
          sliderInput(ns("regions"), "Share of regions visited",
                      min = 0.2, max = 1, value = 0.5, step = 0.1),
          selectInput(ns("clusters"), "Share of each region's clusters fielded",
                      choices = c("15%" = 0.15, "25%" = 0.25, "40%" = 0.4,
                                  "60%" = 0.6, "80%" = 0.8, "100%" = 1.0),
                      selected = 0.4)
        ),
        uiOutput(ns("design_readout")),
        tags$br(),
        h6("The whole grid"),
        reactable::reactableOutput(ns("grid")),
        p(class = "text-muted mt-2",
          "Sorted by cost. 'Cost' is the median share of the full survey's ",
          "sample actually used by that design.")
      )
    ),

    nav_panel(
      title = "Read this first",
      icon = bsicons::bs_icon("exclamation-triangle"),
      div(
        h5("What these numbers are, and are not"),
        tags$ul(
          tags$li(tags$strong("This is the survey-only baseline, not the model. "),
            "Every figure describes one estimator: assign each district its ",
            "region's design-based survey mean. It tells you what a survey of a ",
            "given size delivers on its own. The covariate model's performance ",
            "is on the Benchmarks tab."),
          tags$li(tags$strong("The regional mean is jackknifed. "),
            "A district's regional anchor is computed from the region's ",
            tags$em("other"), " districts only. This matters: the headline ",
            "anchoring result reported earlier in this project was withdrawn ",
            "because its anchor included the scored district's own respondents. ",
            "This analysis was written afterwards and does not repeat that."),
          tags$li(tags$strong("Zero regional error at 100% is definitional. "),
            "At a full survey the regional estimate is being compared with ",
            "itself, so its error is 0 by construction, not by merit. Read the ",
            "regional curve for its shape, not its endpoint."),
          tags$li(tags$strong("Whole clusters are dropped, not individuals. "),
            "A survey planner buys clusters. Retaining a fraction of people ",
            "within every cluster would cost the same as the full survey."),
          tags$li(tags$strong("Pooled across countries and outcomes. "),
            "A specific country with few regions will do worse than the pooled ",
            "curve suggests. Country region counts range from 6 to 27.")
        ),
        h5("Provenance"),
        uiOutput(ns("prov"))
      )
    )
  )
}

mod_survey_planning_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    bundle <- reactive({
      p <- file.path("data", "survey_planning.rds")
      if (!file.exists(p)) return(NULL)
      readRDS(p)
    })

    # Smallest and largest designs on the grid, for the value boxes.
    #
    # The region-level endpoint needs care. At region_share = 1 AND
    # fraction_clusters = 1 the regional estimate is compared with itself, so
    # its error is 0 by construction. Quoting "6.8 -> 0.0 pp" in a value box
    # would advertise a definitional artefact as a result, so the region box
    # reports the largest design that is NOT the full survey. The district box
    # keeps the full survey, where 9.2 pp is a real residual error: even with
    # the whole survey in hand, giving every district its region mean is 9.2 pp
    # out, which is the entire point of the tab.
    ends <- reactive({
      b <- bundle(); req(b)
      p <- b$pooled
      full <- abs(p$region_share - 1) < 1e-6 & abs(p$fraction_clusters - 1) < 1e-6
      sub <- p[!full, ]
      list(lo = p[which.min(p$pct_survey), ],
           hi = p[which.max(p$pct_survey), ],
           hi_sub = sub[which.max(sub$pct_survey), ])
    })

    output$vb_a2 <- renderText({
      e <- ends(); req(e)
      sprintf("%.1f → %.1f pp", e$lo$mae_a2, e$hi$mae_a2)
    })
    output$vb_a1 <- renderText({
      e <- ends(); req(e)
      sprintf("%.1f → %.1f pp", e$lo$mae_a1, e$hi_sub$mae_a1)
    })

    output$curve <- renderPlotly({
      b <- bundle(); req(b)
      p <- b$pooled[order(b$pooled$pct_survey), ]
      # The full-survey point is a self-comparison for the REGION series only
      # (error 0 by construction), so it is dropped from that trace. It is a
      # genuine measurement for the district series and stays there.
      pr <- p[!(abs(p$region_share - 1) < 1e-6 &
                  abs(p$fraction_clusters - 1) < 1e-6), ]
      plot_ly() |>
        add_markers(data = p, x = ~pct_survey, y = ~mae_a2,
                    name = "District (admin-2)",
                    marker = list(size = 9, color = "#d9822b"),
                    hovertemplate = paste0(
                      "%{y:.1f} pp district error<br>",
                      "%{x:.0f}% of a full survey<extra></extra>")) |>
        add_markers(data = pr, x = ~pct_survey, y = ~mae_a1,
                    name = "Region (admin-1)",
                    marker = list(size = 9, color = "#2c7fb8"),
                    hovertemplate = paste0(
                      "%{y:.1f} pp region error<br>",
                      "%{x:.0f}% of a full survey<extra></extra>")) |>
        layout(
          xaxis = list(title = "Share of a full survey fielded (%)"),
          yaxis = list(title = "Mean absolute error (percentage points)",
                       rangemode = "tozero"),
          legend = list(orientation = "h", x = 0, y = 1.12),
          hovermode = "closest")
    })

    picked <- reactive({
      b <- bundle(); req(b)
      p <- b$pooled
      rs <- round(as.numeric(input$regions), 1)
      fc <- as.numeric(input$clusters)
      r <- p[abs(p$region_share - rs) < 1e-6 &
               abs(p$fraction_clusters - fc) < 1e-6, ]
      if (!nrow(r)) NULL else r[1, ]
    })

    output$design_readout <- renderUI({
      b <- bundle(); req(b)
      r <- picked()
      if (is.null(r))
        return(div(class = "alert alert-secondary",
                   "That combination is not on the measured grid. Pick another."))
      full <- b$baseline_mae
      div(class = "alert alert-light border",
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          div(h6(class = "text-muted mb-1", "Cost"),
              h4(sprintf("%.0f%%", r$pct_survey)),
              tags$small(class = "text-muted", "of a full survey")),
          div(h6(class = "text-muted mb-1", "District error"),
              h4(sprintf("%.1f pp", r$mae_a2)),
              tags$small(class = "text-muted",
                         sprintf("full survey: %.1f pp", full))),
          div(h6(class = "text-muted mb-1", "Region error"),
              h4(sprintf("%.1f pp", r$mae_a1)),
              tags$small(class = "text-muted", "vs. full-survey regional value")),
          div(h6(class = "text-muted mb-1", "District bias"),
              h4(sprintf("%+.1f pp", r$bias_a2)),
              tags$small(class = "text-muted", "negative = under-estimates"))
        ),
        tags$hr(class = "my-2"),
        tags$small(sprintf(
          paste("This design uses %.0f%% of the survey and lands within %.1f pp",
                "of what the full survey's district estimates give you (%.1f pp).",
                "Spending the remaining %.0f%% of the budget buys %.1f pp of",
                "district accuracy."),
          r$pct_survey, abs(r$mae_a2 - full), full,
          100 - r$pct_survey, max(r$mae_a2 - full, 0))))
    })

    output$grid <- reactable::renderReactable({
      b <- bundle(); req(b)
      d <- b$pooled[order(b$pooled$pct_survey), ]
      d$regions  <- paste0(round(100 * d$region_share), "%")
      d$clusters <- paste0(round(100 * d$fraction_clusters), "%")
      reactable::reactable(
        d[, c("regions", "clusters", "pct_survey", "mae_a2", "mae_a1", "bias_a2")],
        columns = list(
          regions    = reactable::colDef(name = "Regions visited", width = 130,
                                         align = "center"),
          clusters   = reactable::colDef(name = "Clusters fielded", width = 130,
                                         align = "center"),
          pct_survey = reactable::colDef(name = "Cost (% of survey)", width = 150,
                         align = "center",
                         format = reactable::colFormat(digits = 1)),
          mae_a2     = reactable::colDef(name = "District error (pp)",
                         format = reactable::colFormat(digits = 2)),
          mae_a1     = reactable::colDef(name = "Region error (pp)",
                         format = reactable::colFormat(digits = 2)),
          bias_a2    = reactable::colDef(name = "District bias (pp)",
                         format = reactable::colFormat(digits = 2))),
        striped = TRUE, highlight = TRUE, defaultPageSize = 12,
        searchable = FALSE)
    })

    output$prov <- renderUI({
      b <- bundle(); req(b)
      div(class = "small text-muted",
        p(sprintf("Source: %s (%d replicates per design).", b$source, b$n_reps)),
        p(sprintf("Built %s.", format(b$build_time, "%Y-%m-%d %H:%M"))),
        p(b$note))
    })
  })
}
