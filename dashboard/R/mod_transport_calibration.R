# =============================================================================
# Module: Transport calibration
# =============================================================================
# Interactive predicted-vs-true prevalence calibration under leave-one-country-
# out (LOCO) holdout. The user selects across three axes:
#   (1) Outcome              — micronutrient × population
#   (2) Aggregation level    — National / Admin-1 / Admin-2
#   (3) Model / approach     — depends on level (individual vs area SL; loss)
#
# A point is one held-out unit (country at national level; district/region at
# sub-national levels). Distance from the 45-degree line is the calibration
# error; the fitted slope (<1 = compression toward the training mean, <0 =
# inverted ordering) summarises transportability. Reads the optional
# `transport_cal` object built by 01_prepare_dashboard_data.R; renders an
# empty state when a source is missing.

# Approaches available at each aggregation level (label -> internal key).
.tc_approaches <- function(level) {
  if (level == "National") {
    c("Individual-level SuperLearner" = "indiv",
      "Area-level SL (continuous loss)" = "area_continuous",
      "Area-level SL (logit loss)"      = "area_logit")
  } else {
    c("Area-level SuperLearner" = "area")
  }
}

.tc_pretty_country <- function(x) gsub("([a-z])([A-Z])", "\\1 \\2", as.character(x))

.tc_empty <- function(msg) {
  plotly::plotly_empty(type = "scatter", mode = "markers") |>
    plotly::layout(annotations = list(text = msg, showarrow = FALSE,
                                      font = list(size = 14)))
}

mod_transport_calibration_ui <- function(id) {
  ns <- NS(id)
  div(
    p("Predicted versus true deficiency prevalence when the model is ",
      "transported to a held-out country (leave-one-country-out). Points on ",
      "the dashed 45° line are perfectly calibrated; a shallow or negative ",
      "fitted slope means predictions collapse toward the cross-country ",
      "average. Choose an outcome, the spatial level at which prevalence is ",
      "compared, and the modelling approach."),
    layout_columns(
      col_widths = c(4, 4, 4),
      selectInput(ns("outcome"), "Outcome", choices = NULL),
      radioButtons(ns("level"), "Aggregation level",
                   choices = c("National", "Admin-1", "Admin-2"),
                   selected = "National", inline = TRUE),
      selectInput(ns("approach"), "Model / approach", choices = NULL)
    ),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box(title = "Calibration slope", value = textOutput(ns("vb_slope")),
                theme = "primary", showcase = bsicons::bs_icon("rulers")),
      value_box(title = "Mean abs. error", value = textOutput(ns("vb_mae")),
                theme = "secondary", showcase = bsicons::bs_icon("dash-circle")),
      value_box(title = "Held-out units", value = textOutput(ns("vb_n")),
                theme = "secondary", showcase = bsicons::bs_icon("geo-alt")),
      value_box(title = "Perfect = 1.00", value = "reference",
                theme = "light", showcase = bsicons::bs_icon("bullseye"))
    ),
    plotlyOutput(ns("scatter"), height = "460px"),
    methods_note(
      "National points use the pooled ", tags$b("individual-level"),
      " SuperLearner (one point per held-out country); admin-1 and admin-2 use ",
      "the ", tags$b("area-level"), " SuperLearner (one point per district / ",
      "region). Admin-1 is the population-weighted aggregate of the admin-2 ",
      "predictions. Predictors are proxy-only. The fitted slope is the OLS ",
      "slope of predicted on true prevalence across the held-out units; MAE is ",
      "the mean absolute calibration gap in percentage points. ",
      tags$br(), tags$br(),
      "Folate and B12 are available at the national level only (women, three ",
      "countries); zinc is Malawi-only and cannot be transported. Additional ",
      "transport methods (penalized, spatial+soil, CORAL) will populate the ",
      "approach selector as those results are finalised."
    )
  )
}


mod_transport_calibration_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    tc <- transport_cal

    # Outcomes available for the current (level, approach) selection.
    available_outcomes <- reactive({
      lvl <- input$level %||% "National"
      app <- input$approach %||% "indiv"
      ocs <- if (lvl == "National" && app == "indiv") {
        if (!is.null(tc$nat_indiv)) unique(tc$nat_indiv$outcome) else character(0)
      } else if (lvl == "National") {
        if (!is.null(tc$nat_area)) unique(tc$nat_area$outcome) else character(0)
      } else {
        if (!is.null(tc$adm2_area)) unique(tc$adm2_area$outcome) else character(0)
      }
      # Order by the canonical outcome list, label nicely.
      ord <- intersect(names(meta$outcome_labels), ocs)
      setNames(ord, meta$outcome_labels[ord])
    })

    # Keep the approach selector in sync with the level.
    observeEvent(input$level, {
      apps <- .tc_approaches(input$level)
      updateSelectInput(session, "approach", choices = apps,
                        selected = unname(apps)[1])
    }, ignoreNULL = TRUE)

    # Keep the outcome selector in sync with (level, approach).
    observe({
      ocs <- available_outcomes()
      cur <- isolate(input$outcome)
      sel <- if (!is.null(cur) && cur %in% ocs) cur else unname(ocs)[1]
      updateSelectInput(session, "outcome", choices = ocs, selected = sel)
    })

    # Assemble the calibration points for the current selection.
    plot_data <- reactive({
      lvl <- input$level %||% "National"
      app <- input$approach %||% "indiv"
      oc  <- input$outcome
      req(oc)

      if (lvl == "National" && app == "indiv") {
        d <- tc$nat_indiv
        if (is.null(d)) return(NULL)
        d <- d[d$outcome == oc, , drop = FALSE]
        if (!nrow(d)) return(NULL)
        data.frame(x = d$true_prev * 100, y = d$pred_prev * 100,
                   group = .tc_pretty_country(d$held_out),
                   label = .tc_pretty_country(d$held_out),
                   stringsAsFactors = FALSE)

      } else if (lvl == "National") {
        mt <- sub("^area_", "", app)  # continuous / logit
        d <- tc$nat_area
        if (is.null(d)) return(NULL)
        d <- d[d$outcome == oc & d$model_type == mt, , drop = FALSE]
        if (!nrow(d)) return(NULL)
        data.frame(x = d$true_prev * 100, y = d$pred_prev * 100,
                   group = .tc_pretty_country(d$held_out),
                   label = .tc_pretty_country(d$held_out),
                   stringsAsFactors = FALSE)

      } else {
        d <- tc$adm2_area
        if (is.null(d)) return(NULL)
        d <- d[d$outcome == oc, , drop = FALSE]
        if (!nrow(d)) return(NULL)
        d$country_label <- unname(meta$countries[d$country])
        if (lvl == "Admin-2") {
          return(data.frame(x = d$true_prev * 100, y = d$pred_prev * 100,
                            group = d$country_label, label = d$Admin2,
                            stringsAsFactors = FALSE))
        }
        # Admin-1: population-weighted aggregate of admin-2 predictions.
        a1 <- unique(admin2_pred[, c("country", "Admin1", "Admin2")])
        pop <- admin2_pop[, c("country", "Admin2", "population")]
        m <- merge(d, a1, by.x = c("country_label", "Admin2"),
                   by.y = c("country", "Admin2"), all.x = TRUE)
        m <- merge(m, pop, by.x = c("country_label", "Admin2"),
                   by.y = c("country", "Admin2"), all.x = TRUE)
        m$w <- ifelse(is.na(m$population) | m$population <= 0, m$n_svy, m$population)
        m$w[is.na(m$w) | m$w <= 0] <- 1
        m <- m[!is.na(m$Admin1), , drop = FALSE]
        if (!nrow(m)) return(NULL)
        wm <- function(v, w) sum(v * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
        agg <- do.call(rbind, lapply(split(m, list(m$country_label, m$Admin1), drop = TRUE),
          function(g) data.frame(
            x = wm(g$true_prev, g$w) * 100, y = wm(g$pred_prev, g$w) * 100,
            group = g$country_label[1], label = g$Admin1[1],
            stringsAsFactors = FALSE)))
        agg
      }
    })

    metrics <- reactive({
      d <- plot_data()
      if (is.null(d) || nrow(d) < 1) return(list(slope = NA, mae = NA, n = 0))
      slope <- if (nrow(d) >= 2 && stats::sd(d$x) > 0)
                 unname(stats::coef(stats::lm(y ~ x, data = d))[2]) else NA_real_
      list(slope = slope, mae = mean(abs(d$y - d$x), na.rm = TRUE), n = nrow(d))
    })

    output$vb_slope <- renderText({
      s <- metrics()$slope; if (is.na(s)) "—" else sprintf("%.2f", s) })
    output$vb_mae   <- renderText({
      m <- metrics()$mae; if (is.na(m)) "—" else sprintf("%.1f pp", m) })
    output$vb_n     <- renderText({ as.character(metrics()$n) })

    output$scatter <- renderPlotly({
      d <- plot_data()
      if (is.null(d) || nrow(d) == 0)
        return(.tc_empty("No transport results for this selection yet."))
      rng <- range(c(d$x, d$y), na.rm = TRUE); pad <- max(diff(rng) * 0.08, 1)
      lim <- c(max(0, rng[1] - pad), rng[2] + pad)
      many <- nrow(d) > 12

      p <- plot_ly(d) |>
        add_markers(x = ~x, y = ~y, color = ~group,
                    size = I(if (many) 60 else 110), hoverinfo = "text",
                    text = ~sprintf("%s<br>true: %.1f%%<br>pred: %.1f%%", label, x, y)) |>
        add_segments(x = lim[1], xend = lim[2], y = lim[1], yend = lim[2],
                     line = list(dash = "dot", color = "#999"),
                     showlegend = FALSE, hoverinfo = "skip", inherit = FALSE)
      if (!many)  # label the few held-out countries directly
        p <- p |> add_text(x = ~x, y = ~y, text = ~label, textposition = "top center",
                           showlegend = FALSE, hoverinfo = "skip", inherit = FALSE,
                           textfont = list(size = 10, color = "#444"))
      p |> layout(
        xaxis = list(title = "True prevalence (%)", range = lim),
        yaxis = list(title = "Predicted prevalence (%)", range = lim),
        legend = list(title = list(text = if (many) "Country" else "Held-out")))
    })
  })
}
