# =============================================================================
# Module: Decision value  (stakeholder translation)
# =============================================================================
# Answers the question a decision-maker actually has — "IF / WHEN / HOW is this
# model useful to ME?" — instead of MSE/AUC. Three panels, all computed from
# EXISTING dashboard results (no model refit):
#
#   (1) Targeting accuracy   — if you act on the model to pick the highest-burden
#       districts, what share of the truly-deficient population do you reach,
#       vs. having no sub-national data (national average → can only pick by
#       chance)? A burden-concentration ("enrichment") curve + plain headline.
#   (2) Status quo vs model  — per outcome/country, the targeting "lift" and
#       top-quintile hit-rate of model-guided vs. average-based prioritisation.
#   (3) Fitness for purpose  — a traffic-light scorecard (when to trust the model
#       for targeting, by outcome × country, from the OUT-OF-SAMPLE transport
#       diagnostics) plus a plain-language calibration readout
#       ("when the model says 30%, the truth averages X%").
#
# Sources: transport_cal$adm2_area (out-of-sample LOCO predictions vs survey
# truth — the honest "transported to a new place" case), admin2_pred (within-
# country), admin2_pop (burden weights), benchmarks_data$admin2_error,
# model_diag$calibration. Each panel renders an empty state when its source is
# missing.

.dv_empty <- function(msg) {
  plotly::plotly_empty(type = "scatter", mode = "markers") |>
    plotly::layout(annotations = list(text = msg, showarrow = FALSE,
                                      font = list(size = 14)))
}

.dv_pretty_outcome <- function(oc) {
  lab <- meta$outcome_labels[oc]
  ifelse(is.na(lab), oc, lab)
}

# Pull the analysis frame (Admin2, true_prev, pred_prev, country label) for a
# given setting. "loco" = out-of-sample transported predictions (realistic);
# "within" = within-country fitted predictions (optimistic upper bound).
.dv_frame <- function(setting) {
  if (setting == "within") {
    ap <- admin2_pred
    if (is.null(ap) || nrow(ap) == 0) return(NULL)
    data.frame(outcome = ap$outcome,
               country_label = ap$country,
               Admin2 = ap$Admin2,
               true_prev = ap$obs_prev,
               pred_prev = ap$pred_prev,
               stringsAsFactors = FALSE)
  } else {
    d <- transport_cal$adm2_area
    if (is.null(d) || nrow(d) == 0) return(NULL)
    data.frame(outcome = d$outcome,
               country_label = unname(meta$countries[d$country]),
               Admin2 = d$Admin2,
               true_prev = d$true_prev,
               pred_prev = d$pred_prev,
               stringsAsFactors = FALSE)
  }
}

# Attach the subgroup population for each district (child vs women), falling
# back to equal weights when population is unavailable.
.dv_add_weight <- function(df, oc) {
  df$w <- NA_real_
  if (!is.null(admin2_pop) && nrow(admin2_pop) > 0) {
    is_child <- startsWith(oc, "child_")
    pcol <- if (is_child) "pop_child" else "pop_women"
    if (!pcol %in% colnames(admin2_pop)) pcol <- "population"
    pop <- admin2_pop[, c("country", "Admin2", pcol)]
    colnames(pop) <- c("country_label", "Admin2", "w_pop")
    df <- merge(df, pop, by = c("country_label", "Admin2"),
                all.x = TRUE, sort = FALSE)
    df$w <- df$w_pop
  }
  df$w[!is.finite(df$w) | df$w <= 0] <- NA_real_
  if (all(is.na(df$w))) df$w <- 1 else
    df$w[is.na(df$w)] <- stats::median(df$w, na.rm = TRUE)
  df
}

# Burden-concentration curves: cumulative share of true deficient people
# captured (y) as districts are added in priority order (x = cumulative share
# of the subgroup population reached).
.dv_curves <- function(df) {
  df <- df[is.finite(df$true_prev) & is.finite(df$pred_prev), , drop = FALSE]
  if (nrow(df) < 2) return(NULL)
  burden <- df$true_prev * df$w
  totB <- sum(burden, na.rm = TRUE); totW <- sum(df$w, na.rm = TRUE)
  if (!is.finite(totB) || totB <= 0) return(NULL)
  cum <- function(ord) data.frame(
    x = c(0, cumsum(df$w[ord]) / totW),
    y = c(0, cumsum(burden[ord]) / totB))
  list(model   = cum(order(df$pred_prev, decreasing = TRUE)),
       perfect = cum(order(df$true_prev, decreasing = TRUE)),
       n = nrow(df))
}

.dv_capture_at <- function(curve, f) {
  if (is.null(curve)) return(NA_real_)
  stats::approx(curve$x, curve$y, xout = f, ties = "ordered", rule = 2)$y
}

# Top-fraction hit-rate: of the districts the model flags as the top f, how many
# are in the true top f (set overlap).
.dv_hit_rate <- function(df, f) {
  df <- df[is.finite(df$true_prev) & is.finite(df$pred_prev), , drop = FALSE]
  n <- nrow(df); if (n < 2) return(NA_real_)
  k <- max(1L, round(f * n))
  top_m <- order(df$pred_prev, decreasing = TRUE)[seq_len(k)]
  top_t <- order(df$true_prev, decreasing = TRUE)[seq_len(k)]
  length(intersect(top_m, top_t)) / k
}


mod_decision_value_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "Targeting accuracy",
      icon = bsicons::bs_icon("bullseye"),
      div(
        p(strong("The decision: "), "you can only act in some districts first. ",
          "If you rank districts by the model and start with the worst, how much ",
          "of the truly-deficient population do you reach — compared with having ",
          "no sub-national data (where you can only pick districts by chance / ",
          "population size)?"),
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput(ns("ta_outcome"), "Outcome", choices = outcome_choices,
                      selected = "women_iron"),
          selectInput(ns("ta_country"), "Country",
                      choices = c("All countries (pooled)" = "__all__",
                                  country_choices),
                      selected = "__all__"),
          radioButtons(ns("ta_setting"), "Model setting",
                       choices = c("Transported to a new place (realistic)" = "loco",
                                   "Within-country fit (optimistic)" = "within"),
                       selected = "loco")
        ),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          value_box(title = "Reach (act in worst 20% of pop)",
                    value = textOutput(ns("vb_capture")),
                    theme = "primary", showcase = bsicons::bs_icon("people")),
          value_box(title = "Lift vs no targeting",
                    value = textOutput(ns("vb_lift")),
                    theme = "secondary", showcase = bsicons::bs_icon("arrow-up-right")),
          value_box(title = "Top-quintile hit rate",
                    value = textOutput(ns("vb_hit")),
                    theme = "secondary", showcase = bsicons::bs_icon("check2-square")),
          value_box(title = "Districts evaluated",
                    value = textOutput(ns("vb_n")),
                    theme = "light", showcase = bsicons::bs_icon("geo-alt"))
        ),
        plotlyOutput(ns("ta_curve"), height = "440px"),
        methods_note(
          "The curve shows the cumulative share of truly-deficient people ",
          "reached (y) as districts are added worst-first (x = cumulative share ",
          "of the relevant subgroup population). ", tags$b("Model"),
          " orders districts by predicted prevalence; ", tags$b("Perfect"),
          " by the (unknown) true prevalence; ", tags$b("No targeting"),
          " is the 45° line — the best you can do with only a national average, ",
          "i.e. reaching people in proportion to where they live. The gap ",
          "between Model and the diagonal is the decision value the model adds. ",
          tags$br(), tags$br(),
          "“Truth” here is the survey-observed district prevalence (itself ",
          "subject to sampling error). The ", tags$b("realistic"), " setting ",
          "uses leave-one-country-out predictions (the model never saw that ",
          "country) — the honest analogue of deploying to a new area; the ",
          tags$b("optimistic"), " setting uses within-country fitted values and ",
          "is an upper bound."
        )
      )
    ),

    nav_panel(
      title = "Status quo vs model",
      icon = bsicons::bs_icon("arrow-left-right"),
      div(
        p("For every outcome × country, how much better is model-guided ",
          "prioritisation than acting on the national average? ",
          strong("Lift"), " is the prevalence in the model's worst-quintile of ",
          "districts divided by the overall prevalence (1.0 = no better than ",
          "average). ", strong("Hit rate"), " is the share of the truly worst ",
          "quintile the model correctly flags."),
        radioButtons(ns("sq_setting"), "Model setting",
                     choices = c("Transported (realistic)" = "loco",
                                 "Within-country (optimistic)" = "within"),
                     selected = "loco", inline = TRUE),
        reactableOutput(ns("sq_table")),
        methods_note(
          "Computed on the surveyed districts where a survey-observed prevalence ",
          "exists. ", tags$b("Reach@20%"), " is the share of deficient people ",
          "reached by acting in the model's worst 20% of the population (vs 20% ",
          "for no targeting). A lift near 1 or a low hit rate means the model ",
          "adds little targeting value for that outcome/country and a blanket or ",
          "survey-based approach is preferable."
        )
      )
    ),

    nav_panel(
      title = "Fitness for purpose",
      icon = bsicons::bs_icon("clipboard-check"),
      div(
        h5("When should you trust the model for targeting?"),
        p("A traffic-light verdict per outcome × country, based on the ",
          strong("out-of-sample (transported)"), " agreement between predicted ",
          "and survey-observed district prevalence — the regime that matters ",
          "when deploying to areas without a survey."),
        reactableOutput(ns("fp_table")),
        br(),
        h5("Plain-language calibration"),
        p("What a model number actually means on the ground, by prediction band:"),
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("cal_country"), "Country", choices = country_choices,
                      selected = "ghana"),
          selectInput(ns("cal_outcome"), "Outcome", choices = outcome_choices,
                      selected = "women_iron")
        ),
        uiOutput(ns("cal_text")),
        methods_note(
          "Verdict thresholds (transported): ",
          tags$b("Good for targeting"), " = rank correlation r ≥ 0.6 and mean ",
          "absolute error ≤ 7 pp; ", tags$b("Use with caution"), " = r ≥ 0.4; ",
          tags$b("Not reliable for targeting"), " otherwise. These are ",
          "pragmatic cut-offs for ranking districts, not formal tests. The ",
          "within-country correlation (last column) is shown for contrast and is ",
          "always optimistic relative to deployment. Calibration statements come ",
          "from binned out-of-fold predictions; bands with few observations are ",
          "noisier."
        )
      )
    )
  )
}


mod_decision_value_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Panel 1: targeting accuracy ──────────────────────────────────────
    ta_data <- reactive({
      setting <- input$ta_setting %||% "loco"
      oc <- input$ta_outcome %||% "women_iron"
      f <- .dv_frame(setting)
      if (is.null(f)) return(NULL)
      f <- f[f$outcome == oc, , drop = FALSE]
      if (!is.null(input$ta_country) && input$ta_country != "__all__") {
        lbl <- unname(meta$countries[input$ta_country])
        f <- f[f$country_label == lbl, , drop = FALSE]
      }
      if (nrow(f) < 2) return(NULL)
      .dv_add_weight(f, oc)
    })

    ta_curves <- reactive({ d <- ta_data(); if (is.null(d)) NULL else .dv_curves(d) })

    output$vb_capture <- renderText({
      cu <- ta_curves(); if (is.null(cu)) return("—")
      sprintf("%.0f%%", .dv_capture_at(cu$model, 0.20) * 100)
    })
    output$vb_lift <- renderText({
      cu <- ta_curves(); if (is.null(cu)) return("—")
      cap <- .dv_capture_at(cu$model, 0.20)
      if (!is.finite(cap) || cap <= 0) return("—")
      sprintf("%.1f×", cap / 0.20)
    })
    output$vb_hit <- renderText({
      d <- ta_data(); if (is.null(d)) return("—")
      hr <- .dv_hit_rate(d, 0.20); if (is.na(hr)) return("—")
      sprintf("%.0f%%", hr * 100)
    })
    output$vb_n <- renderText({
      cu <- ta_curves(); if (is.null(cu)) "—" else as.character(cu$n)
    })

    output$ta_curve <- renderPlotly({
      cu <- ta_curves()
      if (is.null(cu))
        return(.dv_empty("Not enough districts with truth + prediction for this selection."))
      plot_ly() |>
        add_trace(x = cu$perfect$x, y = cu$perfect$y, type = "scatter",
                  mode = "lines", name = "Perfect targeting",
                  line = list(color = "#1a9850", width = 2, dash = "dash")) |>
        add_trace(x = cu$model$x, y = cu$model$y, type = "scatter",
                  mode = "lines", name = "Model targeting",
                  line = list(color = "#2c7bb6", width = 3)) |>
        add_segments(x = 0, xend = 1, y = 0, yend = 1, name = "No targeting",
                     line = list(color = "#999", dash = "dot"),
                     inherit = FALSE) |>
        layout(xaxis = list(title = "Share of subgroup population reached",
                            tickformat = ".0%", range = c(0, 1)),
               yaxis = list(title = "Share of deficient people reached",
                            tickformat = ".0%", range = c(0, 1)),
               legend = list(x = 0.55, y = 0.15))
    })

    # ── Panel 2: status quo vs model (summary table) ─────────────────────
    output$sq_table <- renderReactable({
      setting <- input$sq_setting %||% "loco"
      f <- .dv_frame(setting)
      validate(need(!is.null(f),
                    "Targeting source not available. Re-run the data-prep step."))
      rows <- list()
      combos <- unique(f[, c("outcome", "country_label")])
      for (i in seq_len(nrow(combos))) {
        oc <- combos$outcome[i]; cl <- combos$country_label[i]
        sub <- f[f$outcome == oc & f$country_label == cl, , drop = FALSE]
        sub <- .dv_add_weight(sub, oc)
        cu <- .dv_curves(sub)
        if (is.null(cu)) next
        ss <- sub[is.finite(sub$true_prev) & is.finite(sub$pred_prev), ]
        overall <- stats::weighted.mean(ss$true_prev, ss$w, na.rm = TRUE)
        k <- max(1L, round(0.20 * nrow(ss)))
        sel <- order(ss$pred_prev, decreasing = TRUE)[seq_len(k)]
        sel_prev <- stats::weighted.mean(ss$true_prev[sel], ss$w[sel], na.rm = TRUE)
        rows[[length(rows) + 1]] <- data.frame(
          Outcome = .dv_pretty_outcome(oc), Country = cl,
          n = cu$n,
          reach20 = .dv_capture_at(cu$model, 0.20),
          lift = if (is.finite(overall) && overall > 0) sel_prev / overall else NA_real_,
          hit20 = .dv_hit_rate(sub, 0.20),
          stringsAsFactors = FALSE)
      }
      validate(need(length(rows) > 0, "No targeting rows for this setting."))
      d <- do.call(rbind, rows)
      reactable(
        d, compact = TRUE, striped = TRUE, searchable = TRUE,
        defaultPageSize = 16, defaultSorted = list(lift = "desc"),
        columns = list(
          Outcome = colDef(minWidth = 170),
          n = colDef(name = "Districts", format = colFormat(digits = 0)),
          reach20 = colDef(name = "Reach @ worst 20%",
                           format = colFormat(percent = TRUE, digits = 0)),
          lift = colDef(name = "Lift vs average",
            format = colFormat(digits = 2, suffix = "×"),
            style = function(value) {
              if (is.na(value)) return(NULL)
              list(color = if (value >= 1.5) "#1a9850"
                           else if (value >= 1.2) "#fdae61" else "#d7191c",
                   fontWeight = "bold")
            }),
          hit20 = colDef(name = "Top-quintile hit rate",
                         format = colFormat(percent = TRUE, digits = 0))
        )
      )
    })

    # ── Panel 3a: fitness-for-purpose scorecard ──────────────────────────
    output$fp_table <- renderReactable({
      d <- transport_cal$adm2_area
      validate(need(!is.null(d) && nrow(d) > 0,
                    "Out-of-sample transport predictions not available."))
      ae <- benchmarks_data$admin2_error
      rows <- list()
      combos <- unique(d[, c("outcome", "country")])
      for (i in seq_len(nrow(combos))) {
        oc <- combos$outcome[i]; ck <- combos$country[i]
        cl <- unname(meta$countries[ck])
        sub <- d[d$outcome == oc & d$country == ck, , drop = FALSE]
        sub <- sub[is.finite(sub$true_prev) & is.finite(sub$pred_prev), ]
        if (nrow(sub) < 3) next
        r <- suppressWarnings(stats::cor(sub$pred_prev, sub$true_prev))
        mae <- mean(abs(sub$pred_prev - sub$true_prev), na.rm = TRUE) * 100
        verdict <- if (!is.na(r) && r >= 0.6 && mae <= 7) "Good for targeting"
                   else if (!is.na(r) && r >= 0.4) "Use with caution"
                   else "Not reliable for targeting"
        within_r <- NA_real_
        if (!is.null(ae)) {
          m <- ae[ae$outcome == oc & ae$country == cl, , drop = FALSE]
          if (nrow(m) > 0) within_r <- m$pearson_r[1]
        }
        rows[[length(rows) + 1]] <- data.frame(
          Outcome = .dv_pretty_outcome(oc), Country = cl, n = nrow(sub),
          transport_r = round(r, 2), mae_pp = round(mae, 1),
          verdict = verdict, within_r = within_r, stringsAsFactors = FALSE)
      }
      validate(need(length(rows) > 0, "No fitness-for-purpose rows."))
      tab <- do.call(rbind, rows)
      reactable(
        tab, compact = TRUE, striped = TRUE, searchable = TRUE,
        defaultPageSize = 16,
        columns = list(
          Outcome = colDef(minWidth = 170),
          n = colDef(name = "Districts", format = colFormat(digits = 0)),
          transport_r = colDef(name = "Transported r"),
          mae_pp = colDef(name = "MAE (pp)"),
          within_r = colDef(name = "Within-country r",
                            format = colFormat(digits = 2)),
          verdict = colDef(name = "Verdict", minWidth = 170,
            style = function(value) {
              col <- switch(value,
                "Good for targeting" = "#1a9850",
                "Use with caution" = "#fd8d3c",
                "#d7191c")
              list(color = col, fontWeight = "bold")
            })
        )
      )
    })

    # ── Panel 3b: plain-language calibration ─────────────────────────────
    output$cal_text <- renderUI({
      cal <- model_diag$calibration
      if (is.null(cal) || nrow(cal) == 0)
        return(p(em("Calibration data not available."), style = "color:#888;"))
      cl <- unname(meta$countries[input$cal_country])
      oc <- input$cal_outcome
      d <- cal[cal$country == cl & cal$outcome == oc, , drop = FALSE]
      if (nrow(d) == 0)
        return(p(em(sprintf("No calibration bins for %s × %s.",
                            cl, .dv_pretty_outcome(oc))), style = "color:#888;"))
      d <- d[is.finite(d$mean_pred) & is.finite(d$mean_obs), , drop = FALSE]
      if (nrow(d) == 0)
        return(p(em("Calibration bins for this selection have no usable values."),
                 style = "color:#888;"))
      d <- d[order(d$mean_pred), ]
      items <- lapply(seq_len(nrow(d)), function(i) {
        mp <- d$mean_pred[i] * 100; mo <- d$mean_obs[i] * 100
        n <- if (is.finite(d$n[i])) d$n[i] else 0L
        dirn <- if (mo > mp + 2) "higher than predicted"
                else if (mo < mp - 2) "lower than predicted"
                else "about as predicted"
        tags$li(sprintf(
          "When the model predicts around %.0f%%, the survey-observed rate averages %.0f%% (%s; n = %s).",
          mp, mo, dirn, formatC(n, big.mark = ",", format = "d")))
      })
      tags$ul(items)
    })

  })
}
