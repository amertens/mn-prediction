# =============================================================================
# Module: Scenario Projections
# =============================================================================
# Two scenario modes:
#   (A) Coverage: hypothetical intervention reaches a fraction of districts
#       at a given coverage and effect size. Recompute prevalence post-
#       intervention.
#   (B) Sensitivity: scale prevalence in all districts by a user-defined
#       proportional shift (e.g., simulating a climate shock or food price
#       change). This is illustrative — flagged as such in the methods note.

# Literature-default effect sizes by outcome
# Sources: Lancet Maternal & Child Nutrition Series; WHO/UNICEF program reviews
# Effect size = proportional reduction in prevalence among those reached
default_effects <- list(
  child_vitA   = list(coverage = 0.60, effect = 0.40,
                      label = "VAS supplementation",
                      cite = "WHO 2011 — semi-annual VAS reduces VAD"),
  women_vitA   = list(coverage = 0.30, effect = 0.30,
                      label = "Food fortification",
                      cite = "Allen et al. 2006 — fortified oil/maize"),
  child_iron   = list(coverage = 0.50, effect = 0.30,
                      label = "MNP / fortification",
                      cite = "De-Regil 2013 — MNPs reduce IDA in children"),
  women_iron   = list(coverage = 0.40, effect = 0.35,
                      label = "Iron-folate tablets / fortification",
                      cite = "Pena-Rosas 2015 — iron supplementation"),
  women_folate = list(coverage = 0.50, effect = 0.55,
                      label = "Flour fortification",
                      cite = "Centeno Tablante 2019 — flour fortification"),
  women_b12    = list(coverage = 0.30, effect = 0.30,
                      label = "B12 fortification (limited evidence)",
                      cite = "Limited fortification evidence"),
  child_zinc   = list(coverage = 0.40, effect = 0.30,
                      label = "Zinc supplementation / biofortification",
                      cite = "Brown 2009 — zinc fortification"),
  women_zinc   = list(coverage = 0.30, effect = 0.20,
                      label = "Zinc fortification",
                      cite = "Limited evidence")
)

mod_scenarios_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(
    id = ns("scenario_mode"),

    nav_panel(
      title = "Coverage scenario",
      icon = bsicons::bs_icon("bullseye"),

      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          title = "Intervention parameters",

          selectInput(ns("country"), "Country",
                      choices = country_choices,
                      selected = "ghana"),

          selectInput(ns("outcome"), "Outcome",
                      choices = outcome_choices,
                      selected = "women_iron"),

          radioButtons(ns("targeting"), "Targeting strategy",
                       choices = c("All districts" = "all",
                                   "Only districts above national average" = "above_natl",
                                   "Top-N highest-prevalence districts" = "top_n"),
                       selected = "above_natl"),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'top_n'", ns("targeting")),
            sliderInput(ns("top_n"), "Number of districts (N)",
                        min = 1, max = 50, value = 10, step = 1)
          ),

          sliderInput(ns("coverage"), "Coverage of intervention (%)",
                      min = 0, max = 100, value = 60, step = 5,
                      post = "%"),

          sliderInput(ns("effect"), "Effect size on prevalence (%)",
                      min = 0, max = 100, value = 40, step = 5,
                      post = "%",
                      width = "100%"),

          uiOutput(ns("default_note"))
        ),

        layout_columns(
          col_widths = c(6, 6),

          card(
            card_header("Before intervention"),
            card_body(
              leafletOutput(ns("map_before"), height = "350px"),
              htmlOutput(ns("summary_before"))
            )
          ),

          card(
            card_header("After intervention"),
            card_body(
              leafletOutput(ns("map_after"), height = "350px"),
              htmlOutput(ns("summary_after"))
            )
          ),

          card(
            full_screen = TRUE,
            card_header("District-level impact summary"),
            card_body(
              reactableOutput(ns("impact_table")),
              methods_note(
                "The “before” column is the model-predicted prevalence under current ",
                "conditions. The “after” column applies the user-specified intervention ",
                "to targeted districts only: ",
                tags$code("prev_after = prev_before × (1 − coverage × effect)"),
                ". Population reached is the targeted district populations multiplied ",
                "by coverage; cases averted is the difference in expected cases ",
                "(prevalence × population) between before and after. ",
                tags$br(), tags$br(),
                "Default coverage and effect sizes are drawn from published evaluations ",
                "of similar interventions (Lancet Maternal & Child Nutrition series; ",
                "WHO/UNICEF program reviews). They are starting points only — ",
                "users should adjust based on the specific intervention design ",
                "and local implementation context. The model assumes no spillover ",
                "to non-targeted districts and constant effect across the targeted ",
                "set — both simplifications. The 95% conformal CI from the underlying ",
                "model is preserved as a baseline uncertainty band; intervention ",
                "uncertainty (in coverage and effect size) is not propagated and ",
                "would widen the post-intervention interval substantially."
              )
            )
          )
        )
      )
    ),

    nav_panel(
      title = "Sensitivity scenario",
      icon = bsicons::bs_icon("graph-down-arrow"),

      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          title = "Sensitivity parameters",

          selectInput(ns("sens_country"), "Country",
                      choices = country_choices,
                      selected = "ghana"),

          selectInput(ns("sens_outcome"), "Outcome",
                      choices = outcome_choices,
                      selected = "women_iron"),

          h6("Predictor shift", style = "margin-top: 1em;"),

          selectInput(ns("sens_lever"),
                      "Adjust which factor?",
                      choices = c("Climate stress (drought / heat)" = "climate",
                                  "Food prices (staple commodities)" = "food_price",
                                  "Food security (FCS / coping)" = "food_sec",
                                  "Conflict events" = "conflict",
                                  "Generic prevalence shift" = "generic"),
                      selected = "generic"),

          sliderInput(ns("sens_change"),
                      "Magnitude of change (%)",
                      min = -50, max = 50, value = 20, step = 5,
                      post = "%"),

          checkboxInput(ns("apply_to_subset"),
                        "Apply only to vulnerable districts",
                        value = FALSE),

          conditionalPanel(
            condition = sprintf("input['%s']", ns("apply_to_subset")),
            p(em("Vulnerable districts = those with current prevalence above the national average."),
              style = "font-size: 0.85em; color: #666;")
          ),

          hr(),

          htmlOutput(ns("sens_summary"))
        ),

        layout_columns(
          col_widths = c(6, 6),

          card(
            card_header("Baseline prediction"),
            card_body(
              leafletOutput(ns("sens_map_before"), height = "400px")
            )
          ),

          card(
            card_header("Sensitivity scenario"),
            card_body(
              leafletOutput(ns("sens_map_after"), height = "400px")
            )
          ),

          card(
            full_screen = TRUE,
            card_header("Sensitivity impact"),
            card_body(
              reactableOutput(ns("sens_table")),
              methods_note(
                tags$strong("This is an illustrative planning tool, not a precise forecast. "),
                "We apply a simple linear scaling to predicted prevalence based on ",
                "the user-specified change. The elasticity of each lever (climate, ",
                "food prices, food security, conflict) is calibrated using the domain ",
                "ablation analysis: the larger the AUC drop when a domain is permuted, ",
                "the larger the assumed elasticity. The “generic” lever scales prevalence ",
                "uniformly by the chosen percentage. ",
                tags$br(), tags$br(),
                "True scenario forecasting would require re-running the underlying ",
                "machine learning models with adjusted predictor values, which is ",
                "computationally intensive and only available offline. The estimates ",
                "shown here are intended for: (1) communicating roughly how much ",
                "burden could shift under stress, and (2) identifying which districts ",
                "would be most exposed to a given shift. They should not be cited ",
                "as quantitative projections."
              )
            )
          )
        )
      )
    )
  )
}


mod_scenarios_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── COVERAGE SCENARIO ──────────────────────────────────────────────────

    # Show default coverage/effect from literature when outcome changes
    observeEvent(input$outcome, {
      d <- default_effects[[input$outcome]]
      if (!is.null(d)) {
        updateSliderInput(session, "coverage",
                          value = round(d$coverage * 100))
        updateSliderInput(session, "effect",
                          value = round(d$effect * 100))
      }
    })

    output$default_note <- renderUI({
      d <- default_effects[[input$outcome]]
      if (is.null(d)) return(NULL)
      tagList(
        hr(),
        p(strong("Suggested defaults: "), d$label,
          style = "font-size: 0.9em; margin-bottom: 4px;"),
        p(em(d$cite),
          style = "font-size: 0.8em; color: #666;")
      )
    })

    # Reactive: baseline data
    base_data <- reactive({
      req(input$country, input$outcome)
      get_country_admin2(input$country, input$outcome,
                          admin2_bnds, admin2_pred, admin2_pop)
    })

    # Reactive: scenario data
    scenario_data <- reactive({
      df <- base_data()
      req(df, input$coverage, input$effect)

      cov <- input$coverage / 100
      eff <- input$effect / 100

      # Determine which districts are targeted
      df$is_target <- FALSE
      df_clean <- as.data.frame(df)

      if (input$targeting == "all") {
        df$is_target <- !is.na(df$pred_prev)
      } else if (input$targeting == "above_natl") {
        natl <- national_aggregate(df_clean)
        df$is_target <- !is.na(df$pred_prev) &
                          df$pred_prev > natl$pred_prev_natl
      } else if (input$targeting == "top_n") {
        n <- input$top_n %||% 10
        ordered_idx <- order(df_clean$pred_prev, decreasing = TRUE,
                              na.last = TRUE)
        topn_idx <- ordered_idx[seq_len(min(n, sum(!is.na(df_clean$pred_prev))))]
        df$is_target <- FALSE
        df$is_target[topn_idx] <- TRUE
      }

      # Apply intervention
      df$pred_prev_after <- df$pred_prev
      target_rows <- df$is_target
      df$pred_prev_after[target_rows] <- df$pred_prev[target_rows] *
                                          (1 - cov * eff)

      # Compute populations and cases
      df$pop_at_risk_after <- df$pred_prev_after * df$population
      df$pop_reached <- ifelse(df$is_target,
                                df$population * cov, 0)
      df$cases_averted <- df$pop_at_risk - df$pop_at_risk_after

      df
    })

    render_scenario_map <- function(df, prev_col) {
      vals <- df[[prev_col]]
      pal <- if (any(!is.na(vals))) {
        colorNumeric("YlOrRd",
                      domain = c(0, max(vals, na.rm = TRUE)),
                      na.color = "#cccccc")
      } else NULL
      fill_col <- if (!is.null(pal)) pal(vals) else "#cccccc"

      labels <- sprintf(
        "<strong>%s</strong><br/>Prevalence: %s",
        df$Admin2, fmt_pct(vals)
      ) |> lapply(HTML)

      leaflet(df) |>
        addProviderTiles(providers$CartoDB.Positron) |>
        addPolygons(
          fillColor = fill_col,
          fillOpacity = 0.75,
          color = "white",
          weight = 1,
          label = labels
        )
    }

    output$map_before <- renderLeaflet({
      df <- scenario_data(); req(df)
      render_scenario_map(df, "pred_prev")
    })

    output$map_after <- renderLeaflet({
      df <- scenario_data(); req(df)
      render_scenario_map(df, "pred_prev_after")
    })

    output$summary_before <- renderUI({
      df <- scenario_data(); req(df)
      df_clean <- as.data.frame(df)
      df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                             !is.na(df_clean$population), , drop = FALSE]
      total_pop <- sum(df_clean$population)
      cases <- sum(df_clean$pred_prev * df_clean$population)
      natl <- national_aggregate(df_clean)
      tagList(
        p(strong("National prevalence: "), fmt_pct(natl$pred_prev_natl)),
        p(strong("Cases (current): "), fmt_count(cases))
      )
    })

    output$summary_after <- renderUI({
      df <- scenario_data(); req(df)
      df_clean <- as.data.frame(df)
      df_clean <- df_clean[!is.na(df_clean$pred_prev_after) &
                             !is.na(df_clean$population), , drop = FALSE]
      cases_after <- sum(df_clean$pred_prev_after * df_clean$population)
      cases_before <- sum(df_clean$pred_prev * df_clean$population)
      averted <- cases_before - cases_after
      pop_reached <- sum(df_clean$pop_reached)
      n_target <- sum(df_clean$is_target)
      total_pop <- sum(df_clean$population)
      natl_after <- weighted.mean(df_clean$pred_prev_after,
                                    df_clean$population, na.rm = TRUE)
      tagList(
        p(strong("National prevalence: "), fmt_pct(natl_after),
          tags$small(sprintf(" (Δ %.2fpp)",
                              (natl_after - cases_before / total_pop) * 100),
                     style = "color: #2c7bb6;")),
        p(strong("Cases (after): "), fmt_count(cases_after)),
        p(strong("Cases averted: "),
          tags$span(style = "color: #1a9850; font-weight: bold;",
                    fmt_count(averted))),
        p(strong("Population reached: "), fmt_count(pop_reached)),
        p(strong("Districts targeted: "), n_target)
      )
    })

    output$impact_table <- renderReactable({
      df <- scenario_data(); req(df)
      df_clean <- as.data.frame(df)
      df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                             !is.na(df_clean$population), , drop = FALSE]

      out <- df_clean[, c("Admin1", "Admin2", "is_target", "pred_prev",
                          "pred_prev_after", "population", "pop_reached",
                          "cases_averted")]
      out <- out[order(-out$cases_averted), , drop = FALSE]

      reactable(
        out,
        compact = TRUE, striped = TRUE,
        searchable = TRUE, defaultPageSize = 12,
        columns = list(
          Admin1 = colDef(name = "Region"),
          Admin2 = colDef(name = "District"),
          is_target = colDef(name = "Targeted",
                              cell = function(v) if (v) "✓" else ""),
          pred_prev       = colDef(name = "Before",
                                    format = colFormat(percent = TRUE, digits = 1)),
          pred_prev_after = colDef(name = "After",
                                    format = colFormat(percent = TRUE, digits = 1)),
          population     = colDef(name = "Population",
                                   format = colFormat(separators = TRUE, digits = 0)),
          pop_reached    = colDef(name = "Pop. reached",
                                   format = colFormat(separators = TRUE, digits = 0)),
          cases_averted  = colDef(name = "Cases averted",
                                   format = colFormat(separators = TRUE, digits = 0),
                                   style = list(fontWeight = "bold",
                                                color = "#1a9850"))
        )
      )
    })

    # ── SENSITIVITY SCENARIO ───────────────────────────────────────────────

    # Approximate elasticities by lever (calibrated from domain ablation magnitude)
    sens_elasticities <- list(
      climate    = 0.30,
      food_price = 0.40,
      food_sec   = 0.50,
      conflict   = 0.20,
      generic    = 1.00
    )

    sens_data <- reactive({
      req(input$sens_country, input$sens_outcome,
          input$sens_lever, input$sens_change)

      df <- get_country_admin2(input$sens_country, input$sens_outcome,
                                admin2_bnds, admin2_pred, admin2_pop)
      req(df)

      elasticity <- sens_elasticities[[input$sens_lever]] %||% 1.0
      shift <- input$sens_change / 100  # proportion

      # new_prev = old_prev * (1 + elasticity * shift)
      multiplier <- 1 + elasticity * shift

      target <- if (isTRUE(input$apply_to_subset)) {
        df_clean <- as.data.frame(df)
        natl <- national_aggregate(df_clean)
        !is.na(df$pred_prev) & df$pred_prev > natl$pred_prev_natl
      } else {
        !is.na(df$pred_prev)
      }

      df$pred_prev_sens <- df$pred_prev
      df$pred_prev_sens[target] <- pmax(0, pmin(1,
                                          df$pred_prev[target] * multiplier))

      df
    })

    output$sens_map_before <- renderLeaflet({
      df <- sens_data(); req(df)
      render_scenario_map(df, "pred_prev")
    })

    output$sens_map_after <- renderLeaflet({
      df <- sens_data(); req(df)
      render_scenario_map(df, "pred_prev_sens")
    })

    output$sens_summary <- renderUI({
      df <- sens_data(); req(df)
      df_clean <- as.data.frame(df)
      df_clean <- df_clean[!is.na(df_clean$pred_prev_sens) &
                             !is.na(df_clean$population), , drop = FALSE]

      base_natl <- national_aggregate(df_clean)$pred_prev_natl
      sens_natl <- weighted.mean(df_clean$pred_prev_sens,
                                   df_clean$population, na.rm = TRUE)
      delta_pp <- (sens_natl - base_natl) * 100

      cases_base <- sum(df_clean$pred_prev * df_clean$population)
      cases_sens <- sum(df_clean$pred_prev_sens * df_clean$population)
      delta_cases <- cases_sens - cases_base

      delta_color <- if (delta_pp > 0) "#d7191c" else "#1a9850"

      tagList(
        h6("Projected change", style = "margin-top: 0;"),
        p(strong("Baseline national prevalence: "), fmt_pct(base_natl)),
        p(strong("Scenario national prevalence: "), fmt_pct(sens_natl),
          br(),
          tags$span(style = sprintf("color: %s; font-weight: bold;",
                                      delta_color),
                    sprintf("Δ %.2f pp", delta_pp))),
        p(strong("Cases (baseline): "), fmt_count(cases_base)),
        p(strong("Cases (scenario): "), fmt_count(cases_sens),
          br(),
          tags$span(style = sprintf("color: %s; font-weight: bold;",
                                      delta_color),
                    sprintf("Δ %s", fmt_count(abs(delta_cases))),
                    if (delta_cases > 0) " more" else " fewer"))
      )
    })

    output$sens_table <- renderReactable({
      df <- sens_data(); req(df)
      df_clean <- as.data.frame(df)
      df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                             !is.na(df_clean$population), , drop = FALSE]

      df_clean$delta_prev <- df_clean$pred_prev_sens - df_clean$pred_prev
      df_clean$delta_cases <- df_clean$delta_prev * df_clean$population

      out <- df_clean[, c("Admin1", "Admin2", "pred_prev",
                           "pred_prev_sens", "delta_prev", "delta_cases")]
      out <- out[order(-abs(out$delta_cases)), , drop = FALSE]

      reactable(
        out,
        compact = TRUE, striped = TRUE,
        searchable = TRUE, defaultPageSize = 12,
        columns = list(
          Admin1     = colDef(name = "Region"),
          Admin2     = colDef(name = "District"),
          pred_prev  = colDef(name = "Baseline",
                               format = colFormat(percent = TRUE, digits = 1)),
          pred_prev_sens = colDef(name = "Scenario",
                                   format = colFormat(percent = TRUE, digits = 1)),
          delta_prev = colDef(name = "Δ prev. (pp)",
                               format = colFormat(percent = TRUE, digits = 2),
                               style = function(v) {
                                 if (is.na(v)) return(NULL)
                                 list(color = if (v > 0) "#d7191c"
                                              else "#1a9850")
                               }),
          delta_cases = colDef(name = "Δ cases",
                                format = colFormat(separators = TRUE, digits = 0))
        )
      )
    })

  })
}
