# =============================================================================
# Module: Predictor Importance
# =============================================================================
# Surfaces three views of which factors drive predicted prevalence:
#   (1) Domain ablation heatmap — AUC drop when each conceptual domain is
#       permuted, by country × outcome.
#   (2) Top single-variable importance ranking — the predictors that lose
#       the most AUC when individually shuffled.
#   (3) Per-district top SHAP factors — for a specific district, which
#       predictors push the prediction up or down.

mod_importance_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "Domain importance heatmap",
      icon = bsicons::bs_icon("grid-3x3"),
      div(
        p("Domain importance is measured by permutation: each thematic ",
          "group of predictors (Climate, Food Security, Demographics, etc.) ",
          "is shuffled in turn, and we measure the resulting drop in model ",
          "performance (ROC-AUC). Larger drops indicate the domain ",
          "contributes more signal."),
        plotlyOutput(ns("ablation_heatmap"), height = "550px"),
        methods_note(
          "The heatmap shows, for each country × outcome model, the AUC drop ",
          "when each thematic predictor domain is removed (permuted). Darker ",
          "cells indicate domains that the model relies on more heavily for ",
          "that specific outcome and country. ",
          tags$br(), tags$br(),
          "This is a population-level importance measure, not an explanation ",
          "of any single district's prediction. For per-district explanations ",
          "see the “District factors” sub-tab. Domains with near-zero drops ",
          "are either redundant (their signal is captured by other domains) ",
          "or genuinely uninformative for the outcome. ",
          tags$br(), tags$br(),
          "Note that permutation importance can underestimate the contribution ",
          "of correlated predictors. If two domains carry overlapping ",
          "information, removing either one alone may not hurt performance ",
          "much, even though the combined signal is essential."
        )
      )
    ),

    nav_panel(
      title = "Top variables",
      icon = bsicons::bs_icon("list-ol"),
      div(
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("var_country"), "Country",
                      choices = country_choices, selected = "ghana"),
          selectInput(ns("var_outcome"), "Outcome",
                      choices = outcome_choices, selected = "women_iron")
        ),
        plotlyOutput(ns("varimp_plot"), height = "500px"),
        methods_note(
          "Single-variable permutation importance: each individual predictor ",
          "is shuffled three times, the model re-predicts, and we compute the ",
          "average AUC drop. Larger bars indicate predictors the model relies ",
          "on more heavily. The top 30 variables are shown. ",
          tags$br(), tags$br(),
          "Variable names follow these prefixes: ", tags$code("chirps_*"),
          " (rainfall), ", tags$code("gee_*"), " (satellite-derived), ",
          tags$code("MAP_*"), " (Malaria Atlas Project), ",
          tags$code("worldpop_*"), " (population), ",
          tags$code("wfp_*"), " (food security or food prices), ",
          tags$code("dhs*_*"), " (DHS Admin-2 indicators), ",
          tags$code("ihme_*"), " (IHME modeled health indicators), ",
          tags$code("soil_*"), " (soil properties)."
        )
      )
    ),

    nav_panel(
      title = "District factors",
      icon = bsicons::bs_icon("crosshair"),
      div(
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput(ns("shap_country"), "Country",
                      choices = country_choices, selected = "ghana"),
          selectInput(ns("shap_outcome"), "Outcome",
                      choices = outcome_choices, selected = "women_iron"),
          uiOutput(ns("shap_district_picker"))
        ),
        plotlyOutput(ns("shap_plot"), height = "450px"),
        methods_note(
          "SHAP (SHapley Additive exPlanations) values quantify how much ",
          "each predictor contributes to a single district's prediction. ",
          "Positive values (red) indicate predictors that push the prediction ",
          "above the model's average; negative values (blue) push it below. ",
          tags$br(), tags$br(),
          "These attributions are extracted from the xgboost base learner ",
          "within the SuperLearner ensemble, which means they reflect ",
          "xgboost's view of feature contributions rather than the full ",
          "ensemble. They are accurate in direction and approximately right ",
          "in magnitude when xgboost has substantial weight in the ensemble. ",
          tags$br(), tags$br(),
          "These per-district explanations are particularly useful for ",
          "ground-truthing model predictions: if the “top factors” for a ",
          "district align with what local health officials know about the ",
          "area, that supports the model. If they do not, that prompts ",
          "further investigation."
        )
      )
    )
  )
}


mod_importance_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Domain ablation heatmap ──────────────────────────────────────────
    output$ablation_heatmap <- renderPlotly({
      abl <- importance_data$ablation
      if (is.null(abl) || nrow(abl) == 0) {
        return(plotly_empty(type = "scatter", mode = "markers") |>
          layout(annotations = list(text = "Domain ablation data not yet built. Re-run pipeline.",
                                      showarrow = FALSE,
                                      font = list(size = 14))))
      }

      # Compute baseline AUC per country×outcome (domain_removed == "none")
      base <- abl[abl$domain_removed == "none", c("country", "outcome", "auc")]
      colnames(base)[3] <- "auc_baseline"
      d <- merge(abl[abl$domain_removed != "none", ], base,
                  by = c("country", "outcome"), all.x = TRUE)
      d$drop <- d$auc_baseline - d$auc

      d$label <- paste0(d$country, " — ", outcome_labels[d$outcome])

      plot_ly(d,
              x = ~domain_removed, y = ~label, z = ~drop,
              type = "heatmap",
              colorscale = "YlOrRd",
              hoverinfo = "text",
              text = ~sprintf("%s × %s<br>Domain: %s<br>AUC drop: %.3f",
                                country, outcome, domain_removed, drop)) |>
        layout(xaxis = list(title = "Domain", tickangle = -45),
               yaxis = list(title = ""),
               margin = list(l = 250, b = 100))
    })

    # ── Single-variable importance plot ──────────────────────────────────
    output$varimp_plot <- renderPlotly({
      vi <- importance_data$varimp
      ctry_label <- meta$countries[input$var_country]
      if (is.null(vi) || nrow(vi) == 0) {
        return(plotly_empty(type = "scatter", mode = "markers") |>
          layout(annotations = list(
            text = "Per-variable importance data not yet built. Re-run pipeline with single_var_ablation targets.",
            showarrow = FALSE, font = list(size = 14))))
      }

      d <- vi[vi$country == ctry_label & vi$outcome == input$var_outcome, ,
              drop = FALSE]
      if (nrow(d) == 0) {
        return(plotly_empty(type = "scatter", mode = "markers") |>
          layout(annotations = list(
            text = sprintf("No variable importance data for %s × %s",
                            ctry_label, input$var_outcome),
            showarrow = FALSE, font = list(size = 14))))
      }

      d$variable <- factor(d$variable, levels = rev(d$variable[order(d$drop_mean)]))

      plot_ly(d,
              x = ~drop_mean, y = ~variable,
              type = "bar", orientation = "h",
              marker = list(color = "#2c7bb6"),
              hoverinfo = "text",
              text = ~sprintf("%s<br>AUC drop: %.4f", variable, drop_mean)) |>
        layout(xaxis = list(title = "AUC drop when permuted"),
               yaxis = list(title = ""),
               margin = list(l = 250))
    })

    # ── Per-district SHAP picker ─────────────────────────────────────────
    output$shap_district_picker <- renderUI({
      sh <- importance_data$shap
      ctry_label <- meta$countries[input$shap_country]
      if (is.null(sh) || nrow(sh) == 0) {
        return(p(em("SHAP data not yet available."), style = "color: #888;"))
      }
      ds <- sort(unique(sh$Admin2[sh$country == ctry_label &
                                     sh$outcome == input$shap_outcome]))
      if (length(ds) == 0) {
        return(p(em("No districts with SHAP data."), style = "color: #888;"))
      }
      selectInput(ns("shap_district"), "District",
                  choices = ds, selected = ds[1])
    })

    output$shap_plot <- renderPlotly({
      sh <- importance_data$shap
      ctry_label <- meta$countries[input$shap_country]
      req(input$shap_district, sh)

      d <- sh[sh$country == ctry_label &
              sh$outcome == input$shap_outcome &
              sh$Admin2 == input$shap_district, , drop = FALSE]
      if (nrow(d) == 0) {
        return(plotly_empty(type = "scatter", mode = "markers"))
      }

      d <- d[order(d$rank), ]
      d$variable <- factor(d$variable, levels = rev(d$variable))
      d$bar_color <- ifelse(d$shap > 0, "#d7191c", "#2c7bb6")

      plot_ly(d,
              x = ~shap, y = ~variable,
              type = "bar", orientation = "h",
              marker = list(color = ~bar_color),
              hoverinfo = "text",
              text = ~sprintf("%s<br>SHAP: %+.4f<br>%s",
                                variable, shap, direction)) |>
        layout(xaxis = list(title = "SHAP value (effect on predicted prevalence)",
                             zeroline = TRUE, zerolinecolor = "#333",
                             zerolinewidth = 1),
               yaxis = list(title = ""),
               margin = list(l = 250))
    })

  })
}
