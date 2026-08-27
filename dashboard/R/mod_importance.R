# =============================================================================
# Module: Predictor Importance
# =============================================================================
# Three views of what drives a predicted prevalence, ordered so the one a
# program decision actually turns on comes first:
#   (1) Per-district SHAP factors — why THIS district's estimate is where it is,
#       and the panel to sense-check against local knowledge.
#   (2) Domain ablation heatmap — AUC drop when each predictor family is
#       permuted, by country × outcome.
#   (3) Top single-variable importance — the predictors that lose the most AUC
#       when individually shuffled.
#
# (1) used to be last. It is the only one framed around a place rather than
# around the model, so it now leads.

mod_importance_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    # "Why is this district high?" comes first. It used to be third, behind two
    # panels about the model rather than about a place.
    nav_panel(
      title = "Why a district is high",
      icon = bsicons::bs_icon("crosshair"),
      div(
        p("Pick a district to see which conditions pushed its estimate up or ",
          "down: malaria burden, rainfall, soil, food prices, and so on. ",
          "Read it as the model showing its working, and check it against what ",
          "you know about the place."),
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
          "Bars are SHAP values: each predictor's contribution to this ",
          "district's estimate. Red pushes the estimate above the country ",
          "average, blue pulls it below, and bar length is how much. ",
          tags$br(), tags$br(),
          "They are computed against the fitted ensemble's own predictions ",
          "using sampling (Monte-Carlo) Shapley values, so they reflect the ",
          "whole model rather than one component. To keep computation ",
          "manageable, they cover each outcome's most important predictors over ",
          "a sample of individuals per district. ",
          tags$br(), tags$br(),
          "This is the most useful panel for sense-checking. If a district's ",
          "top factors match what local health staff know about the area, that ",
          "is a point in the model's favor. If they do not, the estimate is ",
          "worth questioning."
        )
      )
    ),

    nav_panel(
      title = "Which data sources matter",
      icon = bsicons::bs_icon("grid-3x3"),
      div(
        p("How much each family of predictors contributes, by country and ",
          "nutrient. Each family is scrambled in turn and we measure how much ",
          "worse the model gets. Darker cells mean the model leans on that ",
          "family more."),
        plotlyOutput(ns("ablation_heatmap"), height = "550px"),
        methods_note(
          "Cells show the drop in ROC-AUC when a thematic group of predictors ",
          "is permuted, for each country × outcome model. ",
          tags$br(), tags$br(),
          "This describes the model overall, not any one district — for that, ",
          "use the first panel. A near-zero cell means the family is either ",
          "uninformative for that outcome or redundant with another family. ",
          tags$br(), tags$br(),
          "The measure understates predictors that overlap: when two families ",
          "carry the same information, removing either alone changes little, ",
          "even though the information itself matters."
        )
      )
    ),

    nav_panel(
      title = "Individual predictors",
      icon = bsicons::bs_icon("list-ol"),
      div(
        p("The thirty predictors the model relies on most, for one country and ",
          "nutrient."),
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("var_country"), "Country",
                      choices = country_choices, selected = "ghana"),
          selectInput(ns("var_outcome"), "Outcome",
                      choices = outcome_choices, selected = "women_iron")
        ),
        plotlyOutput(ns("varimp_plot"), height = "500px"),
        methods_note(
          "Each predictor is scrambled three times, the model re-predicts, and ",
          "the bar is the average drop in ROC-AUC. Longer bars mean heavier ",
          "reliance. ",
          tags$br(), tags$br(),
          "Name prefixes: ", tags$code("chirps_*"), " rainfall, ",
          tags$code("gee_*"), " satellite-derived, ",
          tags$code("MAP_*"), " Malaria Atlas Project, ",
          tags$code("worldpop_*"), " population, ",
          tags$code("wfp_*"), " food security and food prices, ",
          tags$code("dhs*_*"), " DHS district indicators, ",
          tags$code("ihme_*"), " IHME modeled health indicators, ",
          tags$code("soil_*"), " soil properties."
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

      d$label <- paste0(d$country, " — ", meta$outcome_labels[d$outcome])

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
