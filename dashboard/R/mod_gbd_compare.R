# =============================================================================
# Module: GBD / IHME Comparison
# =============================================================================
# Compares pipeline national-level estimates against the Global Burden of
# Disease (GBD) study estimates from IHME. The framework supports either
# placeholder or actual GBD Results Tool data; an RA task is surfaced when
# placeholders are in use.

mod_gbd_compare_ui <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = 12,

    # RA notice (only shown if data is placeholder)
    uiOutput(ns("placeholder_warning")),

    card(
      card_header("Pipeline vs. GBD national estimates"),
      card_body(
        p("This view compares national prevalence estimates from this project's ",
          "subnational modelling pipeline against the Global Burden of Disease ",
          "(GBD) study estimates published by the Institute for Health Metrics ",
          "and Evaluation (IHME). Discrepancies highlight where subnational ",
          "modelling reveals patterns that may differ from GBD's country-level ",
          "approach."),

        selectInput(ns("year_filter"), "Reference year",
                    choices = c("Survey year (matched per country)" = "survey",
                                "2023 projection" = "2023"),
                    selected = "survey"),

        reactableOutput(ns("compare_table")),
        methods_note(
          "GBD prevalence estimates are sourced from the IHME Global Burden of Disease ",
          "study and reflect a country-year point estimate with uncertainty bounds. ",
          "Pipeline estimates are population-weighted national averages of district-level ",
          "model predictions, with conformal 95% prediction intervals propagated by ",
          "the same population weights. ",
          tags$br(), tags$br(),
          "When the two estimates agree, that corroborates the national-level ",
          "estimate. When they disagree, the gap highlights where one method may be ",
          "missing important signal — e.g., pipeline estimates may capture subnational ",
          "heterogeneity that population-level GBD modelling smooths over, while GBD ",
          "estimates may incorporate broader temporal trends and external indicators ",
          "that the pipeline does not. ",
          tags$br(), tags$br(),
          "The “Diff (pp)” column is pipeline minus GBD in percentage points; positive ",
          "values mean the pipeline reports higher prevalence than GBD."
        )
      )
    ),

    card(
      card_header("Hidden burden vs. GBD"),
      card_body(
        p("A useful framing for funders: how many people live in districts where ",
          "the pipeline's local estimate substantially diverges from GBD's national ",
          "average? These are individuals whose burden may be mischaracterized by ",
          "country-level estimates alone."),
        sliderInput(ns("hidden_threshold"),
                    "Threshold for ‘substantial divergence’ (percentage points)",
                    min = 1, max = 20, value = 5, step = 1),
        reactableOutput(ns("hidden_gbd_table")),
        methods_note(
          "For each country × outcome, we identify districts where the pipeline's ",
          "local prediction differs from the GBD country-level estimate by more than ",
          "the user-set threshold. The “mischaracterized population” is the sum of ",
          "the populations of those districts. ",
          tags$br(), tags$br(),
          "This metric is intended to quantify, in absolute terms, the value of subnational ",
          "data for resource allocation. A high mischaracterized population suggests that ",
          "country-level estimates would lead to substantial misallocation — either ",
          "underservicing high-burden districts or overservicing low-burden ones. ",
          tags$br(), tags$br(),
          "Note that the pipeline and GBD use different methodologies; some divergence ",
          "is expected. The threshold should be chosen with the policy context in mind: ",
          "5 pp is a common rule of thumb for distinguishing meaningful from random ",
          "variation, but for some interventions (e.g., severe deficiency targeting), ",
          "even smaller divergences may matter."
        )
      )
    ),

    card(
      card_header("Methodological bridge: integrating with GBD"),
      card_body(
        p("This dashboard's framework is designed to be complementary to GBD, not ",
          "competing. Two integration pathways are envisioned:"),
        tags$ol(
          tags$li(strong("Pipeline as input to GBD: "),
                  "Subnational predictions could provide finer-grained ",
                  "biomarker-grounded data points to refine GBD's country-level ",
                  "estimates, particularly for countries where direct survey data is sparse."),
          tags$li(strong("GBD as prior for pipeline: "),
                  "GBD estimates could serve as a national-level Bayesian prior, ",
                  "with subnational predictions adjusting away from the prior where ",
                  "evidence supports doing so. This would tighten uncertainty bounds ",
                  "and ensure cross-country comparability.")
        ),
        p("Engaging IHME on either pathway requires shared data formats and ",
          "documented uncertainty quantification, both of which are in scope for ",
          "this project."),
        methods_note(
          "GBD methodology relies on a Bayesian meta-regression framework (DisMod-MR) ",
          "that integrates surveys, administrative data, and modelled covariates to produce ",
          "country-level estimates with uncertainty. Our pipeline differs by operating at ",
          "Admin-2 resolution and using machine learning ensembles (SuperLearner) on a ",
          "broader proxy predictor set. The two approaches complement each other: ",
          "GBD provides comparable cross-country ",
          "estimates with consistent methodology, while subnational modelling reveals ",
          "within-country heterogeneity. Joint use is the intended long-term integration ",
          "pathway, with the pipeline contributing local detail and GBD contributing ",
          "global comparability and temporal trends."
        )
      )
    )
  )
}


mod_gbd_compare_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$placeholder_warning <- renderUI({
      if (!isTRUE(gbd_meta$using_placeholder)) return(NULL)

      div(
        class = "alert alert-warning",
        style = "margin-bottom: 1em;",
        h5(bsicons::bs_icon("exclamation-triangle"),
           " GBD comparison data is currently placeholder",
           style = "margin-top: 0;"),
        p("The values shown in these tables are illustrative and based on rough ",
          "regional estimates from the published literature, not on actual GBD ",
          "Results Tool exports. To replace with real data:"),
        tags$ol(
          tags$li("Visit the ", tags$a(href = gbd_meta$source_url,
                                       target = "_blank",
                                       "IHME GBD Results Tool"),
                  " and download prevalence estimates for the relevant ",
                  "country-outcome-year combinations."),
          tags$li("Save the file as ", tags$code("dashboard/data-raw/gbd_estimates.csv"),
                  " with columns: country, outcome, year, gbd_prev, gbd_lo, gbd_hi."),
          tags$li("Re-run ", tags$code("dashboard/data-raw/02_gbd_placeholder.R"),
                  " to refresh the dashboard data.")
        ),
        p(em("Until then, this tab demonstrates the framework but should not be ",
              "cited or shared externally."))
      )
    })

    # Build the comparison table
    compare_data <- reactive({
      gbd <- gbd_data$data
      yr_filter <- input$year_filter %||% "survey"

      rows <- list()
      for (ctry_label in unique(gbd$country)) {
        ctry_key <- names(meta$countries)[meta$countries == ctry_label]
        if (length(ctry_key) == 0) next

        if (yr_filter == "survey") {
          target_year <- meta$survey_years[ctry_key]
        } else {
          target_year <- 2023
        }

        gbd_subset <- gbd[gbd$country == ctry_label &
                            gbd$year == target_year, , drop = FALSE]

        for (oc in unique(gbd_subset$outcome)) {
          gbd_row <- gbd_subset[gbd_subset$outcome == oc, , drop = FALSE]
          if (nrow(gbd_row) == 0) next

          # Pipeline estimate from admin2_pred + admin2_pop
          df <- get_country_admin2(ctry_key, oc,
                                    admin2_bnds, admin2_pred, admin2_pop)
          if (is.null(df)) next
          natl <- national_aggregate(as.data.frame(df))

          rows[[paste(ctry_label, oc, target_year, sep = "_")]] <- data.frame(
            Country     = ctry_label,
            Outcome     = meta$outcome_labels[oc],
            Year        = target_year,
            Pipeline    = natl$pred_prev_natl,
            GBD         = gbd_row$gbd_prev[1],
            GBD_lo      = gbd_row$gbd_lo[1],
            GBD_hi      = gbd_row$gbd_hi[1],
            Diff_pp     = (natl$pred_prev_natl - gbd_row$gbd_prev[1]) * 100,
            stringsAsFactors = FALSE
          )
        }
      }
      bind_rows(rows)
    })

    output$compare_table <- renderReactable({
      df <- compare_data()
      req(nrow(df) > 0)

      reactable(
        df,
        compact = TRUE, striped = TRUE,
        searchable = TRUE, defaultPageSize = 15,
        columns = list(
          Country  = colDef(name = "Country"),
          Outcome  = colDef(name = "Outcome", minWidth = 200),
          Year     = colDef(name = "Year"),
          Pipeline = colDef(name = "Pipeline est.",
                             format = colFormat(percent = TRUE, digits = 1)),
          GBD      = colDef(name = "GBD est.",
                             format = colFormat(percent = TRUE, digits = 1)),
          GBD_lo   = colDef(name = "GBD low",
                             format = colFormat(percent = TRUE, digits = 1)),
          GBD_hi   = colDef(name = "GBD high",
                             format = colFormat(percent = TRUE, digits = 1)),
          Diff_pp  = colDef(name = "Diff (pp)",
                             format = colFormat(digits = 1),
                             style = function(v) {
                               if (is.na(v)) return(NULL)
                               color <- if (abs(v) > 5) "#d7191c"
                                 else if (abs(v) > 2) "#fdae61"
                                 else "#1a9850"
                               list(color = color, fontWeight = "bold")
                             })
        )
      )
    })

    # Hidden burden table — districts that diverge from GBD
    output$hidden_gbd_table <- renderReactable({
      thr <- (input$hidden_threshold %||% 5) / 100
      gbd <- gbd_data$data

      rows <- list()
      for (ctry_label in unique(gbd$country)) {
        ctry_key <- names(meta$countries)[meta$countries == ctry_label]
        if (length(ctry_key) == 0) next
        target_year <- meta$survey_years[ctry_key]

        gbd_subset <- gbd[gbd$country == ctry_label &
                            gbd$year == target_year, , drop = FALSE]

        for (oc in unique(gbd_subset$outcome)) {
          gbd_row <- gbd_subset[gbd_subset$outcome == oc, , drop = FALSE]
          if (nrow(gbd_row) == 0) next

          df <- get_country_admin2(ctry_key, oc,
                                    admin2_bnds, admin2_pred, admin2_pop)
          if (is.null(df)) next
          df_clean <- as.data.frame(df)
          df_clean <- df_clean[!is.na(df_clean$pred_prev) &
                                 !is.na(df_clean$population), , drop = FALSE]

          gbd_natl <- gbd_row$gbd_prev[1]
          diverging <- abs(df_clean$pred_prev - gbd_natl) > thr
          n_div <- sum(diverging, na.rm = TRUE)
          pop_div <- sum(df_clean$population[diverging], na.rm = TRUE)
          total_pop <- sum(df_clean$population, na.rm = TRUE)

          rows[[paste(ctry_label, oc, sep = "_")]] <- data.frame(
            Country = ctry_label,
            Outcome = meta$outcome_labels[oc],
            GBD_natl = gbd_natl,
            N_diverging = n_div,
            Pop_diverging = pop_div,
            Pct = pop_div / total_pop,
            stringsAsFactors = FALSE
          )
        }
      }
      df <- bind_rows(rows)
      req(nrow(df) > 0)

      reactable(
        df,
        compact = TRUE, striped = TRUE,
        searchable = TRUE, defaultPageSize = 15,
        columns = list(
          Country       = colDef(name = "Country"),
          Outcome       = colDef(name = "Outcome", minWidth = 200),
          GBD_natl      = colDef(name = "GBD national",
                                  format = colFormat(percent = TRUE, digits = 1)),
          N_diverging   = colDef(name = "Districts diverging"),
          Pop_diverging = colDef(name = "Pop. diverging",
                                  format = colFormat(separators = TRUE, digits = 0)),
          Pct           = colDef(name = "% pop. mischaracterized",
                                  format = colFormat(percent = TRUE, digits = 1),
                                  style = function(v) {
                                    if (is.na(v)) return(NULL)
                                    color <- if (v > 0.4) "#d7191c"
                                      else if (v > 0.2) "#fdae61"
                                      else "#1a9850"
                                    list(color = color, fontWeight = "bold")
                                  })
        )
      )
    })

  })
}
