# =============================================================================
# Module: Methods Explainer
# =============================================================================
# Plain-language description of the modelling approach, validation, and
# limitations. Targeted at technical stakeholders (ministry of health
# advisors, academic collaborators) who need to assess credibility.

mod_methods_ui <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = 12,

    card(
      card_header("How these estimates are produced"),
      card_body(
        h5("In short"),
        p("We trained machine learning models on individual-level biomarker data ",
          "from nationally representative micronutrient surveys in four ",
          "sub-Saharan African countries, linked to a wide range of routinely ",
          "available proxy data (satellite imagery, climate reanalysis, modeled ",
          "disease burden, food security indicators, household survey ",
          "indicators, and food prices). The trained models then predict ",
          "deficiency prevalence at the district level, including for districts ",
          "without direct survey coverage."),

        h5("Why proxy-only?"),
        p("The models exclude all variables from the original biomarker survey ",
          "as predictors. This is intentional: it ensures the predictions can ",
          "be applied to areas where no biomarker survey exists, using only ",
          "publicly or routinely available data sources."),

        h5("How predictions are validated"),
        p("Models are evaluated using two complementary cross-validation ",
          "designs:"),
        tags$ul(
          tags$li(strong("Cluster-blocked cross-validation: "),
                  "individuals from the same survey cluster are assigned to ",
                  "the same fold, so the model is always tested on individuals ",
                  "spatially separated from those used to train it. This ",
                  "guards against optimistic estimates that would result from ",
                  "splitting clusters across folds."),
          tags$li(strong("Leave-one-country-out (LOCO) cross-validation: "),
                  "the model is trained on three countries and tested on the ",
                  "fourth. This estimates how well the model would perform if ",
                  "applied to a new country with no biomarker data — a key ",
                  "consideration for scaling up.")
        )
      )
    ),

    card(
      card_header("Model performance summary"),
      card_body(
        p("Cross-validated discrimination (ROC-AUC) for binary deficiency ",
          "models. AUC ranges from 0.5 (no discrimination, random) to 1.0 ",
          "(perfect discrimination)."),
        reactableOutput(ns("perf_table")),
        methods_note(
          "ROC-AUC is computed on cluster-blocked cross-validated predictions, ",
          "so it reflects performance on held-out clusters (not held-out individuals ",
          "within clusters). PR-AUC (precision-recall area under curve) is more ",
          "informative than ROC-AUC for rare outcomes — when prevalence is below ",
          "10%, ROC-AUC can appear high simply because the model correctly ",
          "predicts most non-events. ",
          tags$br(), tags$br(),
          "Brier skill score compares the model's mean squared prediction error ",
          "against a baseline that always predicts the prevalence. Positive values ",
          "mean the model improves on the baseline; values near zero mean the ",
          "model adds little signal beyond knowing the overall prevalence. ",
          tags$br(), tags$br(),
          "Sample sizes (N) reflect the number of individuals with non-missing ",
          "outcome and at least some predictor data after preprocessing. Models ",
          "with N below approximately 500 may be unstable; treat their estimates ",
          "as preliminary."
        )
      )
    ),

    card(
      card_header("Confidence intervals"),
      card_body(
        h5("Conformal prediction intervals"),
        p("Uncertainty around predicted prevalence is quantified using ",
          tags$em("split-conformal prediction"), " — a distribution-free ",
          "approach that uses cross-validation residuals to construct ",
          "intervals with provable coverage guarantees, without assuming any ",
          "particular noise distribution."),
        p("Our intervals target 95% marginal coverage, meaning that if we ",
          "repeat this exercise across many districts, approximately 95% of ",
          "the constructed intervals should contain the true prevalence. ",
          "Intervals are computed at the Admin-1 level and broadcast to the ",
          "constituent Admin-2 districts."),
        p("Wider intervals indicate higher predictor variability or sparser ",
          "training data for that area; narrow intervals indicate a more ",
          "confident prediction.")
      )
    ),

    card(
      card_header("Limitations"),
      card_body(
        tags$ul(
          tags$li(strong("Cross-country generalizability is limited. "),
                  "When models are trained in three countries and applied to a ",
                  "fourth, performance drops substantially (LOCO AUC 0.50–0.73). ",
                  "Predictions for unsurveyed countries should be treated as ",
                  "screening estimates rather than definitive."),
          tags$li(strong("Survey-vintage assumption. "),
                  "Each prediction reflects conditions at the time of the original ",
                  "biomarker survey. We assume the relationship between proxy ",
                  "indicators and deficiency status is stable over time, which ",
                  "may not hold under rapid food system or climate change."),
          tags$li(strong("Population denominators are uncertain. "),
                  "Population-affected counts use WorldPop estimates, which ",
                  "themselves carry uncertainty (typically 5–15% at Admin-2 level). ",
                  "Counts should be interpreted as planning estimates."),
          tags$li(strong("Some outcomes are very rare. "),
                  "When deficiency prevalence is below approximately 5%, models ",
                  "have few positive examples to learn from and predictions may ",
                  "be unstable.")
        )
      )
    ),

    card(
      card_header("Data sources"),
      card_body(
        reactableOutput(ns("sources_table")),
        methods_note(
          "All external data sources are publicly available or accessible through ",
          "common APIs. Variables from each source are linked to individual survey ",
          "respondents through their cluster's geographic location (latitude/longitude) ",
          "or by matching cluster Admin-2 districts to externally aggregated values. ",
          tags$br(), tags$br(),
          "Where multiple temporal versions of a source exist, we select the ",
          "version closest to the survey fieldwork year. For climate variables ",
          "(rainfall, temperature) the survey year is used directly; for slowly ",
          "changing indicators (soil properties) the most recent global version is used."
        )
      )
    )
  )
}


mod_methods_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    output$perf_table <- renderReactable({
      df <- cv_perf
      df <- df[df$model_type == "binary" & !is.na(df$auc), , drop = FALSE]

      df$Outcome <- meta$outcome_labels[df$outcome]
      df_show <- df[, c("country", "Outcome", "n", "auc", "brier")]
      colnames(df_show) <- c("Country", "Outcome", "N", "AUC", "Brier")

      reactable(
        df_show,
        compact = TRUE, striped = TRUE,
        searchable = TRUE, defaultPageSize = 12,
        columns = list(
          AUC = colDef(name = "ROC-AUC",
                       format = colFormat(digits = 2),
                       style = function(value) {
                         if (is.na(value)) return(NULL)
                         color <- if (value >= 0.80) "#1a9850"
                           else if (value >= 0.70) "#91cf60"
                           else if (value >= 0.60) "#fdae61"
                           else "#d7191c"
                         list(color = color, fontWeight = "bold")
                       }),
          Brier = colDef(format = colFormat(digits = 3))
        )
      )
    })

    output$sources_table <- renderReactable({
      sources <- data.frame(
        Source = c(
          "Biomarker surveys (GMNS, GMS, SLMS, MNS)",
          "CHIRPS rainfall",
          "WorldPop population density and age structure",
          "VIIRS / NASA Black Marble nighttime lights",
          "Malaria Atlas Project (MAP)",
          "ISRIC SoilGrids",
          "Global Data Lab Subnational HDI",
          "WFP HungerMap (food security)",
          "WFP food prices",
          "IPC / Cadre Harmonisé food security",
          "ACLED conflict events",
          "HarvestStat Africa crop statistics",
          "DHS Admin-2 indicators (smoothed)",
          "MODIS land surface temperature, NDVI",
          "FLDAS climate reanalysis",
          "GHSL settlement, GPW grasslands"
        ),
        Domain = c(
          "Outcome",
          "Climate",
          "Demography",
          "Infrastructure",
          "Disease burden",
          "Soil",
          "Development",
          "Food security",
          "Food prices",
          "Food security",
          "Conflict",
          "Agriculture",
          "Health system / sociodemographic",
          "Environment",
          "Climate",
          "Settlement / land use"
        ),
        Resolution = c(
          "Individual (cluster-georef.)",
          "5 km, monthly",
          "100 m, annual",
          "500 m, annual",
          "5 km, annual",
          "250 m, static",
          "Admin-1, annual",
          "Admin-1, current",
          "Market, monthly",
          "Admin-2, annual",
          "Event-level",
          "Admin-1, annual",
          "Admin-2, smoothed",
          "500 m – 1 km, annual",
          "10 km, monthly",
          "100 m – 1 km, annual"
        ),
        stringsAsFactors = FALSE
      )

      reactable(sources, compact = TRUE, striped = TRUE,
               defaultPageSize = 20)
    })

  })
}
