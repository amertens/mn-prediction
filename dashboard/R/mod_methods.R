# =============================================================================
# Module: Methods
# =============================================================================
# Plain-language account of the model, how it is tested and how it does, the
# surveys, and the limits. The audit history and the biomarker notes are folded
# away. Numbers come from the bundles.

mod_methods_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_header("The model"),
         card_body(
           p(sprintf(paste("Every district is described by %d public data layers in %d groups. Inside each country, each layer is",
                           "replaced by the district's rank on it, so that surveys on different scales can be pooled. Each group is",
                           "condensed to a few summary scores (principal components covering 80 percent of its variation, at most 12),",
                           "each score is weighted by how strongly it tracked deficiency in the training districts, and the weighted",
                           "scores are added up. Nothing is tuned; there is no setting a modeller chooses."), Q$n_predictors, Q$n_domains)),
           p("For the maps, the model is fitted on all of a country's surveyed districts and applied to every district. Because",
             " the score is a weighted sum, each district's score splits exactly into the part each layer contributes, which is",
             " what the district profiles and What drives the estimate show. The planning prevalence anchors the ranking to the",
             " survey's national figure: the order is the model's, the level is the survey's."))),
    card(card_header("How it is tested, and how it does"),
         card_body(
           tags$ul(
             tags$li("Three tests, reported separately: a district hidden inside a surveyed country, a whole region hidden, and a",
                     " whole country hidden."),
             tags$li("Every comparison method sees exactly what the model sees. The survey's own regional average is computed without",
                     " the scored district, and a neighbour smoother runs on the same splits."),
             tags$li(sprintf(paste("Every random split is repeated ten times, and chance is measured by shuffling the outcome: a ranking",
                                   "with no information reaches %s across districts."), fmt_num(Q$null_d))),
             tags$li(sprintf(paste("The survey's own noise is reported. Most districts hold one or two survey clusters, so even a perfect",
                                   "predictor could reach only about %s on prevalence and %s on the biomarker level."),
                             fmt_num(Q$ceiling_prev), fmt_num(Q$ceiling_level)))
           ),
           reactableOutput(ns("perf")),
           methods_note("Mean ranking accuracy over country and outcome combinations, biomarker level, by test and method. Ranking",
                        " accuracy is the Spearman correlation between the model's order of districts and the survey's (0 = no",
                        " better than guessing, 1 = a perfect match). Differences under 0.03 are ties."))),
    card(card_header("The surveys"),
         card_body(reactableOutput(ns("outcomes")),
                   p(style = "font-size:0.9em; color:#555;", "Vitamin A is measured by retinol-binding protein, converted to retinol",
                     " with each survey's own calibration before the 0.70 micromol per litre cut-off; iron is inflammation-adjusted",
                     " by each survey's method. No variable from the biomarker survey itself is used as a predictor, and any",
                     " household-survey indicator that names a target nutrient is excluded. The full list of data layers, with",
                     " sources and weights, is the Predictor catalogue."))),
    card(card_header("Limits"),
         card_body(
           tags$ul(
             tags$li("The model ranks. For a publishable prevalence figure in a surveyed district, a geostatistical model with an",
                     " interval is the better tool."),
             tags$li("The ranking crosses borders; the prevalence level does not. A country with no survey gets a ranking and needs",
                     " one national figure to turn it into percentages."),
             tags$li("Four countries, three of them West African, bound the result for new countries. The next survey added is the",
                     " real test."),
             tags$li("Where a survey exists, smoothing between neighbouring districts does almost as well; the model's value is",
                     " where there is no survey."),
             tags$li("Rare outcomes (women's vitamin A is under 3 percent everywhere) and Malawi-only zinc have low ceilings or no",
                     " cross-border test.")
           ))),
    accordion(
      id = ns("more"), open = FALSE,
      accordion_panel("What changed since the earlier version", icon = bsicons::bs_icon("clock-history"),
        p("An audit of the first version of this work found that evaluation choices had decided its headlines, in both directions:"),
        tags$ul(
          tags$li("A within-country accuracy of 0.06 was one random split; the replicated median was 0.22, and the corrected index reaches ",
                  fmt_num(Q$infill), "."),
          tags$li("A baseline the models appeared to lose to had read the scored district's own respondents; recomputed without them, it loses to the index."),
          tags$li("A regional transport headline of 0.50 to 0.56 rested on three of four countries because of a spelling mismatch in a population file; it is withdrawn."),
          tags$li("A reliability ceiling used to argue that districts were unresolvable was biased low by a factor of about five."),
          tags$li("The anchoring gain shown on an earlier version of this dashboard (0.16 to 0.41) did not survive a matched control and is withdrawn."),
          tags$li("Person-level prediction of who is deficient, the earlier headline estimator, has an AUC of about ", fmt_num(Q$il_auc),
                  " under honest folds and is no longer offered.")
        ),
        p("Everything on this dashboard comes from the corrected protocol. The findings notes, the protocol scripts and the",
          " pre-registration for the next countries are in the project repository.")),
      accordion_panel("Notes on each biomarker", icon = bsicons::bs_icon("droplet"),
        tags$ul(lapply(names(biomarker_caveats), function(k) tags$li(strong(meta$outcome_labels[[k]], ": "), biomarker_caveats[[k]]))))
    )
  )
}

mod_methods_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$perf <- renderReactable({
      B <- EV$benchmarks_summary; validate(need(!is.null(B), "Benchmark summary not built."))
      d <- B[B$target == "level" & B$arm %in% names(arm_label), ] |> select(arm, estimand, mean_spearman) |>
        pivot_wider(names_from = estimand, values_from = mean_spearman)
      d$Method <- arm_label[d$arm]; d <- d[order(-d$infill), ]
      t <- data.frame(Method = d$Method, `District hidden` = round(d$infill, 2), `Region hidden` = round(d$region, 2),
                      `Country hidden` = round(d$country, 2), check.names = FALSE)
      t <- rbind(t, data.frame(Method = "Chance (a shuffled outcome, 95th percentile)", `District hidden` = NA, `Region hidden` = NA,
                               `Country hidden` = round(Q$null_d, 2), check.names = FALSE))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE,
                defaultColDef = colDef(format = colFormat(digits = 2), na = ""))
    })
    output$outcomes <- renderReactable({
      t <- data.frame(
        Survey = c("The Gambia 2018", "Ghana 2017", "Sierra Leone 2013", "Malawi 2015 to 2016"),
        Fieldwork = c("January to April 2018", "April to June 2017", "November to December 2013", "December 2015 to February 2016"),
        Districts = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) g1(idx_national$n_districts[idx_national$country == c]), numeric(1)),
        Surveyed = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) g1(idx_national$n_surveyed[idx_national$country == c]), numeric(1)),
        Outcomes = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) paste(outcome_short[idx_national$outcome[idx_national$country == c]], collapse = "; "), character(1)),
        stringsAsFactors = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE, columns = list(Outcomes = colDef(minWidth = 300)))
    })
  })
}
