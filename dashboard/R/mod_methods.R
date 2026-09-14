# =============================================================================
# Module: Methods
# =============================================================================
# Plain-language account of the data, the model, the protocol it is scored
# under, what the audit changed, and the limits. Numbers come from the bundles.

mod_methods_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_header("The model"),
         card_body(
           p(sprintf(paste("To every district we attach %d public data layers in %d groups. Inside each country every layer is",
                           "replaced by its rank across districts, so that surveys on different scales can be pooled. Each group is",
                           "reduced to a few summary axes (principal components to 80 percent of its variance, at most 12), each axis is",
                           "weighted by how strongly it tracked deficiency in the training districts (a Fisher-z of its rank",
                           "correlation), and the weighted axes are summed. Nothing is tuned; the index has no setting a modeller",
                           "chooses. The sum is rescaled to the training districts' mean and spread and ranked."), Q$n_predictors, Q$n_domains)),
           p("Because the index is linear in the ranked layers, its weights project back exactly onto the original",
             " layers, which is what What drives the estimate and the district decompositions show. The same projection",
             " defines the twenty-layer composite a country could compute for itself."),
           p("For the maps, the index is fitted on all of a country's surveyed districts and applied to every district. The",
             " planning prevalence anchors the resulting ranking to the survey's national prevalence: the order is the model's,",
             " the level the survey's."))),
    card(card_header("How it is tested"),
         card_body(
           tags$ol(
             tags$li(strong("Three questions, answered separately. "), "A district hidden inside a surveyed country (five folds over districts,",
                     " ten random draws); a whole region hidden; a whole country hidden. They have different answers and are reported as such."),
             tags$li(strong("Every comparator sees what the model sees. "), "The survey's own regional average is computed without the scored",
                     " district; a covariate-free neighbour smoother runs on the same folds."),
             tags$li(strong("Every random split is repeated. "), "One split of 14 to 87 districts can move a result by 0.1; every number is the",
                     " median of ten draws."),
             tags$li(strong("Chance is measured. "), sprintf("Shuffling the outcome across countries shows what a ranking reaches with no information: %s across districts, %s across regions.",
                                                           fmt_num(Q$null_d), fmt_num(Q$null_a1))),
             tags$li(strong("The target's own reliability is reported. "), sprintf("Most districts hold one or two survey clusters; a variance model gives the ceiling a perfect predictor could reach (about %s on prevalence, %s on the biomarker level).",
                                                                                   fmt_num(Q$ceiling_prev), fmt_num(Q$ceiling_level)))
           ),
           p("Districts are weighted by an effective sample size from measured design effects (median 2.4). Two targets are carried",
             " for every country and outcome: the survey-weighted district prevalence, and the district mean of the log biomarker",
             " concentration, which keeps information a cut-off discards."))),
    card(card_header("Performance, in one table"),
         card_body(reactableOutput(ns("perf")),
                   methods_note("Mean ranking accuracy over country-outcome combinations, biomarker level, by test and method.",
                                " Differences under 0.03 are ties. The chance level is the 95th percentile of the permutation null."))),
    card(card_header("The surveys and outcomes"),
         card_body(reactableOutput(ns("outcomes")),
                   p(style = "font-size:0.9em; color:#555;", "Vitamin A uses each survey's own calibration of retinol-binding protein to",
                     " retinol before the 0.70 micromol per litre cut-off; iron is inflammation-adjusted by each survey's method. Fieldwork",
                     " dates were recovered for every cluster so that time-varying layers match the survey year."))),
    card(card_header("The data"),
         card_body(reactableOutput(ns("sources")),
                   p(style = "font-size:0.9em; color:#555;", "Any household-survey indicator that names a target nutrient or its biomarker",
                     " is excluded as leakage. National values enter pooled models only. The full list, with definitions, is the",
                     " Predictor catalogue."))),
    card(card_header("What changed since the earlier version"),
         card_body(
           p("An audit of the first version of this work found that evaluation choices had decided its headlines, in both directions:"),
           tags$ul(
             tags$li("A within-country accuracy of 0.06 was one random split; the replicated median was 0.22, and the corrected index reaches ",
                     fmt_num(Q$infill), "."),
             tags$li("A baseline the models appeared to lose to had read the scored district's own respondents; recomputed without them it loses to the index."),
             tags$li("A regional transport headline of 0.50 to 0.56 rested on three of four countries because of a spelling mismatch in a population file; it is withdrawn."),
             tags$li("A reliability ceiling used to argue that districts were unresolvable was biased low by a factor of about five."),
             tags$li("The anchoring gain shown on an earlier version of this dashboard (0.16 to 0.41) did not survive a matched control and is withdrawn."),
             tags$li("Person-level prediction of who is deficient, the earlier headline estimator, has an AUC of about ", fmt_num(Q$il_auc),
                     " under honest folds and is no longer offered.")
           ),
           p("Everything on this dashboard comes from the corrected protocol. The findings notes, the protocol scripts and the",
             " pre-registration for the next countries are in the project repository."))),
    card(card_header("Limits"),
         card_body(
           tags$ul(
             tags$li("The model ranks. For a publishable prevalence figure in a surveyed district, the geostatistical model is the better tool."),
             tags$li("Levels do not cross borders; rankings do. A country with no survey gets a ranking and needs one national number to turn it into prevalence."),
             tags$li("Four countries, three of them West African, bound the transport result. The next survey added is the real test."),
             tags$li("Where a survey exists, smoothing between neighbours does almost as well; the model's value is where there is no survey."),
             tags$li("Household-survey inputs help inside a country and hurt in a new one; the cross-border model leans on the physical environment."),
             tags$li("Rare outcomes (women's vitamin A under 3 percent everywhere) and Malawi-only zinc have low ceilings or no cross-border test.")
           ),
           h6("Biomarker notes"),
           tags$ul(lapply(names(biomarker_caveats), function(k) tags$li(strong(meta$outcome_labels[[k]], ": "), biomarker_caveats[[k]])))))
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
      t <- rbind(t, data.frame(Method = "Chance level (permutation null, 95th percentile)", `District hidden` = NA, `Region hidden` = NA,
                               `Country hidden` = round(Q$null_d, 2), check.names = FALSE))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE)
    })
    output$outcomes <- renderReactable({
      t <- data.frame(
        Survey = c("The Gambia 2018", "Ghana 2017", "Sierra Leone 2013", "Malawi 2015 to 2016"),
        Fieldwork = c("January to April 2018", "April to June 2017", "November to December 2013", "December 2015 to February 2016"),
        Districts = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) g1(idx_national$n_districts[idx_national$country == c]), numeric(1)),
        Surveyed = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) g1(idx_national$n_surveyed[idx_national$country == c]), numeric(1)),
        Outcomes = vapply(c("Gambia", "Ghana", "Sierra Leone", "Malawi"), function(c) paste(outcome_short[idx_national$outcome[idx_national$country == c]], collapse = "; "), character(1)),
        stringsAsFactors = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, columns = list(Outcomes = colDef(minWidth = 300)))
    })
    output$sources <- renderReactable({
      S <- CAT$sources; validate(need(!is.null(S), "Source table not built."))
      t <- S[order(-S$n_columns), ]
      names(t) <- c("Source", "Layers", "Domains", "With a definition")
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, columns = list(Domains = colDef(minWidth = 320)))
    })
  })
}
