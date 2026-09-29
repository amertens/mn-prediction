# =============================================================================
# Module: Methods
# =============================================================================
# Plain-language account of the model, how it is tested and how it does, the
# surveys, the limits, and where the survey-planning ideas come from. Project
# history and the notes on each biomarker are folded away. Numbers come from
# the bundles.

mod_methods_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_header("The model"),
         card_body(
           p(sprintf(paste("Each district is described by %d public data layers in %d groups, and the model uses %s of them.",
                           "It leaves out the DHS survey summaries and layers that are the same across a whole country. Within",
                           "each country, each layer is replaced by the district's rank on that layer, so that surveys measured",
                           "on different scales can be combined. Each group of layers is summarised by a few principal",
                           "components (enough to cover 80%% of the group's variation, at most 12). Each component is weighted by",
                           "how strongly it followed deficiency in the training districts, and the weighted components are added",
                           "together. No setting was adjusted to improve the results."),
                     Q$n_predictors, Q$n_domains, if (is.finite(Q$n_in_model)) Q$n_in_model else "383")),
           p("For the maps, the model is built on all of a country's surveyed districts and applied to every district.",
             " Because the score is a weighted sum, each district's score can be split exactly into the contribution of",
             " each data layer, which is what the district profiles show. The estimated prevalence combines the ranking with",
             " the national survey figure: the order of districts comes from the model and the overall level from the survey."))),
    card(card_header("How it is tested, and how it does"),
         card_body(
           tags$ul(
             tags$li("Three tests, reported separately: a district hidden inside a surveyed country, a whole region hidden,",
                     " and a whole country left out."),
             tags$li("Every comparison method gets the same information as the model. For example, the survey's regional",
                     " average leaves out the district being predicted, and neighbour averaging uses the same splits."),
             tags$li(sprintf("Every random split is repeated ten times. A random ranking, made by shuffling the outcome, stays below %s in 95%% of tries.",
                             fmt_num(Q$null_d))),
             tags$li(sprintf(paste("The survey's own noise is measured. Most districts have one or two survey clusters, so even a",
                                   "perfect model could reach only about %s on prevalence and %s on biomarker levels."),
                             fmt_num(Q$ceiling_prev), fmt_num(Q$ceiling_level))),
             tags$li(sprintf(paste("Percentage ranges are checked on districts the model had not seen: 90%% ranges contained the",
                                   "survey figure %s of the time. Ranges made by re-estimating the model on resampled districts",
                                   "contained it only %s of the time, so they are used only to show how firmly districts are ranked."),
                             pc(Q$cal_cov), pc(Q$stabprev_cov)))
           ),
           reactableOutput(ns("perf")),
           methods_note("Average ranking accuracy across country-outcome pairs (biomarker levels), by test and method. Ranking",
                        " accuracy is the Spearman correlation between a method's order of districts and the survey's (0 = no",
                        " better than random, 1 = the same order). Differences under 0.03 are ties."))),
    card(card_header("The surveys"),
         card_body(reactableOutput(ns("outcomes")),
                   p(style = "font-size:0.9em; color:#555;",
                     "Vitamin A deficiency is retinol-binding protein, adjusted for inflammation with the BRINDA method and",
                     " converted to a retinol value with each survey's published conversion, below 0.70 micromoles per litre.",
                     " Iron deficiency is each survey's inflammation-adjusted ferritin below 12 micrograms per litre (children)",
                     " or 15 (women), not iron-deficiency anaemia. Children are 6 to 59 months old; women are 15 to 49 and not",
                     " pregnant, as in the survey reports. Sierra Leone's child data as supplied cover only the children found",
                     " to be anaemic. No variable from the biomarker surveys is used as a predictor, and household-survey",
                     " indicators that name a target nutrient are excluded. The full list of data layers is in the data",
                     " layer catalogue."))),
    card(card_header("Limits"),
         card_body(
           tags$ul(
             tags$li("The model ranks districts. For a published prevalence figure in a surveyed district, a geostatistical",
                     " model with an interval is the better tool, and the national prevalence, which is the main result of a",
                     " biomarker survey, still needs a probability sample."),
             tags$li("The ranking carries over to a new country; prevalence levels do not. A country without a survey gets",
                     " a ranking and needs at least a national figure to turn it into percentages."),
             tags$li("Four training countries, three of them in West Africa, limit what can be said about new countries.",
                     " The next national survey will be the real test."),
             tags$li("Where a survey exists, averaging neighbouring districts does almost as well as the model; the model",
                     " adds most where there is no survey."),
             tags$li("Rare outcomes (vitamin A deficiency in women is under 3% everywhere) and outcomes measured in Malawi",
                     " only (zinc, selenium, iodine) either have little to rank or cannot be checked in another country.")
           ))),
    card(card_header("Where the survey-planning ideas come from"),
         card_body(
           p(style = "font-size:0.9em;",
             "The Plan a survey tools adapt ideas used in survey design for disease control programmes:"),
           tags$ul(style = "font-size:0.85em; color:#555;",
             tags$li("Choosing the next sampling locations where predictions are uncertain or close to a decision",
                     " threshold (adaptive geostatistical design): Chipeta, Terlouw, Phiri and Diggle (2016), Spatial",
                     " Statistics; used in repeated malaria surveys in Malawi by Kabaghe and colleagues (2017), PLOS ONE."),
             tags$li("Choosing batches of locations to find areas above a prevalence threshold with fewer samples:",
                     " Andrade-Pacheco and colleagues (2020), Scientific Reports."),
             tags$li("Using a model to classify areas against elimination thresholds in neglected tropical disease",
                     " programmes: Fronterre, Amoah, Giorgi, Stanton and Diggle (2020), Journal of Infectious Diseases;",
                     " Diggle and colleagues (2021), Transactions of the Royal Society of Tropical Medicine and Hygiene;",
                     " Amoah and colleagues (2022), International Journal of Epidemiology."),
             tags$li("This dashboard applies their logic for choosing districts to micronutrient surveys. It does not",
                     " implement their cluster-level geostatistical methods, and its planning tools have been tested on",
                     " four existing surveys but not in the field.")))),
    accordion(
      id = ns("more"), open = FALSE,
      accordion_panel("What changed since the earlier version", icon = bsicons::bs_icon("clock-history"),
        p("An audit of the first version of this work found that choices in how results were evaluated had driven its",
          " main findings, in both directions:"),
        tags$ul(
          tags$li("A within-country accuracy of 0.06 came from one random split; repeated splits gave 0.22, and the corrected",
                  " model reaches ", fmt_num(Q$infill), "."),
          tags$li("A comparison method that the models appeared to lose to had used the predicted district's own survey",
                  " respondents. Recomputed without them, it does worse than the model."),
          tags$li("A regional result of 0.50 to 0.56 for new countries rested on three of the four countries, because a",
                  " spelling mismatch in a population file dropped one. It was withdrawn."),
          tags$li("An estimate of the best achievable score, used to argue that districts could not be told apart, was about",
                  " five times too low."),
          tags$li("A gain from adding a national survey figure (0.16 to 0.41) shown on an earlier version of this dashboard",
                  " did not hold up against a matched comparison and was withdrawn."),
          tags$li("Predicting which individual people are deficient, the main method of the earlier version, has an AUC of",
                  " about ", fmt_num(Q$il_auc), " when tested properly, and is no longer offered.")
        ),
        p("Everything on this dashboard comes from the corrected analysis. The findings notes, the analysis scripts and",
          " the pre-registration for the next countries are in the project repository.")),
      accordion_panel("Notes on each biomarker", icon = bsicons::bs_icon("droplet"),
        tags$ul(lapply(names(biomarker_caveats), function(k) tags$li(strong(meta$outcome_labels[[k]] %||% k, ": "), biomarker_caveats[[k]]))))
    )
  )
}

mod_methods_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$perf <- renderReactable({
      B <- EV$benchmarks_summary; validate(need(!is.null(B), "Results summary not built."))
      d <- B[B$target == "level" & B$arm %in% names(arm_label), ] |> select(arm, estimand, mean_spearman) |>
        pivot_wider(names_from = estimand, values_from = mean_spearman)
      d$Method <- arm_label[d$arm]; d <- d[order(-d$infill), ]
      t <- data.frame(Method = d$Method, `District hidden` = round(d$infill, 2), `Region hidden` = round(d$region, 2),
                      `Country left out` = round(d$country, 2), check.names = FALSE)
      t <- rbind(t, data.frame(Method = "Random ranking (95% of tries fall below)", `District hidden` = NA, `Region hidden` = NA,
                               `Country left out` = round(Q$null_d, 2), check.names = FALSE))
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE,
                defaultColDef = colDef(format = colFormat(digits = 2), na = ""))
    })
    output$outcomes <- renderReactable({
      cn <- c("Gambia", "Ghana", "Sierra Leone", "Malawi")
      fw <- vapply(cn, function(c) gsub(" – ", " to ", meta$fieldwork_dates[[c]] %||% ""), character(1))
      t <- data.frame(
        Survey = c("The Gambia 2018", "Ghana 2017", "Sierra Leone 2013", "Malawi 2015 to 2016"),
        Fieldwork = unname(fw),
        Districts = vapply(cn, function(c) g1(idx_national$n_districts[idx_national$country == c]), numeric(1)),
        Surveyed = vapply(cn, function(c) g1(idx_national$n_surveyed[idx_national$country == c]), numeric(1)),
        Outcomes = vapply(cn, function(c) paste(outcome_short[idx_national$outcome[idx_national$country == c]], collapse = "; "), character(1)),
        stringsAsFactors = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE, columns = list(Outcomes = colDef(minWidth = 300)))
    })
  })
}
