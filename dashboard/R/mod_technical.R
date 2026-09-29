# =============================================================================
# Module: Technical notes
# =============================================================================
# The appendix of the concise dashboard: every explanation that the pages
# themselves now give in one line, in full, grouped by page. Also the method,
# the tests, the surveys, the limits, the sources of the survey-planning
# ideas, the project history and the notes on each biomarker.

tech_section <- function(title, ...) accordion_panel(title, ...)

mod_technical_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_body(
      p(class = "lead", style = "margin-bottom:0;", "The detail behind every page, for readers who want it."),
      p(style = "color:#555; font-size:0.9em; margin-bottom:0;",
        "Each section below matches a page of the dashboard. Open the one you need."))),
    accordion(
      id = ns("tech"), open = FALSE, multiple = TRUE,

      tech_section("How to read the numbers",
        tags$dl(
          tags$dt("Ranking accuracy"),
          tags$dd("How closely the model's order of districts matches the survey's order, from 0 (no better than a random",
                  " order) to 1 (the same order). It is a Spearman correlation. Differences smaller than 0.03 are ties."),
          tags$dt("Random ranking"),
          tags$dd(sprintf("The score a random order of districts reaches in 95%% of tries: %s for districts and %s for regions.",
                          f2(Q$null_d), f2(Q$null_a1))),
          tags$dt("Best achievable score"),
          tags$dd("The highest ranking accuracy any model could reach against the survey's district figures. It is below 1",
                  " because most districts have only one or two survey clusters, so the survey figures are themselves uncertain."),
          tags$dt("Estimated prevalence"),
          tags$dd("The model's ranking converted into a percentage using the country's national survey figure. The order of",
                  " districts comes from the model and the overall level from the survey."),
          tags$dt("Reliability of the district percentages"),
          tags$dd("How well the model's district percentages matched the survey's district figures in districts the model had",
                  " not seen (a correlation from 0 to 1): below 0.10 not informative, 0.10 to 0.30 weak, 0.30 to 0.50",
                  " moderate, above 0.50 good. Where it is not informative, every district is shown close to the national figure."),
          tags$dt("Checked 90% range"),
          tags$dd(calibrated_note()),
          tags$dt("Rank range when re-estimated"),
          tags$dd(paste0(stability_note(), ".")),
          tags$dt("Placed in the worst fifth (share of runs)"),
          tags$dd(sprintf(paste("How often a surveyed district was placed in its country's worst fifth when the model was",
                                "re-estimated 40 times with that district hidden. It points in the right direction but",
                                "overstates certainty: districts placed there in at least 80%% of runs were in the survey's worst",
                                "fifth %s of the time, and those placed there in at most 20%% of runs %s of the time, against",
                                "20%% by chance. Part of the gap is because the survey's own worst fifth is uncertain."),
                          pc(Q$cal_top), pc(Q$cal_low))),
          tags$dt("Cross-hatching on the map"),
          tags$dd("A district is cross-hatched when its rank range spans more than half of the country's list."))),

      tech_section("The model",
        p(sprintf(paste("Each district is described by %d public data layers in %d groups, including satellite imagery, climate,",
                        "soil, crops, livestock, malaria and other diseases, food prices, and summaries of public household",
                        "surveys (MICS and household budget surveys). The model uses %s of them. It leaves out the DHS survey",
                        "summaries, so that it relies only on data a country without a recent DHS survey could also collect,",
                        "and layers that are the same across a whole country. Within each country, each layer is replaced by",
                        "the district's rank on that layer, so that surveys measured on different scales can be combined. Each",
                        "group of layers is summarised by a few principal components (enough to cover 80%% of the group's",
                        "variation, at most 12). Each component is weighted by how strongly it followed deficiency in the",
                        "training districts, and the weighted components are added together. No setting was adjusted to",
                        "improve the results."),
                  Q$n_predictors, Q$n_domains, if (is.finite(Q$n_in_model)) Q$n_in_model else "383")),
        p("For the maps, the model is built on all of a country's surveyed districts and applied to every district, so",
          " for surveyed districts the map reflects their own survey data. Because the score is a weighted sum, each",
          " district's score can be split exactly into the contribution of each data layer, which is what the district",
          " profiles show. These contributions show where deficiency is likely, not what causes it.")),

      tech_section("How it is tested",
        tags$ul(
          tags$li("Three tests, reported separately: a district hidden inside a surveyed country, a whole region hidden,",
                  " and a whole country left out (the situation of a country with no survey)."),
          tags$li("Every comparison method gets the same information as the model. For example, the survey's regional",
                  " average leaves out the district being predicted, and neighbour averaging uses the same splits."),
          tags$li(sprintf("Every random split is repeated ten times. A random ranking, made by shuffling the outcome, stays below %s in 95%% of tries.",
                          fmt_num(Q$null_d))),
          tags$li(sprintf(paste("Most districts have one or two survey clusters, so even a perfect model could reach only about %s",
                                "on prevalence and %s on biomarker levels. The model reaches %s and %s. More survey clusters per",
                                "district would raise the best achievable score."),
                          fmt_num(Q$ceiling_prev), fmt_num(Q$ceiling_level), fmt_num(Q$infill_prev), fmt_num(Q$infill))),
          tags$li(sprintf(paste("Accuracy in a country left out of training rose from %s with one training survey to %s with",
                                "three, about %s for each survey added; with three points we cannot tell where it levels off."),
                          fmt_num(Q$lc1), fmt_num(Q$lc_max), fmt_num(Q$lc_step))),
          tags$li(sprintf(paste("With a country left out, the model scores %s with all data layers and %s with climate and soil",
                                "only. The climate-and-soil version was chosen after seeing these results, so its result on a",
                                "fifth country is a prediction that has not yet been tested."),
                          fmt_num(Q$tr), fmt_num(Q$cs))),
          tags$li(sprintf(paste("Percentage ranges were checked on districts the model had not seen: 90%% ranges contained the",
                                "survey figure %s of the time. Ranges made by re-estimating the model on resampled districts",
                                "contained it only %s of the time, so they are used only to show how firmly districts are ranked."),
                          pc(Q$cal_cov), pc(Q$stabprev_cov)))),
        reactableOutput(ns("perf")),
        p(style = "font-size:0.85em; color:#555;", "Average ranking accuracy across country-outcome pairs (biomarker levels), by test and method.")),

      tech_section("Compared with a geostatistical model",
        p(sprintf(paste("A geostatistical model, the kind the DHS Program uses, combines a smooth map of location with data",
                        "layers and is fitted to the locations of survey clusters. It ranks districts less accurately than this",
                        "model (%s against %s inside a surveyed country) but estimates the percentage more accurately (%s against",
                        "%s percentage points of error). It needs survey clusters inside the country, so it cannot be used for a",
                        "country without a survey. For a published prevalence in a surveyed district it is the better tool; for",
                        "deciding which districts to look at first, or for a country without a survey, the ranking is more useful."),
                  fmt_num(Q$mbg_rank), fmt_num(Q$mbg_rank_index), fmt_num(Q$mbg_err, 1), fmt_num(Q$mbg_err_index, 1)))),

      tech_section("What else was tried",
        tags$ul(
          tags$li(sprintf(paste("Machine-learning ensembles (SuperLearner, with 12 to 16 learners) at best matched this model, and",
                                "did worse when tuned to minimise squared error. With 14 to 87 districts per country, methods that",
                                "tune their own settings did not do better. Predicting which individual people are deficient reached",
                                "an AUC of %s, close to a coin toss."), fmt_num(Q$il_auc))),
          tags$li("Other ways of weighting the same data groups (ridge regression, lasso and elastic net) did slightly worse",
                  " within a country. Two variants did about 0.03 better in new countries and will be tested on the next country."),
          tags$li("Fitting the model to individual survey clusters instead of districts did not do better in any test."),
          tags$li("Adding further data (livestock, distance to water and coast, intestinal worms, updated health maps, and",
                  " food prices and temperatures at the time of fieldwork) changed accuracy by less than 0.01 each.")),
        uiOutput(ns("weights"))),

      tech_section("The six-country check",
        p(sprintf(paste("In September 2026 the climate-and-soil ranking was compared with survey results that other teams had",
                        "submitted to the WHO micronutrient database, for Zambia, Ethiopia, Sudan and Nigeria and for Pakistan",
                        "and India. None of these surveys was used to build the model. In the African countries the ranking",
                        "matched at %s for biomarker levels (above zero in %s of %s tests; a random ranking stays below %s in",
                        "95%% of tries when a country's outcomes are shuffled together) and %s for prevalence; in Pakistan and",
                        "India it reached %s for prevalence. The same climate-and-soil test inside the four training countries,",
                        "at the regional level, scores %s. A global soil layer worked as well as the Africa-only one (%s",
                        "against %s)."),
                  f2(Q$xv_level), Q$xv_level_pos, Q$xv_level_cells, f2(Q$xv_level_null), f2(Q$xv_prev), f2(Q$xv_off),
                  f2(Q$cs_a1), f2(Q$xv_sg_level), f2(Q$xv_level))),
        tags$ul(
          tags$li("These surveys report results for regions or provinces, not districts."),
          tags$li("Only rankings are compared: cut-offs and laboratory methods differ between surveys."),
          tags$li("Nigeria's results cover only six zones, so Nigeria counts only in the combined result."),
          tags$li("Vitamin A did poorly in the African surveys but was the strongest outcome in South Asia; the reason is not known."),
          tags$li("This does not replace the planned test on a new country's own survey data."))),

      tech_section("Targeting and WHO bands",
        p(sprintf(paste("A programme directed to the fifth of districts the model ranks worst would reach %s of a country's",
                        "deficient people, against %s for the survey's regional averages, %s for a random fifth and %s for the",
                        "true worst fifth. Deficiency is spread across many districts, which is why even the true worst fifth",
                        "holds under half of it. In a country left out of the model this gain disappears."),
                  fmt_pct(Q$cap_index, 0), fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_null, 0), fmt_pct(Q$cap_oracle, 0))),
        p(sprintf(paste("Using the WHO public health bands for vitamin A (under 2%%, 2%% to 10%%, 10%% to 20%%, 20%% or more), the",
                        "model places %s of districts in the correct band and %s within one band; the survey's own regional",
                        "averages do slightly better (%s and %s). Malawi's districts all fall in the lowest band and are excluded."),
                  fmt_pct(Q$band_exact, 0), fmt_pct(Q$band_w1, 0), fmt_pct(Q$band_exact_jk, 0), fmt_pct(Q$band_w1_jk, 0)))),

      tech_section("Survey planning",
        p("The planner scores each district from five parts: the model's estimated burden, how far its rank moves when the",
          " model is re-estimated, how close its estimated prevalence is to a WHO threshold, its population, and whether a",
          " previous survey reached it. The chosen aim sets how much each part counts. It helps choose districts; choosing",
          " clusters within districts, setting sample sizes and weighting the results still need a survey statistician."),
        p(sprintf(paste("Small surveys: the 'national sample + ranking' design uses a small national sample only to estimate the",
                        "national prevalence and orders districts with the ranking from a model built without that country. With",
                        "5%% of a full survey's sample, its district estimates are off by a median of %s points, against %s for a",
                        "regional survey and %s for a district survey of the same size; the district survey needs about %s of the",
                        "full sample to do as well. A district survey of any size reaches more of the deficient population, and",
                        "blending the model into a small survey's own district estimates did not improve them."),
                  fmt_num(Q$ar_a1, 1), fmt_num(Q$ar_c, 1), fmt_num(Q$ar_b, 1), fmt_pct(Q$ar_b_match, 0))),
        p(sprintf(paste("Choosing districts: a test on the four surveys, with its design written down before it was run,",
                        "treated only some surveyed districts as visited and scored the model's ranking of the rest (40 repeats,",
                        "%s country-outcome pairs). At half the districts, spreading visits across the model's ranking scored %s",
                        "against %s for random choice, and choosing by population %s. The model from the other countries alone",
                        "scores %s, which beats a model built on a third of the districts or fewer. For the national estimate,",
                        "random choice, or spread choice weighted by population, stays unbiased; visiting only the highest- and",
                        "lowest-ranked districts biases it by about 1.4 points."),
                  g1(Q$plan_cells), f2(Q$plan_spread), f2(Q$plan_random), f2(Q$plan_pps), f2(Q$plan_transport))),
        p(style = "font-size:0.88em; color:#555;",
          "Sources of the ideas: adaptive geostatistical design (Chipeta, Terlouw, Phiri and Diggle, 2016, Spatial Statistics;",
          " Kabaghe and colleagues, 2017, PLOS ONE); batch sampling to find areas above a threshold (Andrade-Pacheco and",
          " colleagues, 2020, Scientific Reports); model-based classification against thresholds in neglected tropical disease",
          " programmes (Fronterre and colleagues, 2020, Journal of Infectious Diseases; Diggle and colleagues, 2021, Transactions",
          " of the Royal Society of Tropical Medicine and Hygiene; Amoah and colleagues, 2022, International Journal of",
          " Epidemiology). This dashboard applies their logic for choosing districts, not their cluster-level methods.")),

      tech_section("Data layers and their weights",
        p("Weight: how much a layer moves a district's score for each step up the country's ranking of that layer; positive",
          " means more deficiency. Share: the part of the model's differences between districts that comes from the layer.",
          " Range: where the weight fell in 90% of re-estimates. Same direction: in how many countries' own models, and in how",
          " many models with one country left out, the layer pushes the same way."),
        p(sprintf(paste("Satellite, climate and soil data account for %s to %s of the model; climate and soil help most in a",
                        "new country, while groups built from household surveys help within a country but not in a new one. A",
                        "simple equal-weight score from the 20 largest-weight layers ranks a new country at %s, against %s for",
                        "the full model. More livestock goes with more iron and B12 deficiency at district level because the",
                        "livestock-keeping areas are the drier and poorer ones; the weights show where deficiency is likely, not",
                        "what to change."),
                  fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0), fmt_num(Q$sparse20_tr), fmt_num(Q$tr))),
        p("Some layer names were reviewed by a person; the others are direct translations of the layer's code, which is always shown.")),

      tech_section("The surveys",
        reactableOutput(ns("outcomes")),
        p(style = "font-size:0.9em; color:#555;",
          "Vitamin A deficiency is retinol-binding protein, adjusted for inflammation with the BRINDA method and converted to",
          " a retinol value with each survey's published conversion, below 0.70 micromoles per litre. Iron deficiency is each",
          " survey's inflammation-adjusted ferritin below 12 micrograms per litre (children) or 15 (women), not iron-deficiency",
          " anaemia. Children are 6 to 59 months old; women are 15 to 49 and not pregnant. Sierra Leone's child data as supplied",
          " cover only the children found to be anaemic. No variable from the biomarker surveys is used as a predictor."),
        p(style = "font-size:0.9em; color:#555;", GENERAL_CAVEAT)),

      tech_section("Notes on each biomarker",
        tags$ul(lapply(names(biomarker_caveats_long), function(k) tags$li(strong(meta$outcome_labels[[k]] %||% k, ": "), biomarker_caveats_long[[k]])))),

      tech_section("Limits",
        tags$ul(
          tags$li("The model ranks districts. For a published prevalence in a surveyed district a geostatistical model is",
                  " better, and the national prevalence still needs a probability sample."),
          tags$li(sprintf(paste("Prevalence levels do not carry over to a new country; a new country needs at least a national figure.",
                                "A separate model of national prevalence from national indicators (World Bank and similar), scored",
                                "on WHO's national survey database with each country left out in turn, was typically %s to %s points off",
                                "(%s to %s countries per nutrient), and for %s of %s nutrient groups no better than the average of the",
                                "other countries. Predicting vitamin A for our four countries this way missed by up to %s points."),
                          fmt_num(Q$natpred_lo, 0), fmt_num(Q$natpred_hi, 0), Q$natpred_ctry_lo, Q$natpred_ctry_hi,
                          Q$natpred_nobetter, Q$natpred_panels, fmt_num(Q$natpred_own_hi, 0))),
          tags$li("Four training countries, three in West Africa, limit what can be said about new countries."),
          tags$li("Where a survey exists, averaging neighbouring districts does almost as well as the model."),
          tags$li("Rare outcomes and outcomes measured only in Malawi have little to rank or cannot be checked elsewhere."))),

      tech_section("What changed since the earlier version",
        tags$ul(
          tags$li("A within-country accuracy of 0.06 came from one random split; repeated splits gave 0.22, and the corrected",
                  " model reaches ", fmt_num(Q$infill), "."),
          tags$li("A comparison method that the models appeared to lose to had used the predicted district's own respondents."),
          tags$li("A regional result of 0.50 to 0.56 for new countries rested on three of four countries and was withdrawn."),
          tags$li("An estimate of the best achievable score was about five times too low."),
          tags$li("A gain from adding a national survey figure (0.16 to 0.41) did not hold up and was withdrawn."),
          tags$li("Predicting which individuals are deficient has an AUC of about ", fmt_num(Q$il_auc), " and is no longer offered.")))
    )
  )
}

mod_technical_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$perf <- renderReactable({
      B <- EV$benchmarks_summary; validate(need(!is.null(B), "Results summary not built."))
      d <- B[B$target == "level" & B$arm %in% names(arm_label), ] |> select(arm, estimand, mean_spearman) |>
        pivot_wider(names_from = estimand, values_from = mean_spearman)
      d$Method <- arm_label[d$arm]; d <- d[order(-d$infill), ]
      t <- data.frame(Method = d$Method, `District hidden` = round(d$infill, 2), `Region hidden` = round(d$region, 2),
                      `Country left out` = round(d$country, 2), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, rownames = FALSE,
                defaultColDef = colDef(format = colFormat(digits = 2), na = ""))
    })
    output$weights <- renderUI({
      WS <- EV$weight_sources; if (is.null(WS)) return(NULL)
      lab <- c(domain_index = "This model (weights from rank correlations)", ridge_min = "Ridge regression",
               domain_enet = "Elastic net", lasso_min = "Lasso",
               index_decor = "Weights adjusted for overlap between groups", index_soft1 = "Small weights set to zero",
               sparse20 = "Equal-weight score from 20 layers")
      w <- WS[WS$arm %in% names(lab) & WS$target == "level", ] |> select(arm, estimand, mean_spearman) |>
        pivot_wider(names_from = estimand, values_from = mean_spearman)
      w$Method <- lab[w$arm]; w <- w[order(match(w$arm, names(lab))), ]
      tags$table(class = "table table-sm", style = "font-size:0.9em;",
                 tags$thead(tags$tr(tags$th("Ways of weighting the same data groups"), tags$th("District hidden"), tags$th("Region hidden"), tags$th("Country left out"))),
                 tags$tbody(lapply(seq_len(nrow(w)), function(i) tags$tr(tags$td(w$Method[i]), tags$td(fmt_num(w$infill[i])), tags$td(fmt_num(w$region[i])), tags$td(fmt_num(w$country[i]))))))
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
