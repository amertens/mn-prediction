# =============================================================================
# Module: Start here
# =============================================================================
# The entry point for a reader with no statistics: what the dashboard is, how
# the model was built and checked, the headline numbers, a short path through
# the site, how to use the map, and what the model can and cannot do. The
# worked example is folded away. Every number is computed from the bundles.

mod_start_here_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,

    card(
      card_header("What this is"),
      card_body(
        p(sprintf(paste("Micronutrient deficiencies are measured by national blood surveys, which are expensive and cover",
                        "only some districts. This dashboard ranks all %d districts of The Gambia, Ghana, Sierra Leone and",
                        "Malawi by how likely they are to be among the worst affected. It uses public data, so it also covers",
                        "districts the surveys did not reach."), Q$n_districts)),
        p(sprintf(paste("The model learned from each country's national biomarker survey (The Gambia 2018, Ghana 2017,",
                        "Sierra Leone 2013, Malawi 2015 to 2016), which together collected blood samples in %d of the districts.",
                        "For every district we gathered %d public data layers, such as climate, soil, satellite imagery,",
                        "crops, livestock, disease burden, food prices and public household surveys; the model uses %s of them.",
                        "It combines the layers into one score per district and ranks districts by that score. The result",
                        "is a ranking, not a measurement."),
                  Q$n_surveyed, Q$n_predictors, if (is.finite(Q$n_in_model)) Q$n_in_model else "383")),
        p("Every accuracy figure on this site comes from districts, regions or whole countries that were hidden from the",
          " model while it was built. Simpler alternatives, such as the survey's own regional averages and a random",
          " ranking, are tested the same way, so each figure has a comparison. The ranking is also applied to",
          " Cote d'Ivoire, which has not had a national biomarker survey.")
      )
    ),

    uiOutput(ns("hero")),
    uiOutput(ns("hero2")),

    card(
      card_header("Short on time?"),
      card_body(
        div(style = "display:flex; gap:10px; flex-wrap:wrap;",
          actionButton(ns("go5_map"), "1. Your country's map", icon = shiny::icon("map"), class = "btn-sm btn-primary"),
          actionButton(ns("go5_trust"), "2. How accurate it is", icon = shiny::icon("gauge"), class = "btn-sm btn-outline-primary"),
          actionButton(ns("go5_plan"), "3. Planning the next survey", icon = shiny::icon("route"), class = "btn-sm btn-outline-primary")),
        p(style = "font-size:0.85em; color:#666; margin:8px 0 0;",
          "These three pages show where deficiency is likely, how far to trust the estimates, and what they mean",
          " for the next survey. The other pages give the supporting detail."))),

    card(
      card_header("How to use the map"),
      card_body(
        tags$ul(
          tags$li("Choose a country and an outcome. Darker districts rank higher on the list (priority score 100 = ranked",
                  " worst in the country)."),
          tags$li("A dark outline marks districts where the survey collected data; their survey estimate is shown next to",
                  " the model's. Districts with a thin outline had no survey data, so the model is the only estimate there."),
          tags$li("Other views show how firmly each district is ranked, and an estimated prevalence based on the national",
                  " survey figure, with a checked range."),
          tags$li("Click a district to see its figures and the data layers that move it up or down the list.")
        ),
        p(em("Rely on the ranking more than on the percentages. Which districts are worse off is the more reliable result;",
             " the percentages take their overall level from the national survey and are best read as rough bands."))
      )
    ),

    card(
      card_header("What the model can and cannot do"),
      card_body(uiOutput(ns("checklist")),
                p(style = "font-size:0.9em; color:#555; margin-top:6px;", GENERAL_CAVEAT))
    ),

    accordion(
      id = ns("more"), open = FALSE,
      accordion_panel("Worked example: children's vitamin A in Ghana", icon = bsicons::bs_icon("geo-alt"),
                      uiOutput(ns("example")))
    )
  )
}

mod_start_here_server <- function(id, go_to = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    box <- function(num, lab, sub, colour) div(class = "col-md-4",
      div(class = "card h-100", style = sprintf("border-left:5px solid %s;", colour),
          div(class = "card-body", style = "padding:14px 18px;",
              div(style = sprintf("font-size:2.0em; font-weight:700; color:%s; line-height:1.1;", colour), num),
              p(style = "margin:6px 0 0; font-size:1.0em;", lab),
              p(style = "margin:4px 0 0; color:#666; font-size:0.85em;", sub))))

    output$hero <- renderUI({
      div(class = "row g-3", style = "margin-bottom:1em;",
          box(fmt_num(Q$infill), "how closely the model's ranking matches the survey's, within a surveyed country",
              sprintf("On a scale where 0 is a random order and 1 is a perfect match. The survey's own regional averages reach %s.", fmt_num(Q$infill_jk)),
              PROXY_COL),
          box(fmt_num(Q$tr), "the same measure in a country the model was not trained on",
              sprintf(paste("Using all data layers. With only climate and soil layers it reaches %s; that combination was chosen",
                            "after seeing these results and is now being tested on new surveys. A random ranking stays below %s in 95%% of tries."),
                      fmt_num(Q$cs), fmt_num(Q$null_d)),
              PROXY_COL),
          box(fmt_pct(Q$cap_index, 0), "of a country's deficient people live in the fifth of districts the model ranks worst",
              sprintf("%s for the fifth chosen with the survey's regional averages, %s for a random fifth, and %s for the true worst fifth.",
                      fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_null, 0), fmt_pct(Q$cap_oracle, 0)),
              PROXY_COL))
    })

    output$hero2 <- renderUI({
      div(class = "row g-3", style = "margin-bottom:1em;",
          box(sprintf("%s countries", if (is.finite(Q$xv_countries)) Q$xv_countries else 6),
              "outside the training data where the ranking was checked against survey results held by the WHO",
              sprintf(paste("The ranking matched the surveys' regional results at %s (biomarker levels) and %s (prevalence) in four",
                            "African countries, and %s (prevalence) in Pakistan and India."),
                      f2(Q$xv_level), f2(Q$xv_prev), f2(Q$xv_off)),
              SURVEY_COL),
          box(f2(Q$strong_mean), sprintf("average ranking accuracy in the model's best %s of %s within-country tests", Q$strong_n, Q$infill_cells),
              sprintf(paste("These are vitamin A and iron in The Gambia, B12 in Ghana and Malawi, and children's iron in Ghana.",
                            "The other %s tests average %s."), Q$weak_n, f2(Q$weak_mean)),
              SURVEY_COL),
          box(pc(Q$cal_cov), "of the checked 90% prevalence ranges contained the survey's own district figure",
              sprintf(paste("Tested on districts the model had not seen. The ranges are wide (a median of %s percentage points",
                            "either side), partly because each district's survey figure is itself uncertain."),
                      fmt_num(Q$cal_half_med, 0)),
              SURVEY_COL))
    })

    output$example <- renderUI({
      ck <- "ghana"; oc <- "child_vitA"
      d <- idx_districts[idx_districts$country_key == ck & idx_districts$outcome == oc, ]
      nat <- idx_national[idx_national$country_key == ck & idx_national$outcome == oc, ]
      if (!nrow(d)) return(p(em("Ghana children's vitamin A is not in the data build.")))
      n <- nrow(d)
      firm <- sum(d$p_worst_fifth >= 0.8, na.rm = TRUE); unlikely <- sum(d$p_worst_fifth <= 0.2, na.rm = TRUE)
      scored <- sum(is.finite(d$p_worst_fifth))
      worst <- d[order(d$rank_worst), ][seq_len(min(5, n)), ]
      BC <- EV$benchmarks_cells; TCc <- EV$targeting_cells
      rho <- if (!is.null(BC)) g1(BC$spearman[BC$country == "Ghana" & BC$outcome == oc & BC$estimand == "infill" & BC$arm == "domain_index" & BC$target == "prev"]) else NA
      cap <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "domain_index"], na.rm = TRUE)) else NA
      cap_jk <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "region_mean_jk"], na.rm = TRUE)) else NA
      tagList(
        p(sprintf("Ghana's 2017 survey measured children's vitamin A in %d of its %d districts and found a national prevalence of %s. The other %d districts have no measurement.",
                  nat$n_surveyed, n, fmt_pct(nat$national_prev), n - nat$n_surveyed)),
        p("The model, built on the surveyed districts and applied to all of them, ranks these five highest: ",
          strong(paste0(paste(worst$Admin2, collapse = ", "), ".")),
          if (is.finite(rho)) sprintf(" With each surveyed district hidden in turn, its ranking matched the survey's at %s.", fmt_num(rho)) else ""),
        p(sprintf(paste("In 40 test runs with the district hidden, %d of the %d surveyed districts were placed in the worst fifth in at",
                        "least 80%% of runs, and %d in at most 20%%. The rest are uncertain."), firm, scored, unlikely)),
        if (is.finite(cap)) p(sprintf(paste("A programme directed to the fifth of districts the model ranks worst would reach about %s of Ghana's",
                                             "children with vitamin A deficiency; one guided by the survey's regional averages would reach %s."),
                                       fmt_pct(cap, 1), fmt_pct(cap_jk, 1))),
        div(style = "margin-top:10px;",
            actionButton(ns("go_map"), "See this on the map", icon = shiny::icon("map"), class = "btn-sm btn-primary"),
            actionButton(ns("go_targeting"), "See the targeting test", icon = shiny::icon("bullseye"),
                         class = "btn-sm btn-outline-primary", style = "margin-left:6px;"))
      )
    })

    output$checklist <- renderUI({
      row <- function(q, a, ev) tags$tr(tags$td(q), tags$td(strong(a)), tags$td(style = "color:#555;", ev))
      tags$table(class = "table table-sm", style = "font-size:0.95em;",
        tags$thead(tags$tr(tags$th("Can it..."), tags$th("Answer"), tags$th("Evidence"))),
        tags$tbody(
          row("Rank districts inside a surveyed country?", "Yes",
              sprintf("Ranking accuracy %s, compared with %s for the survey's own regional averages.", fmt_num(Q$infill), fmt_num(Q$infill_jk))),
          row("Rank districts in a country with no survey?", "Yes, less accurately",
              sprintf("%s using all data layers (%s with climate and soil only); above zero in %s of %s country-outcome tests.",
                      fmt_num(Q$tr), fmt_num(Q$cs), Q$tr_pos, Q$tr_n)),
          row("Work in countries beyond these four?", "Yes, in a first check",
              sprintf("Compared with survey results held by the WHO for %s more countries, the regional rankings matched at %s to %s.",
                      if (is.finite(Q$xv_countries)) Q$xv_countries else 6, fmt_num(Q$xv_prev), fmt_num(Q$xv_level))),
          row("Show which data matter?", "Yes",
              sprintf("Climate and soil carry most of what transfers to a new country. Satellite, climate and soil data together make up %s to %s of the model.",
                      fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0))),
          row("Give a prevalence for a country with no survey at all?", "No",
              "The model ranks districts. A percentage needs at least a national survey figure."),
          row("Give district percentages if a small national survey is added?", "Yes, roughly",
              sprintf("With a national sample 5%% the size of a full survey, district estimates were off by a median of %s percentage points.", fmt_num(Q$ar_a1, 1))),
          row("Replace biomarker surveys?", "No",
              "The national prevalence still needs a probability survey, and the model can only be as good as the survey data it learns from."),
          row("Help plan the next survey?", "Yes, as a design aid",
              "Plan a survey shows tests on these four surveys. There has been no field trial yet."),
          row("What would improve it?", "More and better surveys",
              sprintf(paste("Accuracy in countries left out of training rose with each survey added (%s with one, %s with three), and",
                            "more survey clusters per district would make the survey figures the model learns from less noisy."),
                      fmt_num(Q$lc1), fmt_num(Q$lc_max)))
        ))
    })

    if (!is.null(go_to)) {
      observeEvent(input$go_map, go_to("Map explorer", "ghana", "child_vitA"))
      observeEvent(input$go_targeting, go_to("What the ranking buys", NULL, "child_vitA"))
      observeEvent(input$go5_map, go_to("Map explorer"))
      observeEvent(input$go5_trust, go_to("How well it works"))
      observeEvent(input$go5_plan, go_to("Plan a survey"))
    }
  })
}
