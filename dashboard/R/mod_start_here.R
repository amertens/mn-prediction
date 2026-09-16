# =============================================================================
# Module: Start here
# =============================================================================
# The entry point for a reader with no statistics: what the project is, how the
# model was built and checked (including the four training surveys), the three
# headline numbers, how to use the map, and the checklist of what the model can
# and cannot do yet. The worked example is folded away under the checklist.
# Every number is computed from the bundles.

mod_start_here_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,

    card(
      card_header("What this is"),
      card_body(
        p(sprintf(paste("Vitamin A, iron and other micronutrient deficiencies are measured by national blood surveys, which are",
                        "costly and reach only some districts. This dashboard ranks every one of the %d districts of The Gambia,",
                        "Ghana, Sierra Leone and Malawi by how likely it is to be among the worst affected, using public data alone,",
                        "so that programmes and the next survey know where to look first."), Q$n_districts)),
        p(sprintf(paste("The model was trained on the four countries' national biomarker surveys (The Gambia 2018, Ghana 2017,",
                        "Sierra Leone 2013, Malawi 2015 to 2016), which between them measured deficiency in %d of those districts.",
                        "For every district, surveyed or not, we assembled %d public data layers in %d groups: satellite imagery,",
                        "climate, soil, crops, livestock, malaria and other disease, market prices and summaries from the public",
                        "MICS and household budget surveys; the DHS aggregates are kept as a check, not in the model, so every",
                        "layer is one a country without a recent DHS can rebuild. Each group is condensed to a few summary",
                        "scores, each score is weighted by how well it tracked",
                        "deficiency in the surveyed districts, and the weighted scores are added up. No setting is chosen by hand.",
                        "The result is a ranking of districts, not a measurement."), Q$n_surveyed, Q$n_predictors, Q$n_domains)),
        p("Every result is scored on districts the model had not seen: districts hidden one at a time inside a country, whole",
          " regions hidden, and whole countries hidden. The survey's own regional averages and a shuffled outcome are scored the",
          " same way, so each number has a comparison. The ranking is also applied to Cote d'Ivoire, which has never had a",
          " biomarker survey.")
      )
    ),

    uiOutput(ns("hero")),

    card(
      card_header("How to use the map"),
      card_body(
        tags$ul(
          tags$li("Choose a country and an outcome. Darker districts are ranked nearer the top of the list; the colour is a",
                  " priority score from 100 (ranked worst in the country) down."),
          tags$li("Districts with a dark outline had their own survey clusters, so the survey's estimate is shown beside the",
                  " model's. Thin outlines are districts the survey never reached, which is where the model earns its place."),
          tags$li("Switch the view to see how sure the model is about each surveyed district, or a planning prevalence that turns",
                  " the ranking into a percentage anchored to the national survey."),
          tags$li("Click a district for its numbers and for the public data that push it up or down the list.")
        ),
        p(em("Use the ranking more than the percentage. Which districts are worst off is the part that holds up; the percentage",
             " borrows its level from the national survey and should be read in bands, not to the decimal."))
      )
    ),

    card(
      card_header("What the model can and cannot do yet"),
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

    output$hero <- renderUI({
      box <- function(num, lab, sub) div(class = "col-md-4",
        div(class = "card h-100", style = "border-left:5px solid #0F7B8A;",
            div(class = "card-body", style = "padding:14px 18px;",
                div(style = "font-size:2.0em; font-weight:700; color:#0F7B8A; line-height:1.1;", num),
                p(style = "margin:6px 0 0; font-size:1.0em;", lab),
                p(style = "margin:4px 0 0; color:#666; font-size:0.85em;", sub))))
      div(class = "row g-3", style = "margin-bottom:1em;",
          box(fmt_num(Q$infill), "match between the model's ranking and the survey's, inside a surveyed country",
              sprintf("on a scale where 0 is guessing and 1 is a perfect match; the survey's own regional averages reach %s", fmt_num(Q$infill_jk))),
          box(fmt_num(Q$cs), "the same match in a country the model has never seen",
              sprintf("from climate and soil layers alone (%s with every layer); guessing rarely beats %s", fmt_num(Q$tr), fmt_num(Q$null_d))),
          box(fmt_pct(Q$cap_index, 0), "of a country's deficient people live in the fifth of districts the model ranks worst",
              sprintf("%s if a programme used the survey's regional averages instead, %s if it knew the truth", fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_oracle, 0))))
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
      # per-cell numbers from the protocol's benchmark and targeting tables
      BC <- EV$benchmarks_cells; TCc <- EV$targeting_cells
      rho <- if (!is.null(BC)) g1(BC$spearman[BC$country == "Ghana" & BC$outcome == oc & BC$estimand == "infill" & BC$arm == "domain_index" & BC$target == "prev"]) else NA
      cap <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "domain_index"], na.rm = TRUE)) else NA
      cap_jk <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "region_mean_jk"], na.rm = TRUE)) else NA
      tagList(
        p(sprintf("Ghana's 2017 survey measured children's vitamin A in %d of its %d districts and found a national prevalence of %s. The other %d districts have no measurement.",
                  nat$n_surveyed, n, fmt_pct(nat$national_prev), n - nat$n_surveyed)),
        p("The model, fitted on the surveyed districts and applied to all of them, puts these five at the top of the list: ",
          strong(paste0(paste(worst$Admin2, collapse = ", "), ".")),
          if (is.finite(rho)) sprintf(" Scored with each surveyed district hidden in turn, it matches the survey's ranking at %s.", fmt_num(rho)) else ""),
        p(sprintf(paste("Refitted 40 times with each district hidden, %d of the %d surveyed districts land in the worst fifth at least 80 percent",
                        "of the time and %d at most 20 percent of the time. The rest are the uncertain middle."), firm, scored, unlikely)),
        if (is.finite(cap)) p(sprintf(paste("A programme sent to the fifth of districts the model ranks worst would reach about %s of Ghana's vitamin A",
                                             "deficient children; one guided by the survey's regional averages would reach %s."),
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
              sprintf("Matches the survey's ranking at %s (1 = perfect); the survey's own regional averages reach %s.", fmt_num(Q$infill), fmt_num(Q$infill_jk))),
          row("Rank districts in a country with no survey?", "Yes, more roughly",
              sprintf("%s with every layer, %s from climate and soil alone; positive in %s of %s country and outcome tests.", fmt_num(Q$tr), fmt_num(Q$cs), Q$tr_pos, Q$tr_n)),
          row("Say which public data carry the signal?", "Yes",
              sprintf("Satellite imagery, climate and soil make up %s to %s of the model and are what works across borders.", fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0))),
          row("Give a prevalence figure for a new country from the model alone?", "Not yet",
              "The model says which districts are worse, not how bad. A percentage needs a national survey figure to anchor it."),
          row("Give a prevalence figure with one small national blood sample added?", "Yes",
              sprintf("A sample 5 percent the size of a full survey puts district figures within about %s percentage points.", fmt_num(Q$ar_a1, 1))),
          row("Replace biomarker surveys?", "No",
              sprintf("The survey's own noise sets a ceiling the model cannot pass; it is at that ceiling in %s of %s cases.", Q$vc_at, Q$vc_n)),
          row("Tell the next survey where to sample?", "Yes, in design",
              "The anchor-and-rank design is a result on these four surveys, not yet a pilot.")
        ))
    })

    if (!is.null(go_to)) {
      observeEvent(input$go_map, go_to("Map explorer", "ghana", "child_vitA"))
      observeEvent(input$go_targeting, go_to("What the ranking buys", NULL, "child_vitA"))
    }
  })
}
