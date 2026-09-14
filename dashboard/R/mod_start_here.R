# =============================================================================
# Module: Start here
# =============================================================================
# The plain-language entry point: what the dashboard shows, how to read the
# map, one worked example with honest numbers, and the checklist of what the
# models can and cannot do yet. Every number is computed from the bundles.

mod_start_here_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    uiOutput(ns("hero")),

    card(
      card_header("What this dashboard shows"),
      card_body(
        p("A ranking of districts by likely micronutrient deficiency, for every district of four countries,",
          " built from public data anyone can download and scored against the countries' own biomarker surveys.",
          " Most districts have never had a blood test; this says which ones to look at first."),
        tags$ul(
          tags$li(strong("Who it is for: "), "ministries of health, funders and survey planners deciding where",
                  " programmes and the next survey should go."),
          tags$li(strong("What it gives: "), "a priority score per district, how sure the model is about each,",
                  " a planning prevalence anchored to the national survey, and the survey's own figure where one exists."),
          tags$li(strong("What it is not: "), "a measurement. The model ranks; a survey measures. It says nothing",
                  " about causes or about which programme to run.")
        )
      )
    ),

    card(
      card_header("How to read the map"),
      card_body(
        tags$ul(
          tags$li(strong("Colour is the priority score, "), "from 100 (ranked worst in its country) down. Darker means",
                  " the model puts the district nearer the top of the list."),
          tags$li(strong("A dark outline "), "means the district had its own survey clusters; the survey's estimate",
                  " is shown alongside the model's. Thin outlines are districts the survey never reached, which is",
                  " where the model earns its place."),
          tags$li(strong("How sure "), "shows, for surveyed districts, how often the district landed in the worst fifth",
                  " when the model was refitted 40 times with that district hidden. Firm at the top and bottom,",
                  " uncertain in the middle."),
          tags$li(strong("Planning prevalence "), "turns the ranking into a percentage by anchoring it to the",
                  " country's national survey figure. The order is the model's; the level is the survey's."),
          tags$li(strong("Click a district "), "for its numbers and for the predictors that push it up or down the list.")
        ),
        p(em("Differences in ranking accuracy under 0.03 are ties. Read percentages in bands, not to the decimal."))
      )
    ),

    card(
      card_header("Worked example: children's vitamin A in Ghana"),
      card_body(uiOutput(ns("example")),
                methods_note("The model's ranking is scored with each district hidden from the fit, so the",
                             " accuracy quoted is for districts the model had not seen. The burden figures are the",
                             " protocol's own targeting test, averaged over all 24 country and outcome combinations;",
                             " the Ghana-specific line uses the same test for this one cell."))
    ),

    card(
      card_header("What the models can and cannot do yet"),
      card_body(uiOutput(ns("checklist")))
    ),

    card(
      card_header("What to trust, and what not to"),
      card_body(
        tags$ul(
          tags$li(strong("Use the ranking, not the percentage. "), "Which districts are worst off is the part that",
                  " holds up. The planning prevalence borrows its level from the national survey."),
          tags$li(strong("Where a survey exists, smoothing between neighbours does almost as well. "),
                  "The model's value is in the districts, regions and countries no survey reached."),
          tags$li(strong("The survey's own district values are noisy. "),
                  uiOutput(ns("ceiling_line"), inline = TRUE)),
          tags$li(strong("Some nutrients are measured better than others. "), "Vitamin A in women and B12 rest on",
                  " weaker markers, zinc is Malawi only, and rare outcomes have low ceilings. The note under the",
                  " outcome selector says which."),
          tags$li(GENERAL_CAVEAT)
        )
      )
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
          box(fmt_num(Q$infill), "ranking accuracy inside a surveyed country",
              sprintf("against %s for the survey's own regional averages and %s for chance", fmt_num(Q$infill_jk), fmt_num(Q$null_d))),
          box(fmt_num(Q$cs), "in a country the model has never seen",
              sprintf("from climate and soil layers alone (%s with everything); chance is %s", fmt_num(Q$tr), fmt_num(Q$null_d))),
          box(fmt_pct(Q$cap_index, 0), "of a country's deficient people are in the fifth of districts the model ranks worst",
              sprintf("against %s for the survey's regional averages and %s with perfect knowledge", fmt_pct(Q$cap_jk, 0), fmt_pct(Q$cap_oracle, 0))))
    })

    output$ceiling_line <- renderUI({
      span(sprintf(paste("A perfect predictor could reach about %s on the biomarker level and %s on prevalence;",
                         "the model reaches %s and %s, about two-thirds of the way, and is at the ceiling in %s of %s cells."),
                   fmt_num(Q$ceiling_level), fmt_num(Q$ceiling_prev), fmt_num(Q$infill), fmt_num(Q$infill_prev), Q$vc_at, Q$vc_n))
    })

    output$example <- renderUI({
      ck <- "ghana"; oc <- "child_vitA"
      d <- idx_districts[idx_districts$country_key == ck & idx_districts$outcome == oc, ]
      nat <- idx_national[idx_national$country_key == ck & idx_national$outcome == oc, ]
      if (!nrow(d)) return(p(em("Ghana children's vitamin A is not in the data build.")))
      n <- nrow(d); k <- ceiling(n / 5)
      firm <- sum(d$p_worst_fifth >= 0.8, na.rm = TRUE); unlikely <- sum(d$p_worst_fifth <= 0.2, na.rm = TRUE)
      scored <- sum(is.finite(d$p_worst_fifth))
      worst <- d[order(d$rank_worst), ][seq_len(min(5, n)), ]
      # honest per-cell numbers from the protocol's benchmark and targeting tables
      BC <- EV$benchmarks_cells; TCc <- EV$targeting_cells
      rho <- if (!is.null(BC)) g1(BC$spearman[BC$country == "Ghana" & BC$outcome == oc & BC$estimand == "infill" & BC$arm == "domain_index" & BC$target == "prev"]) else NA
      cap <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "domain_index"], na.rm = TRUE)) else NA
      cap_jk <- if (!is.null(TCc)) g1(mean(TCc$capture_top20[TCc$country == "Ghana" & TCc$outcome == oc & TCc$estimand == "infill" & TCc$arm == "region_mean_jk"], na.rm = TRUE)) else NA
      tagList(
        p(sprintf("Ghana's 2017 survey measured children's vitamin A in %d of its %d districts and found a national prevalence of %s. The other %d districts have no measurement.",
                  nat$n_surveyed, n, fmt_pct(nat$national_prev), n - nat$n_surveyed)),
        p("The model, fitted on the surveyed districts and applied to all of them, puts these five at the top of the list: ",
          strong(paste(worst$Admin2, collapse = ", ")), ". ",
          if (is.finite(rho)) sprintf("Scored with each surveyed district hidden in turn, its ranking accuracy for this outcome is %s against %s for chance.", fmt_num(rho), fmt_num(Q$null_d)) else ""),
        p(sprintf(paste("Refitted 40 times with each district hidden, %d of the %d surveyed districts land in the worst fifth at least 80 percent of the time and %d",
                        "at most 20 percent of the time. The rest are the uncertain middle, and a programme should treat them as such."), firm, scored, unlikely)),
        if (is.finite(cap)) p(sprintf(paste("Sent to the fifth of districts the model ranks worst, a programme would reach about %s of Ghana's vitamin A deficient",
                                             "children; the survey's own regional averages reach %s."), fmt_pct(cap, 0), fmt_pct(cap_jk, 0))),
        div(style = "margin-top:10px;",
            actionButton(ns("go_map"), "See this on the map", icon = shiny::icon("map"), class = "btn-sm btn-primary"),
            actionButton(ns("go_targeting"), "See the targeting test", icon = shiny::icon("bullseye"),
                         class = "btn-sm btn-outline-primary", style = "margin-left:6px;"))
      )
    })

    output$checklist <- renderUI({
      row <- function(q, a, ev) tags$tr(tags$td(q), tags$td(strong(a)), tags$td(style = "color:#555;", ev))
      tags$table(class = "table table-sm", style = "font-size:0.95em;",
        tags$thead(tags$tr(tags$th("Can we..."), tags$th("Answer"), tags$th("Evidence"))),
        tags$tbody(
          row("Rank districts inside a surveyed country?", "Yes",
              sprintf("%s against %s for the survey's regional averages; chance %s.", fmt_num(Q$infill), fmt_num(Q$infill_jk), fmt_num(Q$null_d))),
          row("Rank districts in a country with no survey?", "Yes, more roughly",
              sprintf("%s with everything, %s from climate and soil, positive in %s of %s country-outcome pairs.", fmt_num(Q$tr), fmt_num(Q$cs), Q$tr_pos, Q$tr_n)),
          row("Say which public data carry the signal?", "Yes",
              sprintf("Satellite imagery, climate and soil carry %s to %s of the model and are what travels across borders.", fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0))),
          row("Give a prevalence figure for a new country from the model alone?", "Not yet", "Levels do not cross borders; rankings do."),
          row("Give a prevalence figure with one national blood sample added?", "Yes",
              sprintf("District levels within about %s points at 5 percent of a full survey's sample.", fmt_num(Q$ar_a1, 1))),
          row("Replace biomarker surveys?", "No", sprintf("The ceiling is the survey's; the model is at it in %s of %s cells.", Q$vc_at, Q$vc_n)),
          row("Tell the next survey where to sample?", "Yes, in design", "The anchor-and-rank design is a result on these four surveys, not yet a pilot.")
        ))
    })

    if (!is.null(go_to)) {
      observeEvent(input$go_map, go_to("Map explorer", "ghana", "child_vitA"))
      observeEvent(input$go_targeting, go_to("What the ranking buys", NULL, "child_vitA"))
    }
  })
}
