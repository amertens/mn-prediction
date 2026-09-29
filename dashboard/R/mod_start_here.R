# =============================================================================
# Module: Start here (concise version)
# =============================================================================
# One card per question a reader brings: the national prevalence, the
# prevalence in each district, and the ranking of districts. Each card gives the
# answer, the tested number behind it, and a button to the page that shows it.

mod_start_here_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,
    card(card_header("What this is"),
      card_body(
        p(sprintf(paste("This dashboard estimates micronutrient deficiency in all %d districts of The Gambia, Ghana, Sierra Leone",
                        "and Malawi. It uses public data, such as climate, soil, crops and disease maps, so it also covers",
                        "districts that national blood surveys did not reach."), Q$n_districts)),
        p("The model learned from each country's national biomarker survey and was tested on districts, regions and whole",
          "countries it had not seen. How far it can be trusted depends on the question."))),
    uiOutput(ns("answers")),
    p(style = "font-size:0.85em; color:#666; margin-top:0.5em;",
      "Each page explains its figures briefly. ",
      actionLink(ns("go_tech"), "Technical notes"), " has the full detail.")
  )
}

mod_start_here_server <- function(id, go_to = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    answer <- function(title, verdict, num, num_lab, body, btn_id, btn_lab, btn_icon, col = PROXY_COL) {
      card(style = sprintf("border-top:5px solid %s;", col),
        card_header(title),
        card_body(
          p(style = "font-weight:600; margin-bottom:0.4em;", verdict),
          div(style = sprintf("font-size:1.7em; font-weight:700; color:%s; line-height:1.1;", col), num),
          p(style = "margin:4px 0 0.8em; font-size:0.88em; color:#555;", num_lab),
          lapply(body, function(t) p(style = "font-size:0.92em;", t))),
        card_footer(actionButton(ns(btn_id), btn_lab, icon = shiny::icon(btn_icon), class = "btn-sm btn-primary")))
    }

    output$answers <- renderUI({
      xv <- c(Q$xv_prev, Q$xv_level, Q$xv_off); xv <- xv[is.finite(xv)]
      n_xv <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")[Q$xv_countries] %||% Q$xv_countries
      layout_columns(
        col_widths = c(4, 4, 4),
        answer("National prevalence",
          "Needs a survey. The model can help plan it.",
          sprintf("%s to %s points", fmt_num(Q$natpred_lo, 0), fmt_num(Q$natpred_hi, 0)),
          sprintf("typical error when a country's national prevalence was predicted from national indicators, without its own survey (WHO survey data, %s to %s countries per nutrient).",
                  Q$natpred_ctry_lo, Q$natpred_ctry_hi),
          list(
            sprintf("For vitamin A in our four countries the miss reached %s points. The district model learns which districts are worse off, not how common deficiency is overall.",
                    fmt_num(Q$natpred_own_hi, 0)),
            sprintf(paste("In past surveys, keeping only half the districts, chosen at random or spread across the model's ranking,",
                          "still gave the full survey's national figure on average (off by %s points; margin of error ±%s to ±%s).",
                          "Keeping only the highest- and lowest-ranked districts pulled it off by %s points. Not yet tried in a new survey."),
                    fmt_num(mean(c(Q$nat_bias_random, Q$nat_bias_strat)), 1), fmt_num(Q$nat_ci_strat, 1),
                    fmt_num(Q$nat_ci_random, 1), fmt_num(Q$nat_bias_ext, 1))),
          "go_plan", "Plan a survey", "route", col = SURVEY_COL),
        answer("Prevalence in each district",
          "A rough estimate, with a checked range.",
          sprintf("about %s points", fmt_num(Q$prev_err, 0)),
          "typical distance from the survey's own figure, for districts the model did not see.",
          list(
            sprintf("Each district gets a 90%% range. In testing the ranges held the survey's figure %s of the time, but they are wide: a median of ±%s points.",
                    pc(Q$cal_cov), fmt_num(Q$cal_half_med, 0)),
            sprintf(paste("District estimates are scaled to a national figure, so they need one, but a small sample will do: with 5%% of a",
                          "full survey (about %s to %s people per nutrient) they were %s points off, against %s for a district survey of the same size."),
                    fmt_count(Q$anchor_n_lo), fmt_count(Q$anchor_n_hi), fmt_num(Q$ar_a1, 0), fmt_num(Q$ar_b, 0))),
          "go_map", "See the district map", "map"),
        answer("Ranking districts",
          "The model's strongest use.",
          fmt_num(Q$infill),
          sprintf("agreement with the survey's own ranking of districts inside a surveyed country (0 = random, 1 = perfect). The survey's regional averages scored %s.",
                  fmt_num(Q$infill_jk)),
          list(
            sprintf("In a country the model was not trained on it scored %s. Against regional results from %s more countries' surveys, %s to %s.",
                    fmt_num(Q$tr), n_xv, fmt_num(min(xv)), fmt_num(max(xv))),
            "It ranks some nutrients better than others; each map says how reliable its ranking is."),
          "go_trust", "How accurate it is", "gauge"))
    })

    if (!is.null(go_to)) {
      observeEvent(input$go_plan, go_to("Plan a survey"))
      observeEvent(input$go_map, go_to("Map explorer"))
      observeEvent(input$go_trust, go_to("How well it works"))
      observeEvent(input$go_tech, go_to("Technical notes"))
    }
  })
}
