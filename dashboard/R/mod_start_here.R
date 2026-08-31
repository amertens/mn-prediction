# =============================================================================
# Module: Start here  (plain-language guide + worked use case)
# =============================================================================
# A non-technical entry point: what the dashboard shows, how to read the map,
# and a concrete worked example (Ghana iron) showing how subnational estimates
# would be used. Numbers in the example are computed live from the data.

mod_start_here_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = 12,

    # The single number worth leading with. It comes from the Decision value
    # tab's targeting analysis, which was five tabs in and doing the persuading
    # where nobody arriving would see it.
    uiOutput(ns("hero")),

    card(
      card_header("What this dashboard shows"),
      card_body(
        p("It maps the estimated ", strong("prevalence of micronutrient deficiency"),
          " (e.g. iron, vitamin A) at the ", strong("district level"),
          " across four countries — including districts that were never surveyed —",
          " so that limited programs can be aimed where the need is greatest."),
        tags$ul(
          tags$li(strong("Who it's for: "), "ministries of health, funders, and ",
                  "researchers planning surveys or interventions."),
          tags$li(strong("What's new: "), "national surveys give one number per ",
                  "country; this fills in the ", em("within-country"), " picture."),
          tags$li(strong("What it is not: "), "it predicts where deficiency is high; ",
                  "it does not recommend specific programs or doses.")
        )
      )
    ),

    card(
      card_header("How to read the map (the “Map explorer” tab)"),
      card_body(
        tags$ul(
          tags$li(strong("Colour = predicted prevalence. "), "Darker = higher ",
                  "deficiency. The legend shows the scale."),
          tags$li(strong("Border = data source. "), "A dark, thick border means the ",
                  "district had an actual survey; a thin grey border means the value ",
                  "is model-predicted (no local survey)."),
          tags$li(strong("Prediction model "), "(left panel): the ",
                  em("area-level"), " models cover every district; the ",
                  em("Fay-Herriot"), " and ", em("BYM2"), " options add a 95% ",
                  "uncertainty range; the ", em("SuperLearner"),
                  " option is a sensitivity analysis (surveyed districts only)."),
          tags$li(strong("Display layer: "), "switch between predicted prevalence, ",
                  "the survey-observed value, the gap between them, the WHO severity ",
                  "class, and ", strong("WHO-class certainty"),
                  " (whether the uncertainty is wide enough that the district's ",
                  "severity category could change)."),
          tags$li(strong("Click a district "), "for its detailed numbers; the ",
                  "“District profiles” tab profiles one place (or the whole country) ",
                  "across all outcomes.")
        ),
        p(em("Tip: the percentages are planning estimates — read them in broad ",
             "bands (e.g. above/below the 20% “severe” line), not to the decimal."))
      )
    ),

    card(
      card_header("Worked example: is subnational targeting worth it? (Ghana, iron in women)"),
      card_body(
        uiOutput(ns("ghana_example")),
        methods_note(
          "This illustrates the model's intended use — identifying whether ",
          "deficiency is concentrated in hotspots (so targeting helps) versus ",
          "spread evenly (so a national approach is better). It is not advice on ",
          "which program to run."
        )
      )
    ),

    card(
      card_header("What to trust, and what not to"),
      card_body(
        tags$ul(
          tags$li(GENERAL_CAVEAT),
          tags$li(strong("Some nutrients are measured better than others. "),
                  "Vitamin A in women and B12 rest on weaker biomarkers. The ",
                  "note under the outcome selector says which, and the Methods ",
                  "tab explains why."),
          # Was a repeat of the old "not for citation" banner. That framing is
          # gone from the header, and leaving it here would have been the one
          # place the app still told people not to use it.
          tags$li(strong("Use the ranking more than the number. "),
                  "These are working estimates. Which districts are worst off ",
                  "is the part that holds up; the exact percentage is a planning ",
                  "figure, and how firm it is varies a lot by country and ",
                  "nutrient. Decision value shows where it holds and where it ",
                  "does not."),
          # Added 2026-08-30. The district maps are anchored to design-based
          # regional totals, and that anchoring is what makes them usable --
          # it more than doubles rank correlation (0.164 -> 0.413) and cuts
          # mean absolute bias from 3.2 to 1.6 pp, better in 20 of 24 cells.
          # Saying so is more honest than presenting the map as if the
          # covariates alone produced it.
          tags$li(strong("The survey sets the level; the model sets the pattern. "),
                  "District estimates are anchored to each region's ",
                  "design-based survey total. That anchoring is doing much of ",
                  "the work: without it the same model's district ranking is ",
                  "less than half as accurate. Predictions for a country with ",
                  "no survey of its own inherit the pattern but not the level, ",
                  "and should be read as a ranking only. Scored arm by arm ",
                  "under ", tags$em("Resolution & anchoring"), "."),
          tags$li(strong("Geography explains most of the district pattern. "),
                  "A model using only where a district is -- no proxy ",
                  "variables at all -- matches the full 294-predictor model at ",
                  "district level. Tested individually against district ",
                  "prevalence, no single proxy survives multiple-comparison ",
                  "correction once spatial structure is accounted for. The ",
                  "proxies earn their place at region level, not district ",
                  "level.")
        )
      )
    )
  )
}

#' @param go_to optional callback, `function(tab, country, outcome)`, supplied by
#'   the app server. When present the worked example gets buttons that switch
#'   tabs with the country and outcome already set, instead of telling the
#'   reader to go and find them.
mod_start_here_server <- function(id, go_to = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Share of a country's deficient population you reach by treating the
    # highest-prevalence districts covering `frac` of the population, versus
    # reaching the same share without knowing where deficiency is concentrated.
    targeting_gain <- function(ck, oc, frac = 0.20) {
      tryCatch({
        pd <- pred_model_data(DEFAULT_PRED_MODEL)
        df <- as.data.frame(get_country_admin2(ck, oc, admin2_bnds, pd, admin2_pop))
        df <- df[is.finite(df$pred_prev) & is.finite(df$population) & df$population > 0, ]
        if (!nrow(df)) return(NULL)
        w <- df$population; burden <- df$pred_prev * w
        ord <- order(df$pred_prev, decreasing = TRUE)
        cx <- cumsum(w[ord]) / sum(w); cy <- cumsum(burden[ord]) / sum(burden)
        reach <- stats::approx(c(0, cx), c(0, cy), xout = frac, rule = 2)$y
        list(reach = reach, lift = reach / frac, frac = frac)
      }, error = function(e) NULL)
    }

    output$hero <- renderUI({
      g <- targeting_gain("ghana", "women_iron")
      if (is.null(g) || !is.finite(g$lift)) return(NULL)
      div(
        class = "card",
        style = "border-left:5px solid #2c7bb6; margin-bottom:1em;",
        div(
          class = "card-body",
          style = "padding:14px 18px;",
          div(style = "font-size:2.1em; font-weight:700; color:#2c7bb6; line-height:1.1;",
              sprintf("%.1f×", g$lift)),
          p(style = "margin:6px 0 0; font-size:1.02em;",
            "more deficient women reached, for the same program size."),
          p(style = "margin:4px 0 0; color:#555; font-size:0.9em;",
            sprintf(paste("Reaching the worst-affected %.0f%% of Ghana's population",
                          "using these district estimates covers about %.0f%% of all",
                          "women with iron deficiency. Without district data you",
                          "would reach roughly %.0f%%."),
                    100 * g$frac, 100 * g$reach, 100 * g$frac)),
          p(style = "margin:8px 0 0; color:#777; font-size:0.82em;",
            "Ghana, iron deficiency in women. The gain varies by country and ",
            "nutrient, and is much smaller for some. Decision value has the rest.")
        )
      )
    })

    # Wire the worked example's buttons, when the app supplied a handler.
    if (!is.null(go_to)) {
      observeEvent(input$go_map, go_to("Map explorer", "ghana", "women_iron"))
      observeEvent(input$go_decision, go_to("Decision value", "ghana", "women_iron"))
    }

    output$ghana_example <- renderUI({
      pd <- pred_model_data(DEFAULT_PRED_MODEL)
      ck <- "ghana"; oc <- "women_iron"

      natl <- tryCatch({
        d <- get_country_admin2(ck, oc, admin2_bnds, pd, admin2_pop)
        national_aggregate(d)$pred_prev_natl
      }, error = function(e) NA_real_)

      reg <- tryCatch({
        a1 <- as.data.frame(get_country_admin1(ck, oc, admin1_bnds, admin2_bnds, pd, admin2_pop))
        a1 <- a1[is.finite(a1$pred_prev), , drop = FALSE]
        list(hi = a1[which.max(a1$pred_prev), ], lo = a1[which.min(a1$pred_prev), ])
      }, error = function(e) NULL)

      targ <- tryCatch({
        df <- as.data.frame(get_country_admin2(ck, oc, admin2_bnds, pd, admin2_pop))
        df <- df[is.finite(df$pred_prev) & is.finite(df$population) & df$population > 0, ]
        w <- df$population; burden <- df$pred_prev * w
        ord <- order(df$pred_prev, decreasing = TRUE)
        cx <- cumsum(w[ord]) / sum(w); cy <- cumsum(burden[ord]) / sum(burden)
        reach20 <- stats::approx(c(0, cx), c(0, cy), xout = 0.20, rule = 2)$y
        list(reach20 = reach20, lift = reach20 / 0.20)
      }, error = function(e) NULL)

      pct <- function(x) if (is.na(x)) "—" else sprintf("%.0f%%", x * 100)
      tagList(
        p("A single national figure for Ghana — about ",
          strong(pct(natl)), " of women with iron deficiency — hides large ",
          "differences between districts."),
        if (!is.null(reg)) p(
          "The model estimates roughly ", strong(pct(reg$hi$pred_prev)),
          " in ", strong(reg$hi$Admin1), " versus ", strong(pct(reg$lo$pred_prev)),
          " in ", strong(reg$lo$Admin1), " — the kind of north–south gap a ",
          "national average would mask."),
        if (!is.null(targ)) p(
          "If a program could only reach the worst ", strong("20%"),
          " of the population, using these estimates to pick those districts would reach about ",
          strong(pct(targ$reach20)), " of all deficient women — roughly ",
          strong(sprintf("%.1f×", targ$lift)),
          " more than reaching 20% without district-level data."),
        div(
          style = "margin-top:10px;",
          # shiny::icon(), not bsicons::bs_icon() — actionButton validates its
          # icon argument and rejects a raw SVG.
          actionButton(ns("go_map"), "See this on the map",
                       icon = shiny::icon("map"),
                       class = "btn-sm btn-primary"),
          actionButton(ns("go_decision"), "See the full targeting analysis",
                       icon = shiny::icon("bullseye"),
                       class = "btn-sm btn-outline-primary",
                       style = "margin-left:6px;")
        )
      )
    })
  })
}
