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
      card_header("What to trust — and what not to"),
      card_body(
        tags$ul(
          tags$li(GENERAL_CAVEAT),
          tags$li(strong("Biomarker caveats vary by nutrient "),
                  "(e.g. vitamin A in women, B12) — see the note under the outcome ",
                  "selector and the Methods tab."),
          tags$li(strong("In development: "), "preliminary results for internal ",
                  "review; not for citation or external distribution.")
        )
      )
    )
  )
}

mod_start_here_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    output$ghana_example <- renderUI({
      pd <- if (exists("admin2_area_pred") && !is.null(admin2_area_pred))
              admin2_area_pred else admin2_pred
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
          " of the population, using the model to pick those districts would reach about ",
          strong(pct(targ$reach20)), " of all deficient women — roughly ",
          strong(sprintf("%.1f×", targ$lift)),
          " more than reaching 20% without subnational data."),
        p(em("Open the Map explorer (Ghana, iron in women) and the Decision value ",
             "tab to explore this interactively."),
          style = "color:#555;")
      )
    })
  })
}
