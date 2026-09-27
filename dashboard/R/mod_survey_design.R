# =============================================================================
# Module: Plan a survey
# =============================================================================
# Three questions a survey planner asks, in order: WHERE should the next
# survey's district visits go (an interactive prioritiser borrowing the logic
# of adaptive geostatistical design), HOW SMALL can the survey be for a given
# job (the anchor-and-rank design, AR-01), and DOES CHOOSING WELL MATTER
# (SP-01, a pre-registered retrospective test on the four surveys).
#
# Honesty rules for this tab: the prioritiser is a screening aid, not a
# sampling frame; stability ranges are stability, not truth (VZ-01); the
# augmentation negatives are shown, not hidden.

PLANNER_PRESETS <- list(
  confirm = c(burden = 45, uncert = 10, thresh = 10, popn = 20, nosvy = 15),
  shrink  = c(burden = 10, uncert = 45, thresh = 25, popn = 20, nosvy = 0),
  cover   = c(burden = 15, uncert = 15, thresh = 0,  popn = 25, nosvy = 45)
)

mod_survey_design_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    id = ns("tabs"),
    # ── 1. Where to sample next ─────────────────────────────────────────────
    nav_panel(
      title = "Where to sample next", icon = bsicons::bs_icon("geo-alt"),
      layout_sidebar(
        sidebar = sidebar(width = 330, title = "What should the next survey learn?",
          selectInput(ns("country"), "Country", choices = country_choices, selected = "ghana"),
          selectInput(ns("outcome"), "Outcome", choices = NULL),
          sliderInput(ns("k"), "District visits the budget allows", min = 5, max = 40, value = 15, step = 1),
          radioButtons(ns("preset"), "Objective",
            choices = c("Confirm the worst-off districts" = "confirm",
                        "Firm up what the model cannot place" = "shrink",
                        "Reach the never-surveyed" = "cover",
                        "Custom weights" = "custom"), selected = "shrink"),
          conditionalPanel(sprintf("input['%s'] == 'custom'", ns("preset")),
            sliderInput(ns("w_burden"), "Predicted burden", 0, 100, 20, 5),
            sliderInput(ns("w_uncert"), "Ranking instability", 0, 100, 30, 5),
            sliderInput(ns("w_thresh"), "Near a WHO threshold", 0, 100, 20, 5),
            sliderInput(ns("w_popn"),   "Population", 0, 100, 20, 5),
            sliderInput(ns("w_nosvy"),  "Never surveyed", 0, 100, 10, 5)),
          hr(),
          uiOutput(ns("plan_summary"))),
        layout_columns(col_widths = c(7, 5),
          card(full_screen = TRUE,
            card_header("The districts this objective sends the survey to"),
            card_body(padding = 0, leafletOutput(ns("plan_map"), height = "560px")),
            card_footer(tags$small(class = "text-muted",
              "Outlined districts are the proposed visits. Fill: the priority each district gets under the current",
              " objective. A visit here means enough clusters in that district to give it its own estimate."))),
          card(card_header("Why each district is on the list"),
            card_body(reactableOutput(ns("plan_table")),
              methods_note("Each district is scored on five ingredients, all shown on hover: the model's predicted burden,",
                " how far its rank moves across refits (instability; a stability measure, not a confidence interval),",
                " how close its planning prevalence sits to the WHO action threshold (where one exists), its population,",
                " and whether the last survey reached it. The objective sets the weights; the list is the top of the",
                " weighted score. This is a screening aid for choosing DISTRICTS - cluster selection inside districts,",
                " sample sizes and weighting still need a survey statistician."))))
        ),
      uiOutput(ns("borrowed"))
    ),
    # ── 2. How small can it be ──────────────────────────────────────────────
    nav_panel(
      title = "How small can a survey be", icon = bsicons::bs_icon("cash-coin"),
      layout_sidebar(
        sidebar = sidebar(width = 330, title = "Design choices",
          selectInput(ns("rank_from"), "Ranking learned from", choices = NULL),
          selectInput(ns("set"), "Predictors in the ranking", choices = NULL),
          selectInput(ns("size_country"), "Read the sample as", choices = country_choices, selected = "ghana"),
          sliderInput(ns("fraction"), "Share of a full survey's sample", min = 0.05, max = 1, value = 0.05, step = 0.05),
          hr(),
          uiOutput(ns("at_fraction"))),
        layout_columns(col_widths = c(12, 12),
          card(card_header("A national anchor plus the ranking beats a small direct survey on district error"),
            card_body(plotlyOutput(ns("mae"), height = "360px"),
              methods_note(sprintf(paste("Median district error in percentage points over 22 country-outcome combinations.",
                "The anchored design measures only a national prevalence from the given share of a full survey's",
                "respondents and orders districts by the transported ranking; the district and regional surveys spend",
                "the same sample on direct estimates. At 5 percent of the sample the anchored design gives %s points,",
                "a regional survey %s and a district survey %s; the district survey needs about %s of the full sample",
                "to match."), fmt_num(Q$ar_a1, 1), fmt_num(Q$ar_c, 1), fmt_num(Q$ar_b, 1), fmt_pct(Q$ar_b_match, 0))))),
          card(card_header("What the cheap design gives up, and what does not work"),
            card_body(plotlyOutput(ns("capture"), height = "300px"),
              methods_note("Burden captured is where the anchored design loses: a district survey of any size reaches more",
                " of the deficient population, because burden sits in populous districts that surveys measure precisely.",
                " Two more honest negatives from the same tests: blending the model into a small survey's own district",
                " estimates does not beat those estimates where clusters exist, and the model does not substitute for",
                " sample size. The anchored design buys a level for the ranking; it does not buy precision."))))
      )
    ),
    # ── 3. Does choosing well matter ────────────────────────────────────────
    nav_panel(
      title = "Does choosing districts well matter?", icon = bsicons::bs_icon("diagram-3"),
      layout_columns(col_widths = c(8, 4),
        card(card_header("Surveying k districts: what the resulting map reaches, by how the k were chosen"),
          card_body(plotlyOutput(ns("planner_curve"), height = "360px"),
            methods_note(sprintf(paste("Pre-registered test (SP-01) on the four surveys: pretend only a share of the surveyed",
              "districts get biomarkers, choose them by each rule, fit the index on those, and score the map on the",
              "districts held back; 40 draws, mean over 22 country-outcome cells. Spreading the visits across the",
              "TRANSPORTED model's gradient beats choosing at random in %s of %s cells at half the districts (mean",
              "gain %s), and choosing districts by population is WORSE than random. The gain from choosing well is",
              "real but modest - the big lever is how many districts, not only which."),
              g1(Q$plan_better), g1(Q$plan_cells), f2(Q$plan_delta, 3))))),
        card(card_header("How to read it"),
          card_body(
            p(style = "font-size:0.9em;", "Half the districts, chosen sensibly, keep about ",
              strong(sprintf("%s of the full survey's ranking accuracy", if (is.finite(Q$plan_spread) && is.finite(Q$infill)) fmt_pct(Q$plan_spread / Q$infill, 0) else "—")),
              sprintf(" (%s against %s).", f2(Q$plan_spread), f2(Q$infill))),
            p(style = "font-size:0.9em;", "The dashed line is what the transported ranking gives with ",
              strong("no in-country biomarkers at all"), " - the floor any survey improves on."),
            p(style = "font-size:0.85em; color:#555;", "Design and readings were written before the numbers:",
              " docs/findings/SP-01_SURVEY_PLANNER_DESIGN_2026-09-27.md.")))),
      layout_columns(col_widths = c(8, 4),
        card(card_header("And the national prevalence - the survey's primary product?"),
          card_body(plotlyOutput(ns("national"), height = "300px"),
            methods_note("The district-selection component of the national estimate: each line is the 95% spread of the",
              " population-weighted national figure across draws of which districts are visited (within-district sampling",
              " noise comes on top of this). Choosing districts at random, or spreading them across the model's gradient",
              " and weighting by stratum, keeps the national number unbiased; targeting the model's extremes is an",
              " informative sample and shifts it. That is why this dashboard proposes model-guided choice as a stratified",
              " probability design, never as 'survey only where the model says it is bad'."))),
        card(card_header("The bias rule"),
          card_body(uiOutput(ns("national_bias")))))
    )
  )
}

mod_survey_design_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # ── 1. the prioritiser ────────────────────────────────────────────────
    observeEvent(input$country, {
      ch <- outcomes_for(input$country)
      sel <- if ((input$outcome %||% "") %in% ch) input$outcome else ch[[1]]
      updateSelectInput(session, "outcome", choices = ch, selected = sel)
      n <- g1(idx_national$n_districts[idx_national$country_key == input$country])
      if (is.finite(n)) updateSliderInput(session, "k", max = max(10, min(60, n - 2)), value = min(15, max(5, round(n / 5))))
    })

    weights_now <- reactive({
      if ((input$preset %||% "shrink") != "custom") PLANNER_PRESETS[[input$preset]]
      else c(burden = input$w_burden, uncert = input$w_uncert, thresh = input$w_thresh,
             popn = input$w_popn, nosvy = input$w_nosvy)
    })

    plan_data <- reactive({
      req(input$country, input$outcome)
      df <- get_country_admin2(input$country, input$outcome)
      req(df, nrow(df) > 0)
      w <- weights_now(); w <- w / max(sum(w), 1)
      rk01 <- function(x) { r <- rank(x, ties.method = "average", na.last = "keep"); (r - 1) / max(sum(is.finite(x)) - 1, 1) }
      comp <- data.frame(
        burden = ifelse(is.finite(df$priority), df$priority / 100, 0),
        uncert = ifelse(is.finite(df$width_share), pmin(df$width_share, 1), 0),
        thresh = { th <- g1(df$th_moderate_plus[is.finite(df$th_moderate_plus)])
                   if (is.finite(th)) ifelse(is.finite(df$prev_anchored), pmax(0, 1 - abs(df$prev_anchored - th) / 0.15), 0)
                   else 0 },
        popn   = ifelse(is.finite(df$population), rk01(df$population), 0),
        nosvy  = ifelse(is.na(df$surveyed) | !df$surveyed, 1, 0))
      df$plan_score <- as.numeric(as.matrix(comp) %*% w[c("burden", "uncert", "thresh", "popn", "nosvy")])
      for (nm in names(comp)) df[[paste0("c_", nm)]] <- comp[[nm]] * w[[nm]]
      df$selected <- rank(-df$plan_score, ties.method = "first") <= min(input$k, nrow(df))
      df
    })

    output$plan_map <- renderLeaflet({
      df <- plan_data(); req(nrow(df) > 0)
      pal <- colorNumeric("YlGnBu", domain = range(df$plan_score, na.rm = TRUE), na.color = "#d9d9d9")
      labels <- sprintf(paste0("<strong>%s</strong>%s<br/>Priority under this objective: %s<br/>",
                               "burden %s | instability %s | near threshold %s | population %s | never surveyed %s"),
                        df$Admin2, ifelse(df$selected, " — PROPOSED VISIT", ""), fmt_num(df$plan_score, 2),
                        fmt_num(df$c_burden, 2), fmt_num(df$c_uncert, 2), fmt_num(df$c_thresh, 2),
                        fmt_num(df$c_popn, 2), fmt_num(df$c_nosvy, 2)) |> lapply(HTML)
      leaflet(df) |> addProviderTiles(providers$Esri.WorldGrayCanvas) |>
        addPolygons(fillColor = pal(df$plan_score), fillOpacity = 0.75,
                    color = ifelse(df$selected, "#C8641E", "#9aa0a6"),
                    weight = ifelse(df$selected, 3, 0.5), opacity = 1,
                    label = labels, labelOptions = labelOptions(textsize = "13px"),
                    highlightOptions = highlightOptions(weight = 4, color = "#333", bringToFront = TRUE)) |>
        addLegend(pal = pal, values = df$plan_score[is.finite(df$plan_score)], opacity = 0.75,
                  title = "Visit priority", position = "bottomright")
    })

    output$plan_table <- renderReactable({
      df <- sf::st_drop_geometry(plan_data()); req(nrow(df) > 0)
      s <- df[df$selected, ]
      s <- s[order(-s$plan_score), ]
      why <- vapply(seq_len(nrow(s)), function(i) {
        parts <- c(burden = s$c_burden[i], instability = s$c_uncert[i], `near threshold` = s$c_thresh[i],
                   population = s$c_popn[i], `never surveyed` = s$c_nosvy[i])
        paste(names(sort(parts, decreasing = TRUE))[1:2], collapse = " + ")
      }, "")
      t <- data.frame(District = s$Admin2, Region = s$Admin1,
                      `Mostly because` = why,
                      `Model rank` = ifelse(is.finite(s$rank_worst), sprintf("%d of %d", s$rank_worst, s$n_districts), "—"),
                      `Rank moves across refits` = ifelse(is.finite(s$rank_width), sprintf("±%.0f places", s$rank_width / 2), "—"),
                      `Surveyed before` = ifelse(!is.na(s$surveyed) & s$surveyed, "yes", "no"),
                      check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 15, searchable = FALSE)
    })

    output$plan_summary <- renderUI({
      df <- sf::st_drop_geometry(plan_data()); req(nrow(df) > 0)
      s <- df[df$selected, ]
      tagList(
        p(style = "font-size:0.9em;", sprintf("%d visits proposed: %d never-surveyed districts, %d of the model's worst fifth, about %s people in the target group.",
          nrow(s), sum(!s$surveyed | is.na(s$surveyed)), sum(s$worst_fifth, na.rm = TRUE), fmt_count(sum(s$population, na.rm = TRUE)))),
        p(style = "font-size:0.8em; color:#777;",
          "Instability is how far a district's rank moves when the model is refitted on resampled training data - ",
          "a stability measure (see How well it works), not a guarantee."))
    })

    output$borrowed <- renderUI({
      methods_note(
        "Where these ideas come from: choosing the next survey locations by prediction uncertainty and threshold",
        " proximity is adaptive geostatistical design (Chipeta et al. 2016, Spatial Statistics; used for rolling",
        " malaria surveys in Malawi, Kabaghe et al. 2017, PLOS ONE); batch selection against a prevalence threshold",
        " is the hotspot-finding design of Andrade-Pacheco et al. 2020 (Scientific Reports); classifying areas",
        " against elimination thresholds with model support is standard in neglected-tropical-disease programmes",
        " (Fronterre et al. 2020, J Infect Dis; Diggle et al. 2021, Trans R Soc Trop Med Hyg; Amoah et al. 2022,",
        " Int J Epidemiol). This tab adapts their district-selection logic; it does not implement their",
        " cluster-level geostatistical machinery.")
    })

    # ── 2. the size story (AR-01, reframed) ───────────────────────────────
    AR <- EV$design_summary
    design_lab <- c(A1_anchor_rank = "National anchor + ranking", A2_region_anchor_rank = "Regional anchors + ranking",
                    B_district_survey = "District survey", C_regional_survey = "Regional survey")
    design_col <- c("National anchor + ranking" = PROXY_COL, "Regional anchors + ranking" = "#7fb3bd",
                    "District survey" = SURVEY_COL, "Regional survey" = "#e0a878")
    observe({
      req(AR)
      rf <- unique(AR$rank_from); st <- unique(AR$set)
      updateSelectInput(session, "rank_from", choices = setNames(rf, c(prev = "Prevalence", level = "Biomarker level")[rf]),
                        selected = if ("prev" %in% rf) "prev" else rf[1])
      updateSelectInput(session, "set", choices = setNames(st, c(climate_soil = "Climate and soil (pre-registered)", full = "Full index")[st] %||% st),
                        selected = if ("climate_soil" %in% st) "climate_soil" else st[1])
      fr <- sort(unique(AR$fraction))
      updateSliderInput(session, "fraction", min = min(fr), max = max(fr), value = min(fr), step = min(diff(fr)))
    })
    dat <- reactive({
      req(AR, input$rank_from, input$set)
      d <- AR[AR$rank_from == input$rank_from & AR$set == input$set & AR$design %in% names(design_lab), ]
      d$Design <- design_lab[d$design]; d
    })
    output$mae <- renderPlotly({
      d <- dat(); validate(need(nrow(d) > 0, "No design rows for this choice."))
      plot_ly(d, x = ~fraction, y = ~mae, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%% of the sample: %.1f points", Design, 100 * fraction, mae), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the full survey sample", tickformat = ".0%"),
               yaxis = list(title = "Median district error (percentage points)", rangemode = "tozero"),
               shapes = list(list(type = "line", x0 = input$fraction, x1 = input$fraction, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(dash = "dot", color = "#888"))),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$capture <- renderPlotly({
      d <- dat(); validate(need(nrow(d) > 0 && "capture" %in% names(d), "No capture rows."))
      d <- d[is.finite(d$capture), ]
      plot_ly(d, x = ~fraction, y = ~capture, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%%: %.0f%% of burden in the picked fifth", Design, 100 * fraction, 100 * capture), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the full survey sample", tickformat = ".0%"),
               yaxis = list(title = "Burden in the worst-ranked fifth", tickformat = ".0%"),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$at_fraction <- renderUI({
      d <- dat(); req(nrow(d) > 0)
      fr <- d$fraction[which.min(abs(d$fraction - input$fraction))]
      s <- d[d$fraction == fr, ]
      # translate the fraction into people and clusters for a chosen country
      SM <- EV$survey_design_meta
      trans <- NULL
      if (!is.null(SM)) {
        sm <- SM[SM$country == meta$countries[[input$size_country %||% "ghana"]], ]
        if (nrow(sm)) {
          nr <- stats::median(sm$n_raw, na.rm = TRUE); np <- stats::median(sm$n_psu, na.rm = TRUE)
          deff <- stats::median(sm$deff_binary, na.rm = TRUE); if (!is.finite(deff)) deff <- 2.4
          # the national prevalence is the survey's primary product: what does
          # this sample size do to ITS precision? SE = sqrt(p(1-p) deff / n).
          pn <- idx_national$national_prev[idx_national$country_key == (input$size_country %||% "ghana")]
          pn <- stats::median(pn[is.finite(pn) & pn > 0.02], na.rm = TRUE); if (!is.finite(pn)) pn <- 0.2
          half <- 1.96 * sqrt(pn * (1 - pn) * deff / max(fr * nr, 1))
          trans <- tagList(
            p(style = "font-size:0.85em; color:#555;",
              sprintf("For a survey like %s's: about %s respondents per outcome group in roughly %s clusters (a full survey there ran %s in %s clusters; design effect near %s).",
                      meta$countries[[input$size_country %||% "ghana"]],
                      fmt_count(fr * nr), fmt_count(max(1, round(fr * np))), fmt_count(nr), fmt_count(np), fmt_num(deff, 1)))
            ,
            p(style = "font-size:0.85em; color:#555;",
              strong("The national prevalence - the survey's primary product - "),
              sprintf("keeps its validity at any of these sizes as long as the sample stays a national probability sample (the anchored design's sample is one), but its precision shrinks: at this size a prevalence near %s carries a 95%% interval of about ±%s points, against ±%s for the full survey.",
                      fmt_pct(pn, 0), fmt_num(100 * half, 1),
                      fmt_num(100 * 1.96 * sqrt(pn * (1 - pn) * deff / nr), 1))))
        }
      }
      tagList(
        h6(sprintf("At %s of a full survey's sample", fmt_pct(fr, 0))),
        tags$table(class = "table table-sm", style = "font-size:0.88em;",
                   tags$thead(tags$tr(tags$th("Design"), tags$th("Error, points"), tags$th("Ranking accuracy"), tags$th("Burden reached"))),
                   tags$tbody(lapply(seq_len(nrow(s)), function(i) tags$tr(tags$td(s$Design[i]), tags$td(fmt_num(s$mae[i], 1)), tags$td(fmt_num(s$spearman[i])), tags$td(fmt_pct(s$capture[i], 0)))))),
        trans)
    })

    # ── 3. SP-01 ───────────────────────────────────────────────────────────
    output$planner_curve <- renderPlotly({
      PV <- EV$planner_validation; validate(need(!is.null(PV), "Planner validation not built yet (script 64)."))
      lab <- c(random = "Chosen at random", pps = "Chosen by population", spread_model = "Spread across the model's gradient",
               extremes_model = "The model's extremes")
      cols <- c("Chosen at random" = "#8c8c8c", "Chosen by population" = "#e0a878",
                "Spread across the model's gradient" = PROXY_COL, "The model's extremes" = "#7fb3bd")
      d <- PV |> group_by(fraction, arm) |> summarise(m = mean(spearman, na.rm = TRUE), .groups = "drop")
      d$Rule <- lab[d$arm]
      t0 <- mean(PV$spearman_transport_only, na.rm = TRUE)
      plot_ly(d, x = ~fraction, y = ~m, color = ~Rule, colors = cols, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s, %.0f%% of districts surveyed: %.2f", Rule, 100 * fraction, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of districts the survey visits", tickformat = ".0%"),
               yaxis = list(title = "Ranking accuracy on the districts left out", rangemode = "tozero"),
               shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = t0, y1 = t0,
                                  line = list(dash = "dash", color = "#666"))),
               annotations = list(list(x = 0.02, y = t0, xref = "paper", yref = "y", xanchor = "left", yanchor = "bottom",
                                       text = "no in-country biomarkers at all (transported ranking)", showarrow = FALSE,
                                       font = list(size = 11, color = "#666"))),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$national <- renderPlotly({
      PV <- EV$planner_validation
      validate(need(!is.null(PV) && "nat_ci_pp" %in% names(PV), "National-estimate metrics not built yet (script 64 amendment)."))
      lab <- c(random = "Chosen at random", pps = "Chosen by population", spread_model = "Model-stratified (naive weights)",
               extremes_model = "The model's extremes")
      cols <- c("Chosen at random" = "#8c8c8c", "Chosen by population" = "#e0a878",
                "Model-stratified (naive weights)" = PROXY_COL, "Model-stratified (design weights)" = "#083D45",
                "The model's extremes" = "#7fb3bd")
      d <- PV |> group_by(fraction, arm) |>
        summarise(ci = mean(nat_ci_pp, na.rm = TRUE), bias = mean(abs(nat_bias_pp), na.rm = TRUE), .groups = "drop")
      d$Rule <- lab[d$arm]
      ds <- PV |> filter(arm == "spread_model") |> group_by(fraction) |>
        summarise(ci = mean(nat_strat_ci_pp, na.rm = TRUE), bias = mean(abs(nat_strat_bias_pp), na.rm = TRUE), .groups = "drop") |>
        mutate(Rule = "Model-stratified (design weights)", arm = "spread_strat")
      d <- bind_rows(d, ds)
      plot_ly(d, x = ~fraction, y = ~ci, color = ~Rule, colors = cols, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s, %.0f%% of districts:<br>national figure moves ±%.1f points across draws<br>typical bias %.1f points", Rule, 100 * fraction, ci / 2, bias),
              hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of districts the survey visits", tickformat = ".0%"),
               yaxis = list(title = "95% spread of the national figure (pp)", rangemode = "tozero"),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$national_bias <- renderUI({
      PV <- EV$planner_validation
      if (is.null(PV) || !"nat_bias_pp" %in% names(PV)) return(empty_state("Rebuild the data to see the national-estimate metrics."))
      at <- function(a, col) { v <- PV[[col]][PV$arm == a & PV$fraction == 0.5]; mean(abs(v), na.rm = TRUE) }
      tagList(
        p(style = "font-size:0.9em;", sprintf("At half the districts, the typical shift in the national figure is %s points when districts are chosen at random, %s under model-stratified choice with design weights, and %s when only the model's extremes are visited.",
          fmt_num(at("random", "nat_bias_pp"), 1), fmt_num(at("spread_model", "nat_strat_bias_pp"), 1), fmt_num(at("extremes_model", "nat_bias_pp"), 1))),
        p(style = "font-size:0.9em;", strong("Rule: "), "the model may decide WHICH districts anchor each stratum; it must never decide which districts count. Keep the sample a probability sample and the national number keeps its meaning."),
        p(style = "font-size:0.82em; color:#777;", "These are the between-district components; respondent-level noise adds the usual survey interval on top (see How small can a survey be)."))
    })
  })
}
