# =============================================================================
# Module: Plan a survey (concise version)
# =============================================================================

PLANNER_PRESETS <- list(
  confirm = c(burden = 45, uncert = 10, thresh = 10, popn = 20, nosvy = 15),
  shrink  = c(burden = 10, uncert = 45, thresh = 25, popn = 20, nosvy = 0),
  cover   = c(burden = 15, uncert = 15, thresh = 0,  popn = 25, nosvy = 45)
)

mod_survey_design_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    id = ns("tabs"),
    nav_panel(
      title = "Where to sample next", icon = bsicons::bs_icon("geo-alt"),
      layout_sidebar(
        sidebar = sidebar(width = 320, title = "What should the next survey learn?",
          selectInput(ns("country"), "Country", choices = country_choices, selected = "ghana"),
          selectInput(ns("outcome"), "Outcome", choices = NULL),
          sliderInput(ns("k"), "Districts the survey can visit", min = 5, max = 40, value = 15, step = 1),
          radioButtons(ns("preset"), "Aim",
            choices = c("Confirm the worst-off districts" = "confirm",
                        "Improve estimates where the model is unsure" = "shrink",
                        "Reach districts never surveyed" = "cover",
                        "Set my own weights" = "custom"), selected = "shrink"),
          conditionalPanel(sprintf("input['%s'] == 'custom'", ns("preset")),
            sliderInput(ns("w_burden"), "Estimated burden", 0, 100, 20, 5),
            sliderInput(ns("w_uncert"), "Rank uncertainty", 0, 100, 30, 5),
            sliderInput(ns("w_thresh"), "Close to a WHO threshold", 0, 100, 20, 5),
            sliderInput(ns("w_popn"),   "Population", 0, 100, 20, 5),
            sliderInput(ns("w_nosvy"),  "Never surveyed", 0, 100, 10, 5)),
          hr(),
          uiOutput(ns("plan_summary"))),
        layout_columns(col_widths = c(7, 5),
          card(full_screen = TRUE,
            card_header("Districts proposed for a survey visit"),
            card_body(padding = 0, leafletOutput(ns("plan_map"), height = "560px"))),
          card(card_header("Why each district was chosen"),
            card_body(reactableOutput(ns("plan_table")),
              p(style = "font-size:0.85em; color:#666;",
                "Choosing clusters and sample sizes still needs a survey statistician. How the scores work: Technical notes."))))
        )
    ),
    nav_panel(
      title = "How small can a survey be", icon = bsicons::bs_icon("cash-coin"),
      layout_sidebar(
        sidebar = sidebar(width = 320, title = "Design choices",
          selectInput(ns("rank_from"), "Ranking learned from", choices = NULL),
          selectInput(ns("set"), "Data layers used", choices = NULL),
          selectInput(ns("size_country"), "Show sample sizes for", choices = country_choices, selected = "ghana"),
          sliderInput(ns("fraction"), "Share of a full survey's sample", min = 0.05, max = 1, value = 0.05, step = 0.05),
          hr(),
          uiOutput(ns("at_fraction"))),
        layout_columns(col_widths = c(12, 12),
          card(card_header("A small national sample plus the ranking gives lower district error than a small direct survey"),
            card_body(plotlyOutput(ns("mae"), height = "340px"),
              p(style = "font-size:0.85em; color:#666;",
                sprintf("At 5%% of a full survey's sample: off by %s points, against %s for a district survey of the same size.",
                        fmt_num(Q$ar_a1, 1), fmt_num(Q$ar_b, 1))))),
          card(card_header("But a district survey reaches more deficient people"),
            card_body(plotlyOutput(ns("capture"), height = "280px"))))
      )
    ),
    nav_panel(
      title = "Does choosing districts well matter?", icon = bsicons::bs_icon("diagram-3"),
      layout_columns(col_widths = c(8, 4),
        card(card_header("Ranking of the districts not visited, by how the visited districts were chosen"),
          card_body(plotlyOutput(ns("planner_curve"), height = "340px"))),
        card(card_header("What it shows"),
          card_body(
            p(style = "font-size:0.9em;", sprintf("Spreading visits across the model's ranking did about as well as random choice (%s against %s at half the districts). Choosing by population did worse (%s).",
                                                  f2(Q$plan_spread), f2(Q$plan_random), f2(Q$plan_pps))),
            p(style = "font-size:0.9em;", sprintf("With a third of districts or fewer, the model from other countries (dashed line, %s) ranks the rest better.",
                                                  f2(Q$plan_transport)))))),
      layout_columns(col_widths = c(8, 4),
        card(card_header("What happens to the national prevalence"),
          card_body(plotlyOutput(ns("national"), height = "280px"))),
        card(card_header("Keeping the national estimate valid"),
          card_body(uiOutput(ns("national_bias")))))
    )
  )
}

mod_survey_design_server <- function(id) {
  moduleServer(id, function(input, output, session) {
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
      labels <- sprintf("<strong>%s</strong>%s<br/>Priority under this aim: %s",
                        df$Admin2, ifelse(df$selected, " (proposed visit)", ""), fmt_num(df$plan_score, 2)) |> lapply(HTML)
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
        parts <- c(burden = s$c_burden[i], `rank uncertainty` = s$c_uncert[i], `near threshold` = s$c_thresh[i],
                   population = s$c_popn[i], `never surveyed` = s$c_nosvy[i])
        paste(names(sort(parts, decreasing = TRUE))[1:2], collapse = " + ")
      }, "")
      t <- data.frame(District = s$Admin2, `Mainly because` = why,
                      `Model rank` = ifelse(is.finite(s$rank_worst), sprintf("%d of %d", s$rank_worst, s$n_districts), "—"),
                      `Surveyed before` = ifelse(!is.na(s$surveyed) & s$surveyed, "yes", "no"),
                      check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 15, searchable = FALSE)
    })

    output$plan_summary <- renderUI({
      df <- sf::st_drop_geometry(plan_data()); req(nrow(df) > 0)
      s <- df[df$selected, ]
      p(style = "font-size:0.9em;", sprintf("%d visits: %d never surveyed, %d in the worst fifth.",
        nrow(s), sum(!s$surveyed | is.na(s$surveyed)), sum(s$worst_fifth, na.rm = TRUE)))
    })

    AR <- EV$design_summary
    design_lab <- c(A1_anchor_rank = "National sample + ranking", A2_region_anchor_rank = "Regional samples + ranking",
                    B_district_survey = "District survey", C_regional_survey = "Regional survey")
    design_col <- c("National sample + ranking" = PROXY_COL, "Regional samples + ranking" = "#7fb3bd",
                    "District survey" = SURVEY_COL, "Regional survey" = "#e0a878")
    observe({
      req(AR)
      rf <- unique(AR$rank_from); st <- unique(AR$set)
      updateSelectInput(session, "rank_from", choices = setNames(rf, c(prev = "Prevalence", level = "Biomarker level")[rf]),
                        selected = if ("prev" %in% rf) "prev" else rf[1])
      updateSelectInput(session, "set", choices = setNames(st, c(climate_soil = "Climate and soil only", full = "All data layers")[st] %||% st),
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
      d <- dat(); validate(need(nrow(d) > 0, "No results for this choice."))
      plot_ly(d, x = ~fraction, y = ~mae, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%% of the sample: %.1f points", Design, 100 * fraction, mae), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of a full survey's sample", tickformat = ".0%"),
               yaxis = list(title = "District error (percentage points)", rangemode = "tozero"),
               shapes = list(list(type = "line", x0 = input$fraction, x1 = input$fraction, y0 = 0, y1 = 1, yref = "paper",
                                  line = list(dash = "dot", color = "#888"))),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$capture <- renderPlotly({
      d <- dat(); validate(need(nrow(d) > 0 && "capture" %in% names(d), "No results."))
      d <- d[is.finite(d$capture), ]
      plot_ly(d, x = ~fraction, y = ~capture, color = ~Design, colors = design_col, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s at %.0f%%: %.0f%%", Design, 100 * fraction, 100 * capture), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of a full survey's sample", tickformat = ".0%"),
               yaxis = list(title = "Deficient people in the worst-ranked fifth", tickformat = ".0%"),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })
    output$at_fraction <- renderUI({
      d <- dat(); req(nrow(d) > 0)
      fr <- d$fraction[which.min(abs(d$fraction - input$fraction))]
      s <- d[d$fraction == fr, ]
      SM <- EV$survey_design_meta
      trans <- NULL
      if (!is.null(SM)) {
        sm <- SM[SM$country == meta$countries[[input$size_country %||% "ghana"]], ]
        if (nrow(sm)) {
          nr <- stats::median(sm$n_raw, na.rm = TRUE); np <- stats::median(sm$n_psu, na.rm = TRUE)
          deff <- stats::median(sm$deff_binary, na.rm = TRUE); if (!is.finite(deff)) deff <- 2.4
          pn <- idx_national$national_prev[idx_national$country_key == (input$size_country %||% "ghana")]
          pn <- stats::median(pn[is.finite(pn) & pn > 0.02], na.rm = TRUE); if (!is.finite(pn)) pn <- 0.2
          half <- 1.96 * sqrt(pn * (1 - pn) * deff / max(fr * nr, 1))
          trans <- p(style = "font-size:0.85em; color:#555;",
            sprintf("About %s people in %s clusters. National prevalence: ±%s points (full survey: ±%s).",
                    fmt_count(fr * nr), fmt_count(max(1, round(fr * np))), fmt_num(100 * half, 1),
                    fmt_num(100 * 1.96 * sqrt(pn * (1 - pn) * deff / nr), 1)))
        }
      }
      tagList(
        h6(sprintf("At %s of a full survey's sample", fmt_pct(fr, 0))),
        tags$table(class = "table table-sm", style = "font-size:0.88em;",
                   tags$thead(tags$tr(tags$th("Design"), tags$th("Error (points)"), tags$th("Deficient people reached"))),
                   tags$tbody(lapply(seq_len(nrow(s)), function(i) tags$tr(tags$td(s$Design[i]), tags$td(fmt_num(s$mae[i], 1)), tags$td(fmt_pct(s$capture[i], 0)))))),
        trans)
    })

    output$planner_curve <- renderPlotly({
      PV <- EV$planner_validation; validate(need(!is.null(PV), "District-choice test not built yet."))
      lab <- c(random = "Chosen at random", pps = "Most populous", spread_model = "Spread across the model's ranking",
               extremes_model = "Highest- and lowest-ranked only")
      cols <- c("Chosen at random" = "#8c8c8c", "Most populous" = "#e0a878",
                "Spread across the model's ranking" = PROXY_COL, "Highest- and lowest-ranked only" = "#7fb3bd")
      d <- PV |> group_by(fraction, arm) |> summarise(m = mean(spearman, na.rm = TRUE), .groups = "drop")
      d$Rule <- lab[d$arm]
      t0 <- mean(PV$spearman_transport_only, na.rm = TRUE)
      plot_ly(d, x = ~fraction, y = ~m, color = ~Rule, colors = cols, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s, %.0f%% visited: %.2f", Rule, 100 * fraction, m), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of districts visited", tickformat = ".0%"),
               yaxis = list(title = "Ranking accuracy, districts not visited", rangemode = "tozero"),
               shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = t0, y1 = t0,
                                  line = list(dash = "dash", color = "#666"))),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$national <- renderPlotly({
      PV <- EV$planner_validation
      validate(need(!is.null(PV) && "nat_ci_pp" %in% names(PV), "National-estimate results not built yet."))
      lab <- c(random = "Chosen at random", pps = "Most populous", spread_model = "Spread across the ranking, unweighted",
               extremes_model = "Highest- and lowest-ranked only")
      cols <- c("Chosen at random" = "#8c8c8c", "Most populous" = "#e0a878",
                "Spread across the ranking, unweighted" = "#7fb3bd", "Spread across the ranking, weighted" = PROXY_COL,
                "Highest- and lowest-ranked only" = "#b0b0b0")
      d <- PV |> group_by(fraction, arm) |>
        summarise(ci = mean(nat_ci_pp, na.rm = TRUE), bias = mean(abs(nat_bias_pp), na.rm = TRUE), .groups = "drop")
      d$Rule <- lab[d$arm]
      ds <- PV |> filter(arm == "spread_model") |> group_by(fraction) |>
        summarise(ci = mean(nat_strat_ci_pp, na.rm = TRUE), bias = mean(abs(nat_strat_bias_pp), na.rm = TRUE), .groups = "drop") |>
        mutate(Rule = "Spread across the ranking, weighted", arm = "spread_strat")
      d <- bind_rows(d, ds)
      plot_ly(d, x = ~fraction, y = ~ci, color = ~Rule, colors = cols, type = "scatter", mode = "lines+markers",
              text = ~sprintf("%s, %.0f%% visited: ±%.1f points, bias %.1f", Rule, 100 * fraction, ci, bias), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of districts visited", tickformat = ".0%"),
               yaxis = list(title = "± points, national estimate", rangemode = "tozero"),
               legend = list(orientation = "h", y = -0.25), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$national_bias <- renderUI({
      PV <- EV$planner_validation
      if (is.null(PV) || !"nat_bias_pp" %in% names(PV)) return(empty_state("Rebuild the data to see the national-estimate results."))
      at <- function(a, col) { v <- PV[[col]][PV$arm == a & PV$fraction == 0.5]; mean(abs(v), na.rm = TRUE) }
      tagList(
        p(style = "font-size:0.9em;", sprintf("At half the districts the national estimate is off by %s points with random choice and %s when only the extremes are visited.",
          fmt_num(at("random", "nat_bias_pp"), 1), fmt_num(at("extremes_model", "nat_bias_pp"), 1))),
        p(style = "font-size:0.9em;", "Let the model guide choices within a probability sample, not replace it."))
    })
  })
}
