# =============================================================================
# Module: What drives the estimate
# =============================================================================
# The model is a weighted sum of the rank-transformed data layers, so its
# weights can be read directly for each layer. This tab shows those weights
# with their uncertainty: the leading layers per outcome with their range
# across re-estimates, a searchable table of every weight in every model, which
# data groups carry the model and which help in a new country, the twenty-layer
# score, and the pattern shared across outcomes.

mod_importance_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Leading data layers", icon = bsicons::bs_icon("bar-chart-line"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome"), "Outcome", choices = NULL),
            radioButtons(ns("target"), "Outcome measured as", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
            p(style = "font-size:0.85em; color:#555;",
              "Each dot shows how strongly a data layer moves a district's score in the model built on all four",
              " countries. Dots to the right push a district toward more deficiency. The grey bar shows where the value",
              " fell in 90% of re-estimates on different samples of districts. Colour shows the data source. A star",
              " marks a layer whose direction differs in at least one country's own model.")),
        plotlyOutput(ns("top"), height = "460px")),
      methods_note("These layers show where deficiency is likely, not what to change. For example, more livestock goes",
                   " with more iron and B12 deficiency at district level because in these four countries the",
                   " livestock-keeping areas are the drier and poorer ones. No single layer accounts for more than about",
                   " 2% of the model, and a layer whose bar crosses zero could easily be replaced by another.")
    ),
    nav_panel(
      title = "Search all layers", icon = bsicons::bs_icon("search"),
      layout_columns(col_widths = c(3, 9),
        div(
          selectInput(ns("sc_scope"), "Model built on", choices = c("All four countries" = "pooled",
                                                                   "One country only" = "country",
                                                                   "All countries except one" = "loco")),
          conditionalPanel(sprintf("input['%s'] != 'pooled'", ns("sc_scope")),
                           selectInput(ns("sc_country"), "Country", choices = NULL)),
          selectInput(ns("sc_outcome"), "Outcome", choices = NULL),
          radioButtons(ns("sc_target"), "Outcome measured as", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
          selectInput(ns("sc_domain"), "Data group", choices = NULL),
          downloadButton(ns("sc_download"), "Download this table", class = "btn-sm btn-outline-primary"),
          p(style = "font-size:0.8em; color:#777; margin-top:8px;",
            "Search by name, code or source. The table covers every model used in the tests: the model built on all",
            " four countries, each country's own model, and the models with one country left out.")),
        div(reactableOutput(ns("search_table")),
            methods_note("Weight: how much the layer moves a district's score for each step up the country's ranking of",
                         " that layer; positive means more deficiency. Share: the part of the model's differences between",
                         " districts that comes from this layer. Range: where the weight fell in 90% of re-estimates",
                         " (four-country model only). Same direction: in how many countries' own models, and in how many",
                         " models with one country left out, the layer pushes the same way. That is the most direct check",
                         " of whether a layer's role is consistent.")))
    ),
    nav_panel(
      title = "Which data groups matter", icon = bsicons::bs_icon("layers"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Importance within a country and in a new country"),
             plotlyOutput(ns("domain_scatter"), height = "460px")),
        card(card_header("Effect of dropping a whole data source"),
             plotlyOutput(ns("source_ablation"), height = "460px"))),
      methods_note(sprintf(paste("Horizontal axis: the share of the model that a data group accounts for within surveyed",
                                 "countries (average across outcomes, with its typical range across re-estimates). Vertical",
                                 "axis: how much accuracy a country left out of the model loses when that group is removed.",
                                 "Satellite, climate and soil data account for %s to %s of the model; climate and soil are",
                                 "the groups that help most in a new country. Groups built from household surveys, such as",
                                 "water and sanitation or education, help within a country but not in a new one. The DHS",
                                 "survey summaries showed that pattern more strongly, which is one reason they are not used."),
                           fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0)))
    ),
    nav_panel(
      title = "Twenty public layers", icon = bsicons::bs_icon("list-ol"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome20"), "Outcome", choices = NULL),
            p(style = "font-size:0.85em; color:#555;",
              sprintf(paste("A simple score that gives equal weight to the 20 layers with the largest weights ranks a new",
                            "country at %s, compared with %s for the full model. It is less accurate, but a country could",
                            "assemble these 20 layers itself, and none of them needs a blood sample."),
                      fmt_num(Q$sparse20_tr), fmt_num(Q$tr)))),
        reactableOutput(ns("twenty")))
    ),
    nav_panel(
      title = "A shared pattern", icon = bsicons::bs_icon("arrow-down-up"),
      p(class = "lead", "The same kinds of district rank worst for every deficiency."),
      p("The table lists data layers in the top 20 for three or more of the six outcomes, with their direction for each.",
        " Districts that are drier and grassier, less productive, more dependent on livestock, poorer and with more child",
        " stunting rank worst on every deficiency measured. Patterns specific to single nutrients sit on top of this shared one."),
      reactableOutput(ns("patterns"))
    )
  )
}

mod_importance_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    IT <- EV$importance_top
    IA <- EV$importance_all
    PW <- if (length(UE)) UE$pooled_weights else NULL
    observe({
      req(IT)
      ocs <- intersect(names(outcome_short), unique(IT$outcome))
      updateSelectInput(session, "outcome", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
      updateSelectInput(session, "outcome20", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
    })
    observe({
      req(IA)
      ocs <- intersect(names(outcome_short), unique(IA$outcome))
      updateSelectInput(session, "sc_outcome", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
      cts <- sort(unique(IA$country[IA$scope %in% c("country", "loco") & nzchar(IA$country %||% "")]))
      lbl <- ifelse(cts == "SierraLeone", "Sierra Leone", cts)
      updateSelectInput(session, "sc_country", choices = setNames(cts, lbl))
      dms <- sort(unique(IA$domain))
      updateSelectInput(session, "sc_domain", choices = c("All groups" = "", setNames(dms, dom_disp(dms))))
    })
    src_pal <- c("DHS household surveys" = "#C8641E", "MICS household surveys (public microdata)" = "#e08214", "MICS via WHO HEAT" = "#fdb863",
                 "Household budget surveys (public microdata)" = "#b35806", "World Bank RTFP market prices" = "#fee0b6",
                 "IHME modelled surfaces" = "#6a3d9a", "Malaria Atlas Project" = "#cab2d6",
                 "Earth Engine (climate, land, built environment)" = "#0F7B8A", "AlphaEarth satellite embedding" = "#7fcdbb",
                 "SoilGrids / iSDA soil" = "#8c510a", "MapSPAM crops" = "#33a02c", "Gridded Livestock of the World" = "#b15928",
                 "WHO ESPEN helminths" = "#9e9ac8", "WFP market prices" = "#fdbf6f")
    col_of <- function(s) { c <- src_pal[s]; c[is.na(c)] <- "#8c8c8c"; unname(c) }

    output$top <- renderPlotly({
      req(IT, input$outcome, input$target)
      d <- IT[IT$outcome == input$outcome & IT$target == input$target & IT$rank <= 10, ]
      validate(need(nrow(d) > 0, "No weights for this outcome."))
      d$source <- pred_source(d$column); d$label <- unique_labels(pred_label(d$column), d$column)
      d$label <- ifelse(d$incountry_sign_agree == d$incountry_fits, d$label, paste(d$label, "*"))
      d$lo <- d$hi <- d$med <- NA_real_
      if (!is.null(PW)) {
        u <- PW[PW$outcome == input$outcome & PW$target == input$target, ]
        j <- match(d$column, u$column)
        d$lo <- u$beta_lo[j]; d$hi <- u$beta_hi[j]; d$med <- u$beta_med[j]
      }
      use_ue <- all(is.finite(d$med))
      d$x <- if (use_ue) d$med else d$beta_std
      d <- d[order(d$x), ]; d$label <- factor(d$label, levels = d$label)
      p <- plot_ly(d)
      if (use_ue) p <- p |> add_segments(x = ~lo, xend = ~hi, y = ~label, yend = ~label,
                                         line = list(color = "#c9c9c9", width = 7), hoverinfo = "none", showlegend = FALSE)
      p |> add_markers(x = ~x, y = ~label, marker = list(color = col_of(d$source), size = 11),
              text = ~sprintf("%s<br>%s<br>%s<br>weight %+.2f%s, share of model %.1f%%<br>same direction in %d of %d countries' own models",
                              column, source, dom_disp(domain), x,
                              ifelse(is.finite(lo), sprintf(" (%+.2f to %+.2f when re-estimated)", lo, hi), ""),
                              100 * share, incountry_sign_agree, incountry_fits), hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Weight in the model (right = more deficiency; bar = 90% range when re-estimated)", zeroline = TRUE),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })

    search_data <- reactive({
      req(IA, input$sc_scope, input$sc_outcome, input$sc_target)
      d <- IA[IA$scope == input$sc_scope & IA$outcome == input$sc_outcome & IA$target == input$sc_target, ]
      if (input$sc_scope != "pooled" && nzchar(input$sc_country %||% "")) d <- d[d$country == input$sc_country, ]
      if (nzchar(input$sc_domain %||% "")) d <- d[d$domain == input$sc_domain, ]
      d <- d[is.finite(d$beta_std) & d$beta_std != 0, ]
      d <- d[order(d$rank), ]
      d$lo <- d$hi <- NA_real_
      if (!is.null(PW) && input$sc_scope == "pooled") {
        u <- PW[PW$outcome == input$sc_outcome & PW$target == input$sc_target, ]
        j <- match(d$column, u$column); d$lo <- u$beta_lo[j]; d$hi <- u$beta_hi[j]
      }
      d
    })

    output$search_table <- renderReactable({
      d <- search_data()
      validate(need(nrow(d) > 0, "No weights for this choice."))
      t <- data.frame(Rank = d$rank, `Data layer` = pred_label(d$column), Code = d$column,
                      Group = dom_disp(d$domain), Source = pred_source(d$column),
                      Weight = round(d$beta_std, 3),
                      `Range when re-estimated` = ifelse(is.finite(d$lo), sprintf("%+.2f to %+.2f", d$lo, d$hi), "—"),
                      Share = sprintf("%.1f%%", 100 * d$share),
                      `Same direction: countries` = if ("incountry_sign_agree" %in% names(d)) sprintf("%d of %d", d$incountry_sign_agree, d$incountry_fits) else "—",
                      `Same direction: one country left out` = if ("loco_sign_agree" %in% names(d)) sprintf("%d of %d", d$loco_sign_agree, d$loco_fits) else "—",
                      check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, searchable = TRUE, defaultPageSize = 15,
                columns = list(Code = colDef(show = FALSE),
                               Weight = colDef(style = function(v) if (is.numeric(v) && is.finite(v)) list(color = if (v > 0) "#b2182b" else "#2166ac"))),
                details = function(i) div(style = "padding:6px 24px; font-size:0.85em; color:#555;", t$Code[i]))
    })

    output$sc_download <- downloadHandler(
      filename = function() sprintf("model_weights_%s_%s_%s_%s.csv", input$sc_scope, input$sc_outcome, input$sc_target, Sys.Date()),
      content = function(file) write.csv(search_data(), file, row.names = FALSE))

    output$domain_scatter <- renderPlotly({
      D <- CAT$domains; validate(need(!is.null(D) && "share_mean" %in% names(D), "Data group table not built."))
      D <- D[is.finite(D$share_mean), ]
      D$disp <- dom_disp(D$domain)
      D$sh_lo <- D$sh_hi <- NA_real_
      if (length(UE) && !is.null(UE$pooled_domains)) {
        u <- UE$pooled_domains[UE$pooled_domains$target == "level", ] |>
          group_by(domain) |> summarise(lo = mean(share_lo, na.rm = TRUE), hi = mean(share_hi, na.rm = TRUE), .groups = "drop")
        j <- match(D$domain, u$domain); D$sh_lo <- u$lo[j]; D$sh_hi <- u$hi[j]
      }
      D$carry <- ifelse(is.finite(D$transport_cost) & D$transport_cost > 0.008, "helps in a new country",
                        ifelse(is.finite(D$transport_cost) & D$transport_cost < -0.004, "hurts in a new country", "little effect"))
      plot_ly(D, x = ~share_mean, y = ~transport_cost, type = "scatter", mode = "markers+text", text = ~disp, textposition = "top center",
              textfont = list(size = 10), color = ~carry, colors = c("helps in a new country" = PROXY_COL, "hurts in a new country" = SURVEY_COL, "little effect" = "#9a9a9a"),
              error_x = ~list(array = sh_hi - share_mean, arrayminus = share_mean - sh_lo, thickness = 1, color = "#c9c9c9"),
              marker = list(size = 11), hovertext = ~sprintf("%s<br>%d layers, %d summary scores kept<br>share of model %.1f%% (typically %.1f to %.1f%% when re-estimated)<br>accuracy lost in a new country if dropped: %+.3f",
                                                            disp, n_columns, n_axes, 100 * share_mean, 100 * sh_lo, 100 * sh_hi, transport_cost), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the model within surveyed countries", tickformat = ".0%"),
               yaxis = list(title = "Accuracy a new country loses if the group is dropped", zeroline = TRUE),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$source_ablation <- renderPlotly({
      S <- EV$source_ablation; validate(need(!is.null(S), "Source table not built."))
      S <- S[S$target == "level" & S$n_cols >= 4, ]; S$source <- sub(" [(;].*$", "", S$source)
      S <- S[order(S$delta_drop), ]; S$source <- factor(S$source, levels = S$source)
      plot_ly(S, x = ~delta_drop, y = ~source, type = "bar", orientation = "h",
              marker = list(color = ifelse(S$delta_drop > 0, PROXY_COL, SURVEY_COL)),
              text = ~sprintf("%s: %d layers<br>accuracy in a new country %.3f with it, %.3f without (%+.3f)", source, n_cols, full, drop, delta_drop), hoverinfo = "text") |>
        layout(xaxis = list(title = "Accuracy lost in a new country when the source is dropped (negative = better without it)"),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })

    output$twenty <- renderReactable({
      req(IT, input$outcome20)
      d <- IT[IT$outcome == input$outcome20 & IT$target == "level" & IT$rank <= 20, ]
      d <- d[order(d$rank), ]
      d$lo <- d$hi <- NA_real_
      if (!is.null(PW)) {
        u <- PW[PW$outcome == input$outcome20 & PW$target == "level", ]
        j <- match(d$column, u$column); d$lo <- u$beta_lo[j]; d$hi <- u$beta_hi[j]
      }
      t <- data.frame(Rank = d$rank, `Data layer` = pred_label(d$column), Column = d$column,
                      Direction = ifelse(d$beta > 0, "more deficiency", "less deficiency"),
                      `Weight range when re-estimated` = ifelse(is.finite(d$lo), sprintf("%+.2f to %+.2f", d$lo, d$hi), "—"),
                      Source = pred_source(d$column), Group = dom_disp(d$domain),
                      `Same direction with each country left out` = sprintf("%d of %d", d$loco_sign_agree, d$loco_fits),
                      `Same direction in each country's own model` = sprintf("%d of %d", d$incountry_sign_agree, d$incountry_fits), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, searchable = TRUE,
                columns = list(Column = colDef(show = FALSE)),
                details = function(i) div(style = "padding:6px 24px; font-size:0.85em; color:#555;", t$Column[i]))
    })

    output$patterns <- renderReactable({
      IP <- EV$importance_patterns; validate(need(!is.null(IP), "Pattern table not built."))
      d <- IP[IP$target == "level" & IP$outcomes_in_top20 >= 3, ]
      d <- d[order(-d$outcomes_in_top20, d$mean_rank), ]
      t <- data.frame(`Data layer` = pred_label(d$column), Group = dom_disp(d$domain), `Outcomes (in top 20)` = d$outcomes_in_top20,
                      Directions = d$signs, `Average rank` = round(d$mean_rank, 1), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 15)
    })
  })
}
