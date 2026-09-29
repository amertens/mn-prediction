# =============================================================================
# Module: What drives the estimate (concise version)
# =============================================================================

mod_importance_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Leading data layers", icon = bsicons::bs_icon("bar-chart-line"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome"), "Outcome", choices = NULL),
            radioButtons(ns("target"), "Outcome measured as", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
            p(style = "font-size:0.85em; color:#555;",
              "Right = pushes toward more deficiency. Grey bars: range when the model is re-estimated.",
              " * = direction differs in at least one country.")),
        plotlyOutput(ns("top"), height = "460px")),
      p(style = "font-size:0.85em; color:#666;", "These show where deficiency is likely, not what to change.")
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
          downloadButton(ns("sc_download"), "Download this table", class = "btn-sm btn-outline-primary")),
        div(reactableOutput(ns("search_table")),
            p(style = "font-size:0.85em; color:#666;", "Positive weight = more deficiency. Definitions in Technical notes.")))
    ),
    nav_panel(
      title = "Which data groups matter", icon = bsicons::bs_icon("layers"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Importance within a country and in a new country"),
             plotlyOutput(ns("domain_scatter"), height = "460px")),
        card(card_header("Effect of dropping a whole data source"),
             plotlyOutput(ns("source_ablation"), height = "460px"))),
      p(style = "font-size:0.85em; color:#666;", "Climate and soil help most in a new country; household-survey data help only within a country.")
    ),
    nav_panel(
      title = "Twenty public layers", icon = bsicons::bs_icon("list-ol"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome20"), "Outcome", choices = NULL),
            p(style = "font-size:0.85em; color:#555;",
              sprintf("A simple score from these 20 layers ranks a new country at %s (full model: %s). None needs a blood sample.",
                      fmt_num(Q$sparse20_tr), fmt_num(Q$tr)))),
        reactableOutput(ns("twenty")))
    ),
    nav_panel(
      title = "A shared pattern", icon = bsicons::bs_icon("arrow-down-up"),
      p(class = "lead", "Drier, poorer and more livestock-dependent districts with more child stunting rank worst for every deficiency."),
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
              text = ~sprintf("%s<br>%s<br>weight %+.2f", column, source, x), hoverinfo = "text", showlegend = FALSE) |>
        layout(xaxis = list(title = "Weight in the model (right = more deficiency)", zeroline = TRUE),
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
                      Group = dom_disp(d$domain),
                      Weight = round(d$beta_std, 3),
                      Range = ifelse(is.finite(d$lo), sprintf("%+.2f to %+.2f", d$lo, d$hi), "—"),
                      `Same direction (countries)` = if ("incountry_sign_agree" %in% names(d)) sprintf("%d of %d", d$incountry_sign_agree, d$incountry_fits) else "—",
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
      D$carry <- ifelse(is.finite(D$transport_cost) & D$transport_cost > 0.008, "helps in a new country",
                        ifelse(is.finite(D$transport_cost) & D$transport_cost < -0.004, "hurts in a new country", "little effect"))
      plot_ly(D, x = ~share_mean, y = ~transport_cost, type = "scatter", mode = "markers+text", text = ~disp, textposition = "top center",
              textfont = list(size = 10), color = ~carry, colors = c("helps in a new country" = PROXY_COL, "hurts in a new country" = SURVEY_COL, "little effect" = "#9a9a9a"),
              marker = list(size = 11), hovertext = ~sprintf("%s<br>share of model %.1f%%<br>accuracy lost in a new country if dropped: %+.3f",
                                                            disp, 100 * share_mean, transport_cost), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the model within surveyed countries", tickformat = ".0%"),
               yaxis = list(title = "Accuracy a new country loses if dropped", zeroline = TRUE),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$source_ablation <- renderPlotly({
      S <- EV$source_ablation; validate(need(!is.null(S), "Source table not built."))
      S <- S[S$target == "level" & S$n_cols >= 4, ]; S$source <- sub(" [(;].*$", "", S$source)
      S <- S[order(S$delta_drop), ]; S$source <- factor(S$source, levels = S$source)
      plot_ly(S, x = ~delta_drop, y = ~source, type = "bar", orientation = "h",
              marker = list(color = ifelse(S$delta_drop > 0, PROXY_COL, SURVEY_COL)),
              text = ~sprintf("%s: %+.3f", source, delta_drop), hoverinfo = "text") |>
        layout(xaxis = list(title = "Accuracy lost in a new country when dropped"),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })

    output$twenty <- renderReactable({
      req(IT, input$outcome20)
      d <- IT[IT$outcome == input$outcome20 & IT$target == "level" & IT$rank <= 20, ]
      d <- d[order(d$rank), ]
      t <- data.frame(Rank = d$rank, `Data layer` = pred_label(d$column), Column = d$column,
                      Direction = ifelse(d$beta > 0, "more deficiency", "less deficiency"),
                      Group = dom_disp(d$domain), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, searchable = TRUE,
                columns = list(Column = colDef(show = FALSE)),
                details = function(i) div(style = "padding:6px 24px; font-size:0.85em; color:#555;", t$Column[i]))
    })

    output$patterns <- renderReactable({
      IP <- EV$importance_patterns; validate(need(!is.null(IP), "Pattern table not built."))
      d <- IP[IP$target == "level" & IP$outcomes_in_top20 >= 3, ]
      d <- d[order(-d$outcomes_in_top20, d$mean_rank), ]
      t <- data.frame(`Data layer` = pred_label(d$column), Group = dom_disp(d$domain), `Outcomes (in top 20)` = d$outcomes_in_top20,
                      Directions = d$signs, check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 15)
    })
  })
}
