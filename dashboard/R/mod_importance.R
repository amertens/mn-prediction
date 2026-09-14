# =============================================================================
# Module: What drives the estimate
# =============================================================================
# The index is linear in the rank-normalised predictors, so its weights project
# back exactly onto the columns. This tab shows that projection: the leading
# predictors per outcome, which kinds of data carry the model and which travel
# to a new country, the twenty-layer composite, and the gradient that recurs.

mod_importance_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    nav_panel(
      title = "Leading predictors", icon = bsicons::bs_icon("bar-chart-line"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome"), "Outcome", choices = NULL),
            radioButtons(ns("target"), "Target", choices = c("Biomarker level" = "level", "Prevalence" = "prev"), selected = "level"),
            p(style = "font-size:0.85em; color:#555;",
              "Bar length is the predictor's exact weight in the four-country fit, per unit of its within-country rank.",
              " Right means the district ranks worse. Colour is the kind of data. A star marks a sign that is not",
              " reproduced in every country's own fit.")),
        plotlyOutput(ns("top"), height = "460px")),
      methods_note("These weights say where deficiency is, not what to change. Livestock density weighting toward",
                   " more iron and B12 deficiency is ecological: in these four countries the pastoral zones are the",
                   " dry, poor zones. No single predictor carries more than about 2 percent of the model.")
    ),
    nav_panel(
      title = "Which data carry it", icon = bsicons::bs_icon("layers"),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Inside a country against across borders"),
             plotlyOutput(ns("domain_scatter"), height = "460px")),
        card(card_header("Dropping a whole data source"),
             plotlyOutput(ns("source_ablation"), height = "460px"))),
      methods_note(sprintf(paste("Horizontal axis: the share of the model a data group carries inside a surveyed country (pooled fit,",
                                 "mean over outcomes). Vertical: what a country the model has never seen loses when that group is",
                                 "removed. Satellite imagery, climate and soil carry %s to %s of the model; climate and soil are",
                                 "what a new country needs, and the household-survey aggregates sit below zero: they help inside",
                                 "a country and hurt in a new one."), fmt_pct(Q$env_lo, 0), fmt_pct(Q$env_hi, 0)))
    ),
    nav_panel(
      title = "Twenty public layers", icon = bsicons::bs_icon("list-ol"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput(ns("outcome20"), "Outcome", choices = NULL),
            p(style = "font-size:0.85em; color:#555;",
              sprintf(paste("An equal-weight composite of the twenty predictors with the largest weights, membership chosen",
                            "inside each training fold, ranks a new country at %s against %s for the full %d-column index.",
                            "A country can assemble this list itself; none of it needs a blood sample."),
                      fmt_num(Q$sparse20_tr), fmt_num(Q$tr), Q$n_predictors))),
        reactableOutput(ns("twenty")))
    ),
    nav_panel(
      title = "One gradient", icon = bsicons::bs_icon("arrow-down-up"),
      p(class = "lead", "The same kind of district ranks worst on every deficiency."),
      p("Predictors in the leading twenty of three or more of the six outcomes, with the sign in each. Districts that are",
        " drier and grassier, less productive, more pastoral, poorer and with more stunted children rank worst on every",
        " deficiency measured. Outcome-specific signals sit beneath that gradient."),
      reactableOutput(ns("patterns"))
    )
  )
}

mod_importance_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    IT <- EV$importance_top
    observe({
      req(IT)
      ocs <- intersect(names(outcome_short), unique(IT$outcome))
      updateSelectInput(session, "outcome", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
      updateSelectInput(session, "outcome20", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
    })
    src_pal <- c("DHS / MICS household surveys" = "#C8641E", "IHME modelled surfaces" = "#6a3d9a", "Malaria Atlas Project" = "#cab2d6",
                 "Earth Engine (climate, land, built environment)" = "#0F7B8A", "AlphaEarth satellite embedding" = "#7fcdbb",
                 "SoilGrids / iSDA soil" = "#8c510a", "MapSPAM crops" = "#33a02c", "Gridded Livestock of the World" = "#b15928",
                 "WHO ESPEN helminths" = "#9e9ac8", "WFP market prices" = "#fdbf6f")
    col_of <- function(s) { c <- src_pal[s]; c[is.na(c)] <- "#8c8c8c"; unname(c) }

    output$top <- renderPlotly({
      req(IT, input$outcome, input$target)
      d <- IT[IT$outcome == input$outcome & IT$target == input$target & IT$rank <= 10, ]
      validate(need(nrow(d) > 0, "No importance rows for this outcome."))
      d$source <- pred_source(d$column); d$label <- unique_labels(pred_label(d$column), d$column)
      d$label <- ifelse(d$incountry_sign_agree == d$incountry_fits, d$label, paste(d$label, "*"))
      d <- d[order(d$beta_std), ]; d$label <- factor(d$label, levels = d$label)
      plot_ly(d, x = ~beta_std, y = ~label, type = "bar", orientation = "h", marker = list(color = col_of(d$source)),
              text = ~sprintf("%s<br>%s<br>%s<br>weight %+.2f, share of model %.1f%%<br>same sign in %d of %d countries' own fits",
                              column, source, domain, beta_std, 100 * share, incountry_sign_agree, incountry_fits), hoverinfo = "text") |>
        layout(xaxis = list(title = "Weight in the index (right = more deficiency)"), yaxis = list(title = ""),
               margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })

    output$domain_scatter <- renderPlotly({
      D <- CAT$domains; validate(need(!is.null(D) && "share_mean" %in% names(D), "Domain table not built."))
      D <- D[is.finite(D$share_mean), ]
      D$carry <- ifelse(is.finite(D$transport_cost) & D$transport_cost > 0.008, "carries a new country",
                        ifelse(is.finite(D$transport_cost) & D$transport_cost < -0.004, "holds a new country back", "little effect"))
      plot_ly(D, x = ~share_mean, y = ~transport_cost, type = "scatter", mode = "markers+text", text = ~domain, textposition = "top center",
              textfont = list(size = 10), color = ~carry, colors = c("carries a new country" = PROXY_COL, "holds a new country back" = SURVEY_COL, "little effect" = "#9a9a9a"),
              marker = list(size = 11), hovertext = ~sprintf("%s<br>%d columns, %d axes kept<br>share of model %.1f%%<br>transport cost when dropped %+.3f",
                                                            domain, n_columns, n_axes, 100 * share_mean, transport_cost), hoverinfo = "text") |>
        layout(xaxis = list(title = "Share of the model inside a surveyed country", tickformat = ".0%"),
               yaxis = list(title = "Accuracy a new country loses if the group is dropped", zeroline = TRUE),
               legend = list(orientation = "h", y = -0.2), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$source_ablation <- renderPlotly({
      S <- EV$source_ablation; validate(need(!is.null(S), "Source ablation not built."))
      S <- S[S$target == "level" & S$n_cols >= 4, ]; S$source <- sub(" [(;].*$", "", S$source)
      S <- S[order(S$delta_drop), ]; S$source <- factor(S$source, levels = S$source)
      plot_ly(S, x = ~delta_drop, y = ~source, type = "bar", orientation = "h",
              marker = list(color = ifelse(S$delta_drop > 0, PROXY_COL, SURVEY_COL)),
              text = ~sprintf("%s: %d columns<br>transport %.3f with, %.3f without (%+.3f)", source, n_cols, full, drop, delta_drop), hoverinfo = "text") |>
        layout(xaxis = list(title = "Accuracy lost in a new country when the source is dropped (negative = it helps to drop)"),
               yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 50)) |> config(displayModeBar = FALSE)
    })

    output$twenty <- renderReactable({
      req(IT, input$outcome20)
      d <- IT[IT$outcome == input$outcome20 & IT$target == "level" & IT$rank <= 20, ]
      d <- d[order(d$rank), ]
      t <- data.frame(Rank = d$rank, Predictor = pred_label(d$column), Column = d$column,
                      Direction = ifelse(d$beta > 0, "more deficiency", "less deficiency"),
                      Source = pred_source(d$column), Domain = d$domain,
                      `Same sign, every held-out country` = sprintf("%d of %d", d$loco_sign_agree, d$loco_fits),
                      `Same sign, each country alone` = sprintf("%d of %d", d$incountry_sign_agree, d$incountry_fits), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, pagination = FALSE, searchable = TRUE,
                columns = list(Column = colDef(show = FALSE)),
                details = function(i) div(style = "padding:6px 24px; font-size:0.85em; color:#555;", t$Column[i]))
    })

    output$patterns <- renderReactable({
      IP <- EV$importance_patterns; validate(need(!is.null(IP), "Pattern table not built."))
      d <- IP[IP$target == "level" & IP$outcomes_in_top20 >= 3, ]
      d <- d[order(-d$outcomes_in_top20, d$mean_rank), ]
      t <- data.frame(Predictor = pred_label(d$column), Domain = d$domain, `Outcomes (in top 20)` = d$outcomes_in_top20,
                      Signs = d$signs, `Mean rank` = round(d$mean_rank, 1), check.names = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, defaultPageSize = 15)
    })
  })
}
