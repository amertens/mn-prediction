# =============================================================================
# Module: Predictor catalogue
# =============================================================================
# Every predictor the model can see, browsed by data source or by conceptual
# domain: what it is, where it comes from, which countries have it, how much
# weight the index gives it for each outcome, whether its district association
# points the same way in every country, and a small map of its values.
# Definitions come from the variable annotation sheet; where the sheet has none
# yet the row says so, which makes this tab the worklist for finishing it.

mod_catalogue_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330, title = "Browse the predictors",
      radioButtons(ns("by"), "Browse by", choices = c("Conceptual domain" = "domain", "Data source" = "source"), inline = TRUE),
      selectizeInput(ns("pick"), "Show", choices = NULL, multiple = TRUE, options = list(placeholder = "all")),
      selectInput(ns("outcome"), "Weight shown for", choices = NULL),
      checkboxInput(ns("only_defined"), "Only predictors with a plain-language name", FALSE),
      checkboxInput(ns("only_composite"), "Only members of a twenty-layer composite", FALSE),
      checkboxInput(ns("only_travel"), "Only climate and soil (the layers that travel)", FALSE),
      hr(),
      uiOutput(ns("selection_summary")),
      downloadButton(ns("download"), "Download this table", class = "btn-sm btn-outline-primary")
    ),
    layout_columns(
      col_widths = c(12, 12),
      card(card_header(textOutput(ns("table_title"), inline = TRUE)),
           card_body(reactableOutput(ns("table")),
                     methods_note("Weight is the predictor's exact contribution per unit of its within-country rank in the",
                                  " four-country fit for the chosen outcome; positive means more deficiency. Replicated counts",
                                  " the outcomes for which the predictor's district association carries the same sign in every",
                                  " country that measured them. Mechanism is the annotation sheet's template for the predictor's",
                                  " group, not a per-variable definition; predictors without a plain-language name show a",
                                  " cleaned code. Click a row for detail and a map."))),
      card(card_header("Selected predictor"), card_body(uiOutput(ns("detail"))))
    )
  )
}

mod_catalogue_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    V <- if (!is.null(CAT)) CAT$variables else NULL

    observe({
      req(V)
      ch <- if (input$by == "domain") sort(unique(V$domain)) else sort(unique(V$source_label))
      updateSelectizeInput(session, "pick", choices = ch, selected = character(0))
    })
    observe({
      req(CAT$weights)
      ocs <- intersect(names(outcome_short), unique(CAT$weights$outcome))
      updateSelectInput(session, "outcome", choices = setNames(ocs, outcome_short[ocs]), selected = "child_iron")
    })

    filtered <- reactive({
      req(V)
      d <- V
      if (length(input$pick)) d <- d[(if (input$by == "domain") d$domain else d$source_label) %in% input$pick, ]
      if (isTRUE(input$only_defined)) d <- d[!is.na(d$plain_name), ]
      if (isTRUE(input$only_composite)) d <- d[!is.na(d$composite_outcomes), ]
      if (isTRUE(input$only_travel)) d <- d[d$climate_soil, ]
      if (!is.null(CAT$weights) && nzchar(input$outcome %||% "")) {
        w <- CAT$weights[CAT$weights$outcome == input$outcome & CAT$weights$target == "level", ]
        i <- match(d$column, w$column); d$weight <- w$beta_std[i]; d$weight_rank <- w$rank[i]
      } else { d$weight <- NA_real_; d$weight_rank <- NA_integer_ }
      d[order(-abs(d$weight), d$column), ]
    })

    output$table_title <- renderText({
      d <- filtered(); sprintf("%d predictors%s", nrow(d), if (length(input$pick)) paste0(" in ", paste(input$pick, collapse = ", ")) else "")
    })

    output$selection_summary <- renderUI({
      d <- filtered(); req(nrow(d) > 0)
      lines <- list(p(sprintf("%d predictors from %d sources across %d domains; %d have a plain-language name, the rest show a cleaned code.",
                              nrow(d), length(unique(d$source_label)), length(unique(d$domain)), sum(!is.na(d$plain_name))), style = "font-size:0.9em;"))
      if (input$by == "domain" && length(input$pick) && !is.null(CAT$domains)) {
        dm <- CAT$domains[CAT$domains$domain %in% input$pick, ]
        if (nrow(dm) && "share_mean" %in% names(dm))
          lines <- c(lines, list(p(sprintf("Share of the model inside a country: %s (mean over outcomes). Accuracy a new country loses if dropped: %s.",
                                           fmt_pct(sum(dm$share_mean, na.rm = TRUE), 0),
                                           if (any(is.finite(dm$transport_cost))) fmt_num(sum(dm$transport_cost, na.rm = TRUE), 3) else "not measured"),
                                   style = "font-size:0.9em;")))
      }
      lines <- c(lines, list(p(sprintf("%d are members of a twenty-layer composite; %d are climate or soil layers.",
                                       sum(!is.na(d$composite_outcomes)), sum(d$climate_soil)), style = "font-size:0.9em;")))
      tagList(lines)
    })

    output$table <- renderReactable({
      d <- filtered(); req(nrow(d) > 0)
      t <- data.frame(Predictor = d$label, Code = d$column,
                      Mechanism = ifelse(is.na(d$mechanism), "", d$mechanism),
                      Unit = ifelse(is.na(d$unit), "", d$unit),
                      Domain = d$domain, Source = d$source_label,
                      Countries = d$n_countries, Complete = d$completeness,
                      Weight = d$weight, Rank = d$weight_rank,
                      Replicated = ifelse(is.na(d$n_replicated), 0L, d$n_replicated),
                      Composite = ifelse(is.na(d$composite_outcomes), "", gsub("_", " ", d$composite_outcomes)),
                      Travels = ifelse(d$climate_soil, "yes", ""), check.names = FALSE, stringsAsFactors = FALSE)
      reactable(t, compact = TRUE, striped = TRUE, searchable = TRUE, filterable = TRUE, defaultPageSize = 15,
                selection = "single", onClick = "select", highlight = TRUE,
                columns = list(
                  Predictor = colDef(minWidth = 170, style = function(v, i) if (is.na(d$plain_name[i])) list(color = "#666", fontStyle = "italic") else NULL),
                  Code = colDef(minWidth = 150, style = list(fontFamily = "monospace", fontSize = "0.8em")),
                  Mechanism = colDef(minWidth = 240, style = list(color = "#666", fontSize = "0.9em")),
                  Unit = colDef(width = 90), Countries = colDef(width = 90, align = "center"),
                  Complete = colDef(width = 90, format = colFormat(percent = TRUE, digits = 0)),
                  Weight = colDef(width = 90, format = colFormat(digits = 2),
                                  style = function(v) if (is.na(v)) NULL else list(color = if (v > 0) "#b2182b" else "#2166ac")),
                  Rank = colDef(width = 70), Replicated = colDef(width = 100, align = "center"),
                  Composite = colDef(minWidth = 150), Travels = colDef(width = 80, align = "center")))
    })

    selected_col <- reactive({
      s <- getReactableState("table", "selected", session); d <- filtered()
      if (is.null(s) || !length(s) || s > nrow(d)) return(NULL)
      d$column[s]
    })

    output$detail <- renderUI({
      col <- selected_col()
      if (is.null(col)) return(p(em("Click a row above to see the predictor's definition, its weight for every outcome, and a map of its values."), style = "color:#888;"))
      v <- V[V$column == col, ][1, ]
      tagList(
        h5(v$label, tags$small(style = "color:#888; font-family:monospace; margin-left:8px;", col), style = "margin-top:0;"),
        if (is.na(v$plain_name)) p(em("No plain-language name yet; the code is shown cleaned. The mechanism below is the annotation sheet's template for this predictor's group.")),
        if (!is.na(v$mechanism)) p(strong("Mechanism (group template): "), v$mechanism),
        tags$dl(class = "row", style = "font-size:0.9em;",
                tags$dt(class = "col-sm-3", "Domain"), tags$dd(class = "col-sm-9", v$domain),
                tags$dt(class = "col-sm-3", "Source"), tags$dd(class = "col-sm-9", sprintf("%s (%s)", v$source_label, v$source)),
                tags$dt(class = "col-sm-3", "Unit"), tags$dd(class = "col-sm-9", ifelse(is.na(v$unit), "not recorded", v$unit)),
                tags$dt(class = "col-sm-3", "Time"), tags$dd(class = "col-sm-9", ifelse(is.na(v$temporal_kind), "not recorded", v$temporal_kind)),
                tags$dt(class = "col-sm-3", "Countries"), tags$dd(class = "col-sm-9", sprintf("%s (%s of districts have a value)", v$countries, fmt_pct(v$completeness, 0))),
                tags$dt(class = "col-sm-3", "Range"), tags$dd(class = "col-sm-9", ifelse(is.na(v$value_range), "", v$value_range)),
                if (!is.na(v$coverage_note)) tagList(tags$dt(class = "col-sm-3", "Coverage"), tags$dd(class = "col-sm-9", v$coverage_note)),
                if (!is.na(v$source_note)) tagList(tags$dt(class = "col-sm-3", "Note"), tags$dd(class = "col-sm-9", v$source_note)),
                if (!is.na(v$composite_outcomes)) tagList(tags$dt(class = "col-sm-3", "Composite"), tags$dd(class = "col-sm-9", paste("Member of the twenty-layer composite for", gsub("_", " ", v$composite_outcomes))))),
        layout_columns(col_widths = c(6, 6),
                       card(card_header("Weight in the index, by outcome"), plotlyOutput(ns("weights_plot"), height = "260px")),
                       card(card_header("District association, by outcome"), plotlyOutput(ns("signal_plot"), height = "260px"))),
        layout_columns(col_widths = c(3, 9),
                       selectInput(ns("map_country"), "Map the predictor in", choices = country_choices, selected = "ghana"),
                       leafletOutput(ns("mini_map"), height = "360px")),
        tags$small(style = "color:#777;", "The map shows the value the model sees: the predictor's rank within the country, on a normal scale, darker = higher. Grey = no value.")
      )
    })

    output$weights_plot <- renderPlotly({
      col <- selected_col(); req(col, CAT$weights)
      w <- CAT$weights[CAT$weights$column == col & CAT$weights$target == "level", ]
      validate(need(nrow(w) > 0, "No weight for this predictor."))
      w$label <- outcome_short[w$outcome]
      plot_ly(w, x = ~beta_std, y = ~label, type = "bar", orientation = "h", marker = list(color = ifelse(w$beta_std > 0, "#b2182b", "#2166ac")),
              text = ~sprintf("%s<br>weight %+.2f, rank %d of %d<br>share of model %.2f%%", label, beta_std, rank, Q$n_predictors, 100 * share), hoverinfo = "text") |>
        layout(xaxis = list(title = "Weight (right = more deficiency)"), yaxis = list(title = ""), margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$signal_plot <- renderPlotly({
      col <- selected_col(); req(col)
      s <- CAT$signal; validate(need(!is.null(s), "Signal table not built."))
      s <- s[s$column == col, ]; validate(need(nrow(s) > 0, "This predictor was not in the association scan."))
      s$label <- outcome_short[s$outcome]; s$rep <- sprintf("%d of %d countries agree", s$sign_agree, s$k_countries)
      plot_ly(s, x = ~meta_z, y = ~label, type = "bar", orientation = "h",
              marker = list(color = ifelse(s$sign_agree >= s$k_countries & s$k_countries >= 3, ifelse(s$meta_z > 0, "#b2182b", "#2166ac"), "#bdbdbd")),
              text = ~sprintf("%s<br>pooled z %+.2f<br>%s", label, meta_z, rep), hoverinfo = "text") |>
        layout(xaxis = list(title = "Pooled association (right = more deficiency; grey = signs disagree)"), yaxis = list(title = ""),
               margin = list(l = 10, r = 10, t = 10, b = 40)) |> config(displayModeBar = FALSE)
    })

    output$mini_map <- renderLeaflet({
      col <- selected_col(); req(col, input$map_country)
      xr <- CAT$xr[[input$map_country]]; bnd <- admin2_bnds[[input$map_country]]
      validate(need(!is.null(xr) && col %in% colnames(xr), "This predictor has no usable values in this country (below the 70 percent coverage floor)."))
      bnd <- bnd[!is_water(bnd$Admin2), ]
      vals <- xr[match(.key(bnd$Admin1, bnd$Admin2), rownames(xr)), col]
      pal <- colorNumeric("YlGnBu", domain = range(vals, na.rm = TRUE), na.color = "#d9d9d9")
      leaflet(bnd) |> addProviderTiles(providers$CartoDB.Positron) |>
        addPolygons(fillColor = pal(vals), fillOpacity = 0.8, color = "#777", weight = 0.4,
                    label = sprintf("%s (%s): %s", bnd$Admin2, bnd$Admin1, fmt_num(vals, 2))) |>
        addLegend(pal = pal, values = vals[is.finite(vals)], title = "Within-country rank (normal scale)", position = "bottomright")
    })

    output$download <- downloadHandler(
      filename = function() sprintf("predictor_catalogue_%s.csv", Sys.Date()),
      content = function(file) write.csv(filtered(), file, row.names = FALSE))
  })
}
