# =============================================================================
# Module: Map explorer
# =============================================================================
# One estimator, several views: the priority ranking, how consistently each
# surveyed district is placed in the worst fifth, how far each rank moves when
# the model is re-estimated (optionally cross-hatched), the estimated
# prevalence with its checked range and WHO-threshold chance, and the survey's
# own estimate. A click gives the district's numbers and the data layers that
# push it up or down the list.

HATCH_SHARE <- 0.5   # cross-hatch when the rank range spans more than half the list

mod_map_explorer_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 340, title = "Map controls",
      selectInput(ns("country"), "Country", choices = country_choices, selected = "ghana"),
      selectInput(ns("outcome"), "Outcome", choices = outcome_choices, selected = "child_vitA"),
      uiOutput(ns("outcome_caveat")),
      radioButtons(ns("admin_level"), "Level", choices = c("District" = "admin2", "Region" = "admin1"),
                   selected = "admin2", inline = TRUE),
      selectInput(ns("layer"), "What to show",
                  choices = list(
                    "The ranking" = c("Priority score (100 = ranked worst)" = "priority",
                                      "How often placed in the worst fifth (surveyed districts)" = "p_worst_fifth",
                                      "Rank range when re-estimated" = "rank_width"),
                    "Percentages" = c("Estimated prevalence" = "prev_anchored",
                                      "Chance of being at or above the WHO 'moderate' level" = "p_modplus_cal",
                                      "Survey estimate (surveyed districts only)" = "survey_prev",
                                      "WHO severity class (estimated prevalence)" = "who_class",
                                      "People affected (estimated prevalence x population)" = "people_affected")),
                  selected = "priority"),
      checkboxInput(ns("hatch_unstable"), "Cross-hatch districts the model cannot rank firmly", value = FALSE),
      hr(),
      uiOutput(ns("headline")),
      hr(),
      uiOutput(ns("district_detail")),
      hr(),
      downloadButton(ns("brief"), "One-page country brief", class = "btn-sm btn-outline-secondary")
    ),
    card(
      full_screen = TRUE,
      card_header("District ranking of micronutrient deficiency",
                  downloadButton(ns("download"), "Download CSV", class = "btn-sm btn-outline-primary float-end")),
      card_body(padding = 0, leafletOutput(ns("map"), height = "650px")),
      card_footer(textOutput(ns("caption")),
                  tags$small(class = "text-muted",
                             "The model is built on this country's surveyed districts and then applied to every district,",
                             " so for surveyed districts the map reflects their own survey data. The accuracy figures under",
                             " How well it works come from tests in which each district was hidden from the model."))
    )
  )
}

mod_map_explorer_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$country, {
      ch <- outcomes_for(input$country)
      sel <- if ((input$outcome %||% "") %in% ch) input$outcome else ch[[1]]
      updateSelectInput(session, "outcome", choices = ch, selected = sel)
    })

    map_data <- reactive({
      req(input$country, input$outcome)
      if (!(input$outcome %in% outcomes_for(input$country))) return(NULL)
      if (input$admin_level == "admin1") get_country_admin1(input$country, input$outcome)
      else get_country_admin2(input$country, input$outcome)
    })
    area_col <- reactive(if (input$admin_level == "admin1") "Admin1" else "Admin2")

    output$outcome_caveat <- renderUI({
      cv <- biomarker_caveats[[input$outcome %||% ""]]
      if (is.null(cv)) return(NULL)
      div(style = "font-size:0.78em; color:#8a6d3b; background:#fcf8e3; border-left:3px solid #e0c97f; padding:5px 8px; border-radius:3px; margin:-2px 0 6px;",
          bsicons::bs_icon("info-circle"), " ", cv)
    })

    output$headline <- renderUI({
      req(input$country, input$outcome)
      nat <- idx_national[idx_national$country_key == input$country & idx_national$outcome == input$outcome, ]
      if (!nrow(nat)) return(empty_state("This country's survey did not measure this outcome."))
      d <- idx_districts[idx_districts$country_key == input$country & idx_districts$outcome == input$outcome, ]
      k <- ceiling(nat$n_districts / 5)
      worst <- d[order(d$rank_worst), ][seq_len(k), ]
      firm <- sum(d$p_worst_fifth >= 0.8, na.rm = TRUE); scored <- sum(is.finite(d$p_worst_fifth))
      tagList(
        h5(meta$countries[[input$country]], style = "margin-top:0;"),
        p(em(meta$outcome_labels[[input$outcome]])),
        p("National prevalence (survey): ", tags$span(style = "font-size:1.2em; color:#C8641E;", strong(fmt_pct(nat$national_prev)))),
        p(sprintf("%d of %d districts were surveyed. The model ranks all %d.", nat$n_surveyed, nat$n_districts, nat$n_districts),
          style = "font-size:0.9em; color:#555;"),
        p(level_skill_badge(nat$rho_train), style = "margin-bottom:2px;"),
        p(level_skill_text(nat$rho_train), style = "font-size:0.82em; color:#555;"),
        p(strong(sprintf("Worst fifth (%d districts): ", k)),
          paste(head(worst$Admin2, 8), collapse = ", "), if (k > 8) sprintf(" and %d more", k - 8) else "",
          style = "font-size:0.9em;"),
        if (scored > 0) p(sprintf("Of %d surveyed districts, %d were placed in the worst fifth in at least 80%% of 40 test runs.", scored, firm),
                          style = "font-size:0.85em; color:#555;")
        else p("This check needs at least 15 surveyed districts, so it is not available here.", style = "font-size:0.85em; color:#555;")
      )
    })

    output$map <- renderLeaflet({
      df <- map_data(); req(df, nrow(df) > 0)
      layer <- input$layer
      if (!layer %in% names(df) || !any(is.finite(suppressWarnings(as.numeric(df[[layer]])))) && layer != "who_class") layer <- "priority"
      vals <- df[[layer]]
      legend_title <- "Priority score"; pal <- NULL; fill <- rep("#d9d9d9", nrow(df)); lab_fmt <- labelFormat()
      if (layer == "who_class") {
        pal <- colorFactor(unname(who_colors), levels = names(who_colors), na.color = "#cccccc")
        fill <- pal(df$who_class); legend_title <- "WHO class (estimated prevalence)"
      } else if (layer == "priority") {
        pal <- colorNumeric("YlOrRd", domain = c(0, 100), na.color = "#d9d9d9"); fill <- pal(vals)
      } else if (layer == "p_worst_fifth") {
        pal <- colorNumeric(c("#f7fbff", "#9ecae1", "#08519c"), domain = c(0, 1), na.color = "#d9d9d9"); fill <- pal(vals)
        legend_title <- "Placed in worst fifth (share of runs)"; lab_fmt <- labelFormat(transform = function(x) round(100 * x), suffix = "%")
      } else if (layer == "rank_width") {
        nz <- vals[is.finite(vals)]
        if (length(nz)) { pal <- colorNumeric(c("#252525", "#bdbdbd", "#f7f7f7"), domain = range(nz), na.color = "#d9d9d9"); fill <- pal(vals) }
        legend_title <- "Places the rank moves (dark = firm)"
      } else if (layer == "p_modplus_cal") {
        pal <- colorNumeric(c("#f7fbff", "#fdae61", "#d7191c"), domain = c(0, 1), na.color = "#d9d9d9"); fill <- pal(vals)
        legend_title <- "Chance at or above 'moderate'"; lab_fmt <- labelFormat(transform = function(x) round(100 * x), suffix = "%")
      } else if (layer == "people_affected") {
        nz <- vals[is.finite(vals) & vals > 0]
        if (length(nz)) { pal <- colorNumeric("YlOrRd", domain = log10(range(pmax(nz, 1))), na.color = "#d9d9d9")
          fill <- ifelse(is.finite(vals) & vals > 0, pal(log10(pmax(vals, 1))), "#d9d9d9") }
        legend_title <- "People affected"; lab_fmt <- labelFormat(transform = function(x) round(10^x))
      } else {
        nz <- vals[is.finite(vals)]
        if (length(nz)) { pal <- colorNumeric("YlOrRd", domain = range(nz), na.color = "#d9d9d9"); fill <- pal(vals) }
        legend_title <- if (layer == "prev_anchored") "Estimated prevalence" else "Survey estimate"
        lab_fmt <- labelFormat(transform = function(x) round(100 * x, 1), suffix = "%")
      }
      surveyed <- isTRUE(input$admin_level == "admin2") & !is.na(df$surveyed) & df$surveyed
      if (input$admin_level == "admin1") surveyed <- df$n_surveyed > 0 & !is.na(df$n_surveyed)
      name <- df[[area_col()]]
      sub <- if (area_col() == "Admin2") ifelse(is.na(df$Admin1), "", df$Admin1) else rep("", nrow(df))
      rng <- if ("rank_lo" %in% names(df)) ifelse(is.finite(df$rank_lo), sprintf("<br/>Rank range when re-estimated: %s to %s", round(df$rank_lo), round(df$rank_hi)), "") else ""
      labels <- sprintf(paste0("<strong>%s</strong><br/>%s<br/>Priority score: %s (rank %s of %s)%s<br/>Estimated prevalence: %s<br/>",
                               "Survey estimate: %s<br/>Placed in worst fifth: %s<br/>Population: %s"),
                        name, sub, fmt_num(df$priority, 0), ifelse(is.finite(df$rank_worst), df$rank_worst, "—"),
                        if (area_col() == "Admin2") df$n_districts else nrow(df), rng,
                        fmt_pct(df$prev_anchored), ifelse(is.finite(df$survey_prev), fmt_pct(df$survey_prev), "not surveyed"),
                        ifelse(is.finite(df$p_worst_fifth), sprintf("%s of runs", fmt_pct(df$p_worst_fifth, 0)), "—"),
                        fmt_count(df$population)) |> lapply(HTML)
      m <- leaflet(df) |> addProviderTiles(providers$Esri.WorldGrayCanvas) |>
        addPolygons(fillColor = fill, fillOpacity = 0.78, color = ifelse(surveyed, "#1a1a1a", "#9aa0a6"),
                    weight = ifelse(surveyed, 1.6, 0.5), opacity = 1,
                    highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE),
                    label = labels, labelOptions = labelOptions(textsize = "13px", direction = "auto"),
                    layerId = name)
      hatch_key <- ""
      if (isTRUE(input$hatch_unstable) && input$admin_level == "admin2" && "width_share" %in% names(df)) {
        flag <- is.finite(df$width_share) & df$width_share > HATCH_SHARE
        hl <- if (any(flag)) hatch_lines(df[flag, ]) else NULL
        if (!is.null(hl)) m <- m |> addPolylines(data = hl, color = "#1a1a1a", weight = 0.9, opacity = 0.8, options = pathOptions(interactive = FALSE))
        hatch_key <- paste0("<br/><span style='display:inline-block;width:16px;height:11px;vertical-align:middle;",
                            "background:repeating-linear-gradient(45deg,#222 0 1px,transparent 1px 4px),",
                            "repeating-linear-gradient(-45deg,#222 0 1px,transparent 1px 4px);'></span> rank not firmly placed")
      }
      m <- m |> addControl(position = "bottomleft", html = HTML(paste0(
          "<div style='background:rgba(255,255,255,0.88);padding:4px 8px;border-radius:4px;font-size:11px;line-height:1.5;'>",
          "<span style='display:inline-block;width:16px;border-top:3px solid #1a1a1a;vertical-align:middle;'></span> surveyed<br/>",
          "<span style='display:inline-block;width:16px;border-top:2px solid #9aa0a6;vertical-align:middle;'></span> no survey clusters",
          hatch_key, "</div>")))
      if (layer == "who_class") m <- m |> addLegend(colors = unname(who_colors), labels = names(who_colors), opacity = 0.78, title = legend_title, position = "bottomright")
      else if (!is.null(pal)) {
        lv <- if (layer == "people_affected") log10(pmax(vals[is.finite(vals) & vals > 0], 1)) else if (layer == "priority") c(0, 100) else if (layer %in% c("p_worst_fifth", "p_modplus_cal")) c(0, 1) else vals[is.finite(vals)]
        m <- m |> addLegend(pal = pal, values = lv, opacity = 0.78, title = legend_title, position = "bottomright", labFormat = lab_fmt)
      }
      m
    })

    clicked <- reactiveVal(NULL)
    observeEvent(input$map_shape_click, clicked(input$map_shape_click$id))
    observeEvent(list(input$admin_level, input$country, input$outcome), clicked(NULL), ignoreInit = TRUE)

    output$district_detail <- renderUI({
      area <- clicked(); df <- map_data()
      if (is.null(area) || is.null(df)) return(p(em("Click a district to see its figures and what drives them."), style = "color:#888;"))
      row <- df[df[[area_col()]] == area, , drop = FALSE]; if (!nrow(row)) return(NULL)
      row <- row[1, ]
      dec <- if (area_col() == "Admin2") decompose_district(input$country, input$outcome, row$Admin1, row$Admin2) else NULL
      tagList(
        h5(area, style = "margin-top:0;"),
        if (area_col() == "Admin2" && !is.na(row$Admin1)) p(em(row$Admin1)),
        p(strong("Priority score: "), fmt_num(row$priority, 0),
          if (is.finite(row$rank_worst)) sprintf(" (ranked %d of %d; 1 = worst)", row$rank_worst, if (area_col() == "Admin2") row$n_districts else nrow(df)) else ""),
        if (is.finite(row$p_worst_fifth)) p(strong("Placed in the worst fifth: "), sprintf("%s of 40 test runs", fmt_pct(row$p_worst_fifth, 0)),
                                             tags$small(" (the district was hidden in each run)")),
        if (is.finite(row$rank_lo)) p(strong("Rank range when re-estimated: "),
                                      sprintf("%d to %d", round(row$rank_lo), round(row$rank_hi)),
                                      tags$small(" (how firmly the model places it; not a confidence interval)")),
        p(strong("Estimated prevalence: "), fmt_pct(row$prev_anchored),
          if (is.finite(row$prev_cal_lo)) tags$small(sprintf(" (checked 90%% range %s to %s: where a survey of this district would be expected to fall 9 times in 10)",
                                                             fmt_pct(row$prev_cal_lo), fmt_pct(row$prev_cal_hi)))
          else tags$small(" (the ranking converted using the national survey figure)")),
        if (is.finite(row$p_modplus_cal)) p(strong("Chance of being at or above the WHO 'moderate' level: "), fmt_pct(row$p_modplus_cal, 0),
                                            if (is.finite(row$p_severe_cal) && row$p_severe_cal > 0.005) tags$small(sprintf(" ('severe' level: %s)", fmt_pct(row$p_severe_cal, 0)))),
        if (isTRUE(row$surveyed) || (area_col() == "Admin1" && isTRUE(row$n_surveyed > 0)))
          p(strong("Survey estimate: "), fmt_pct(row$survey_prev),
            if (area_col() == "Admin2" && is.finite(row$survey_lo)) sprintf(" (95%% range %s to %s; %s respondents in %s clusters)",
                                                                             fmt_pct(row$survey_lo, 0), fmt_pct(row$survey_hi, 0), fmt_count(row$n_resp), row$n_clusters) else "")
        else p(em("The survey had no clusters here, so the model is the only estimate."), style = "color:#888;"),
        p(strong("WHO class (estimated prevalence): "), row$who_class),
        if (is.finite(row$population)) p(strong("Population in this group: "), fmt_count(row$population),
                                         if (is.finite(row$people_affected)) sprintf(", about %s affected", fmt_count(row$people_affected)) else ""),
        if (!is.null(dec) && nrow(dec)) {
          top <- head(dec, 6)
          tagList(
            h6("What moves this district up or down the list", style = "margin-top:10px;"),
            tags$table(class = "table table-sm", style = "font-size:0.8em;",
                       tags$tbody(lapply(seq_len(nrow(top)), function(i) tags$tr(
                         tags$td(style = if (top$contribution[i] > 0) "color:#b2182b;" else "color:#2166ac;",
                                 if (top$contribution[i] > 0) "▲" else "▼"),
                         tags$td(top$label[i]), tags$td(style = "color:#777;", top$source[i]))))),
            tags$small(style = "color:#777;", "▲ pushes the district toward more deficiency, ▼ toward less. These are parts of the model's score, not causes of deficiency.")
          )
        }
      )
    })

    output$caption <- renderText({
      nat <- idx_national[idx_national$country_key == input$country & idx_national$outcome == input$outcome, ]
      skill <- if (nrow(nat) && input$layer %in% c("prev_anchored", "who_class", "people_affected")) {
        b <- level_skill(nat$rho_train)$band
        sprintf(" Reliability of the district percentages here: %s (%s).%s", skill_word[[b]], fmt_num(g1(nat$rho_train), 2),
                if (b %in% c("none", "weak")) " Use the ranking view instead." else "")
      } else ""
      hatch <- if (isTRUE(input$hatch_unstable) && input$admin_level == "admin2")
        " Cross-hatched districts: the rank range when re-estimated spans more than half of the country's list." else ""
      sprintf("%s, survey year %s, %s. %s Grey areas have no prediction.%s%s",
              meta$countries[[input$country]], meta$survey_years[[input$country]],
              if (input$admin_level == "admin1") "regions (population-weighted from districts)" else "districts",
              switch(input$layer,
                     priority = "Colour shows the priority score: 100 = ranked worst in the country.",
                     p_worst_fifth = paste("Colour shows how often the district was placed in the worst fifth when the model was",
                                           "re-estimated 40 times with that district hidden (surveyed districts only). Districts placed",
                                           sprintf("there in at least 80%% of runs were in the survey's worst fifth %s of the time, against 20%% by chance.", pc(Q$cal_top))),
                     rank_width = paste("Colour shows how many places the district's rank moves when the model is re-estimated;",
                                        "dark districts are firmly placed.", paste0(stability_note(), ".")),
                     prev_anchored = sprintf("Colour shows the estimated prevalence: the model's ranking converted into a percentage using the national survey figure of %s.", fmt_pct(g1(nat$national_prev))),
                     p_modplus_cal = paste("Colour shows the chance that the district's prevalence is at or above the WHO 'moderate",
                                           "public health problem' level.", calibrated_note()),
                     survey_prev = "Colour shows the survey's own district estimate. Blank where the survey had no clusters.",
                     who_class = "Colour shows the WHO severity class of the estimated prevalence.",
                     people_affected = "Colour shows the estimated prevalence multiplied by the population of the group (log scale).", ""),
              skill, hatch)
    })

    output$download <- downloadHandler(
      filename = function() sprintf("ranking_%s_%s_%s_%s.csv", input$country, input$outcome, input$admin_level, Sys.Date()),
      content = function(file) { df <- map_data(); write.csv(sf::st_drop_geometry(df), file, row.names = FALSE) })

    output$brief <- downloadHandler(
      filename = function() sprintf("brief_%s_%s.html", input$country, Sys.Date()),
      content = function(file) {
        src <- file.path("briefs", sprintf("brief_%s.html", input$country))
        validate(need(file.exists(src), "Brief not built: run dashboard/data-raw/07_build_country_briefs.R"))
        file.copy(src, file, overwrite = TRUE)
      })
  })
}
