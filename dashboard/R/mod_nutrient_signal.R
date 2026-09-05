# =============================================================================
# Module: What tracks which nutrient
# =============================================================================
# The question a nutrition programme actually brings, and the one the dashboard
# could not previously answer. Every other tab is about WHERE deficiency is
# predicted or HOW the model works; this one is about WHICH CONDITIONS TRAVEL
# WITH WHICH NUTRIENT.
#
# It leads with a plain-language table because that is what the audience reads,
# and puts the indicator-level evidence behind it for anyone who wants it.
#
# EVIDENCE STANDARD. The headline number is CROSS-COUNTRY SIGN AGREEMENT, not a
# p-value. Under correction for having examined 451 indicators at 14-87
# districts per country, few clear conventional thresholds - the per-cell screen
# has median power 0.005 - while the same indicator with the same sign in four
# independent national surveys is far harder to explain away. The family-wise p
# is carried as a secondary column so nobody thinks it is being hidden.

mod_nutrient_signal_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(

    nav_panel(
      title = "By nutrient",
      icon = bsicons::bs_icon("clipboard2-pulse"),
      div(
        p(class = "lead",
          "Different nutrients have different geographies. These are the ",
          "conditions that travel with each one across all four surveyed ",
          "countries."),
        p("Read the last column first. It says in how many of the four ",
          "countries the indicator pointed the same way. Four out of four means ",
          "the pattern held in countries with different diets, different ",
          "farming and different survey teams, which is a stronger claim than ",
          "any single significance test at these sample sizes."),
        reactable::reactableOutput(ns("headline")),
        tags$br(),
        div(class = "alert alert-warning border",
          tags$strong("How to read these — important. "),
          "These are district-level associations, and several run ",
          tags$strong("opposite"), " to what individual-level nutrition would ",
          "predict. Districts where more households eat legumes or keep cattle ",
          "have ", tags$em("more"), " deficiency, not less. That is because ",
          "legume-eating and cattle-keeping mark rural subsistence districts, ",
          "which are poorer and more deficient — it is not evidence that ",
          "legumes or cattle harm anyone. Treating these as dietary effects ",
          "would be an ecological fallacy. Use them to decide ",
          tags$em("where"), " to look, never ", tags$em("what"),
          " to change.")
      )
    ),

    nav_panel(
      title = "Indicator families",
      icon = bsicons::bs_icon("grid-3x3-gap"),
      div(
        p("Each of the 19 indicator families, scored against each nutrient. ",
          "Positive strength means the family tracks ", tags$em("more"),
          " deficiency."),
        selectInput(ns("dom_scale"), "Outcome measured as",
                    choices = c("Deficiency prevalence",
                                "Biomarker concentration"), width = "320px"),
        reactable::reactableOutput(ns("domains")),
        p(class = "text-muted mt-2",
          "No single family is decisive: the signal is spread across many ",
          "indicators each carrying a little. The practical consequence is that ",
          "there is no shortcut indicator to collect instead of a survey.")
      )
    ),

    nav_panel(
      title = "Individual indicators",
      icon = bsicons::bs_icon("list-ul"),
      div(
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput(ns("ind_outcome"), "Nutrient", choices = NULL),
          selectInput(ns("ind_scale"), "Outcome measured as",
                      choices = c("Deficiency prevalence",
                                  "Biomarker concentration")),
          checkboxInput(ns("ind_agree"), "Only those agreeing in all countries",
                        value = TRUE)
        ),
        reactable::reactableOutput(ns("indicators")),
        p(class = "text-muted mt-2",
          "Sorted by strength of association. Positive strength means the ",
          "indicator tracks MORE deficiency. No malaria burden indicator ",
          "reaches the top of any list, which is a real result: it was tested ",
          "for and not found. The strongest malaria-domain signal is indoor ",
          "residual spraying coverage, which marks where control programmes ",
          "operate rather than where transmission is high.")
      )
    )
  )
}

mod_nutrient_signal_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    bundle <- reactive({
      p <- file.path("data", "nutrient_signal.rds")
      if (!file.exists(p)) return(NULL)
      readRDS(p)
    })

    observe({
      b <- bundle(); req(b)
      outs <- sort(unique(b$predictors$outcome))
      updateSelectInput(session, "ind_outcome", choices = outs,
                        selected = outs[1])
    })

    output$headline <- reactable::renderReactable({
      b <- bundle(); req(b)
      reactable::reactable(
        b$headline,
        columns = list(
          outcome   = reactable::colDef(name = "Nutrient", minWidth = 130),
          tracks    = reactable::colDef(name = "What tracks it", minWidth = 220,
                                        html = TRUE),
          direction = reactable::colDef(name = "Direction", minWidth = 170),
          agreeing  = reactable::colDef(name = "Countries agreeing", width = 130,
                                        align = "center"),
          reading   = reactable::colDef(name = "How to read it", minWidth = 260)
        ),
        striped = TRUE, highlight = TRUE, defaultPageSize = 6, wrap = TRUE)
    })

    output$domains <- reactable::renderReactable({
      b <- bundle(); req(b)
      d <- b$domains[b$domains$scale == input$dom_scale, ]
      d <- d[order(d$outcome, -abs(d$meta_z)), ]
      reactable::reactable(
        d[, c("outcome", "domain_name", "meta_z", "direction", "agree", "p_fwer")],
        columns = list(
          outcome     = reactable::colDef(name = "Nutrient", minWidth = 130),
          domain_name = reactable::colDef(name = "Indicator family", minWidth = 200),
          meta_z      = reactable::colDef(name = "Strength", width = 90,
                          format = reactable::colFormat(digits = 2)),
          direction   = reactable::colDef(name = "Direction", minWidth = 130),
          agree       = reactable::colDef(name = "Countries agreeing", width = 130,
                                          align = "center"),
          p_fwer      = reactable::colDef(name = "p (corrected)", width = 110,
                          format = reactable::colFormat(digits = 3))
        ),
        groupBy = "outcome", striped = TRUE, highlight = TRUE,
        defaultPageSize = 12)
    })

    output$indicators <- reactable::renderReactable({
      b <- bundle(); req(b)
      d <- b$predictors
      d <- d[d$outcome == input$ind_outcome & d$scale == input$ind_scale, ]
      if (isTRUE(input$ind_agree)) d <- d[d$agree >= d$countries, ]
      d <- d[order(-abs(d$meta_z)), ]
      d <- utils::head(d, 40)
      reactable::reactable(
        d[, c("indicator", "domain", "meta_z", "direction", "agree", "p_fwer")],
        columns = list(
          indicator = reactable::colDef(name = "Indicator", minWidth = 200),
          domain    = reactable::colDef(name = "Family", minWidth = 180),
          meta_z    = reactable::colDef(name = "Strength", width = 90,
                        format = reactable::colFormat(digits = 2)),
          direction = reactable::colDef(name = "Direction", minWidth = 130),
          agree     = reactable::colDef(name = "Countries agreeing", width = 130,
                                        align = "center"),
          p_fwer    = reactable::colDef(name = "p (corrected)", width = 110,
                        format = reactable::colFormat(digits = 3))
        ),
        striped = TRUE, highlight = TRUE, defaultPageSize = 15)
    })
  })
}
