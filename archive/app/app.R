# =============================================================================
# Shiny App: MN-Prediction Results Explorer
#
# Loads targets pipeline results and displays:
#   1. 2x2 map grid (observed, predicted, absolute error, relative error)
#   2. Performance metrics table
#   3. Plain-language metric explanations
# =============================================================================

library(shiny)
library(bslib)
library(targets)
library(sf)
library(ggplot2)
library(patchwork)
library(dplyr)
library(geodata)
library(DT)
library(scales)

# ── Globals -------------------------------------------------------------------

STORE <- here::here("_targets")

COUNTRIES <- c(
  "Gambia"       = "gambia",
  "Ghana"        = "ghana",
  "Sierra Leone" = "sierraleone",
  "Malawi"       = "malawi"
)

GADM_CODES <- c(
  gambia      = "GMB",
  ghana       = "GHA",
  sierraleone = "SLE",
  malawi      = "MWI"
)

OUTCOMES <- c(
  "Child Vitamin A Deficiency" = "child_vitA",
  "Women Vitamin A Deficiency" = "women_vitA",
  "Child Iron Deficiency"      = "child_iron",
  "Women Iron Deficiency"      = "women_iron"
)

# Cache shapefiles in-session
shp_cache <- new.env(parent = emptyenv())

get_admin1_shp <- function(country_key) {
  if (!is.null(shp_cache[[country_key]])) return(shp_cache[[country_key]])
  gadm_code <- GADM_CODES[[country_key]]
  shp <- tryCatch({
    raw <- geodata::gadm(gadm_code, level = 1,
                         path = here::here("data", "gadm"))
    sf::st_as_sf(raw)
  }, error = function(e) NULL)
  shp_cache[[country_key]] <- shp
  shp
}

# Safe target reader
safe_tar_read <- function(name) {
  tryCatch(
    targets::tar_read_raw(name, store = STORE),
    error = function(e) NULL
  )
}

# ── UI ------------------------------------------------------------------------

ui <- page_sidebar(
  title = "MN-Prediction Results Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 280,
    selectInput("country", "Country",
                choices = COUNTRIES, selected = "gambia"),
    selectInput("outcome", "Outcome",
                choices = OUTCOMES, selected = "child_vitA"),
    hr(),
    p(class = "text-muted small",
      "Select a country and outcome to view Admin-1 level maps ",
      "comparing observed survey prevalence to model predictions."),
    hr(),
    actionButton("refresh", "Reload Data", class = "btn-outline-primary btn-sm")
  ),

  layout_columns(
    col_widths = c(12),

    # ── Row 1: Map grid (2x2) ──
    card(
      card_header(textOutput("map_title")),
      card_body(
        layout_columns(
          col_widths = c(6, 6),
          plotOutput("map_obs",  height = "450px"),
          plotOutput("map_pred", height = "450px"),
          plotOutput("map_abs",  height = "450px"),
          plotOutput("map_rel",  height = "450px")
        )
      )
    ),

    # ── Row 2: Performance metrics ──
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Binary Model Diagnostics"),
        card_body(DTOutput("metrics_binary"))
      ),
      card(
        card_header("Continuous Model Diagnostics"),
        card_body(DTOutput("metrics_continuous"))
      )
    ),

    # ── Row 3: Metric explanations ──
    card(
      card_header("How to Judge Model Quality"),
      card_body(
        class = "small",
        h5("Discrimination — Can the model rank who is deficient?"),
        tags$dl(
          tags$dt("ROC-AUC"),
          tags$dd("Area under the Receiver Operating Characteristic curve. ",
                   "Ranges from 0.5 (coin flip) to 1.0 (perfect). ",
                   "Measures how well the model separates deficient from non-deficient ",
                   "individuals, regardless of the chosen threshold. ",
                   "An AUC of 0.70 means that a randomly chosen deficient person ",
                   "scores higher than a non-deficient person 70% of the time."),
          tags$dt("PR-AUC"),
          tags$dd("Area under the Precision-Recall curve. ",
                   "More informative than ROC-AUC when the outcome is rare (< 10% prevalence). ",
                   "Baseline equals the prevalence rate (e.g., 0.04 for a 4% outcome), ",
                   "so a PR-AUC of 0.15 when prevalence is 4% is actually strong. ",
                   "Values near baseline indicate no discrimination.")
        ),
        h5("Calibration — Are the predicted probabilities accurate?"),
        tags$dl(
          tags$dt("Calibration Intercept"),
          tags$dd("From a logistic recalibration model: log-odds(Y) ~ α + β × log-odds(ŷ). ",
                   "Ideal = 0. Positive values mean the model underestimates risk; ",
                   "negative values mean overestimation. Large magnitude = systematic bias."),
          tags$dt("Calibration Slope"),
          tags$dd("From the same recalibration model. Ideal = 1. ",
                   "Values < 1 indicate overfitting (predictions are too extreme); ",
                   "values > 1 suggest the model is too conservative (predictions are compressed).")
        ),
        h5("Overall Accuracy"),
        tags$dl(
          tags$dt("Brier Score"),
          tags$dd("Mean squared difference between predicted probability and actual outcome (0/1). ",
                   "Ranges from 0 (perfect) to 1 (worst). ",
                   "Decomposable into calibration and discrimination components."),
          tags$dt("Brier Skill Score"),
          tags$dd("Improvement over a naive model that always predicts the prevalence. ",
                   "Positive = better than naive; 0 = same as naive; negative = worse. ",
                   "A skill score of 0.10 means 10% improvement in squared error."),
          tags$dt("R² / Correlation (continuous)"),
          tags$dd("R² is the proportion of variance explained by the model. ",
                   "Pearson r measures linear agreement between predicted and observed. ",
                   "For continuous outcomes like hemoglobin or retinol concentration."),
          tags$dt("RMSE / MAE (continuous)"),
          tags$dd("Root Mean Squared Error and Mean Absolute Error — ",
                   "average prediction error in the original units (e.g., µg/dL). ",
                   "Lower is better. RMSE penalises large errors more heavily.")
        ),
        h5("Prevalence Recovery — Does the model get the population rate right?"),
        tags$dl(
          tags$dt("Observed vs Predicted Prevalence"),
          tags$dd("Comparing the mean predicted deficiency rate (from thresholded continuous ",
                   "predictions or binary probabilities) against the survey-measured rate. ",
                   "Close agreement means the model recovers population-level burden ",
                   "even if individual predictions are noisy."),
          tags$dt("Admin-1 Scatter"),
          tags$dd("Comparing prevalence at the subnational level. Points near the 1:1 line ",
                   "mean the model captures geographic variation. The maps above show this spatially.")
        ),
        h5("Evaluating Rare Outcomes (prevalence < 10%)"),
        p("When deficiency is rare, accuracy metrics can be misleading. A model that ",
          "predicts 'not deficient' for everyone achieves 96% accuracy when prevalence is 4%, ",
          "but is useless for identifying at-risk populations. For rare outcomes:"),
        tags$ul(
          tags$li(tags$strong("PR-AUC"), " is more informative than ROC-AUC"),
          tags$li(tags$strong("Brier Skill Score"), " shows improvement over naive baseline"),
          tags$li(tags$strong("Prevalence recovery"), " at Admin-1 level matters most for policy"),
          tags$li("Even modest AUC (0.60–0.70) with good calibration can produce useful ",
                  "subnational prevalence maps")
        )
      )
    )
  )
)

# ── Server --------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive data loader
  data <- reactiveValues()

  observeEvent(list(input$country, input$outcome, input$refresh), {
    req(input$country, input$outcome)
    ctry <- input$country
    oc   <- input$outcome
    tag  <- paste0(ctry, "_", oc)

    data$admin1     <- safe_tar_read(paste0("admin1_prev_", tag))
    data$cv_perf    <- safe_tar_read(paste0("cv_perf_", tag))
    data$diagnostics <- safe_tar_read(paste0("diagnostics_", tag))
    data$bootstrap  <- safe_tar_read(paste0("bootstrap_ci_", tag))
    data$shp        <- get_admin1_shp(ctry)
  }, ignoreNULL = FALSE)

  # ── Shared map data reactive ──
  map_prepared <- reactive({
    req(data$admin1, data$shp)

    admin1 <- data$admin1
    shp    <- data$shp

    shp$NAME_1_clean <- toupper(trimws(shp$NAME_1))
    admin1$Admin1_clean <- toupper(trimws(admin1$Admin1))

    map_data <- shp %>%
      dplyr::left_join(admin1, by = c("NAME_1_clean" = "Admin1_clean"))

    has_obs <- "obs_prev" %in% colnames(map_data) &&
               any(!is.na(map_data$obs_prev))

    if (has_obs) {
      map_data <- map_data %>%
        dplyr::mutate(
          abs_error = abs(sl_prev - obs_prev),
          rel_error = dplyr::if_else(
            obs_prev > 0.005,
            abs(sl_prev - obs_prev) / obs_prev * 100,
            NA_real_
          )
        )
    }

    prev_range <- range(
      c(map_data$sl_prev, if (has_obs) map_data$obs_prev),
      na.rm = TRUE
    )
    prev_range[1] <- max(0, prev_range[1] - 0.02)
    prev_range[2] <- min(1, prev_range[2] + 0.02)

    bbox <- sf::st_bbox(shp)

    list(
      map_data   = map_data,
      has_obs    = has_obs,
      prev_range = prev_range,
      bbox       = bbox
    )
  })

  # Shared theme and coord helper
  map_theme_fn <- function() {
    theme_minimal(base_size = 13) +
      theme(
        axis.text      = element_blank(),
        axis.ticks     = element_blank(),
        axis.title     = element_blank(),
        panel.grid     = element_blank(),
        legend.position  = "bottom",
        legend.key.width  = unit(1.8, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.text    = element_text(size = 10),
        legend.title   = element_text(size = 11),
        plot.title     = element_text(face = "bold", size = 14),
        plot.margin    = margin(5, 5, 5, 5)
      )
  }

  coord_fn <- function(bbox) {
    coord_sf(
      xlim   = c(bbox["xmin"], bbox["xmax"]),
      ylim   = c(bbox["ymin"], bbox["ymax"]),
      expand = FALSE
    )
  }

  # ── Card title ──
  output$map_title <- renderText({
    ctry_label <- names(COUNTRIES)[COUNTRIES == input$country]
    oc_label   <- names(OUTCOMES)[OUTCOMES == input$outcome]
    paste0(ctry_label, " \u2014 ", oc_label, ": Admin-1 Prevalence Maps")
  })

  # ── Panel 1: Observed ──
  output$map_obs <- renderPlot({
    mp <- map_prepared()
    if (mp$has_obs) {
      ggplot(mp$map_data) +
        geom_sf(aes(fill = obs_prev), color = "white", linewidth = 0.4) +
        scale_fill_viridis_c(
          "Prevalence", limits = mp$prev_range,
          labels = percent_format(accuracy = 1),
          option = "mako", direction = -1, na.value = "grey85"
        ) +
        labs(title = "Observed (Survey)") +
        coord_fn(mp$bbox) + map_theme_fn()
    } else {
      ggplot(mp$map_data) +
        geom_sf(fill = "grey85", color = "white", linewidth = 0.4) +
        labs(title = "Observed (Survey)", subtitle = "Not available") +
        coord_fn(mp$bbox) + map_theme_fn()
    }
  }, res = 96)

  # ── Panel 2: Predicted ──
  output$map_pred <- renderPlot({
    mp <- map_prepared()
    ggplot(mp$map_data) +
      geom_sf(aes(fill = sl_prev), color = "white", linewidth = 0.4) +
      scale_fill_viridis_c(
        "Prevalence", limits = mp$prev_range,
        labels = percent_format(accuracy = 1),
        option = "mako", direction = -1, na.value = "grey85"
      ) +
      labs(title = "Predicted (SuperLearner)") +
      coord_fn(mp$bbox) + map_theme_fn()
  }, res = 96)

  # ── Panel 3: Absolute Error ──
  output$map_abs <- renderPlot({
    mp <- map_prepared()
    if (mp$has_obs) {
      ggplot(mp$map_data) +
        geom_sf(aes(fill = abs_error), color = "white", linewidth = 0.4) +
        scale_fill_distiller(
          "Abs. Error", palette = "YlOrRd", direction = 1,
          labels = percent_format(accuracy = 0.1),
          na.value = "grey85"
        ) +
        labs(title = "Absolute Error |Pred \u2212 Obs|") +
        coord_fn(mp$bbox) + map_theme_fn()
    } else {
      ggplot(mp$map_data) +
        geom_sf(fill = "grey85", color = "white", linewidth = 0.4) +
        labs(title = "Absolute Error", subtitle = "Not available") +
        coord_fn(mp$bbox) + map_theme_fn()
    }
  }, res = 96)

  # ── Panel 4: Relative Error ──
  output$map_rel <- renderPlot({
    mp <- map_prepared()
    if (mp$has_obs) {
      ggplot(mp$map_data) +
        geom_sf(aes(fill = rel_error), color = "white", linewidth = 0.4) +
        scale_fill_distiller(
          "Rel. Error (%)", palette = "YlOrRd", direction = 1,
          na.value = "grey85"
        ) +
        labs(title = "Relative Error (%)") +
        coord_fn(mp$bbox) + map_theme_fn()
    } else {
      ggplot(mp$map_data) +
        geom_sf(fill = "grey85", color = "white", linewidth = 0.4) +
        labs(title = "Relative Error (%)", subtitle = "Not available") +
        coord_fn(mp$bbox) + map_theme_fn()
    }
  }, res = 96)

  # ── Metrics tables ──

  output$metrics_binary <- renderDT({
    diag <- data$diagnostics
    if (is.null(diag)) return(NULL)

    bin <- diag$binary_metrics
    if (is.null(bin) || nrow(bin) == 0) return(NULL)

    display <- bin %>%
      dplyr::transmute(
        Metric = "Binary Model",
        N         = n,
        Events    = n_events,
        Prevalence = sprintf("%.1f%%", prevalence * 100),
        `ROC-AUC`  = round(roc_auc, 3),
        `PR-AUC`   = round(pr_auc, 3),
        Brier      = round(brier, 4),
        `Brier Skill` = round(brier_skill, 3),
        `Calib. Int.`   = round(calib_int, 3),
        `Calib. Slope`  = round(calib_slope, 3)
      )

    datatable(
      display, rownames = FALSE,
      options = list(
        dom = "t", paging = FALSE, searching = FALSE,
        initComplete = htmlwidgets::JS(
          "function(settings, json) {",
          "  $(this.api().table().container()).css({'font-size': '11px'});",
          "}"
        )
      ),
      class = "compact stripe"
    )
  })

  output$metrics_continuous <- renderDT({
    diag <- data$diagnostics
    if (is.null(diag)) return(NULL)

    cont <- diag$continuous_metrics
    if (is.null(cont) || nrow(cont) == 0) return(NULL)

    display <- cont %>%
      dplyr::transmute(
        Metric = "Continuous Model",
        N    = n,
        RMSE = round(rmse, 4),
        MAE  = round(mae, 4),
        `R²` = round(r2, 3),
        r    = round(r, 3),
        `Prev. Obs`  = if ("prev_obs" %in% names(.))
          sprintf("%.1f%%", prev_obs * 100) else NA_character_,
        `Prev. Pred` = if ("prev_pred" %in% names(.))
          sprintf("%.1f%%", prev_pred * 100) else NA_character_
      )

    datatable(
      display, rownames = FALSE,
      options = list(
        dom = "t", paging = FALSE, searching = FALSE,
        initComplete = htmlwidgets::JS(
          "function(settings, json) {",
          "  $(this.api().table().container()).css({'font-size': '11px'});",
          "}"
        )
      ),
      class = "compact stripe"
    )
  })
}

shinyApp(ui, server)
