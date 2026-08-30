# =============================================================================
# R/plotting.R
#
# Functions for generating pipeline figures.
# Each function takes data + config and returns a ggplot (or saves a file).
# =============================================================================

#' Plot Admin1 scatter: predicted vs observed prevalence
#'
#' @param admin1_prev Data frame from aggregate_admin1()
#' @param oc Outcome config
#' @param cc Country config
#' @return ggplot object
plot_admin1_scatter <- function(admin1_prev, oc, cc) {
  if (is.null(admin1_prev) || !"obs_prev" %in% colnames(admin1_prev)) return(NULL)

  df <- admin1_prev %>% dplyr::filter(!is.na(sl_prev), !is.na(obs_prev))
  if (nrow(df) < 2) return(NULL)

  mae <- mean(abs(df$sl_prev - df$obs_prev))
  r   <- cor(df$sl_prev, df$obs_prev)

  ggplot2::ggplot(df, ggplot2::aes(x = obs_prev, y = sl_prev)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = "steelblue", alpha = 0.8) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Admin1), size = 2.5,
                              max.overlaps = 15) +
    ggplot2::annotate("text", x = min(df$obs_prev), y = max(df$sl_prev),
                       label = sprintf("MAE = %.1f pp\nr = %.2f", mae * 100, r),
                       hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
    ggplot2::scale_size_area(name = "n", max_size = 6) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x        = "Observed prevalence",
      y        = "Predicted prevalence",
      title    = sprintf("Admin1: %s", oc$label),
      subtitle = cc$country
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Admin2 combined coverage map (surveyed + unsurveyed)
#'
#' @param area_model_result Output from fit_area_level_model()
#' @param oc Outcome config
#' @param cc Country config
#' @return ggplot object
plot_admin2_coverage <- function(area_model_result, oc, cc) {
  if (is.null(area_model_result$area_preds)) return(NULL)

  polys <- area_model_result$polygons
  preds <- area_model_result$area_preds

  # Pair key where available: joining 256 Malawi polygons to 256 predictions
  # on Admin2 alone is a many-to-many join (dplyr warned), which paints
  # same-named districts with each other's prediction and inflates the frame.
  map_data <- polys %>%
    dplyr::left_join(
      preds %>% dplyr::select(dplyr::any_of(c("Admin1", "Admin2")),
                              area_pred_prev, has_survey),
      by = admin2_join_by(polys, preds)
    )

  p <- ggplot2::ggplot(map_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = area_pred_prev),
                      colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_sf(
      data = map_data %>% dplyr::filter(!has_survey),
      fill = NA, colour = "red", linewidth = 1.0
    ) +
    ggplot2::scale_fill_viridis_c(
      name   = "Predicted\nprevalence",
      labels = scales::percent_format(accuracy = 1),
      na.value = "grey85", option = "plasma"
    ) +
    ggplot2::labs(
      title    = sprintf("Area-Level Model: %s", oc$label),
      subtitle = sprintf("%s | Red borders = unsurveyed", cc$country)
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  p
}


#' Plot Admin2 forest plot with CIs
#'
#' @param area_model_result Output from fit_area_level_model()
#' @param oc Outcome config
#' @param cc Country config
#' @return ggplot object
plot_admin2_forest <- function(area_model_result, oc, cc) {
  preds <- area_model_result$area_preds
  if (is.null(preds) || !"ci_lo" %in% colnames(preds)) return(NULL)

  # Deduplicate Admin2 names (e.g. Malawi has duplicate TA names across districts)
  preds$Admin2 <- make.unique(as.character(preds$Admin2), sep = " — ")

  df <- preds %>%
    dplyr::arrange(boot_mean) %>%
    dplyr::mutate(
      Admin2 = factor(Admin2, levels = Admin2),
      source = ifelse(has_survey, "Surveyed", "Unsurveyed")
    )

  ggplot2::ggplot(df, ggplot2::aes(x = boot_mean, y = Admin2,
                                     xmin = ci_lo, xmax = ci_hi,
                                     colour = source)) +
    ggplot2::geom_errorbarh(height = 0.4, linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(has_survey),
      ggplot2::aes(x = svy_prev), shape = 4, size = 2.5,
      colour = "black", stroke = 0.8
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_colour_manual(
      name = "Admin2 type",
      values = c("Surveyed" = "steelblue", "Unsurveyed" = "firebrick")
    ) +
    ggplot2::labs(
      x     = "Predicted prevalence (bootstrap mean + 95% CI)",
      y     = NULL,
      title = sprintf("All Admin2: %s", oc$label),
      subtitle = sprintf("%s | x = survey-weighted estimate", cc$country)
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7)
    )
}
