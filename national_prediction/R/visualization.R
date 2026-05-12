# =============================================================================
# National Prediction Pipeline — Visualization
# =============================================================================

#' VMNIS data availability heatmap
#'
#' Dot-plot heatmap showing micronutrient data coverage across LMIC countries,
#' color-coded by scaled mean value, sized by number of surveys.
#' Ported from src/vmnis_heatmaps.R.
#'
#' @param heatmap_data data.frame from clean_vmnis_for_heatmap()
#' @param population Character: population filter (e.g., "Pregnant women")
#' @param exclude_haemoglobin Logical: drop haemoglobin from plot?
#' @return ggplot object
plot_vmnis_heatmap <- function(heatmap_data,
                               population = "Pregnant women",
                               exclude_haemoglobin = FALSE) {
  plotdf <- heatmap_data %>%
    dplyr::filter(Population == population)

  if (exclude_haemoglobin) {
    plotdf <- plotdf %>%
      dplyr::filter(mn_group != "Haemoglobin") %>%
      droplevels()
  }

  ggplot2::ggplot(plotdf, ggplot2::aes(y = mn_group, x = Country, fill = mean)) +
    ggplot2::geom_point(ggplot2::aes(size = n_surveys), shape = 21, colour = "black") +
    ggplot2::scale_size_continuous(range = c(1, 6), name = "Number of surveys") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(
      low = "navy", high = "skyblue",
      na.value = "grey90", name = "Scaled mean\n(latest survey)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title    = paste0("VMNIS data availability: ", population),
      subtitle = "Most recent national survey\ncovering both urban & rural areas",
      x = "Country", y = "Micronutrient"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}


#' VMNIS heatmap faceted by all populations
#'
#' @param heatmap_data data.frame from clean_vmnis_for_heatmap()
#' @return ggplot object
plot_vmnis_heatmap_all <- function(heatmap_data) {
  ggplot2::ggplot(heatmap_data,
                  ggplot2::aes(y = mn_group, x = Country, fill = mean)) +
    ggplot2::geom_point(ggplot2::aes(size = n_surveys), shape = 21, colour = "black") +
    ggplot2::scale_size_continuous(range = c(1, 6), name = "Number of surveys") +
    ggplot2::facet_wrap(~Population) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(
      low = "navy", high = "skyblue",
      na.value = "grey90", name = "Scaled mean\n(latest survey)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(
      title    = "VMNIS data availability by population",
      subtitle = "Most recent national survey covering both urban & rural areas",
      x = "Country", y = "Micronutrient"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}


#' Performance comparison dot plot
#'
#' Shows key metrics (R-squared, correlation, RMSE) across all outcomes,
#' colored by data source (BRINDA vs VMNIS).
#'
#' @param perf_df data.frame: combined performance from both sources,
#'   must have columns: outcome, source, R2_vs_mean, correlation, rmse
#' @return ggplot object
plot_performance_comparison <- function(perf_df) {
  # Pivot to long format for faceted display
  plot_data <- perf_df %>%
    dplyr::select(outcome, source, R2_vs_mean, correlation, rmse) %>%
    tidyr::pivot_longer(
      cols      = c(R2_vs_mean, correlation, rmse),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric == "R2_vs_mean"  ~ "R-squared (vs. mean)",
        metric == "correlation" ~ "Correlation",
        metric == "rmse"        ~ "RMSE"
      )
    )

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = value, y = reorder(outcome, value),
                               color = source)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c("BRINDA" = "#E64B35", "VMNIS" = "#4DBBD5"),
      name   = "Data source"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Model performance across outcomes",
      x     = "Metric value",
      y     = NULL
    )
}


#' Observed vs predicted scatter plot
#'
#' @param fit_results List of fit results from run_national_sl()
#' @param outcome_labels Named character vector mapping outcome names to labels
#' @return ggplot object
plot_obs_vs_pred <- function(fit_results, outcome_labels = NULL) {
  # Combine res data.frames from all fit results
  all_res <- purrr::map2_dfr(fit_results, names(fit_results), function(res, nm) {
    if (inherits(res, "try-error") || is.null(res)) return(NULL)
    res$res %>%
      dplyr::mutate(outcome_name = nm)
  })

  if (nrow(all_res) == 0) {
    warning("No valid results to plot")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Apply labels if provided
  if (!is.null(outcome_labels)) {
    all_res$outcome_label <- outcome_labels[all_res$outcome_name]
  } else {
    all_res$outcome_label <- all_res$outcome_name
  }

  ggplot2::ggplot(all_res, ggplot2::aes(x = Y, y = yhat_full)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~outcome_label, scales = "free") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Observed vs. predicted national prevalence",
      x     = "Observed",
      y     = "Predicted"
    )
}


#' Variable importance bar chart
#'
#' @param varimp_top data.frame: top N variables from calc_importance()
#' @param outcome_label Character: outcome label for plot title
#' @return ggplot object
plot_varimp_bar <- function(varimp_top, outcome_label = "") {
  var_col   <- colnames(varimp_top)[1]
  value_col <- colnames(varimp_top)[2]

  ggplot2::ggplot(
    varimp_top,
    ggplot2::aes(
      x = stats::reorder(.data[[var_col]], .data[[value_col]]),
      y = .data[[value_col]]
    )
  ) +
    ggplot2::geom_col(fill = "#4DBBD5") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = paste0("Top predictors: ", outcome_label),
      x     = NULL,
      y     = "Importance (permutation ratio)"
    )
}
