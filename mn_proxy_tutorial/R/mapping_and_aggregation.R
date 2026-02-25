# =============================================================================
# R/mapping_and_aggregation.R
#
# Functions for:
#   - Predicting on a spatial grid.
#   - Aggregating individual predictions to Admin1 prevalence.
#   - Four prevalence estimation approaches (plug-in, cluster binomial,
#     calibrated, continuous threshold).
#   - Comparison table across approaches.
# =============================================================================

# ── Grid prediction ───────────────────────────────────────────────────────────

#' Predict biomarker / deficiency probability on a covariate grid.
#'
#' @param sl_fit     A fitted SuperLearner object.
#' @param grid_df    Data frame with the same predictor columns used in training.
#' @param X_cols     Predictor column names (must match those used in fitting).
#' @param type       "continuous" or "binary" (determines output interpretation).
#' @return The input grid_df with columns `pred_Y` (continuous) and / or
#'         `pred_prob` (binary probability) appended.
predict_grid <- function(sl_fit, grid_df, X_cols) {
  # Align columns to those used in fitting
  fit_cols <- colnames(sl_fit$X)
  available <- intersect(fit_cols, colnames(grid_df))
  missing_cols <- setdiff(fit_cols, colnames(grid_df))

  newX <- grid_df[, available, drop = FALSE]
  for (col in missing_cols) newX[[col]] <- 0   # fill missing with 0

  newX <- newX[, fit_cols, drop = FALSE]        # exact column order

  preds <- predict(sl_fit, newdata = newX)$pred
  grid_df$sl_pred <- as.numeric(preds)
  grid_df
}


# ── Aggregate to Admin1 ───────────────────────────────────────────────────────

#' Aggregate individual-level predictions to Admin1 prevalence.
#'
#' @param pred_df    Data frame with at minimum: pred_col and admin1_col.
#' @param admin1_col Column name for Admin1 identifier.
#' @param pred_col   Column name of predictions (probabilities or binary).
#' @param weight_col Optional column name for survey weights.
#' @return Data frame: admin1, n, pred_prevalence.
aggregate_to_admin1 <- function(pred_df,
                                 admin1_col = "admin1_name",
                                 pred_col   = "sl_pred",
                                 weight_col = NULL) {
  if (!is.null(weight_col) && weight_col %in% colnames(pred_df)) {
    pred_df %>%
      dplyr::group_by(.data[[admin1_col]]) %>%
      dplyr::summarise(
        n               = dplyr::n(),
        pred_prevalence = stats::weighted.mean(.data[[pred_col]],
                                               .data[[weight_col]], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    pred_df %>%
      dplyr::group_by(.data[[admin1_col]]) %>%
      dplyr::summarise(
        n               = dplyr::n(),
        pred_prevalence = mean(.data[[pred_col]], na.rm = TRUE),
        .groups         = "drop"
      )
  }
}


# ── Four prevalence estimation approaches ────────────────────────────────────

# Approach 1 ──────────────────────────────────────────────────────────────────

#' Approach 1: Plug-in mean of predicted probabilities.
#'
#' p̂_a^plug = (1/N_a) Σ_{i: A_i = a} P̂(D_i = 1 | X_i)
#'
#' @param pred_df    Data frame with admin1_col and a predicted probability column.
#' @param admin1_col Admin1 column name.
#' @param pred_col   Column of predicted probabilities.
#' @return Named numeric vector: Admin1 → prevalence estimate.
prevalence_plugin <- function(pred_df,
                               admin1_col = "admin1_name",
                               pred_col   = "pred_bin") {
  tbl <- pred_df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(prev = mean(.data[[pred_col]], na.rm = TRUE), .groups = "drop")
  setNames(tbl$prev, tbl[[admin1_col]])
}


# Approach 2 ──────────────────────────────────────────────────────────────────

#' Approach 2: Aggregate to cluster level, fit binomial GLM, predict to Admin1.
#'
#' Step 1: Compute y_c = Σ D_i, n_c for each cluster.
#' Step 2: Fit logistic regression y_c ~ X_c_mean (cluster-mean covariates).
#' Step 3: Predict p_c for each cluster, average within Admin1.
#'
#' @param df          Individual-level data frame.
#' @param admin1_col  Admin1 column name.
#' @param cluster_col Cluster column name.
#' @param outcome     Binary outcome column name.
#' @param X_cols      Predictor columns.
#' @return Named numeric vector: Admin1 → prevalence estimate.
prevalence_cluster_binomial <- function(df,
                                         admin1_col  = "admin1_name",
                                         cluster_col = "cluster_id",
                                         outcome     = "D",
                                         X_cols      = NULL) {
  if (is.null(X_cols)) X_cols <- get_predictor_cols(df)
  X_cols <- intersect(X_cols, colnames(df))
  # Keep only numeric cluster-level or near-constant variables
  X_num <- X_cols[sapply(df[, X_cols, drop = FALSE], is.numeric)]

  # Aggregate to cluster level
  clust_df <- df %>%
    dplyr::group_by(.data[[cluster_col]], .data[[admin1_col]]) %>%
    dplyr::summarise(
      y_c = sum(.data[[outcome]], na.rm = TRUE),
      n_c = dplyr::n(),
      dplyr::across(dplyr::all_of(X_num), mean, na.rm = TRUE),
      .groups = "drop"
    )

  X_clust <- clust_df[, X_num, drop = FALSE]

  # Drop near-zero variance columns
  nzv <- which(sapply(X_clust, function(x) sd(x, na.rm = TRUE) < 1e-8))
  if (length(nzv) > 0) X_clust <- X_clust[, -nzv, drop = FALSE]

  formula_str <- paste("cbind(y_c, n_c - y_c) ~",
                        paste(names(X_clust), collapse = " + "))

  fit <- tryCatch(
    glm(as.formula(formula_str),
         data   = cbind(clust_df[, c("y_c", "n_c", admin1_col)], X_clust),
         family = binomial()),
    error = function(e) {
      glm(cbind(y_c, n_c - y_c) ~ 1, data = clust_df, family = binomial())
    }
  )

  clust_df$pred_p <- predict(fit, type = "response")

  tbl <- clust_df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(prev = mean(pred_p), .groups = "drop")
  setNames(tbl$prev, tbl[[admin1_col]])
}


# Approach 3 ──────────────────────────────────────────────────────────────────

#' Approach 3: Calibration-corrected Admin1 prevalence.
#'
#' Fits logit(P(D=1)) = β₀ + β₁ · logit(p̂) on held-out CV predictions,
#' recalibrates all predictions, then re-aggregates to Admin1.
#'
#' @param df          Data frame with cv_pred and outcome columns.
#' @param admin1_col  Admin1 column name.
#' @param obs_col     Observed binary outcome column.
#' @param cv_pred_col Cross-validated predicted probability column.
#' @return Named numeric vector: Admin1 → calibrated prevalence estimate.
prevalence_calibrated <- function(df,
                                   admin1_col  = "admin1_name",
                                   obs_col     = "D",
                                   cv_pred_col = "pred_bin") {
  lp   <- qlogis(pmax(1e-6, pmin(1 - 1e-6, df[[cv_pred_col]])))
  m    <- glm(df[[obs_col]] ~ lp, family = binomial())
  p_cal <- plogis(coef(m)[1] + coef(m)[2] * lp)

  df$pred_cal <- p_cal
  tbl <- df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(prev = mean(pred_cal, na.rm = TRUE), .groups = "drop")
  setNames(tbl$prev, tbl[[admin1_col]])
}


# Approach 4 ──────────────────────────────────────────────────────────────────

#' Approach 4: Threshold continuous biomarker prediction → Admin1 prevalence.
#'
#' Converts predicted continuous values to binary using the clinical cutoff,
#' then averages within Admin1.
#'
#' NOTE: This approach can be biased when residual variance is ignored.
#' A probabilistic version accounting for residual SD is also provided.
#'
#' @param df           Data frame.
#' @param admin1_col   Admin1 column name.
#' @param cont_pred_col Column of predicted continuous biomarker.
#' @param cutoff       Classification threshold (default -0.5).
#' @param direction    "less" → D = 1 if Y < cutoff.
#' @param residual_sd  Optional: if supplied, use Gaussian tail probability
#'                     P(Y < cutoff | Ŷ, σ) instead of hard threshold.
#' @return Named numeric vector: Admin1 → prevalence estimate.
prevalence_threshold_continuous <- function(df,
                                             admin1_col    = "admin1_name",
                                             cont_pred_col = "pred_cont",
                                             cutoff        = -0.5,
                                             direction     = "less",
                                             residual_sd   = NULL) {
  yhat <- df[[cont_pred_col]]

  if (!is.null(residual_sd)) {
    # Soft thresholding: P(Y < cutoff | Ŷ, σ̂)
    prob <- if (direction == "less") pnorm((cutoff - yhat) / residual_sd)
            else                     pnorm((yhat - cutoff) / residual_sd)
  } else {
    prob <- apply_threshold(yhat, cutoff, direction)
  }

  df$pred_thresh <- prob
  tbl <- df %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(prev = mean(pred_thresh, na.rm = TRUE), .groups = "drop")
  setNames(tbl$prev, tbl[[admin1_col]])
}


# ── Comparison ────────────────────────────────────────────────────────────────

#' Compare multiple prevalence estimation approaches against true values.
#'
#' @param true_prev  Named numeric vector of true Admin1 prevalences.
#' @param est_list   Named list of named numeric vectors (one per approach).
#' @return Data frame: approach × admin1 × true × est × bias × abs_error.
compare_prevalence_approaches <- function(true_prev, est_list) {
  do.call(rbind, lapply(names(est_list), function(approach) {
    ev <- eval_prevalence(true_prev, est_list[[approach]])
    ev$approach <- approach
    ev
  })) %>%
    dplyr::arrange(approach, admin1)
}

#' Summary table (MAE, RMSE, mean bias) per approach.
#' @param comparison_df Output of \code{compare_prevalence_approaches}.
#' @return Data frame: one row per approach with scalar error metrics.
approach_error_summary <- function(comparison_df) {
  comparison_df %>%
    dplyr::group_by(approach) %>%
    dplyr::summarise(
      MAE       = round(mean(abs_error) * 100, 2),
      RMSE      = round(sqrt(mean(bias^2)) * 100, 2),
      mean_bias = round(mean(bias) * 100, 2),
      .groups   = "drop"
    )
}


# ── Plotting helpers ──────────────────────────────────────────────────────────

#' Dot plot comparing true vs estimated Admin1 prevalence.
#' @param comparison_df Output of \code{compare_prevalence_approaches}.
#' @return ggplot object.
plot_prevalence_comparison <- function(comparison_df) {
  ggplot2::ggplot(comparison_df,
                   ggplot2::aes(x = true * 100, y = est * 100,
                                colour = approach)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                          colour = "grey50") +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::facet_wrap(~ approach, ncol = 2) +
    ggplot2::labs(x     = "True Admin1 prevalence (%)",
                   y     = "Estimated prevalence (%)",
                   title = "Prevalence estimation approaches: true vs estimated") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}

#' Choropleth map of Admin1 predicted prevalence.
#'
#' @param admin1_sf    sf object with Admin1 polygons; must have a column
#'                     matching admin1_col.
#' @param pred_tbl     Data frame with admin1_col and pred_col.
#' @param admin1_col   Join key.
#' @param pred_col     Column to map.
#' @param title        Plot title.
#' @return ggplot object.
plot_admin1_map <- function(admin1_sf, pred_tbl, admin1_col = "admin1_name",
                             pred_col = "pred_prevalence", title = "") {
  map_df <- admin1_sf %>%
    dplyr::left_join(pred_tbl, by = admin1_col)

  ggplot2::ggplot(map_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[pred_col]] * 100),
                     colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(
      name     = "Prevalence (%)",
      limits   = c(0, 100),
      na.value = "grey80"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
