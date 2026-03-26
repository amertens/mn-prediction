# =============================================================================
# R/admin1_analysis.R
#
# Functions for Admin1-level prediction, aggregation, and CV performance.
# =============================================================================

#' Apply a cutoff threshold to continuous predictions
#'
#' @param pred Numeric vector of continuous predictions
#' @param cutoff Threshold value
#' @param direction "less" (deficient if pred < cutoff) or "greater"
#' @return Numeric 0/1 vector
apply_threshold <- function(pred, cutoff, direction = "less") {
  if (direction == "less") as.numeric(pred < cutoff)
  else                     as.numeric(pred > cutoff)
}


#' Extract CV performance metrics from a fitted SL object
#'
#' @param sl_result Output from DHS_SL_clustered (one element of fit_sl_models)
#' @param oc Outcome config
#' @param model_type "continuous" or "binary"
#' @return data.frame with one row of metrics
extract_cv_performance <- function(sl_result, oc, model_type) {
  if (is.null(sl_result)) return(NULL)

  res <- sl_result$res
  Y    <- res$Y
  yhat <- res$yhat_full

  metrics <- data.frame(
    outcome    = oc$tag,
    model_type = model_type,
    n          = length(Y),
    stringsAsFactors = FALSE
  )

  if (model_type == "continuous") {
    metrics$rmse <- sqrt(mean((Y - yhat)^2, na.rm = TRUE))
    metrics$mae  <- mean(abs(Y - yhat), na.rm = TRUE)
    metrics$r2   <- 1 - sum((Y - yhat)^2) / sum((Y - mean(Y))^2)
    metrics$r    <- cor(Y, yhat, use = "complete.obs")
  }

  if (model_type == "binary") {
    # AUC
    metrics$auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat, quiet = TRUE))),
      error = function(e) NA_real_
    )
    # Brier score
    metrics$brier <- mean((Y - yhat)^2, na.rm = TRUE)
  }

  # Prevalence via threshold (for continuous models)
  if (model_type == "continuous" && !is.null(oc$cutoff)) {
    pred_bin <- apply_threshold(yhat, oc$cutoff, oc$cutoff_dir)
    metrics$prev_pred <- mean(pred_bin, na.rm = TRUE)
    metrics$prev_obs  <- if (!is.null(oc$binary)) mean(res[[oc$binary]], na.rm = TRUE) else NA
  }

  metrics
}


#' Aggregate individual-level predictions to Admin1 prevalence
#'
#' @param sl_result Output from DHS_SL_clustered
#' @param outcome_data Output from build_outcome_dataset
#' @param cc Country config
#' @param oc Outcome config
#' @param model_type "binary_prob" or "cont_threshold"
#' @return data.frame with Admin1-level prevalence
aggregate_admin1 <- function(sl_result, outcome_data, cc, oc, model_type = "binary_prob") {
  if (is.null(sl_result)) return(NULL)

  d   <- outcome_data$data
  res <- sl_result$res

  # Match predictions back to data by row
  if (model_type == "binary_prob") {
    d$pred_deficient <- res$yhat_full
  } else {
    d$pred_deficient <- apply_threshold(res$yhat_full, oc$cutoff, oc$cutoff_dir)
  }

  admin1_col <- cc$admin1_col

  # Use survey weights for prevalence estimation if available.
  # All MN surveys (Ghana GMS 2017, Gambia GMNS 2018, SL SLMS 2013, Malawi)
  # use stratified cluster sampling with unequal probability of selection —
  # unweighted means can be biased.
  wt_col <- cc$weight_col
  has_weights <- !is.null(wt_col) && wt_col %in% colnames(d) &&
                 any(!is.na(d[[wt_col]]) & d[[wt_col]] > 0)

  if (has_weights) {
    d$`.wt` <- as.numeric(d[[wt_col]])
    d$`.wt`[is.na(d$`.wt`) | d$`.wt` <= 0] <- 1
  } else {
    d$`.wt` <- 1
  }

  admin1_prev <- d %>%
    dplyr::filter(!is.na(.data[[admin1_col]])) %>%
    dplyr::group_by(.data[[admin1_col]]) %>%
    dplyr::summarise(
      n          = dplyr::n(),
      sl_prev    = stats::weighted.mean(pred_deficient, `.wt`, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::rename(Admin1 = dplyr::all_of(admin1_col))

  # Attach observed binary prevalence (survey-weighted) where available
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    obs_prev <- d %>%
      dplyr::filter(!is.na(.data[[admin1_col]])) %>%
      dplyr::group_by(.data[[admin1_col]]) %>%
      dplyr::summarise(
        obs_prev = stats::weighted.mean(.data[[oc$binary]], `.wt`, na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      dplyr::rename(Admin1 = dplyr::all_of(admin1_col))

    admin1_prev <- dplyr::left_join(admin1_prev, obs_prev, by = "Admin1")
  }

  # Compute error metrics
  if ("obs_prev" %in% colnames(admin1_prev)) {
    admin1_prev <- admin1_prev %>%
      dplyr::mutate(
        abs_error     = abs(sl_prev - obs_prev),
        pct_error     = dplyr::if_else(
          obs_prev > 0.01,
          100 * abs(sl_prev - obs_prev) / obs_prev,
          NA_real_
        )
      )
  }

  admin1_prev$outcome    <- oc$tag
  admin1_prev$model_type <- model_type
  admin1_prev
}
