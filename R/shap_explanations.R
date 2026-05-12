# =============================================================================
# R/shap_explanations.R
#
# SHAP-based feature attributions for SuperLearner predictions.
#
# We extract feature contributions from the tree-based base learners
# (xgboost and ranger) within the ensemble. SHAP values for the full
# SuperLearner are approximated by weighting each base learner's SHAP
# contributions by their NNLS metalearner weight (when available) or
# falling back to xgboost-only attributions.
#
# Two outputs:
#   - District-level top factors: for each Admin-2 polygon, the top N
#     predictors driving the predicted prevalence above (or below) the
#     country mean. Used by the dashboard's District Profile cards.
#   - Global feature importance: across-district mean |SHAP| per predictor.
# =============================================================================


#' Compute SHAP-based feature attributions for an SL fit
#'
#' Uses the xgboost base learner inside the SuperLearner if available
#' (it's the most efficient SHAP-supported learner). For each observation,
#' returns the top predictors by absolute SHAP contribution.
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param sl_fit Output from fit_mlr3_models()
#' @param cc Country config
#' @param oc Outcome config
#' @param n_top Number of top predictors to return per observation
#' @return list with district_factors (Admin-2 × top-N table) and
#'   global_importance (predictor × mean |SHAP| table).
compute_shap_explanations <- function(outcome_data, sl_fit, cc, oc, n_top = 5L) {

  fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
  if (is.null(fit)) {
    return(list(district_factors = NULL, global_importance = NULL,
                method = "unavailable"))
  }

  d <- outcome_data$data
  admin2_col <- cc$admin2_col %||% "Admin2"
  admin1_col <- cc$admin1_col %||% "Admin1"

  # ── Find the xgboost learner inside the SL ──
  xgb_model <- NULL
  if (!is.null(fit$sl_fit) && !is.null(fit$sl_fit$fits)) {
    learners <- fit$sl_fit$fits$models %||% fit$sl_fit$fits
    for (lrn_name in names(learners)) {
      if (grepl("xgb", lrn_name, ignore.case = TRUE)) {
        candidate <- learners[[lrn_name]]
        if (inherits(candidate, c("xgb.Booster", "Booster"))) {
          xgb_model <- candidate
          break
        }
        # mlr3 wraps the booster — check $model and $state$model
        if (!is.null(candidate$model) && inherits(candidate$model, "xgb.Booster")) {
          xgb_model <- candidate$model
          break
        }
      }
    }
  }

  if (is.null(xgb_model)) {
    cat(sprintf("[shap] %s — %s | xgboost learner not available, skipping\n",
                cc$country, oc$tag))
    return(list(district_factors = NULL, global_importance = NULL,
                method = "no_xgb_learner"))
  }

  # ── Recreate the predictor matrix used by the SL ──
  Xvars <- fit$Xvars %||% outcome_data$Xvars
  if (is.null(Xvars)) {
    return(list(district_factors = NULL, global_importance = NULL,
                method = "no_xvars"))
  }

  X <- as.matrix(d[, intersect(Xvars, colnames(d)), drop = FALSE])
  X[!is.finite(X)] <- NA

  # ── Compute SHAP values ──
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    return(list(district_factors = NULL, global_importance = NULL,
                method = "xgboost_not_installed"))
  }

  shap_values <- tryCatch({
    suppressWarnings(predict(xgb_model, X, predcontrib = TRUE,
                              approxcontrib = FALSE))
  }, error = function(e) {
    cat(sprintf("[shap] predict(predcontrib=TRUE) failed: %s\n", e$message))
    NULL
  })

  if (is.null(shap_values)) {
    return(list(district_factors = NULL, global_importance = NULL,
                method = "shap_extract_failed"))
  }

  # Drop the BIAS column (last column in xgboost SHAP output)
  bias_col <- grep("^BIAS$", colnames(shap_values))
  if (length(bias_col) > 0) shap_values <- shap_values[, -bias_col, drop = FALSE]

  # ── Aggregate to Admin-2 ──
  admin2_vec <- d[[admin2_col]]
  admin1_vec <- if (admin1_col %in% colnames(d)) d[[admin1_col]] else NA_character_

  district_rows <- list()
  for (a2 in unique(admin2_vec[!is.na(admin2_vec)])) {
    idx <- which(admin2_vec == a2)
    if (length(idx) == 0) next
    mean_shap <- colMeans(shap_values[idx, , drop = FALSE], na.rm = TRUE)
    # Order by absolute mean SHAP, take top n
    abs_order <- order(abs(mean_shap), decreasing = TRUE)
    top_idx <- abs_order[seq_len(min(n_top, length(abs_order)))]

    a1 <- if (!all(is.na(admin1_vec))) admin1_vec[idx[1]] else NA_character_

    district_rows[[a2]] <- data.frame(
      country  = cc$country,
      outcome  = oc$tag,
      Admin1   = a1,
      Admin2   = a2,
      rank     = seq_along(top_idx),
      variable = colnames(shap_values)[top_idx],
      shap     = mean_shap[top_idx],
      direction = ifelse(mean_shap[top_idx] > 0, "increases risk",
                          "decreases risk"),
      stringsAsFactors = FALSE
    )
  }
  district_factors <- if (length(district_rows) > 0)
    do.call(rbind, district_rows) else NULL

  # ── Global importance (across all observations) ──
  global_imp <- data.frame(
    variable = colnames(shap_values),
    mean_abs_shap = colMeans(abs(shap_values), na.rm = TRUE),
    mean_shap = colMeans(shap_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  global_imp <- global_imp[order(-global_imp$mean_abs_shap), ]
  global_imp$rank <- seq_len(nrow(global_imp))
  global_imp$country <- cc$country
  global_imp$outcome <- oc$tag

  cat(sprintf("[shap] %s — %s | %d districts × top-%d factors extracted\n",
              cc$country, oc$tag, nrow(district_factors %||% data.frame()),
              n_top))

  list(
    district_factors  = district_factors,
    global_importance = global_imp,
    method = "xgboost_treeshap"
  )
}
