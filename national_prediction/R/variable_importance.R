# =============================================================================
# National Prediction Pipeline — Variable Importance
# =============================================================================

#' Compute permutation-based variable importance from an sl3 fit
#'
#' Uses sl3's importance() function to compute permutation importance.
#' Ported from src/0-functions.R calc_importance().
#'
#' @param sl_fit Trained sl3 Lrnr_sl object
#' @param eval_fun Loss function (default: loss_squared_error)
#' @param importance_metric Character: "ratio" or "difference"
#' @param n_vars Integer: number of top variables to return (default 20)
#' @param covariate_groups Optional named list of covariate groups
#' @return List with varimp (full table), varimp_top (top N), plot
calc_importance <- function(sl_fit,
                            eval_fun = sl3::loss_squared_error,
                            importance_metric = "ratio",
                            n_vars = 20L,
                            covariate_groups = NULL) {
  set.seed(983)

  varimp <- sl3::importance(
    fit               = sl_fit,
    eval_fun          = eval_fun,
    type              = "permute",
    importance_metric = importance_metric,
    covariate_groups  = covariate_groups
  )

  varimp_top <- varimp %>%
    dplyr::arrange(dplyr::desc(dplyr::across(2))) %>%
    utils::head(n = n_vars)

  p <- sl3::importance_plot(x = varimp_top)

  list(
    varimp     = varimp,
    varimp_top = varimp_top,
    plot       = p
  )
}


#' Safely compute variable importance, returning NULL on error
#'
#' @param fit_result Result from run_national_sl()
#' @param outcome Character: outcome label (for warning messages)
#' @return List from calc_importance(), or NULL on failure
safe_calc_importance <- function(fit_result, outcome) {
  if (inherits(fit_result, "try-error") || is.null(fit_result)) {
    warning("Skipping variable importance for ", outcome, ": model fit failed")
    return(NULL)
  }

  tryCatch(
    calc_importance(fit_result$sl_fit),
    error = function(e) {
      warning("Variable importance failed for ", outcome, ": ", e$message)
      NULL
    }
  )
}
