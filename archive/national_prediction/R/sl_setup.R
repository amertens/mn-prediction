# =============================================================================
# National Prediction Pipeline — SuperLearner Setup (sl3)
# =============================================================================

#' Build the full SuperLearner stack
#'
#' 7-learner stack: mean, lasso, ridge, ranger, xgboost, plus screened pipelines.
#' Matches the legacy slfull from src/0-SL-setup.R.
#' @return sl3 Lrnr_sl object with squared error loss
build_sl_full <- function() {
  # Base learners
  lrn_mean    <- sl3::make_learner(sl3::Lrnr_mean)
  lrn_lasso   <- sl3::Lrnr_glmnet$new(alpha = 1)
  lrn_ridge   <- sl3::Lrnr_glmnet$new(alpha = 0)
  lrn_rf      <- sl3::Lrnr_ranger$new()
  lrn_xgb     <- sl3::Lrnr_xgboost$new()

  # Screeners
  rf_importance    <- sl3::Lrnr_ranger$new(importance = "impurity_corrected")
  screen_rf_top20  <- sl3::Lrnr_screener_importance$new(
    learner = rf_importance, num_screen = 20
  )
  screen_lasso     <- sl3::Lrnr_screener_coefs$new(
    learner = lrn_lasso, threshold = 0
  )
  screen_cor       <- sl3::Lrnr_screener_correlation$new(
    type = "threshold", pvalue_threshold = 0.1, min_screen = 2
  )

  # Screened pipelines
  lrn_hal <- sl3::Lrnr_hal9001$new(max_degree = 2, num_knots = c(3, 2), nfolds = 5)
  screen_rf_top20_hal_stack <- sl3::Pipeline$new(
    screen_rf_top20,
    sl3::Stack$new(lrn_hal, lrn_rf)
  )

  lrn_glm      <- sl3::Lrnr_glm$new()
  lrn_bayes    <- sl3::Lrnr_bayesglm$new()
  lrn_gam      <- sl3::Lrnr_gam$new()
  lrn_polspline <- sl3::Lrnr_polspline$new()
  ranger_lrnr  <- sl3::Lrnr_ranger$new()

  screen_lasso_stack <- sl3::Pipeline$new(
    screen_lasso,
    sl3::Stack$new(lrn_glm, lrn_rf, lrn_bayes, lrn_gam, ranger_lrnr,
                   lrn_xgb, lrn_polspline)
  )
  screen_cor_stack <- sl3::Pipeline$new(
    screen_cor,
    sl3::Stack$new(lrn_glm, ranger_lrnr, lrn_xgb)
  )

  # Full stack
  full_stack <- sl3::Stack$new(
    lrn_mean,
    lrn_lasso,
    lrn_ridge,
    lrn_rf,
    screen_rf_top20_hal_stack,
    screen_lasso_stack,
    screen_cor_stack
  )

  # Metalearner

  metalearner <- sl3::make_learner(sl3::Lrnr_nnls)

  sl3::make_learner(
    sl3::Lrnr_sl,
    learners      = full_stack,
    loss_function = sl3::loss_squared_error,
    metalearner   = metalearner
  )
}


#' Build the simple SuperLearner stack
#'
#' 2-learner stack: mean + lasso. Used as fallback for outcomes where the full
#' stack fails (e.g., Vitamin B12, Vitamin D due to small sample sizes).
#' @return sl3 Lrnr_sl object with squared error loss
build_sl_simple <- function() {
  lrn_mean  <- sl3::make_learner(sl3::Lrnr_mean)
  lrn_lasso <- sl3::Lrnr_glmnet$new()

  stack <- sl3::make_learner(sl3::Stack, lrn_mean, lrn_lasso)
  metalearner <- sl3::make_learner(sl3::Lrnr_nnls)

  sl3::make_learner(
    sl3::Lrnr_sl,
    learners      = stack,
    loss_function = sl3::loss_squared_error,
    metalearner   = metalearner
  )
}


#' Get the appropriate SL stack by name
#'
#' @param stack_name Either "full" or "simple"
#' @return sl3 Lrnr_sl object
get_sl_stack <- function(stack_name = "full") {
  switch(stack_name,
    full   = build_sl_full(),
    simple = build_sl_simple(),
    stop("Unknown stack name: ", stack_name, ". Use 'full' or 'simple'.")
  )
}
