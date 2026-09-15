# =============================================================================
# National Prediction Pipeline — Model Performance Evaluation
# =============================================================================

#' Evaluate SuperLearner model performance
#'
#' Computes standard regression metrics from an sl3 fit object: MSE, RMSE, MAE,
#' R-squared (variance-based and vs. mean predictor), correlation, and bias.
#' Ported from src/0-functions.R evaluate_sl_performance().
#'
#' @param res List returned by run_national_sl() — must contain $Y, $yhat_full,
#'   and $cv_risk_w_sl_revere
#' @param outcome Character: label for the outcome
#' @return data.frame (1 row) with 12 performance metrics
evaluate_sl_performance <- function(res, outcome) {
  Y    <- as.numeric(res$Y)
  Yhat <- as.numeric(res$yhat_full)

  # Baseline MSE (mean predictor)
  mse_mean <- as.numeric(
    res$cv_risk_w_sl_revere[res$cv_risk_w_sl_revere$learner == "Lrnr_mean", "MSE"]
  )

  # SuperLearner MSE
  mse <- as.numeric(
    res$cv_risk_w_sl_revere[res$cv_risk_w_sl_revere$learner == "SuperLearner", "MSE"]
  )

  rmse <- sqrt(mse)
  mae  <- mean(abs(Y - Yhat))
  varY <- var(Y)

  # R-squared: variance-based and relative to mean predictor
 R2          <- 1 - mse / varY
  R2_baseline <- 1 - mse / mse_mean
  rel_improve <- (1 - mse / mse_mean) * 100

  # Correlation
  ct        <- cor.test(Y, Yhat)
  corr      <- as.numeric(ct$estimate)
  corr_pval <- as.numeric(ct$p.value)

  # Bias
  bias <- mean(Yhat - Y)

  data.frame(
    outcome            = outcome,
    mse_model          = mse,
    mse_mean           = mse_mean,
    rmse               = rmse,
    mae                = mae,
    outcome_variance   = varY,
    R2_variance_based  = R2,
    R2_vs_mean         = R2_baseline,
    relative_improvement = rel_improve,
    correlation        = corr,
    correlation_pvalue = corr_pval,
    mean_bias          = bias,
    stringsAsFactors   = FALSE
  )
}


#' Safely evaluate performance, returning NA row on error
#'
#' Wraps evaluate_sl_performance() in tryCatch so that failed model fits
#' don't break the pipeline.
#' @param res Result from run_national_sl() (may be a try-error)
#' @param outcome Character: outcome label
#' @return data.frame (1 row)
safe_evaluate_performance <- function(res, outcome) {
  if (inherits(res, "try-error") || is.null(res)) {
    return(data.frame(
      outcome = outcome,
      mse_model = NA_real_, mse_mean = NA_real_, rmse = NA_real_,
      mae = NA_real_, outcome_variance = NA_real_,
      R2_variance_based = NA_real_, R2_vs_mean = NA_real_,
      relative_improvement = NA_real_, correlation = NA_real_,
      correlation_pvalue = NA_real_, mean_bias = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  tryCatch(
    evaluate_sl_performance(res, outcome),
    error = function(e) {
      warning("Performance evaluation failed for ", outcome, ": ", e$message)
      data.frame(
        outcome = outcome,
        mse_model = NA_real_, mse_mean = NA_real_, rmse = NA_real_,
        mae = NA_real_, outcome_variance = NA_real_,
        R2_variance_based = NA_real_, R2_vs_mean = NA_real_,
        relative_improvement = NA_real_, correlation = NA_real_,
        correlation_pvalue = NA_real_, mean_bias = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
}
