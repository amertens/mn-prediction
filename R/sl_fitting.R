# =============================================================================
# R/sl_fitting.R
#
# Functions for fitting SuperLearner models (continuous + binary).
# These wrap the existing DHS_SL_clustered() machinery as pure functions.
# =============================================================================

#' Set up the SuperLearner learner stacks
#'
#' Creates the sl3 learner objects used by all model-fitting functions.
#'
#' Both stacks are designed for the "many weak predictors, rare outcome" setting:
#' ~200+ proxy predictors (IHME, MAP, GEE, DHS, MICS) with binary outcomes that
#' can have prevalence as low as 3%. Key design choices:
#'
#' - Lasso (alpha=1) for aggressive variable selection
#' - Prescreened RF/XGBoost to avoid splitting on noise
#' - Ridge to capture distributed weak signal across correlated predictors
#' - Conservative min.node.size / min_child_weight for rare outcomes
#'
#' In "fast" mode, uses a 4-learner stack (~5-10 min per outcome).
#' In "full" mode, uses an 8-learner stack with multiple prescreening
#' strategies (~30-60 min per outcome).
#'
#' @param params Pipeline parameters (from get_pipeline_params())
#' @return list with slmod (continuous stack) and slmod2_bin (binary stack)
setup_sl_learners <- function(params) {

  stack_mode <- params$sl_stack %||% "fast"
  metalearner <- sl3::make_learner(sl3::Lrnr_nnls)

  if (stack_mode == "fast") {
    # ── Fast stack: 3 learners ─────────────────────────────────────────
    # Focused on regularization for many weak predictors.
    # All learners are robust to high p>>n settings (no screener pipelines
    # that can fail when lasso drops all coefficients or ranger hits NA).
    cat("[sl_learners] Using FAST stack (mean + lasso + ranger)\n")

    lrnr_mean  <- sl3::make_learner(sl3::Lrnr_mean)

    # Pure lasso — aggressive variable selection, best single-model selector.
    # Handles p >> n natively via regularization.
    lrnr_lasso <- sl3::Lrnr_glmnet$new(alpha = 1)

    # Ranger with conservative settings for high-p.
    # min.node.size=10 prevents overfitting; no prescreening needed since
    # ranger handles high-p natively and screener pipelines can fail when
    # p >> n (e.g., Lrnr_screener_importance errors on NA importance).
    lrnr_ranger <- sl3::Lrnr_ranger$new(
      num.trees     = 500,
      min.node.size = 10,
      importance    = "none"
    )

    stack_fast <- sl3::make_learner(
      sl3::Stack,
      lrnr_mean,
      lrnr_lasso,
      lrnr_ranger
    )

    slmod <- sl3::make_learner(
      sl3::Lrnr_sl,
      learners      = stack_fast,
      loss_function = sl3::loss_squared_error,
      metalearner   = metalearner
    )

    slmod_bin <- sl3::make_learner(
      sl3::Lrnr_sl,
      learners      = stack_fast,
      loss_function = sl3::loss_loglik_binomial,
      metalearner   = metalearner
    )

    list(slmod = slmod, slmod2_bin = slmod_bin)

  } else {
    # ── Full stack: 7 learners ─────────────────────────────────────────
    # Designed for many weak predictors + rare binary outcomes.
    # Avoids screener pipelines (Lrnr_screener_coefs, Lrnr_screener_importance)
    # that crash when p >> n and lasso drops all coefficients or importance
    # values contain NA. Instead, all learners handle high-p natively.
    cat("[sl_learners] Using FULL stack (7 learners, all high-p robust)\n")

    lrnr_mean <- sl3::make_learner(sl3::Lrnr_mean)

    # ── Regularized linear models (different alpha) ──
    lrnr_lasso    <- sl3::Lrnr_glmnet$new(alpha = 1)
    lrnr_enet_075 <- sl3::Lrnr_glmnet$new(alpha = 0.75)
    lrnr_ridge    <- sl3::Lrnr_glmnet$new(alpha = 0)

    # ── Random forests with conservative settings ──
    # min.node.size prevents overfitting on rare outcomes
    rf_conservative <- sl3::Lrnr_ranger$new(
      num.trees     = 1000,
      min.node.size = 10,
      importance    = "none"
    )

    rf_deeper <- sl3::Lrnr_ranger$new(
      num.trees     = 1000,
      min.node.size = 5,
      importance    = "none"
    )

    # ── XGBoost with conservative regularization ──
    lrnr_xgb <- sl3::Lrnr_xgboost$new(
      max_depth        = 3,
      eta              = 0.05,
      nrounds          = 300,
      min_child_weight = 20,
      subsample        = 0.8,
      colsample_bytree = 0.5,
      eval_metric      = "logloss"
    )

    stack_full <- sl3::make_learner(
      sl3::Stack,
      lrnr_mean,
      lrnr_lasso,
      lrnr_enet_075,
      lrnr_ridge,
      rf_conservative,
      rf_deeper,
      lrnr_xgb
    )

    slmod <- sl3::make_learner(
      sl3::Lrnr_sl,
      learners      = stack_full,
      loss_function = sl3::loss_squared_error,
      metalearner   = metalearner
    )

    slmod_bin <- sl3::make_learner(
      sl3::Lrnr_sl,
      learners      = stack_full,
      loss_function = sl3::loss_loglik_binomial,
      metalearner   = metalearner
    )

    list(slmod = slmod, slmod2_bin = slmod_bin)
  }
}


#' Fit continuous + binary SL for one outcome
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @return list with cont_fit and bin_fit (each a DHS_SL_clustered result or NULL)
fit_sl_models <- function(outcome_data, cc, oc, sl_learners, params) {

  # Source helper functions (DHS_SL_clustered, etc.)
  source(here::here("src", "analysis", "sl_helpers.R"), local = FALSE)
  source(here::here("src", "0-functions.R"), local = FALSE)

  # Make cfg available for sl_helpers (which references cfg$seed, cfg$admin1, etc.)
  cfg <- list(
    seed       = params$seed,
    admin1     = cc$admin1_col,
    cluster_id = cc$cluster_id,
    K          = params$K
  )
  assign("cfg", cfg, envir = globalenv())

  d     <- outcome_data$data

  # PRIMARY: proxy-only predictors (no gw_ survey variables).
  # These are the predictors available everywhere (DHS, MICS, GEE, IHME, etc.)
  # and are the basis for extrapolation to unsurveyed areas.
  # gw_ variables are from the same survey as the outcome — using them would be

  # circular (predicting survey results from survey results) and inflates AUC.
  Xvars <- outcome_data$Xvars

  cat(sprintf("\n[fit_sl] %s — %s | n = %d | p = %d (proxy-only, no gw_)\n",
              cc$country, oc$tag, nrow(d), length(Xvars)))

  # Bail early if no data (e.g., outcome not collected in this country)
  if (nrow(d) == 0 || length(Xvars) == 0) {
    cat("  SKIPPING: no observations or no predictors available.\n")
    return(list(cont_fit = NULL, bin_fit = NULL,
                outcome = oc$tag, country = cc$country))
  }

  # ── Continuous model ──
  cat("  Fitting continuous SL...\n")
  t0 <- proc.time()
  cont_fit <- DHS_SL_clustered(
    d          = d,
    Xvars      = Xvars,
    outcome    = oc$continuous,
    population = oc$population,
    id         = cc$cluster_id,
    folds      = params$K,
    sl         = sl_learners$slmod
  )
  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  Continuous done in %.1f min\n", elapsed / 60))

  # ── Binary model ──
  bin_fit <- NULL
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    y_bin <- d[[oc$binary]]
    y_bin_nona <- y_bin[!is.na(y_bin)]
    is_binary <- length(unique(y_bin_nona)) >= 2 &&
                 all(y_bin_nona >= 0 & y_bin_nona <= 1)
    if (!is_binary) {
      cat(sprintf("  SKIPPING binary: values not in [0,1] (range: %.2f-%.2f)\n",
                  min(y_bin_nona), max(y_bin_nona)))
    }
    if (is_binary) {
      cat("  Fitting binary SL...\n")
      t0 <- proc.time()
      bin_fit <- DHS_SL_clustered(
        d          = d,
        Xvars      = Xvars,
        outcome    = oc$binary,
        population = oc$population,
        id         = cc$cluster_id,
        folds      = params$K,
        sl         = sl_learners$slmod2_bin
      )
      elapsed <- (proc.time() - t0)["elapsed"]
      cat(sprintf("  Binary done in %.1f min\n", elapsed / 60))
    }
  }

  list(
    cont_fit = cont_fit,
    bin_fit  = bin_fit,
    outcome  = oc$tag,
    country  = cc$country
  )
}
