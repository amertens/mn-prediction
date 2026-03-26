# =============================================================================
# R/domain_ablation.R
#
# Permutation-based domain importance using the fitted SL model.
# Shuffles each domain's columns and measures performance degradation,
# WITHOUT refitting the model. Much cheaper than refit-based ablation.
#
# The refit-based version is saved in R/domain_ablation_refit.R.bak
# =============================================================================

#' Compute performance metrics from predictions
#'
#' @param Y Observed outcome
#' @param yhat Predicted values
#' @param use_binary Logical
#' @return 1-row data.frame with auc/brier or rmse/mae/r2
compute_perf_metrics <- function(Y, yhat, use_binary) {
  if (use_binary) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat, quiet = TRUE))),
      error = function(e) NA_real_
    )
    data.frame(auc = auc_val,
               brier = mean((Y - yhat)^2, na.rm = TRUE),
               stringsAsFactors = FALSE)
  } else {
    data.frame(
      rmse = sqrt(mean((Y - yhat)^2, na.rm = TRUE)),
      mae  = mean(abs(Y - yhat), na.rm = TRUE),
      r2   = 1 - sum((Y - yhat)^2, na.rm = TRUE) /
                  sum((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}


#' Permute one domain and predict using the existing fitted model
#'
#' Shuffles the columns belonging to one domain in the SL task data,
#' generates new predictions from the already-fitted model, and
#' computes performance metrics.
#'
#' @param sl_fit_obj The trained SL object (sl_fit$cont_fit or sl_fit$bin_fit)
#' @param domain_cols Character vector of column names to permute
#' @param domain_name Name of the domain (for labelling)
#' @param use_binary Logical
#' @param n_perm Number of permutation replicates (default 5)
#' @param seed Random seed
#' @return 1-row data.frame with domain_removed, n_removed, and mean metrics
permute_domain <- function(sl_fit_obj, domain_cols, domain_name,
                           use_binary, n_perm = 5, seed = 12345L) {

  task     <- sl_fit_obj$task
  fit      <- sl_fit_obj$sl_fit
  Y        <- sl_fit_obj$res$Y
  covars   <- task$nodes$covariates

  # Which domain columns are actually in the final task covariates?
  perm_cols <- intersect(domain_cols, covars)
  if (length(perm_cols) == 0) return(NULL)

  # Get the underlying data.table
  task_dt <- data.table::copy(task$data)

  perm_results <- list()

  for (r in seq_len(n_perm)) {
    set.seed(seed + r)
    dt_perm <- data.table::copy(task_dt)

    # Shuffle each domain column independently (breaks the relationship
    # between these predictors and the outcome while preserving marginals)
    n_rows <- nrow(dt_perm)
    for (col in perm_cols) {
      data.table::set(dt_perm, j = col, value = dt_perm[[col]][sample.int(n_rows)])
    }

    # Build a new task with permuted data and predict
    perm_task <- tryCatch(
      sl3::sl3_Task$new(
        data       = dt_perm,
        covariates = covars,
        outcome    = "Y",
        id         = task$nodes$id,
        folds      = task$folds
      ),
      error = function(e) NULL
    )
    if (is.null(perm_task)) next

    yhat_perm <- tryCatch(
      as.numeric(fit$predict(perm_task)),
      error = function(e) NULL
    )
    if (is.null(yhat_perm) || length(yhat_perm) != length(Y)) next

    perm_results[[r]] <- compute_perf_metrics(Y, yhat_perm, use_binary)
  }

  if (length(perm_results) == 0) return(NULL)

  # Average metrics across permutation replicates
  perm_df <- dplyr::bind_rows(perm_results)
  avg <- as.data.frame(lapply(perm_df, mean, na.rm = TRUE))
  avg$domain_removed <- domain_name
  avg$n_removed      <- length(perm_cols)
  avg$n_perm_valid   <- length(perm_results)

  avg
}


#' Run permutation-based domain importance for one outcome
#'
#' For each predictor domain, shuffles that domain's columns in the
#' fitted model's task data and measures performance degradation.
#' Uses the EXISTING fitted model — no refitting.
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param sl_fit Output from fit_sl_models()
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learners Not used (kept for interface compatibility)
#' @param params Pipeline parameters (for seed, n_perm)
#' @return data.frame with one row per domain (including "none" baseline)
run_domain_ablation <- function(outcome_data, sl_fit, cc, oc, sl_learners, params) {

  use_binary <- !is.null(sl_fit$bin_fit)
  fit_obj    <- if (use_binary) sl_fit$bin_fit else sl_fit$cont_fit

  if (is.null(fit_obj)) {
    cat(sprintf("[permutation] %s — %s | No fitted model, skipping\n",
                cc$country, oc$tag))
    return(data.frame())
  }

  domains <- names(outcome_data$domain_vars)
  domains <- setdiff(domains, "GW")  # GW not used in primary model

  # Only check domains with columns in the final fitted model's covariates
  final_covars <- fit_obj$task$nodes$covariates
  domains <- domains[sapply(domains, function(nm) {
    length(intersect(outcome_data$domain_vars[[nm]], final_covars)) > 0
  })]

  n_domains <- length(domains)
  n_perm <- params$n_perm %||% 5L

  cat(sprintf("\n[permutation] %s — %s | %d domains, %d permutations each\n",
              cc$country, oc$tag, n_domains, n_perm))

  # Baseline metrics (unpermuted CV predictions)
  baseline <- compute_perf_metrics(fit_obj$res$Y, fit_obj$res$yhat_full, use_binary)
  baseline$domain_removed <- "none"
  baseline$n_removed      <- 0L
  baseline$n_perm_valid   <- NA_integer_

  # Permute each domain
  perm_list <- lapply(domains, function(dom) {
    dom_cols <- outcome_data$domain_vars[[dom]]
    cat(sprintf("  Permuting %s (%d vars in model)\n", dom,
                length(intersect(dom_cols, final_covars))))
    permute_domain(fit_obj, dom_cols, dom, use_binary,
                   n_perm = n_perm, seed = params$seed)
  })

  all_results <- dplyr::bind_rows(c(list(baseline), Filter(Negate(is.null), perm_list)))
  all_results$outcome <- oc$tag
  all_results$country <- cc$country

  all_results
}
