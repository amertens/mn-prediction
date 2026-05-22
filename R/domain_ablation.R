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

  # Extract components — handle both sl3 and mlr3 wrapper formats
  fit      <- sl_fit_obj$sl_fit
  Y        <- as.numeric(sl_fit_obj$res$Y)

  # Get covariates and task data from whichever format is available
  if (!is.null(sl_fit_obj$task) && !is.null(sl_fit_obj$task$nodes)) {
    # sl3 format
    covars  <- sl_fit_obj$task$nodes$covariates
    task_dt <- data.table::copy(sl_fit_obj$task$data)
  } else {
    # mlr3 wrapper format — use stored train_data
    covars  <- sl_fit_obj$Xvars
    task_dt <- sl_fit_obj$train_data

    if (is.null(task_dt)) {
      cat(sprintf("  [permute] No train_data stored for %s, skipping\n", domain_name))
      return(NULL)
    }
    task_dt <- data.table::as.data.table(task_dt)
  }

  # Which domain columns are actually in the final task covariates?
  perm_cols <- intersect(domain_cols, covars)
  if (length(perm_cols) == 0) return(NULL)

  perm_results <- list()

  for (r in seq_len(n_perm)) {
    set.seed(seed + r)
    dt_perm <- data.table::copy(task_dt)

    # Shuffle each domain column independently (breaks the relationship
    # between these predictors and the outcome while preserving marginals)
    n_rows <- nrow(dt_perm)
    for (col in perm_cols) {
      if (col %in% colnames(dt_perm)) {
        data.table::set(dt_perm, j = col, value = dt_perm[[col]][sample.int(n_rows)])
      }
    }

    # Predict on permuted data
    yhat_perm <- tryCatch({
      perm_df <- data.frame(dt_perm)
      p <- as.numeric(fit$predict(perm_df))
      # Detect degenerate predictions (all identical = predict failed silently)
      if (length(unique(round(p, 6))) <= 1) {
        cat(sprintf("    [permute] WARNING: degenerate predictions for %s rep %d (all %.4f)\n",
                    domain_name, r, p[1]))
      }
      p
    }, error = function(e) {
      cat(sprintf("    [permute] Predict failed for %s rep %d: %s\n",
                  domain_name, r, e$message))
      NULL
    })
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

  # Get the model's actual covariates (what it was trained on)
  final_covars <- fit_obj$Xvars
  if (is.null(final_covars) || length(final_covars) == 0) {
    final_covars <- tryCatch(fit_obj$task$nodes$covariates,
                             error = function(e) character())
  }
  if (length(final_covars) == 0) {
    cat(sprintf("[permutation] %s — %s | No covariates found, skipping\n",
                cc$country, oc$tag))
    return(data.frame())
  }

  # Classify variables into conceptual domains using the mapping file
  source_classify <- tryCatch(
    classify_variables(final_covars),
    error = function(e) {
      cat(sprintf("  [permutation] classify_variables failed: %s\n", e$message))
      setNames(rep("Unknown", length(final_covars)), final_covars)
    }
  )

  # Build domain_vars list: domain name -> vector of variable names
  domain_vars <- split(names(source_classify), source_classify)
  domain_vars <- domain_vars[names(domain_vars) != "Unknown"]
  domain_vars <- domain_vars[names(domain_vars) != "Missingness Indicators"]

  domains <- names(domain_vars)
  domains <- domains[sapply(domains, function(nm) {
    length(intersect(domain_vars[[nm]], final_covars)) > 0
  })]

  n_domains <- length(domains)
  n_perm <- params$n_perm %||% 5L

  cat(sprintf("\n[permutation] %s — %s | %d domains, %d permutations each\n",
              cc$country, oc$tag, n_domains, n_perm))

  # Baseline metrics. Permutation importance compares the SAME fitted model's
  # predictions with vs without a domain shuffled; the permuted predictions
  # below are in-sample (fit$predict on permuted data), so the baseline must
  # also be in-sample. yhat_full is now out-of-fold, so use yhat_insample
  # (falls back to yhat_full for models fit before that column existed).
  yhat_base_is <- fit_obj$res$yhat_insample %||% fit_obj$res$yhat_full
  baseline <- compute_perf_metrics(fit_obj$res$Y, yhat_base_is, use_binary)
  baseline$domain_removed <- "none"
  baseline$n_removed      <- 0L
  baseline$n_perm_valid   <- NA_integer_

  # Permute each domain
  perm_list <- lapply(domains, function(dom) {
    dom_cols <- domain_vars[[dom]]
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
