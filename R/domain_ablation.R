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
    # Use the underlying mlr3superlearner object directly — the R6 wrapper
    # predict method is more reliable than the closure-based wrapper after
    # serialization (R6 objects can lose state when saved/loaded via targets).
    yhat_perm <- tryCatch({
      perm_df <- data.frame(dt_perm)
      # Try direct mlr3superlearner predict first
      mlr3_obj <- fit$mlr3_fit
      if (!is.null(mlr3_obj) && inherits(mlr3_obj, "mlr3superlearner")) {
        as.numeric(predict(mlr3_obj, perm_df))
      } else {
        # Fallback to wrapper predict
        as.numeric(fit$predict(perm_df))
      }
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
#' For each predictor domain, shuffles that domain's columns and re-evaluates.
#' Because mlr3 R6 learner objects don't survive RDS serialization (targets
#' saves/loads objects between steps), we refit a lightweight glmnet model
#' on the stored training data for permutation predictions. This is fast
#' (~seconds) and avoids the stale-object bug.
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

  # Get the model's actual covariates and training data
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

  train_data <- fit_obj$train_data
  if (is.null(train_data)) {
    cat(sprintf("[permutation] %s — %s | No train_data, skipping\n",
                cc$country, oc$tag))
    return(data.frame())
  }
  train_data <- data.frame(train_data)

  Y <- as.numeric(train_data$Y)

  # ── Refit a lightweight glmnet for permutation predictions ───────────────
  # mlr3 R6 objects lose trained state after RDS serialization. Instead of
  # relying on the (potentially stale) SL object, refit a cv.glmnet on the
  # same training data. This is fast and gives a reasonable proxy for the
  # full SL's feature usage.
  X_mat <- as.matrix(train_data[, final_covars, drop = FALSE])
  # Replace NaN/Inf
  X_mat[!is.finite(X_mat)] <- 0

  family <- if (use_binary) "binomial" else "gaussian"
  perm_model <- tryCatch({
    set.seed(params$seed)
    glmnet::cv.glmnet(x = X_mat, y = Y, family = family, alpha = 0.5,
                       nfolds = min(10, max(3, floor(nrow(X_mat) / 20))))
  }, error = function(e) {
    cat(sprintf("[permutation] glmnet refit failed: %s\n", e$message))
    NULL
  })
  if (is.null(perm_model)) return(data.frame())

  # Baseline: predict on unpermuted data with the refit model
  yhat_baseline <- as.numeric(predict(perm_model, X_mat, s = "lambda.min",
                                       type = "response"))
  baseline <- compute_perf_metrics(Y, yhat_baseline, use_binary)
  baseline$domain_removed <- "none"
  baseline$n_removed      <- 0L
  baseline$n_perm_valid   <- NA_integer_

  cat(sprintf("[permutation] Baseline AUC (glmnet refit): %.3f\n", baseline$auc))

  # Classify variables into conceptual domains
  source_classify <- tryCatch(
    classify_variables(final_covars),
    error = function(e) {
      cat(sprintf("  [permutation] classify_variables failed: %s\n", e$message))
      setNames(rep("Unknown", length(final_covars)), final_covars)
    }
  )

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

  # Permute each domain using the refit glmnet
  perm_list <- lapply(domains, function(dom) {
    dom_cols <- intersect(domain_vars[[dom]], final_covars)
    if (length(dom_cols) == 0) return(NULL)

    cat(sprintf("  Permuting %s (%d vars)\n", dom, length(dom_cols)))

    perm_results <- list()
    for (r in seq_len(n_perm)) {
      set.seed(params$seed + r)
      X_perm <- X_mat
      for (col in dom_cols) {
        col_idx <- match(col, final_covars)
        if (!is.na(col_idx)) {
          X_perm[, col_idx] <- X_perm[sample.int(nrow(X_perm)), col_idx]
        }
      }
      yhat_perm <- tryCatch(
        as.numeric(predict(perm_model, X_perm, s = "lambda.min", type = "response")),
        error = function(e) NULL
      )
      if (!is.null(yhat_perm) && length(yhat_perm) == length(Y)) {
        perm_results[[r]] <- compute_perf_metrics(Y, yhat_perm, use_binary)
      }
    }

    if (length(perm_results) == 0) return(NULL)
    avg <- as.data.frame(lapply(dplyr::bind_rows(perm_results), mean, na.rm = TRUE))
    avg$domain_removed <- dom
    avg$n_removed      <- length(dom_cols)
    avg$n_perm_valid   <- length(perm_results)
    avg
  })

  all_results <- dplyr::bind_rows(c(list(baseline), Filter(Negate(is.null), perm_list)))
  all_results$outcome <- oc$tag
  all_results$country <- cc$country

  all_results
}
