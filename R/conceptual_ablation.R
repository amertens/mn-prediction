# =============================================================================
# R/conceptual_ablation.R
#
# Permutation-based variable importance grouped by conceptual domain.
# Uses the fitted SL model (no refitting) — shuffles predictor columns within
# each conceptual domain and measures the drop in performance.
#
# This is much cheaper than the data-source domain ablation (which refits the
# model for each domain), making it feasible for many subgroups.
# =============================================================================

#' Load the conceptual domain mapping from metadata CSV
#'
#' @return data.frame with variable_pattern, data_source, level1_domain,
#'   level2_domain, description, match_type
load_conceptual_domains <- function() {
  path <- here::here("metadata", "variable_conceptual_domains.csv")
  if (!file.exists(path)) stop("Conceptual domain mapping not found: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}


#' Load the DHS indicator lookup (IndicatorId -> Level1, Level2)
#'
#' Built from metadata/dhs_indicators.csv. Contains exact IndicatorId matches
#' and 2-part prefix matches (e.g., CH_VACC -> Child Health / Vaccinations).
#'
#' @return list with `exact` (IndicatorId -> l1, l2) and `prefix` (2-part -> l1, l2)
load_dhs_lookup <- function() {
  exact_path  <- here::here("metadata", "dhs_indicator_lookup.csv")
  prefix_path <- here::here("metadata", "dhs_prefix_lookup.csv")

  exact <- if (file.exists(exact_path)) {
    read.csv(exact_path, stringsAsFactors = FALSE)
  } else data.frame()

  prefix <- if (file.exists(prefix_path)) {
    read.csv(prefix_path, stringsAsFactors = FALSE)
  } else data.frame()

  list(exact = exact, prefix = prefix)
}


#' Extract the DHS IndicatorId from a variable name
#'
#' e.g., "dhs2014_FE_FRTR_W_A10" -> "FE_FRTR_W_A10"
#'       "dhs2016_CH_VACC_C_OP2" -> "CH_VACC_C_OP2"
#' @param varname character
#' @return character (indicator ID without the dhs####_ prefix), or NA
extract_dhs_indicator_id <- function(varname) {
  m <- regmatches(varname, regexpr("^dhs\\d{4}_(.+)$", varname, perl = TRUE))
  if (length(m) == 0) return(NA_character_)
  sub("^dhs\\d{4}_", "", m)
}


#' Classify predictor variables into conceptual Level-1 domains
#'
#' Uses a multi-pass strategy:
#'   1. DHS variables: match IndicatorId exactly against dhs_indicators.csv
#'   2. DHS variables: match 2-part prefix (e.g., CH_VACC) against prefix lookup
#'   3. All variables: fall back to pattern matching from variable_conceptual_domains.csv
#'
#' @param xvars Character vector of predictor variable names
#' @param mapping Output from load_conceptual_domains() (or NULL to auto-load)
#' @return Named character vector: names = xvars, values = level1_domain
classify_variables <- function(xvars, mapping = NULL) {
  if (is.null(mapping)) mapping <- load_conceptual_domains()
  dhs_lk <- load_dhs_lookup()

  result <- setNames(rep("Unknown", length(xvars)), xvars)

  # Pass 1: DHS exact IndicatorId match
  if (nrow(dhs_lk$exact) > 0) {
    exact_l1 <- setNames(dhs_lk$exact$level1, dhs_lk$exact$indicator_id)
    for (i in seq_along(xvars)) {
      iid <- extract_dhs_indicator_id(xvars[i])
      if (!is.na(iid) && iid %in% names(exact_l1)) {
        result[i] <- exact_l1[[iid]]
      }
    }
  }

  # Pass 2: DHS 2-part prefix match (e.g., CH_VACC -> Child Health)
  if (nrow(dhs_lk$prefix) > 0) {
    prefix_l1 <- setNames(dhs_lk$prefix$level1, dhs_lk$prefix$prefix)
    for (i in seq_along(xvars)) {
      if (result[i] != "Unknown") next
      iid <- extract_dhs_indicator_id(xvars[i])
      if (is.na(iid)) next
      parts <- strsplit(iid, "_")[[1]]
      if (length(parts) >= 2) {
        pfx <- paste(parts[1], parts[2], sep = "_")
        if (pfx %in% names(prefix_l1)) {
          result[i] <- prefix_l1[[pfx]]
        }
      }
    }
  }

  # Pass 3: General pattern matching (non-DHS variables, or remaining unknowns)
  for (i in seq_len(nrow(mapping))) {
    pat    <- mapping$variable_pattern[i]
    domain <- mapping$level1_domain[i]
    mtype  <- mapping$match_type[i]

    if (mtype == "exact") {
      idx <- which(xvars == pat)
    } else {
      regex_pat <- gsub("####", "\\\\d{4}", pat, fixed = TRUE)
      regex_pat <- paste0("^", regex_pat)
      idx <- grep(regex_pat, xvars)
    }

    for (j in idx) {
      if (result[j] == "Unknown") result[j] <- domain
    }
  }

  result
}


#' Classify predictor variables into conceptual Level-2 domains
#'
#' Same multi-pass strategy as classify_variables() but returns Level-2.
#'
#' @param xvars Character vector of predictor variable names
#' @param mapping Output from load_conceptual_domains() (or NULL to auto-load)
#' @return Named character vector: names = xvars, values = level2_domain
classify_variables_l2 <- function(xvars, mapping = NULL) {
  if (is.null(mapping)) mapping <- load_conceptual_domains()
  dhs_lk <- load_dhs_lookup()

  result <- setNames(rep("Unknown", length(xvars)), xvars)

  # Pass 1: DHS exact IndicatorId match
  if (nrow(dhs_lk$exact) > 0) {
    exact_l2 <- setNames(dhs_lk$exact$level2, dhs_lk$exact$indicator_id)
    for (i in seq_along(xvars)) {
      iid <- extract_dhs_indicator_id(xvars[i])
      if (!is.na(iid) && iid %in% names(exact_l2)) {
        val <- exact_l2[[iid]]
        if (!is.na(val) && nchar(val) > 0) result[i] <- val
      }
    }
  }

  # Pass 2: DHS 2-part prefix match
  if (nrow(dhs_lk$prefix) > 0) {
    prefix_l2 <- setNames(dhs_lk$prefix$level2, dhs_lk$prefix$prefix)
    for (i in seq_along(xvars)) {
      if (result[i] != "Unknown") next
      iid <- extract_dhs_indicator_id(xvars[i])
      if (is.na(iid)) next
      parts <- strsplit(iid, "_")[[1]]
      if (length(parts) >= 2) {
        pfx <- paste(parts[1], parts[2], sep = "_")
        if (pfx %in% names(prefix_l2)) {
          result[i] <- prefix_l2[[pfx]]
        }
      }
    }
  }

  # Pass 3: General pattern matching
  for (i in seq_len(nrow(mapping))) {
    pat    <- mapping$variable_pattern[i]
    domain <- mapping$level2_domain[i]
    mtype  <- mapping$match_type[i]

    if (mtype == "exact") {
      idx <- which(xvars == pat)
    } else {
      regex_pat <- gsub("####", "\\\\d{4}", pat, fixed = TRUE)
      regex_pat <- paste0("^", regex_pat)
      idx <- grep(regex_pat, xvars)
    }

    for (j in idx) {
      if (result[j] == "Unknown") result[j] <- domain
    }
  }

  result
}


#' Run conceptual domain permutation importance
#'
#' For each conceptual domain, shuffles the columns belonging to that domain
#' and measures the change in a performance metric. Repeats n_perm times
#' and returns the mean drop.
#'
#' @param sl_fit Output from fit_sl_models() — must contain cont_fit and/or bin_fit
#' @param outcome_data Output from build_outcome_dataset()
#' @param oc Outcome config
#' @param level "level1" or "level2" for domain granularity
#' @param n_perm Number of permutations per domain (default 10)
#' @param seed Random seed
#' @return data.frame with domain, n_vars, baseline metric, permuted metric,
#'   importance (drop), and CI
run_conceptual_permutation <- function(sl_fit, outcome_data, oc,
                                       level = "level1", n_perm = 10,
                                       seed = 12345) {
  set.seed(seed)

  mapping <- load_conceptual_domains()

  # Determine which fit to use
  use_binary <- !is.null(sl_fit$bin_fit) && !is.null(oc$binary)
  fit_obj    <- if (use_binary) sl_fit$bin_fit else sl_fit$cont_fit
  if (is.null(fit_obj)) {
    warning("No fitted model available")
    return(data.frame())
  }

  # Get model details
  sl_model   <- fit_obj$sl_fit
  task       <- fit_obj$task
  xvars      <- fit_obj$Xvars
  Y          <- as.double(unclass(fit_obj$res$Y))
  yhat_base  <- as.double(unclass(fit_obj$res$yhat_full))

  # Classify variables
  if (level == "level1") {
    var_domains <- classify_variables(xvars, mapping)
  } else {
    var_domains <- classify_variables_l2(xvars, mapping)
  }

  # Get unique domains (excluding Unknown if too small)
  domains <- sort(unique(var_domains))

  cat(sprintf("[conceptual_perm] %d variables across %d %s domains, %d permutations\n",
              length(xvars), length(domains), level, n_perm))

  # Baseline performance
  if (use_binary) {
    baseline_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat_base, quiet = TRUE))),
      error = function(e) NA_real_
    )
    baseline_brier <- mean((Y - yhat_base)^2, na.rm = TRUE)
    cat(sprintf("  Baseline: AUC = %.3f, Brier = %.4f\n", baseline_auc, baseline_brier))
  } else {
    baseline_r2  <- 1 - sum((Y - yhat_base)^2) / sum((Y - mean(Y))^2)
    baseline_mae <- mean(abs(Y - yhat_base), na.rm = TRUE)
    cat(sprintf("  Baseline: R2 = %.3f, MAE = %.4f\n", baseline_r2, baseline_mae))
  }

  # Get task data as data.table, stripping haven labels to avoid arithmetic errors
  task_data <- data.table::as.data.table(task$data)
  for (col in names(task_data)) {
    if (inherits(task_data[[col]], "haven_labelled")) {
      data.table::set(task_data, j = col, value = as.double(unclass(task_data[[col]])))
    }
  }

  results <- list()

  for (dom in domains) {
    dom_vars <- names(var_domains[var_domains == dom])
    # Only include vars that are actually in the task covariates
    dom_vars <- intersect(dom_vars, xvars)
    if (length(dom_vars) == 0) next

    perm_metrics <- list()

    for (p in seq_len(n_perm)) {
      # Copy task data and shuffle domain columns
      shuffled_data <- data.table::copy(task_data)
      perm_idx <- sample(nrow(shuffled_data))
      for (v in dom_vars) {
        data.table::set(shuffled_data, j = v, value = shuffled_data[[v]][perm_idx])
      }

      # Create new task with shuffled data
      perm_task <- tryCatch(
        sl3::sl3_Task$new(
          data       = shuffled_data,
          covariates = xvars,
          outcome    = "Y",
          id         = "id",
          folds      = task$folds
        ),
        error = function(e) NULL
      )
      if (is.null(perm_task)) next

      # Use predict_fold "validation" to get cross-validated predictions,
      # matching how the baseline AUC was computed. predict() would use
      # the full model on training data, inflating performance.
      yhat_perm <- tryCatch(
        as.numeric(sl_model$predict_fold(perm_task, "validation")),
        error = function(e) {
          # Fall back to regular predict if predict_fold fails
          tryCatch(as.numeric(sl_model$predict(perm_task)),
                   error = function(e2) NULL)
        }
      )
      if (is.null(yhat_perm) || length(yhat_perm) != length(Y)) next

      if (use_binary) {
        perm_auc <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(Y, yhat_perm, quiet = TRUE))),
          error = function(e) NA_real_
        )
        perm_brier <- mean((Y - yhat_perm)^2, na.rm = TRUE)
        perm_metrics[[p]] <- data.frame(auc = perm_auc, brier = perm_brier)
      } else {
        perm_r2  <- 1 - sum((Y - yhat_perm)^2) / sum((Y - mean(Y))^2)
        perm_mae <- mean(abs(Y - yhat_perm), na.rm = TRUE)
        perm_metrics[[p]] <- data.frame(r2 = perm_r2, mae = perm_mae)
      }
    }

    if (length(perm_metrics) == 0) next
    perm_df <- do.call(rbind, perm_metrics)

    if (use_binary) {
      results[[dom]] <- data.frame(
        domain       = dom,
        n_vars       = length(dom_vars),
        n_perm_valid = nrow(perm_df),
        baseline_auc = baseline_auc,
        perm_auc     = mean(perm_df$auc, na.rm = TRUE),
        auc_drop     = baseline_auc - mean(perm_df$auc, na.rm = TRUE),
        auc_drop_lo  = baseline_auc - quantile(perm_df$auc, 0.975, na.rm = TRUE),
        auc_drop_hi  = baseline_auc - quantile(perm_df$auc, 0.025, na.rm = TRUE),
        baseline_brier = baseline_brier,
        perm_brier   = mean(perm_df$brier, na.rm = TRUE),
        brier_increase = mean(perm_df$brier, na.rm = TRUE) - baseline_brier,
        stringsAsFactors = FALSE
      )
    } else {
      results[[dom]] <- data.frame(
        domain       = dom,
        n_vars       = length(dom_vars),
        n_perm_valid = nrow(perm_df),
        baseline_r2  = baseline_r2,
        perm_r2      = mean(perm_df$r2, na.rm = TRUE),
        r2_drop      = baseline_r2 - mean(perm_df$r2, na.rm = TRUE),
        r2_drop_lo   = baseline_r2 - quantile(perm_df$r2, 0.975, na.rm = TRUE),
        r2_drop_hi   = baseline_r2 - quantile(perm_df$r2, 0.025, na.rm = TRUE),
        baseline_mae = baseline_mae,
        perm_mae     = mean(perm_df$mae, na.rm = TRUE),
        mae_increase = mean(perm_df$mae, na.rm = TRUE) - baseline_mae,
        stringsAsFactors = FALSE
      )
    }

    cat(sprintf("  %s: %d vars, %s\n", dom, length(dom_vars),
                if (use_binary) sprintf("AUC drop = %.4f", results[[dom]]$auc_drop)
                else sprintf("R2 drop = %.4f", results[[dom]]$r2_drop)))
  }

  if (length(results) == 0) return(data.frame())
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Run "leave-one-in" conceptual domain permutation importance
#'
#' The complement of run_conceptual_permutation(): for each domain, permute
#' ALL OTHER domains and leave only that domain intact. This measures each
#' domain's standalone predictive contribution (how well it predicts alone),
#' rather than its marginal contribution on top of everything else.
#'
#' Together with "permute-one" (run_conceptual_permutation), these two analyses
#' bracket a domain's importance:
#'   - permute-one: marginal contribution (redundant domains score low)
#'   - leave-one-in: standalone signal (even redundant domains can score high)
#'
#' @param sl_fit Output from fit_sl_models()
#' @param outcome_data Output from build_outcome_dataset()
#' @param oc Outcome config
#' @param level "level1" or "level2"
#' @param n_perm Number of permutations per domain (default 10)
#' @param seed Random seed
#' @return data.frame with domain, n_vars, baseline metric, solo metric, and
#'   recovery (how much of baseline performance is retained by this domain alone)
run_leave_one_in <- function(sl_fit, outcome_data, oc,
                             level = "level1", n_perm = 10,
                             seed = 12345) {
  set.seed(seed)

  mapping <- load_conceptual_domains()

  use_binary <- !is.null(sl_fit$bin_fit) && !is.null(oc$binary)
  fit_obj    <- if (use_binary) sl_fit$bin_fit else sl_fit$cont_fit
  if (is.null(fit_obj)) {
    warning("No fitted model available")
    return(data.frame())
  }

  sl_model   <- fit_obj$sl_fit
  task       <- fit_obj$task
  xvars      <- fit_obj$Xvars
  Y          <- as.double(unclass(fit_obj$res$Y))
  yhat_base  <- as.double(unclass(fit_obj$res$yhat_full))

  if (level == "level1") {
    var_domains <- classify_variables(xvars, mapping)
  } else {
    var_domains <- classify_variables_l2(xvars, mapping)
  }

  domains <- sort(unique(var_domains))

  cat(sprintf("[leave_one_in] %d variables across %d %s domains, %d permutations\n",
              length(xvars), length(domains), level, n_perm))

  # Baseline performance
  if (use_binary) {
    baseline_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat_base, quiet = TRUE))),
      error = function(e) NA_real_
    )
    baseline_brier <- mean((Y - yhat_base)^2, na.rm = TRUE)
    cat(sprintf("  Baseline: AUC = %.3f, Brier = %.4f\n", baseline_auc, baseline_brier))
  } else {
    baseline_r2  <- 1 - sum((Y - yhat_base)^2) / sum((Y - mean(Y))^2)
    baseline_mae <- mean(abs(Y - yhat_base), na.rm = TRUE)
    cat(sprintf("  Baseline: R2 = %.3f, MAE = %.4f\n", baseline_r2, baseline_mae))
  }

  # Get task data, strip haven labels
  task_data <- data.table::as.data.table(task$data)
  for (col in names(task_data)) {
    if (inherits(task_data[[col]], "haven_labelled")) {
      data.table::set(task_data, j = col, value = as.double(unclass(task_data[[col]])))
    }
  }

  results <- list()

  for (dom in domains) {
    # Variables to KEEP intact (this domain only)
    keep_vars <- names(var_domains[var_domains == dom])
    keep_vars <- intersect(keep_vars, xvars)
    if (length(keep_vars) == 0) next

    # Variables to SHUFFLE (everything else)
    shuffle_vars <- setdiff(xvars, keep_vars)
    if (length(shuffle_vars) == 0) {
      # This domain contains ALL variables — solo = baseline
      if (use_binary) {
        results[[dom]] <- data.frame(
          domain = dom, n_vars = length(keep_vars), n_perm_valid = 0,
          baseline_auc = baseline_auc, solo_auc = baseline_auc,
          recovery_auc = 1.0, baseline_brier = baseline_brier,
          solo_brier = baseline_brier, stringsAsFactors = FALSE
        )
      } else {
        results[[dom]] <- data.frame(
          domain = dom, n_vars = length(keep_vars), n_perm_valid = 0,
          baseline_r2 = baseline_r2, solo_r2 = baseline_r2,
          recovery_r2 = 1.0, baseline_mae = baseline_mae,
          solo_mae = baseline_mae, stringsAsFactors = FALSE
        )
      }
      next
    }

    perm_metrics <- list()

    for (p in seq_len(n_perm)) {
      shuffled_data <- data.table::copy(task_data)
      perm_idx <- sample(nrow(shuffled_data))

      # Shuffle ALL variables EXCEPT the domain of interest
      for (v in shuffle_vars) {
        data.table::set(shuffled_data, j = v, value = shuffled_data[[v]][perm_idx])
      }

      perm_task <- tryCatch(
        sl3::sl3_Task$new(
          data       = shuffled_data,
          covariates = xvars,
          outcome    = "Y",
          id         = "id",
          folds      = task$folds
        ),
        error = function(e) NULL
      )
      if (is.null(perm_task)) next

      yhat_perm <- tryCatch(
        as.numeric(sl_model$predict_fold(perm_task, "validation")),
        error = function(e) {
          tryCatch(as.numeric(sl_model$predict(perm_task)),
                   error = function(e2) NULL)
        }
      )
      if (is.null(yhat_perm) || length(yhat_perm) != length(Y)) next

      if (use_binary) {
        perm_auc <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(Y, yhat_perm, quiet = TRUE))),
          error = function(e) NA_real_
        )
        perm_brier <- mean((Y - yhat_perm)^2, na.rm = TRUE)
        perm_metrics[[p]] <- data.frame(auc = perm_auc, brier = perm_brier)
      } else {
        perm_r2  <- 1 - sum((Y - yhat_perm)^2) / sum((Y - mean(Y))^2)
        perm_mae <- mean(abs(Y - yhat_perm), na.rm = TRUE)
        perm_metrics[[p]] <- data.frame(r2 = perm_r2, mae = perm_mae)
      }
    }

    if (length(perm_metrics) == 0) next
    perm_df <- do.call(rbind, perm_metrics)

    if (use_binary) {
      solo_auc <- mean(perm_df$auc, na.rm = TRUE)
      # Recovery: what fraction of baseline AUC does this domain alone achieve?
      # 0.5 is chance-level AUC, so recovery = (solo - 0.5) / (baseline - 0.5)
      recovery <- if (!is.na(baseline_auc) && baseline_auc > 0.5)
        (solo_auc - 0.5) / (baseline_auc - 0.5) else NA_real_

      results[[dom]] <- data.frame(
        domain       = dom,
        n_vars       = length(keep_vars),
        n_perm_valid = nrow(perm_df),
        baseline_auc = baseline_auc,
        solo_auc     = solo_auc,
        recovery_auc = recovery,
        solo_auc_lo  = quantile(perm_df$auc, 0.025, na.rm = TRUE),
        solo_auc_hi  = quantile(perm_df$auc, 0.975, na.rm = TRUE),
        baseline_brier = baseline_brier,
        solo_brier   = mean(perm_df$brier, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    } else {
      solo_r2 <- mean(perm_df$r2, na.rm = TRUE)
      recovery <- if (!is.na(baseline_r2) && baseline_r2 > 0)
        solo_r2 / baseline_r2 else NA_real_

      results[[dom]] <- data.frame(
        domain       = dom,
        n_vars       = length(keep_vars),
        n_perm_valid = nrow(perm_df),
        baseline_r2  = baseline_r2,
        solo_r2      = solo_r2,
        recovery_r2  = recovery,
        solo_r2_lo   = quantile(perm_df$r2, 0.025, na.rm = TRUE),
        solo_r2_hi   = quantile(perm_df$r2, 0.975, na.rm = TRUE),
        baseline_mae = baseline_mae,
        solo_mae     = mean(perm_df$mae, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }

    cat(sprintf("  %s: %d vars kept, %s\n", dom, length(keep_vars),
                if (use_binary) sprintf("solo AUC = %.3f (%.0f%% recovery)",
                                        results[[dom]]$solo_auc,
                                        results[[dom]]$recovery_auc * 100)
                else sprintf("solo R2 = %.3f (%.0f%% recovery)",
                             results[[dom]]$solo_r2,
                             results[[dom]]$recovery_r2 * 100)))
  }

  if (length(results) == 0) return(data.frame())
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Plot forest plot of "leave-one-in" domain importance
#'
#' Shows each domain's solo AUC (or R2) as a dot, with the baseline as a
#' reference line. Domains that retain more of the baseline performance
#' when used alone are more independently informative.
#'
#' @param loi_results Output from run_leave_one_in()
#' @param metric "auc" or "r2"
#' @param title Optional plot title
#' @return ggplot object
plot_leave_one_in_forest <- function(loi_results, metric = NULL, title = NULL) {
  if (nrow(loi_results) == 0) return(NULL)

  if (is.null(metric)) {
    metric <- if ("solo_auc" %in% colnames(loi_results)) "auc" else "r2"
  }

  if (metric == "auc") {
    loi_results$solo  <- loi_results$solo_auc
    loi_results$lo    <- if ("solo_auc_lo" %in% names(loi_results))
      loi_results$solo_auc_lo else loi_results$solo
    loi_results$hi    <- if ("solo_auc_hi" %in% names(loi_results))
      loi_results$solo_auc_hi else loi_results$solo
    baseline_val <- loi_results$baseline_auc[1]
    x_label <- "Solo AUC (domain alone, all others permuted)"
    chance_val <- 0.5
  } else {
    loi_results$solo <- loi_results$solo_r2
    loi_results$lo   <- if ("solo_r2_lo" %in% names(loi_results))
      loi_results$solo_r2_lo else loi_results$solo
    loi_results$hi   <- if ("solo_r2_hi" %in% names(loi_results))
      loi_results$solo_r2_hi else loi_results$solo
    baseline_val <- loi_results$baseline_r2[1]
    x_label <- expression("Solo" ~ R^2 ~ "(domain alone, all others permuted)")
    chance_val <- 0
  }

  loi_results <- loi_results[order(loi_results$solo), ]
  loi_results$label <- sprintf("%s (n=%d)", loi_results$domain, loi_results$n_vars)
  loi_results$label <- factor(loi_results$label, levels = loi_results$label)

  if (is.null(title)) title <- "Domain Standalone Contribution (Leave-One-In)"

  ggplot2::ggplot(loi_results, ggplot2::aes(x = solo, y = label)) +
    ggplot2::geom_point(size = 3, color = "#d7191c") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lo, xmax = hi),
      height = 0.3, color = "#d7191c", linewidth = 0.7
    ) +
    ggplot2::geom_vline(xintercept = baseline_val, linetype = "solid",
                        color = "#2c7bb6", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = chance_val, linetype = "dashed",
                        color = "grey50") +
    ggplot2::annotate("text", x = baseline_val, y = Inf,
                      label = "Full model", hjust = -0.1, vjust = 2,
                      color = "#2c7bb6", size = 3.5) +
    ggplot2::labs(title = title, x = x_label, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Plot forest plot of conceptual domain importance
#'
#' @param perm_results Output from run_conceptual_permutation()
#' @param metric "auc" or "r2" or "brier" — which metric to plot
#' @param title Optional plot title
#' @return ggplot object
plot_conceptual_forest <- function(perm_results, metric = NULL, title = NULL) {
  if (nrow(perm_results) == 0) return(NULL)

  # Auto-detect metric
  if (is.null(metric)) {
    metric <- if ("auc_drop" %in% colnames(perm_results)) "auc" else "r2"
  }

  if (metric == "auc") {
    perm_results$importance <- perm_results$auc_drop
    perm_results$lo         <- perm_results$auc_drop_lo
    perm_results$hi         <- perm_results$auc_drop_hi
    x_label <- "AUC Drop (higher = more important)"
  } else {
    perm_results$importance <- perm_results$r2_drop
    perm_results$lo         <- perm_results$r2_drop_lo
    perm_results$hi         <- perm_results$r2_drop_hi
    x_label <- expression(R^2 ~ "Drop (higher = more important)")
  }

  # Sort by importance and create domain label with var count
  perm_results <- perm_results[order(perm_results$importance), ]
  perm_results$label <- sprintf("%s (n=%d)", perm_results$domain, perm_results$n_vars)
  perm_results$label <- factor(perm_results$label, levels = perm_results$label)

  if (is.null(title)) title <- "Conceptual Domain Importance (Permutation)"

  ggplot2::ggplot(perm_results, ggplot2::aes(x = importance, y = label)) +
    ggplot2::geom_point(size = 3, color = "#2c7bb6") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lo, xmax = hi),
      height = 0.3, color = "#2c7bb6", linewidth = 0.7
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = title,
      x     = x_label,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      panel.grid.minor = ggplot2::element_blank()
    )
}
