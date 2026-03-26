# =============================================================================
# R/domain_ablation.R
#
# Functions for domain ablation (leave-one-domain-out feature importance).
# =============================================================================

#' Fit a reduced model with one domain removed
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @param domain_to_remove Name of domain to drop (e.g., "GEE")
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @param use_binary If TRUE, fit binary SL
#' @return Metrics data.frame or NULL
fit_reduced_domain <- function(outcome_data, cc, oc, domain_to_remove,
                                sl_learners, params, use_binary = TRUE) {

  source(here::here("src", "analysis", "sl_helpers.R"), local = FALSE)

  d           <- outcome_data$data
  Xvars_proxy <- outcome_data$Xvars    # proxy-only (no gw_)
  domain_vars <- outcome_data$domain_vars

  # Remove the target domain columns from the proxy-only set
  drop_cols <- domain_vars[[domain_to_remove]]
  X_reduced <- setdiff(Xvars_proxy, drop_cols)

  if (length(X_reduced) < 2) return(NULL)

  Y_name     <- if (use_binary) oc$binary else oc$continuous
  sl_learner <- if (use_binary) sl_learners$slmod2_bin else sl_learners$slmod

  cat(sprintf("  Ablating %s (%d vars removed, %d remaining)\n",
              domain_to_remove, length(drop_cols), length(X_reduced)))

  fit_r <- tryCatch(
    DHS_SL_clustered(
      d          = d,
      Xvars      = X_reduced,
      outcome    = Y_name,
      population = oc$population,
      id         = cc$cluster_id,
      sl         = sl_learner
    ),
    error = function(e) {
      warning(sprintf("  Ablation %s failed: %s", domain_to_remove, e$message))
      NULL
    }
  )
  if (is.null(fit_r)) return(NULL)

  res  <- fit_r$res
  Y    <- res$Y
  yhat <- res$yhat_full

  metrics <- data.frame(
    domain_removed = domain_to_remove,
    n_removed      = length(drop_cols),
    n_remaining    = length(X_reduced),
    stringsAsFactors = FALSE
  )

  if (use_binary) {
    metrics$auc   <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat, quiet = TRUE))),
      error = function(e) NA_real_
    )
    metrics$brier <- mean((Y - yhat)^2, na.rm = TRUE)
  } else {
    metrics$rmse <- sqrt(mean((Y - yhat)^2, na.rm = TRUE))
    metrics$mae  <- mean(abs(Y - yhat), na.rm = TRUE)
    metrics$r2   <- 1 - sum((Y - yhat)^2) / sum((Y - mean(Y))^2)
  }

  metrics
}


#' Run full domain ablation for one outcome
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param sl_fit Output from fit_sl_models() (for baseline metrics)
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @return data.frame with one row per domain (including "none" baseline)
run_domain_ablation <- function(outcome_data, sl_fit, cc, oc, sl_learners, params) {

  # Make cfg available for sl_helpers
  cfg <- list(
    seed = params$seed, admin1 = cc$admin1_col,
    cluster_id = cc$cluster_id, K = params$K
  )
  assign("cfg", cfg, envir = globalenv())

  # Increase future globals size limit — SL learner objects can be >500 MB
  old_maxsize <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 3 * 1024^3)  # 3 GB
  on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)

  use_binary <- !is.null(sl_fit$bin_fit)
  domains <- names(outcome_data$domain_vars)
  # Only ablate proxy domains (skip GW — not used in primary model)
  domains <- setdiff(domains, "GW")
  # Only ablate domains that have columns in the proxy predictor set
  proxy_vars <- outcome_data$Xvars
  domains <- domains[sapply(domains, function(nm) {
    dom_cols <- outcome_data$domain_vars[[nm]]
    length(intersect(dom_cols, proxy_vars)) > 0
  })]

  n_domains <- length(domains)
  cat(sprintf("\n[ablation] %s — %s | %d domains to ablate\n",
              cc$country, oc$tag, n_domains))

  # Baseline metrics (full model)
  baseline_fit <- if (use_binary) sl_fit$bin_fit else sl_fit$cont_fit
  baseline <- extract_cv_performance(baseline_fit, oc,
                                      if (use_binary) "binary" else "continuous")
  baseline$domain_removed <- "none"

  # ── Parallel domain ablation via future.apply ────────────────────────────
  # Each domain fit is independent. Scale workers to number of domains.
  # Use parallelly::availableCores() which respects constraints set by
  # tar_make_future() (mc.cores=1 inside a worker). If only 1 core is
  # available (i.e., we're inside a targets worker), run sequentially.
  avail <- parallelly::availableCores()
  n_abl_workers <- max(1L, min(n_domains, as.integer(
    Sys.getenv("ABLATION_WORKERS", max(1L, floor(avail / 2)))
  )))

  if (n_abl_workers > 1 && n_domains > 1 && avail > 1) {
    cat(sprintf("  Using %d parallel workers for %d domains\n",
                n_abl_workers, n_domains))

    old_plan <- future::plan()
    future::plan(future::multisession, workers = n_abl_workers)
    on.exit(future::plan(old_plan), add = TRUE)

    reduced_list <- future.apply::future_lapply(domains, function(dom) {
      assign("cfg", cfg, envir = globalenv())
      fit_reduced_domain(outcome_data, cc, oc, dom, sl_learners, params, use_binary)
    }, future.seed = TRUE)
  } else {
    cat("  Running ablation sequentially\n")
    reduced_list <- lapply(domains, function(dom) {
      fit_reduced_domain(outcome_data, cc, oc, dom, sl_learners, params, use_binary)
    })
  }

  all_results <- dplyr::bind_rows(c(list(baseline), Filter(Negate(is.null), reduced_list)))
  all_results$outcome <- oc$tag
  all_results$country <- cc$country

  all_results
}
