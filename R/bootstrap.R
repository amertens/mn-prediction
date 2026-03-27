# =============================================================================
# R/bootstrap.R
#
# Functions for cluster-level bootstrap uncertainty estimation.
# =============================================================================

#' Run one bootstrap replicate: resample clusters, refit, predict, aggregate
#'
#' @param b Replicate index (used for progress only)
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learner sl3 learner object (slmod or slmod2_bin)
#' @param params Pipeline parameters
#' @param use_binary If TRUE, fit binary model; otherwise continuous
#' @return list(admin1 = data.frame, national = scalar) or NULL on failure
one_bootstrap_rep <- function(b, outcome_data, cc, oc, sl_learner, params,
                              use_binary = TRUE) {

  # Source helpers (need DHS_SL_clustered)
  source(here::here("src", "analysis", "sl_helpers.R"), local = TRUE)

  d         <- outcome_data$data
  Xvars     <- outcome_data$Xvars      # proxy-only (no gw_)
  cluster_col <- cc$cluster_id
  admin1_col  <- cc$admin1_col

  Y_name <- if (use_binary && !is.null(oc$binary)) oc$binary else oc$continuous

  # Use replicate-specific seed so each bootstrap sample is different
  # even if DHS_SL_clustered() resets set.seed() internally.
  set.seed(params$seed + b * 1000L)

  # Resample clusters with replacement
  clusters <- unique(d[[cluster_col]])
  boot_clusters <- sample(clusters, replace = TRUE)

  # Build bootstrap dataset (with duplicates for resampled clusters)
  boot_list <- lapply(seq_along(boot_clusters), function(i) {
    rows <- d[d[[cluster_col]] == boot_clusters[i], ]
    rows[[cluster_col]] <- paste0(boot_clusters[i], "_b", i)
    rows
  })
  d_boot <- dplyr::bind_rows(boot_list)

  # Fit SL on bootstrap sample
  fit_b <- tryCatch(
    DHS_SL_clustered(
      d          = d_boot,
      Xvars      = Xvars,
      outcome    = Y_name,
      population = oc$population,
      id         = cluster_col,
      sl         = sl_learner
    ),
    error = function(e) {
      message(sprintf("    Boot %d: SL fit failed: %s", b, e$message))
      NULL
    }
  )
  if (is.null(fit_b)) return(NULL)

  # Predict on original data using the saved recipe for consistent preprocessing.
  # The bootstrap model's DHS_SL_clustered() returns auto_recipe and pre_recipe_cols
  # which define the exact transforms applied during training. We must apply the
  # same chain to the prediction data to ensure column count + scale match.
  pred <- tryCatch({
    # Step 1: select raw predictor columns, unlabel, basic cleaning
    X0   <- d %>% dplyr::select(dplyr::any_of(Xvars)) %>% as.data.frame()
    cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
    cov0 <- cov0[, !sapply(cov0, function(x) all(is.na(x))), drop = FALSE]

    # NZV removal
    cov0 <- cov0 %>%
      dplyr::select(dplyr::where(~{
        non_na <- .x[!is.na(.x)]
        length(non_na) > 0 && length(unique(non_na)) > 1
      }))
    nzv0 <- caret::nearZeroVar(cov0)
    if (length(nzv0) > 0) cov0 <- cov0[, -nzv0, drop = FALSE]

    # Impute
    cov0 <- as.data.frame(
      ck37r::impute_missing_values(cov0, type = "standard",
                                   add_indicators = TRUE,
                                   prefix = "missing_")$data
    )
    nzv0 <- caret::nearZeroVar(cov0)
    if (length(nzv0) > 0) cov0 <- cov0[, -nzv0, drop = FALSE]

    # Align to the pre-recipe columns the bootstrap model expects
    pre_cols <- fit_b$pre_recipe_cols
    for (col in setdiff(pre_cols, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, pre_cols, drop = FALSE]

    # Apply the saved recipe (step_corr + step_normalize)
    cov0 <- recipes::bake(fit_b$auto_recipe, new_data = cov0)
    cov0 <- data.frame(cov0)

    # Final alignment to the exact post-recipe columns
    final_covars <- fit_b$Xvars
    for (col in setdiff(final_covars, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, final_covars, drop = FALSE]

    # Build sl3 task and predict
    pred_task <- sl3::sl3_Task$new(
      data       = data.table::data.table(cov0),
      covariates = final_covars,
      outcome    = NULL
    )
    as.numeric(fit_b$sl_fit$predict(pred_task))
  }, error = function(e) {
    message(sprintf("    Boot %d: prediction failed: %s", b, e$message))
    NULL
  })

  if (is.null(pred) || length(pred) != nrow(d)) return(NULL)

  # Convert to prevalence
  if (!use_binary && !is.null(oc$cutoff)) {
    prev <- apply_threshold(pred, oc$cutoff, oc$cutoff_dir)
  } else {
    prev <- pred
  }

  # Aggregate to Admin1
  admin1 <- data.frame(
    Admin1 = d[[admin1_col]],
    prev   = prev,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(prev = mean(prev, na.rm = TRUE), .groups = "drop")

  list(
    admin1   = admin1,
    national = mean(prev, na.rm = TRUE)
  )
}


#' Run the full bootstrap and summarise CIs
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learners Output from setup_sl_learners()
#' @param params Pipeline parameters
#' @return list with admin1_ci (data.frame) and national_ci (data.frame)
run_bootstrap_ci <- function(outcome_data, cc, oc, sl_learners, params) {

  B <- params$B_boot
  use_binary <- !is.null(oc$binary) && oc$binary %in% colnames(outcome_data$data)
  sl_learner <- if (use_binary) sl_learners$slmod2_bin else sl_learners$slmod

  # Make cfg available for sl_helpers (needed in one_bootstrap_rep -> DHS_SL_clustered)
  cfg <- list(
    seed = params$seed, admin1 = cc$admin1_col,
    cluster_id = cc$cluster_id, K = params$K
  )
  assign("cfg", cfg, envir = globalenv())

  # Increase future globals size limit — SL learner objects can be >500 MB
  old_maxsize <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 3 * 1024^3)  # 3 GB
  on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)

  cat(sprintf("\n[bootstrap] %s — %s | B = %d\n", cc$country, oc$tag, B))

  # ── Parallel bootstrap via future.apply ──────────────────────────────────
  # Each replicate is independent: resample clusters, refit SL, predict.
  # Use parallelly::availableCores() which respects constraints set by
  # tar_make_future() (mc.cores=1 inside a worker). If only 1 core is
  # available (i.e., we're inside a targets worker), run sequentially.
  avail <- parallelly::availableCores()
  n_boot_workers <- max(1L, min(B, as.integer(
    Sys.getenv("BOOT_WORKERS", max(1L, floor(avail / 2)))
  )))

  if (n_boot_workers > 1 && B > 1 && avail > 1) {
    cat(sprintf("  Using %d parallel workers for %d replicates\n",
                n_boot_workers, B))

    # Save and restore the outer future plan (targets may set its own)
    old_plan <- future::plan()
    future::plan(future::multisession, workers = n_boot_workers)
    on.exit(future::plan(old_plan), add = TRUE)

    set.seed(params$seed)
    results <- future.apply::future_lapply(seq_len(B), function(b) {
      assign("cfg", cfg, envir = globalenv())
      one_bootstrap_rep(b, outcome_data, cc, oc, sl_learner, params, use_binary)
    }, future.seed = TRUE)
  } else {
    cat("  Running bootstrap sequentially\n")
    set.seed(params$seed)
    results <- lapply(seq_len(B), function(b) {
      if (b %% 10 == 0 || b == 1) cat(sprintf("  Replicate %d/%d\n", b, B))
      one_bootstrap_rep(b, outcome_data, cc, oc, sl_learner, params, use_binary)
    })
  }

  valid <- Filter(Negate(is.null), results)
  cat(sprintf("  Valid replicates: %d/%d\n", length(valid), B))

  if (length(valid) < 3) {
    warning("Too few valid bootstrap replicates for CIs")
    return(list(admin1_ci = NULL, national_ci = NULL))
  }

  # Admin1 CIs
  all_admin1 <- dplyr::bind_rows(lapply(valid, `[[`, "admin1"))
  admin1_ci <- all_admin1 %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n_boot    = dplyr::n(),
      boot_mean = mean(prev, na.rm = TRUE),
      ci_lo     = quantile(prev, 0.025, na.rm = TRUE),
      ci_hi     = quantile(prev, 0.975, na.rm = TRUE),
      ci_width  = quantile(prev, 0.975, na.rm = TRUE) - quantile(prev, 0.025, na.rm = TRUE),
      .groups   = "drop"
    )

  # National CI
  nat_vec <- sapply(valid, `[[`, "national")
  national_ci <- data.frame(
    outcome   = oc$tag,
    country   = cc$country,
    n_reps    = length(nat_vec),
    boot_mean = mean(nat_vec, na.rm = TRUE),
    ci_lo     = quantile(nat_vec, 0.025, na.rm = TRUE),
    ci_hi     = quantile(nat_vec, 0.975, na.rm = TRUE)
  )

  list(admin1_ci = admin1_ci, national_ci = national_ci)
}
