# =============================================================================
# src/analysis/04_bootstrap_uncertainty.R
#
# Cluster-level bootstrap to produce uncertainty intervals (95% CIs) for
# Admin1 and national deficiency prevalence estimates.
#
# Algorithm:
#   For each of B bootstrap replicates:
#     1. Sample gw_cnum clusters WITH REPLACEMENT (cluster bootstrap)
#     2. Reconstruct the resampled individual-level dataset
#     3. Refit the SL model (binary if available, else continuous+threshold)
#     4. Predict deficiency for every individual in the ORIGINAL data
#        (out-of-bag in bootstrap sense; but we predict on full original data
#        for stable Admin1 aggregation)
#     5. Aggregate to Admin1 and national prevalence
#   Summarise: mean, 2.5th and 97.5th percentiles across B replicates.
#
# Parallelisation: future::plan(multicore) as established by the runner script.
# Memory safety  : future.globals.maxSize set to 5 GB by runner script.
#
# Outputs:
#   results/tables/gambia_admin1_ci_{tag}.csv
#   results/tables/gambia_national_ci_{tag}.csv
#
# Inputs  : cfg, gw_data_list, sl_results_cont, sl_results_bin,
#           DHS_SL_clustered (defined in 02_fit_sl_models.R)
# =============================================================================

#temp for speed:
cfg$B=10
cfg$K=2

cat("\n[04] Bootstrap uncertainty...\n")
cat(sprintf("  B = %d replicates, K = %d folds, seed = %d\n",
            cfg$B, cfg$K, cfg$seed))

library(future.apply)

dir.create(cfg$out_tables, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: one bootstrap iteration ----------------------------------------
# Returns a named list: Admin1 prevalence vector and national prevalence scalar.

one_bootstrap <- function(b, d_boot_orig, Xvars_b, outcome_b, population_b,
                          id_col, K, sl_obj,
                          d_predict, cutoff, cutoff_dir,
                          binary_outcome = FALSE,
                          seed_base = 12345L) {

  set.seed(seed_base + b)

  # 1. Resample clusters with replacement
  clusters     <- unique(d_boot_orig[[id_col]])
  boot_clusters <- sample(clusters, size = length(clusters), replace = TRUE)

  # Build resampled dataset (rows may be duplicated)
  d_b <- do.call(rbind, lapply(boot_clusters, function(cl) {
    d_boot_orig[d_boot_orig[[id_col]] == cl, , drop = FALSE]
  }))

  # 2. Refit SL on the resampled data
  fit_b <- tryCatch(
    DHS_SL_clustered(
      d          = d_b,
      Xvars      = Xvars_b[Xvars_b %in% colnames(d_b)],
      outcome    = outcome_b,
      population = population_b,
      id         = id_col,
      folds      = K,
      CV         = FALSE,
      prescreen  = TRUE,
      sl         = sl_obj
    ),
    error = function(e) NULL
  )

  if (is.null(fit_b)) return(NULL)

  # 3. Predict on the ORIGINAL full dataset
  #    We build a minimal task from the original covariates using the
  #    column names that survived preprocessing in the bootstrap fit.
  final_covars <- fit_b$Xvars   # preprocessed column names

  # The original dataset may not have these exactly (e.g. imputed indicators);
  # we re-preprocess the original data through the SAME recipe path.
  # For simplicity, we re-run the minimal preprocessing on d_predict and
  # predict using the bootstrap-fitted model.

  X_pred <- tryCatch({
    X0  <- d_predict %>%
      dplyr::select(dplyr::any_of(Xvars_b)) %>%
      as.data.frame()
    cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
    # Impute
    cov0 <- cov0 %>%
      do(ck37r::impute_missing_values(., type = "standard",
                                      add_indicators = TRUE,
                                      prefix = "missing_")$data) %>%
      as.data.frame()
    # Keep only columns the model was trained on
    keep <- intersect(final_covars, colnames(cov0))
    if (length(keep) == 0) return(NULL)
    # Add any missing columns as 0
    for (col in setdiff(final_covars, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, final_covars, drop = FALSE]
    data.table::data.table(cov0)
  }, error = function(e) NULL)

  if (is.null(X_pred)) return(NULL)

  pred_task <- tryCatch(
    sl3::sl3_Task$new(
      data       = X_pred,
      covariates = final_covars,
      outcome    = NULL
    ),
    error = function(e) NULL
  )

  if (is.null(pred_task)) return(NULL)

  yhat <- tryCatch(
    as.numeric(fit_b$sl_fit$predict(pred_task)),
    error = function(e) NULL
  )

  if (is.null(yhat) || length(yhat) != nrow(d_predict)) return(NULL)

  # 4. Convert to deficiency indicator / probability
  if (binary_outcome) {
    # binary model: yhat = predicted probability
    deficient_pred <- yhat
  } else {
    # continuous model: threshold at cutoff
    deficient_pred <- as.numeric(apply_threshold(yhat, cutoff, cutoff_dir))
  }

  # 5. Aggregate
  admin1_col <- cfg$admin1
  out_df <- data.frame(
    Admin1         = d_predict[[admin1_col]],
    deficient_pred = deficient_pred
  )
  out_df <- out_df[!is.na(out_df$Admin1), ]

  a1_prev <- out_df %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(prev = mean(deficient_pred, na.rm = TRUE), .groups = "drop")

  nat_prev <- mean(deficient_pred, na.rm = TRUE)

  list(admin1 = a1_prev, national = nat_prev)
}


# ---- Main loop over outcomes -----------------------------------------------

for (tag in names(cfg$outcomes)) {

  outcome_cfg <- cfg$outcomes[[tag]]
  d_orig      <- gw_data_list[[tag]]
  res_bin     <- sl_results_bin[[tag]]
  res_cont    <- sl_results_cont[[tag]]

  if (is.null(res_cont) && is.null(res_bin)) {
    cat(sprintf("  [skip] no model for %s\n", tag))
    next
  }

  cat(sprintf("\n  Bootstrap for: %s\n", tag))

  # Decide which model type to use for bootstrap refitting
  use_binary   <- !is.null(res_bin)
  outcome_col  <- if (use_binary) outcome_cfg$binary else outcome_cfg$continuous
  sl_obj       <- if (use_binary) slmod2_bin else slmod
  model_type   <- if (use_binary) "binary_prob" else "continuous_threshold"

  cat(sprintf("    Model type: %s  |  outcome: %s\n", model_type, outcome_col))

  if (!outcome_col %in% colnames(d_orig)) {
    cat(sprintf("    [skip] column '%s' not found in dataset\n", outcome_col))
    next
  }

  d_fit <- d_orig[!is.na(d_orig[[outcome_col]]), ]

  Xvars_b <- Xvars_full[Xvars_full %in% colnames(d_fit)]

  # Run B bootstrap replicates in parallel
  boot_results <- future.apply::future_lapply(
    seq_len(cfg$B),
    FUN       = one_bootstrap,
    d_boot_orig    = d_fit,
    Xvars_b        = Xvars_b,
    outcome_b      = outcome_col,
    population_b   = outcome_cfg$population,
    id_col         = cfg$cluster_id,
    K              = cfg$K,
    sl_obj         = sl_obj,
    d_predict      = d_orig,
    cutoff         = outcome_cfg$cutoff,
    cutoff_dir     = outcome_cfg$cutoff_dir,
    binary_outcome = use_binary,
    seed_base      = cfg$seed,
    future.seed    = TRUE,
    future.globals = TRUE
  )

  # Filter out NULL results
  boot_results <- Filter(Negate(is.null), boot_results)
  n_valid      <- length(boot_results)
  cat(sprintf("    Valid replicates: %d / %d\n", n_valid, cfg$B))

  if (n_valid < 10) {
    cat("    [warn] Too few valid replicates; skipping CI output.\n")
    next
  }

  # ---- Admin1 CIs -----------------------------------------------------------
  # Collect prevalence per Admin1 across replicates
  all_a1_wide <- do.call(rbind, lapply(seq_along(boot_results), function(i) {
    df_i <- boot_results[[i]]$admin1
    df_i$rep <- i
    df_i
  }))

  admin1_ci <- all_a1_wide %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(
      n_reps     = dplyr::n(),
      boot_mean  = mean(prev,            na.rm = TRUE),
      ci_lo      = quantile(prev, 0.025, na.rm = TRUE),
      ci_hi      = quantile(prev, 0.975, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(outcome = tag, model_type = model_type)

  # Merge observed prevalence from the main results
  obs_tbl <- admin1_prev_list[[tag]]
  if (!is.null(obs_tbl)) {
    admin1_ci <- admin1_ci %>%
      dplyr::left_join(
        obs_tbl %>% dplyr::select(Admin1, obs_prevalence, pred_prevalence),
        by = "Admin1"
      )
  }

  a1_csv <- file.path(
    cfg$out_tables,
    sprintf("gambia_admin1_ci_%s.csv", outcome_cfg$table_tag)
  )
  readr::write_csv(admin1_ci, a1_csv)
  cat(sprintf("    Saved: %s\n", basename(a1_csv)))

  # ---- National CIs ---------------------------------------------------------
  nat_prev_vec <- sapply(boot_results, function(x) x$national)

  national_ci <- data.frame(
    outcome    = tag,
    model_type = model_type,
    n_reps     = length(nat_prev_vec),
    boot_mean  = mean(nat_prev_vec,            na.rm = TRUE),
    ci_lo      = quantile(nat_prev_vec, 0.025, na.rm = TRUE),
    ci_hi      = quantile(nat_prev_vec, 0.975, na.rm = TRUE)
  )

  # Observed national prevalence (from admin1_prev_list if available)
  obs_nat <- if (!is.null(admin1_prev_list[[tag]])) {
    obs_src <- d_orig[!is.na(d_orig[[cfg$admin1]]), ]
    bin_col <- outcome_cfg$binary
    if (!is.null(bin_col) && bin_col %in% colnames(obs_src)) {
      mean(obs_src[[bin_col]], na.rm = TRUE)
    } else {
      mean(apply_threshold(obs_src[[outcome_cfg$continuous]],
                           outcome_cfg$cutoff, outcome_cfg$cutoff_dir),
           na.rm = TRUE)
    }
  } else NA_real_

  national_ci$obs_prevalence <- obs_nat

  nat_csv <- file.path(
    cfg$out_tables,
    sprintf("gambia_national_ci_%s.csv", outcome_cfg$table_tag)
  )
  readr::write_csv(national_ci, nat_csv)
  cat(sprintf("    Saved: %s\n", basename(nat_csv)))

  # Print summary
  cat(sprintf(
    "    National prevalence: %.1f%% [%.1f%%, %.1f%%]\n",
    national_ci$boot_mean * 100,
    national_ci$ci_lo * 100,
    national_ci$ci_hi * 100
  ))
}

cat("[04] Done.\n\n")
