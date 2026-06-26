# =============================================================================
# src/analysis/sl_helpers.R
#
# Defines two helper functions used across the pipeline:
#
#   DHS_SL_clustered  – cluster-blocked SuperLearner wrapper
#                       (02_fit_sl_models.R and 04_bootstrap_uncertainty.R)
#
#   one_bootstrap     – single bootstrap iteration for CI estimation
#                       (04_bootstrap_uncertainty.R)
#
# Sourcing this file updates both function definitions in the current session
# without re-running any model fitting.
# =============================================================================

DHS_SL_clustered <- function(d, Xvars, outcome = "mod_sev_anemia",
                              population, id = "gw_cnum",
                              folds = 5L, CV = FALSE,
                              prescreen = TRUE, sl,
                              seed = NULL,
                              normalize_method = "zscore") {
  # normalize_method: "zscore" (legacy step_normalize) or "rank"
  #   (step_percentile, outside="both"). Rank/quantile normalisation is robust to
  #   the heavy tails of environmental predictors and transfers better across
  #   countries; validated in sandbox_fe/ (rank >= z-score on the full SL stack).
  #   step_percentile is bake-aware, so the bootstrap/prediction path in
  #   one_bootstrap() reproduces it automatically via the saved recipe. We pass
  #   outside="both" so new-country values beyond the training range clamp to
  #   [0,1] instead of becoming NA (which would break glmnet/ranger prediction).
  normalize_method <- match.arg(normalize_method, c("zscore", "rank"))
  # Resolve seed: explicit arg > cfg global > default
  if (is.null(seed)) {
    seed <- tryCatch(cfg$seed, error = function(e) 12345L)
  }

  # Ensure future globals size limit is large enough for SL learner objects
  options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

  X      <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  cov    <- labelled::unlabelled(X, user_na_to_na = TRUE)
  Y      <- d[[outcome]]
  id_vec <- d[[id]]
  dataid <- if ("dataid" %in% colnames(d)) d$dataid else seq_len(nrow(d))

  # Drop all-NA columns
  cov <- cov[, !sapply(cov, function(x) all(is.na(x))), drop = FALSE]

  # Drop zero/near-zero variance columns
  cov <- cov %>%
    dplyr::select(dplyr::where(~{
      non_na <- .x[!is.na(.x)]
      length(non_na) > 0 && length(unique(non_na)) > 1
    }))
  nzv_idx <- caret::nearZeroVar(cov)
  if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

  # Impute missing values (ck37r)
  # Note: avoid dplyr::do() which is deprecated and fragile with edge cases
  cov <- as.data.frame(
    ck37r::impute_missing_values(cov, type = "standard",
                                 add_indicators = TRUE,
                                 prefix = "missing_")$data
  )

  # Second NZV pass (new indicator columns may be zero-variance)
  nzv_idx <- caret::nearZeroVar(cov)
  if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

  covars <- colnames(cov)

  # Screening (skip for binary; washb_prescreen uses gaussian family)
  if (prescreen) {
    family_screen <- if (length(unique(Y[!is.na(Y)])) == 2) "binomial" else "gaussian"
    Wvars <- washb::washb_prescreen(Y = Y, Ws = cov,
                                    family = family_screen,
                                    pval = 0.2, print = FALSE)
    cov <- cov %>% dplyr::select(dplyr::all_of(Wvars))
  }

  # Save the column names that go INTO the recipe (pre-recipe columns).
  # This is the column set AFTER imputation + NZV + prescreening, but BEFORE
  # the recipe applies step_corr / step_normalize.  Needed by one_bootstrap()
  # to prepare prediction data identically.
  pre_recipe_cols <- colnames(cov)

  # Recipes-based final cleaning (correlation removal, normalisation).
  # Final scaling step is configurable: legacy z-score (step_normalize) or
  # rank/quantile (step_percentile). Both are bake-aware -> prediction path
  # reproduces them via the saved recipe.
  auto_recipe <- recipes::recipe(~ ., data = cov) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_numeric(), threshold = 0.85)
  auto_recipe <- if (normalize_method == "rank") {
    auto_recipe %>% recipes::step_percentile(recipes::all_numeric(), outside = "both")
  } else {
    auto_recipe %>% recipes::step_normalize(recipes::all_numeric())
  }
  auto_recipe <- recipes::prep(auto_recipe)
  cov    <- recipes::bake(auto_recipe, new_data = cov)
  cov    <- data.frame(cov)
  covars <- colnames(cov)   # final covariate names (excludes Y and id)

  dat <- data.table::data.table(Y = Y, id = id_vec, cov)

  # Cluster-blocked folds
  set.seed(seed)
  fold_obj <- origami::make_folds(cluster_ids = id_vec, V = folds)

  SL_task <- sl3::make_sl3_Task(
    data       = dat,
    covariates = covars,
    outcome    = "Y",
    id         = "id",
    folds      = fold_obj
  )

  if (CV) {
    suppressMessages(cv_sl   <- sl3::make_learner(sl3::Lrnr_cv, sl, full_fit = TRUE))
    suppressMessages(sl_fit  <- cv_sl$train(SL_task))
  } else {
    suppressMessages(sl_fit  <- sl$train(SL_task))
  }

  cv_risk <- sl_fit$cv_risk(
    eval_fun = sl3::loss_squared_error,
    get_sl_revere_risk = TRUE
  )

  yhat_full <- sl_fit$predict_fold(SL_task, "validation")

  res <- data.frame(
    dataid     = dataid,
    clusterid  = id_vec,
    outcome    = outcome,
    population = population,
    Y          = Y,
    yhat_full  = yhat_full
  )

  # Return the prepped recipe so that prediction on new data can use

  # recipes::bake(auto_recipe, new_data = ...) for identical transforms.
  list(
    sl_fit               = sl_fit,
    res                  = res,
    cv_risk_w_sl_revere  = cv_risk,
    task                 = SL_task,
    Xvars                = covars,
    auto_recipe          = auto_recipe,
    pre_recipe_cols      = pre_recipe_cols
  )
}


# =============================================================================
# one_bootstrap
# =============================================================================
# Single bootstrap iteration for 04_bootstrap_uncertainty.R.
# Returns a named list: Admin1 prevalence (data.frame) and national prevalence
# (scalar), or NULL if the replicate failed.
#
# Key fix: use fit_b$task$nodes$covariates (the covariate names the SL task
# was actually built with) rather than fit_b$Xvars.  The latter historically
# included the internally-created Y and id columns from the data.table passed
# to sl3, making the prediction matrix 2 columns wider than what glmnet and
# ranger were trained on and causing dimension-mismatch errors.

one_bootstrap <- function(b, d_boot_orig, Xvars_b, outcome_b, population_b,
                          id_col, K, sl_obj,
                          d_predict, cutoff, cutoff_dir,
                          binary_outcome = FALSE,
                          seed_base = 12345L,
                          admin1_col = NULL,
                          normalize_method = "zscore") {
  # normalize_method threads the recipe choice to the bootstrap refit so CIs are
  # computed under the same preprocessing as the point estimate. The prediction
  # step below bakes the saved recipe, so it adapts automatically regardless.
  # Resolve admin1_col: explicit arg > cfg global > default
  if (is.null(admin1_col)) {
    admin1_col <- tryCatch(cfg$admin1, error = function(e) "Admin1")
  }

  set.seed(seed_base + b)

  # 1. Resample clusters with replacement
  clusters      <- unique(d_boot_orig[[id_col]])
  boot_clusters <- sample(clusters, size = length(clusters), replace = TRUE)

  d_b <- do.call(rbind, lapply(boot_clusters, function(cl) {
    d_boot_orig[d_boot_orig[[id_col]] == cl, , drop = FALSE]
  }))

  # 2. Refit SL on the resampled data
  fit_b <- tryCatch(
    DHS_SL_clustered(
      d                = d_b,
      Xvars            = Xvars_b[Xvars_b %in% colnames(d_b)],
      outcome          = outcome_b,
      population       = population_b,
      id               = id_col,
      folds            = K,
      CV               = FALSE,
      prescreen        = TRUE,
      sl               = sl_obj,
      normalize_method = normalize_method
    ),
    error = function(e) NULL
  )

  if (is.null(fit_b)) return(NULL)

  # 3. Predict on the ORIGINAL full dataset.
  #    Must reproduce the EXACT same preprocessing chain as DHS_SL_clustered():
  #      a) Select raw columns → unlabel → drop all-NA → NZV → impute → NZV
  #      b) Subset to prescreened columns (pre_recipe_cols)
  #      c) Bake with the saved recipe (step_corr, then step_normalize OR
  #         step_percentile depending on normalize_method)
  #    This guarantees column count + scale match the trained model.
  final_covars <- fit_b$Xvars

  X_pred <- tryCatch({
    # a) Raw preprocessing (mirrors lines 21-48 of DHS_SL_clustered)
    X0   <- d_predict %>%
      dplyr::select(dplyr::any_of(Xvars_b)) %>%
      as.data.frame()
    cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
    cov0 <- cov0[, !sapply(cov0, function(x) all(is.na(x))), drop = FALSE]

    # NZV removal (same as training)
    cov0 <- cov0 %>%
      dplyr::select(dplyr::where(~{
        non_na <- .x[!is.na(.x)]
        length(non_na) > 0 && length(unique(non_na)) > 1
      }))
    nzv0 <- caret::nearZeroVar(cov0)
    if (length(nzv0) > 0) cov0 <- cov0[, -nzv0, drop = FALSE]

    # Impute (same as training)
    cov0 <- as.data.frame(
      ck37r::impute_missing_values(cov0, type = "standard",
                                   add_indicators = TRUE,
                                   prefix = "missing_")$data
    )
    nzv0 <- caret::nearZeroVar(cov0)
    if (length(nzv0) > 0) cov0 <- cov0[, -nzv0, drop = FALSE]

    # b) Subset to the prescreened columns the recipe expects
    pre_cols <- fit_b$pre_recipe_cols
    # Add any missing columns as 0 (e.g., missing_ indicators not created)
    for (col in setdiff(pre_cols, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, pre_cols, drop = FALSE]

    # c) Apply the saved recipe (step_corr + step_normalize)
    cov0 <- recipes::bake(fit_b$auto_recipe, new_data = cov0)
    cov0 <- data.frame(cov0)

    # Final alignment: ensure exact column match
    for (col in setdiff(final_covars, colnames(cov0))) cov0[[col]] <- 0
    cov0 <- cov0[, final_covars, drop = FALSE]

    data.table::data.table(cov0)
  }, error = function(e) {
    message(sprintf("    Boot %d: prediction prep failed: %s", b, e$message))
    NULL
  })

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
  deficient_pred <- if (binary_outcome) {
    yhat
  } else {
    as.numeric(apply_threshold(yhat, cutoff, cutoff_dir))
  }

  # 5. Aggregate to Admin1 and national
  # admin1_col resolved in function header (explicit arg > cfg > default)
  out_df <- data.frame(Admin1 = d_predict[[admin1_col]],
                       deficient_pred = deficient_pred)
  out_df <- out_df[!is.na(out_df$Admin1), ]

  a1_prev <- out_df %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(prev = mean(deficient_pred, na.rm = TRUE), .groups = "drop")

  list(admin1 = a1_prev, national = mean(deficient_pred, na.rm = TRUE))
}
