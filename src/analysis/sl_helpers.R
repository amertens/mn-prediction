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
                              prescreen = TRUE, sl) {

  X      <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  cov    <- labelled::unlabelled(X, user_na_to_na = TRUE)
  Y      <- d[[outcome]]
  id_vec <- d[[id]]
  dataid <- d$dataid

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
  cov <- cov %>%
    do(ck37r::impute_missing_values(., type = "standard",
                                    add_indicators = TRUE,
                                    prefix = "missing_")$data) %>%
    as.data.frame()

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

  # Recipes-based final cleaning (correlation removal, normalisation)
  auto_recipe <- recipes::recipe(~ ., data = cov) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_numeric(), threshold = 0.9) %>%
    recipes::step_normalize(recipes::all_numeric()) %>%
    recipes::prep()
  cov    <- recipes::bake(auto_recipe, new_data = cov)
  cov    <- data.frame(cov)
  covars <- colnames(cov)   # final covariate names (excludes Y and id)

  dat <- data.table::data.table(Y = Y, id = id_vec, cov)

  # Cluster-blocked folds
  set.seed(cfg$seed)
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

  # NOTE: return covars (not SL_task$column_names).
  # SL_task$column_names includes the internally-created Y and id columns,
  # which are 2 wider than what the learners were trained on.  Returning
  # covars ensures the bootstrap prediction task has the correct dimensions.
  list(
    sl_fit               = sl_fit,
    res                  = res,
    cv_risk_w_sl_revere  = cv_risk,
    task                 = SL_task,
    Xvars                = covars
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
                          seed_base = 12345L) {

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

  # 3. Predict on the ORIGINAL full dataset.
  #    Use task$nodes$covariates — the covariate names the sl3 task was built
  #    with — NOT fit_b$Xvars, which in older session objects may include the
  #    Y and id columns added to the data.table inside DHS_SL_clustered.
  final_covars <- fit_b$task$nodes$covariates

  X_pred <- tryCatch({
    X0   <- d_predict %>%
      dplyr::select(dplyr::any_of(Xvars_b)) %>%
      as.data.frame()
    cov0 <- labelled::unlabelled(X0, user_na_to_na = TRUE)
    cov0 <- cov0 %>%
      do(ck37r::impute_missing_values(., type = "standard",
                                      add_indicators = TRUE,
                                      prefix = "missing_")$data) %>%
      as.data.frame()
    # Align to the columns the model was trained on
    keep <- intersect(final_covars, colnames(cov0))
    if (length(keep) == 0) return(NULL)
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
  deficient_pred <- if (binary_outcome) {
    yhat
  } else {
    as.numeric(apply_threshold(yhat, cutoff, cutoff_dir))
  }

  # 5. Aggregate to Admin1 and national
  admin1_col <- cfg$admin1
  out_df <- data.frame(Admin1 = d_predict[[admin1_col]],
                       deficient_pred = deficient_pred)
  out_df <- out_df[!is.na(out_df$Admin1), ]

  a1_prev <- out_df %>%
    dplyr::group_by(Admin1) %>%
    dplyr::summarise(prev = mean(deficient_pred, na.rm = TRUE), .groups = "drop")

  list(admin1 = a1_prev, national = mean(deficient_pred, na.rm = TRUE))
}
