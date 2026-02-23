# =============================================================================
# src/analysis/sl_helpers.R
#
# Defines DHS_SL_clustered — the cluster-blocked SuperLearner wrapper used
# by 02_fit_sl_models.R (fitting) and 04_bootstrap_uncertainty.R (bootstrap).
#
# Sourcing this file updates the function definition in the current session
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
