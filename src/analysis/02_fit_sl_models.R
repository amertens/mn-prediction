# =============================================================================
# src/analysis/02_fit_sl_models.R
#
# Fits SuperLearner models for each (population × micronutrient) outcome.
#
# Deliverable A – Continuous SL fits (original behaviour, standardised):
#   Uses `slmod`  (loss_squared_error) from src/0-SL-setup.R.
#   Saved to results/models/ with the canonical naming convention.
#
# Deliverable B/C/D – Binary SL fits (direct deficiency probability):
#   Uses `slmod2_bin` (loss_loglik_binomial) from src/0-SL-setup.R.
#   Requires the pre-computed binary deficiency columns in each dataset.
#   Saved to results/models/ with a `res_bin_` prefix.
#
# Both model types use cluster-blocked K-fold CV (folds on gw_cnum clusters)
# so that no cluster is split across folds.
#
# Inputs  : cfg, gw_data_list, Xvars_full (from 01_load_and_construct.R)
#           slmod, slmod2_bin              (from src/0-SL-setup.R)
# Outputs : results/models/res_GW_Gambia_SL_*.rds  (continuous)
#           results/models/res_bin_GW_Gambia_SL_*.rds (binary)
#           sl_results_cont (named list of continuous fit objects in environment)
#           sl_results_bin  (named list of binary   fit objects in environment)
# =============================================================================

cat("\n[02] Fitting SuperLearner models...\n")

# ---- Helper: DHS_SL with cluster-blocked folds ----------------------------
#
# Identical to DHS_SL() in src/0-SL-setup.R, with one change:
#   make_folds(cluster_ids = id_vec, V = folds)
# This ensures that all observations from the same gw_cnum cluster fall in
# the same CV fold, preventing cluster-level leakage.
#
# Parameters match DHS_SL exactly so the two functions are interchangeable.

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
  covars <- colnames(cov)

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

  list(
    sl_fit               = sl_fit,
    res                  = res,
    cv_risk_w_sl_revere  = cv_risk,
    task                 = SL_task,
    Xvars                = covars
  )
}


# ---- Convenience wrapper with error handling --------------------------------
safe_fit <- function(label, d, outcome, population, Xvars, id,
                     folds, sl, outfile) {
  cat(sprintf("  Fitting %-35s ... ", label))
  t0  <- proc.time()
  res <- tryCatch(
    DHS_SL_clustered(
      d          = d,
      Xvars      = Xvars,
      outcome    = outcome,
      population = population,
      id         = id,
      folds      = folds,
      CV         = FALSE,
      prescreen  = TRUE,
      sl         = sl
    ),
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )
  elapsed <- round((proc.time() - t0)["elapsed"], 1)
  if (!is.null(res)) {
    saveRDS(res, file = outfile)
    cat(sprintf("done (%gs)  → %s\n", elapsed, basename(outfile)))
  }
  res
}

# ---- Ensure output directory exists ----------------------------------------
dir.create(cfg$out_models, showWarnings = FALSE, recursive = TRUE)

id    <- cfg$cluster_id
K     <- cfg$K

# ============================================================================
# A. Continuous SL fits  (slmod, loss_squared_error)
# ============================================================================
cat("\n  -- Continuous models (slmod) --\n")

sl_results_cont <- list()

sl_results_cont$child_vitA <- safe_fit(
  label      = "children / Vitamin A (continuous)",
  d          = gw_data_list$child_vitA,
  outcome    = cfg$outcomes$child_vitA$continuous,
  population = cfg$outcomes$child_vitA$population,
  Xvars      = Xvars_full[Xvars_full %in% colnames(gw_data_list$child_vitA)],
  id         = id,
  folds      = K,
  sl         = slmod,
  outfile    = file.path(cfg$out_models, cfg$outcomes$child_vitA$model_file_cont)
)

sl_results_cont$women_vitA <- safe_fit(
  label      = "women / Vitamin A (continuous)",
  d          = gw_data_list$women_vitA,
  outcome    = cfg$outcomes$women_vitA$continuous,
  population = cfg$outcomes$women_vitA$population,
  Xvars      = Xvars_full[Xvars_full %in% colnames(gw_data_list$women_vitA)],
  id         = id,
  folds      = K,
  sl         = slmod,
  outfile    = file.path(cfg$out_models, cfg$outcomes$women_vitA$model_file_cont)
)

sl_results_cont$child_iron <- safe_fit(
  label      = "children / Iron (continuous)",
  d          = gw_data_list$child_iron,
  outcome    = cfg$outcomes$child_iron$continuous,
  population = cfg$outcomes$child_iron$population,
  Xvars      = Xvars_full[Xvars_full %in% colnames(gw_data_list$child_iron)],
  id         = id,
  folds      = K,
  sl         = slmod,
  outfile    = file.path(cfg$out_models, cfg$outcomes$child_iron$model_file_cont)
)

sl_results_cont$women_iron <- safe_fit(
  label      = "women / Iron (continuous)",
  d          = gw_data_list$women_iron,
  outcome    = cfg$outcomes$women_iron$continuous,
  population = cfg$outcomes$women_iron$population,
  Xvars      = Xvars_full[Xvars_full %in% colnames(gw_data_list$women_iron)],
  id         = id,
  folds      = K,
  sl         = slmod,
  outfile    = file.path(cfg$out_models, cfg$outcomes$women_iron$model_file_cont)
)

# ============================================================================
# B. Binary SL fits  (slmod2_bin, loss_loglik_binomial)
#
# Uses the pre-computed binary deficiency columns.  These exist in the dataset
# (gw_cVAD_Thurn, gw_wVAD_Thurn, gw_cIDA_Brinda, gw_wIDA_Brinda) and represent
# the clinical standard definitions (Thurnham-adjusted VAD; Brinda-adjusted IDA).
# Direct binary modeling of deficiency probability is preferred over post-hoc
# thresholding of continuous predictions.
# ============================================================================
cat("\n  -- Binary models (slmod2_bin) --\n")

# Helper: return dataset only if binary outcome column is present & non-trivial
get_bin_data <- function(d, bin_col) {
  if (!bin_col %in% colnames(d)) {
    cat(sprintf("    [skip] binary column '%s' not found in dataset\n", bin_col))
    return(NULL)
  }
  d_bin <- d[!is.na(d[[bin_col]]), ]
  tab   <- table(d_bin[[bin_col]])
  if (length(tab) < 2) {
    cat(sprintf("    [skip] '%s' has only one level\n", bin_col))
    return(NULL)
  }
  d_bin
}

sl_results_bin <- list()

# -- children VitA binary --
d_tmp <- get_bin_data(gw_data_list$child_vitA, cfg$outcomes$child_vitA$binary)
if (!is.null(d_tmp)) {
  sl_results_bin$child_vitA <- safe_fit(
    label      = "children / Vitamin A (binary)",
    d          = d_tmp,
    outcome    = cfg$outcomes$child_vitA$binary,
    population = cfg$outcomes$child_vitA$population,
    Xvars      = Xvars_full[Xvars_full %in% colnames(d_tmp)],
    id         = id,
    folds      = K,
    sl         = slmod2_bin,
    outfile    = file.path(cfg$out_models, cfg$outcomes$child_vitA$model_file_bin)
  )
}

# -- women VitA binary --
d_tmp <- get_bin_data(gw_data_list$women_vitA, cfg$outcomes$women_vitA$binary)
if (!is.null(d_tmp)) {
  sl_results_bin$women_vitA <- safe_fit(
    label      = "women / Vitamin A (binary)",
    d          = d_tmp,
    outcome    = cfg$outcomes$women_vitA$binary,
    population = cfg$outcomes$women_vitA$population,
    Xvars      = Xvars_full[Xvars_full %in% colnames(d_tmp)],
    id         = id,
    folds      = K,
    sl         = slmod2_bin,
    outfile    = file.path(cfg$out_models, cfg$outcomes$women_vitA$model_file_bin)
  )
}

# -- children iron binary --
d_tmp <- get_bin_data(gw_data_list$child_iron, cfg$outcomes$child_iron$binary)
if (!is.null(d_tmp)) {
  sl_results_bin$child_iron <- safe_fit(
    label      = "children / Iron (binary)",
    d          = d_tmp,
    outcome    = cfg$outcomes$child_iron$binary,
    population = cfg$outcomes$child_iron$population,
    Xvars      = Xvars_full[Xvars_full %in% colnames(d_tmp)],
    id         = id,
    folds      = K,
    sl         = slmod2_bin,
    outfile    = file.path(cfg$out_models, cfg$outcomes$child_iron$model_file_bin)
  )
}

# -- women iron binary --
d_tmp <- get_bin_data(gw_data_list$women_iron, cfg$outcomes$women_iron$binary)
if (!is.null(d_tmp)) {
  sl_results_bin$women_iron <- safe_fit(
    label      = "women / Iron (binary)",
    d          = d_tmp,
    outcome    = cfg$outcomes$women_iron$binary,
    population = cfg$outcomes$women_iron$population,
    Xvars      = Xvars_full[Xvars_full %in% colnames(d_tmp)],
    id         = id,
    folds      = K,
    sl         = slmod2_bin,
    outfile    = file.path(cfg$out_models, cfg$outcomes$women_iron$model_file_bin)
  )
}

cat("[02] Done.\n\n")
