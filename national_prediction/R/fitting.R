# =============================================================================
# National Prediction Pipeline — Model Fitting
# =============================================================================

#' Fit a national-level SuperLearner model
#'
#' Prepares data (NZV removal, subsetting), then calls fit_sl() for the core
#' sl3-based fitting workflow. This is the main entry point for each outcome.
#'
#' @param d Merged data.frame (outcomes + predictors)
#' @param outcome Character: outcome column name
#' @param Xvars_full Character vector: candidate predictor names
#' @param id Character: grouping variable for CV folds (default "country")
#' @param folds Integer: number of CV folds (default 2)
#' @param CV Logical: use Lrnr_cv wrapper? (default FALSE)
#' @param sl sl3 Lrnr_sl object
#' @return List with sl_fit, res (data.frame), cv_risk, task, Xvars
run_national_sl <- function(d, outcome, Xvars_full, id = "country",
                            folds = 2L, CV = FALSE, sl) {
  df <- d %>%
    dplyr::select(dplyr::all_of(c(outcome, Xvars_full, id))) %>%
    dplyr::mutate(id = factor(.data[[id]])) %>%
    as.data.frame()

  df <- df[!is.na(df[[outcome]]), ]

  # Remove near-zero variance columns
  nzv_idx <- caret::nearZeroVar(df)
  if (length(nzv_idx) > 0) {
    cols_to_drop <- colnames(df)[nzv_idx]
    df <- df[, !(colnames(df) %in% cols_to_drop)]
  }

  Xvars_full <- Xvars_full[Xvars_full %in% colnames(df)]

  res <- try(fit_sl(
    d = df, outcome = outcome, Xvars = Xvars_full,
    id = id, folds = folds, CV = CV, sl = sl
  ))
  res
}


#' Core SuperLearner fitting workflow
#'
#' Handles: drop all-NA columns, drop zero/near-zero variance, impute missing,
#' prescreen (washb_prescreen), recipe-based cleaning, then sl3 task + training.
#' Ported from src/0-SL-setup.R DHS_SL().
#'
#' @param d data.frame with outcome + predictors
#' @param Xvars Character vector: predictor column names
#' @param outcome Character: outcome column name
#' @param id Character: ID column for fold stratification
#' @param folds Integer: number of CV folds
#' @param CV Logical: wrap SL in Lrnr_cv?
#' @param sl sl3 Lrnr_sl object
#' @param prescreen Logical: run washb_prescreen? (default TRUE)
#' @param prescreen_pval Numeric: p-value threshold for prescreening
#' @return List: sl_fit, res, cv_risk_w_sl_revere, task, Xvars
fit_sl <- function(d, Xvars, outcome, id = "country",
                   folds = 2L, CV = FALSE, sl,
                   prescreen = TRUE, prescreen_pval = 0.2) {

  X <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  cov <- labelled::unlabelled(X, user_na_to_na = TRUE)
  Y   <- d[[outcome]]
  id_vec <- d[[id]]

  # Drop columns that are entirely NA
  cov <- cov[, !sapply(cov, function(x) all(is.na(x))), drop = FALSE]

  # Drop columns with no variation
  cov <- cov %>%
    dplyr::select(dplyr::where(~{
      non_na <- .x[!is.na(.x)]
      length(non_na) > 0 && length(unique(non_na)) > 1
    }))

  # Near-zero variance removal
  nzv_cols <- caret::nearZeroVar(cov)
  if (length(nzv_cols) > 0) {
    cov <- cov[, -nzv_cols, drop = FALSE]
  }

  # Impute missing values
  cov <- cov %>%
    do(ck37r::impute_missing_values(
      ., type = "standard", add_indicators = TRUE, prefix = "missing_"
    )$data) %>%
    as.data.frame()

  # Second NZV pass after imputation
  nzv_cols <- caret::nearZeroVar(cov)
  if (length(nzv_cols) > 0) {
    cov <- cov[, -nzv_cols, drop = FALSE]
  }

  # Prescreen: univariate association filter
  if (prescreen) {
    Wvars <- washb::washb_prescreen(
      Y = Y, Ws = cov, family = "gaussian",
      pval = prescreen_pval, print = FALSE
    )
    cov <- cov %>% dplyr::select(dplyr::all_of(Wvars))
  }

  # Recipe-based cleaning
  auto_recipe <- recipes::recipe(~ ., data = cov) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_numeric(), threshold = 0.9) %>%
    recipes::step_normalize(recipes::all_numeric()) %>%
    recipes::prep()

  cov <- recipes::bake(auto_recipe, new_data = cov)
  cov <- data.frame(cov)
  covars <- colnames(cov)

  dat <- data.table::data.table(Y = Y, id = id_vec, cov)

  # Create sl3 task
  set.seed(12345)
  SL_task <- sl3::make_sl3_Task(
    data       = dat,
    covariates = covars,
    outcome    = "Y",
    id         = "id",
    folds      = origami::make_folds(dat, fold_fun = origami::folds_vfold, V = folds)
  )

  # Train
  if (CV) {
    cv_sl <- sl3::make_learner(sl3::Lrnr_cv, sl, full_fit = TRUE)
    sl_fit <- suppressMessages(cv_sl$train(SL_task))
  } else {
    sl_fit <- suppressMessages(sl$train(SL_task))
  }

  # CV risk
  cv_risk_w_sl_revere <- sl_fit$cv_risk(
    eval_fun = sl3::loss_squared_error, get_sl_revere_risk = TRUE
  )

  # Predictions
  yhat_full <- sl_fit$predict_fold(SL_task, "validation")

  res <- data.frame(
    id        = id_vec,
    outcome   = outcome,
    Y         = Y,
    yhat_full = yhat_full
  )

  list(
    sl_fit              = sl_fit,
    res                 = res,
    cv_risk_w_sl_revere = cv_risk_w_sl_revere,
    task                = SL_task,
    Xvars               = SL_task$column_names,
    Y                   = Y,
    yhat_full           = yhat_full
  )
}
