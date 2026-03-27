# =============================================================================
# explore_sl_library.R
#
# Standalone script to explore a wide range of SuperLearner algorithms
# on a single country×outcome combination using TWO implementations:
#   1. sl3 (tlverse) — the current pipeline's framework
#   2. mlr3superlearner — modern mlr3-based implementation
#
# Compares learner rankings and ensemble performance across frameworks
# to check if implementation matters.
#
# Usage:
#   source("scripts/explore_sl_library.R")
#
# Modify COUNTRY, OUTCOME, RUN_SL3, RUN_MLR3 below.
# =============================================================================

library(here)
library(dplyr)
library(data.table)

# ── Configuration ─────────────────────────────────────────────────────────────
COUNTRY  <- "ghana"         # lowercase: gambia, ghana, sierraleone, malawi
OUTCOME  <- "child_iron"    # child_vitA, women_vitA, child_iron, women_iron
USE_BIN  <- TRUE            # TRUE = binary outcome, FALSE = continuous
K_FOLDS  <- 5               # number of CV folds
SEED     <- 12345

RUN_SL3  <- TRUE            # run sl3 library exploration
RUN_MLR3 <- TRUE            # run mlr3superlearner comparison

set.seed(SEED)

# ── Output directories ────────────────────────────────────────────────────────
out_dir   <- here("results", "sl_exploration")
model_dir <- here("results", "sl_exploration", "models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data from targets ────────────────────────────────────────────────────
cat("Loading data from targets store...\n")
store <- here("_targets")
if (!dir.exists(store)) store <- here("_targets_fast_backup")

outcome_data <- targets::tar_read_raw(
  paste0("outcome_data_", COUNTRY, "_", OUTCOME), store = store
)

d     <- outcome_data$data
Xvars <- outcome_data$Xvars

# Source helpers
source(here("src", "analysis", "sl_helpers.R"), local = FALSE)
source(here("src", "0-functions.R"), local = FALSE)

# Get country config
source(here("R", "config.R"))
all_configs <- get_country_configs()
cc_name <- names(all_configs)[tolower(names(all_configs)) == COUNTRY]
cc <- all_configs[[cc_name]]
oc <- cc$outcomes[[OUTCOME]]

cat(sprintf("\n=== Exploring SL Library ===\n"))
cat(sprintf("Country: %s | Outcome: %s | n = %d | p = %d\n",
            COUNTRY, OUTCOME, nrow(d), length(Xvars)))

# Strip haven labels from ALL columns to avoid <haven_labelled> arithmetic errors
for (col in colnames(d)) {
  if (inherits(d[[col]], "haven_labelled")) {
    d[[col]] <- as.double(unclass(d[[col]]))
  }
}

# Determine outcome column
if (USE_BIN && !is.null(oc$binary) && oc$binary %in% colnames(d)) {
  outcome_col <- oc$binary
  task_type   <- "binomial"
  cat(sprintf("Using binary outcome: %s (prevalence: %.1f%%)\n",
              outcome_col, 100 * mean(d[[outcome_col]], na.rm = TRUE)))
} else {
  outcome_col <- oc$continuous
  task_type   <- "continuous"
  cat(sprintf("Using continuous outcome: %s\n", outcome_col))
}

# ── Preprocess (mirrors DHS_SL_clustered) ─────────────────────────────────────
cat("\nPreprocessing...\n")
X <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
cov <- labelled::unlabelled(X, user_na_to_na = TRUE)
Y <- as.numeric(d[[outcome_col]])
id_vec <- d[[cc$cluster_id]]

# Remove all-NA and constant columns
cov <- cov[, !sapply(cov, function(x) all(is.na(x))), drop = FALSE]
cov <- cov %>% dplyr::select(dplyr::where(~{
  non_na <- .x[!is.na(.x)]
  length(non_na) > 0 && length(unique(non_na)) > 1
}))

# NZV
nzv_idx <- caret::nearZeroVar(cov)
if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

# Impute
cov <- as.data.frame(
  ck37r::impute_missing_values(cov, type = "standard",
                                add_indicators = TRUE, prefix = "missing_")$data
)
nzv_idx <- caret::nearZeroVar(cov)
if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

# Prescreen
family_screen <- if (length(unique(Y[!is.na(Y)])) == 2) "binomial" else "gaussian"
Wvars <- washb::washb_prescreen(Y = Y, Ws = cov, family = family_screen,
                                 pval = 0.2, print = FALSE)
cov <- cov %>% dplyr::select(dplyr::all_of(Wvars))

# Recipe
auto_recipe <- recipes::recipe(~ ., data = cov) %>%
  recipes::step_zv(recipes::all_predictors()) %>%
  recipes::step_nzv(recipes::all_predictors()) %>%
  recipes::step_corr(recipes::all_numeric(), threshold = 0.85) %>%
  recipes::step_normalize(recipes::all_numeric()) %>%
  recipes::prep()
cov <- data.frame(recipes::bake(auto_recipe, new_data = cov))
covars <- colnames(cov)

cat(sprintf("After preprocessing: n = %d, p = %d\n", nrow(cov), length(covars)))

# Combined data frame for both frameworks
analysis_df <- data.frame(Y = Y, cov)

# Strip any remaining haven labels
for (col in colnames(analysis_df)) {
  if (inherits(analysis_df[[col]], "haven_labelled")) {
    analysis_df[[col]] <- as.double(unclass(analysis_df[[col]]))
  }
}

# #############################################################################
#
#                          PART 1: sl3 EXPLORATION
#
# #############################################################################

sl3_results <- NULL

if (RUN_SL3) {

  library(sl3)
  library(origami)

  cat("\n\n")
  cat("###################################################################\n")
  cat("#                      sl3 EXPLORATION                           #\n")
  cat("###################################################################\n")

  # Loss function
  loss_fn   <- if (task_type == "binomial") sl3::loss_loglik_binomial else sl3::loss_squared_error
  loss_name <- if (task_type == "binomial") "loglik_binomial" else "squared_error"

  # ── Build sl3 task ──────────────────────────────────────────────────────
  dat <- data.table(Y = Y, id = id_vec, cov)
  fold_obj <- origami::make_folds(cluster_ids = id_vec, V = K_FOLDS)
  task <- sl3::make_sl3_Task(
    data = dat, covariates = covars, outcome = "Y",
    id = "id", folds = fold_obj
  )

  # ── Build learner library ───────────────────────────────────────────────
  cat("\nBuilding sl3 learner library...\n")

  learners <- list()

  # 1. Baseline
  learners$mean <- Lrnr_mean$new()

  # 2. Regularized linear
  learners$lasso       <- Lrnr_glmnet$new(alpha = 1)
  learners$enet_025    <- Lrnr_glmnet$new(alpha = 0.25)
  learners$enet_050    <- Lrnr_glmnet$new(alpha = 0.50)
  learners$enet_075    <- Lrnr_glmnet$new(alpha = 0.75)
  learners$ridge       <- Lrnr_glmnet$new(alpha = 0)

  # 3. GLM
  learners$glm_fast <- tryCatch(Lrnr_glm_fast$new(), error = function(e) NULL)

  # 4. Earth / MARS
  learners$earth_d1 <- tryCatch(Lrnr_earth$new(degree = 1), error = function(e) NULL)
  learners$earth_d2 <- tryCatch(Lrnr_earth$new(degree = 2), error = function(e) NULL)

  # 5. GAM
  learners$gam <- tryCatch(Lrnr_gam$new(), error = function(e) NULL)

  # 6. Polspline
  learners$polspline <- tryCatch(Lrnr_polspline$new(), error = function(e) NULL)

  # 7. Random forests (multiple configs)
  learners$ranger_500_mn5   <- Lrnr_ranger$new(num.trees = 500, min.node.size = 5, importance = "none")
  learners$ranger_500_mn10  <- Lrnr_ranger$new(num.trees = 500, min.node.size = 10, importance = "none")
  learners$ranger_1000_mn10 <- Lrnr_ranger$new(num.trees = 1000, min.node.size = 10, importance = "none")
  learners$ranger_1000_mn20 <- Lrnr_ranger$new(num.trees = 1000, min.node.size = 20, importance = "none")
  learners$ranger_mtry_sqrt <- Lrnr_ranger$new(
    num.trees = 500, min.node.size = 10,
    mtry = max(floor(sqrt(length(covars))), 1), importance = "none"
  )
  learners$ranger_mtry_log <- Lrnr_ranger$new(
    num.trees = 500, min.node.size = 10,
    mtry = max(floor(log2(length(covars))), 1), importance = "none"
  )

  # 8. XGBoost (multiple configs)
  learners$xgb_conservative <- Lrnr_xgboost$new(
    max_depth = 3, eta = 0.05, nrounds = 300,
    min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
    eval_metric = "logloss"
  )
  learners$xgb_deeper <- Lrnr_xgboost$new(
    max_depth = 5, eta = 0.05, nrounds = 300,
    min_child_weight = 10, subsample = 0.8, colsample_bytree = 0.5,
    eval_metric = "logloss"
  )
  learners$xgb_shallow <- Lrnr_xgboost$new(
    max_depth = 2, eta = 0.1, nrounds = 200,
    min_child_weight = 30, subsample = 0.7, colsample_bytree = 0.4,
    eval_metric = "logloss"
  )
  learners$xgb_fast <- Lrnr_xgboost$new(
    max_depth = 4, eta = 0.1, nrounds = 150,
    min_child_weight = 15, subsample = 0.8, colsample_bytree = 0.6,
    eval_metric = "logloss"
  )
  learners$xgb_deep_reg <- Lrnr_xgboost$new(
    max_depth = 6, eta = 0.03, nrounds = 500,
    min_child_weight = 20, subsample = 0.7, colsample_bytree = 0.4,
    lambda = 1, alpha = 0.5, eval_metric = "logloss"
  )

  # 9. BART
  learners$dbarts <- tryCatch(
    Lrnr_dbarts$new(ndpost = 200, ntree = 100), error = function(e) NULL
  )

  # 10. Neural nets — need MaxNWts scaled to p, and prescreening for large p
  #     nnet fails when size * p > MaxNWts (default 1000). With p=363, size=5 needs 1826.
  nnet_maxwts <- max((length(covars) + 1) * 25 + 25, 5000)
  learners$nnet_5  <- tryCatch(Lrnr_nnet$new(size = 5, decay = 0.1, maxit = 200, linout = !USE_BIN, MaxNWts = nnet_maxwts), error = function(e) NULL)
  learners$nnet_10 <- tryCatch(Lrnr_nnet$new(size = 10, decay = 0.01, maxit = 200, linout = !USE_BIN, MaxNWts = nnet_maxwts), error = function(e) NULL)
  learners$nnet_20 <- tryCatch(Lrnr_nnet$new(size = 20, decay = 0.001, maxit = 300, linout = !USE_BIN, MaxNWts = nnet_maxwts), error = function(e) NULL)

  # 11. SVM
  learners$svm_radial <- tryCatch(Lrnr_svm$new(kernel = "radial", cost = 1), error = function(e) NULL)
  learners$svm_linear <- tryCatch(Lrnr_svm$new(kernel = "linear", cost = 1), error = function(e) NULL)

  # 12. RPART — sl3's Lrnr_rpart has a bug with binary outcomes (expects factor).
  #     Only include for continuous tasks.
  if (task_type != "binomial") {
    learners$rpart      <- tryCatch(Lrnr_rpart$new(cp = 0.01), error = function(e) NULL)
    learners$rpart_deep <- tryCatch(Lrnr_rpart$new(cp = 0.001, maxdepth = 10), error = function(e) NULL)
  }

  # 13. HAL
  learners$hal_fast <- tryCatch(
    Lrnr_hal9001$new(max_degree = 1, num_knots = c(5), smoothness_orders = 0),
    error = function(e) NULL
  )

  # 14. GBM
  learners$gbm <- tryCatch(
    Lrnr_gbm$new(n.trees = 500, interaction.depth = 3, shrinkage = 0.05,
                  n.minobsinnode = 20),
    error = function(e) NULL
  )

  # ── Screener pipelines ─────────────────────────────────────────────────
  cat("Building screener pipelines...\n")

  safe_pipeline <- function(screener, learner, name) {
    tryCatch(make_learner(Pipeline, screener, learner),
             error = function(e) { cat(sprintf("  [skip] %s: %s\n", name, e$message)); NULL })
  }

  # Lasso-based screening
  lasso_screen <- tryCatch(
    Lrnr_screener_coefs$new(learner = Lrnr_glmnet$new(alpha = 1), threshold = 0),
    error = function(e) NULL
  )
  if (!is.null(lasso_screen)) {
    learners$lasso_then_glm    <- safe_pipeline(lasso_screen, Lrnr_glm_fast$new(), "lasso→glm")
    learners$lasso_then_ranger <- safe_pipeline(lasso_screen, Lrnr_ranger$new(num.trees = 500, min.node.size = 10, importance = "none"), "lasso→ranger")
    learners$lasso_then_xgb   <- safe_pipeline(lasso_screen, Lrnr_xgboost$new(max_depth = 3, eta = 0.05, nrounds = 200, min_child_weight = 15, eval_metric = "logloss"), "lasso→xgb")
    learners$lasso_then_earth  <- safe_pipeline(lasso_screen, Lrnr_earth$new(degree = 2), "lasso→earth")
  }

  # Correlation screening
  for (k in c(10, 30, 50)) {
    scr <- tryCatch(Lrnr_screener_correlation$new(type = "rank", num_screen = k), error = function(e) NULL)
    if (!is.null(scr)) {
      learners[[paste0("cor", k, "_ranger")]] <- safe_pipeline(scr, Lrnr_ranger$new(num.trees = 500, min.node.size = 10, importance = "none"), paste0("cor", k, "→ranger"))
      learners[[paste0("cor", k, "_xgb")]]    <- safe_pipeline(scr, Lrnr_xgboost$new(max_depth = 4, eta = 0.05, nrounds = 200, min_child_weight = 10, eval_metric = "logloss"), paste0("cor", k, "→xgb"))
      if (k <= 30) {
        learners[[paste0("cor", k, "_earth")]] <- safe_pipeline(scr, Lrnr_earth$new(degree = 2), paste0("cor", k, "→earth"))
        learners[[paste0("cor", k, "_nnet")]]  <- safe_pipeline(scr, Lrnr_nnet$new(size = 10, decay = 0.01, maxit = 200, linout = !USE_BIN, MaxNWts = nnet_maxwts), paste0("cor", k, "→nnet"))
      }
    }
  }

  # RF importance screening
  for (k in c(20, 40)) {
    scr <- tryCatch(
      Lrnr_screener_importance$new(learner = Lrnr_ranger$new(num.trees = 500, importance = "impurity"), num_screen = k),
      error = function(e) NULL
    )
    if (!is.null(scr)) {
      learners[[paste0("rfimp", k, "_ranger")]] <- safe_pipeline(scr, Lrnr_ranger$new(num.trees = 1000, min.node.size = 5, importance = "none"), paste0("rfimp", k, "→ranger"))
      learners[[paste0("rfimp", k, "_xgb")]]    <- safe_pipeline(scr, Lrnr_xgboost$new(max_depth = 4, eta = 0.05, nrounds = 300, min_child_weight = 10, eval_metric = "logloss"), paste0("rfimp", k, "→xgb"))
      learners[[paste0("rfimp", k, "_earth")]]   <- safe_pipeline(scr, Lrnr_earth$new(degree = 2), paste0("rfimp", k, "→earth"))
      learners[[paste0("rfimp", k, "_gam")]]     <- safe_pipeline(scr, Lrnr_gam$new(), paste0("rfimp", k, "→gam"))
    }
  }

  # Remove NULLs
  learners <- learners[!sapply(learners, is.null)]
  cat(sprintf("\nsl3 library: %d learners\n", length(learners)))

  # ── Build stack and SL ──────────────────────────────────────────────────
  stack <- make_learner(Stack, unname(learners))
  metalearner <- make_learner(Lrnr_nnls)
  sl <- make_learner(Lrnr_sl, learners = stack, loss_function = loss_fn, metalearner = metalearner)

  # ── Fit ─────────────────────────────────────────────────────────────────
  cat(sprintf("\nFitting sl3 SuperLearner (%d learners × %d folds)...\n",
              length(learners), K_FOLDS))
  t0_sl3 <- proc.time()
  sl_fit <- sl$train(task)
  time_sl3 <- (proc.time() - t0_sl3)["elapsed"]
  cat(sprintf("sl3 done in %.1f min\n", time_sl3 / 60))

  # ── Extract results ─────────────────────────────────────────────────────
  cv_risk <- sl_fit$cv_risk(eval_fun = loss_fn)

  # Map learner names
  learner_name_map <- setNames(names(learners), sapply(learners, function(x) x$name))
  cv_risk$learner_name <- learner_name_map[cv_risk$learner]
  cv_risk$learner_name[is.na(cv_risk$learner_name)] <- cv_risk$learner[is.na(cv_risk$learner_name)]

  # Ensemble AUC
  sl3_ensemble_auc <- NA
  if (task_type == "binomial") {
    yhat_cv <- sl_fit$predict_fold(task, "validation")
    sl3_ensemble_auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(Y, yhat_cv, quiet = TRUE))),
      error = function(e) NA
    )
    cat(sprintf("\nsl3 ensemble AUC: %.4f\n", sl3_ensemble_auc))
  }

  # Format results
  sl3_results <- cv_risk %>%
    dplyr::transmute(
      framework   = "sl3",
      learner     = learner_name,
      cv_risk     = risk,
      nnls_weight = coefficients
    ) %>%
    dplyr::arrange(cv_risk)

  sl3_results$rank <- seq_len(nrow(sl3_results))

  cat("\n=== sl3 Top 15 Learners ===\n")
  top15 <- head(sl3_results, 15)
  for (i in seq_len(nrow(top15))) {
    wt_str <- if (top15$nnls_weight[i] > 0.001) sprintf(" [wt=%.3f]", top15$nnls_weight[i]) else ""
    cat(sprintf("  %2d. risk=%.6f  %s%s\n", i, top15$cv_risk[i], top15$learner[i], wt_str))
  }

  cat(sprintf("\nNNLS selected %d / %d learners\n",
              sum(cv_risk$coefficients > 0.001), nrow(cv_risk)))

  # ── Save sl3 fit and results IMMEDIATELY ─────────────────────────────
  # (Don't wait for mlr3 part — if that errors, we lose everything)
  sl3_save_path <- file.path(model_dir, sprintf("sl3_fit_%s_%s.rds", COUNTRY, OUTCOME))
  saveRDS(sl_fit, sl3_save_path)
  cat(sprintf("sl3 fit saved to: %s\n", sl3_save_path))

  # Save sl3 results table
  sl3_csv_path <- file.path(dirname(model_dir), sprintf("sl3_results_%s_%s.csv", COUNTRY, OUTCOME))
  readr::write_csv(sl3_results, sl3_csv_path)
  cat(sprintf("sl3 results saved to: %s\n", sl3_csv_path))

  # Save ensemble metrics
  sl3_meta <- data.frame(
    framework    = "sl3",
    country      = COUNTRY,
    outcome      = OUTCOME,
    task_type    = task_type,
    n_learners   = length(learners),
    n_selected   = sum(sl3_results$nnls_weight > 0.001),
    ensemble_auc = if (exists("sl3_ensemble_auc")) sl3_ensemble_auc else NA_real_,
    best_learner = sl3_results$learner[1],
    best_risk    = sl3_results$cv_risk[1],
    time_min     = time_sl3 / 60,
    stringsAsFactors = FALSE
  )
  sl3_meta_path <- file.path(dirname(model_dir), sprintf("sl3_meta_%s_%s.csv", COUNTRY, OUTCOME))
  readr::write_csv(sl3_meta, sl3_meta_path)
  cat(sprintf("sl3 meta saved to: %s\n", sl3_meta_path))
}


# #############################################################################
#
#                     PART 2: mlr3superlearner EXPLORATION
#
# #############################################################################

mlr3_results <- NULL

if (RUN_MLR3) {

  # Raise future globals limit — sl3 fit object can be >1 GB and if it's
  # still in the environment, future will try to serialize it to workers.
  # Also explicitly remove the sl3 fit to free memory before mlr3 runs.
  options(future.globals.maxSize = 5 * 1024^3)  # 5 GB
  if (exists("sl_fit")) {
    cat("Freeing sl3 fit object from memory before mlr3...\n")
    rm(sl_fit); gc()
  }

  library(mlr3superlearner)
  library(mlr3extralearners)

  cat("\n\n")
  cat("###################################################################\n")
  cat("#                  mlr3superlearner EXPLORATION                   #\n")
  cat("###################################################################\n")

  # mlr3superlearner takes a data.frame with outcome column
  mlr3_df <- analysis_df

  # ── Define learner library ──────────────────────────────────────────────
  # mlr3superlearner uses string names + optional hyperparameters
  mlr3_library <- list(
    # Baselines
    "mean",
    "glm",

    # Regularized linear
    "glmnet",
    "cv_glmnet",

    # Tree-based
    "rpart",
    list("ranger", num.trees = 500, min.node.size = 5, id = "ranger_500_mn5"),
    list("ranger", num.trees = 500, min.node.size = 10, id = "ranger_500_mn10"),
    list("ranger", num.trees = 1000, min.node.size = 10, id = "ranger_1000_mn10"),
    list("ranger", num.trees = 1000, min.node.size = 20, id = "ranger_1000_mn20"),

    # XGBoost configs
    list("xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
         min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
         id = "xgb_conservative"),
    list("xgboost", max_depth = 5, eta = 0.05, nrounds = 300,
         min_child_weight = 10, subsample = 0.8, colsample_bytree = 0.5,
         id = "xgb_deeper"),
    list("xgboost", max_depth = 2, eta = 0.1, nrounds = 200,
         min_child_weight = 30, subsample = 0.7, colsample_bytree = 0.4,
         id = "xgb_shallow"),
    list("xgboost", max_depth = 6, eta = 0.03, nrounds = 500,
         min_child_weight = 20, subsample = 0.7, colsample_bytree = 0.4,
         id = "xgb_deep_reg"),

    # Splines
    "earth",

    # Neural net
    list("nnet", size = 5, decay = 0.1, maxit = 200, id = "nnet_5"),
    list("nnet", size = 10, decay = 0.01, maxit = 200, id = "nnet_10"),
    list("nnet", size = 20, decay = 0.001, maxit = 300, id = "nnet_20"),

    # SVM
    "svm",

    # BART
    "bart",

    # GAM
    "gam",

    # Gaussian process
    "gaussianprocess",

    # Boosting (mboost)
    "glmboost"
  )

  mlr3_type <- if (task_type == "binomial") "binomial" else "continuous"

  cat(sprintf("\nFitting mlr3superlearner (%d learners, type = %s)...\n",
              length(mlr3_library), mlr3_type))

  t0_mlr3 <- proc.time()
  mlr3_fit <- tryCatch(
    mlr3superlearner::mlr3superlearner(
      data      = mlr3_df,
      target    = "Y",
      library   = mlr3_library,
      outcome_type = mlr3_type,
      folds     = K_FOLDS
    ),
    error = function(e) {
      cat(sprintf("mlr3superlearner FAILED: %s\n", e$message))
      NULL
    }
  )
  time_mlr3 <- (proc.time() - t0_mlr3)["elapsed"]
  cat(sprintf("mlr3superlearner done in %.1f min\n", time_mlr3 / 60))

  if (!is.null(mlr3_fit)) {
    # Print the fit summary
    cat("\n=== mlr3superlearner Results ===\n")
    print(mlr3_fit)

    # Extract risk and weights from the fit object
    mlr3_summary <- tryCatch({
      # The print output has Risk and Coefficients
      fit_data <- mlr3_fit$fits
      risks <- mlr3_fit$risk
      coefs <- mlr3_fit$coef

      data.frame(
        framework   = "mlr3",
        learner     = names(risks),
        cv_risk     = as.numeric(risks),
        nnls_weight = as.numeric(coefs),
        stringsAsFactors = FALSE
      ) %>% dplyr::arrange(cv_risk)
    }, error = function(e) {
      cat(sprintf("  Could not extract mlr3 summary: %s\n", e$message))
      # Try alternative extraction
      tryCatch({
        # Capture print output
        out <- capture.output(print(mlr3_fit))
        cat("  Raw output:\n")
        cat(paste(out, collapse = "\n"), "\n")
        NULL
      }, error = function(e2) NULL)
    })

    if (!is.null(mlr3_summary)) {
      mlr3_summary$rank <- seq_len(nrow(mlr3_summary))
      mlr3_results <- mlr3_summary

      cat("\n=== mlr3 Top 15 Learners ===\n")
      top15m <- head(mlr3_results, 15)
      for (i in seq_len(nrow(top15m))) {
        wt_str <- if (top15m$nnls_weight[i] > 0.001) sprintf(" [wt=%.3f]", top15m$nnls_weight[i]) else ""
        cat(sprintf("  %2d. risk=%.6f  %s%s\n", i, top15m$cv_risk[i], top15m$learner[i], wt_str))
      }
    }

    # Ensemble AUC
    mlr3_ensemble_auc <- NA
    if (task_type == "binomial") {
      preds_mlr3 <- tryCatch(predict(mlr3_fit, mlr3_df), error = function(e) NULL)
      if (!is.null(preds_mlr3)) {
        mlr3_ensemble_auc <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(Y, preds_mlr3, quiet = TRUE))),
          error = function(e) NA
        )
        cat(sprintf("\nmlr3 ensemble AUC (resubstitution): %.4f\n", mlr3_ensemble_auc))
      }
    }

    # Save fit
    mlr3_save_path <- file.path(model_dir, sprintf("mlr3_fit_%s_%s.rds", COUNTRY, OUTCOME))
    saveRDS(mlr3_fit, mlr3_save_path)
    cat(sprintf("mlr3 fit saved to: %s\n", mlr3_save_path))
  }
}


# #############################################################################
#
#                       PART 3: COMPARISON
#
# #############################################################################

cat("\n\n")
cat("###################################################################\n")
cat("#                       COMPARISON                               #\n")
cat("###################################################################\n")

# ── Combine results ───────────────────────────────────────────────────────────
all_results <- dplyr::bind_rows(sl3_results, mlr3_results)

if (nrow(all_results) > 0) {
  all_results$country <- COUNTRY
  all_results$outcome <- OUTCOME

  # Save combined table
  comparison_file <- file.path(out_dir, sprintf("sl_comparison_%s_%s.csv", COUNTRY, OUTCOME))
  readr::write_csv(all_results, comparison_file)
  cat(sprintf("\nComparison table saved to: %s\n", comparison_file))

  # Print side-by-side summary
  cat("\n=== Framework Comparison ===\n")
  cat(sprintf("%-20s  %10s  %10s\n", "Metric", "sl3", "mlr3"))
  cat(sprintf("%-20s  %10s  %10s\n", "----", "----", "----"))

  if (!is.null(sl3_results))
    cat(sprintf("%-20s  %10d  %10s\n", "Learners tested",
                nrow(sl3_results),
                if (!is.null(mlr3_results)) nrow(mlr3_results) else "N/A"))

  if (!is.null(sl3_results))
    cat(sprintf("%-20s  %10d  %10s\n", "NNLS selected",
                sum(sl3_results$nnls_weight > 0.001),
                if (!is.null(mlr3_results)) sum(mlr3_results$nnls_weight > 0.001) else "N/A"))

  if (RUN_SL3)
    cat(sprintf("%-20s  %10.1f  %10s\n", "Time (min)",
                time_sl3 / 60,
                if (RUN_MLR3) sprintf("%.1f", time_mlr3 / 60) else "N/A"))

  if (exists("sl3_ensemble_auc") && exists("mlr3_ensemble_auc")) {
    cat(sprintf("%-20s  %10.4f  %10s\n", "Ensemble AUC",
                sl3_ensemble_auc,
                if (!is.na(mlr3_ensemble_auc)) sprintf("%.4f*", mlr3_ensemble_auc) else "N/A"))
    cat("* mlr3 AUC is resubstitution (not CV) — expect it to be higher\n")
  }

  # Best learner per framework
  if (!is.null(sl3_results))
    cat(sprintf("\nsl3 best:  %s (risk = %.6f)\n", sl3_results$learner[1], sl3_results$cv_risk[1]))
  if (!is.null(mlr3_results))
    cat(sprintf("mlr3 best: %s (risk = %.6f)\n", mlr3_results$learner[1], mlr3_results$cv_risk[1]))

  # Shared learner comparison
  if (!is.null(sl3_results) && !is.null(mlr3_results)) {
    cat("\n=== Shared Learner Rankings ===\n")
    # Map common names (approximate)
    common_map <- data.frame(
      sl3_pattern  = c("^lasso$", "^ridge$", "^ranger_500_mn10$", "^xgb_conservative$",
                        "^earth_d1$", "^nnet_10$", "^mean$", "^dbarts$", "^gam$", "^rpart$"),
      mlr3_pattern = c("glmnet", "glmnet", "ranger_500_mn10", "xgb_conservative",
                        "earth", "nnet_10", "featureless|mean", "bart", "gam", "rpart"),
      label        = c("Lasso/glmnet", "Ridge", "Ranger(500,mn10)", "XGB(conservative)",
                        "Earth", "NNet(10)", "Mean", "BART", "GAM", "Rpart"),
      stringsAsFactors = FALSE
    )
    cat(sprintf("%-20s  %8s %8s  %8s %8s\n", "Learner", "sl3_rank", "sl3_risk", "mlr3_rank", "mlr3_risk"))
    for (i in seq_len(nrow(common_map))) {
      sl3_row  <- sl3_results[grepl(common_map$sl3_pattern[i], sl3_results$learner, ignore.case = TRUE), ]
      mlr3_row <- mlr3_results[grepl(common_map$mlr3_pattern[i], mlr3_results$learner, ignore.case = TRUE), ]
      if (nrow(sl3_row) > 0 || nrow(mlr3_row) > 0) {
        cat(sprintf("%-20s  %8s %8s  %8s %8s\n",
                    common_map$label[i],
                    if (nrow(sl3_row) > 0) sl3_row$rank[1] else "—",
                    if (nrow(sl3_row) > 0) sprintf("%.5f", sl3_row$cv_risk[1]) else "—",
                    if (nrow(mlr3_row) > 0) mlr3_row$rank[1] else "—",
                    if (nrow(mlr3_row) > 0) sprintf("%.5f", mlr3_row$cv_risk[1]) else "—"))
      }
    }
  }
}


# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n\n=== FILES SAVED ===\n")
cat(sprintf("  Comparison table: %s\n", file.path(out_dir, sprintf("sl_comparison_%s_%s.csv", COUNTRY, OUTCOME))))
if (RUN_SL3)  cat(sprintf("  sl3 model:  %s\n", file.path(model_dir, sprintf("sl3_fit_%s_%s.rds", COUNTRY, OUTCOME))))
if (RUN_MLR3) cat(sprintf("  mlr3 model: %s\n", file.path(model_dir, sprintf("mlr3_fit_%s_%s.rds", COUNTRY, OUTCOME))))

cat("\nDone!\n")
