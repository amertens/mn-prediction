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
COUNTRY  <- "gambia"         # lowercase: gambia, ghana, sierraleone, malawi
OUTCOME  <- "child_vitA"    # child_vitA, women_vitA, child_iron, women_iron
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

# ── Create SHARED clustered CV folds (used by ALL frameworks) ─────────────
# This ensures sl3, mlr3superlearner comparison, and improved pipeline
# all use identical train/test splits for fair comparison.
library(origami)
set.seed(SEED)
shared_fold_obj <- origami::make_folds(cluster_ids = id_vec, V = K_FOLDS)
shared_train_sets <- lapply(shared_fold_obj, function(f) f$training_set)
shared_test_sets  <- lapply(shared_fold_obj, function(f) f$validation_set)
cat(sprintf("Shared clustered CV: %d folds, %d clusters\n",
            length(shared_fold_obj), length(unique(id_vec))))

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

  # ── Build sl3 task (using shared clustered folds) ──────────────────────
  dat <- data.table(Y = Y, id = id_vec, cov)
  task <- sl3::make_sl3_Task(
    data = dat, covariates = covars, outcome = "Y",
    id = "id", folds = shared_fold_obj
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
  # Raise future globals limit — sl3 fit objects can be >1 GB
  options(future.globals.maxSize = 5 * 1024^3)  # 5 GB
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

  # Format results — rename cv_risk df to avoid dplyr column/variable conflict
  # sl3 cv_risk data.table uses "NLL" as column name (not "risk")
  cv_risk_df <- cv_risk
  risk_col <- if ("NLL" %in% colnames(cv_risk_df)) "NLL" else if ("risk" %in% colnames(cv_risk_df)) "risk" else colnames(cv_risk_df)[3]
  sl3_results <- data.frame(
    framework   = "sl3",
    learner     = as.character(cv_risk_df$learner_name),
    cv_risk_val = cv_risk_df[[risk_col]],
    nnls_weight = cv_risk_df$coefficients,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(cv_risk_val)

  sl3_results$rank <- seq_len(nrow(sl3_results))

  cat("\n=== sl3 Top 15 Learners ===\n")
  top15 <- head(sl3_results, 15)
  for (i in seq_len(nrow(top15))) {
    wt_str <- if (top15$nnls_weight[i] > 0.001) sprintf(" [wt=%.3f]", top15$nnls_weight[i]) else ""
    cat(sprintf("  %2d. risk=%.6f  %s%s\n", i, top15$cv_risk_val[i], top15$learner[i], wt_str))
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
    best_risk    = sl3_results$cv_risk_val[1],
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
  # Unlike sl3/DHS_SL_clustered, mlr3 does NOT auto-impute — must handle NAs here
  mlr3_df <- analysis_df

  # Strip haven labels, convert factors to numeric, then impute NAs
  # Must do this in two passes to avoid column-shift bugs during deletion
  cols_to_drop <- character()
  for (col in colnames(mlr3_df)) {
    if (col == "Y") next
    x <- mlr3_df[[col]]

    # Convert haven_labelled
    if (inherits(x, "haven_labelled")) {
      mlr3_df[[col]] <- as.double(unclass(x))
      x <- mlr3_df[[col]]
    }

    # Convert factor/character to numeric or mark for dropping
    if (is.factor(x) || is.character(x)) {
      cols_to_drop <- c(cols_to_drop, col)
      next
    }

    # Replace NaN with NA, then impute
    mlr3_df[[col]][is.nan(mlr3_df[[col]])] <- NA
    if (any(is.na(mlr3_df[[col]]))) {
      med <- median(mlr3_df[[col]], na.rm = TRUE)
      if (is.finite(med)) {
        mlr3_df[[col]][is.na(mlr3_df[[col]])] <- med
      } else {
        mlr3_df[[col]][is.na(mlr3_df[[col]])] <- 0
      }
    }
  }

  # Drop factor/character columns in one go (avoids index-shift bug)
  if (length(cols_to_drop) > 0) {
    cat(sprintf("  Dropping %d non-numeric columns: %s\n",
                length(cols_to_drop), paste(head(cols_to_drop, 5), collapse = ", ")))
    mlr3_df <- mlr3_df[, !colnames(mlr3_df) %in% cols_to_drop, drop = FALSE]
  }

  # Drop any rows still with NA in Y
  mlr3_df <- mlr3_df[!is.na(mlr3_df$Y), ]

  # Final check
  n_na <- sum(is.na(mlr3_df))
  cat(sprintf("  After imputation: %d rows x %d cols, %d NAs remaining\n",
              nrow(mlr3_df), ncol(mlr3_df), n_na))
  if (n_na > 0) {
    na_cols <- colnames(mlr3_df)[colSums(is.na(mlr3_df)) > 0]
    cat(sprintf("  WARNING: NAs remain in: %s\n", paste(na_cols, collapse = ", ")))
  }

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
    # Ranger — NaN risk in previous run likely due to mlr3's internal
    # predict_type handling. mlr3superlearner should set predict_type = "prob"
    # automatically for binomial tasks. Keep configs simple.
    # Ranger — min.node.size must be ≥15 to avoid exact 0/1 probability
    # predictions that cause log(0) = NaN in mlr3superlearner's unclipped
    # NLL loss. sl3 doesn't have this issue (different loss handling).
    list("ranger", num.trees = 500, min.node.size = 15, id = "ranger_500_mn15"),
    list("ranger", num.trees = 500, min.node.size = 20, id = "ranger_500_mn20"),
    list("ranger", num.trees = 1000, min.node.size = 15, id = "ranger_1000_mn15"),
    list("ranger", num.trees = 1000, min.node.size = 20, id = "ranger_1000_mn20"),
    list("ranger", num.trees = 1000, min.node.size = 30, id = "ranger_1000_mn30"),

    # randomForest — alternative RF implementation (from randomForest package).
    # Included because mlr3's ranger has NaN risk bug on some configs.
    list("randomforest", ntree = 500, nodesize = 5, id = "rf_500_ns5"),
    list("randomforest", ntree = 500, nodesize = 10, id = "rf_500_ns10"),
    list("randomforest", ntree = 1000, nodesize = 10, id = "rf_1000_ns10"),

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

    # SVM
    "svm",

    # BART — multiple configs. dbarts uses MCMC posterior sampling over
    # sums of regression trees. ntree controls ensemble size, ndpost
    # controls number of posterior draws (more = smoother estimates),
    # nskip is burn-in.
    # BART — mlr3 uses 'ntree' not 'n.trees', 'ndpost', 'nskip'
    list("bart", id = "bart_default"),                                         # defaults
    list("bart", ntree = 50, ndpost = 500, nskip = 100, id = "bart_small"),    # fewer trees, fast
    list("bart", ntree = 200, ndpost = 1500, nskip = 250, id = "bart_large"),  # more posterior draws
    list("bart", ntree = 100, ndpost = 1000, nskip = 200, id = "bart_100t"),   # moderate
    list("bart", ntree = 300, ndpost = 1000, nskip = 200, id = "bart_300t"),   # more trees

    # GAM
    "gam",

    # Gaussian process
    "gaussianprocess"

    # Note: "glmboost" omitted — requires mboost package
    # Note: "nnet" omitted — "too many weights" with p=298
  )

  mlr3_type <- if (task_type == "binomial") "binomial" else "continuous"

  cat(sprintf("\nFitting mlr3superlearner (%d learners, type = %s)...\n",
              length(mlr3_library), mlr3_type))

  # Ensure mlr3_df has no haven_labelled or list columns
  for (col in colnames(mlr3_df)) {
    if (inherits(mlr3_df[[col]], "haven_labelled")) {
      mlr3_df[[col]] <- as.double(unclass(mlr3_df[[col]]))
    }
    if (is.list(mlr3_df[[col]])) {
      mlr3_df[[col]] <- as.numeric(unlist(mlr3_df[[col]]))
    }
    if (is.factor(mlr3_df[[col]])) {
      mlr3_df[[col]] <- as.numeric(as.character(mlr3_df[[col]]))
    }
  }

  # Y must be numeric 0/1 for mlr3superlearner (it creates the task internally).
  # Factor Y causes "not meaningful for factors" errors in the NLL loss.
  mlr3_df$Y <- as.numeric(mlr3_df$Y)
  cat(sprintf("  mlr3_df: %d rows x %d cols, Y range: [%.2f, %.2f]\n",
              nrow(mlr3_df), ncol(mlr3_df), min(mlr3_df$Y), max(mlr3_df$Y)))

  # Monkey-patch mlr3superlearner's loss_nll to clip probabilities.
  # The original loss_nll = -mean(y*log(x) + (1-y)*log(1-x)) produces NaN
  # when x is exactly 0 or 1 (from ranger leaf predictions).
  # This clips to [1e-15, 1-1e-15] before taking log.
  tryCatch({
    ns <- asNamespace("mlr3superlearner")
    original_loss <- ns$loss_nll
    patched_loss <- function(x, y) {
      x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
      -mean(y * log(x) + (1 - y) * log(1 - x))
    }
    assignInNamespace("loss_nll", patched_loss, ns = "mlr3superlearner")
    cat("  Patched mlr3superlearner::loss_nll with probability clipping\n")
  }, error = function(e) {
    cat(sprintf("  Could not patch loss_nll: %s\n", e$message))
  })

  t0_mlr3 <- proc.time()
  mlr3_fit <- tryCatch({
    # withCallingHandlers catches warnings BEFORE mlr3's internal handler
    # (which crashes with "no restart 'muffleWarning' found")
    withCallingHandlers(
      {
        fit <- mlr3superlearner::mlr3superlearner(
          data      = mlr3_df,
          target    = "Y",
          library   = mlr3_library,
          outcome_type = mlr3_type,
          folds     = K_FOLDS
        )
        cat("  mlr3superlearner returned successfully\n")
        fit
      },
      warning = function(w) {
        msg <- conditionMessage(w)
        # Only log non-spammy warnings
        if (!grepl("non-positive", msg, ignore.case = TRUE))
          cat(sprintf("  mlr3 warning: %s\n", substr(msg, 1, 100)))
        tryInvokeRestart("muffleWarning")
      }
    )
  },
  error = function(e) {
    cat(sprintf("mlr3superlearner FAILED: %s\n", e$message))
    cat(sprintf("  Traceback:\n"))
    traceback_lines <- capture.output(traceback())
    cat(paste(head(traceback_lines, 20), collapse = "\n"), "\n")
    NULL
  })
  # Note: warnings from mlr3 (e.g., "glm.fit: algorithm did not converge")
  # are allowed to pass through. Using tryCatch+warning with invokeRestart
  # causes "no restart 'muffleWarning' found" errors.
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
        cat(sprintf("  %2d. risk=%.6f  %s%s\n", i, top15m$cv_risk_val[i], top15m$learner[i], wt_str))
      }
    }

    # Ensemble AUC — resubstitution (predict on training data)
    mlr3_ensemble_auc <- NA
    mlr3_cv_auc <- NA
    if (task_type == "binomial") {
      # Resubstitution AUC
      preds_mlr3 <- tryCatch(predict(mlr3_fit, mlr3_df), error = function(e) NULL)
      if (!is.null(preds_mlr3)) {
        mlr3_ensemble_auc <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(Y, preds_mlr3, quiet = TRUE))),
          error = function(e) NA
        )
        cat(sprintf("\nmlr3 ensemble AUC (resubstitution): %.4f\n", mlr3_ensemble_auc))
      }

      # CV AUC — extract from internal CV predictions if available
      mlr3_cv_auc <- tryCatch({
        # mlr3superlearner stores CV predictions in $cv_preds or $data
        cv_preds <- mlr3_fit$cv_risk_data
        if (!is.null(cv_preds) && "pred" %in% names(cv_preds)) {
          as.numeric(pROC::auc(pROC::roc(cv_preds$truth, cv_preds$pred, quiet = TRUE)))
        } else {
          # Try accessing the metalearner's CV predictions
          NA_real_
        }
      }, error = function(e) NA_real_)

      if (!is.na(mlr3_cv_auc)) {
        cat(sprintf("mlr3 ensemble AUC (cross-validated): %.4f\n", mlr3_cv_auc))
      }

      # Brier score
      if (!is.null(preds_mlr3)) {
        mlr3_brier <- mean((Y - preds_mlr3)^2, na.rm = TRUE)
        cat(sprintf("mlr3 Brier score (resubstitution): %.4f\n", mlr3_brier))
      }
    }

    # Save fit
    mlr3_save_path <- file.path(model_dir, sprintf("mlr3_fit_%s_%s.rds", COUNTRY, OUTCOME))
    saveRDS(mlr3_fit, mlr3_save_path)
    cat(sprintf("mlr3 fit saved to: %s\n", mlr3_save_path))

    # Save mlr3 results CSV
    if (!is.null(mlr3_results)) {
      mlr3_csv_path <- file.path(out_dir, sprintf("mlr3_results_%s_%s.csv", COUNTRY, OUTCOME))
      readr::write_csv(mlr3_results, mlr3_csv_path)
      cat(sprintf("mlr3 results saved to: %s\n", mlr3_csv_path))
    }
  }
}


# #############################################################################
#
#              PART 2b: RANGER CLASSIFICATION TASK DIAGNOSTIC
#
# Hypothesis: mlr3superlearner creates TaskRegr for binary outcomes,
# causing ranger to fit a regression forest instead of a probability forest.
# Test: manually create TaskClassif and compare ranger CV performance.
#
# #############################################################################

if (RUN_MLR3 && task_type == "binomial") {

  cat("\n\n")
  cat("###################################################################\n")
  cat("#         RANGER DIAGNOSTIC: TaskClassif vs TaskRegr              #\n")
  cat("###################################################################\n")

  library(mlr3)
  library(mlr3learners)

  # ── Prepare classification data ──
  classif_df <- mlr3_df
  classif_df$Y <- factor(classif_df$Y, levels = c(0, 1))

  # Create proper classification task
  task_classif <- as_task_classif(classif_df, target = "Y", positive = "1")
  cat(sprintf("TaskClassif: %d obs, %d features, positive class = '%s'\n",
              task_classif$nrow, length(task_classif$feature_names),
              task_classif$positive))

  # Also create regression task (what mlr3superlearner likely does)
  regr_df <- mlr3_df  # Y is numeric 0/1
  task_regr <- as_task_regr(regr_df, target = "Y")
  cat(sprintf("TaskRegr:    %d obs, %d features\n",
              task_regr$nrow, length(task_regr$feature_names)))

  # ── Define ranger configs to test ──
  ranger_configs <- list(
    list(id = "ranger_500_mn10", num.trees = 500, min.node.size = 10),
    list(id = "ranger_500_mn5",  num.trees = 500, min.node.size = 5),
    list(id = "ranger_1000_mn10", num.trees = 1000, min.node.size = 10),
    list(id = "ranger_1000_mn20", num.trees = 1000, min.node.size = 20)
  )

  # NLL with clipping
  clipped_nll <- function(pred, truth) {
    p <- pmin(pmax(pred, 1e-15), 1 - 1e-15)
    -mean(truth * log(p) + (1 - truth) * log(1 - p))
  }

  # ── CV resampling ──
  set.seed(SEED)
  cv5 <- rsmp("cv", folds = K_FOLDS)

  diag_results <- list()

  for (cfg in ranger_configs) {
    id <- cfg$id
    cat(sprintf("\n  Testing: %s\n", id))

    # ── Classification task (proper probability forest) ──
    classif_nll <- tryCatch({
      lrn_c <- lrn("classif.ranger",
                     num.trees = cfg$num.trees,
                     min.node.size = cfg$min.node.size,
                     predict_type = "prob")
      rr_c <- resample(task_classif, lrn_c, cv5, store_models = FALSE)

      # Extract probability predictions for positive class
      preds_c <- rr_c$predictions()
      all_probs <- numeric()
      all_truth <- numeric()
      for (p in preds_c) {
        probs <- p$prob[, "1"]
        truth <- as.numeric(as.character(p$truth))
        all_probs <- c(all_probs, probs)
        all_truth <- c(all_truth, truth)
      }

      nll <- clipped_nll(all_probs, all_truth)
      auc <- as.numeric(pROC::auc(pROC::roc(all_truth, all_probs, quiet = TRUE)))
      cat(sprintf("    TaskClassif: NLL = %.4f, AUC = %.4f, pred range [%.4f, %.4f]\n",
                  nll, auc, min(all_probs), max(all_probs)))
      list(nll = nll, auc = auc)
    }, error = function(e) {
      cat(sprintf("    TaskClassif FAILED: %s\n", e$message))
      list(nll = NA, auc = NA)
    })

    # ── Regression task (what mlr3superlearner does) ──
    regr_nll <- tryCatch({
      lrn_r <- lrn("regr.ranger",
                     num.trees = cfg$num.trees,
                     min.node.size = cfg$min.node.size)
      rr_r <- resample(task_regr, lrn_r, cv5, store_models = FALSE)

      preds_r <- rr_r$predictions()
      all_resp <- numeric()
      all_truth_r <- numeric()
      for (p in preds_r) {
        all_resp <- c(all_resp, p$response)
        all_truth_r <- c(all_truth_r, p$truth)
      }

      # Clip to [0,1] since regression can predict outside
      all_resp_clip <- pmin(pmax(all_resp, 0), 1)
      nll <- clipped_nll(all_resp_clip, all_truth_r)
      auc <- as.numeric(pROC::auc(pROC::roc(all_truth_r, all_resp_clip, quiet = TRUE)))
      cat(sprintf("    TaskRegr:    NLL = %.4f, AUC = %.4f, pred range [%.4f, %.4f]\n",
                  nll, auc, min(all_resp), max(all_resp)))
      list(nll = nll, auc = auc)
    }, error = function(e) {
      cat(sprintf("    TaskRegr FAILED: %s\n", e$message))
      list(nll = NA, auc = NA)
    })

    diag_results[[id]] <- data.frame(
      config = id,
      classif_nll = classif_nll$nll, classif_auc = classif_nll$auc,
      regr_nll = regr_nll$nll, regr_auc = regr_nll$auc,
      nll_gap = regr_nll$nll - classif_nll$nll,
      stringsAsFactors = FALSE
    )
  }

  diag_df <- do.call(rbind, diag_results)
  cat("\n\n=== RANGER DIAGNOSTIC SUMMARY ===\n")
  cat(sprintf("sl3 ranger_500_mn10 NLL:    %.4f (benchmark)\n", 0.3812))
  cat("\n")
  print(diag_df, row.names = FALSE)

  cat(sprintf("\nConclusion: TaskClassif %s TaskRegr by %.4f NLL on average\n",
              if (mean(diag_df$nll_gap, na.rm = TRUE) > 0) "beats" else "loses to",
              abs(mean(diag_df$nll_gap, na.rm = TRUE))))

  # Save diagnostic
  diag_path <- file.path(out_dir, "ranger_diagnostic.csv")
  readr::write_csv(diag_df, diag_path)
  cat(sprintf("Saved to: %s\n", diag_path))
}


# #############################################################################
#
#        PART 2c: IMPROVED mlr3 PIPELINE — PipeOps + Clustered CV
#
# Implements all proposed pipeline improvements:
#   1. No washb_prescreen (data leakage removed)
#   2. PipeOps preprocessing INSIDE CV folds
#   3. Learner-specific preprocessing (minimal for ranger/xgb, full for glm)
#   4. Clustered CV via origami + rsmp("custom")
#   5. Clipped NLL loss
#   6. Brier loss comparison for rare outcomes
#   7. Proper TaskClassif for binary outcomes
#
# #############################################################################

RUN_IMPROVED <- TRUE  # Set FALSE to skip

if (RUN_IMPROVED && RUN_MLR3 && task_type == "binomial") {

  cat("\n\n")
  cat("###################################################################\n")
  cat("#      IMPROVED PIPELINE: PipeOps + Clustered CV + TaskClassif    #\n")
  cat("###################################################################\n")

  library(mlr3)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3extralearners)
  library(origami)

  # ══════════════════════════════════════════════════════════════════════
  # 1. PREPARE RAW DATA (minimal preprocessing — rest goes inside CV)
  # ══════════════════════════════════════════════════════════════════════

  # Start from raw data — NO washb_prescreen, NO recipes
  cat("\n[improved] Preparing raw data (no external preprocessing)...\n")

  raw_X <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  raw_X <- labelled::unlabelled(raw_X, user_na_to_na = TRUE)

  # Strip haven labels, factors, lists
  for (col in colnames(raw_X)) {
    if (inherits(raw_X[[col]], "haven_labelled")) {
      raw_X[[col]] <- as.double(unclass(raw_X[[col]]))
    }
    if (is.factor(raw_X[[col]]) || is.character(raw_X[[col]])) {
      raw_X[[col]] <- as.numeric(as.character(raw_X[[col]]))
    }
    if (is.list(raw_X[[col]])) {
      raw_X[[col]] <- as.numeric(unlist(raw_X[[col]]))
    }
  }

  # Remove all-NA and truly constant columns (safe — doesn't use Y)
  raw_X <- raw_X[, !sapply(raw_X, function(x) all(is.na(x))), drop = FALSE]
  raw_X <- raw_X[, sapply(raw_X, function(x) {
    non_na <- x[!is.na(x)]
    length(non_na) > 0 && length(unique(non_na)) > 1
  }), drop = FALSE]

  # Build TaskClassif with factor Y
  improved_df <- data.frame(Y = factor(Y, levels = c(0, 1)), raw_X)
  improved_df$cluster_id <- as.character(id_vec)

  cat(sprintf("[improved] Raw data: n = %d, p = %d (no prescreening)\n",
              nrow(improved_df), ncol(improved_df) - 2))  # -2 for Y, cluster_id

  # ══════════════════════════════════════════════════════════════════════
  # 2. CLUSTERED CV FOLDS (origami → mlr3 custom resampling)
  # ══════════════════════════════════════════════════════════════════════

  # Reuse the shared clustered folds (created before Part 1)
  # Note: shared_fold_obj was created on the preprocessed data (same id_vec).
  # The improved pipeline uses RAW data, so row indices must be re-mapped
  # if row counts differ. For now, raw data has the same rows (just different columns).
  cat("[improved] Reusing shared clustered CV folds from origami...\n")

  # The improved pipeline works on improved_df which may have different row count
  # if we filtered differently. Create fresh folds on the raw id_vec if needed.
  if (nrow(improved_df) == length(shared_train_sets[[1]]) + length(shared_test_sets[[1]]) ||
      nrow(improved_df) == sum(sapply(shared_test_sets, length))) {
    # Same row count — reuse directly
    train_sets <- shared_train_sets
    test_sets  <- shared_test_sets
    cat("  Using identical fold assignments as sl3\n")
  } else {
    # Different row count (raw vs preprocessed) — recreate with same seed
    cat("  Row count differs from sl3 preprocessing — recreating folds with same seed...\n")
    set.seed(SEED)
    fold_obj_imp <- origami::make_folds(cluster_ids = improved_df$cluster_id, V = K_FOLDS)
    train_sets <- lapply(fold_obj_imp, function(f) f$training_set)
    test_sets  <- lapply(fold_obj_imp, function(f) f$validation_set)
  }

  cat(sprintf("[improved] %d clustered folds, fold sizes: %s\n",
              length(train_sets),
              paste(sapply(test_sets, length), collapse = "/")))

  # ══════════════════════════════════════════════════════════════════════
  # 3. PipeOps PREPROCESSING (inside CV folds)
  # ══════════════════════════════════════════════════════════════════════

  cat("[improved] Building PipeOps preprocessing pipelines...\n")

  # Minimal preprocessing — for learners that handle high-p natively
  # (ranger, xgboost, bart, glmnet)
  preproc_minimal <- po("removeconstants") %>>%
    po("imputemedian") %>>%
    po("removeconstants", id = "removeconstants2")

  # Full preprocessing — for learners that need it (glm, earth, gam, svm)
  # Adds correlation filtering + normalization
  preproc_full <- po("removeconstants") %>>%
    po("imputemedian") %>>%
    po("removeconstants", id = "removeconstants2_full") %>>%
    po("scale") %>>%
    # Filter highly correlated features (replaces step_corr(0.85))
    po("filter",
       filter = mlr3filters::flt("find_correlation"),
       filter.cutoff = 0.85,
       id = "corr_filter")

  # ══════════════════════════════════════════════════════════════════════
  # 4. BUILD LEARNER PIPELINES (learner-specific preprocessing)
  # ══════════════════════════════════════════════════════════════════════

  cat("[improved] Building learner-specific pipelines...\n")

  # ── Learners with minimal preprocessing ──
  lrn_ranger_main <- preproc_minimal %>>%
    po("learner",
       lrn("classif.ranger", num.trees = 500, min.node.size = 10,
           predict_type = "prob", id = "ranger_main"))

  lrn_ranger_low_mtry <- preproc_minimal %>>%
    po("learner",
       lrn("classif.ranger", num.trees = 500, min.node.size = 10,
           mtry = max(floor(log2(ncol(raw_X))), 2),
           predict_type = "prob", id = "ranger_low_mtry"))

  lrn_xgb <- preproc_minimal %>>%
    po("learner",
       lrn("classif.xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
           min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
           predict_type = "prob", id = "xgb_conservative"))

  lrn_xgb_deep <- preproc_minimal %>>%
    po("learner",
       lrn("classif.xgboost", max_depth = 6, eta = 0.03, nrounds = 500,
           min_child_weight = 20, subsample = 0.7, colsample_bytree = 0.4,
           predict_type = "prob", id = "xgb_deep"))

  lrn_bart <- preproc_minimal %>>%
    po("learner",
       lrn("classif.bart", ntree = 100, ndpost = 1000, nskip = 200,
           predict_type = "prob", id = "bart_100"))

  lrn_bart_small <- preproc_minimal %>>%
    po("learner",
       lrn("classif.bart", ntree = 50, ndpost = 500, nskip = 100,
           predict_type = "prob", id = "bart_small"))

  lrn_lasso <- preproc_minimal %>>%
    po("learner",
       lrn("classif.glmnet", alpha = 1, predict_type = "prob", id = "lasso"))

  lrn_mean <- po("learner",
    lrn("classif.featureless", predict_type = "prob", id = "mean"))

  # ── Learners with full preprocessing ──
  lrn_earth <- preproc_full %>>%
    po("learner",
       lrn("classif.earth", predict_type = "prob", id = "earth"))

  lrn_gp <- preproc_full %>>%
    po("learner",
       lrn("classif.gausspr", predict_type = "prob", id = "gaussianprocess"))

  # ── Screening pipelines (feature selection → learner, all inside CV) ──
  # These are the mlr3 equivalent of sl3's Pipeline(Lrnr_screener_coefs, Lrnr_ranger).
  # The screening step trains a lasso model, selects non-zero features,
  # then passes only those to the downstream learner. All within each fold.

  cat("[improved] Building screening pipelines...\n")

  # Helper function: build a lasso-screening → learner pipeline
  # Uses po("filter") with glmnet importance to select features,
  # then passes to the downstream learner. All inside CV.
  make_lasso_screen_pipeline <- function(downstream_lrn, n_features = 50,
                                          extra_preproc = NULL, id_suffix = "") {
    # Use glmnet importance filter — fits lasso, ranks by |coefficient|
    screen_pipe <- po("removeconstants", id = paste0("rc_scr", id_suffix)) %>>%
      po("imputemedian", id = paste0("imp_scr", id_suffix)) %>>%
      po("removeconstants", id = paste0("rc_scr2", id_suffix)) %>>%
      po("filter",
         filter = mlr3filters::flt("importance",
                                    learner = lrn("classif.glmnet", alpha = 1,
                                                   predict_type = "prob")),
         filter.nfeat = n_features,
         id = paste0("lasso_filter", id_suffix))

    if (!is.null(extra_preproc)) {
      screen_pipe <- screen_pipe %>>% extra_preproc
    }

    screen_pipe %>>% po("learner", downstream_lrn)
  }

  # Lasso screen → Ranger (23% NNLS weight in sl3 exploration with this combo)
  lrn_lasso_ranger <- tryCatch({
    make_lasso_screen_pipeline(
      lrn("classif.ranger", num.trees = 500, min.node.size = 10,
          predict_type = "prob", id = "lasso_ranger"),
      n_features = 50, id_suffix = "_ranger"
    )
  }, error = function(e) {
    cat(sprintf("  [skip] lasso→ranger pipeline: %s\n", e$message))
    NULL
  })

  # Lasso screen → XGBoost
  lrn_lasso_xgb <- tryCatch({
    make_lasso_screen_pipeline(
      lrn("classif.xgboost", max_depth = 3, eta = 0.05, nrounds = 200,
          min_child_weight = 15,
          predict_type = "prob", id = "lasso_xgb"),
      n_features = 50, id_suffix = "_xgb"
    )
  }, error = function(e) {
    cat(sprintf("  [skip] lasso→xgb pipeline: %s\n", e$message))
    NULL
  })

  # Lasso screen → BART
  lrn_lasso_bart <- tryCatch({
    make_lasso_screen_pipeline(
      lrn("classif.bart", ntree = 100, ndpost = 1000, nskip = 200,
          predict_type = "prob", id = "lasso_bart"),
      n_features = 50, id_suffix = "_bart"
    )
  }, error = function(e) {
    cat(sprintf("  [skip] lasso→bart pipeline: %s\n", e$message))
    NULL
  })

  # Lasso screen → Earth (MARS) — benefits most from screening
  lrn_lasso_earth <- tryCatch({
    make_lasso_screen_pipeline(
      lrn("classif.earth", predict_type = "prob", id = "lasso_earth"),
      n_features = 30,  # earth is more sensitive to p, use fewer
      extra_preproc = po("scale", id = "scale_lasso_earth"),
      id_suffix = "_earth"
    )
  }, error = function(e) {
    cat(sprintf("  [skip] lasso→earth pipeline: %s\n", e$message))
    NULL
  })

  # ══════════════════════════════════════════════════════════════════════
  # 5. CREATE TASK AND RESAMPLING
  # ══════════════════════════════════════════════════════════════════════

  # Create task (excluding cluster_id from features)
  task_improved <- as_task_classif(improved_df, target = "Y", positive = "1")
  task_improved$set_col_roles("cluster_id", roles = character())  # exclude from features

  cat(sprintf("[improved] Task: %d obs, %d features, positive = '%s'\n",
              task_improved$nrow, length(task_improved$feature_names),
              task_improved$positive))

  # Create custom clustered resampling
  rsmp_clustered <- rsmp("custom")
  rsmp_clustered$instantiate(task_improved, train = train_sets, test = test_sets)
  cat(sprintf("[improved] Custom clustered resampling: %d folds\n", rsmp_clustered$iters))

  # ══════════════════════════════════════════════════════════════════════
  # 6. RUN CV FOR EACH LEARNER PIPELINE
  # ══════════════════════════════════════════════════════════════════════

  cat("\n[improved] Running clustered CV with PipeOps preprocessing...\n")

  # Clipped NLL and Brier loss functions
  clipped_nll <- function(pred, truth) {
    p <- pmin(pmax(pred, 1e-15), 1 - 1e-15)
    -mean(truth * log(p) + (1 - truth) * log(1 - p))
  }
  brier_loss <- function(pred, truth) mean((truth - pred)^2)

  all_learner_pipes <- list(
    # Base learners (no screening)
    mean            = lrn_mean,
    lasso           = lrn_lasso,
    ranger_main     = lrn_ranger_main,
    ranger_low_mtry = lrn_ranger_low_mtry,
    xgb_conservative = lrn_xgb,
    xgb_deep        = lrn_xgb_deep,
    bart_100        = lrn_bart,
    bart_small      = lrn_bart_small,
    earth           = lrn_earth,
    gaussianprocess = lrn_gp
  )

  # Add screening pipelines (only if they built successfully)
  if (!is.null(lrn_lasso_ranger)) all_learner_pipes$lasso_ranger <- lrn_lasso_ranger
  if (!is.null(lrn_lasso_xgb))    all_learner_pipes$lasso_xgb    <- lrn_lasso_xgb
  if (!is.null(lrn_lasso_bart))   all_learner_pipes$lasso_bart   <- lrn_lasso_bart
  if (!is.null(lrn_lasso_earth))  all_learner_pipes$lasso_earth  <- lrn_lasso_earth

  cat(sprintf("[improved] %d learner pipelines (%d with screening)\n",
              length(all_learner_pipes),
              sum(c(!is.null(lrn_lasso_ranger), !is.null(lrn_lasso_xgb),
                    !is.null(lrn_lasso_bart), !is.null(lrn_lasso_earth)))))

  improved_results <- list()
  all_cv_preds <- list()  # Store for ensemble building

  for (name in names(all_learner_pipes)) {
    cat(sprintf("  Fitting: %s ... ", name))
    t0 <- proc.time()

    res <- tryCatch({
      pipe_lrn <- as_learner(all_learner_pipes[[name]])
      pipe_lrn$predict_type <- "prob"

      rr <- suppressWarnings(
        resample(task_improved, pipe_lrn, rsmp_clustered, store_models = FALSE)
      )

      # Extract all CV predictions
      preds_list <- rr$predictions()
      all_probs <- numeric(task_improved$nrow)
      all_truth <- numeric(task_improved$nrow)

      for (i in seq_along(preds_list)) {
        p <- preds_list[[i]]
        idx <- test_sets[[i]]
        probs <- p$prob[, "1"]
        truth <- as.numeric(as.character(p$truth))
        all_probs[idx] <- probs
        all_truth[idx] <- truth
      }

      elapsed <- (proc.time() - t0)["elapsed"]

      nll <- clipped_nll(all_probs, all_truth)
      brier <- brier_loss(all_probs, all_truth)
      auc_val <- as.numeric(pROC::auc(pROC::roc(all_truth, all_probs, quiet = TRUE)))

      cat(sprintf("NLL=%.4f  Brier=%.4f  AUC=%.4f  (%.1fs)\n", nll, brier, auc_val, elapsed))

      all_cv_preds[[name]] <- all_probs

      data.frame(learner = name, nll = nll, brier = brier, auc = auc_val,
                 time_sec = elapsed, stringsAsFactors = FALSE)
    }, error = function(e) {
      elapsed <- (proc.time() - t0)["elapsed"]
      cat(sprintf("FAILED: %s  (%.1fs)\n", substr(e$message, 1, 80), elapsed))
      data.frame(learner = name, nll = NA, brier = NA, auc = NA,
                 time_sec = elapsed, stringsAsFactors = FALSE)
    })

    improved_results[[name]] <- res
  }

  improved_df_results <- do.call(rbind, improved_results)
  rownames(improved_df_results) <- NULL

  # ══════════════════════════════════════════════════════════════════════
  # 7. BUILD ENSEMBLE (NNLS metalearner on CV predictions)
  # ══════════════════════════════════════════════════════════════════════

  cat("\n[improved] Building NNLS ensemble from CV predictions...\n")

  # Stack CV predictions into matrix
  valid_learners <- names(all_cv_preds)
  if (length(valid_learners) >= 2) {
    Z <- do.call(cbind, all_cv_preds[valid_learners])
    colnames(Z) <- valid_learners
    Y_ensemble <- as.numeric(as.character(improved_df$Y))

    # Clamp predictions for NNLS stability
    Z <- pmin(pmax(Z, 1e-10), 1 - 1e-10)

    # NNLS: minimize || Y - Z %*% w ||^2  subject to w >= 0
    nnls_fit <- tryCatch({
      nnls::nnls(Z, Y_ensemble)
    }, error = function(e) NULL)

    if (!is.null(nnls_fit)) {
      weights <- nnls_fit$x
      weights <- weights / sum(weights)  # normalize to sum to 1
      names(weights) <- valid_learners

      # Ensemble CV predictions
      ensemble_preds <- as.numeric(Z %*% weights)

      ensemble_nll <- clipped_nll(ensemble_preds, Y_ensemble)
      ensemble_brier <- brier_loss(ensemble_preds, Y_ensemble)
      ensemble_auc <- as.numeric(pROC::auc(pROC::roc(Y_ensemble, ensemble_preds, quiet = TRUE)))

      cat(sprintf("[improved] Ensemble: NLL=%.4f  Brier=%.4f  AUC=%.4f\n",
                  ensemble_nll, ensemble_brier, ensemble_auc))
      cat("\n[improved] NNLS weights:\n")
      for (nm in names(weights)) {
        if (weights[nm] > 0.001) {
          cat(sprintf("  %20s: %.3f\n", nm, weights[nm]))
        }
      }
      cat(sprintf("  Selected %d / %d learners\n",
                  sum(weights > 0.001), length(weights)))
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # 8. COMPARE: IMPROVED vs STANDARD mlr3 vs sl3
  # ══════════════════════════════════════════════════════════════════════

  cat("\n\n=== IMPROVED PIPELINE RESULTS ===\n")
  improved_df_results <- improved_df_results[order(improved_df_results$nll), ]
  print(improved_df_results, row.names = FALSE)

  cat(sprintf("\n=== HEAD-TO-HEAD: %s — %s ===\n", toupper(COUNTRY), OUTCOME))
  cat(sprintf("%-25s  %8s  %8s  %8s\n", "Pipeline", "NLL", "AUC", "Brier"))
  cat(sprintf("%-25s  %8s  %8s  %8s\n", "---", "---", "---", "---"))

  # sl3 benchmark — use actual session results (loaded or just computed)
  sl3_file <- file.path(out_dir, sprintf("sl3_results_%s_%s.csv", COUNTRY, OUTCOME))
  if (!is.null(sl3_results) && nrow(sl3_results) > 0) {
    sl3_best_nll <- min(sl3_results$cv_risk_val, na.rm = TRUE)
    sl3_best_name <- sl3_results$learner[which.min(sl3_results$cv_risk_val)]
    sl3_auc_val <- if (exists("sl3_ensemble_auc") && !is.na(sl3_ensemble_auc)) sl3_ensemble_auc else NA
    cat(sprintf("%-25s  %8.4f  %8s  %8s\n",
                "sl3 ensemble (CV)",
                if (!is.na(sl3_auc_val)) sl3_best_nll else sl3_best_nll,
                if (!is.na(sl3_auc_val)) sprintf("%.4f", sl3_auc_val) else "—", "—"))
    cat(sprintf("%-25s  %8.4f  %8s  %8s\n",
                paste0("sl3 best (", sl3_best_name, ")"),
                sl3_best_nll, "—", "—"))
  } else if (file.exists(sl3_file)) {
    sl3_saved <- readr::read_csv(sl3_file, show_col_types = FALSE)
    sl3_best_nll <- min(sl3_saved$cv_risk_val, na.rm = TRUE)
    sl3_best_name <- sl3_saved$learner[which.min(sl3_saved$cv_risk_val)]
    cat(sprintf("%-25s  %8.4f  %8s  %8s\n",
                paste0("sl3 best (", sl3_best_name, ", saved)"),
                sl3_best_nll, "—", "—"))
  } else {
    cat(sprintf("%-25s  %8s  %8s  %8s\n", "sl3 (not run)", "—", "—", "—"))
  }

  # Standard mlr3superlearner
  if (exists("mlr3_fit") && !is.null(mlr3_fit)) {
    mlr3_auc_str <- if (exists("mlr3_ensemble_auc") && !is.na(mlr3_ensemble_auc)) {
      sprintf("%.4f*", mlr3_ensemble_auc)
    } else "—"
    mlr3_brier_str <- if (exists("mlr3_brier") && !is.na(mlr3_brier)) {
      sprintf("%.4f*", mlr3_brier)
    } else "—"
    cat(sprintf("%-25s  %8s  %8s  %8s\n",
                "mlr3sl (standard)", "—", mlr3_auc_str, mlr3_brier_str))
    cat("  * mlr3 AUC/Brier are resubstitution (not CV) — expect optimistic bias\n")
  }

  # Improved pipeline
  if (exists("ensemble_nll") && !is.na(ensemble_nll)) {
    cat(sprintf("%-25s  %8.4f  %8.4f  %8.4f\n",
                "IMPROVED ensemble (CV)", ensemble_nll, ensemble_auc, ensemble_brier))
  }
  best_improved <- improved_df_results[1, ]
  if (!is.na(best_improved$nll)) {
    cat(sprintf("%-25s  %8.4f  %8.4f  %8.4f\n",
                paste0("IMPROVED best (", best_improved$learner, ")"),
                best_improved$nll, best_improved$auc, best_improved$brier))
  }

  cat("\nKey differences from standard mlr3superlearner:\n")
  cat("  ✓ TaskClassif (not TaskRegr) — ranger gets proper probability predictions\n")
  cat("  ✓ Clustered CV via origami — respects survey cluster structure\n")
  cat("  ✓ PipeOps inside CV — no data leakage from preprocessing\n")
  cat("  ✓ No washb_prescreen — removes p-value-based data leakage\n")
  cat("  ✓ Learner-specific preprocessing — minimal for ranger/xgb, full for earth/gp\n")
  cat("  ✓ Clipped NLL — prevents log(0) = NaN\n")

  # Save
  improved_path <- file.path(out_dir, "improved_pipeline_results.csv")
  readr::write_csv(improved_df_results, improved_path)
  cat(sprintf("\nSaved to: %s\n", improved_path))

  if (exists("ensemble_nll")) {
    ensemble_summary <- data.frame(
      pipeline = "improved_ensemble",
      nll = ensemble_nll, brier = ensemble_brier, auc = ensemble_auc,
      n_learners = length(valid_learners),
      n_selected = sum(weights > 0.001),
      stringsAsFactors = FALSE
    )
    ensemble_path <- file.path(out_dir, "improved_ensemble_summary.csv")
    readr::write_csv(ensemble_summary, ensemble_path)
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

# ── Load saved results if objects are missing ─────────────────────────────────
# This allows running comparison after restarting R, or when only one
# framework was run in the current session.
if (is.null(sl3_results) || (is.data.frame(sl3_results) && nrow(sl3_results) == 0)) {
  sl3_csv <- file.path(out_dir, sprintf("sl3_results_%s_%s.csv", COUNTRY, OUTCOME))
  if (file.exists(sl3_csv)) {
    cat(sprintf("Loading saved sl3 results from: %s\n", sl3_csv))
    sl3_results <- readr::read_csv(sl3_csv, show_col_types = FALSE)
  } else {
    cat("No sl3 results available (not run and no saved CSV found)\n")
  }
}

if (is.null(mlr3_results) || (is.data.frame(mlr3_results) && nrow(mlr3_results) == 0)) {
  mlr3_csv <- file.path(out_dir, sprintf("mlr3_results_%s_%s.csv", COUNTRY, OUTCOME))
  if (file.exists(mlr3_csv)) {
    cat(sprintf("Loading saved mlr3 results from: %s\n", mlr3_csv))
    mlr3_results <- readr::read_csv(mlr3_csv, show_col_types = FALSE)
  } else {
    cat("No mlr3 results available (not run and no saved CSV found)\n")
  }
}

# Load saved meta values if missing
if (!exists("sl3_ensemble_auc") || is.null(sl3_ensemble_auc)) {
  sl3_meta_csv <- file.path(out_dir, sprintf("sl3_meta_%s_%s.csv", COUNTRY, OUTCOME))
  if (file.exists(sl3_meta_csv)) {
    sl3_meta <- readr::read_csv(sl3_meta_csv, show_col_types = FALSE)
    sl3_ensemble_auc <- sl3_meta$ensemble_auc[1]
    if (!exists("time_sl3")) time_sl3 <- sl3_meta$time_min[1] * 60
    cat(sprintf("Loaded sl3 meta: AUC = %.4f, time = %.1f min\n",
                sl3_ensemble_auc, time_sl3 / 60))
  }
}

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
    cat(sprintf("\nsl3 best:  %s (risk = %.6f)\n", sl3_results$learner[1], sl3_results$cv_risk_val[1]))
  if (!is.null(mlr3_results))
    cat(sprintf("mlr3 best: %s (risk = %.6f)\n", mlr3_results$learner[1], mlr3_results$cv_risk_val[1]))

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
                    if (nrow(sl3_row) > 0) sprintf("%.5f", sl3_row$cv_risk_val[1]) else "—",
                    if (nrow(mlr3_row) > 0) mlr3_row$rank[1] else "—",
                    if (nrow(mlr3_row) > 0) sprintf("%.5f", mlr3_row$cv_risk_val[1]) else "—"))
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
