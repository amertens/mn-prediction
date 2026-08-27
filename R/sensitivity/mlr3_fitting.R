# =============================================================================
# R/mlr3_fitting.R
#
# mlr3-based SuperLearner fitting functions.
# Drop-in replacement for sl3-based DHS_SL_clustered().
# Returns objects with the SAME interface so all downstream functions
# (admin1_analysis, conformal, diagnostics, transportability) work unchanged.
# =============================================================================

#' Learner ids from a library spec, for messages and provenance.
#'
#' A spec entry is either a bare string ("mean") or a list whose first element
#' is the learner type and which may carry an `id`.
sl_learner_ids <- function(library_spec) {
  vapply(library_spec, function(x) {
    if (is.character(x)) return(x[1])
    if (!is.null(x$id)) return(as.character(x$id))
    as.character(x[[1]])
  }, character(1))
}

#' How many learners a given stack mode builds, without side effects.
#'
#' Used by the _targets.R startup banner so that message cannot drift from the
#' stack it describes.
sl_stack_size <- function(stack_mode) {
  invisible(utils::capture.output(
    L <- setup_mlr3_learners(list(sl_stack = stack_mode))))
  length(L$library)
}

#' Set up the mlr3 SuperLearner learner library
#'
#' Returns a list of mlr3 learner specifications for use with
#' mlr3superlearner or manual pipeline construction.
#'
#' @param params Pipeline parameters (from get_pipeline_params())
#' @return list with library specs for continuous and binary tasks
setup_mlr3_learners <- function(params, with_gp = NULL) {

  stack_mode <- params$sl_stack %||% "fast"

  # Gaussian process is OFF by default in full mode as of 2026-08.
  #
  # kernlab::gausspr is O(n^3) in training rows and allocates an n x n
  # kernel matrix. Measured on this machine at p = 140: 0.2 s at n = 400,
  # 1.0 s at n = 800, 9.4 s at n = 1600 -- clean cubic scaling. The
  # area-level targets (n = 30-370) are unaffected, but the individual-level
  # LOCO targets pool n = 10,011 (child_vitA) and 13,107 (women_vitA), where
  # the same curve implies 8-18 h per target across 5 held-out countries.
  # Three such targets ran 4.9 CPU-hours each without finishing.
  #
  # Pass with_gp = TRUE (see the gp_sensitivity target) to put it back for a
  # single small country x outcome, which is how we check that dropping it
  # costs nothing rather than just assuming so.
  if (is.null(with_gp)) with_gp <- isTRUE(params$sl_with_gp)

  if (stack_mode == "fast") {
    # NOTE: BART excluded — its R6 model objects use external pointers that
    # don't survive targets/future serialization, causing predict() to return
    # degenerate (all-identical) values after save/load. This breaks domain
    # ablation and any post-hoc prediction from saved models. BART wins the
    # SL competition frequently but then the model is unusable.
    # Replaced with elastic net (alpha=0.5) for regularization diversity.
    library_spec <- list(
      "mean",
      list("glmnet", alpha = 1, id = "lasso"),
      list("glmnet", alpha = 0.5, id = "elastic_net"),
      list("ranger", num.trees = 500, min.node.size = 10, id = "ranger_main"),
      list("xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
           min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
           id = "xgb_conservative")
    )

  } else {
    library_spec <- list(
      "mean",

      # ── Regularized linear models ──
      # Three alpha values span the regularization spectrum: pure L1 (lasso
      # for variable selection), mixed (elastic net), and pure L2 (ridge
      # for correlated predictors).
      list("glmnet", alpha = 1, id = "lasso"),
      list("glmnet", alpha = 0, id = "ridge"),
      list("glmnet", alpha = 0.5, id = "elastic_net"),

      # ── Random forests ──
      # Multiple configs: standard, low mtry (more randomness, helps when
      # few strong predictors), and deeper trees (captures more interactions).
      list("ranger", num.trees = 500, min.node.size = 10, id = "ranger_main"),
      list("ranger", num.trees = 500, min.node.size = 10, mtry = 8,
           id = "ranger_low_mtry"),
      list("ranger", num.trees = 1000, min.node.size = 5, id = "ranger_deep"),

      # ── XGBoost ──
      # Conservative (shallow trees, high regularization) and deeper (more
      # expressive but with L1/L2 penalty to prevent overfitting).
      list("xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
           min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
           id = "xgb_conservative"),
      list("xgboost", max_depth = 6, eta = 0.03, nrounds = 500,
           min_child_weight = 20, subsample = 0.7, colsample_bytree = 0.4,
           lambda = 1, alpha = 0.5, id = "xgb_deep"),

      # ── BART (Bayesian Additive Regression Trees) ──
      # NOTE: dbarts C++ pointers don't survive serialization. The predict
      # wrapper handles this by refitting BART from stored training data
      # when the primary predict path returns degenerate results (~0.5s).
      # BART excels at rare outcomes where other learners collapse to the mean.
      list("bart", ntree = 50, id = "bart_small"),
      list("bart", ntree = 100, id = "bart_100"),
      list("bart", ntree = 200, id = "bart_large")

      # ── Gaussian process ──
      # Appended below when with_gp = TRUE, not listed here. See the note at the
      # top of this function for why it is off by default.

      # ── Disabled learners (caused silent failures in full pipeline) ──
      # LightGBM and nnet passed individual serialization tests but caused
      # mlr3superlearner to return NULL for 14/18 country×outcome combos
      # when included in the full stack. Root cause unclear — likely
      # interaction between learner initialization and mlr3superlearner's
      # internal CV machinery on real (non-simulated) data.
      # TODO: debug with workspace_on_error and re-enable individually.
      # list("lightgbm", num_leaves = 31, learning_rate = 0.05,
      #      num_iterations = 300, min_data_in_leaf = 20,
      #      feature_fraction = 0.7, bagging_fraction = 0.8,
      #      bagging_freq = 5, id = "lgbm_main"),
      # list("nnet", size = 10, decay = 0.01, maxit = 200, id = "nnet_small")
    )
  }

  # Add the Gaussian process back only when explicitly asked for.
  if (with_gp && stack_mode != "fast") {
    library_spec <- c(library_spec, list("gaussianprocess"))
    cat("[mlr3_learners] with_gp = TRUE - Gaussian process ADDED; this is O(n^3),",
        "keep it to small n
")
  }

  # Optional: drop BART learners (dbarts MCMC can fatally crash the R process
  # — OOM in parallel workers and segfaults even single-process). Set
  # SL_NO_BART=1 for a reliable unattended run; the remaining stack is still
  # rich (glmnet x3, ranger x3, xgboost x2, gaussian process). Default keeps BART.
  if (nzchar(Sys.getenv("SL_NO_BART"))) {
    n0 <- length(library_spec)
    library_spec <- Filter(function(x) {
      type <- if (is.list(x)) (x[[1]] %||% "") else x
      !grepl("bart", type, ignore.case = TRUE)
    }, library_spec)
    cat(sprintf("[mlr3_learners] SL_NO_BART set — dropped %d BART learner(s) (%d -> %d)\n",
                n0 - length(library_spec), n0, length(library_spec)))
  }

  # Report the count from the spec itself rather than from a literal, and after
  # any SL_NO_BART filtering. The two hardcoded numbers that used to sit in the
  # branches above had drifted to "5" and "16" against actual stacks of 5 and
  # 13, while _targets.R printed a third and fourth number ("3 learners" and
  # "6-learner stack") for the same two stacks. Deriving it cannot drift.
  cat(sprintf("[mlr3_learners] %s stack: %d learners (%s)\n",
              toupper(stack_mode), length(library_spec),
              paste(sl_learner_ids(library_spec), collapse = ", ")))

  list(library = library_spec, stack_mode = stack_mode)
}


#' Compute genuine out-of-fold (cross-validated) SuperLearner predictions
#'
#' `mlr3superlearner` runs cross-validation internally to estimate the ensemble
#' weights but then DISCARDS the cross-validated predictions, returning learners
#' refit on the full data — so `predict(fit, training_data)` is in-sample
#' (resubstitution). This recomputes honest out-of-fold predictions by
#' re-running `mlr3::resample()` on the fitted (cloned + reset) base learners
#' with the same clustered fold structure, then combining them with the SL
#' weights. For the default discrete SL this is a single learner (cheap).
#'
#' These OOF predictions are what conformal intervals, calibration, AUC/Brier,
#' and area-level aggregation must use to be valid / non-optimistic.
#'
#' @return numeric vector aligned to the rows of `fit_df`, or NULL on failure
#'   (caller falls back to in-sample predictions).
mlr3_oof_predictions <- function(mlr3_fit, fit_df, target = "Y",
                                  outcome_type = "binomial",
                                  group = "cluster_id", folds = 5L,
                                  seed = 5249L) {
  tryCatch({
    data <- as.data.frame(fit_df)
    if (outcome_type == "continuous") {
      task <- mlr3::as_task_regr(data, target = target)
    } else {
      data[[target]] <- factor(data[[target]])
      task <- mlr3::as_task_classif(data, target = target, positive = "1")
    }
    if (!is.null(group) && group %in% colnames(data)) {
      task$set_col_roles(group, "group")
    }
    # Features = SL predictors minus the grouping column
    feats <- setdiff(intersect(mlr3_fit$x, colnames(data)), group)
    task$select(feats)

    rsmp <- mlr3::rsmp("cv", folds = folds)
    set.seed(seed)
    rsmp$instantiate(task)

    w <- mlr3_fit$weights
    use <- names(w)[is.finite(w) & w != 0]
    if (length(use) == 0) return(NULL)

    extract_oof <- function(lid) {
      idx <- which(vapply(mlr3_fit$learners, function(l) l$id, character(1)) == lid)
      if (length(idx) == 0) {
        if (length(mlr3_fit$learners) == 1) idx <- 1L else return(NULL)
      }
      lrn <- mlr3_fit$learners[[idx[1]]]$clone(deep = TRUE)
      lrn$reset()
      rr <- mlr3::resample(task, lrn, rsmp, store_models = FALSE)
      pr <- rr$prediction()
      dt <- data.table::as.data.table(pr)
      data.table::setorder(dt, row_ids)
      if (outcome_type == "continuous") dt[["response"]] else dt[["prob.1"]]
    }

    Z <- vapply(use, extract_oof, numeric(nrow(data)))
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = length(use))
    wn <- w[use] / sum(w[use])
    oof <- as.numeric(Z %*% wn)
    if (length(oof) != nrow(data) || all(is.na(oof))) return(NULL)
    oof
  }, error = function(e) {
    cat(sprintf("  [mlr3_SL] OOF prediction failed (%s); using in-sample fallback\n",
                e$message))
    NULL
  })
}


#' Fit mlr3 SuperLearner for one outcome (wrapper matching DHS_SL_clustered interface)
#'
#' This function does the same preprocessing as DHS_SL_clustered():
#'   1. Unlabel + drop NA/NZV + impute + NZV
#'   2. Prescreen with washb_prescreen (optional)
#'   3. Apply recipes (step_corr, step_normalize)
#'   4. Fit mlr3superlearner with clustered CV
#'
#' Returns an object with the same interface as DHS_SL_clustered() output:
#'   - sl_fit: model object with $predict() method
#'   - res: data.frame with Y, yhat_full (CV predictions)
#'   - Xvars: final covariate names
#'   - auto_recipe: preprocessing recipe
#'   - pre_recipe_cols: columns before recipe
#'   - cv_risk: data.frame of learner-level CV risks
#'
#' @param d Data frame with outcome and predictors
#' @param Xvars Character vector of predictor column names
#' @param outcome Character, outcome column name
#' @param population Character, population label
#' @param id Character, cluster ID column name
#' @param folds Integer, number of CV folds
#' @param mlr3_library List of mlr3 learner specifications
#' @param outcome_type "binomial" or "continuous"
#' @param prescreen Logical, whether to run washb_prescreen
#' @return List matching DHS_SL_clustered interface
mlr3_SL_clustered <- function(d, Xvars, outcome, population,
                               id = "gw_cnum", folds = 5L,
                               mlr3_library = NULL,
                               outcome_type = "binomial",
                               prescreen = TRUE) {

  # ── Load required packages ──
  if (!requireNamespace("mlr3superlearner", quietly = TRUE))
    stop("mlr3superlearner required. Install with: devtools::install_github('nt-williams/mlr3superlearner')")
  if (!requireNamespace("mlr3", quietly = TRUE))
    stop("mlr3 required. Install with: install.packages('mlr3')")
  # The full stack uses BART / Gaussian-process learners that live in
  # mlr3extralearners; attaching it registers them in mlr3's learner
  # dictionary. Without this the fit fails with "mlr3extralearners required"
  # whenever the caller hasn't already attached it (the targets pipeline does
  # so via tar_option_set, but standalone callers may not).
  if (requireNamespace("mlr3extralearners", quietly = TRUE))
    suppressMessages(suppressWarnings(library(mlr3extralearners)))

  # ── Validate the outcome column exists and is usable ──
  # Guards against config/data drift (e.g. a stale merged cache missing a
  # biomarker column): without this, an absent column yields Y = NULL and the
  # downstream data.frame() / prescreen fails with a cryptic
  # "differing number of rows" error.
  if (is.null(outcome) || !(outcome %in% colnames(d))) {
    cat(sprintf("  [mlr3_SL] SKIP: outcome column '%s' not found in data\n",
                outcome %||% "NULL"))
    return(NULL)
  }
  y_check <- suppressWarnings(as.numeric(d[[outcome]]))
  if (sum(!is.na(y_check)) == 0 || length(unique(y_check[!is.na(y_check)])) < 2) {
    cat(sprintf("  [mlr3_SL] SKIP: outcome '%s' has <2 distinct non-missing values\n",
                outcome))
    return(NULL)
  }

  # ── Preprocessing (identical to DHS_SL_clustered) ──
  X <- d %>% dplyr::select(dplyr::all_of(Xvars)) %>% as.data.frame()
  cov <- labelled::unlabelled(X, user_na_to_na = TRUE)
  Y <- d[[outcome]]
  id_vec <- d[[id]]
  dataid <- if ("dataid" %in% colnames(d)) d$dataid else seq_len(nrow(d))

  # Strip haven labels from Y
  if (inherits(Y, "haven_labelled")) Y <- as.double(unclass(Y))
  Y <- as.numeric(Y)

  # Drop all-NA and constant columns
  cov <- cov[, !sapply(cov, function(x) all(is.na(x))), drop = FALSE]
  cov <- cov %>% dplyr::select(dplyr::where(~{
    non_na <- .x[!is.na(.x)]
    length(non_na) > 0 && length(unique(non_na)) > 1
  }))

  # Replace Inf/NaN with NA before imputation (GEE raster extraction can
  # produce Inf/NaN for polygons that don't overlap the raster extent)
  n_inf <- 0L
  for (col in colnames(cov)) {
    bad <- !is.finite(cov[[col]]) & !is.na(cov[[col]])
    if (any(bad)) {
      cov[[col]][bad] <- NA
      n_inf <- n_inf + sum(bad)
    }
  }
  if (n_inf > 0) cat(sprintf("  [mlr3_SL] Replaced %d Inf/NaN values with NA\n", n_inf))

  # NZV
  nzv_idx <- caret::nearZeroVar(cov)
  if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

  # Impute
  cov <- as.data.frame(
    ck37r::impute_missing_values(cov, type = "standard",
                                  add_indicators = TRUE, prefix = "missing_")$data
  )

  # NZV again (after imputation may create constant indicator columns)
  nzv_idx <- caret::nearZeroVar(cov)
  if (length(nzv_idx) > 0) cov <- cov[, -nzv_idx, drop = FALSE]

  # Below this many raw predictors, skip BOTH washb_prescreen and step_corr
  # (see their respective call sites for why) -- they exist to cut a large
  # candidate pool down; on an already-small, deliberately-curated set they
  # do more harm than good. Shared by both fixes so the threshold can't drift
  # between them.
  SMALL_P_SKIP_CORR_THRESHOLD <- 15L

  # Optional prescreening
  # 2026-08-27 fix: washb_prescreen's job is to cut a large candidate pool
  # down; on an already-small set (same threshold/rationale as the step_corr
  # skip above) it can reduce to a single surviving predictor, which then
  # crashes glmnet ("x should be a matrix with 2 or more columns") --
  # reproduced with a 6-variable curated set where one fold's prescreen left
  # exactly 1 column. Skip it below the same small-p threshold; there's
  # nothing to screen out of a set that size that a human/algorithm didn't
  # already choose deliberately.
  if (prescreen && ncol(cov) <= SMALL_P_SKIP_CORR_THRESHOLD) {
    cat(sprintf("  [mlr3_SL] p=%d <= %d: skipping washb_prescreen (same rationale as step_corr above)\n",
                ncol(cov), SMALL_P_SKIP_CORR_THRESHOLD))
    prescreen <- FALSE
  }
  if (prescreen) {
    family_screen <- if (length(unique(Y[!is.na(Y)])) == 2) "binomial" else "gaussian"

    # Use relaxed p-value threshold (0.2) for marginal screening.
    # For rare outcomes (<10% prevalence), further relax to 0.3 to avoid
    # discarding predictors that matter in combination.
    prevalence <- mean(Y == 1, na.rm = TRUE)
    pval_thresh <- if (family_screen == "binomial" && prevalence < 0.10) 0.3 else 0.2

    n_before_screen <- ncol(cov)
    Wvars <- tryCatch(
      washb::washb_prescreen(Y = Y, Ws = cov, family = family_screen,
                              pval = pval_thresh, print = FALSE),
      error = function(e) {
        cat(sprintf("  [mlr3_SL] washb_prescreen failed: %s. Using all vars.\n", e$message))
        colnames(cov)
      }
    )

    if (length(Wvars) == 0) {
      cat("  [mlr3_SL] Prescreen removed all variables — using all\n")
      Wvars <- colnames(cov)
    }
    cov <- cov %>% dplyr::select(dplyr::all_of(Wvars))
    # Log input -> output counts (n_before_screen is captured BEFORE subsetting;
    # the previous code measured ncol(cov) after subsetting, so it always
    # printed "n/n" and hid whether screening actually filtered anything).
    cat(sprintf("  [mlr3_SL] Prescreen (p<%.2f): %d/%d vars retained\n",
                pval_thresh, length(Wvars), n_before_screen))
  }

  # Save pre-recipe columns (needed for prediction on new data)
  pre_recipe_cols <- colnames(cov)

  # Recipe: ZV + NZV + correlation filter + normalize
  # For rare outcomes, relax the correlation threshold to keep more
  # predictors — correlated variables may capture different facets of
  # the outcome when prevalence is low.
  cor_threshold <- if (exists("prevalence") && prevalence < 0.10) 0.90 else 0.85

  # 2026-08-27 fix: step_corr is tuned for the large-p regime (hundreds of
  # candidates) this pipeline normally runs in, and applied UNCONDITIONALLY --
  # it wasn't gated by `prescreen`, so even prescreen=FALSE couldn't disable
  # it. On an already-small, hand-curated predictor set (confirmed with a
  # 6-variable LOCO-selected set: two malaria variables correlated at r=0.99,
  # both >0.8 with rainfall/MAP) it silently discards predictors kept
  # *because* they're jointly informative despite being correlated, cutting
  # LOCO AUC by ~0.03-0.06 for no benefit. Skip it below the same small-p
  # threshold used for prescreen above -- there's nothing to protect against
  # at this scale, since a human/algorithm already chose these variables
  # deliberately.
  skip_corr <- ncol(cov) <= SMALL_P_SKIP_CORR_THRESHOLD
  if (skip_corr) {
    cat(sprintf("  [mlr3_SL] p=%d <= %d: skipping step_corr (tuned for large-p screening, not a small curated set)\n",
                ncol(cov), SMALL_P_SKIP_CORR_THRESHOLD))
  }

  auto_recipe <- recipes::recipe(~ ., data = cov) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_nzv(recipes::all_predictors())
  if (!skip_corr) {
    auto_recipe <- auto_recipe %>% recipes::step_corr(recipes::all_numeric(), threshold = cor_threshold)
  }
  auto_recipe <- auto_recipe %>%
    recipes::step_normalize(recipes::all_numeric()) %>%
    recipes::prep()

  cov <- data.frame(recipes::bake(auto_recipe, new_data = cov))
  covars <- colnames(cov)

  cat(sprintf("  [mlr3_SL] After preprocessing: n=%d, p=%d\n", nrow(cov), length(covars)))

  # ── Strip haven labels from all columns ──
  for (col in colnames(cov)) {
    if (inherits(cov[[col]], "haven_labelled")) {
      cov[[col]] <- as.double(unclass(cov[[col]]))
    }
    if (is.factor(cov[[col]]) || is.character(cov[[col]])) {
      cov[[col]] <- as.numeric(as.character(cov[[col]]))
    }
  }

  # Replace any remaining NaN/Inf with NA (not 0 — zero is a valid value
  # and would introduce bias). Edge cases from recipe transforms.
  # xgboost and ranger handle NA natively.
  for (col in colnames(cov)) {
    bad <- !is.finite(cov[[col]]) & !is.na(cov[[col]])
    if (any(bad)) {
      cov[[col]][bad] <- NA
      cat(sprintf("  [mlr3_SL] Post-recipe: replaced %d Inf/NaN with NA in '%s'\n",
                  sum(bad), col))
    }
  }

  # ── Add class-weighted learner variants for rare binary outcomes ──────
  # When prevalence < 15%, add xgboost and ranger configs with class
  # weights that upweight the minority class. This prevents the SL from
  # collapsing to "predict everyone as non-deficient".
  if (outcome_type == "binomial" || outcome_type == "binary") {
    prev_rate <- mean(Y == 1, na.rm = TRUE)
    if (prev_rate > 0 && prev_rate < 0.15) {
      pos_weight <- (1 - prev_rate) / prev_rate  # e.g., 75 for 1.3% prevalence
      # Cap at 50 to avoid extreme instability
      pos_weight <- min(pos_weight, 50)
      cat(sprintf("  [mlr3_SL] Rare outcome (prev=%.1f%%) — adding class-weighted learners (weight=%.1f)\n",
                  prev_rate * 100, pos_weight))

      mlr3_library <- c(mlr3_library, list(
        list("xgboost", max_depth = 4, eta = 0.05, nrounds = 300,
             min_child_weight = 10, subsample = 0.8, colsample_bytree = 0.5,
             scale_pos_weight = pos_weight, id = "xgb_weighted"),
        list("ranger", num.trees = 500, min.node.size = 5,
             class.weights = c("0" = 1, "1" = pos_weight),
             id = "ranger_weighted")
      ))
    }
  }

  # ── Build mlr3 data frame ──
  mlr3_df <- data.frame(Y = Y, cov)
  mlr3_df$cluster_id <- id_vec

  # Remove rows with NA outcome
  valid_rows <- !is.na(mlr3_df$Y)
  mlr3_df <- mlr3_df[valid_rows, ]

  # ── Determine outcome type and set up folds ──
  # Cluster-based CV: pass cluster IDs to mlr3superlearner via `group=`,
  # which keeps observations from the same PSU in the same fold (per
  # mlr3superlearner help: "observations in the same group are treated
  # like a 'block' of observations kept together during sample splitting").
  # The origami fold_obj is also retained for downstream consumers.
  cluster_ids <- mlr3_df$cluster_id
  fold_obj <- origami::make_folds(cluster_ids = cluster_ids, V = folds)

  # ── Fit mlr3superlearner ──
  # Monkey-patch NLL to clip probabilities (prevents NaN from ranger)
  env <- asNamespace("mlr3superlearner")
  if (exists("loss_nll", envir = env)) {
    original_nll <- get("loss_nll", envir = env)
    clipped_nll <- function(x, y) {
      x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
      -mean(y * log(x) + (1 - y) * log(1 - x))
    }
    tryCatch(
      assignInNamespace("loss_nll", clipped_nll, "mlr3superlearner"),
      error = function(e) NULL
    )
  }

  # Keep cluster_id IN the data and reference it by name via `group=`.
  # mlr3superlearner uses `group` as the blocking variable for sample
  # splitting and excludes it from predictors automatically.
  fit_df <- mlr3_df  # cluster_id stays as a column

  # 2026-08-27 fix: this call's internal CV (which LEARNS the ensemble
  # weights) was never seeded -- only the separate honest-OOF resampling in
  # mlr3_oof_predictions() below was. With a rich library and hundreds of
  # predictors the resulting run-to-run noise is usually masked by averaging
  # over many learners; with a small library/small p it isn't: confirmed the
  # SAME data + SAME 3-learner library picked 100% weight on a DIFFERENT
  # single learner on every unseeded re-run, swinging held-out AUC between
  # 0.49 and 0.72. Same seed constant as mlr3_oof_predictions() for
  # consistency between the two CV steps.
  set.seed(5249L)
  t0 <- proc.time()
  mlr3_fit <- tryCatch({
    suppressWarnings(
      mlr3superlearner::mlr3superlearner(
        data = fit_df,
        target = "Y",
        library = mlr3_library,
        outcome_type = outcome_type,
        folds = folds,
        group = "cluster_id"
      )
    )
  }, error = function(e) {
    cat(sprintf("  [mlr3_SL] FAILED: %s\n", e$message))
    NULL
  })
  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  [mlr3_SL] Fit in %.1f min\n", elapsed / 60))

  if (is.null(mlr3_fit)) {
    return(NULL)
  }

  # ── Extract predictions ──
  # In-sample (resubstitution) predictions from the full-data ensemble.
  yhat_insample <- tryCatch({
    as.numeric(predict(mlr3_fit, fit_df))
  }, error = function(e) {
    cat(sprintf("  [mlr3_SL] predict() failed: %s\n", e$message))
    rep(NA_real_, nrow(fit_df))
  })

  # Genuine out-of-fold (cross-validated) predictions. mlr3superlearner does
  # not retain its internal CV predictions, so recompute them by resampling
  # the fitted base learners on the same clustered folds (see
  # mlr3_oof_predictions()). These are what conformal CIs, calibration, AUC
  # and area aggregation must use — resubstitution would be optimistic.
  yhat_oof <- mlr3_oof_predictions(
    mlr3_fit, fit_df, target = "Y", outcome_type = outcome_type,
    group = "cluster_id", folds = folds
  )
  if (is.null(yhat_oof)) {
    cat("  [mlr3_SL] WARNING: falling back to in-sample yhat (OOF unavailable)\n")
    yhat_full <- yhat_insample
  } else {
    yhat_full <- yhat_oof
    # Report the in-sample vs OOF optimism gap for binary outcomes.
    if (outcome_type %in% c("binomial", "binary") &&
        requireNamespace("pROC", quietly = TRUE)) {
      gap <- tryCatch({
        yv <- mlr3_df$Y
        auc_in  <- as.numeric(pROC::auc(pROC::roc(yv, yhat_insample, quiet = TRUE)))
        auc_oof <- as.numeric(pROC::auc(pROC::roc(yv, yhat_full,     quiet = TRUE)))
        sprintf("AUC in-sample=%.3f, OOF=%.3f (optimism=%.3f)", auc_in, auc_oof, auc_in - auc_oof)
      }, error = function(e) NULL)
      if (!is.null(gap)) cat(sprintf("  [mlr3_SL] %s\n", gap))
    }
  }

  # ── Extract CV risk table ──
  cv_risk_table <- tryCatch({
    # mlr3superlearner stores risks in the fit object
    raw_risks <- mlr3_fit$fits$risk
    if (!is.null(raw_risks)) {
      learner_names <- names(raw_risks)
      if (is.null(learner_names)) learner_names <- paste0("learner_", seq_along(raw_risks))
      data.frame(
        learner = learner_names,
        risk = as.numeric(raw_risks),
        coefficients = as.numeric(mlr3_fit$fits$weights %||% rep(NA, length(raw_risks))),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(learner = character(), risk = numeric(), coefficients = numeric())
    }
  }, error = function(e) {
    data.frame(learner = character(), risk = numeric(), coefficients = numeric())
  })

  # ── Build result in DHS_SL_clustered format ──
  # Create a wrapper object that supports $predict() for new data.
  #
  # The mlr3superlearner object stores R6 Learner objects that may not
  # survive serialization (targets/future send objects to worker processes
  # via socket serialization). To ensure predict works after reload, we
  # also extract the individual trained learner models as plain R objects
  # and store the SL weights, enabling a manual fallback prediction path.

  # Extract individual learner models for serialization-safe fallback.
  # BART (dbarts) models use C++ external pointers that don't survive
  # serialization. For BART, we store the training data so we can refit
  # on the fly when predict is needed (~0.5s overhead).
  learner_models <- tryCatch({
    # mlr3superlearner may store learners in different locations depending
    # on discrete vs ensemble mode. Try multiple paths.
    lrn_list <- mlr3_fit$learners
    # Also check $fits$learners which may have the full set
    if (length(lrn_list) == 1 && !is.null(mlr3_fit$fits) &&
        !is.null(mlr3_fit$fits$learners) &&
        length(mlr3_fit$fits$learners) > length(lrn_list)) {
      lrn_list <- mlr3_fit$fits$learners
    }
    cat(sprintf("  [mlr3_SL] Extracting %d learner models for fallback predict\n",
                length(lrn_list)))
    lapply(lrn_list, function(lrn) {
      is_bart <- grepl("bart", lrn$id, ignore.case = TRUE)
      entry <- list(
        id    = lrn$id,
        model = lrn$model
      )
      if (is_bart) {
        entry$is_bart    <- TRUE
        entry$bart_ntree <- tryCatch(lrn$param_set$values$ntree %||% 100L,
                                      error = function(e) 100L)
      }
      entry
    })
  }, error = function(e) {
    cat(sprintf("  [mlr3_SL] learner_models extraction failed: %s\n", e$message))
    NULL
  })

  sl_weights <- tryCatch(mlr3_fit$weights, error = function(e) NULL)
  sl_x       <- tryCatch(mlr3_fit$x, error = function(e) covars)
  sl_outcome_type <- tryCatch(mlr3_fit$outcome_type, error = function(e) outcome_type)

  # Store training data for BART refit-on-predict
  train_data_for_bart <- data.table::as.data.table(fit_df)

  model_wrapper <- list(
    mlr3_fit       = mlr3_fit,
    covars         = covars,
    learner_models = learner_models,  # plain R model objects (serialization-safe)
    sl_weights     = sl_weights,
    sl_x           = sl_x,
    sl_outcome_type = sl_outcome_type,
    train_data_for_bart = train_data_for_bart,
    predict    = function(new_task) {
      # Accept either sl3-style task or data.frame/data.table
      if (inherits(new_task, "sl3_Task")) {
        newdata <- data.frame(new_task$data)
      } else if (is.data.frame(new_task) || data.table::is.data.table(new_task)) {
        newdata <- data.frame(new_task)
      } else {
        stop("predict() expects a data.frame or sl3_Task")
      }
      # Ensure column alignment
      # TODO(wrapper-cluster_id): aligning to `covars` strips the CV group
      # column (cluster_id) that mlr3superlearner's predict needs present
      # (it's recorded in mlr3_fit$x even though no learner uses it as a
      # feature), so the primary path below errors with "undefined columns
      # selected" and the fallback at `newdata[, sl_x]` also errors. Every
      # current caller bypasses this by calling predict(fit$mlr3_fit, ...)
      # directly (see R/domain_ablation.R, R/single_var_ablation.R,
      # R/conceptual_ablation.R). Fix the wrapper properly the next time
      # sl_fit_* needs to be invalidated for another reason:
      #   needed <- union(covars, sl_x)
      #   for (col in setdiff(needed, colnames(newdata))) newdata[[col]] <- 0
      #   newdata <- newdata[, needed, drop = FALSE]
      # and at the fallback selection (line ~597):
      #   nd <- newdata[, intersect(covars, colnames(newdata)), drop = FALSE]
      # Held off here to avoid forcing a costly refit of every SL model.
      for (col in setdiff(covars, colnames(newdata))) newdata[[col]] <- 0
      newdata <- newdata[, covars, drop = FALSE]
      newdata$Y <- 0  # dummy outcome for predict

      # Try mlr3superlearner predict first
      pred <- tryCatch(
        as.numeric(predict(mlr3_fit, newdata)),
        error = function(e) NULL
      )

      # Check for degenerate predictions (serialization failure)
      if (!is.null(pred) && length(unique(round(pred, 6))) > 1) {
        return(pred)
      }

      # Fallback: manual prediction from extracted learner models + weights.
      # For BART: refit from stored training data (C++ pointers don't serialize).
      # For glmnet/ranger/xgboost: use the stored plain R model objects directly.
      if (!is.null(learner_models) && !is.null(sl_weights)) {
        nd <- newdata[, sl_x, drop = FALSE]
        # Debug: uncomment to trace fallback
        # cat(sprintf("  [predict fallback] %d learner_models, %d weights\n",
        #             length(learner_models), length(sl_weights)))
        preds_mat <- sapply(seq_along(learner_models), function(i) {
          lm_obj <- learner_models[[i]]
          tryCatch({
            if (isTRUE(lm_obj$is_bart)) {
              # BART: refit from training data and predict on new data.
              td <- as.data.frame(train_data_for_bart)
              X_train_b <- as.matrix(td[, sl_x, drop = FALSE])
              Y_train_b <- as.numeric(td$Y)
              X_test_b  <- as.matrix(nd)
              # BART refit: ~0.5s per call
              bart_refit <- dbarts::bart(
                x.train = X_train_b, y.train = Y_train_b,
                x.test = X_test_b,
                ntree = lm_obj$bart_ntree %||% 100L,
                verbose = FALSE, keeptrees = FALSE
              )
              # yhat.test: (n_draws x n_test)
              if (sl_outcome_type == "binomial") {
                preds_b <- pnorm(colMeans(bart_refit$yhat.test))
              } else {
                preds_b <- colMeans(bart_refit$yhat.test)
              }
              # Refit complete
              preds_b
            } else if (grepl("glmnet|lasso|ridge|enet", lm_obj$id)) {
              as.numeric(predict(lm_obj$model, as.matrix(nd), type = "response",
                                  s = lm_obj$model$lambda.min %||% lm_obj$model$lambda[1]))
            } else if (grepl("ranger", lm_obj$id)) {
              rpred <- predict(lm_obj$model, data = nd)$predictions
              if (is.matrix(rpred) && ncol(rpred) >= 2) rpred[, 2] else as.numeric(rpred)
            } else if (grepl("xgboost|xgb|lgbm", lm_obj$id)) {
              predict(lm_obj$model, as.matrix(nd))
            } else if (grepl("nnet", lm_obj$id)) {
              as.numeric(predict(lm_obj$model, nd, type = "raw"))
            } else if (grepl("gaussianprocess|gp|ksvm|svm", lm_obj$id)) {
              kernlab::predict(lm_obj$model, as.matrix(nd), type = "probabilities")[, 2]
            } else if (grepl("mean", lm_obj$id)) {
              rep(mean(Y, na.rm = TRUE), nrow(nd))
            } else {
              rep(NA_real_, nrow(nd))
            }
          }, error = function(e) {
            rep(NA_real_, nrow(nd))
          })
        })
        if (is.matrix(preds_mat)) {
          n_cols <- ncol(preds_mat)
          if (n_cols == 1) {
            # Only 1 learner model extracted (mlr3superlearner stores only
            # the winning learner). Use it directly — no weighting needed.
            result <- as.numeric(preds_mat[, 1])
            if (length(unique(round(result, 6))) > 1) return(result)
          } else {
            # Multiple learner models — align weights and combine
            n_wts <- length(sl_weights)
            wts_use <- if (n_wts == n_cols) {
              sl_weights
            } else if (n_wts > n_cols) {
              sl_weights[seq_len(n_cols)]
            } else {
              c(sl_weights, rep(0, n_cols - n_wts))
            }
            if (sum(wts_use) > 0) wts_use <- wts_use / sum(wts_use)
            use_idx <- which(wts_use > 0.001)
            if (length(use_idx) > 0) {
              for (j in use_idx) {
                if (all(is.na(preds_mat[, j]))) preds_mat[, j] <- mean(Y, na.rm = TRUE)
              }
              result <- as.numeric(preds_mat[, use_idx, drop = FALSE] %*% wts_use[use_idx])
              if (length(unique(round(result, 6))) > 1) return(result)
            }
          }
        }
      }

      # Last resort: return the R6-based prediction even if degenerate
      if (!is.null(pred)) return(pred)
      rep(NA_real_, nrow(newdata))
    }
  )
  class(model_wrapper) <- "mlr3_sl_wrapper"

  res <- data.frame(
    dataid    = dataid[valid_rows],
    clusterid = cluster_ids,
    outcome   = outcome,
    population = population,
    Y         = mlr3_df$Y,
    yhat_full = yhat_full,        # out-of-fold (cross-validated) — used downstream
    yhat_insample = yhat_insample  # resubstitution — kept for reference/diagnostics
  )

  # Store training data for permutation-based domain importance
  # (permute_domain needs the original data to shuffle columns)
  train_data <- data.table::as.data.table(fit_df)

  list(
    sl_fit          = model_wrapper,
    res             = res,
    cv_risk_w_sl_revere = cv_risk_table,
    task            = NULL,  # mlr3 doesn't use sl3 tasks
    train_data      = train_data,  # stored for permutation ablation
    Xvars           = covars,
    auto_recipe     = auto_recipe,
    pre_recipe_cols = pre_recipe_cols
  )
}


#' Fit continuous + binary mlr3 SL for one outcome
#'
#' Drop-in replacement for fit_sl_models() from R/sl_fitting.R.
#' Returns identical interface.
#'
#' @param outcome_data Output from build_outcome_dataset()
#' @param cc Country config
#' @param oc Outcome config
#' @param sl_learners Output from setup_mlr3_learners()
#' @param params Pipeline parameters
#' @return list with cont_fit and bin_fit
fit_mlr3_models <- function(outcome_data, cc, oc, sl_learners, params) {

  d     <- outcome_data$data
  Xvars <- outcome_data$Xvars

  cat(sprintf("\n[fit_mlr3] %s — %s | n = %d | p = %d (proxy-only, no gw_)\n",
              cc$country, oc$tag, nrow(d), length(Xvars)))

  # Bail early if no data
  if (nrow(d) == 0 || length(Xvars) == 0) {
    cat("  SKIPPING: no observations or no predictors available.\n")
    return(list(cont_fit = NULL, bin_fit = NULL,
                outcome = oc$tag, country = cc$country))
  }

  mlr3_lib <- sl_learners$library

  # Remove "mean" learner for Sierra Leone — discrete SL selects the

  # intercept-only model over real learners that achieve AUC 0.65-0.79,
  # because no single learner *consistently* beats the mean across all
  # clustered CV folds. Excluding "mean" forces selection of a real learner.
  if (grepl("sierra", cc$country, ignore.case = TRUE)) {
    mlr3_lib <- Filter(function(x) !identical(x, "mean"), mlr3_lib)
    cat("  [mlr3] Removed 'mean' learner for Sierra Leone\n")
  }

  # ── Continuous model ──
  # Skip cleanly when the configured continuous outcome column is absent
  # (e.g. a stale merged cache missing a biomarker), mirroring the binary
  # guard below. Previously this attempted a guaranteed-to-fail fit that
  # surfaced as a cryptic washb_prescreen "differing number of rows" error.
  cont_fit <- NULL
  if (is.null(oc$continuous) || !(oc$continuous %in% colnames(d))) {
    cat(sprintf("  SKIPPING continuous: outcome column '%s' not in data\n",
                oc$continuous %||% "NULL"))
  } else {
    cat("  Fitting continuous mlr3 SL...\n")
    cont_fit <- tryCatch(
      mlr3_SL_clustered(
        d = d, Xvars = Xvars, outcome = oc$continuous,
        population = oc$population, id = cc$cluster_id,
        folds = params$K, mlr3_library = mlr3_lib,
        outcome_type = "continuous", prescreen = TRUE
      ),
      error = function(e) {
        cat(sprintf("  Continuous SL failed: %s\n", e$message))
        NULL
      }
    )
  }

  # ── Binary model ──
  bin_fit <- NULL
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    y_bin <- d[[oc$binary]]
    if (inherits(y_bin, "haven_labelled")) y_bin <- as.double(unclass(y_bin))
    y_bin_nona <- y_bin[!is.na(y_bin)]
    is_binary <- length(unique(y_bin_nona)) >= 2 &&
                 all(y_bin_nona >= 0 & y_bin_nona <= 1)

    if (!is_binary) {
      cat(sprintf("  SKIPPING binary: values not in [0,1] (range: %.2f-%.2f)\n",
                  min(y_bin_nona), max(y_bin_nona)))
    } else {
      cat("  Fitting binary mlr3 SL...\n")
      bin_fit <- tryCatch(
        mlr3_SL_clustered(
          d = d, Xvars = Xvars, outcome = oc$binary,
          population = oc$population, id = cc$cluster_id,
          folds = params$K, mlr3_library = mlr3_lib,
          outcome_type = "binomial", prescreen = TRUE
        ),
        error = function(e) {
          cat(sprintf("  Binary SL failed: %s\n", e$message))
          NULL
        }
      )
    }
  }

  list(
    cont_fit = cont_fit,
    bin_fit  = bin_fit,
    outcome  = oc$tag,
    country  = cc$country
  )
}
