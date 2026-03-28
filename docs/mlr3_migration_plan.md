# Plan: Migrate SuperLearner Pipeline from sl3 to mlr3superlearner

## Motivation

- **Speed**: mlr3 runs 40-50x faster (3 min vs 130 min for 52 learners)
- **BART**: mlr3's dbarts implementation dramatically outperforms sl3's (NLL 0.396 vs 0.775)
- **PipeOps**: Preprocessing can be moved inside CV folds (fixes data leakage)
- **Lighter objects**: 1.6 MB vs 215 MB model files
- **Active development**: mlr3 ecosystem is more actively maintained than sl3/tlverse

## Critical Issues to Resolve

### Issue 1: Ranger performance gap (BLOCKER)

**Problem**: Ranger in mlr3 (NLL 0.43-0.48) is far worse than sl3 (0.381).
mlr3superlearner creates a regression task internally for binomial outcomes
(Y is numeric 0/1), so ranger fits a regression forest, not a probability
forest. The predictions are continuous values that may exceed [0,1].

**Fix**: In the forked mlr3superlearner, modify `make_mlr3_task()` to create
a `TaskClassif` (not `TaskRegr`) when `outcome_type == "binomial"`. This
requires Y to be a factor internally, while accepting numeric input.
The task creation should be:
```r
if (outcome_type == "binomial") {
  data[[target]] <- factor(data[[target]], levels = c(0, 1))
  task <- as_task_classif(data, target = target, positive = "1")
} else {
  task <- as_task_regr(data, target = target)
}
```
Then all classif learners will use `predict_type = "prob"` correctly.

### Issue 2: Clustered cross-validation (BLOCKER)

**Problem**: Standard K-fold CV ignores survey cluster structure. Observations
from the same PSU/EA appear in both train and test folds, leaking spatial
information and inflating performance by ~0.02-0.05 AUC.

**Fix**: Modify `make_mlr3_resampling()` to support custom fold assignments:
```r
make_mlr3_resampling <- function(task, folds, outcome_type,
                                  custom_folds = NULL) {
  if (!is.null(custom_folds)) {
    # custom_folds is a list of lists, each with $training_set and $validation_set
    rsmp_custom <- rsmp("custom")
    rsmp_custom$instantiate(task,
      train = lapply(custom_folds, function(f) f$training_set),
      test  = lapply(custom_folds, function(f) f$validation_set))
    return(rsmp_custom)
  }
  # ... existing code for standard CV
}
```

**Integration**: Pre-compute folds using origami:
```r
fold_obj <- origami::make_folds(cluster_ids = cluster_vec, V = K)
custom_folds <- lapply(fold_obj, function(f) {
  list(training_set = f$training_set, validation_set = f$validation_set)
})
```

### Issue 3: Unclipped NLL loss (MODERATE)

**Problem**: `loss_nll` computes `log(0) = -Inf` when predictions are exactly
0 or 1, causing NaN risk for some learners.

**Fix**: Replace `loss_nll` in the fork:
```r
loss_nll <- function(x, y) {
  x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
  -mean(y * log(x) + (1 - y) * log(1 - x))
}
```

## Phase 1: Fork and Fix mlr3superlearner (Week 1)

### 1.1 Fork the repository

```bash
gh repo fork nt-williams/mlr3superlearner --clone
```

### 1.2 Fix task creation (Issue 1)

File: `R/make_mlr3_task.R`

- For `outcome_type == "binomial"`: create `TaskClassif` with factor Y
- For `outcome_type == "continuous"`: keep `TaskRegr` as-is
- Handle the Y conversion (numeric 0/1 -> factor) transparently
- Ensure `predict_type = "prob"` is set for all classif learners

### 1.3 Add clustered CV support (Issue 2)

File: `R/make_mlr3_resampling.R`

- Add `custom_folds` parameter
- Accept origami-style fold objects or plain train/test index lists
- Fall back to standard CV when `custom_folds = NULL`

File: `R/mlr3superlearner.R`

- Add `custom_folds` parameter to main function signature
- Pass through to `make_mlr3_resampling()`
- Document in function help

### 1.4 Clip NLL loss (Issue 3)

File: `R/compute_super_learner_weights.R`

- Add probability clipping to `loss_nll`
- Consider also adding Brier score loss as an option

### 1.5 Validate the fork

- Run on Ghana child_iron (the benchmark case)
- Verify ranger NLL matches sl3 (~0.38)
- Verify clustered folds match origami assignments
- Compare ensemble AUC with sl3 (should be within 0.01)

## Phase 2: Add PipeOps Preprocessing (Week 1-2)

### 2.1 Replace external preprocessing with PipeOps

Current preprocessing (outside CV, causes data leakage):
```
unlabel -> drop all-NA -> NZV -> impute -> NZV -> washb_prescreen(p<0.2) ->
step_zv -> step_nzv -> step_corr(0.85) -> step_normalize
```

New preprocessing (inside CV via PipeOps):
```r
preproc <- po("removeconstants") %>>%       # replaces step_zv + step_nzv
  po("imputemedian") %>>%                    # replaces ck37r::impute_missing_values
  po("removeconstants") %>>%                 # post-imputation NZV
  po("collapseFactors") %>>%                 # handle rare factor levels
  po("scale") %>>%                           # replaces step_normalize
  po("filter", filter = flt("find_correlation"),  # replaces step_corr(0.85)
     filter.cutoff = 0.85)
```

### 2.2 Drop washb_prescreen

**Rationale**:
- Uses full dataset Y for variable selection (data leakage)
- Marginal p-value screening is a weak filter — most variables pass at p<0.2
- Lasso, ranger, and XGBoost all handle variable selection natively
- In the 52-learner comparison, raw ranger (no prescreening) was the best learner
- Prescreened learners (rfimp20_ranger, cor50_ranger) were slightly worse

**If prescreening is desired**: Use `po("filter", filter = flt("importance"),
filter.nfeat = 50)` inside the PipeOp pipeline so it happens within each CV fold.

### 2.3 Create learner-specific pipelines

Some learners benefit from preprocessing, others don't:

```r
# Ranger/XGBoost: only need imputation (handle high-p natively)
ranger_pipe <- preproc_minimal %>>% lrn("classif.ranger", ...)

# GLM/Earth/GAM: need full preprocessing (can't handle high-p)
glm_pipe <- preproc_full %>>% lrn("classif.log_reg")

# BART: needs imputation only
bart_pipe <- preproc_minimal %>>% lrn("classif.bart", ...)
```

This is a major advantage of mlr3 over sl3 — different learners can have
different preprocessing without manually building Pipeline objects.

## Phase 3: Design Production Library (Week 2)

### 3.1 Evidence-based learner selection

Based on the 52-learner sl3 comparison + mlr3 results:

**Fast stack (3-5 learners, ~2-3 min)**:
```r
library_fast <- list(
  "mean",
  list("ranger", num.trees = 500, min.node.size = 10, id = "ranger_main"),
  list("glmnet", alpha = 1, id = "lasso"),
  list("bart", ntree = 100, ndpost = 1000, id = "bart_100")
)
```

**Full stack (8-10 learners, ~5-10 min)**:
```r
library_full <- list(
  # Baseline
  "mean",

  # Regularized linear (complementary signal, 8% NNLS weight in sl3)
  list("glmnet", alpha = 1, id = "lasso"),

  # Random forests (dominant learner, 53% NNLS weight in sl3)
  list("ranger", num.trees = 500, min.node.size = 10, id = "ranger_main"),
  list("ranger", num.trees = 500, min.node.size = 10,
       mtry.ratio = 0.1, id = "ranger_low_mtry"),  # decorrelated

  # BART (best in mlr3, Bayesian uncertainty for free)
  list("bart", ntree = 100, ndpost = 1000, id = "bart_100"),
  list("bart", ntree = 50, ndpost = 500, id = "bart_small"),

  # XGBoost (23% weight when lasso-screened in sl3)
  list("xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
       min_child_weight = 20, subsample = 0.8,
       colsample_bytree = 0.5, id = "xgb_conservative"),
  list("xgboost", max_depth = 6, eta = 0.03, nrounds = 500,
       min_child_weight = 20, subsample = 0.7,
       colsample_bytree = 0.4, id = "xgb_deep"),

  # Gaussian process (competitive in mlr3 at 0.404)
  "gaussianprocess"
)
```

### 3.2 Loss function options

Add support for multiple metalearner loss functions:

```r
# Default: NLL (current, good for calibration)
mlr3superlearner(..., loss = "nll")

# Alternative: Brier score (better for rare outcomes)
mlr3superlearner(..., loss = "brier")

# In the fork, add:
loss_brier <- function(x, y) mean((y - x)^2)
```

For rare outcomes (women_vitA at 2-4%), compare NLL vs Brier ensemble weights.

## Phase 4: Integrate into Targets Pipeline (Week 2-3)

### 4.1 Replace DHS_SL_clustered()

Current flow:
```
DHS_SL_clustered() -> sl3::make_sl3_Task -> sl3::Lrnr_sl$train() -> sl_fit
```

New flow:
```
fit_sl_mlr3() -> origami::make_folds() -> mlr3superlearner() -> mlr3_fit
```

New function in `R/sl_fitting.R`:
```r
fit_sl_mlr3 <- function(d, Xvars, outcome, id, folds, library,
                         outcome_type = "binomial") {
  # 1. Prepare data (strip haven labels, ensure numeric)
  df <- prepare_mlr3_data(d, Xvars, outcome)

  # 2. Create clustered CV folds
  fold_obj <- origami::make_folds(cluster_ids = d[[id]], V = folds)
  custom_folds <- convert_origami_to_mlr3(fold_obj)

  # 3. Fit mlr3superlearner with clustered folds
  fit <- mlr3superlearner(
    data = df, target = "Y", library = library,
    outcome_type = outcome_type,
    custom_folds = custom_folds
  )

  # 4. Extract CV predictions (from internal resampling)
  cv_preds <- extract_cv_predictions(fit)

  # 5. Return in same format as DHS_SL_clustered for compatibility
  list(
    sl_fit = fit,
    res = data.frame(Y = df$Y, yhat_full = cv_preds),
    Xvars = Xvars,
    cv_risk = fit$risk
  )
}
```

### 4.2 Update setup_sl_learners()

Replace sl3 learner construction with mlr3 library lists:
```r
setup_sl_learners <- function(params) {
  if (params$sl_stack == "fast") {
    list(library = library_fast, library_bin = library_fast)
  } else {
    list(library = library_full, library_bin = library_full)
  }
}
```

### 4.3 Update downstream functions

Functions that reference sl3 objects need updating:
- `R/bootstrap.R` -> `R/conformal.R` (already migrated to conformal)
- `R/conceptual_ablation.R` -> update `sl_model$predict()` calls
- `R/transportability.R` -> update `DHS_SL_clustered` calls
- `R/admin1_analysis.R` -> update `sl_result$res` references
- `R/diagnostics.R` -> update prediction extraction

### 4.4 Backward compatibility

Keep `DHS_SL_clustered()` available for comparison but add a switch:
```r
SL_ENGINE <- Sys.getenv("SL_ENGINE", "mlr3")  # or "sl3"
```

## Phase 5: BART-specific Enhancements (Week 3)

### 5.1 BART posterior uncertainty

Since BART uses MCMC, each prediction has a full posterior distribution.
Extract posterior draws for uncertainty quantification:

```r
# dbarts gives ndpost posterior draws per observation
posterior_draws <- predict(bart_fit, newdata, type = "ppd")  # n x ndpost matrix
pred_mean <- rowMeans(posterior_draws)
pred_ci_lo <- apply(posterior_draws, 1, quantile, 0.025)
pred_ci_hi <- apply(posterior_draws, 1, quantile, 0.975)
```

This gives prediction intervals for FREE — no bootstrap or conformal needed
for the BART component. Can be used alongside conformal intervals for the
full ensemble.

### 5.2 BART variable importance

BART provides variable inclusion proportions — how often each variable
appears in the tree sum. This gives a natural importance measure:

```r
var_counts <- colMeans(bart_fit$varcount)  # proportion of splits using each var
```

## Phase 6: Validation (Week 3-4)

### 6.1 Head-to-head comparison

Run both sl3 and mlr3 pipelines on all 4 countries x 4 outcomes:
- Compare CV-AUC, Brier, PR-AUC, calibration
- Verify ranger NLL matches between frameworks (< 0.01 difference)
- Compare ensemble weights
- Compare runtime

### 6.2 Clustered CV validation

Verify that clustered CV gives different (more honest) results:
- Run with standard K-fold and clustered K-fold on same data
- Compare AUC difference (expect ~0.02-0.05 lower with clustering)
- Confirm no cluster leakage

### 6.3 Preprocessing validation

Verify PipeOps preprocessing matches recipes preprocessing:
- Compare number of variables retained after step_corr
- Compare imputed values
- Compare normalized scales

## Timeline

| Week | Task | Deliverable |
|------|------|-------------|
| 1 | Fork mlr3superlearner, fix Issues 1-3 | Working fork with tests |
| 1-2 | PipeOps preprocessing, library design | Validated preprocessing pipeline |
| 2-3 | Integrate into targets pipeline | `fit_sl_mlr3()` + updated targets |
| 3 | BART uncertainty + variable importance | Posterior CIs, var importance |
| 3-4 | Validation against sl3 | Comparison report |

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Ranger still underperforms after fix | Medium | High | Keep sl3 as fallback |
| mlr3superlearner fork diverges from upstream | Low | Medium | Submit PRs upstream |
| PipeOps preprocessing gives different results | Medium | Low | Validate against recipes |
| BART dominates ensemble (low diversity) | Medium | Low | Include diverse learner types |
| Clustered CV breaks some learners | Low | Medium | Use tryCatch per learner |
