# =============================================================================
# Test 2: Reproduce the exact ablation bug using actual pipeline objects
#
# The simple test passed, so the issue must be in how the wrapper's
# predict() interacts with the actual pipeline data.
# =============================================================================

library(here)
library(targets)
library(mlr3superlearner)
library(mlr3)
library(pROC)

cat("\n=== Testing ablation predict on actual pipeline data ===\n\n")

# Load a small SL fit
cat("Loading sl_fit_sierraleone_child_iron...\n")
sl <- tar_read(sl_fit_sierraleone_child_iron, store = "_targets")
bf <- sl$bin_fit

cat(sprintf("Names: %s\n", paste(names(bf), collapse = ", ")))
cat(sprintf("Xvars: %d\n", length(bf$Xvars)))
cat(sprintf("train_data: %d x %d\n", nrow(bf$train_data), ncol(bf$train_data)))

Y <- as.numeric(bf$res$Y)
yhat_cv <- as.numeric(bf$res$yhat_full)
cat(sprintf("Y range: [%.2f, %.2f], prevalence: %.3f\n", min(Y), max(Y), mean(Y)))
cat(sprintf("yhat_cv range: [%.3f, %.3f]\n", min(yhat_cv), max(yhat_cv)))

# ── Test A: Wrapper predict on original train_data ────────────────────────
cat("\n--- Test A: Wrapper predict on train_data ---\n")
td <- data.frame(bf$train_data)
wrapper <- bf$sl_fit

pred_wrapper <- tryCatch({
  p <- as.numeric(wrapper$predict(td))
  cat(sprintf("  Range: [%.4f, %.4f]\n", min(p), max(p)))
  cat(sprintf("  Mean: %.4f\n", mean(p)))
  cat(sprintf("  All 0.5? %s\n", all(abs(p - 0.5) < 0.001)))
  cat(sprintf("  Cor with yhat_cv: %.3f\n", cor(p, yhat_cv)))
  p
}, error = function(e) {
  cat(sprintf("  FAILED: %s\n", e$message))
  NULL
})

# ── Test B: Direct mlr3superlearner predict (bypass wrapper) ──────────────
cat("\n--- Test B: Direct mlr3superlearner predict ---\n")
mlr3_obj <- wrapper$mlr3_fit
cat(sprintf("  mlr3_fit class: %s\n", paste(class(mlr3_obj), collapse = ", ")))

pred_direct <- tryCatch({
  p <- as.numeric(predict(mlr3_obj, td))
  cat(sprintf("  Range: [%.4f, %.4f]\n", min(p), max(p)))
  cat(sprintf("  Mean: %.4f\n", mean(p)))
  cat(sprintf("  All 0.5? %s\n", all(abs(p - 0.5) < 0.001)))
  cat(sprintf("  Cor with yhat_cv: %.3f\n", cor(p, yhat_cv)))
  p
}, error = function(e) {
  cat(sprintf("  FAILED: %s\n", e$message))
  NULL
})

# ── Test C: Check what predict.mlr3superlearner sees ──────────────────────
cat("\n--- Test C: Inspect mlr3 internals ---\n")
cat(sprintf("  object$x: %s\n", paste(head(mlr3_obj$x, 5), collapse = ", ")))
cat(sprintf("  n features in object$x: %d\n", length(mlr3_obj$x)))
cat(sprintf("  n learners: %d\n", length(mlr3_obj$learners)))
cat(sprintf("  weights: %s\n", paste(round(head(mlr3_obj$weights, 5), 3), collapse = ", ")))
cat(sprintf("  outcome_type: %s\n", mlr3_obj$outcome_type))

# Test individual learner predict
cat("\n--- Test D: Individual learner predictions ---\n")
for (i in seq_along(mlr3_obj$learners)) {
  lrn <- mlr3_obj$learners[[i]]
  cat(sprintf("  Learner %d: %s (class: %s)\n", i, lrn$id, paste(class(lrn), collapse=",")))
  pred_i <- tryCatch({
    p <- lrn$predict_newdata(td[1:5, mlr3_obj$x, drop = FALSE])
    if (mlr3_obj$outcome_type == "binomial") {
      probs <- p$prob[, "1"]
    } else {
      probs <- p$response
    }
    cat(sprintf("    Predictions: %s\n", paste(round(probs, 3), collapse = ", ")))
    probs
  }, error = function(e) {
    cat(sprintf("    FAILED: %s\n", e$message))
    NULL
  })
}

# ── Test E: Permutation test ─────────────────────────────────────────────
cat("\n--- Test E: Permutation on first predictor ---\n")
if (!is.null(pred_wrapper) && !all(abs(pred_wrapper - 0.5) < 0.001)) {
  use_pred <- pred_wrapper
  use_label <- "wrapper"
} else if (!is.null(pred_direct) && !all(abs(pred_direct - 0.5) < 0.001)) {
  use_pred <- pred_direct
  use_label <- "direct"
} else {
  cat("  Both predict methods return 0.5 — BUG CONFIRMED\n")
  use_pred <- NULL
}

if (!is.null(use_pred)) {
  td_perm <- td
  first_xvar <- bf$Xvars[1]
  td_perm[[first_xvar]] <- td_perm[[first_xvar]][sample(nrow(td_perm))]

  pred_perm <- tryCatch({
    if (use_label == "wrapper") as.numeric(wrapper$predict(td_perm))
    else as.numeric(predict(mlr3_obj, td_perm))
  }, error = function(e) NULL)

  if (!is.null(pred_perm)) {
    auc_orig <- as.numeric(pROC::auc(pROC::roc(Y, use_pred, quiet = TRUE)))
    auc_perm <- as.numeric(pROC::auc(pROC::roc(Y, pred_perm, quiet = TRUE)))
    cat(sprintf("  AUC original: %.3f\n", auc_orig))
    cat(sprintf("  AUC permuted %s: %.3f\n", first_xvar, auc_perm))
    cat(sprintf("  Delta: %.4f\n", auc_perm - auc_orig))
  }
}

cat("\n=== Done ===\n")
