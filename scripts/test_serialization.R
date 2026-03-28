# =============================================================================
# Test mlr3superlearner serialization: can we save/load and predict?
#
# Tests three approaches:
#   1. Default RDS (current approach — expected to fail)
#   2. Marshal/unmarshal R6 learners before save
#   3. qs format
#
# Run from project root:
#   source("scripts/test_serialization.R")
# =============================================================================

library(here)
library(mlr3superlearner)
library(mlr3)
library(mlr3learners)

cat("\n=== Testing mlr3superlearner serialization ===\n\n")

# ── Create a small test dataset ──────────────────────────────────────────────
set.seed(42)
n <- 200
p <- 20
X <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(X) <- paste0("x", 1:p)
# Binary outcome with some signal
X$Y <- rbinom(n, 1, plogis(0.5 * X$x1 - 0.3 * X$x2 + 0.2 * X$x3))

cat(sprintf("Test data: n=%d, p=%d, prevalence=%.2f\n\n", n, p, mean(X$Y)))

# ── Fit a small SL ──────────────────────────────────────────────────────────
cat("Fitting mlr3superlearner...\n")
sl_fit <- mlr3superlearner(
  data = X,
  target = "Y",
  library = c("mean", "glmnet", "ranger"),
  outcome_type = "binomial",
  folds = 3
)

# ── Test predictions before serialization ────────────────────────────────────
pred_before <- predict(sl_fit, X[1:10, ])
cat(sprintf("Before save — predictions: %s\n", paste(round(pred_before, 3), collapse = ", ")))
cat(sprintf("Before save — range: [%.3f, %.3f]\n\n", min(pred_before), max(pred_before)))

# ── Test 1: Default saveRDS / readRDS ────────────────────────────────────────
cat("--- Test 1: Default RDS ---\n")
tmp1 <- tempfile(fileext = ".rds")
saveRDS(sl_fit, tmp1)
sl_loaded1 <- readRDS(tmp1)

pred_rds <- tryCatch({
  p <- predict(sl_loaded1, X[1:10, ])
  cat(sprintf("  Predictions: %s\n", paste(round(p, 3), collapse = ", ")))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(p), max(p)))
  cat(sprintf("  Matches original: %s\n", all.equal(p, pred_before)))
  p
}, error = function(e) {
  cat(sprintf("  FAILED: %s\n", e$message))
  NULL
})
cat("\n")

# ── Test 2: Marshal before save, unmarshal after load ────────────────────────
cat("--- Test 2: Marshal/Unmarshal ---\n")

# Check if marshal is available on the learner objects
has_marshal <- tryCatch({
  # mlr3superlearner stores learners in $learners
  learner1 <- sl_fit$learners[[1]]
  "marshal" %in% names(learner1) || is.function(learner1$marshal)
}, error = function(e) FALSE)

cat(sprintf("  Marshal method available: %s\n", has_marshal))

if (has_marshal) {
  # Marshal all learners before saving
  sl_marshal <- sl_fit
  for (i in seq_along(sl_marshal$learners)) {
    tryCatch(sl_marshal$learners[[i]]$marshal(), error = function(e)
      cat(sprintf("  Marshal learner %d failed: %s\n", i, e$message)))
  }

  tmp2 <- tempfile(fileext = ".rds")
  saveRDS(sl_marshal, tmp2)
  sl_loaded2 <- readRDS(tmp2)

  # Unmarshal after loading
  for (i in seq_along(sl_loaded2$learners)) {
    tryCatch(sl_loaded2$learners[[i]]$unmarshal(), error = function(e)
      cat(sprintf("  Unmarshal learner %d failed: %s\n", i, e$message)))
  }

  pred_marshal <- tryCatch({
    p <- predict(sl_loaded2, X[1:10, ])
    cat(sprintf("  Predictions: %s\n", paste(round(p, 3), collapse = ", ")))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(p), max(p)))
    cat(sprintf("  Matches original: %s\n", all.equal(p, pred_before)))
    p
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
    NULL
  })
} else {
  cat("  Skipping — no marshal method found\n")
  # Try alternative: save with qs
}
cat("\n")

# ── Test 3: qs format ───────────────────────────────────────────────────────
cat("--- Test 3: qs format ---\n")
if (requireNamespace("qs", quietly = TRUE)) {
  tmp3 <- tempfile(fileext = ".qs")
  qs::qsave(sl_fit, tmp3)
  sl_loaded3 <- qs::qread(tmp3)

  pred_qs <- tryCatch({
    p <- predict(sl_loaded3, X[1:10, ])
    cat(sprintf("  Predictions: %s\n", paste(round(p, 3), collapse = ", ")))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(p), max(p)))
    cat(sprintf("  Matches original: %s\n", all.equal(p, pred_before)))
    p
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
    NULL
  })
} else {
  cat("  qs package not installed — skipping\n")
}
cat("\n")

# ── Test 4: Permutation on loaded model ──────────────────────────────────────
cat("--- Test 4: Permutation prediction (simulating ablation) ---\n")
# Use whichever loaded model worked
test_model <- NULL
test_label <- ""
if (!is.null(pred_rds) && !all(abs(pred_rds - 0.5) < 0.001)) {
  test_model <- sl_loaded1; test_label <- "RDS"
} else if (exists("pred_marshal") && !is.null(pred_marshal) && !all(abs(pred_marshal - 0.5) < 0.001)) {
  test_model <- sl_loaded2; test_label <- "Marshal"
} else if (exists("pred_qs") && !is.null(pred_qs) && !all(abs(pred_qs - 0.5) < 0.001)) {
  test_model <- sl_loaded3; test_label <- "qs"
}

if (!is.null(test_model)) {
  cat(sprintf("  Using %s-loaded model\n", test_label))
  # Permute x1 (which has signal)
  X_perm <- X
  X_perm$x1 <- X_perm$x1[sample(n)]
  pred_perm <- predict(test_model, X_perm[1:10, ])
  cat(sprintf("  Permuted predictions: %s\n", paste(round(pred_perm, 3), collapse = ", ")))
  cat(sprintf("  Different from original: %s\n", !isTRUE(all.equal(pred_perm, pred_before))))

  # Full AUC comparison
  pred_full <- predict(test_model, X)
  pred_full_perm <- predict(test_model, X_perm)
  auc_orig <- as.numeric(pROC::auc(pROC::roc(X$Y, pred_full, quiet = TRUE)))
  auc_perm <- as.numeric(pROC::auc(pROC::roc(X$Y, pred_full_perm, quiet = TRUE)))
  cat(sprintf("  AUC original: %.3f\n", auc_orig))
  cat(sprintf("  AUC permuted x1: %.3f\n", auc_perm))
  cat(sprintf("  Delta AUC: %.3f (should be negative = x1 is important)\n", auc_perm - auc_orig))
} else {
  cat("  No working loaded model available!\n")
}

cat("\n=== Done ===\n")
