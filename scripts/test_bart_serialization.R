# =============================================================================
# Debug BART serialization issue in mlr3superlearner
#
# Goal: understand exactly why BART predictions fail after serialization
# and find a workaround that lets us use BART in the full stack.
#
# Hypothesis: dbarts stores C++ external pointers that become invalid
# after RDS save/load. We need to find what breaks and whether
# re-creating the model object from stored components is feasible.
# =============================================================================

library(mlr3superlearner)
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(pROC)

cat("\n=== BART Serialization Deep-Dive ===\n\n")

# ── Create test data ─────────────────────────────────────────────────────────
set.seed(42)
n <- 500; p <- 30
X <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(X) <- paste0("x", 1:p)
X$Y <- rbinom(n, 1, plogis(0.8 * X$x1 - 0.5 * X$x2 + 0.3 * X$x3 * X$x4))
cat(sprintf("Data: n=%d, p=%d, prevalence=%.2f\n\n", n, p, mean(X$Y)))

# ── Test 1: BART alone in mlr3superlearner ───────────────────────────────────
cat("=== Test 1: BART-only SL ===\n")
fit_bart <- mlr3superlearner(
  data = X, target = "Y",
  library = list(list("bart", ntree = 100, id = "bart_100")),
  outcome_type = "binomial", folds = 3
)

pred_before <- predict(fit_bart, X[1:10, ])
cat(sprintf("Before save: range [%.3f, %.3f]\n", min(pred_before), max(pred_before)))

# Save/load via RDS
tmp_rds <- tempfile(fileext = ".rds")
saveRDS(fit_bart, tmp_rds)
fit_bart_loaded <- readRDS(tmp_rds)

pred_after_rds <- tryCatch(
  predict(fit_bart_loaded, X[1:10, ]),
  error = function(e) { cat("  RDS predict ERROR:", e$message, "\n"); NULL }
)
if (!is.null(pred_after_rds)) {
  cat(sprintf("After RDS:  range [%.3f, %.3f]\n", min(pred_after_rds), max(pred_after_rds)))
  cat(sprintf("Match: %s\n", all.equal(pred_before, pred_after_rds)))
} else {
  cat("After RDS: predict returned NULL\n")
}

# ── Test 2: Examine the BART learner model object ────────────────────────────
cat("\n=== Test 2: Inspect BART internals ===\n")
bart_learner <- fit_bart$learners[[1]]
cat(sprintf("Learner class: %s\n", paste(class(bart_learner), collapse = ", ")))
cat(sprintf("Learner state: %s\n", bart_learner$state$train_task$nrow))

# Check the raw model
bart_model <- bart_learner$model
cat(sprintf("Model class: %s\n", paste(class(bart_model), collapse = ", ")))
cat(sprintf("Model names: %s\n", paste(names(bart_model), collapse = ", ")))

# Can we predict from the raw model directly?
cat("\n=== Test 3: Predict from raw bart model ===\n")
raw_pred <- tryCatch({
  # dbarts predict method
  p <- predict(bart_model, as.matrix(X[1:10, paste0("x", 1:p)]))
  cat(sprintf("Raw predict class: %s, dims: %s\n", class(p), paste(dim(p), collapse = "x")))
  if (is.matrix(p)) {
    # BART returns posterior draws matrix, take column means
    cat(sprintf("Posterior draws: %d\n", ncol(p)))
    p <- pnorm(rowMeans(p))  # probit link for binary
  }
  cat(sprintf("Raw predict range: [%.3f, %.3f]\n", min(p), max(p)))
  p
}, error = function(e) {
  cat(sprintf("Raw predict FAILED: %s\n", e$message))
  NULL
})

# ── Test 4: Save/load just the raw model ─────────────────────────────────────
cat("\n=== Test 4: Save/load raw BART model ===\n")
tmp2 <- tempfile(fileext = ".rds")
saveRDS(bart_model, tmp2)
bart_model_loaded <- readRDS(tmp2)

raw_pred_loaded <- tryCatch({
  p <- predict(bart_model_loaded, as.matrix(X[1:10, paste0("x", 1:p)]))
  if (is.matrix(p)) p <- pnorm(rowMeans(p))
  cat(sprintf("Loaded raw predict range: [%.3f, %.3f]\n", min(p), max(p)))
  cat(sprintf("Match original: %s\n", all.equal(raw_pred, p)))
  p
}, error = function(e) {
  cat(sprintf("Loaded raw predict FAILED: %s\n", e$message))
  NULL
})

# ── Test 5: Check if dbarts has a serialize/deserialize mechanism ────────────
cat("\n=== Test 5: dbarts serialization support ===\n")
if (requireNamespace("dbarts", quietly = TRUE)) {
  # Check for makeModelMatrixFromDataFrame or similar
  cat(sprintf("dbarts version: %s\n", packageVersion("dbarts")))
  # Check if the fit object has a save/load mechanism
  bart_fit_class <- class(bart_model)
  cat(sprintf("Fit class: %s\n", paste(bart_fit_class, collapse = ", ")))

  # Check if it's a dbarts::bart or dbarts::dbartsControl object
  if (inherits(bart_model, "bart")) {
    # Standard bart object — check for posterior samples
    cat(sprintf("Has yhat.train: %s\n", "yhat.train" %in% names(bart_model)))
    cat(sprintf("Has fit: %s\n", "fit" %in% names(bart_model)))
    cat(sprintf("Has call: %s\n", "call" %in% names(bart_model)))
  }

  # Try dbarts::bart2 which may serialize better
  cat("\n--- Testing dbarts::bart2 directly ---\n")
  bart2_fit <- tryCatch({
    dbarts::bart2(
      as.matrix(X[, paste0("x", 1:p)]),
      X$Y,
      ntree = 50,
      keepTrees = TRUE,
      verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("bart2 failed: %s\n", e$message))
    NULL
  })

  if (!is.null(bart2_fit)) {
    pred_b2_before <- predict(bart2_fit, as.matrix(X[1:10, paste0("x", 1:p)]))
    if (is.matrix(pred_b2_before)) pred_b2_before <- colMeans(pnorm(pred_b2_before))
    cat(sprintf("bart2 before save: range [%.3f, %.3f]\n", min(pred_b2_before), max(pred_b2_before)))

    tmp3 <- tempfile(fileext = ".rds")
    saveRDS(bart2_fit, tmp3)
    bart2_loaded <- readRDS(tmp3)

    pred_b2_after <- tryCatch({
      p <- predict(bart2_loaded, as.matrix(X[1:10, paste0("x", 1:p)]))
      if (is.matrix(p)) p <- colMeans(pnorm(p))
      cat(sprintf("bart2 after load: range [%.3f, %.3f]\n", min(p), max(p)))
      cat(sprintf("Match: %s\n", all.equal(pred_b2_before, p)))
      p
    }, error = function(e) {
      cat(sprintf("bart2 loaded predict FAILED: %s\n", e$message))
      NULL
    })
  }
}

# ── Test 6: Socket serialization (simulates targets/future) ──────────────────
cat("\n=== Test 6: Socket serialization (simulates targets) ===\n")
# This is what targets actually does — serialize via connection, not file
raw_bytes <- serialize(fit_bart, NULL)  # serialize to raw vector
fit_bart_socket <- unserialize(raw_bytes)

pred_socket <- tryCatch(
  predict(fit_bart_socket, X[1:10, ]),
  error = function(e) { cat("  Socket predict ERROR:", e$message, "\n"); NULL }
)
if (!is.null(pred_socket)) {
  cat(sprintf("After socket: range [%.3f, %.3f]\n", min(pred_socket), max(pred_socket)))
  cat(sprintf("Match original: %s\n", all.equal(pred_before, pred_socket)))
  cat(sprintf("All same value: %s\n", length(unique(round(pred_socket, 6))) <= 1))
}

# ── Test 7: mlr3 learner marshal/unmarshal ───────────────────────────────────
cat("\n=== Test 7: mlr3 learner marshal/unmarshal ===\n")
# Check if the BART learner supports marshaling
cat(sprintf("Has marshal: %s\n", "marshal" %in% ls(bart_learner)))
cat(sprintf("Has marshaled_model: %s\n", !is.null(bart_learner$marshaled_model)))

# Try marshaling
tryCatch({
  bart_learner$marshal()
  cat("Marshal succeeded\n")

  # Save marshaled version
  tmp_marshal <- tempfile(fileext = ".rds")
  saveRDS(bart_learner, tmp_marshal)
  bart_learner_loaded <- readRDS(tmp_marshal)

  bart_learner_loaded$unmarshal()
  cat("Unmarshal succeeded\n")

  # Test predict on unmarshaled learner
  task <- mlr3::TaskClassif$new("test", X, target = "Y")
  pred_marshal <- bart_learner_loaded$predict_newdata(X[1:10, paste0("x", 1:p)])
  cat(sprintf("Marshaled predict: %s\n", paste(round(pred_marshal$prob[, "1"], 3), collapse = ", ")))
}, error = function(e) {
  cat(sprintf("Marshal/unmarshal failed: %s\n", e$message))
})

cat("\n=== Done ===\n")
