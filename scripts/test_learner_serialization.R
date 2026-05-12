# Test which mlr3 learners survive serialization
library(mlr3superlearner)
library(mlr3)
library(mlr3learners)

set.seed(42)
n <- 300; p <- 15
X <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(X) <- paste0("x", 1:p)
X$Y <- rbinom(n, 1, plogis(0.5 * X$x1 - 0.3 * X$x2))

cat("=== Testing individual learner serialization ===\n\n")

learners <- c("mean", "glmnet", "ranger", "xgboost")

for (lrn_name in learners) {
  cat(sprintf("--- %s ---\n", lrn_name))

  fit <- tryCatch(
    mlr3superlearner(data = X, target = "Y", library = lrn_name,
                     outcome_type = "binomial", folds = 3),
    error = function(e) { cat("  Fit failed:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) next

  # Predict before save
  pred_before <- predict(fit, X[1:5, ])
  cat(sprintf("  Before: %s\n", paste(round(pred_before, 4), collapse = ", ")))

  # Save and reload
  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  pred_after <- tryCatch(
    predict(fit2, X[1:5, ]),
    error = function(e) { cat("  Predict FAILED:", e$message, "\n"); NULL }
  )

  if (!is.null(pred_after)) {
    cat(sprintf("  After:  %s\n", paste(round(pred_after, 4), collapse = ", ")))
    cat(sprintf("  Match:  %s\n", all.equal(pred_before, pred_after)))
    # Also test permuted
    X_perm <- X; X_perm$x1 <- X_perm$x1[sample(n)]
    pred_perm <- predict(fit2, X_perm[1:5, ])
    cat(sprintf("  Permuted: %s\n", paste(round(pred_perm, 4), collapse = ", ")))
    cat(sprintf("  Changed: %s\n", !isTRUE(all.equal(pred_after, pred_perm))))
  }
  cat("\n")
}

# Now test the FULL stack (all learners together)
cat("--- Full stack (mean + glmnet + ranger + xgboost) ---\n")
fit_full <- mlr3superlearner(
  data = X, target = "Y",
  library = c("mean", "glmnet", "ranger", "xgboost"),
  outcome_type = "binomial", folds = 3
)
pred_before <- predict(fit_full, X[1:5, ])
cat(sprintf("  Before: %s\n", paste(round(pred_before, 4), collapse = ", ")))
cat(sprintf("  Weights: %s\n", paste(round(fit_full$weights, 3), collapse = ", ")))

tmp <- tempfile(fileext = ".rds")
saveRDS(fit_full, tmp)
fit_full2 <- readRDS(tmp)

pred_after <- tryCatch(
  predict(fit_full2, X[1:5, ]),
  error = function(e) { cat("  FAILED:", e$message, "\n"); NULL }
)
if (!is.null(pred_after)) {
  cat(sprintf("  After:  %s\n", paste(round(pred_after, 4), collapse = ", ")))
  cat(sprintf("  Match:  %s\n", all.equal(pred_before, pred_after)))
}

# Test with BART too
cat("\n--- With BART ---\n")
fit_bart <- tryCatch(
  mlr3superlearner(data = X, target = "Y",
                   library = c("mean", "glmnet", "ranger", "bart"),
                   outcome_type = "binomial", folds = 3),
  error = function(e) { cat("  Fit failed:", e$message, "\n"); NULL }
)
if (!is.null(fit_bart)) {
  pred_before <- predict(fit_bart, X[1:5, ])
  cat(sprintf("  Before: %s\n", paste(round(pred_before, 4), collapse = ", ")))
  cat(sprintf("  Weights: %s\n", paste(round(fit_bart$weights, 3), collapse = ", ")))

  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit_bart, tmp)
  fit_bart2 <- readRDS(tmp)

  pred_after <- tryCatch(
    predict(fit_bart2, X[1:5, ]),
    error = function(e) { cat("  FAILED:", e$message, "\n"); NULL }
  )
  if (!is.null(pred_after)) {
    cat(sprintf("  After:  %s\n", paste(round(pred_after, 4), collapse = ", ")))
    cat(sprintf("  Match:  %s\n", all.equal(pred_before, pred_after)))
  }
}

cat("\n=== Done ===\n")
