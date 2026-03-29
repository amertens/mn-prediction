# Test lightgbm in mlr3superlearner: fit, serialize, predict, permute
library(mlr3superlearner)
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(pROC)

set.seed(42)
n <- 400; p <- 20
d <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(d) <- paste0("x", 1:p)
d$Y <- rbinom(n, 1, plogis(0.8 * d$x1 - 0.5 * d$x2 + 0.3 * d$x3))

cat("=== LightGBM serialization test ===\n\n")

# Test 1: Individual fit
cat("--- Individual lightgbm ---\n")
fit <- tryCatch(
  mlr3superlearner(data = d, target = "Y",
                   library = list(list("lightgbm", num_leaves = 31,
                                       learning_rate = 0.05,
                                       num_iterations = 100,
                                       min_data_in_leaf = 20,
                                       id = "lgbm_test")),
                   outcome_type = "binomial", folds = 3),
  error = function(e) { cat("Fit FAILED:", e$message, "\n"); NULL }
)

if (!is.null(fit)) {
  pred_before <- predict(fit, d[1:10, ])
  cat(sprintf("Before: range [%.3f, %.3f]\n", min(pred_before), max(pred_before)))

  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  pred_after <- tryCatch(predict(fit2, d[1:10, ]), error = function(e) NULL)
  if (!is.null(pred_after)) {
    degen <- length(unique(round(pred_after, 6))) <= 1
    cat(sprintf("After RDS: range [%.3f, %.3f], match=%s, degen=%s\n",
                min(pred_after), max(pred_after),
                isTRUE(all.equal(pred_before, pred_after)), degen))
  } else {
    cat("After RDS: predict FAILED\n")
  }

  # Socket serialization
  fit3 <- unserialize(serialize(fit, NULL))
  pred_socket <- tryCatch(predict(fit3, d[1:10, ]), error = function(e) NULL)
  if (!is.null(pred_socket)) {
    degen <- length(unique(round(pred_socket, 6))) <= 1
    cat(sprintf("After socket: range [%.3f, %.3f], degen=%s\n",
                min(pred_socket), max(pred_socket), degen))
  } else {
    cat("After socket: predict FAILED\n")
  }

  # Permutation
  d_perm <- d; d_perm$x1 <- d_perm$x1[sample(n)]
  pred_perm <- tryCatch(predict(fit2, d_perm[1:10, ]), error = function(e) NULL)
  if (!is.null(pred_perm) && !is.null(pred_after)) {
    cat(sprintf("Permutation changes: %s\n", !isTRUE(all.equal(pred_after, pred_perm))))
  }
}

# Test 2: In a mixed stack
cat("\n--- Mixed stack with lightgbm ---\n")
fit_mix <- tryCatch(
  mlr3superlearner(data = d, target = "Y",
                   library = list("mean",
                                  list("glmnet", alpha = 1, id = "lasso"),
                                  list("lightgbm", num_leaves = 31,
                                       learning_rate = 0.05,
                                       num_iterations = 100,
                                       id = "lgbm")),
                   outcome_type = "binomial", folds = 3),
  error = function(e) { cat("Mixed fit FAILED:", e$message, "\n"); NULL }
)

if (!is.null(fit_mix)) {
  cat(sprintf("Weights: %s\n", paste(round(fit_mix$weights, 3), collapse = ", ")))
  auc_val <- as.numeric(auc(roc(d$Y, predict(fit_mix, d), quiet = TRUE)))
  cat(sprintf("AUC: %.3f\n", auc_val))

  # Continuous test
  d$Y_cont <- 2 * d$x1 - d$x2 + rnorm(n, 0, 0.5)
  fit_cont <- tryCatch(
    mlr3superlearner(data = d[, c(paste0("x", 1:p), "Y_cont")],
                     target = "Y_cont",
                     library = list(list("lightgbm", num_leaves = 31,
                                         learning_rate = 0.05,
                                         num_iterations = 100,
                                         id = "lgbm_reg")),
                     outcome_type = "continuous", folds = 3),
    error = function(e) { cat("Continuous fit FAILED:", e$message, "\n"); NULL }
  )
  if (!is.null(fit_cont)) {
    p_cont <- predict(fit_cont, d[1:5, ])
    cat(sprintf("Continuous predictions: %s\n", paste(round(p_cont, 3), collapse = ", ")))

    fit_cont2 <- readRDS({tmp2 <- tempfile(); saveRDS(fit_cont, tmp2); tmp2})
    p_cont2 <- tryCatch(predict(fit_cont2, d[1:5, ]), error = function(e) NULL)
    if (!is.null(p_cont2)) {
      cat(sprintf("Continuous after RDS: %s, match=%s\n",
                  paste(round(p_cont2, 3), collapse = ", "),
                  isTRUE(all.equal(p_cont, p_cont2))))
    }
  }
}

cat("\n=== Done ===\n")
