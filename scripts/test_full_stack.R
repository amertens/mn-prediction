# =============================================================================
# Test the FULL learner stack end-to-end:
# Fit → serialize → deserialize → predict → permute
#
# Tests every learner individually AND the full stack combined.
# Checks both binary and continuous outcome types.
# =============================================================================

library(here)
library(mlr3superlearner)
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(pROC)
library(dplyr)
library(caret)
library(ck37r)
library(washb)
library(recipes)
library(data.table)
library(origami)
library(dbarts)

source(here("R", "mlr3_fitting.R"))

cat("\n=== Full Stack Serialization Test ===\n\n")

# ── Test data ────────────────────────────────────────────────────────────────
set.seed(42)
n <- 400; p <- 25
d <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(d) <- paste0("x", 1:p)
d$Y_bin <- rbinom(n, 1, plogis(0.8 * d$x1 - 0.5 * d$x2 + 0.3 * d$x3))
d$Y_cont <- 2 * d$x1 - 1.5 * d$x2 + 0.5 * d$x3 + rnorm(n, 0, 0.5)
d$cluster_id <- rep(1:(n / 10), each = 10)

cat(sprintf("Data: n=%d, p=%d, binary prev=%.2f\n\n", n, p, mean(d$Y_bin)))

# ── Get the full stack spec ──────────────────────────────────────────────────
params <- list(sl_stack = "full", seed = 12345, K = 3, n_perm = 3)
learner_setup <- setup_mlr3_learners(params)
full_lib <- learner_setup$library

cat(sprintf("Full stack: %d learners\n", length(full_lib)))
for (i in seq_along(full_lib)) {
  lrn <- full_lib[[i]]
  if (is.character(lrn)) {
    cat(sprintf("  %d. %s\n", i, lrn))
  } else {
    cat(sprintf("  %d. %s (id=%s)\n", i, lrn[[1]],
                if (!is.null(lrn$id)) lrn$id else "auto"))
  }
}

# ── Test 1: Individual learner serialization ─────────────────────────────────
cat("\n=== Test 1: Individual learners (binary) ===\n")

# Extract just the learner names for individual testing
learner_names <- sapply(full_lib, function(x) {
  if (is.character(x)) x else x[[1]]
})

individual_results <- list()
for (i in seq_along(full_lib)) {
  lrn_spec <- full_lib[[i]]
  lrn_name <- if (is.character(lrn_spec)) lrn_spec else lrn_spec[[1]]
  lrn_id <- if (is.list(lrn_spec) && !is.null(lrn_spec$id)) lrn_spec$id else lrn_name

  cat(sprintf("\n--- %s ---\n", lrn_id))

  fit <- tryCatch(
    mlr3superlearner(data = d[, c(paste0("x", 1:p), "Y_bin")],
                     target = "Y_bin", library = list(lrn_spec),
                     outcome_type = "binomial", folds = 3),
    error = function(e) { cat("  Fit FAILED:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) {
    individual_results[[lrn_id]] <- "FIT_FAIL"
    next
  }

  pred_before <- tryCatch(predict(fit, d[1:5, ]), error = function(e) NULL)
  if (is.null(pred_before)) {
    cat("  Predict before save FAILED\n")
    individual_results[[lrn_id]] <- "PREDICT_FAIL_BEFORE"
    next
  }

  # Serialize round-trip
  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  pred_after <- tryCatch(predict(fit2, d[1:5, ]), error = function(e) NULL)
  if (is.null(pred_after)) {
    cat("  Predict after load FAILED\n")
    individual_results[[lrn_id]] <- "PREDICT_FAIL_AFTER"
    next
  }

  all_same <- length(unique(round(pred_after, 6))) <= 1
  matches <- isTRUE(all.equal(pred_before, pred_after, tolerance = 0.01))

  cat(sprintf("  Before: [%.3f, %.3f], After: [%.3f, %.3f], match=%s, degen=%s\n",
              min(pred_before), max(pred_before),
              min(pred_after), max(pred_after),
              matches, all_same))

  if (all_same) {
    individual_results[[lrn_id]] <- "DEGENERATE"
  } else if (!matches) {
    individual_results[[lrn_id]] <- "CHANGED"
  } else {
    individual_results[[lrn_id]] <- "OK"
  }

  file.remove(tmp)
}

cat("\n=== Individual learner summary ===\n")
for (nm in names(individual_results)) {
  status <- individual_results[[nm]]
  symbol <- switch(status, OK = "OK", DEGENERATE = "DEGEN", CHANGED = "CHANGED",
                   FIT_FAIL = "FIT_FAIL", "??")
  cat(sprintf("  [%s] %s\n", symbol, nm))
}

# ── Test 2: Full stack via our wrapper (binary) ─────────────────────────────
cat("\n=== Test 2: Full stack via mlr3_SL_clustered (binary) ===\n")

result_bin <- tryCatch(
  mlr3_SL_clustered(
    d = d, Xvars = paste0("x", 1:p), outcome = "Y_bin",
    population = "all", id = "cluster_id",
    folds = 3, mlr3_library = full_lib,
    outcome_type = "binomial", prescreen = FALSE
  ),
  error = function(e) { cat("Full stack fit FAILED:", e$message, "\n"); NULL }
)

if (!is.null(result_bin)) {
  Y <- result_bin$res$Y
  yhat <- result_bin$res$yhat_full
  auc_cv <- as.numeric(auc(roc(Y, yhat, quiet = TRUE)))
  cat(sprintf("CV AUC: %.3f\n", auc_cv))

  # Weights
  wts <- result_bin$sl_fit$sl_weights
  lms <- result_bin$sl_fit$learner_models
  if (!is.null(wts) && !is.null(lms)) {
    n_lms <- min(length(wts), length(lms))
    for (idx in seq_len(n_lms)) {
      if (wts[idx] > 0.001) {
        cat(sprintf("  Weight %.3f -> %s (is_bart=%s)\n",
                    wts[idx], lms[[idx]]$id, isTRUE(lms[[idx]]$is_bart)))
      }
    }
  }

  # Serialize round-trip
  cat("\nSerialize → predict test...\n")
  tmp <- tempfile(fileext = ".rds")
  saveRDS(result_bin, tmp)
  loaded <- readRDS(tmp)
  wrapper <- loaded$sl_fit

  pred_after <- tryCatch(
    as.numeric(wrapper$predict(d[1:10, ])),
    error = function(e) { cat("  Predict FAILED:", e$message, "\n"); NULL }
  )

  if (!is.null(pred_after)) {
    degen <- length(unique(round(pred_after, 6))) <= 1
    cat(sprintf("  After load: range [%.3f, %.3f], unique=%d, degen=%s\n",
                min(pred_after), max(pred_after),
                length(unique(round(pred_after, 6))), degen))

    if (!degen) {
      # Permutation test
      d_perm <- d; d_perm$x1 <- d_perm$x1[sample(n)]
      pred_perm <- as.numeric(wrapper$predict(d_perm[1:10, ]))
      changed <- !isTRUE(all.equal(pred_after, pred_perm))
      cat(sprintf("  Permutation changes predictions: %s\n", changed))

      # Full AUC comparison
      pred_full <- as.numeric(wrapper$predict(d))
      pred_full_perm <- as.numeric(wrapper$predict(d_perm))
      auc_full <- as.numeric(auc(roc(d$Y_bin, pred_full, quiet = TRUE)))
      auc_perm <- as.numeric(auc(roc(d$Y_bin, pred_full_perm, quiet = TRUE)))
      cat(sprintf("  Full AUC: %.3f, Permuted: %.3f, Delta: %.4f\n",
                  auc_full, auc_perm, auc_perm - auc_full))
      cat("  BINARY FULL STACK: PASS\n")
    } else {
      cat("  BINARY FULL STACK: FAIL (degenerate predictions)\n")
    }
  }
  file.remove(tmp)
}

# ── Test 3: Full stack continuous ────────────────────────────────────────────
cat("\n=== Test 3: Full stack via mlr3_SL_clustered (continuous) ===\n")

result_cont <- tryCatch(
  mlr3_SL_clustered(
    d = d, Xvars = paste0("x", 1:p), outcome = "Y_cont",
    population = "all", id = "cluster_id",
    folds = 3, mlr3_library = full_lib,
    outcome_type = "continuous", prescreen = FALSE
  ),
  error = function(e) { cat("Full stack cont fit FAILED:", e$message, "\n"); NULL }
)

if (!is.null(result_cont)) {
  Y <- result_cont$res$Y
  yhat <- result_cont$res$yhat_full
  r2 <- 1 - sum((Y - yhat)^2) / sum((Y - mean(Y))^2)
  cat(sprintf("CV R2: %.3f\n", r2))

  # Weights
  wts <- result_cont$sl_fit$sl_weights
  lms <- result_cont$sl_fit$learner_models
  if (!is.null(wts) && !is.null(lms)) {
    n_lms <- min(length(wts), length(lms))
    for (idx in seq_len(n_lms)) {
      if (wts[idx] > 0.001) {
        cat(sprintf("  Weight %.3f -> %s (is_bart=%s)\n",
                    wts[idx], lms[[idx]]$id, isTRUE(lms[[idx]]$is_bart)))
      }
    }
  }

  # Serialize round-trip
  cat("\nSerialize → predict test...\n")
  tmp <- tempfile(fileext = ".rds")
  saveRDS(result_cont, tmp)
  loaded <- readRDS(tmp)
  wrapper <- loaded$sl_fit

  pred_after <- tryCatch(
    as.numeric(wrapper$predict(d[1:10, ])),
    error = function(e) { cat("  Predict FAILED:", e$message, "\n"); NULL }
  )

  if (!is.null(pred_after)) {
    degen <- length(unique(round(pred_after, 6))) <= 1
    cat(sprintf("  After load: range [%.3f, %.3f], unique=%d, degen=%s\n",
                min(pred_after), max(pred_after),
                length(unique(round(pred_after, 6))), degen))

    if (!degen) {
      d_perm <- d; d_perm$x1 <- d_perm$x1[sample(n)]
      pred_perm <- as.numeric(wrapper$predict(d_perm[1:10, ]))
      changed <- !isTRUE(all.equal(pred_after, pred_perm))
      cat(sprintf("  Permutation changes predictions: %s\n", changed))
      cat("  CONTINUOUS FULL STACK: PASS\n")
    } else {
      cat("  CONTINUOUS FULL STACK: FAIL (degenerate predictions)\n")
    }
  }
  file.remove(tmp)
}

cat("\n=== All tests complete ===\n")
