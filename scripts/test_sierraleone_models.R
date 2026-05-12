# =============================================================================
# scripts/test_sierraleone_models.R
#
# Test alternative model configs for Sierra Leone AUC=0.5 outcomes.
# Uses the ALREADY PREPROCESSED data from sl_fit objects (avoids reprocessing).
# =============================================================================

library(targets)
library(here)
library(dplyr)
library(mlr3superlearner)
library(mlr3extralearners)
library(pROC)
library(origami)

TARGETS_STORE <- here("_targets_full")

# ── Monkey-patch NLL ──────────────────────────────────────────────────────
env <- asNamespace("mlr3superlearner")
if (exists("loss_nll", envir = env)) {
  clipped_nll <- function(x, y) {
    x <- pmin(pmax(x, 1e-15), 1 - 1e-15)
    -mean(y * log(x) + (1 - y) * log(1 - x))
  }
  tryCatch(assignInNamespace("loss_nll", clipped_nll, "mlr3superlearner"),
           error = function(e) NULL)
}

# ── Load preprocessed data from existing sl_fit objects ───────────────────
test_outcomes <- c("child_vitA", "child_iron", "women_vitA")

cat("Loading preprocessed data from sl_fit objects...\n")
prepped <- list()
for (oc in test_outcomes) {
  tname <- paste0("sl_fit_sierraleone_", oc)
  sl <- tryCatch(tar_read_raw(tname, store = TARGETS_STORE), error = function(e) NULL)
  if (is.null(sl)) { cat(sprintf("  %s: target not found\n", oc)); next }

  bf <- sl[["bin_fit"]]
  if (is.null(bf)) { cat(sprintf("  %s: no bin_fit\n", oc)); next }

  Y <- as.numeric(bf[["res"]][["Y"]])
  train_data <- bf[["train_data"]]
  covars <- bf[["Xvars"]]

  if (is.null(train_data) || length(Y) == 0) {
    cat(sprintf("  %s: no train_data or Y\n", oc))
    next
  }

  X <- as.data.frame(train_data)[, covars, drop = FALSE]

  # Clean any remaining issues
  for (col in colnames(X)) {
    bad <- !is.finite(X[[col]])
    if (any(bad)) X[[col]][bad] <- 0
  }

  n_events <- sum(Y == 1)
  prev <- mean(Y == 1)
  cat(sprintf("  %s: n=%d, p=%d, events=%d, prev=%.1f%%\n",
              oc, length(Y), ncol(X), n_events, prev * 100))

  # Get cluster IDs for fold creation
  # Try to extract from the original data
  cluster_id <- tryCatch({
    od <- tar_read_raw(paste0("outcome_data_sierraleone_", oc), store = TARGETS_STORE)
    od[["data"]][["gw_cnum"]]
  }, error = function(e) seq_len(length(Y)))

  if (length(cluster_id) != length(Y)) {
    cat(sprintf("    cluster_id length mismatch (%d vs %d), using row indices\n",
                length(cluster_id), length(Y)))
    cluster_id <- seq_len(length(Y))
  }

  prepped[[oc]] <- list(X = X, Y = Y, cluster_id = cluster_id,
                         n_events = n_events, prev = prev)
}


# ── Learner libraries ────────────────────────────────────────────────────

# Lib 1: No mean, standard learners
lib_no_mean <- list(
  list("glmnet", alpha = 1, id = "lasso"),
  list("glmnet", alpha = 0, id = "ridge"),
  list("glmnet", alpha = 0.5, id = "elastic_net"),
  list("ranger", num.trees = 500, min.node.size = 10, id = "ranger"),
  list("ranger", num.trees = 1000, min.node.size = 5, id = "ranger_deep"),
  list("xgboost", max_depth = 3, eta = 0.05, nrounds = 300,
       min_child_weight = 20, subsample = 0.8, colsample_bytree = 0.5,
       id = "xgb"),
  list("bart", ntree = 100, id = "bart")
)

# Lib 2: Class-weighted for rare outcomes
make_lib_weighted <- function(prev_rate) {
  pw <- min((1 - prev_rate) / prev_rate, 50)
  list(
    list("xgboost", max_depth = 2, eta = 0.05, nrounds = 500,
         min_child_weight = 5, subsample = 0.8, colsample_bytree = 0.6,
         scale_pos_weight = pw, id = "xgb_wt_shallow"),
    list("xgboost", max_depth = 4, eta = 0.03, nrounds = 500,
         min_child_weight = 10, subsample = 0.8, colsample_bytree = 0.5,
         scale_pos_weight = pw, id = "xgb_wt_deep"),
    list("ranger", num.trees = 1000, min.node.size = 3,
         class.weights = c("0" = 1, "1" = pw), id = "ranger_wt"),
    list("bart", ntree = 200, id = "bart_large"),
    list("bart", ntree = 50, id = "bart_small"),
    list("glmnet", alpha = 1, id = "lasso"),
    list("glmnet", alpha = 0.5, id = "enet")
  )
}

# Lib 3: BART-only (often best for rare outcomes)
lib_bart_only <- list(
  list("bart", ntree = 50, id = "bart_50"),
  list("bart", ntree = 100, id = "bart_100"),
  list("bart", ntree = 200, id = "bart_200"),
  list("glmnet", alpha = 0.5, id = "enet")
)

# Lib 4: Simple (minimal overfitting)
lib_simple <- list(
  list("glmnet", alpha = 1, id = "lasso"),
  list("glmnet", alpha = 0.5, id = "enet"),
  list("ranger", num.trees = 500, min.node.size = 10, id = "ranger"),
  list("bart", ntree = 100, id = "bart")
)


# ── Test function ─────────────────────────────────────────────────────────
test_library <- function(X, Y, cluster_id, lib, lib_name, n_folds = 5) {
  cat(sprintf("\n  --- %s (folds=%d) ---\n", lib_name, n_folds))

  if (sum(Y == 1) < 3) {
    cat("    Too few events\n")
    return(NULL)
  }

  fit_df <- data.frame(Y = Y, X)

  t0 <- proc.time()
  sl_fit <- tryCatch(
    mlr3superlearner::mlr3superlearner(
      data = fit_df, target = "Y", library = lib,
      outcome_type = "binomial", folds = n_folds
    ),
    error = function(e) {
      cat(sprintf("    SL fit failed: %s\n", e$message))
      NULL
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(sl_fit)) return(NULL)

  preds <- tryCatch(as.numeric(predict(sl_fit, fit_df)), error = function(e) NULL)
  if (is.null(preds)) { cat("    predict failed\n"); return(NULL) }

  n_unique <- length(unique(round(preds, 6)))
  auc_val <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(Y, preds, quiet = TRUE))),
    error = function(e) NA_real_
  )

  wts <- tryCatch(sl_fit$weights, error = function(e) NULL)
  winner <- if (!is.null(wts)) names(which.max(wts)) else "?"

  degenerate <- n_unique <= 1
  cat(sprintf("    AUC=%.4f | %d unique preds | winner=%s | %.1fs%s\n",
              auc_val, n_unique, winner, elapsed,
              if (degenerate) " ** DEGENERATE **" else ""))

  if (!is.null(wts) && any(wts > 0.01)) {
    active <- wts[wts > 0.01]
    cat(sprintf("    Weights: %s\n",
                paste(sprintf("%s=%.3f", names(active), active), collapse=", ")))
  }

  data.frame(library = lib_name, n_folds = n_folds, auc = auc_val,
             n_unique_preds = n_unique, winner = winner, elapsed = elapsed,
             degenerate = degenerate, stringsAsFactors = FALSE)
}


# ── Run all tests ─────────────────────────────────────────────────────────
results <- list()

for (oc in names(prepped)) {
  p <- prepped[[oc]]
  cat(sprintf("\n\n========== %s (n=%d, events=%d, prev=%.1f%%) ==========\n",
              oc, length(p$Y), p$n_events, p$prev * 100))

  # 5-fold tests
  r <- test_library(p$X, p$Y, p$cluster_id, lib_no_mean, "no_mean_5f", 5)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  r <- test_library(p$X, p$Y, p$cluster_id, make_lib_weighted(p$prev), "weighted_5f", 5)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  r <- test_library(p$X, p$Y, p$cluster_id, lib_bart_only, "bart_only_5f", 5)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  # 3-fold tests (more training data per fold)
  r <- test_library(p$X, p$Y, p$cluster_id, lib_no_mean, "no_mean_3f", 3)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  r <- test_library(p$X, p$Y, p$cluster_id, make_lib_weighted(p$prev), "weighted_3f", 3)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  r <- test_library(p$X, p$Y, p$cluster_id, lib_simple, "simple_3f", 3)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }

  # 10-fold test (standard ML)
  r <- test_library(p$X, p$Y, p$cluster_id, lib_no_mean, "no_mean_10f", 10)
  if (!is.null(r)) { r$outcome <- oc; results[[length(results)+1]] <- r }
}


# ── Summary ───────────────────────────────────────────────────────────────
cat("\n\n==================== SUMMARY ====================\n\n")

if (length(results) > 0) {
  res_df <- bind_rows(results) |> arrange(outcome, desc(auc))

  for (oc in unique(res_df$outcome)) {
    sub <- res_df[res_df$outcome == oc, ]
    cat(sprintf("\n%s (events=%d, prev=%.1f%%):\n",
                oc, prepped[[oc]]$n_events, prepped[[oc]]$prev * 100))
    for (i in seq_len(nrow(sub))) {
      flag <- if (sub$degenerate[i]) " ** DEGENERATE **" else ""
      cat(sprintf("  %-20s AUC=%.4f  uniq=%3d  winner=%-20s %5.1fs%s\n",
                  sub$library[i], sub$auc[i], sub$n_unique_preds[i],
                  sub$winner[i], sub$elapsed[i], flag))
    }
  }

  outfile <- here("results", "tables", "sierraleone_model_tests.csv")
  write.csv(res_df, outfile, row.names = FALSE)
  cat(sprintf("\nSaved to %s\n", outfile))
} else {
  cat("No results generated.\n")
}
