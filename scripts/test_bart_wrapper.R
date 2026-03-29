# Test the BART refit-on-predict wrapper end-to-end
# Simulates what happens in the targets pipeline:
# 1. Fit SL with BART → BART gets weight
# 2. Serialize to RDS (simulates targets save)
# 3. Deserialize (simulates targets load in ablation worker)
# 4. Predict via wrapper (should detect degenerate, fall back to refit)
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

# Source the wrapper
source(here("R", "mlr3_fitting.R"))

cat("\n=== Testing BART Refit-on-Predict Wrapper ===\n\n")

set.seed(42)
n <- 400; p <- 20
d <- data.frame(matrix(rnorm(n * p), ncol = p))
colnames(d) <- paste0("x", 1:p)
d$Y <- rbinom(n, 1, plogis(0.8 * d$x1 - 0.5 * d$x2 + 0.3 * d$x3))
d$cluster_id <- rep(1:(n/10), each = 10)  # 40 clusters of 10

cat(sprintf("Data: n=%d, p=%d, prev=%.2f\n", n, p, mean(d$Y)))

# Fit using our wrapper (which extracts learner_models and stores train_data)
cat("\nFitting mlr3_SL_clustered with BART...\n")
result <- mlr3_SL_clustered(
  d = d, Xvars = paste0("x", 1:p), outcome = "Y",
  population = "all", id = "cluster_id",
  folds = 3,
  mlr3_library = list(
    "mean",
    list("glmnet", alpha = 1, id = "lasso"),
    list("bart", ntree = 50, id = "bart_50")
  ),
  outcome_type = "binomial",
  prescreen = FALSE
)

if (is.null(result)) {
  cat("Fit failed!\n"); q()
}

cat(sprintf("CV AUC: %.3f\n", as.numeric(auc(roc(result$res$Y, result$res$yhat_full, quiet = TRUE)))))

# Check what was extracted
wrapper <- result$sl_fit
cat(sprintf("\nHas learner_models: %s\n", !is.null(wrapper$learner_models)))
cat(sprintf("Has sl_weights: %s\n", !is.null(wrapper$sl_weights)))
cat(sprintf("Has train_data_for_bart: %s\n", !is.null(wrapper$train_data_for_bart)))

if (!is.null(wrapper$sl_weights)) {
  for (i in seq_along(wrapper$learner_models)) {
    lm <- wrapper$learner_models[[i]]
    w <- wrapper$sl_weights[i]
    cat(sprintf("  Learner %d: %s, weight=%.3f, is_bart=%s\n",
                i, lm$id, w, isTRUE(lm$is_bart)))
  }
}

# Test predict before serialization
pred_before <- as.numeric(wrapper$predict(d[1:10, ]))
cat(sprintf("\nBefore save: range [%.3f, %.3f], unique=%d\n",
            min(pred_before), max(pred_before), length(unique(round(pred_before, 6)))))

# Serialize and deserialize (simulate targets)
cat("\nSerializing via RDS...\n")
tmp <- tempfile(fileext = ".rds")
saveRDS(result, tmp)
loaded <- readRDS(tmp)
wrapper2 <- loaded$sl_fit

cat("Predicting from loaded model...\n")
pred_after <- as.numeric(wrapper2$predict(d[1:10, ]))
cat(sprintf("After load: range [%.3f, %.3f], unique=%d\n",
            min(pred_after), max(pred_after), length(unique(round(pred_after, 6)))))

if (length(unique(round(pred_after, 6))) > 1) {
  cat(sprintf("Correlation with before: %.4f\n", cor(pred_before, pred_after)))
  cat("SUCCESS: Predict works after serialization!\n")
} else {
  cat("FAILED: Predictions are degenerate after serialization\n")
}

# Test permutation
cat("\nTesting permutation sensitivity...\n")
d_perm <- d
d_perm$x1 <- d_perm$x1[sample(n)]
pred_perm <- as.numeric(wrapper2$predict(d_perm[1:10, ]))
cat(sprintf("Permuted: range [%.3f, %.3f]\n", min(pred_perm), max(pred_perm)))
cat(sprintf("Changed from unpermuted: %s\n", !isTRUE(all.equal(pred_after, pred_perm))))

# Full AUC
pred_full <- as.numeric(wrapper2$predict(d))
pred_full_perm <- as.numeric(wrapper2$predict(d_perm))
auc_orig <- as.numeric(auc(roc(d$Y, pred_full, quiet = TRUE)))
auc_perm <- as.numeric(auc(roc(d$Y, pred_full_perm, quiet = TRUE)))
cat(sprintf("Full AUC: %.3f, Permuted AUC: %.3f, Delta: %.4f\n",
            auc_orig, auc_perm, auc_perm - auc_orig))

cat("\n=== Done ===\n")
