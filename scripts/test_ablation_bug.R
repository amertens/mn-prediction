# Test the exact ablation flow to find the AUC=0.5 bug
library(here)
library(targets)
library(mlr3superlearner)
library(mlr3)
library(pROC)
source(here("R", "conceptual_ablation.R"))
source(here("R", "domain_ablation.R"))

cat("\n=== Reproducing ablation bug ===\n\n")

sl <- tar_read(sl_fit_sierraleone_child_iron, store = "_targets")
bf <- sl$bin_fit
rm(sl); gc()

covars <- bf$Xvars
Y <- as.numeric(bf$res$Y)
train_data <- data.frame(bf$train_data)
fit <- bf$sl_fit

cat(sprintf("Xvars: %d, n: %d, prevalence: %.3f\n", length(covars), length(Y), mean(Y)))

# Step 1: Verify predict works on original data
pred_orig <- as.numeric(fit$predict(train_data))
cat(sprintf("Predict on original: range [%.4f, %.4f], mean=%.4f\n",
            min(pred_orig), max(pred_orig), mean(pred_orig)))
cat(sprintf("AUC original: %.3f\n", as.numeric(auc(roc(Y, pred_orig, quiet=TRUE)))))

# Step 2: Classify variables
classified <- classify_variables(covars)
domains <- split(names(classified), classified)
domains <- domains[names(domains) != "Unknown"]
domains <- domains[names(domains) != "Missingness Indicators"]
cat(sprintf("\nClassified into %d domains\n", length(domains)))

# Step 3: Try permuting the largest domain
biggest <- names(which.max(sapply(domains, length)))
dom_cols <- intersect(domains[[biggest]], covars)
cat(sprintf("\nPermuting '%s' (%d vars)...\n", biggest, length(dom_cols)))

# Permute
set.seed(12346)
td_perm <- train_data
for (col in dom_cols) {
  td_perm[[col]] <- td_perm[[col]][sample(nrow(td_perm))]
}

pred_perm <- as.numeric(fit$predict(td_perm))
cat(sprintf("Predict permuted: range [%.4f, %.4f], mean=%.4f\n",
            min(pred_perm), max(pred_perm), mean(pred_perm)))
cat(sprintf("All 0.5? %s\n", all(abs(pred_perm - 0.5) < 0.001)))
auc_perm <- as.numeric(auc(roc(Y, pred_perm, quiet=TRUE)))
cat(sprintf("AUC permuted: %.3f (delta: %.4f)\n", auc_perm, auc_perm - as.numeric(auc(roc(Y, pred_orig, quiet=TRUE)))))

# Step 4: Now test using permute_domain() directly
cat("\n--- Testing permute_domain() ---\n")
result <- permute_domain(bf, dom_cols, biggest, use_binary = TRUE, n_perm = 3, seed = 12345L)
if (!is.null(result)) {
  cat(sprintf("permute_domain result: AUC=%.3f, Brier=%.4f\n", result$auc, result$brier))
} else {
  cat("permute_domain returned NULL\n")
}

cat("\n=== Done ===\n")
