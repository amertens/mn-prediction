# Workaround: Extract BART tree structure before serialization,
# use it to refit a fresh BART for prediction after deserialization.
#
# Since dbarts can't predict after serialization, the workaround is:
# 1. After fitting, extract yhat.train (posterior predictions on training data)
# 2. Store the training data X and the BART predictions
# 3. For "predict on new data", use a k-NN approach: find the nearest
#    training points to each new data point and use their BART predictions
#
# OR: simply refit BART on the fly when predict is needed (fast for small data)
library(dbarts)

set.seed(42)
n <- 500; p <- 30
X <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste0("x", 1:p)
Y <- rbinom(n, 1, plogis(0.8 * X[,1] - 0.5 * X[,2]))

cat("=== Workaround: Store BART training predictions + refit ===\n\n")

# Fit BART
fit <- dbarts::bart(X, Y, ntree = 100, verbose = FALSE, keeptrees = TRUE)

# Extract what survives serialization
yhat_train <- fit$yhat.train  # posterior draws matrix (MCMC iterations x observations)
cat(sprintf("yhat.train: %d iterations x %d observations\n", nrow(yhat_train), ncol(yhat_train)))

# The training predictions (column means after probit transform)
pred_train <- pnorm(colMeans(yhat_train))
cat(sprintf("Training pred range: [%.3f, %.3f]\n", min(pred_train), max(pred_train)))

# Store: X_train, Y_train, pred_train — all plain R objects
bart_cache <- list(
  X_train = X,
  Y_train = Y,
  pred_train = pred_train,
  ntree = 100
)

# Save/load the cache
tmp <- tempfile(fileext = ".rds")
saveRDS(bart_cache, tmp)
cache <- readRDS(tmp)

cat(sprintf("Cache file size: %.1f KB\n", file.info(tmp)$size / 1024))

# For prediction: refit BART from cached training data (fast!)
cat("\n=== Refit BART from cache for prediction ===\n")
t0 <- proc.time()
refit <- dbarts::bart(cache$X_train, cache$Y_train,
                       x.test = X[1:10, ],  # predict on new data during fitting
                       ntree = cache$ntree,
                       verbose = FALSE)
elapsed <- (proc.time() - t0)["elapsed"]

pred_refit <- pnorm(colMeans(refit$yhat.test))
cat(sprintf("Refit time: %.2f seconds\n", elapsed))
cat(sprintf("Refit predictions: %s\n", paste(round(pred_refit, 3), collapse = ", ")))

# Compare to original predictions on same data
pred_orig <- pnorm(colMeans(predict(fit, X[1:10, ])))
cat(sprintf("Original predictions: %s\n", paste(round(pred_orig, 3), collapse = ", ")))
cat(sprintf("Correlation: %.4f\n", cor(pred_refit, pred_orig)))

# Test with permuted data
X_perm <- X
X_perm[, 1] <- X_perm[sample(n), 1]
refit_perm <- dbarts::bart(cache$X_train, cache$Y_train,
                           x.test = X_perm[1:10, ],
                           ntree = cache$ntree,
                           verbose = FALSE)
pred_perm <- pnorm(colMeans(refit_perm$yhat.test))
cat(sprintf("Permuted predictions: %s\n", paste(round(pred_perm, 3), collapse = ", ")))
cat(sprintf("Changed from original: %s\n", !isTRUE(all.equal(pred_refit, pred_perm))))

# Full AUC comparison
library(pROC)
refit_full <- dbarts::bart(cache$X_train, cache$Y_train,
                           x.test = X,
                           ntree = cache$ntree,
                           verbose = FALSE)
pred_full <- pnorm(colMeans(refit_full$yhat.test))
auc_val <- as.numeric(auc(roc(Y, pred_full, quiet = TRUE)))
cat(sprintf("\nFull refit AUC: %.3f\n", auc_val))

# Compare to original AUC
pred_full_orig <- pnorm(colMeans(predict(fit, X)))
auc_orig <- as.numeric(auc(roc(Y, pred_full_orig, quiet = TRUE)))
cat(sprintf("Original AUC: %.3f\n", auc_orig))
cat(sprintf("AUC difference: %.4f\n", abs(auc_val - auc_orig)))

cat("\n=== Done ===\n")
