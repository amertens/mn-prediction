# Test if dbarts::bart2 with keepTrees=TRUE survives serialization
library(dbarts)

set.seed(42)
n <- 500; p <- 30
X <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste0("x", 1:p)
Y <- rbinom(n, 1, plogis(0.8 * X[,1] - 0.5 * X[,2]))

cat("=== Testing dbarts::bart2 with keepTrees ===\n\n")

# Fit with keepTrees (stores tree structure in R, not C++ pointer)
fit <- dbarts::bart2(X, Y, n.trees = 100L, keepTrees = TRUE, verbose = FALSE)
cat(sprintf("Class: %s\n", paste(class(fit), collapse = ", ")))

# Predict before save
pred_before <- colMeans(pnorm(predict(fit, X[1:10, ])))
cat(sprintf("Before save: %s\n", paste(round(pred_before, 3), collapse = ", ")))

# RDS save/load
tmp <- tempfile(fileext = ".rds")
saveRDS(fit, tmp)
fit2 <- readRDS(tmp)

pred_after <- tryCatch({
  p <- colMeans(pnorm(predict(fit2, X[1:10, ])))
  cat(sprintf("After RDS:  %s\n", paste(round(p, 3), collapse = ", ")))
  cat(sprintf("Match: %s\n", all.equal(pred_before, p)))
  p
}, error = function(e) {
  cat(sprintf("FAILED: %s\n", e$message))
  NULL
})

# Socket serialization
raw_bytes <- serialize(fit, NULL)
fit3 <- unserialize(raw_bytes)

pred_socket <- tryCatch({
  p <- colMeans(pnorm(predict(fit3, X[1:10, ])))
  cat(sprintf("After socket: %s\n", paste(round(p, 3), collapse = ", ")))
  cat(sprintf("Match: %s\n", all.equal(pred_before, p)))
  p
}, error = function(e) {
  cat(sprintf("Socket FAILED: %s\n", e$message))
  NULL
})

# Permutation test
if (!is.null(pred_after)) {
  X_perm <- X
  X_perm[, 1] <- X_perm[sample(n), 1]
  pred_perm <- colMeans(pnorm(predict(fit2, X_perm[1:10, ])))
  cat(sprintf("\nPermuted: %s\n", paste(round(pred_perm, 3), collapse = ", ")))
  cat(sprintf("Changed: %s\n", !isTRUE(all.equal(pred_after, pred_perm))))
}

# Check file size
cat(sprintf("\nFile size: %.1f KB\n", file.info(tmp)$size / 1024))

cat("\n=== Done ===\n")
