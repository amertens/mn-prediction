# =============================================================================
# sandbox/03_synthetic_control_knn.R
#
# Hypothesis H2: synthetic control via covariate-space kNN. For each held-out
# admin-2 polygon, find the K nearest training admin-2 polygons in
# (PCA-reduced) covariate space, predict prevalence as a kernel-weighted
# mean of their observed prevalences. No regression at all — purely
# nearest-neighbour matching.
#
# Rationale: predictions are bounded by the observed distribution in the
# matched set, so the level-extrapolation problem (large positive bias for
# Gambia child_iron predicted by training-set-mean ≈ 20%) cannot occur.
# A held-out admin-2 with covariates similar to high-iron-deficiency
# training polygons receives a high-prevalence prediction by construction.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

fit_predict_knn_synthetic <- function(train, test, vars,
                                        K = 8, n_pc = 6, kernel = "uniform") {
  vars <- screen_vars(train, vars)
  if (length(vars) < 3) return(NULL)
  imp <- impute_median(train, test, vars)
  Xtr <- scale(as.matrix(imp$train))
  Xte <- scale(as.matrix(imp$test),
                center = attr(Xtr, "scaled:center"),
                scale  = attr(Xtr, "scaled:scale"))
  Xtr[!is.finite(Xtr)] <- 0; Xte[!is.finite(Xte)] <- 0
  # PCA on training; project test.
  pr <- stats::prcomp(Xtr, center = FALSE, scale. = FALSE)
  k_pc <- min(n_pc, ncol(pr$x))
  Ztr <- pr$x[, seq_len(k_pc), drop = FALSE]
  Zte <- Xte %*% pr$rotation[, seq_len(k_pc), drop = FALSE]
  # For each test row, find K nearest training rows by Euclidean distance
  # in PCA space; predict prevalence as kernel-weighted mean.
  Y <- train$svy_prev
  preds <- numeric(nrow(Zte))
  for (i in seq_len(nrow(Zte))) {
    d <- sqrt(colSums((t(Ztr) - Zte[i, ])^2))
    nn_idx <- order(d)[seq_len(min(K, length(d)))]
    w <- switch(kernel,
                 uniform     = rep(1, length(nn_idx)),
                 inverse     = 1 / (d[nn_idx] + 1e-6),
                 gaussian    = exp(-d[nn_idx]^2 / (2 * median(d)^2)))
    preds[i] <- sum(w * Y[nn_idx]) / sum(w)
  }
  pmin(pmax(preds, 0), 1)
}

pooled_all <- load_all_pooled()

res <- list()
configs <- list(
  list(K = 3,  n_pc = 6, kernel = "uniform",  name = "knn_K3_PC6_unif"),
  list(K = 8,  n_pc = 6, kernel = "uniform",  name = "knn_K8_PC6_unif"),
  list(K = 8,  n_pc = 6, kernel = "inverse",  name = "knn_K8_PC6_inv"),
  list(K = 8,  n_pc = 6, kernel = "gaussian", name = "knn_K8_PC6_gauss"),
  list(K = 15, n_pc = 10, kernel = "inverse", name = "knn_K15_PC10_inv")
)
for (cfg in configs) {
  cat(sprintf("\n[H2] %s ...\n", cfg$name))
  r <- loco_evaluate_all(pooled_all,
    predict_fn  = function(train, test, vars)
      fit_predict_knn_synthetic(train, test, vars,
                                  K = cfg$K, n_pc = cfg$n_pc,
                                  kernel = cfg$kernel),
    method_name = cfg$name)
  res[[cfg$name]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "03_knn_synthetic.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "03_knn_synthetic.csv")))
