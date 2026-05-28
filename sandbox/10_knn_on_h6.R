# =============================================================================
# sandbox/10_knn_on_h6.R
#
# Experimental: re-test kNN matching (H2 from sandbox/) but in the LOW-
# DIMENSIONAL feature space selected by the H6 combined filter. H2
# failed at d=152 due to curse-of-dimensionality. With H6 selecting
# the top-12 transportable variables, the nearest-neighbour matching
# might recover useful information.
#
# Not promoted to the main pipeline: kNN is a non-standard
# transportability method that requires the H6 feature selection
# upstream. Documented here as proof-of-concept.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

pooled_all <- load_all_pooled()
cross_pooled <- setNames(lapply(pooled_all, function(x) x$pooled_data),
                          names(pooled_all))

knn_predict_on_h6 <- function(this_outcome, K = 8, kernel = "inverse",
                                top_k_h6 = 12, min_within_ratio = 0.70) {
  function(train, test, vars) {
    # Step 1: get H6-selected variables on this fold's training data.
    out_h6 <- fit_predict_combined_filter(
      train, test, vars,
      cross_outcome_pooled = cross_pooled,
      this_outcome = this_outcome,
      min_within_ratio = min_within_ratio, top_k = top_k_h6)
    if (is.null(out_h6) || is.null(out_h6$selected_vars)) return(NULL)
    selected <- out_h6$selected_vars
    # Step 2: kNN matching in the H6-selected feature space.
    imp <- impute_median(train, test, selected)
    Xtr <- scale(as.matrix(imp$train))
    Xte <- scale(as.matrix(imp$test),
                  center = attr(Xtr, "scaled:center"),
                  scale  = attr(Xtr, "scaled:scale"))
    Xtr[!is.finite(Xtr)] <- 0; Xte[!is.finite(Xte)] <- 0
    Y <- train$svy_prev
    preds <- numeric(nrow(Xte))
    for (i in seq_len(nrow(Xte))) {
      d <- sqrt(colSums((t(Xtr) - Xte[i, ])^2))
      nn_idx <- order(d)[seq_len(min(K, length(d)))]
      w <- switch(kernel,
                   uniform  = rep(1, length(nn_idx)),
                   inverse  = 1 / (d[nn_idx] + 1e-6),
                   gaussian = exp(-d[nn_idx]^2 / (2 * median(d)^2)))
      preds[i] <- sum(w * Y[nn_idx]) / sum(w)
    }
    pmin(pmax(preds, 0), 1)
  }
}

res <- list()
configs <- list(
  list(K = 5,  kernel = "uniform",  name = "knn_on_h6_K5_unif"),
  list(K = 8,  kernel = "inverse",  name = "knn_on_h6_K8_inv"),
  list(K = 12, kernel = "gaussian", name = "knn_on_h6_K12_gauss")
)
for (cfg in configs) {
  cat(sprintf("\n[H6+kNN] %s ...\n", cfg$name))
  per_outcome <- list()
  for (oc in OUTCOMES) {
    fn <- knn_predict_on_h6(oc, K = cfg$K, kernel = cfg$kernel)
    r <- loco_evaluate(pooled_all[[oc]], fn, outcome = oc,
                        method_name = cfg$name)
    per_outcome[[oc]] <- r
  }
  res[[cfg$name]] <- dplyr::bind_rows(per_outcome)
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "10_knn_on_h6.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "10_knn_on_h6.csv")))
