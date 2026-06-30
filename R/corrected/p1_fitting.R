# =============================================================================
# R_corrected/p1_fitting.R  — P1: leakage-free fitting + honest spatial CV
#
# Re-implements the binary deficiency fit so that ALL preprocessing (impute,
# variance filter, SUPERVISED prescreen, correlation filter, standardize) is
# fit INSIDE each training fold. Produces honest out-of-fold (OOF) predictions
# under two deployment-honest fold schemes (cluster-blocked and region/spatial-
# blocked) AND an "optimistic" counterpart that mirrors the production leakage
# (preprocessing fit on the FULL data + naive random CV) so the head-to-head
# delta is computed on the SAME data and SAME learner library.
#
# A compact DISCRETE SuperLearner is used: an inner CV on the training fold
# selects the lowest-loss base learner, which is refit on the full training
# fold and applied to the held-out fold. (The inner selection reuses the
# fold-fitted preprocessing — a deliberate, documented simplification; the
# OUTER honesty that defines the OOF metric is exact because the held-out fold
# never touches preprocessing or fitting.)
# =============================================================================

# One honest OOF pass for a given fold vector (preprocessing fit per train fold)
run_oof_honest <- function(d, predictors, y, fold, library, inner_V = 3L,
                           prescreen_k = 30L, seed = 12345L) {
  n <- nrow(d); oof <- rep(NA_real_, n); chosen <- character(0)
  for (f in sort(unique(fold))) {
    te <- which(fold == f); tr <- which(fold != f)
    if (length(te) == 0 || length(tr) < 20) next
    ytr <- y[tr]
    if (length(unique(ytr)) < 2) { oof[te] <- mean(ytr); next }
    pp <- fit_preproc(d[tr, predictors, drop = FALSE], ytr,
                      prescreen_k = prescreen_k)
    Xtr <- apply_preproc(d[tr, predictors, drop = FALSE], pp)
    Xte <- apply_preproc(d[te, predictors, drop = FALSE], pp)
    best <- select_best_learner(Xtr, ytr, library, inner_V, seed + f)
    chosen <- c(chosen, best)
    m <- fit_learner(best, Xtr, ytr)
    oof[te] <- predict_learner(m, Xte)
  }
  list(oof = data.frame(row = seq_len(n), Y = y, yhat_full = oof, fold = fold),
       chosen = chosen)
}

# Optimistic pass mirroring production leakage: preprocessing fit ONCE on the
# full data (incl. the held-out rows + Y for prescreen), then naive random CV.
run_oof_optimistic <- function(d, predictors, y, V = 5L, library,
                               inner_V = 3L, prescreen_k = 30L, seed = 12345L) {
  n <- nrow(d)
  pp <- fit_preproc(d[, predictors, drop = FALSE], y, prescreen_k = prescreen_k)
  X  <- apply_preproc(d[, predictors, drop = FALSE], pp)   # LEAK: full-data fit
  fold <- make_random_folds(n, V, seed)
  oof <- rep(NA_real_, n); chosen <- character(0)
  for (f in sort(unique(fold))) {
    te <- which(fold == f); tr <- which(fold != f)
    ytr <- y[tr]
    if (length(unique(ytr)) < 2) { oof[te] <- mean(ytr); next }
    best <- select_best_learner(X[tr, , drop = FALSE], ytr, library, inner_V,
                                seed + f)
    chosen <- c(chosen, best)
    m <- fit_learner(best, X[tr, , drop = FALSE], ytr)
    oof[te] <- predict_learner(m, X[te, , drop = FALSE])
  }
  list(oof = data.frame(row = seq_len(n), Y = y, yhat_full = oof, fold = fold),
       chosen = chosen)
}

# Inner-CV discrete SL selection (lowest 1-AUC; ties -> Brier).
select_best_learner <- function(X, y, library, inner_V = 3L, seed = 1L) {
  if (length(library) == 1) return(library)
  inner <- make_random_folds(nrow(X), inner_V, seed)
  loss <- setNames(rep(NA_real_, length(library)), library)
  for (lib in library) {
    pr <- rep(NA_real_, nrow(X))
    for (f in sort(unique(inner))) {
      tr <- which(inner != f); te <- which(inner == f)
      if (length(unique(y[tr])) < 2) { pr[te] <- mean(y[tr]); next }
      m <- fit_learner(lib, X[tr, , drop = FALSE], y[tr])
      pr[te] <- predict_learner(m, X[te, , drop = FALSE])
    }
    a <- auc_manual(y, pr)
    loss[lib] <- if (is.na(a)) mean((pr - y)^2) else 1 - a
  }
  names(which.min(loss))
}

# ── Top-level P1 target: fit all three schemes on one slice ─────────────────
fit_corrected_sl <- function(outcome_data, cc, oc, country_label,
                             library = c("mean", "glmnet", "rpart"),
                             max_pred = 150L, V = 5L, seed = 12345L) {
  d <- outcome_data$data
  ybin <- oc$binary
  if (is.null(ybin) || !ybin %in% colnames(d))
    return(list(error = "binary outcome column missing"))
  y <- as.numeric(d[[ybin]])
  keep <- !is.na(y)
  d <- d[keep, , drop = FALSE]; y <- y[keep]
  predictors <- select_candidates(outcome_data, max_n = max_pred)
  if (length(predictors) < 2) return(list(error = "too few predictors"))

  # Per-individual survey weight, aligned to d (Issue 2): the honest OOF
  # predictions are aggregated to Admin-2 with THESE weights so the predicted
  # district prevalence shares a weighting basis with the design-weighted
  # `svy_prev` truth. Non-finite / non-positive weights fall back to 1.
  wcol  <- cc$weight_col %||% NA_character_
  w_ind <- if (!is.na(wcol) && wcol %in% colnames(d))
    suppressWarnings(as.numeric(d[[wcol]])) else rep(1, nrow(d))
  w_ind[!is.finite(w_ind) | w_ind <= 0] <- 1

  clu <- cc$cluster_id %||% "gw_cnum"
  a1  <- cc$admin1_col %||% "Admin1"
  a2  <- cc$admin2_col %||% "Admin2"
  cluster_vec <- if (clu %in% colnames(d)) d[[clu]] else seq_len(nrow(d))
  region_vec  <- if (a1  %in% colnames(d)) d[[a1]]  else rep("all", nrow(d))

  fold_cluster <- make_group_folds(cluster_vec, V, seed)
  n_reg <- length(unique(region_vec[!is.na(region_vec)]))
  fold_region  <- make_group_folds(region_vec, min(V + 1L, n_reg), seed)

  h_clu <- run_oof_honest(d, predictors, y, fold_cluster, library, seed = seed)
  h_reg <- run_oof_honest(d, predictors, y, fold_region,  library, seed = seed)
  opt   <- run_oof_optimistic(d, predictors, y, V, library, seed = seed)

  attach_meta <- function(x) {
    x$oof$Admin2  <- if (a2 %in% colnames(d)) d[[a2]] else NA_character_
    x$oof$Admin1  <- region_vec
    x$oof$country <- country_label
    x$oof$w       <- w_ind            # per-individual survey weight (Issue 2)
    x
  }
  list(
    country = country_label, outcome = oc$tag %||% oc$label,
    n = nrow(d), n_predictors = length(predictors),
    honest_cluster   = attach_meta(h_clu),
    honest_region    = attach_meta(h_reg),
    optimistic       = attach_meta(opt),
    library = library
  )
}

# Tidy cv_perf_corrected: one row per (scheme) with honest/optimistic flag.
extract_cv_perf_corrected <- function(slfit) {
  if (!is.null(slfit$error)) return(NULL)
  mk <- function(tag, honest, x) {
    m <- binary_metrics(x$oof$Y, x$oof$yhat_full)
    data.frame(country = slfit$country, outcome = slfit$outcome,
               scheme = tag, honest = honest,
               top_learner = paste(sort(unique(x$chosen)), collapse = "+"),
               m, stringsAsFactors = FALSE)
  }
  rbind(
    mk("cluster-block CV", TRUE,  slfit$honest_cluster),
    mk("region/spatial-block CV", TRUE, slfit$honest_region),
    mk("random CV + full-data preprocessing (optimistic)", FALSE, slfit$optimistic)
  )
}
