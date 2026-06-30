# =============================================================================
# R_corrected/00_corrected_utils.R
#
# Shared helpers for the PARALLEL "corrected" pipeline (_targets_corrected.R).
# Nothing here mutates the production pipeline. All function names carry no
# collision risk with R/ (distinct names / _corrected suffix).
#
# P8 (bug fixes) is enforced *by construction* in these helpers:
#   - every survey/area join is keyed on (country, Admin2) — never Admin2 alone
#     (prevents the cross-country cartesian blow-up hit during the review);
#   - any train/test standardization uses TRAINING statistics only (no held-out
#     leakage); out-of-sample imputation uses TRAINING medians, never the
#     target rows' own values.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ── Observability: structured, greppable skip logger ─────────────────────────
# Replaces silent `tryCatch(..., error = function(e) NULL)` traces with a single
# greppable line so a degraded run is visible in the log. Pure stdout (no shared
# file write) so it is safe under parallel tar_make_future(). Rollups
# additionally emit a per-field coverage table (see comparison.R).
.log_skip <- function(country, outcome, stage, reason) {
  cat(sprintf("[SKIP] country=%s outcome=%s stage=%s reason=%s\n",
              country %||% "NA", outcome %||% "NA", stage, reason))
  invisible(NULL)
}

# ── Read an already-built object from the PRODUCTION store ───────────────────
# The corrected pipeline reuses production data-loading outputs (merged data,
# per-outcome datasets, survey aggregates, GEE covariates) instead of
# recomputing the heavy GEE extraction. Production target objects are plain RDS.
PROD_STORE <- Sys.getenv("PROD_STORE", "_targets_full")

read_prod <- function(name, store = PROD_STORE) {
  p <- file.path(store, "objects", name)
  if (!file.exists(p)) {
    warning(sprintf("[corrected] production object not found: %s", p))
    return(NULL)
  }
  tryCatch(readRDS(p), error = function(e) {
    warning(sprintf("[corrected] failed reading %s: %s", name, e$message)); NULL
  })
}

# ── Candidate predictor pool ────────────────────────────────────────────────
# Proxy-only (no gw_ — see feedback_proxy_only). For the lightweight SMOKE we
# cap the pool by UNSUPERVISED variance ranking (no use of Y → not leakage); the
# full run passes max_n = Inf to keep the entire proxy set and let the in-fold
# supervised prescreen do the work.
select_candidates <- function(od, max_n = Inf) {
  d <- od$data
  cand <- intersect(od$Xvars, colnames(d))          # proxy-only
  num <- cand[vapply(cand, function(c) is.numeric(d[[c]]) || is.integer(d[[c]]),
                     logical(1))]
  keep <- num[vapply(num, function(c) {
    x <- d[[c]]
    isTRUE(mean(is.na(x)) < 0.5 && stats::sd(x, na.rm = TRUE) > 1e-8)
  }, logical(1))]
  if (is.finite(max_n) && length(keep) > max_n) {
    v <- vapply(keep, function(c) stats::var(d[[c]], na.rm = TRUE), numeric(1))
    keep <- names(sort(v, decreasing = TRUE))[seq_len(max_n)]
  }
  keep
}

# ── Fold builders (grouped / blocked, never naive-random except optimistic) ──
# Assign whole GROUPS (clusters or regions) to folds so members stay together.
make_group_folds <- function(group, V, seed = 1L) {
  g <- as.character(group)
  ug <- unique(g[!is.na(g)])
  V <- max(2L, min(V, length(ug)))
  set.seed(seed)
  fold_of_group <- setNames(sample(rep_len(seq_len(V), length(ug))), ug)
  out <- fold_of_group[g]
  out[is.na(out)] <- sample(seq_len(V), sum(is.na(out)), replace = TRUE)
  as.integer(out)
}

make_random_folds <- function(n, V, seed = 1L) {
  set.seed(seed); as.integer(sample(rep_len(seq_len(V), n)))
}

# ── Leakage-free preprocessing fit on TRAINING rows only ────────────────────
# Order: drop bad cols -> median impute (train) -> variance filter ->
# SUPERVISED prescreen top-K by |assoc with train Y| -> correlation filter ->
# standardize (train mean/sd). Returned object can transform any X identically.
fit_preproc <- function(Xtr, ytr, prescreen_k = 30L, cor_cut = 0.90) {
  Xtr <- as.data.frame(Xtr)
  cols <- colnames(Xtr)[vapply(Xtr, function(x)
    is.numeric(x) || is.integer(x), logical(1))]
  Xtr <- Xtr[, cols, drop = FALSE]
  # drop all-NA / constant
  ok <- vapply(Xtr, function(x) {
    nn <- x[!is.na(x)]; isTRUE(length(nn) > 1 && stats::sd(nn) > 1e-8)
  }, logical(1))
  Xtr <- Xtr[, ok, drop = FALSE]
  med <- vapply(Xtr, function(x) stats::median(x, na.rm = TRUE), numeric(1))
  for (j in seq_along(Xtr)) Xtr[is.na(Xtr[[j]]), j] <- med[j]
  # variance filter post-impute
  sds <- vapply(Xtr, stats::sd, numeric(1))
  Xtr <- Xtr[, sds > 1e-8, drop = FALSE]; med <- med[colnames(Xtr)]
  # SUPERVISED prescreen (train Y only) — top-K by |point-biserial corr|
  assoc <- vapply(colnames(Xtr), function(c)
    suppressWarnings(abs(stats::cor(Xtr[[c]], ytr, use = "complete.obs"))),
    numeric(1))
  assoc[is.na(assoc)] <- 0
  k <- min(prescreen_k, ncol(Xtr))
  sel <- names(sort(assoc, decreasing = TRUE))[seq_len(k)]
  Xtr <- Xtr[, sel, drop = FALSE]
  # correlation filter (train only)
  if (ncol(Xtr) > 1) {
    cm <- suppressWarnings(stats::cor(Xtr))
    cm[is.na(cm)] <- 0
    drop <- caret_find_correlation(cm, cor_cut)
    if (length(drop)) Xtr <- Xtr[, -drop, drop = FALSE]
  }
  fin <- colnames(Xtr)
  ctr <- vapply(Xtr, mean, numeric(1))
  scl <- vapply(Xtr, stats::sd, numeric(1)); scl[scl < 1e-8] <- 1
  list(final_cols = fin, medians = med[fin], centers = ctr, scales = scl)
}

# Minimal greedy correlation filter (avoids a hard caret dependency).
caret_find_correlation <- function(cm, cutoff) {
  p <- ncol(cm); if (p < 2) return(integer(0))
  cm2 <- abs(cm); diag(cm2) <- 0
  drop <- logical(p)
  avg <- rowMeans(cm2)
  repeat {
    m <- max(cm2[!drop, !drop, drop = FALSE])
    if (m < cutoff) break
    idx <- which(cm2 == m, arr.ind = TRUE)
    idx <- idx[1, , drop = TRUE]
    i <- idx[1]; j <- idx[2]
    rm_i <- if (avg[i] >= avg[j]) i else j
    drop[rm_i] <- TRUE
    cm2[rm_i, ] <- 0; cm2[, rm_i] <- 0
  }
  which(drop)
}

apply_preproc <- function(X, pp) {
  X <- as.data.frame(X)
  out <- matrix(0, nrow = nrow(X), ncol = length(pp$final_cols),
                dimnames = list(NULL, pp$final_cols))
  for (c in pp$final_cols) {
    x <- if (c %in% colnames(X)) X[[c]] else rep(NA_real_, nrow(X))
    x[is.na(x)] <- pp$medians[[c]]                  # TRAIN median (P8: no OOS leakage)
    out[, c] <- (x - pp$centers[[c]]) / pp$scales[[c]]
  }
  out
}

# ── Base learners (kept light; library is configurable) ─────────────────────
fit_learner <- function(type, Xtr, ytr) {
  switch(type,
    mean   = list(type = "mean", p = mean(ytr)),
    glmnet = {
      fit <- tryCatch(
        glmnet::cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = 5),
        error = function(e) NULL)
      list(type = "glmnet", fit = fit, p = mean(ytr))
    },
    rpart = {
      df <- data.frame(.y = factor(ytr), Xtr)
      fit <- tryCatch(rpart::rpart(.y ~ ., data = df, method = "class",
                                   cp = 0.005, minbucket = 20),
                      error = function(e) NULL)
      list(type = "rpart", fit = fit, p = mean(ytr))
    },
    ranger = {
      df <- data.frame(.y = factor(ytr), Xtr)
      fit <- tryCatch(ranger::ranger(.y ~ ., data = df, num.trees = 300,
                                     probability = TRUE, min.node.size = 10),
                      error = function(e) NULL)
      list(type = "ranger", fit = fit, p = mean(ytr))
    },
    stop("unknown learner ", type))
}

predict_learner <- function(m, Xte) {
  if (m$type == "mean") return(rep(m$p, nrow(Xte)))
  if (is.null(m$fit)) return(rep(m$p, nrow(Xte)))
  p <- switch(m$type,
    glmnet = as.numeric(predict(m$fit, newx = Xte, s = "lambda.min",
                                type = "response")),
    rpart  = {
      pr <- predict(m$fit, newdata = data.frame(Xte), type = "prob")
      if ("1" %in% colnames(pr)) pr[, "1"] else pr[, ncol(pr)]
    },
    ranger = {
      pr <- predict(m$fit, data = data.frame(Xte))$predictions
      if ("1" %in% colnames(pr)) pr[, "1"] else pr[, ncol(pr)]
    })
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

# ── Metrics ─────────────────────────────────────────────────────────────────
auc_manual <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

binary_metrics <- function(y, p) {
  y <- as.numeric(y); ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  prev <- mean(y)
  brier <- mean((p - y)^2)
  brier_base <- mean((prev - y)^2)
  data.frame(
    n = length(y), prevalence = round(prev, 4),
    auc = round(auc_manual(y, p), 4),
    brier = round(brier, 4),
    brier_skill = round(1 - brier / brier_base, 4)
  )
}
