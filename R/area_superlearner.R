# =============================================================================
# R/area_superlearner.R
#
# A SuperLearner for AREA-LEVEL rows that honours survey weights and blocked
# cross-validation. Auto-sourced by tar_source("R/").
#
# WHY THIS EXISTS (2026-09-03)
# ----------------------------
# The production area-level SL (fit_predict_sl_prescreened, R/benchmark_models.R)
# and three sites in R/area_level_comparison.R call mlr3superlearner. Reading
# that package's source (0.1.2) shows three things that matter for area-level
# data and are not prominent in its documentation:
#
#   1. No observation-weight argument. Districts carry effective survey n from
#      single digits to several hundred; the SL fitted them all equally.
#   2. make_mlr3_resampling() returns plain rsmp("cv") for a regression task
#      BEFORE it checks the group column, so `group=` is silently ignored for
#      continuous outcomes. Blocking never happened.
#   3. The library is a hard whitelist (available_learners_regr), so a custom
#      learner such as the protocol-v2 domain index cannot be added.
#
# The classic SuperLearner package has obsWeights, custom fold assignment
# (cvControl$validRows) and plain-function wrappers, so this helper wraps it.
# Verified empirically in scripts/protocol_v2/19_sl_with_domain_index.R.
#
# CONTRACT
#   fit_area_superlearner(Y, X, newX, weights, block, library, V, discrete)
#   -> list(pred_new, pred_train, cvRisk, coef, pick, V_used, note) or NULL
#
# BLOCKING. `block` is a vector of group labels, one per training row (e.g.
# paste(country, Admin1)). Whole groups are assigned to folds, largest group
# first into the currently-smallest fold, after a seeded shuffle. Folds that
# end up empty are dropped; if fewer than 2 remain, plain random folds are used
# and the note says so, rather than failing silently.
#
# DISCRETE. The production SL used mlr3superlearner's default, discrete = TRUE,
# so the default here is the discrete pick (lowest CV risk); coef carries the
# NNLS ensemble weights either way so selection can be reported.
# =============================================================================

# ── Library wrappers ─────────────────────────────────────────────────────────
# Defined at top level so SuperLearner can find them by name. Internal CV fold
# counts are capped for small n (glmnet's default nfolds = 10 is more folds
# than some cells have rows).
.asl_nf <- function(n) max(3L, min(5L, floor(n / 5)))
# SuperLearner() passes `id` to every wrapper by name. Wrappers therefore take
# `id` explicitly and do NOT forward `...`, because SL.ranger hands its dots to
# ranger::ranger(), which rejects `id` -- the first version of this file did
# that and every fit failed silently.
SL.asl_lasso <- function(Y, X, newX, family, obsWeights, id, ...)
  SuperLearner::SL.glmnet(Y, X, newX, family, obsWeights, id, alpha = 1,
                          nfolds = .asl_nf(length(Y)))
SL.asl_enet  <- function(Y, X, newX, family, obsWeights, id, ...)
  SuperLearner::SL.glmnet(Y, X, newX, family, obsWeights, id, alpha = 0.5,
                          nfolds = .asl_nf(length(Y)))
SL.asl_ridge <- function(Y, X, newX, family, obsWeights, id, ...)
  SuperLearner::SL.glmnet(Y, X, newX, family, obsWeights, id, alpha = 0,
                          nfolds = .asl_nf(length(Y)))
SL.asl_ranger <- function(Y, X, newX, family, obsWeights, id, ...)
  SuperLearner::SL.ranger(Y, X, newX, family, obsWeights, num.trees = 250,
                          min.node.size = 5, verbose = FALSE)
SL.asl_xgb <- function(Y, X, newX, family, obsWeights, id, ...)
  SuperLearner::SL.xgboost(Y, X, newX, family, obsWeights, id, ntrees = 150,
                           max_depth = 4, shrinkage = 0.05,
                           params = list(subsample = 0.8, colsample_bytree = 0.8))
# SuperLearner() pairs every learner in a character library with the "All"
# screener and finds it by NAME in the calling environment. That works when the
# package is attached and fails with "object 'All' not found" when it is called
# as SuperLearner::SuperLearner() from a fresh session, so the screener is
# bound here next to the wrappers it accompanies.
All <- SuperLearner::All
SL.mean <- SuperLearner::SL.mean          # the one built-in passed by name
.ASL_LIBRARY <- c(mean = "SL.mean", lasso = "SL.asl_lasso", enet = "SL.asl_enet",
                  ridge = "SL.asl_ridge", ranger = "SL.asl_ranger",
                  xgb = "SL.asl_xgb")

# ── Rank-aligned meta-learner ────────────────────────────────────────────────
# method.NNLS minimises cross-validated SQUARED ERROR. At 14-87 noisy districts
# the MSE-optimal combination shrinks hard toward the mean -- in the head-to-head
# of 2026-09-03 the constant learner took the largest NNLS weight and the
# discrete SL picked the best-RANKING learner in 1% of fits. The decision this
# project makes is a ranking (which districts to reach), so this method scores
# each learner by Spearman correlation with the outcome and combines learners
# by non-negative least squares on RANKS. cvRisk = 1 - Spearman, so a discrete
# pick under this method is the best-ranking learner.
.asl_rank_coef <- function(Z, Y) {
  Z <- as.matrix(Z); K <- ncol(Z)
  rho <- vapply(seq_len(K), function(k) {
    z <- Z[, k]
    if (!all(is.finite(z)) || stats::sd(z) == 0) return(-1)
    suppressWarnings(stats::cor(z, Y, method = "spearman"))
  }, numeric(1))
  rho[!is.finite(rho)] <- -1
  coef_discrete <- as.numeric(seq_len(K) == which.max(rho))
  Zr <- apply(Z, 2, rank); Yr <- rank(Y)
  nn <- tryCatch(nnls::nnls(Zr, Yr)$x, error = function(e) coef_discrete)
  coef_nnls <- if (sum(nn) > 0) nn / sum(nn) else coef_discrete
  list(rho = rho, coef_discrete = coef_discrete, coef_nnls = coef_nnls)
}
method.asl_rank <- list(
  require = "nnls",
  computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
    r <- .asl_rank_coef(Z, Y)
    cvRisk <- 1 - r$rho; names(cvRisk) <- libraryNames
    coef <- r$coef_nnls; names(coef) <- libraryNames
    list(cvRisk = cvRisk, coef = coef, optimizer = "rank_nnls")
  },
  computePred = function(predY, coef, control, ...) as.numeric(as.matrix(predY) %*% coef)
)

# ── Population-weighted meta-learners ────────────────────────────────────────
# The project's decision is "which fifth of districts to reach", and districts
# differ in population by two orders of magnitude. Neither MSE nor plain rank
# loss knows that. SuperLearner's computeCoef receives only Z, Y and obsWeights
# (already carrying survey n), so population enters through FACTORIES that
# close over a per-row vector aligned with the training rows.

#' weighted Pearson correlation
.asl_wcor <- function(x, y, w) {
  mx <- stats::weighted.mean(x, w); my <- stats::weighted.mean(y, w)
  d <- sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
  if (!is.finite(d) || d == 0) return(NA_real_)
  sum(w * (x - mx) * (y - my)) / d
}

#' Population-weighted rank loss: weighted Spearman per learner (discrete) and
#' weighted NNLS on ranks (ensemble). w = population, normalised to mean 1.
.asl_wrank_coef <- function(Z, Y, w) {
  Z <- as.matrix(Z); K <- ncol(Z)
  w <- as.numeric(w); w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0], 1)
  w <- w / mean(w)
  Zr <- apply(Z, 2, rank); Yr <- rank(Y)
  rho <- vapply(seq_len(K), function(k)
    if (stats::sd(Z[, k]) == 0) -1 else .asl_wcor(Zr[, k], Yr, w), numeric(1))
  rho[!is.finite(rho)] <- -1
  coef_discrete <- as.numeric(seq_len(K) == which.max(rho))
  nn <- tryCatch(nnls::nnls(sqrt(w) * Zr, sqrt(w) * Yr)$x, error = function(e) coef_discrete)
  coef_nnls <- if (sum(nn) > 0) nn / sum(nn) else coef_discrete
  list(rho = rho, coef_discrete = coef_discrete, coef_nnls = coef_nnls)
}

#' Burden-capture loss: share of total burden (prevalence x population) inside
#' the top `frac` of rows ranked by the prediction. Piecewise constant in the
#' weights, so the ensemble is found by Nelder-Mead over softmax-parameterised
#' weights with restarts from the uniform point and every vertex; a 1e-3 x
#' Spearman term breaks the plateaus toward better overall ordering.
.asl_burden_coef <- function(Z, Y, burden, frac = 0.20) {
  Z <- as.matrix(Z); K <- ncol(Z); n <- nrow(Z)
  burden <- as.numeric(burden); burden[!is.finite(burden) | burden < 0] <- 0
  tot <- sum(burden); k <- max(1L, round(frac * n))
  cap <- function(p) { sel <- order(p, decreasing = TRUE)[seq_len(k)]; if (tot > 0) sum(burden[sel]) / tot else 0 }
  cap_k <- vapply(seq_len(K), function(j) if (stats::sd(Z[, j]) == 0) 0 else cap(Z[, j]), numeric(1))
  coef_discrete <- as.numeric(seq_len(K) == which.max(cap_k))
  obj <- function(theta) {
    b <- exp(theta - max(theta)); b <- b / sum(b); p <- as.numeric(Z %*% b)
    if (stats::sd(p) == 0) return(1)
    rho <- suppressWarnings(stats::cor(p, Y, method = "spearman")); if (!is.finite(rho)) rho <- 0
    -(cap(p) + 1e-3 * rho)
  }
  starts <- c(list(rep(0, K)), lapply(seq_len(K), function(j) { th <- rep(-3, K); th[j] <- 3; th }))
  best <- NULL
  for (s in starts) {
    o <- tryCatch(stats::optim(s, obj, method = "Nelder-Mead", control = list(maxit = 400)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  coef_ens <- if (is.null(best)) coef_discrete else { e <- exp(best$par - max(best$par)); e / sum(e) }
  list(cap = cap_k, coef_discrete = coef_discrete, coef_ens = coef_ens)
}

.asl_pred_lin <- function(predY, coef, control, ...) as.numeric(as.matrix(predY) %*% coef)

#' Factory: weighted-rank method with population `w` aligned to training rows.
make_method_asl_wrank <- function(w) list(
  require = "nnls",
  computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
    r <- .asl_wrank_coef(Z, Y, w)
    list(cvRisk = stats::setNames(1 - r$rho, libraryNames),
         coef = stats::setNames(r$coef_nnls, libraryNames), optimizer = "wrank_nnls")
  },
  computePred = .asl_pred_lin)

#' Factory: burden-capture method with `burden` (prevalence x population)
#' aligned to training rows. cvRisk = 1 - capture, so the discrete pick is the
#' learner whose own ranking captures the most burden.
make_method_asl_burden <- function(burden, frac = 0.20) list(
  require = NULL,
  computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
    r <- .asl_burden_coef(Z, Y, burden, frac)
    list(cvRisk = stats::setNames(1 - r$cap, libraryNames),
         coef = stats::setNames(r$coef_ens, libraryNames), optimizer = "burden_nm")
  },
  computePred = .asl_pred_lin)

#' Assign whole blocks to V folds, balanced by size.
#' @return integer fold id per row, or NULL if blocking cannot give >= 2 folds
.asl_block_folds <- function(block, V) {
  if (is.null(block)) return(NULL)
  block <- as.character(block)
  sizes <- sort(table(block), decreasing = TRUE)
  groups <- names(sizes)
  if (length(groups) < 2) return(NULL)
  # seeded shuffle among equal-sized groups so ties are not alphabetical
  groups <- groups[order(-as.numeric(sizes), stats::runif(length(sizes)))]
  V <- min(V, length(groups))
  load <- integer(V); assign <- integer(length(groups)); names(assign) <- groups
  for (g in groups) {
    k <- which.min(load); assign[g] <- k; load[k] <- load[k] + sizes[[g]]
  }
  fold <- unname(assign[block])
  if (length(unique(fold)) < 2) return(NULL)
  fold
}

#' Weighted, blocked SuperLearner for area-level rows.
#' @param Y numeric outcome (training rows)
#' @param X data.frame of predictors (training rows)
#' @param newX data.frame of predictors to predict (same columns)
#' @param weights per-row weights, e.g. n_svy; NULL for equal
#' @param block per-row group labels for fold blocking; NULL for random folds
#' @param library names from .ASL_LIBRARY, or full SL wrapper names
#' @param V number of CV folds for the meta-learner
#' @param discrete TRUE = predict with the single lowest-risk learner
#' @param meta "mse" = method.NNLS (squared-error meta-learner); "rank" =
#'   method.asl_rank (Spearman / rank-NNLS); "wrank" = population-weighted rank;
#'   "burden" = burden captured in the top fifth. "wrank"/"burden" need `pop`.
#' @param pop per-row population for the training rows (e.g. children under 5);
#'   used by meta = "wrank" (as weights) and "burden" (prevalence x pop). If
#'   missing, those metas fall back to "rank" and say so.
fit_area_superlearner <- function(Y, X, newX, weights = NULL, block = NULL,
                                  library = names(.ASL_LIBRARY), V = 5L,
                                  discrete = TRUE, family = stats::gaussian(),
                                  meta = c("mse", "rank", "wrank", "burden"),
                                  pop = NULL) {
  meta <- match.arg(meta)
  if (meta %in% c("wrank", "burden") && (is.null(pop) || length(pop) != length(Y))) {
    cat("    [area_superlearner] meta =", meta, "needs pop aligned to Y; using meta = rank
")
    meta <- "rank"
  }
  method <- switch(meta,
    mse    = SuperLearner::method.NNLS,
    rank   = method.asl_rank,
    wrank  = make_method_asl_wrank(pop),
    burden = make_method_asl_burden((if (all(Y >= 0 & Y <= 1)) Y else stats::plogis(Y)) * pop))
  .why <- function(msg) { cat("    [area_superlearner] returning NULL:", msg, "
"); NULL }
  if (!requireNamespace("SuperLearner", quietly = TRUE)) return(.why("SuperLearner not installed"))
  n <- length(Y)
  if (n < 6 || ncol(X) < 1) return(.why(sprintf("n=%d, p=%d", n, ncol(X))))
  lib <- ifelse(library %in% names(.ASL_LIBRARY), .ASL_LIBRARY[library], library)
  lib <- unname(lib)
  w <- if (is.null(weights)) rep(1, n) else pmax(as.numeric(weights), 1)
  w[!is.finite(w)] <- 1

  fold <- .asl_block_folds(block, V)
  note_folds <- if (is.null(fold)) {
    V_used <- max(2L, min(V, n)); "random folds"
  } else {
    V_used <- length(unique(fold)); sprintf("blocked folds (%d)", V_used)
  }
  cvc <- if (is.null(fold)) list(V = V_used)
         else list(V = V_used, validRows = unname(split(seq_len(n), fold)))

  # newX = rbind(X, newX) so library.predict covers training rows too
  both <- rbind(X, newX)
  fit <- tryCatch(suppressWarnings(SuperLearner::SuperLearner(
    Y = Y, X = X, newX = both, family = family, SL.library = lib,
    obsWeights = w, cvControl = cvc, method = method)),
    error = function(e) .why(paste("SuperLearner error:", conditionMessage(e))))
  if (is.null(fit)) return(NULL)

  lp <- fit$library.predict
  colnames(lp) <- sub("_All$", "", colnames(lp))
  risk <- fit$cvRisk; names(risk) <- sub("_All$", "", names(risk))
  coef <- fit$coef;   names(coef) <- sub("_All$", "", names(coef))
  ok <- is.finite(risk)
  if (!any(ok)) return(.why("no learner produced a finite CV risk"))
  pick <- names(risk)[ok][which.min(risk[ok])]
  p_all <- if (discrete) lp[, pick] else as.numeric(fit$SL.predict)
  if (!all(is.finite(p_all))) return(.why("non-finite predictions"))
  list(pred_train = p_all[seq_len(n)],
       pred_new   = p_all[n + seq_len(nrow(newX))],
       cvRisk = risk, coef = coef, pick = pick, V_used = V_used,
       meta = meta, Z = fit$Z, Y = Y,
       note = sprintf("%s; weights=%s; meta=%s; pick=%s", note_folds,
                      if (is.null(weights)) "equal" else "supplied", meta, pick))
}
