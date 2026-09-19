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

# ── Domain-PC learners (SL-06, 2026-09-17) ───────────────────────────────────
# The protocol-v2 representation as SuperLearner candidates: each learner
# rebuilds the per-domain principal components (build_domain_pcs_v2, rotations
# from ITS training rows only) from the raw columns it is handed, then fits
#   SL.asl_domain_index     the zero-tuning domain index (arm_domain_index_v2)
#   SL.asl_domain_pc_ridge  weighted ridge on the domain PCs (glmnet, alpha 0)
#   SL.asl_domain_pc_enet   weighted elastic net on the domain PCs (alpha 0.5)
# They need a column -> domain map, which SuperLearner cannot pass through its
# wrapper signature, so fit_area_superlearner(domain_of = ...) parks it in
# .asl_domain_env for the duration of the call. X must already be the
# rank-normalised matrix (prep_predictors_v2), as in every protocol-v2 caller;
# columns without a domain (lon, lat, ...) are dropped by build_domain_pcs_v2.
# Predictions are made through newX at fit time (fit_area_superlearner always
# passes newX = rbind(X, newX)); predict() on the fit object is not supported.
.asl_domain_env <- new.env(parent = emptyenv())
.asl_domain_pcs <- function(X, newX) {
  dom <- .asl_domain_env$domain_of
  if (is.null(dom)) stop("domain-PC learner called without a domain map: pass domain_of= to fit_area_superlearner()")
  Xall <- rbind(as.matrix(X), as.matrix(newX)); storage.mode(Xall) <- "double"
  Xall[!is.finite(Xall)] <- 0
  D <- domain_representation_v2(Xall, dom[colnames(Xall)], sign_rows = seq_len(nrow(X)))
  if (!ncol(D)) stop("domain-PC learner: no column of X has a domain")
  list(tr = D[seq_len(nrow(X)), , drop = FALSE], te = D[nrow(X) + seq_len(nrow(newX)), , drop = FALSE])
}
.asl_domain_fit <- function(pred, what) {
  fit <- list(what = what); class(fit) <- "SL.asl_domain"
  list(pred = as.numeric(pred), fit = fit)
}
predict.SL.asl_domain <- function(object, newdata, ...)
  stop("SL.asl_domain_", object$what, " predicts via newX at fit time")
SL.asl_domain_index <- function(Y, X, newX, family, obsWeights, id, ...) {
  if (!isTRUE(.asl_domain_env$warned_weights) && !is.null(obsWeights) && length(unique(obsWeights)) > 1) {
    cat("    [SL.asl_domain_index] observation weights are ignored (the index has none)\n")
    .asl_domain_env$warned_weights <- TRUE
  }
  D <- .asl_domain_pcs(X, newX); ntr <- nrow(D$tr); nte <- nrow(D$te)
  p <- arm_domain_index_v2(seq_len(ntr), ntr + seq_len(nte), c(Y, rep(NA_real_, nte)),
                           NULL, rbind(D$tr, D$te), NULL)
  .asl_domain_fit(p, "index")
}
.asl_domain_glmnet <- function(Y, X, newX, family, obsWeights, id, alpha, what) {
  D <- .asl_domain_pcs(X, newX)
  if (ncol(D$tr) < 2) return(SuperLearner::SL.mean(Y, X, newX, family, obsWeights, id))
  g <- SuperLearner::SL.glmnet(Y, as.data.frame(D$tr), as.data.frame(D$te), family, obsWeights, id,
                               alpha = alpha, nfolds = .asl_nf(length(Y)))
  .asl_domain_fit(g$pred, what)
}
SL.asl_domain_pc_ridge <- function(Y, X, newX, family, obsWeights, id, ...)
  .asl_domain_glmnet(Y, X, newX, family, obsWeights, id, alpha = 0,   what = "pc_ridge")
SL.asl_domain_pc_enet  <- function(Y, X, newX, family, obsWeights, id, ...)
  .asl_domain_glmnet(Y, X, newX, family, obsWeights, id, alpha = 0.5, what = "pc_enet")
# Ordinary least squares on the first component of every domain: the
# "traditional regression" candidate, ~24 columns, no penalty. Rank-deficient
# fits (inner folds of the smallest cells) drop aliased columns, as lm() does.
SL.asl_domain_pc1_ols <- function(Y, X, newX, family, obsWeights, id, ...) {
  D <- .asl_domain_pcs(X, newX)
  pc1 <- grep("_PC1$", colnames(D$tr), value = TRUE)
  if (length(pc1) < 2) return(SuperLearner::SL.mean(Y, X, newX, family, obsWeights, id))
  dtr <- data.frame(Y = Y, D$tr[, pc1, drop = FALSE]); dte <- as.data.frame(D$te[, pc1, drop = FALSE])
  fit <- stats::lm(Y ~ ., data = dtr, weights = obsWeights)
  p <- suppressWarnings(as.numeric(stats::predict(fit, dte)))
  p[!is.finite(p)] <- mean(Y)
  .asl_domain_fit(p, "pc1_ols")
}
.ASL_DOMAIN_LEARNERS <- c("domain_index", "domain_pc_ridge", "domain_pc_enet", "domain_pc1_ols")

# ── Geostatistical learners (SL-06) ──────────────────────────────────────────
# The protocol's spatial arms as candidates. They read `lon` and `lat` columns
# of X; pair them with screen.asl_coords / All through the `screens` argument
# so the covariate-only learners never see coordinates. Not meaningful under
# transport (no observed outcome inside the held-out country): the caller
# leaves them out of the library there, as arms_for_estimand_v2() does.
screen.asl_coords    <- function(Y, X, family, obsWeights, id, ...) colnames(X) %in% c("lon", "lat")
screen.asl_no_coords <- function(Y, X, family, obsWeights, id, ...) !(colnames(X) %in% c("lon", "lat"))
.asl_coords_of <- function(X, newX) {
  if (!all(c("lon", "lat") %in% colnames(X))) stop("spatial learner needs lon and lat columns in X")
  list(lon = c(X$lon, newX$lon), lat = c(X$lat, newX$lat))
}
SL.asl_spatial_gam <- function(Y, X, newX, family, obsWeights, id, ...) {
  aux <- .asl_coords_of(X, newX); ntr <- nrow(X)
  p <- .v2_spatial_fit(seq_len(ntr), c(Y, rep(NA_real_, nrow(newX))), aux)(ntr + seq_len(nrow(newX)))
  .asl_domain_fit(p, "spatial_gam")
}
SL.asl_spatial_plus_domain <- function(Y, X, newX, family, obsWeights, id, ...) {
  aux <- .asl_coords_of(X, newX); ntr <- nrow(X); nte <- nrow(newX)
  keep <- !(colnames(X) %in% c("lon", "lat"))
  D <- .asl_domain_pcs(X[, keep, drop = FALSE], newX[, keep, drop = FALSE])
  sp <- .v2_spatial_fit(seq_len(ntr), c(Y, rep(NA_real_, nte)), aux)
  add <- .v2_enet(D$tr, Y - sp(seq_len(ntr)), D$te)
  .asl_domain_fit(sp(ntr + seq_len(nte)) + add, "spatial_plus_domain")
}
# SPDE Gaussian field + domain PC1s as fixed effects (R/protocol_v2_mbg.R,
# MB-01), at the rows' own coordinates. INLA: seconds per fit, so opt in.
SL.asl_mbg <- function(Y, X, newX, family, obsWeights, id, ...) {
  if (!exists(".mbg_arm") || !.mbg_ok()) stop("SL.asl_mbg needs INLA, fmesher and R/protocol_v2_mbg.R")
  aux <- .asl_coords_of(X, newX); ntr <- nrow(X); nte <- nrow(newX)
  keep <- !(colnames(X) %in% c("lon", "lat"))
  D <- .asl_domain_pcs(X[, keep, drop = FALSE], newX[, keep, drop = FALSE])
  p <- .mbg_arm(seq_len(ntr), ntr + seq_len(nte), c(Y, rep(NA_real_, nte)), NULL, rbind(D$tr, D$te), aux)
  .asl_domain_fit(p, "mbg")
}
.ASL_SPATIAL_LEARNERS <- c("spatial_gam", "spatial_plus_domain", "mbg")

.ASL_LIBRARY <- c(mean = "SL.mean", lasso = "SL.asl_lasso", enet = "SL.asl_enet",
                  ridge = "SL.asl_ridge", ranger = "SL.asl_ranger",
                  xgb = "SL.asl_xgb",
                  domain_index = "SL.asl_domain_index",
                  domain_pc_ridge = "SL.asl_domain_pc_ridge",
                  domain_pc_enet = "SL.asl_domain_pc_enet",
                  domain_pc1_ols = "SL.asl_domain_pc1_ols",
                  spatial_gam = "SL.asl_spatial_gam",
                  spatial_plus_domain = "SL.asl_spatial_plus_domain",
                  mbg = "SL.asl_mbg",
                  hapc = "SL.hapc", hapc_lasso = "SL.hapc_lasso")   # R/sl_hapc.R
# the six learners the production SL has always used; the domain-PC learners
# are opt-in because they need a domain map
.ASL_LIBRARY_DEFAULT <- c("mean", "lasso", "enet", "ridge", "ranger", "xgb")

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
#' @param screens optional named character vector, library entry -> screener
#'   function name (e.g. c(lasso = "screen.asl_no_coords", spatial_gam =
#'   "screen.asl_coords")); entries not named get "All". Lets coordinate
#'   columns reach only the spatial learners.
#' @param domain_of named character vector, column of X -> domain, required by
#'   the domain-PC learners (domain_index, domain_pc_ridge, domain_pc_enet).
#'   NULL drops those learners from `library` with a note, so callers on a
#'   vocabulary without a domain map (the legacy gee_ set) are unaffected.
#' @return list(pred_train, pred_new, cvRisk, coef, pick, V_used, meta, Z, Y,
#'   library_pred_train, library_pred_new, note) or NULL. The library_pred_*
#'   matrices hold every learner's prediction so other meta-learners can be
#'   derived from one fit (.asl_rank_coef(Z, Y) %*% library_pred_new).
fit_area_superlearner <- function(Y, X, newX, weights = NULL, block = NULL,
                                  library = .ASL_LIBRARY_DEFAULT, V = 5L,
                                  discrete = TRUE, family = stats::gaussian(),
                                  meta = c("mse", "rank", "wrank", "burden"),
                                  pop = NULL, domain_of = NULL, screens = NULL) {
  meta <- match.arg(meta)
  if (is.null(domain_of)) {
    dropped <- intersect(library, c(.ASL_DOMAIN_LEARNERS, "spatial_plus_domain", "mbg"))
    if (length(dropped)) {
      cat("    [area_superlearner] no domain_of map: dropping", paste(dropped, collapse = ", "), "\n")
      library <- setdiff(library, dropped)
    }
  } else {
    old_map <- .asl_domain_env$domain_of
    .asl_domain_env$domain_of <- domain_of
    on.exit(.asl_domain_env$domain_of <- old_map, add = TRUE)
  }
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
  if (!is.null(screens)) {
    sc <- ifelse(library %in% names(screens), screens[library], "All")
    lib <- lapply(seq_along(lib), function(i) c(lib[i], unname(sc[i])))
  }
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

  # library columns are named by the entries the caller passed (short names
  # from .ASL_LIBRARY, or wrapper names), in library order
  lp <- fit$library.predict
  lib_names <- if (ncol(lp) == length(library)) library else sub("_All$", "", colnames(lp))
  colnames(lp) <- lib_names
  risk <- fit$cvRisk; names(risk) <- lib_names
  coef <- fit$coef;   names(coef) <- lib_names
  Z <- fit$Z; if (!is.null(Z) && ncol(Z) == length(lib_names)) colnames(Z) <- lib_names
  ok <- is.finite(risk)
  if (!any(ok)) return(.why("no learner produced a finite CV risk"))
  pick <- names(risk)[ok][which.min(risk[ok])]
  p_all <- if (discrete) lp[, pick] else as.numeric(fit$SL.predict)
  if (!all(is.finite(p_all))) return(.why("non-finite predictions"))
  list(pred_train = p_all[seq_len(n)],
       pred_new   = p_all[n + seq_len(nrow(newX))],
       cvRisk = risk, coef = coef, pick = pick, V_used = V_used,
       meta = meta, Z = Z, Y = Y,
       library_pred_train = lp[seq_len(n), , drop = FALSE],
       library_pred_new   = lp[n + seq_len(nrow(newX)), , drop = FALSE],
       note = sprintf("%s; weights=%s; meta=%s; pick=%s", note_folds,
                      if (is.null(weights)) "equal" else "supplied", meta, pick))
}
