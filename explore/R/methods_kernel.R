# =============================================================================
# explore/R/methods_kernel.R
#
# GENOMIC-PREDICTION METHODS FOR n = 14-87, p = 300+
#
# The project's standing conclusion is that "capacity is a liability at
# n = 14-87": every tuned learner loses to a zero-tuning domain index. That is
# the same regime quantitative genetics lives in (n individuals << p markers),
# and the field's answer there is NOT variable selection. It is to stop
# parameterising by coefficients and parameterise by the SIMILARITY BETWEEN
# UNITS instead: build a relationship matrix K = XX'/p and fit
#
#     y = 1 mu + sum_k u_k + e,    u_k ~ N(0, s2_k K_k),  e ~ N(0, s2_e W^-1)
#
# This is ridge regression in disguise, but with two properties that matter
# here:
#
#   1. The shrinkage is ESTIMATED by REML, not selected by cross-validation.
#      There is no tuning grid to overfit at n = 14, which is exactly the
#      failure mode the project has documented repeatedly.
#   2. With one kernel per domain, the REML variance components ARE a
#      decomposition of how much district-level variation each domain carries
#      - a direct, model-based answer to the question the domain-ablation work
#      has been approaching by deletion.
#
# Everything here is plain base R (no BGLR/rrBLUP dependency).
# =============================================================================

# ── kernels ─────────────────────────────────────────────────────────────────

#' Scale a kernel to mean diagonal 1 so variance components are comparable
k_scale <- function(K) {
  d <- mean(diag(K))
  if (!is.finite(d) || d <= 0) return(K)
  K / d
}

#' Linear (genomic-relationship style) kernel: K = XX' / p
k_linear <- function(X) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- 0
  if (ncol(X) == 0) return(NULL)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  k_scale(tcrossprod(Xc) / ncol(Xc))
}

#' Cosine kernel: the natural similarity for a learned embedding
k_cosine <- function(X) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- 0
  if (ncol(X) == 0) return(NULL)
  nrm <- sqrt(rowSums(X^2)); nrm[nrm == 0] <- 1
  k_scale(tcrossprod(X / nrm))
}

#' Gaussian kernel with the median-heuristic bandwidth
k_gaussian <- function(X, bw = NULL) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- 0
  if (ncol(X) == 0) return(NULL)
  D2 <- as.matrix(stats::dist(X))^2
  if (is.null(bw)) bw <- stats::median(D2[upper.tri(D2)])
  if (!is.finite(bw) || bw <= 0) bw <- 1
  k_scale(exp(-D2 / bw))
}

#' Exponential-decay spatial kernel on great-circle distance between centroids
#' Range defaults to the median pairwise distance (the spatial median heuristic).
k_spatial <- function(lon, lat, range_km = NULL) {
  R <- 6371
  rad <- pi / 180
  n <- length(lon)
  la <- lat * rad; lo <- lon * rad
  cx <- cos(la) * cos(lo); cy <- cos(la) * sin(lo); cz <- sin(la)
  P <- cbind(cx, cy, cz)
  ch <- pmin(pmax(tcrossprod(P), -1), 1)
  D <- R * acos(ch)
  if (is.null(range_km)) range_km <- stats::median(D[upper.tri(D)])
  if (!is.finite(range_km) || range_km <= 0) range_km <- 1
  k_scale(exp(-D / range_km))
}

#' One linear kernel per domain, from a column -> domain map
k_by_domain <- function(X, domain_of, domains = NULL, min_cols = 3L) {
  cols <- colnames(X)
  dm <- domain_of[cols]
  avail <- sort(unique(stats::na.omit(dm)))
  if (!is.null(domains)) avail <- intersect(avail, domains)
  out <- list()
  for (d in avail) {
    cc <- cols[which(dm == d)]
    if (length(cc) < min_cols) next
    K <- k_linear(X[, cc, drop = FALSE])
    if (!is.null(K)) out[[d]] <- K
  }
  out
}

# ── multi-kernel REML ───────────────────────────────────────────────────────

#' EM-REML for y = 1 mu + sum_k u_k + e,  u_k ~ N(0, s2_k K_k)
#'
#' EM updates are slower than AI-REML but cannot leave the parameter space,
#' which matters at n = 14 where AI-REML routinely proposes negative variances.
#' Variances are floored rather than dropped so a domain that carries nothing
#' collapses to ~0 instead of destabilising the fit.
#'
#' @param y numeric outcome (training rows only)
#' @param Klist list of kernels, each already subset to the training rows
#' @param w optional precision weights (n_eff); residual variance is s2_e/w_i
#' @return list(s2 = per-kernel variances, s2e, h2 = share of variance by
#'   kernel, mu, Vi = inverse of V, converged, iters)
reml_em <- function(y, Klist, w = NULL, tol = 1e-7, maxit = 300, floor_v = 1e-8) {
  n <- length(y)
  stopifnot(n > 2, length(Klist) > 0)
  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w); w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0], 1)
  w <- w / mean(w)
  Rm <- diag(1 / w, n, n)            # residual covariance structure

  vy <- stats::var(y)
  if (!is.finite(vy) || vy <= 0) vy <- 1
  nk <- length(Klist)
  s2 <- rep(vy / (nk + 1), nk)
  s2e <- vy / (nk + 1)
  X1 <- matrix(1, n, 1)

  conv <- FALSE; it <- 0L
  for (it in seq_len(maxit)) {
    V <- s2e * Rm
    for (k in seq_len(nk)) V <- V + s2[k] * Klist[[k]]
    Vi <- tryCatch(chol2inv(chol(V)), error = function(e)
      tryCatch(solve(V + diag(1e-6 * mean(diag(V)), n)), error = function(e2) NULL))
    if (is.null(Vi)) break
    XtVi <- crossprod(X1, Vi)
    A <- tryCatch(solve(XtVi %*% X1), error = function(e) NULL)
    if (is.null(A)) break
    P <- Vi - crossprod(XtVi, A %*% XtVi)
    Py <- P %*% y

    s2_new <- s2; s2e_new <- s2e
    for (k in seq_len(nk)) {
      q <- as.numeric(crossprod(Py, Klist[[k]] %*% Py))
      tr <- sum(P * Klist[[k]])                 # tr(P K) without forming P%*%K
      s2_new[k] <- max(floor_v, s2[k] + (s2[k]^2 / n) * (q - tr))
    }
    qe <- as.numeric(crossprod(Py, Rm %*% Py))
    tre <- sum(P * Rm)
    s2e_new <- max(floor_v, s2e + (s2e^2 / n) * (qe - tre))

    delta <- max(abs(c(s2_new, s2e_new) - c(s2, s2e))) / max(vy, 1e-12)
    s2 <- s2_new; s2e <- s2e_new
    if (delta < tol) { conv <- TRUE; break }
  }

  V <- s2e * Rm
  for (k in seq_len(nk)) V <- V + s2[k] * Klist[[k]]
  Vi <- tryCatch(chol2inv(chol(V)), error = function(e)
    solve(V + diag(1e-6 * mean(diag(V)), n)))
  mu <- as.numeric(solve(crossprod(X1, Vi) %*% X1, crossprod(X1, Vi) %*% y))
  tot <- sum(s2) + s2e
  names(s2) <- names(Klist)
  list(s2 = s2, s2e = s2e, h2 = s2 / tot, mu = mu, Vi = Vi,
       converged = conv, iters = it)
}

#' Exact REML for ONE kernel, by eigendecomposition (the EMMA/GEMMA trick)
#'
#' With a single kernel and i.i.d. residuals, V = s2e (delta K + I), so one
#' eigendecomposition of K turns the restricted likelihood into a 1-D function
#' of delta = s2u/s2e that can be optimised exactly. Orders of magnitude faster
#' than EM and not an approximation - which matters because most arms here are
#' single-kernel and get refitted 50 times per cell (5 folds x 10 replicates).
reml_1k <- function(y, K, lower = -10, upper = 10) {
  n <- length(y)
  eg <- eigen(K, symmetric = TRUE)
  d <- pmax(eg$values, 0)
  U <- eg$vectors
  ys <- as.numeric(crossprod(U, y))
  xs <- as.numeric(crossprod(U, rep(1, n)))

  nll <- function(ld) {
    delta <- exp(ld)
    wv <- delta * d + 1
    iw <- 1 / wv
    sxx <- sum(xs^2 * iw)
    if (!is.finite(sxx) || sxx <= 0) return(1e10)
    sxy <- sum(xs * ys * iw)
    mu <- sxy / sxx
    r <- ys - mu * xs
    s2e <- sum(r^2 * iw) / (n - 1)
    if (!is.finite(s2e) || s2e <= 0) return(1e10)
    0.5 * ((n - 1) * log(s2e) + sum(log(wv)) + log(sxx))
  }
  op <- stats::optimize(nll, lower = lower, upper = upper)
  delta <- exp(op$minimum)
  wv <- delta * d + 1
  iw <- 1 / wv
  sxx <- sum(xs^2 * iw); mu <- sum(xs * ys * iw) / sxx
  r <- ys - mu * xs
  s2e <- sum(r^2 * iw) / (n - 1)
  s2u <- delta * s2e
  # V^-1 = U diag(1/(s2e * wv)) U'
  Vi <- U %*% (iw / s2e * t(U))
  list(s2 = c(K = s2u), s2e = s2e, h2 = c(K = s2u / (s2u + s2e)),
       mu = mu, Vi = Vi, delta = delta, converged = TRUE, iters = 1L)
}

#' BLUP prediction for held-out rows
#'
#' pred_te = mu + sum_k s2_k K_k[te, tr] V_tr^-1 (y_tr - mu)
#'
#' Dispatches to the exact single-kernel solver when there is one kernel and no
#' weights, and to EM-REML otherwise.
#'
#' @param Kfull list of kernels on ALL rows (kernels are built from predictors
#'   only, so this uses no held-out outcome information)
reml_blup <- function(y, Kfull, tr, te, w = NULL) {
  Ktr <- lapply(Kfull, function(K) K[tr, tr, drop = FALSE])
  fit <- if (length(Kfull) == 1L && is.null(w))
    reml_1k(y[tr], Ktr[[1]])
  else
    reml_em(y[tr], Ktr, w = if (is.null(w)) NULL else w[tr])
  resid <- y[tr] - fit$mu
  a <- fit$Vi %*% resid
  pred <- rep(fit$mu, length(te))
  s2v <- as.numeric(fit$s2)
  for (k in seq_along(Kfull))
    pred <- pred + s2v[k] * as.numeric(Kfull[[k]][te, tr, drop = FALSE] %*% a)
  list(pred = pred, fit = fit)
}

# ── arms ────────────────────────────────────────────────────────────────────

#' Turn a kernel builder into a protocol-v2 arm
#'
#' @param kfun function(X, D, aux) -> named list of kernels on all rows
#' @param use_weights TRUE passes n_eff as precision weights to REML. Default
#'   FALSE, matching the protocol's own arms (domain_index is unweighted too),
#'   so a BLUP win cannot be an artefact of weighting the comparator lacks.
make_blup_arm <- function(kfun, use_weights = FALSE) {
  function(tr, te, y, X, D, aux) {
    Kf <- kfun(X, D, aux)
    Kf <- Kf[!vapply(Kf, is.null, TRUE)]
    if (!length(Kf)) return(rep(mean(y[tr]), length(te)))
    w <- if (use_weights) aux$w else NULL
    reml_blup(y, Kf, tr, te, w = w)$pred
  }
}

#' Variance components from one fit on the given rows, for reporting
#'
#' Reported on the TRAINING rows of a fold so the numbers are comparable with
#' the predictions they produced.
blup_varcomp <- function(y, Kfull, rows, w = NULL) {
  Ktr <- lapply(Kfull, function(K) K[rows, rows, drop = FALSE])
  fit <- reml_em(y[rows], Ktr, w = if (is.null(w)) NULL else w[rows])
  data.frame(kernel = c(names(fit$s2), "residual"),
             variance = c(as.numeric(fit$s2), fit$s2e),
             share = c(as.numeric(fit$h2), fit$s2e / (sum(fit$s2) + fit$s2e)),
             converged = fit$converged, iters = fit$iters,
             stringsAsFactors = FALSE)
}
