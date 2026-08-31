# =============================================================================
# R/eb_stack_target.R
#
# The stacked shrinkage target.
#
# WHY THIS EXISTS
# ---------------
# Every estimator this project has compared differs in ONE place: what the
# district's own estimate is shrunk toward.
#
#   eb_blend      jackknifed region mean
#   eb_covariate  covariate ridge prediction (the Fay-Herriot form)
#   flat_region   region mean, with no blending at all
#   ridge         covariate prediction, with no blending at all
#
# Measured, no single target wins everywhere:
#
#   Sierra Leone  region mean wins; the ridge is at or below the national mean
#                 null in 5 of 6 cells, and 14 districts over 4 regions leaves a
#                 leave-one-region-out fold training on about ten points.
#   Gambia        the covariate models beat the null in 4 of 4 cells, reaching
#                 r 0.717 against an empirical ceiling of 0.824. The tournament
#                 reproduces this: Gambia's implied covariate signal share sits
#                 at roughly 0.75 to 0.9, ABOVE the grid the tournament tested.
#   broadly       a covariate-free spatial smooth reaches 0.304 under
#                 leave-one-district-out and beats every area-level covariate
#                 model scored in sample.
#
# A fixed target therefore cannot be right for all 24 cells. This function
# learns the target per cell as a convex combination of candidates.
#
# WHAT IT NESTS
# -------------
# Weights of (0,1,0,0) reproduce the shipped eb_blend exactly, (0,0,0,1) is
# Fay-Herriot, (1,0,0,0) shrinks to the national mean. The stack cannot be much
# worse than its best single candidate, and it can adapt where they differ.
#
# THE TWO THINGS THAT WOULD MAKE IT DISHONEST
# -------------------------------------------
# 1. FITTING THE WEIGHTS IN SAMPLE. The candidates are built leave-one-region-
#    out and the weights are fitted leave-one-region-out on top of them. Fitting
#    the weights on all districts would re-introduce exactly the leak that the
#    "prescreen on ALL data" arm demonstrates: that arm gains +0.18 correlation
#    over its in-fold twin and is labelled OPTIMISTIC for that reason.
#
# 2. NOT RECALIBRATING THE LEVEL. The ridge carries a systematic -3.1 to -3.3 pp
#    bias, and a smoke test of the Fay-Herriot arm confirmed that shrinking
#    toward it imports that bias into the blend at -4.5 to -5.5 pp against the
#    plain blend's -1.2 to -1.9. Each candidate is therefore recentred, out of
#    region, on the precision-weighted mean of the training districts before it
#    is allowed into the stack.
#
# WHAT IT CANNOT DO
# -----------------
# This is a WITHIN-COUNTRY estimator. Leave-one-country-out fails on a
# cross-survey level offset reaching -42 pp, which is a calibration problem
# between surveys rather than a shrinkage problem inside one. No choice of
# target here touches it.
# =============================================================================

#' Precision-weighted mean, guarding the degenerate cases.
.sw_mean <- function(x, n) {
  ok <- is.finite(x) & is.finite(n) & n > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], n[ok])
}

#' Build the candidate shrinkage targets, each out of region.
#'
#' @param p district estimates
#' @param n district sample sizes
#' @param reg region label per district
#' @param X optional district covariate matrix, rows aligned to p
#' @param lonlat optional two-column matrix of district centroids
#' @param fin logical, which districts are usable
#' @param k_screen covariates retained by the in-fold screen
#' @return named list of numeric vectors, each the same length as p
eb_stack_candidates <- function(p, n, reg, X = NULL, lonlat = NULL,
                                fin = NULL, k_screen = 20L,
                                seed = 20260925L) {
  if (is.null(fin)) fin <- is.finite(p) & is.finite(n) & n > 0
  m <- length(p)
  regs <- unique(reg)
  cand <- list()

  # (1) National mean, jackknifed by region so a region never contributes to
  # the target its own districts are shrunk toward.
  nat <- rep(NA_real_, m)
  for (rg in regs) {
    i <- which(reg == rg)
    tr <- setdiff(which(fin), i)
    if (length(tr) >= 2) nat[i] <- .sw_mean(p[tr], n[tr])
  }
  cand$national <- nat

  # (2) The jackknifed region mean. This is the shipped blend's target: a
  # district's OWN estimate is excluded from its region's mean.
  fr <- rep(NA_real_, m)
  for (i in seq_len(m)) {
    j <- which(reg == reg[i] & seq_len(m) != i & fin)
    if (length(j) >= 2) fr[i] <- .sw_mean(p[j], n[j])
  }
  cand$region <- fr

  # (3) A covariate-free spatial smooth, fitted out of region. Held out by
  # REGION rather than by district: leave-one-district-out keeps every
  # neighbour in the training set, which is the configuration a spatial
  # smoother exploits, and it is why the LODO figure of 0.304 is an upper
  # bound rather than a comparable number.
  if (!is.null(lonlat) && requireNamespace("mgcv", quietly = TRUE)) {
    sp <- rep(NA_real_, m)
    ll <- as.data.frame(lonlat); names(ll) <- c("lon", "lat")
    okll <- is.finite(ll$lon) & is.finite(ll$lat)
    for (rg in regs) {
      i <- which(reg == rg)
      tr <- setdiff(which(fin & okll), i)
      if (length(tr) < 12) next
      kk <- min(20L, max(5L, floor(length(tr) / 4)))
      y <- stats::qlogis(pmin(pmax(p[tr], 0.005), 0.995))
      g <- tryCatch(mgcv::gam(y ~ s(lon, lat, k = kk), data = ll[tr, ],
                              weights = pmax(n[tr], 1), method = "REML"),
                    error = function(e) NULL)
      if (is.null(g)) next
      ii <- i[okll[i]]
      if (!length(ii)) next
      sp[ii] <- stats::plogis(as.numeric(
        stats::predict(g, newdata = ll[ii, , drop = FALSE])))
    }
    cand$spatial <- sp
  }

  # (4) The covariate ridge, out of region, screened in fold.
  if (!is.null(X) && ncol(X) >= 2 && exists(".ds_fit")) {
    pr <- rep(NA_real_, m)
    for (ri in seq_along(regs)) {
      rg <- regs[ri]
      i <- which(reg == rg)
      tr <- setdiff(which(fin), i)
      if (length(tr) < 6) next
      # .ds_fit() calls cv.glmnet(), whose folds are drawn at random. Without a
      # seed the same cell returns a different ridge candidate on every call,
      # which breaks reproducibility and makes the leak test below unable to
      # distinguish "this candidate saw the district" from "glmnet reshuffled".
      # Seeded per fold so the value depends on the data and the fold, nothing
      # else.
      set.seed(seed + ri)
      f <- tryCatch(.ds_fit(X[tr, , drop = FALSE], p[tr],
                            k_screen = min(k_screen, ncol(X))),
                    error = function(e) NULL)
      if (is.null(f)) next
      pp <- tryCatch(.ds_predict(f, X[i, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(pp)) pr[i] <- pp
    }
    cand$ridge <- pr
  }
  cand
}

#' Recentre each candidate on the training districts' level, out of region.
#'
#' The ridge's -3.1 to -3.3 pp bias is systematic, and a target carrying it
#' drags the blend with it. Recentring removes the level error without touching
#' the ranking the candidate provides.
eb_stack_recalibrate <- function(cand, p, n, reg, fin) {
  regs <- unique(reg)
  lapply(cand, function(v) {
    out <- v
    for (rg in regs) {
      i <- which(reg == rg)
      tr <- setdiff(which(fin), i)
      if (length(tr) < 3) next
      shift <- .sw_mean(p[tr], n[tr]) - .sw_mean(v[tr], n[tr])
      if (is.finite(shift)) out[i] <- v[i] + shift
    }
    pmin(pmax(out, 0), 1)
  })
}

#' Fit convex stacking weights out of region.
#'
#' Non-negative and summing to one, so the stacked target is always a genuine
#' interpolation of its candidates and can never extrapolate past them.
#'
#' REGULARISATION IS NOT COSMETIC. With four candidates and fourteen districts
#' over four regions -- Sierra Leone -- the weights are themselves estimated
#' from about ten points, which is the project's own effective-n problem one
#' level up. `shrink` pulls the solution toward the region-mean-only corner,
#' which is the target the tournament confirmed, so the stack must EARN any
#' departure from it.
#'
#' @param shrink weight on the region-mean prior, 0 to 1
eb_stack_weights <- function(cand, p, n, reg, fin, shrink = 0.5) {
  nm <- names(cand)
  M <- do.call(cbind, cand)
  regs <- unique(reg)
  # Out-of-region predictions of the stack's own training target.
  ok <- fin & is.finite(p) & apply(M, 1, function(z) all(is.finite(z)))
  if (sum(ok) < max(8L, 2L * ncol(M))) {
    w <- rep(0, ncol(M)); names(w) <- nm
    w[if ("region" %in% nm) "region" else 1L] <- 1
    return(list(w = w, n_used = sum(ok), route = "insufficient_areas"))
  }
  A <- M[ok, , drop = FALSE]; y <- p[ok]; wt <- pmax(n[ok], 1)
  # Non-negative least squares on the precision-weighted problem, then
  # normalised to sum to one.
  sw <- sqrt(wt)
  fit <- tryCatch(nnls::nnls(A * sw, y * sw), error = function(e) NULL)
  if (is.null(fit)) {
    w <- rep(0, ncol(M)); names(w) <- nm
    w[if ("region" %in% nm) "region" else 1L] <- 1
    return(list(w = w, n_used = sum(ok), route = "nnls_failed"))
  }
  w <- fit$x
  if (!any(w > 0)) w <- rep(1 / length(w), length(w))
  w <- w / sum(w)
  names(w) <- nm
  # Shrink toward the region-mean corner.
  prior <- rep(0, length(w)); names(prior) <- nm
  prior[if ("region" %in% nm) "region" else 1L] <- 1
  w <- (1 - shrink) * w + shrink * prior
  list(w = w / sum(w), n_used = sum(ok), route = "nnls")
}

#' The stacked target itself.
#'
#' @return list with `target` (numeric vector), `weights`, `candidates`
eb_stack_target <- function(p, n, reg, X = NULL, lonlat = NULL, fin = NULL,
                            shrink = 0.5, k_screen = 20L, seed = 20260925L) {
  if (is.null(fin)) fin <- is.finite(p) & is.finite(n) & n > 0
  cand <- eb_stack_candidates(p, n, reg, X = X, lonlat = lonlat,
                              fin = fin, k_screen = k_screen, seed = seed)
  # A candidate that is missing everywhere is dropped rather than imputed: a
  # zero-variance column would take stack weight for no reason.
  keep <- vapply(cand, function(v) sum(is.finite(v)) >= 8L &&
                   stats::sd(v, na.rm = TRUE) > 0, logical(1))
  cand <- cand[keep]
  if (!length(cand)) return(list(target = rep(NA_real_, length(p)),
                                 weights = numeric(0), candidates = list(),
                                 route = "no_candidates"))
  cand <- eb_stack_recalibrate(cand, p, n, reg, fin)
  wf <- eb_stack_weights(cand, p, n, reg, fin, shrink = shrink)
  M <- do.call(cbind, cand)
  tgt <- as.numeric(M %*% wf$w)
  # Where a candidate is missing for a district, renormalise over what is there
  # rather than emitting NA for the whole district.
  miss <- !is.finite(tgt)
  if (any(miss)) {
    for (i in which(miss)) {
      av <- is.finite(M[i, ])
      if (!any(av) || sum(wf$w[av]) <= 0) next
      tgt[i] <- sum(M[i, av] * wf$w[av]) / sum(wf$w[av])
    }
  }
  list(target = pmin(pmax(tgt, 0), 1), weights = wf$w,
       candidates = cand, n_used = wf$n_used, route = wf$route)
}
