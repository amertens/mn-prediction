# =============================================================================
# R/downscale_admin1.R
#
# Super-resolution downscaling: train at ADMIN-1, predict at ADMIN-2.
#
# THE CASE FOR IT, IN THIS PROJECT'S OWN NUMBERS
# ----------------------------------------------
# Measured across the 24 country x outcome cells with the project's noise model
# (area_sampling_var, deff = 1.5):
#
#   noise ceiling r_max      admin-2  0.132      admin-1  0.664
#   median unit sample size  admin-2  6-54       admin-1  22-214
#
# and after subtracting sampling noise, essentially ALL of the detectable
# between-district variance lies BETWEEN admin-1 units: in 22 of 24 cells the
# observed within-region variance is fully accounted for by sampling error. The
# ordering holds across the whole deff 1.0-2.0 band (5/24, 2/24 and 0/24 cells
# with any detectable within-region signal respectively).
#
# So fitting at admin-2 spends most of its degrees of freedom on noise, while
# fitting at admin-1 gives up almost nothing that these surveys can resolve and
# gets a target that is roughly five times more reliable.
#
# WHAT THIS CANNOT DO -- READ BEFORE QUOTING THE OUTPUT
# -----------------------------------------------------
# "No DETECTABLE within-region signal" is not "no within-region signal". The
# surveys cannot resolve it. Admin-2 predictions from an admin-1-trained model
# are therefore driven entirely by the covariates within a region, and they
# CANNOT BE VALIDATED against these surveys at that resolution. Reporting them
# as validated district estimates would be a claim the data does not support.
# Report: admin-1 performance (leave-one-region-out, where the ceiling is high
# enough for the number to mean something), plus admin-2 agreement with the
# direct estimates shown against r_max, plus the explicit statement that
# sub-regional detail is covariate-driven extrapolation.
#
# HOW THIS DIFFERS FROM THE PAPER THAT MOTIVATED IT
# -------------------------------------------------
# The Planetary Prediction Engine trained ADM1 -> ADM2 on 30 states x 40 months
# = 1,200 rows. We have one cross-section: 6 (Gambia) + 16 (Ghana) + 4 (Sierra
# Leone) + 27 (Malawi) = 53 pooled admin-1 units. Per-country training is
# therefore not viable and is refused below; the arm is POOLED across countries
# by construction, which also makes it a transportability estimator rather than
# a within-country one.
#
# With 53 rows and ~208 candidate covariates the p >> n regime is unchanged, so
# the estimator is a prescreen-then-regularise ridge rather than the full
# SuperLearner. That is a deliberate match to the sample size, not a shortcut.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

#' Aggregate an admin-2 covariate table up to admin-1.
#'
#' @param area_cov data.frame with Admin1, Admin2 and numeric covariates
#' @param weights optional data.frame(Admin2, w); population weights if
#'   available, otherwise districts are averaged equally (recorded in the
#'   returned attribute so the choice is never invisible)
#' @return data.frame with one row per Admin1
aggregate_covariates_to_admin1 <- function(area_cov, weights = NULL) {
  stopifnot(all(c("Admin1", "Admin2") %in% names(area_cov)))
  covs <- setdiff(names(area_cov), c("Admin1", "Admin2"))
  covs <- covs[vapply(covs, function(v) is.numeric(area_cov[[v]]), logical(1))]
  w <- rep(1, nrow(area_cov)); wsrc <- "equal (per district)"
  if (!is.null(weights) && all(c("Admin2", "w") %in% names(weights))) {
    ww <- weights$w[match(area_cov$Admin2, weights$Admin2)]
    if (any(is.finite(ww) & ww > 0)) {
      ww[!is.finite(ww) | ww <= 0] <- stats::median(ww[is.finite(ww) & ww > 0])
      w <- ww; wsrc <- "population"
    }
  }
  g <- trimws(as.character(area_cov$Admin1))
  regions <- sort(unique(g[!is.na(g)]))
  out <- data.frame(Admin1 = regions, stringsAsFactors = FALSE)
  for (v in covs) {
    x <- area_cov[[v]]
    out[[v]] <- vapply(regions, function(r) {
      i <- which(g == r & is.finite(x))
      if (!length(i)) NA_real_ else stats::weighted.mean(x[i], w[i])
    }, numeric(1))
  }
  attr(out, "weighting") <- wsrc
  out
}

#' Prescreen by absolute correlation, then fit a ridge.
#'
#' Screening happens INSIDE whatever fold the caller passes, never on the full
#' data -- selecting features on all rows and then cross-validating them is the
#' optimism the manuscript already documents for the production recipe.
.ds_fit <- function(X, y, k_screen = 20L, alpha = 0) {
  keep <- apply(X, 2, function(c) { v <- stats::var(c, na.rm = TRUE)
                                    is.finite(v) && v > 0 })
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(NULL)
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  r <- abs(suppressWarnings(stats::cor(X, y)))
  r[!is.finite(r)] <- 0
  sel <- colnames(X)[order(r, decreasing = TRUE)][seq_len(min(k_screen, ncol(X)))]
  Xs <- scale(X[, sel, drop = FALSE])
  ctr <- attr(Xs, "scaled:center"); scl <- attr(Xs, "scaled:scale")
  scl[!is.finite(scl) | scl == 0] <- 1
  if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
  ylg <- stats::qlogis(pmin(pmax(y, 0.005), 0.995))
  fit <- tryCatch(
    glmnet::cv.glmnet(Xs, ylg, alpha = alpha,
                      nfolds = max(3L, min(10L, floor(length(ylg) / 3)))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, sel = sel, center = ctr, scale = scl)
}

.ds_predict <- function(m, X) {
  if (is.null(m)) return(NULL)
  Z <- matrix(NA_real_, nrow(X), length(m$sel), dimnames = list(NULL, m$sel))
  for (j in seq_along(m$sel)) {
    v <- m$sel[j]
    x <- if (v %in% colnames(X)) X[, v] else rep(NA_real_, nrow(X))
    x[!is.finite(x)] <- m$center[j]
    Z[, j] <- (x - m$center[j]) / m$scale[j]
  }
  as.numeric(stats::plogis(stats::predict(m$fit, newx = Z, s = "lambda.min")))
}

#' Fit admin-1, predict admin-2.
#'
#' @param a1 data.frame: Admin1, country, prev (admin-1 direct estimate),
#'   n (its sample size) and the aggregated covariates
#' @param a2 data.frame: Admin1, Admin2, country and the admin-2 covariates
#' @param covs covariate names to use (must exist in both)
#' @param k_screen number of covariates kept by the in-fold prescreen
#' @param min_units refuse to fit below this many training regions
#' @return list(pred_admin2, cv_admin1, model)
fit_downscale_admin1 <- function(a1, a2, covs = NULL, k_screen = 20L,
                                 min_units = 30L) {
  covs <- covs %||% intersect(setdiff(names(a1), c("Admin1", "country", "prev", "n")),
                              setdiff(names(a2), c("Admin1", "Admin2", "country")))
  a1 <- a1[is.finite(a1$prev), , drop = FALSE]
  if (nrow(a1) < min_units) {
    cat(sprintf(paste0("[downscale] only %d admin-1 units (need >= %d). ",
                       "Per-country training is not viable at this sample size; ",
                       "pool countries or lower min_units deliberately.\n"),
                nrow(a1), min_units))
    return(NULL)
  }
  X1 <- as.matrix(a1[, covs, drop = FALSE]); y1 <- a1$prev

  # Leave-one-region-out CV, screening re-run inside every fold.
  oof <- rep(NA_real_, nrow(a1))
  for (i in seq_len(nrow(a1))) {
    m <- .ds_fit(X1[-i, , drop = FALSE], y1[-i], k_screen = k_screen)
    p <- .ds_predict(m, X1[i, , drop = FALSE])
    if (!is.null(p)) oof[i] <- p
  }
  ok <- is.finite(oof) & is.finite(y1)
  cv <- data.frame(
    n_units = sum(ok),
    pearson_r = if (sum(ok) > 3) suppressWarnings(stats::cor(y1[ok], oof[ok])) else NA_real_,
    spearman_r = if (sum(ok) > 3) suppressWarnings(stats::cor(y1[ok], oof[ok], method = "spearman")) else NA_real_,
    mae_pp = if (sum(ok) > 0) 100 * mean(abs(y1[ok] - oof[ok])) else NA_real_,
    rmse_pp = if (sum(ok) > 0) 100 * sqrt(mean((y1[ok] - oof[ok])^2)) else NA_real_)

  full <- .ds_fit(X1, y1, k_screen = k_screen)
  pred2 <- .ds_predict(full, as.matrix(a2[, covs, drop = FALSE]))
  out2 <- a2[, intersect(c("country", "Admin1", "Admin2"), names(a2)), drop = FALSE]
  out2$downscale_prev <- pred2

  cat(sprintf("[downscale] trained on %d admin-1 units, %d covariates screened to %d | LORO r = %.3f, MAE = %.2f pp | predicted %d admin-2 areas\n",
              nrow(a1), length(covs), length(full$sel %||% character()),
              cv$pearson_r, cv$mae_pp, sum(is.finite(pred2))))
  list(pred_admin2 = out2, cv_admin1 = cv, model = full,
       oof_admin1 = data.frame(a1[, intersect(c("country", "Admin1"), names(a1))],
                               prev = y1, oof = oof))
}
