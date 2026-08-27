# =============================================================================
# sandbox_parsimony/R/05a_loco_fns.R -- LOCO fitting/metric helpers
# Shared by 05_loco.R and 06_pool_composition.R.
# =============================================================================

center_by <- function(X, grp) {
  for (g in unique(grp)) {
    idx <- grp == g
    X[idx, ] <- sweep(X[idx, , drop = FALSE], 2,
                      colMeans(X[idx, , drop = FALSE]), "-")
  }
  X
}

# One LOCO fit.
#  centering "none" | "train" (as production codes it) | "own"
#  anchor    "train_mean" (production) | "oracle_national"
# oracle_national supplies the held-out country own national prevalence as the
# level. That is realistic in practice -- a national prevalence is almost always
# known from the national survey report even when district estimates are not --
# and it isolates how much LOCO error is level versus spatial pattern.
loco_one <- function(tr, te, vars, feature_fn = NULL, alpha = 0.5,
                     centering = "none", anchor = "train_mean",
                     weight = TRUE, seed = 12345L) {
  pp <- prep_X(tr, te, vars); Xtr <- pp$Xtr; Xte <- pp$Xte; vv <- pp$vars
  ytr <- .logit(tr$svy_prev)
  ytr_fit <- ytr; level <- mean(ytr)

  if (centering == "train") {
    tr_means <- colMeans(Xtr)
    Xtr <- center_by(Xtr, tr$country)
    Xte <- sweep(Xte, 2, tr_means, "-")
    for (g in unique(tr$country)) {
      i <- tr$country == g; ytr_fit[i] <- ytr[i] - mean(ytr[i])
    }
  } else if (centering == "own") {
    Xtr <- center_by(Xtr, tr$country)
    Xte <- sweep(Xte, 2, colMeans(Xte), "-")   # held-out COVARIATE means only
    for (g in unique(tr$country)) {
      i <- tr$country == g; ytr_fit[i] <- ytr[i] - mean(ytr[i])
    }
  }

  if (!is.null(feature_fn) && length(vv) > 1) {
    sel <- tryCatch(feature_fn(Xtr, ytr_fit, vv), error = function(e) vv)
    if (length(sel) >= 1) { Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE] }
  }
  w <- if (weight) pmax(tr$n_svy, 1) else NULL
  if (ncol(Xtr) < 2) return(rep(mean(tr$svy_prev), nrow(te)))
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, ytr_fit, alpha = alpha, weights = w, nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(tr$svy_prev), nrow(te)))
  plog <- as.numeric(stats::predict(cv, newx = Xte, s = "lambda.min"))

  if (centering %in% c("train", "own")) {
    lev <- if (anchor == "oracle_national")
      .logit(stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1))) else level
    plog <- plog + lev
  } else if (anchor == "oracle_national") {
    plog <- plog - mean(plog) +
      .logit(stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1)))
  }
  pmin(pmax(.ilogit(plog), 0), 1)
}

loco_metrics <- function(pred, te, deff = 1.5) {
  ok <- is.finite(pred) & is.finite(te$svy_prev)
  if (sum(ok) < 4) return(NULL)
  p <- pred[ok]; o <- te$svy_prev[ok]
  rel <- reliability(o, te$n_svy[ok], deff)
  dp <- p - mean(p); do <- o - mean(o)
  data.frame(
    n_test = sum(ok),
    r_max = round(rel$r_max, 3),
    pearson = round(suppressWarnings(stats::cor(p, o)), 3),
    spearman = round(suppressWarnings(stats::cor(p, o, method = "spearman")), 3),
    rmse_pp = round(sqrt(mean((p - o)^2)) * 100, 2),
    # RMSE with the level removed from BOTH sides: the error that remains once
    # the national offset is set aside. (Pearson r is already level-invariant,
    # so there is no separate "pattern r".)
    pattern_rmse_pp = round(sqrt(mean((dp - do)^2)) * 100, 2),
    level_bias_pp = round((mean(p) - mean(o)) * 100, 2),
    calib = round(tryCatch(unname(stats::coef(stats::lm(o ~ p))[2]),
                           error = function(e) NA_real_), 2),
    stringsAsFactors = FALSE)
}
