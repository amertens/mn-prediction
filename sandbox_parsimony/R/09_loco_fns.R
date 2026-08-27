# =============================================================================
# sandbox_parsimony/R/09_loco_fns.R -- spatial and z-score LOCO fitters
# Shared by 09_loco_spatial.R and 11_loco_headline.R.
# =============================================================================

#' TPS on lon/lat, optionally plus a centered covariate block
loco_spatial <- function(tr, te, vars, feature_fn = NULL, k_spline = 30L,
                         center_covars = TRUE, anchor = "train_mean") {
  if (!all(is.finite(tr$lon)) || !all(is.finite(te$lon)))
    return(rep(mean(tr$svy_prev), nrow(te)))
  Y  <- .logit(tr$svy_prev)
  wt <- pmax(tr$n_svy, 1)
  ftr <- data.frame(Y = Y, lon = tr$lon, lat = tr$lat)
  fte <- data.frame(lon = te$lon, lat = te$lat)
  rhs <- sprintf("s(lon, lat, k = %d, bs = 'tp')", min(k_spline, nrow(tr) - 5))

  if (!is.null(feature_fn)) {
    pp <- prep_X(tr, te, vars); Xtr <- pp$Xtr; Xte <- pp$Xte; vv <- pp$vars
    if (center_covars) {
      Xtr <- center_by(Xtr, tr$country)
      Xte <- sweep(Xte, 2, colMeans(Xte), "-")
    }
    sel <- tryCatch(feature_fn(Xtr, Y - ave(Y, tr$country), vv), error = function(e) vv)
    sel <- utils::head(sel, 12L)
    if (length(sel) >= 1) {
      Ztr <- as.data.frame(scale(Xtr[, sel, drop = FALSE]))
      ctr <- attr(scale(Xtr[, sel, drop = FALSE]), "scaled:center")
      scl <- attr(scale(Xtr[, sel, drop = FALSE]), "scaled:scale")
      scl[!is.finite(scl) | scl < 1e-8] <- 1
      Zte <- as.data.frame(sweep(sweep(Xte[, sel, drop = FALSE], 2, ctr, "-"), 2, scl, "/"))
      nm <- paste0("z", seq_along(sel)); names(Ztr) <- names(Zte) <- nm
      ftr <- cbind(ftr, Ztr); fte <- cbind(fte, Zte)
      rhs <- paste(c(rhs, nm), collapse = " + ")
    }
  }
  g <- tryCatch(mgcv::gam(stats::as.formula(paste("Y ~", rhs)), data = ftr,
                          weights = wt, method = "REML"), error = function(e) NULL)
  if (is.null(g)) return(rep(mean(tr$svy_prev), nrow(te)))
  pl <- as.numeric(stats::predict(g, newdata = fte))
  if (anchor == "oracle_national")
    pl <- pl - mean(pl) + .logit(stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1)))
  pmin(pmax(.ilogit(pl), 0), 1)
}

#' Train on the within-country z-score of logit prevalence.
#' The model is structurally incapable of learning a between-country level or
#' spread, so it can only transport spatial pattern -- which is what LOCO can
#' realistically deliver. The held-out level and spread come from the anchor.
loco_zscore <- function(tr, te, vars, feature_fn = NULL, alpha = 0.5,
                        anchor = "train_mean") {
  pp <- prep_X(tr, te, vars); Xtr <- pp$Xtr; Xte <- pp$Xte; vv <- pp$vars
  Xtr <- center_by(Xtr, tr$country)
  Xte <- sweep(Xte, 2, colMeans(Xte), "-")
  # scale covariates within country too, so a country with a wider covariate
  # range does not dominate the fit
  for (g in unique(tr$country)) {
    i <- tr$country == g
    s <- apply(Xtr[i, , drop = FALSE], 2, stats::sd); s[!is.finite(s) | s < 1e-8] <- 1
    Xtr[i, ] <- sweep(Xtr[i, , drop = FALSE], 2, s, "/")
  }
  s_te <- apply(Xte, 2, stats::sd); s_te[!is.finite(s_te) | s_te < 1e-8] <- 1
  Xte <- sweep(Xte, 2, s_te, "/")

  y <- .logit(tr$svy_prev); yz <- y
  for (g in unique(tr$country)) {
    i <- tr$country == g
    s <- stats::sd(y[i]); if (!is.finite(s) || s < 1e-8) s <- 1
    yz[i] <- (y[i] - mean(y[i])) / s
  }
  if (!is.null(feature_fn) && length(vv) > 1) {
    sel <- tryCatch(feature_fn(Xtr, yz, vv), error = function(e) vv)
    if (length(sel) >= 1) { Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE] }
  }
  if (ncol(Xtr) < 2) return(rep(mean(tr$svy_prev), nrow(te)))
  set.seed(12345L)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, yz, alpha = alpha,
                                   weights = pmax(tr$n_svy, 1), nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(tr$svy_prev), nrow(te)))
  zhat <- as.numeric(stats::predict(cv, newx = Xte, s = "lambda.min"))

  # Re-scale the predicted z-pattern onto the held-out country logit scale.
  # Level and spread both come from the anchor; the model supplies only shape.
  if (anchor == "oracle_national") {
    lev <- .logit(stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1)))
    spread <- stats::sd(.logit(te$svy_prev))
  } else {
    lev <- mean(y)
    spread <- mean(vapply(unique(tr$country),
                          function(g) stats::sd(y[tr$country == g]), numeric(1)),
                   na.rm = TRUE)
  }
  if (!is.finite(spread) || spread < 1e-8) spread <- 1
  sz <- stats::sd(zhat); if (!is.finite(sz) || sz < 1e-8) sz <- 1
  pmin(pmax(.ilogit(lev + spread * (zhat - mean(zhat)) / sz), 0), 1)
}

