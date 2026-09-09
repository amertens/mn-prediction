# =============================================================================
# R/protocol_v2_mbg.R   [MB-01]
#
# MODEL-BASED GEOSTATISTICS AS A COMPARATOR ARM ON THE CLUSTER TRACK
#
# The DHS Program's Admin-2 estimates are produced by a Bayesian geostatistical
# model fit at the survey cluster: a spatial Gaussian field plus covariates,
# predicted to a fine grid and aggregated to districts with population weights
# (DHS Spatial Analysis Reports 20 and 21). This file adds that estimator as an
# arm with the same (tr, te, y, X, D, aux) interface as every other arm, so it
# is scored under the identical district-held-out folds, the identical
# covariates and the identical aggregation to districts as the domain index.
#
#   mbg      SPDE Gaussian field + the first principal component of every
#            covariate domain as fixed effects (about 24 covariates; DHS uses a
#            comparable number), N(0, 1) priors on the standardised effects
#   mbg_gp   the spatial field alone (pure interpolation of the survey's own
#            clusters), the floor a geostatistical model has without covariates
#
# Likelihood: binomial on the cluster counts (deficient / assayed) for the
# prevalence target, Gaussian for the level target. Predictions are the
# posterior mean of the linear predictor at the held-out clusters, on the
# modelling scale the scorer expects (logit for prevalence). INLA with the
# empirical-Bayes integration strategy for speed (a fit is a few seconds at
# 60-110 clusters). The arm cannot run in the transport estimand (it needs the
# outcome at clusters inside the country), which is the point of the
# comparison.
#
# Requires: INLA, fmesher (installed 2026-09). Registered into ARMS_V2 when this
# file is sourced after R/protocol_v2.R (tar_source sorts by name, and
# "protocol_v2_mbg" sorts after "protocol_v2"). The cluster benchmark passes
# aux$target, aux$n_raw and aux$y_prev for the binomial counts.
# =============================================================================

.mbg_ok <- function() requireNamespace("INLA", quietly = TRUE) && requireNamespace("fmesher", quietly = TRUE)

#' One SPDE fit; returns predictions at `te` on the modelling scale, or NULL on failure.
.mbg_fit_predict <- function(tr, te, y, covs, aux, family = c("gaussian", "binomial"),
                             y_count = NULL, n_trials = NULL) {
  family <- match.arg(family)
  if (!.mbg_ok()) return(NULL)
  lon <- aux$lon; lat <- aux$lat
  ok_tr <- tr[is.finite(lon[tr]) & is.finite(lat[tr]) & is.finite(y[tr])]
  if (length(ok_tr) < 12) return(NULL)
  loc_all <- cbind(lon[c(ok_tr, te)], lat[c(ok_tr, te)])
  loc_all <- loc_all[is.finite(loc_all[, 1]) & is.finite(loc_all[, 2]), , drop = FALSE]
  # mesh in degrees: inner edge ~15 km, outer ~60 km, points closer than ~3 km merged
  mesh <- tryCatch(fmesher::fm_mesh_2d_inla(loc = loc_all, max.edge = c(0.15, 0.6), cutoff = 0.03, offset = c(0.2, 0.8)),
                   error = function(e) NULL)
  if (is.null(mesh)) return(NULL)
  # PC priors: P(range < 0.5 deg ~ 55 km) = 0.5, P(sd > 1) = 0.05
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.5, 0.5), prior.sigma = c(1, 0.05))
  A_tr <- INLA::inla.spde.make.A(mesh, loc = cbind(lon[ok_tr], lat[ok_tr]))
  A_te <- INLA::inla.spde.make.A(mesh, loc = cbind(lon[te], lat[te]))
  idx <- INLA::inla.spde.make.index("field", n.spde = spde$n.spde)
  p <- if (is.null(covs)) 0L else ncol(covs)
  eff_tr <- list(c(idx, list(Intercept = 1)))
  eff_te <- list(c(idx, list(Intercept = 1)))
  A_list_tr <- list(A_tr); A_list_te <- list(A_te)
  if (p > 0) {
    ctr <- as.data.frame(covs[ok_tr, , drop = FALSE]); cte <- as.data.frame(covs[te, , drop = FALSE])
    names(ctr) <- names(cte) <- paste0("v", seq_len(p))
    eff_tr <- list(c(idx, list(Intercept = 1)), ctr); A_list_tr <- list(A_tr, 1)
    eff_te <- list(c(idx, list(Intercept = 1)), cte); A_list_te <- list(A_te, 1)
  }
  ytr <- if (family == "binomial") y_count[ok_tr] else y[ok_tr]
  dat_tr <- if (family == "binomial") list(y = ytr, Ntrials = n_trials[ok_tr]) else list(y = ytr)
  dat_te <- if (family == "binomial") list(y = rep(NA, length(te)), Ntrials = rep(1, length(te))) else list(y = rep(NA, length(te)))
  stk_tr <- INLA::inla.stack(data = dat_tr, A = A_list_tr, effects = eff_tr, tag = "est")
  stk_te <- INLA::inla.stack(data = dat_te, A = A_list_te, effects = eff_te, tag = "pred")
  stk <- INLA::inla.stack(stk_tr, stk_te)
  rhs <- if (p > 0) paste(c("-1 + Intercept", paste0("v", seq_len(p)), "f(field, model = spde)"), collapse = " + ") else "-1 + Intercept + f(field, model = spde)"
  fml <- stats::as.formula(paste("y ~", rhs))
  fit <- tryCatch(suppressWarnings(INLA::inla(
    fml, family = family, data = INLA::inla.stack.data(stk),
    Ntrials = if (family == "binomial") INLA::inla.stack.data(stk)$Ntrials else NULL,
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE, link = 1),
    control.fixed = list(prec = 1, prec.intercept = 0.01),
    control.inla = list(int.strategy = "eb"),
    control.compute = list(config = FALSE, dic = FALSE, waic = FALSE),
    num.threads = "2:1", verbose = FALSE, silent = TRUE)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ip <- INLA::inla.stack.index(stk, "pred")$data
  pred <- fit$summary.linear.predictor[ip, "mean"]
  if (length(pred) != length(te) || !all(is.finite(pred))) return(NULL)
  as.numeric(pred)
}

.mbg_covs <- function(D, tr) {
  if (is.null(D) || !ncol(D)) return(NULL)
  pc1 <- grep("_PC1$", colnames(D), value = TRUE)
  cols <- if (length(pc1) >= 2) pc1 else colnames(D)
  M <- D[, cols, drop = FALSE]
  mu <- colMeans(M[tr, , drop = FALSE]); s <- apply(M[tr, , drop = FALSE], 2, stats::sd); s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(M, 2, mu, "-"), 2, s, "/")
}

.mbg_arm <- function(tr, te, y, X, D, aux, with_covs = TRUE) {
  fallback <- rep(mean(y[tr]), length(te))
  if (is.null(aux$lon) || is.null(aux$lat)) return(fallback)
  covs <- if (with_covs) .mbg_covs(D, tr) else NULL
  if (identical(aux$target, "prev") && !is.null(aux$n_raw) && !is.null(aux$y_prev)) {
    n <- pmax(1, round(aux$n_raw)); k <- pmin(n, pmax(0, round(aux$y_prev * n)))
    p <- .mbg_fit_predict(tr, te, y, covs, aux, family = "binomial", y_count = k, n_trials = n)
  } else {
    p <- .mbg_fit_predict(tr, te, y, covs, aux, family = "gaussian")
  }
  if (is.null(p)) fallback else p
}

# The DHS Program's own covariate shortlist, as far as the cluster-buffer set carries it
# (accessibility, aridity via PET and precipitation, elevation, EVI, night lights,
# population density, precipitation, land-surface and air temperature, and the one
# malaria-intervention layer extracted at the cluster; PfPR and ITN are not in the
# cluster set). Raw rank-normalised columns, not domain axes: the variant the
# collaborator asked for (MB-01, 2026-09-08).
MBG_DHS_SHORTLIST <- c("access_accessibility", "tclim_pet", "elev_elevation", "evi_evi", "ntl_viirs_dmsp2013v1_b1",
                       "popdens_unwpp_adjusted_population_density", "chirps_mean", "lst_night_mean", "tclim_tmmx",
                       "map_interventions_202106_africa_irs_coverage")
.mbg_covs_dhs <- function(X, tr) {
  cols <- intersect(MBG_DHS_SHORTLIST, colnames(X))
  if (length(cols) < 3) return(NULL)
  M <- X[, cols, drop = FALSE]
  mu <- colMeans(M[tr, , drop = FALSE]); s <- apply(M[tr, , drop = FALSE], 2, stats::sd); s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(M, 2, mu, "-"), 2, s, "/")
}
arm_mbg_dhs_v2 <- function(tr, te, y, X, D, aux) {
  fallback <- rep(mean(y[tr]), length(te))
  if (is.null(aux$lon) || is.null(aux$lat)) return(fallback)
  covs <- .mbg_covs_dhs(X, tr); if (is.null(covs)) return(fallback)
  if (identical(aux$target, "prev") && !is.null(aux$n_raw) && !is.null(aux$y_prev)) {
    n <- pmax(1, round(aux$n_raw)); k <- pmin(n, pmax(0, round(aux$y_prev * n)))
    p <- .mbg_fit_predict(tr, te, y, covs, aux, family = "binomial", y_count = k, n_trials = n)
  } else p <- .mbg_fit_predict(tr, te, y, covs, aux, family = "gaussian")
  if (is.null(p)) fallback else p
}

arm_mbg_v2    <- function(tr, te, y, X, D, aux) .mbg_arm(tr, te, y, X, D, aux, with_covs = TRUE)
arm_mbg_gp_v2 <- function(tr, te, y, X, D, aux) .mbg_arm(tr, te, y, X, D, aux, with_covs = FALSE)

if (exists("ARMS_V2") && is.list(ARMS_V2)) {
  ARMS_V2$mbg    <- arm_mbg_v2
  ARMS_V2$mbg_gp <- arm_mbg_gp_v2
  ARMS_V2$mbg_dhs <- arm_mbg_dhs_v2
}
