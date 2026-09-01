# =============================================================================
# R/loco_outlyingness.R
#
# Is a held-out country's covariate profile INSIDE the training countries'
# distribution, or is the model extrapolating?
#
# WHY THIS MATTERS FOR THIS PROJECT
# ---------------------------------
# Leave-one-country-out reports a mean r of 0.184 and a level offset reaching
# -42 pp. Those numbers describe the failure; they do not say whether the model
# was interpolating badly or being asked a question outside its training range.
# Those call for different responses. A model interpolating badly needs a better
# model. A model extrapolating needs to REFUSE, and the honest deliverable is a
# flag rather than a number.
#
# Ported after an independent implementation on another branch found its Sierra
# Leone LOCO predictions to be extrapolations rather than shrinkage. A cheap
# marginal range check here disagreed with that -- it flags Malawi at 100
# percent and Sierra Leone at 0 -- but a marginal check asks whether each
# covariate is separately in range, while the question is whether the JOINT
# point is plausible. A district can pass every marginal range and sit far
# outside the joint distribution. This file answers the joint question.
#
# WHY STAHEL-DONOHO, AND NOT THE TWO OBVIOUS ALTERNATIVES
# -------------------------------------------------------
# HALFSPACE DEPTH saturates. The other implementation verified it returns zero
# for essentially every point at p=16, n=180, against a Gaussian sanity check.
# Our setting is worse, so it is not attempted.
#
# MAHALANOBIS NEEDS AN INVERTIBLE COVARIANCE and we do not have one. The
# training set is 119 to 192 districts against up to 294 predictors, so p > n
# and the sample covariance is singular by construction. The other
# implementation reported condition numbers of 2e5 to 1.8e6 at p=16; ours is not
# merely ill-conditioned but rank-deficient.
#
# STAHEL-DONOHO BY RANDOM PROJECTION WAS TRIED FIRST, AND FAILS HERE. It needs
# no covariance inverse, which is the right instinct at p > n, but it has to be
# computed after a dimension reduction -- random directions in 294 dimensions
# where only ~150 are populated mostly probe noise. And the reduction destroys
# the thing being looked for. Variance-based truncation discards the LOWEST
# variance directions, which is exactly where a point off a correlation ridge
# sits. Measured on a planted joint outlier -- inside every marginal range, two
# standard deviations off a strong ridge -- the retained subspace scored it at
# 0.34 against a training maximum of 5.04, calling the single most anomalous
# point in the set the most typical one. Whitening the components did not fix
# it, because the informative direction had already been dropped.
#
# WHAT IS USED INSTEAD is the standard decomposition for this setting, from
# robust PCA: every point gets TWO distances, and they answer different
# questions.
#
#   SCORE DISTANCE (SD)       how far along the training subspace, robust
#                             Mahalanobis in the retained components. Catches a
#                             point that is an extreme version of the training
#                             pattern.
#
#   ORTHOGONAL DISTANCE (OD)  how far OFF the training subspace. Catches a point
#                             that breaks the correlation structure -- the case
#                             the marginal range check provably cannot see, and
#                             the case the truncated projection method above
#                             threw away.
#
# A district is flagged when it exceeds the training cut-off on EITHER. Both are
# reported, because which one fires says what kind of extrapolation it is.
#
# The subspace, the centre, the scales and both cut-offs are fitted on TRAINING
# ONLY, and the held-out country is projected onto them. Fitting on pooled data
# would let the held-out country help define the space it is judged against,
# which is the leak class this project has found twelve other instances of.
# =============================================================================

#' Rank-deficiency-safe standardisation fitted on the training rows.
.lo_scaler <- function(Xtr) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  sdv <- apply(Xtr, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
  list(mu = mu, sd = sdv, keep = keep)
}
.lo_apply <- function(X, sc) {
  X <- X[, sc$keep, drop = FALSE]
  sweep(sweep(X, 2, sc$mu[sc$keep], "-"), 2, sc$sd[sc$keep], "/")
}

#' Extrapolation distances of test points against a training cloud.
#'
#' @param Xtr training matrix, rows are areas
#' @param Xte test matrix, same columns
#' @param max_pc cap on retained components
#' @param var_target cumulative variance the components should reach
#' @param quant training quantile defining the cut-off
#' @return list of distances, cut-offs and the diagnostics needed to judge
#'   whether the measure is trustworthy in this cell
loco_outlyingness <- function(Xtr, Xte, max_pc = 10L, var_target = 0.95,
                              quant = 0.975, seed = 20260927L, n_dir = NULL) {
  stopifnot(ncol(Xtr) == ncol(Xte))
  sc <- .lo_scaler(Xtr)
  if (sum(sc$keep) < 2) return(NULL)
  Ztr <- .lo_apply(Xtr, sc); Zte <- .lo_apply(Xte, sc)
  Ztr[!is.finite(Ztr)] <- 0; Zte[!is.finite(Zte)] <- 0

  n_tr <- nrow(Ztr)
  k_cap <- max(2L, min(max_pc, floor(n_tr / 5)))
  pc <- tryCatch(stats::prcomp(Ztr, center = FALSE, scale. = FALSE),
                 error = function(e) NULL)
  if (is.null(pc)) return(NULL)
  ev <- pc$sdev^2
  cum <- cumsum(ev) / sum(ev)
  k <- min(k_cap, max(2L, which(cum >= var_target)[1]))
  if (!is.finite(k)) k <- k_cap
  R <- pc$rotation[, seq_len(k), drop = FALSE]
  Ptr <- Ztr %*% R; Pte <- Zte %*% R

  # SCORE DISTANCE: robust Mahalanobis inside the retained subspace. The
  # components are orthogonal by construction, so this is a scaled Euclidean
  # norm and needs no matrix inverse. MAD rather than SD so a handful of
  # extreme training districts cannot inflate the scale and hide everything.
  s_tr <- apply(Ptr, 2, stats::mad)
  bad <- !is.finite(s_tr) | s_tr <= .Machine$double.eps^0.5
  if (any(bad)) s_tr[bad] <- apply(Ptr, 2, stats::sd)[bad]
  bad <- !is.finite(s_tr) | s_tr <= .Machine$double.eps^0.5
  if (any(bad)) s_tr[bad] <- 1
  ctr <- apply(Ptr, 2, stats::median)
  sd_tr <- sqrt(rowSums(sweep(sweep(Ptr, 2, ctr, "-"), 2, s_tr, "/")^2))
  sd_te <- sqrt(rowSums(sweep(sweep(Pte, 2, ctr, "-"), 2, s_tr, "/")^2))

  # ORTHOGONAL DISTANCE: the residual after projecting onto the subspace and
  # back. This is the component the truncated projection method discarded.
  od_tr <- sqrt(rowSums((Ztr - Ptr %*% t(R))^2))
  od_te <- sqrt(rowSums((Zte - Pte %*% t(R))^2))

  cut_sd <- as.numeric(stats::quantile(sd_tr, quant, na.rm = TRUE))
  # ORTHOGONAL DISTANCE NEEDS THE ROBPCA CUT-OFF, NOT A TRAINING QUANTILE.
  #
  # An empirical quantile of the training ODs saturates here and reports every
  # held-out district in every country as extrapolating -- measured at 1.000 for
  # all four. That is uninformative in the same way halfspace depth is, and for
  # a related reason: with ~10 components retained from ~200 covariates the
  # residual outside the subspace is large for everyone, so a cut-off set at the
  # training 97.5th percentile is cleared by any point from a country whose
  # covariance structure differs at all.
  #
  # The standard construction (Hubert et al.) uses the fact that OD^(2/3) is
  # approximately normal, and sets the cut-off from a ROBUST centre and scale of
  # that transform, so it scales with the subspace's own residual size instead
  # of with the training sample's upper tail.
  od23 <- od_tr^(2/3)
  m23 <- stats::median(od23); s23 <- stats::mad(od23)
  if (!is.finite(s23) || s23 <= 0) s23 <- stats::sd(od23)
  cut_od <- if (is.finite(m23) && is.finite(s23) && s23 > 0)
    (m23 + stats::qnorm(quant) * s23)^(3/2) else
    as.numeric(stats::quantile(od_tr, quant, na.rm = TRUE))
  cn <- tryCatch({ sv <- svd(Ptr)$d
    if (length(sv) && min(sv) > 0) max(sv) / min(sv) else Inf },
    error = function(e) NA_real_)

  list(sd_test = as.numeric(sd_te), sd_train = as.numeric(sd_tr),
       od_test = as.numeric(od_te), od_train = as.numeric(od_tr),
       out_test = as.numeric(pmax(sd_te / cut_sd, od_te / cut_od)),
       out_train = as.numeric(pmax(sd_tr / cut_sd, od_tr / cut_od)),
       cut_sd = cut_sd, cut_od = cut_od,
       flag_test = (sd_te > cut_sd) | (od_te > cut_od),
       n_pc = k, n_cov_used = sum(sc$keep),
       var_explained = as.numeric(cum[k]), cond_number = cn,
       thresh_p95 = as.numeric(stats::quantile(pmax(sd_tr / cut_sd,
                                                    od_tr / cut_od), 0.95, na.rm = TRUE)),
       thresh_max = max(pmax(sd_tr / cut_sd, od_tr / cut_od), na.rm = TRUE))
}

#' Summarise one held-out country's extrapolation.
#'
#' The two thresholds answer different questions and both are reported. Above
#' the training 95th percentile means "unusual for a training district". Above
#' the training MAXIMUM means "outside anything the model has seen", which is
#' the one that should gate a published number.
loco_outlyingness_summary <- function(lo, marginal_frac = NA_real_) {
  if (is.null(lo)) return(NULL)
  data.frame(
    n_test = length(lo$out_test),
    n_pc = lo$n_pc, n_cov = lo$n_cov_used,
    var_explained = round(lo$var_explained, 3),
    cond_number = signif(lo$cond_number, 3),
    sd_median = round(stats::median(lo$sd_test, na.rm = TRUE), 2),
    sd_cut = round(lo$cut_sd, 2),
    od_median = round(stats::median(lo$od_test, na.rm = TRUE), 2),
    od_cut = round(lo$cut_od, 2),
    frac_flagged = round(mean(lo$flag_test, na.rm = TRUE), 3),
    # Ratios, because the FRACTIONS above the orthogonal cut-off saturate at
    # 1.000 in every fold: the minimum held-out orthogonal distance exceeds the
    # training cut-off in all four countries. That is a real property of
    # cross-country transport in this covariate set, not a coding artefact --
    # the retained subspace holds about 69 percent of the variance and the
    # residual differs systematically by country -- but a fraction that is
    # always one cannot discriminate, so the magnitude is reported too.
    od_ratio = round(stats::median(lo$od_test, na.rm = TRUE) / lo$cut_od, 2),
    sd_ratio = round(stats::median(lo$sd_test, na.rm = TRUE) / lo$cut_sd, 2),
    frac_score_out = round(mean(lo$sd_test > lo$cut_sd, na.rm = TRUE), 3),
    frac_orth_out = round(mean(lo$od_test > lo$cut_od, na.rm = TRUE), 3),
    frac_above_train_max = round(mean(lo$out_test > lo$thresh_max, na.rm = TRUE), 3),
    frac_marginal_out = round(marginal_frac, 3),
    stringsAsFactors = FALSE)
}

#' The cheap marginal check, kept alongside deliberately.
#'
#' It is the weaker measure and disagreeing with the joint one is informative:
#' marginal-in-range but joint-outlying is exactly the case that motivates
#' computing the joint measure at all.
loco_marginal_out <- function(Xtr, Xte, tol = 0.10) {
  lo <- apply(Xtr, 2, min, na.rm = TRUE)
  hi <- apply(Xtr, 2, max, na.rm = TRUE)
  frac <- rowMeans(sweep(Xte, 2, lo, "<") | sweep(Xte, 2, hi, ">"), na.rm = TRUE)
  mean(frac > tol, na.rm = TRUE)
}
