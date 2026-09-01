# =============================================================================
# R/estimator_tournament.R
#
# WS-B2. Score every candidate district estimator against a KNOWN truth.
#
# WHY AGAINST SIMULATED TRUTH AND NOT AGAINST THE SURVEY
# ------------------------------------------------------
# Every estimator comparison in this project so far has scored predictions
# against the survey's district estimates, which carry sampling error of the
# same order as the differences between estimators. WS1 measured the empirical
# reliability at a median of 0.613, so roughly forty percent of the variance in
# the yardstick is noise. An estimator that shrinks toward a regional mean is
# penalised by that noise in a way an unshrunk estimator is not, which is
# precisely backwards.
#
# Here the truth is generated, so it is observable, and every estimator is
# scored against it.
#
# THE COVARIATE-LINKED TRUTH, AND WHY IT IS NECESSARY
# ---------------------------------------------------
# The WS1c simulator drew district truth independently of the covariates. Under
# that generator a covariate model is scored against noise by construction and
# loses trivially, which would make the tournament uninformative for exactly the
# estimators the shipping decision is about. The truth here is
#
#   logit(pi_d) = mu + sqrt(rho) * z_d + sqrt(1 - rho) * e_d
#
# where z_d is the standardised first principal component of the real district
# covariates and e_d is independent noise, both scaled to a chosen
# between-district standard deviation. rho is the share of true between-district
# variance that the covariates can in principle explain, and it is SWEPT rather
# than assumed, because the shipping decision depends on it and the measured
# value is uncertain. The measured anchor is roughly 0.33: WS4b puts achieved
# district correlation at a median of 0.200 against an empirical ceiling of
# 0.613.
#
# WHAT IS SCORED
#   direct        the survey's own district estimate
#   flat_region   the region mean, jackknifed so a district never sees itself
#   eb_blend      empirical Bayes between the two, lambda_d = tau2/(tau2 + v_d),
#                 with tau2 the RESIDUAL variance after the region mean
#   eb_blend_totvar  the same, but with tau2 from TOTAL between-district
#                 variance. This is the shipped estimator's moment fallback,
#                 kept scored so the cost of that route is measured.
#   ridge         the covariate ridge under leave-one-region-out
#   region_tilt   the jackknifed region mean plus a covariate fit to
#                 within-region residuals
#   eb_covariate  the Fay-Herriot form: shrink the district's own estimate
#                 toward the COVARIATE prediction rather than the region mean,
#                 with tau2 the residual between-district variance after the
#                 covariate fit. The first tournament shrank toward the weaker
#                 of the two available targets: the ridge beat the flat regional
#                 mean as a predictor (0.324 against 0.143 at rho 0.35), so the
#                 better shrinkage target had not been tried.
# Each is scored against TRUTH and, for reference, against the observed estimate.
# =============================================================================

#' One replicate of a covariate-linked simulated world.
#' @keywords internal
.sim_world <- function(a2, cl, w, z, lev, sd_logit, rho, icc, base_p) {
  e <- as.numeric(scale(stats::rnorm(length(lev))))
  lin <- sqrt(rho) * z + sqrt(max(1 - rho, 0)) * e
  mu <- stats::qlogis(min(max(base_p, 0.01), 0.99))
  pi_d <- stats::plogis(mu + sd_logit * as.numeric(scale(lin)))
  names(pi_d) <- lev
  pbar <- pi_d[a2]
  sd_c <- sqrt(pmax(icc, 0) * pbar * (1 - pbar))
  cu <- unique(paste(a2, cl, sep = "|::|"))
  eps <- stats::rnorm(length(cu), 0, 1); names(eps) <- cu
  p_i <- pmin(pmax(pbar + sd_c * eps[paste(a2, cl, sep = "|::|")], 0.001), 0.999)
  list(y = stats::rbinom(length(p_i), 1L, p_i), pi_d = pi_d)
}

#' Run the estimator tournament for one cell.
#'
#' @param d individual data supplying the real survey structure
#' @param X_d district covariate matrix, rows aligned to sorted district labels
#' @param region district-to-region map, aligned to the same labels
#' @param rho grid of covariate-signal shares
#' @param R replicates per setting
#' @return long data.frame: one row per estimator per setting
run_tournament <- function(d, a2_col, cl_col, w_col, y_col, X_d, region,
                           rho = c(0, 0.2, 0.35, 0.6), sd_logit = 0.8,
                           R = 200L, seed = 20260922L, k_screen = 20L,
                           lonlat = NULL, stack_shrink = 0.5) {
  y0 <- suppressWarnings(as.numeric(haven::zap_labels(d[[y_col]])))
  a2 <- as.character(d[[a2_col]]); cl <- as.character(d[[cl_col]])
  w  <- suppressWarnings(as.numeric(haven::zap_labels(d[[w_col]])))
  w[!is.finite(w) | w <= 0] <- 1
  ok <- is.finite(y0) & !is.na(a2) & !is.na(cl)
  y0 <- y0[ok]; a2 <- a2[ok]; cl <- cl[ok]; w <- w[ok]
  lev <- sort(unique(a2))
  keep <- lev %in% rownames(X_d)
  lev <- lev[keep]
  sel <- a2 %in% lev
  y0 <- y0[sel]; a2 <- a2[sel]; cl <- cl[sel]; w <- w[sel]
  if (length(lev) < 10) return(NULL)
  X_d <- X_d[lev, , drop = FALSE]
  reg <- region[lev]
  ic <- estimate_icc(y0, cl, a2); icc <- if (is.finite(ic$icc)) ic$icc else 0
  base_p <- mean(y0)

  # The covariate axis is a property of the covariates, not of a replicate, so
  # it is computed ONCE. The first version recomputed prcomp inside the replicate
  # loop, which cost a principal-component decomposition per draw and also let
  # the truth axis drift between replicates for no reason.
  z_axis <- if (ncol(X_d) < 2) rep(0, length(lev)) else {
    pc <- tryCatch(stats::prcomp(X_d, center = TRUE, scale. = TRUE)$x[, 1],
                   error = function(e) rep(0, nrow(X_d)))
    as.numeric(scale(pc))
  }
  z_axis[!is.finite(z_axis)] <- 0

  set.seed(seed)
  out <- list()
  for (rh in rho) {
    acc <- list()
    for (r in seq_len(R)) {
      sm <- .sim_world(a2, cl, w, z_axis, lev, sd_logit, rh, icc, base_p)
      ys <- sm$y; truth <- sm$pi_d
      g <- factor(a2, levels = lev)
      p_obs <- .wprev(ys, w, g, lev)
      n_obs <- as.integer(table(g))
      fin <- is.finite(p_obs) & n_obs >= 2
      if (sum(fin) < 10) next

      # direct
      est <- list(direct = p_obs)
      # jackknifed region mean
      fr <- rep(NA_real_, length(lev))
      for (i in seq_along(lev)) {
        j <- which(reg == reg[i] & seq_along(lev) != i & fin)
        if (length(j) >= 2) fr[i] <- stats::weighted.mean(p_obs[j], pmax(n_obs[j], 1))
      }
      est$flat_region <- fr
      # empirical Bayes toward the jackknifed region mean.
      #
      # tau2 IS THE RESIDUAL VARIANCE AFTER THE TARGET, NOT THE TOTAL.
      # For a shrinkage p_d -> m_d the model is theta_d = m_d + u_d with
      # u_d ~ N(0, tau2), so tau2 = Var(p - m) - vbar. An earlier version used
      # Var(p) - vbar, the TOTAL between-district variance, which overstates
      # tau2 by the between-REGION component -- about 40 percent of district
      # variance by WS-E2's decomposition. That inflated lambda, made eb_blend
      # shrink less than its own model says it should, and, because the
      # eb_covariate and eb_stack arms below were already residual-based, it
      # confounded "better target" with "better tau2 route" in the one
      # comparison this tournament exists to make.
      v_d <- area_sampling_var_shrunk(p_obs, pmax(n_obs, 1), deff = 1.5)
      okb <- fin & is.finite(fr)
      tau2 <- if (sum(okb) >= 10)
        max(stats::var(p_obs[okb] - fr[okb]) - mean(v_d[okb]), 1e-8) else
        max(stats::var(p_obs[fin]) - mean(v_d[fin]), 1e-8)
      lam <- tau2 / (tau2 + v_d)
      lam[!is.finite(lam)] <- 0
      est$eb_blend <- ifelse(is.finite(fr), lam * p_obs + (1 - lam) * fr, p_obs)
      # The total-variance route, kept as a SCORED ARM rather than deleted:
      # it is what the shipped estimator's moment fallback computes, so its
      # cost needs to be measured rather than assumed.
      tau2_tot <- max(stats::var(p_obs[fin]) - mean(v_d[fin]), 1e-8)
      lam_tot <- tau2_tot / (tau2_tot + v_d)
      lam_tot[!is.finite(lam_tot)] <- 0
      est$eb_blend_totvar <- ifelse(is.finite(fr),
                                    lam_tot * p_obs + (1 - lam_tot) * fr, p_obs)
      # covariate ridge, leave-one-region-out
      pr <- rep(NA_real_, length(lev))
      for (rg in unique(reg)) {
        i <- which(reg == rg); tr <- setdiff(which(fin), i)
        if (length(tr) < 6) next
        f <- .ds_fit(X_d[tr, , drop = FALSE], p_obs[tr],
                     k_screen = min(k_screen, ncol(X_d)))
        pp <- .ds_predict(f, X_d[i, , drop = FALSE])
        if (!is.null(pp)) pr[i] <- pp
      }
      est$ridge <- pr
      # region mean plus covariate tilt on within-region residuals
      tilt <- rep(NA_real_, length(lev))
      resid <- p_obs - fr
      for (rg in unique(reg)) {
        i <- which(reg == rg); tr <- setdiff(which(fin & is.finite(resid)), i)
        if (length(tr) < 6) next
        f <- .ds_fit(X_d[tr, , drop = FALSE], resid[tr],
                     k_screen = min(k_screen, ncol(X_d)))
        pp <- .ds_predict(f, X_d[i, , drop = FALSE])
        if (!is.null(pp)) tilt[i] <- pp
      }
      est$region_tilt <- ifelse(is.finite(fr) & is.finite(tilt),
                                pmin(pmax(fr + tilt, 0), 1), fr)
      # Fay-Herriot: shrink toward the covariate prediction. tau2 is the
      # RESIDUAL between-district variance after the covariate fit, so a
      # covariate model that explains more leaves less to shrink and the
      # district's own estimate gets more weight, not less.
      okr <- fin & is.finite(pr)
      if (sum(okr) >= 10) {
        tau2_r <- max(stats::var(p_obs[okr] - pr[okr]) - mean(v_d[okr]), 1e-8)
        lam_r <- tau2_r / (tau2_r + v_d)
        lam_r[!is.finite(lam_r)] <- 1
        est$eb_covariate <- ifelse(is.finite(pr),
          pmin(pmax(lam_r * p_obs + (1 - lam_r) * pr, 0), 1), est$eb_blend)
      }
      # The stacked target: the same empirical Bayes blend, but shrinking
      # toward a convex combination of national mean, jackknifed region mean,
      # covariate-free spatial smooth and screened ridge, each built out of
      # region and each RECENTRED out of region. It nests eb_blend exactly at
      # weights (0,1,0,0), so it can only lose by paying for flexibility it
      # does not use -- which is the thing this arm exists to measure.
      if (exists("eb_stack_target")) {
        # `pr` is the ridge arm's own out-of-region prediction, computed
        # above. Reusing it keeps the stack's covariate candidate identical to
        # the arm it is scored against, and halves the glmnet cost.
        st <- tryCatch(eb_stack_target(p_obs, pmax(n_obs, 1), reg, X = X_d,
                                       lonlat = lonlat, fin = fin,
                                       shrink = stack_shrink, k_screen = k_screen,
                                       seed = seed + r, ridge_pred = pr),
                       error = function(e) NULL)
        if (!is.null(st) && any(is.finite(st$target))) {
          tg <- st$target
          okt <- fin & is.finite(tg)
          # tau2 is the residual variance after the stacked target, so a better
          # target leaves less to shrink and the district keeps more weight.
          tau2_s <- if (sum(okt) >= 10)
            max(stats::var(p_obs[okt] - tg[okt]) - mean(v_d[okt]), 1e-8) else tau2
          lam_s <- tau2_s / (tau2_s + v_d)
          lam_s[!is.finite(lam_s)] <- 1
          est$eb_stack <- ifelse(is.finite(tg),
            pmin(pmax(lam_s * p_obs + (1 - lam_s) * tg, 0), 1), est$eb_blend)
        }
      }

      for (nm in names(est)) {
        p <- est[[nm]]; f2 <- fin & is.finite(p)
        if (sum(f2) < 10 || stats::sd(truth[f2]) == 0) next
        # TARGETING LIFT AGAINST TRUTH.
        #
        # The tournament scored only error and correlation, and on that basis
        # eb_stack looked like the winner: best MAE across the whole measured
        # range of covariate signal. But it never won correlation, and the
        # project's endpoint is now RANKING areas for programme allocation, not
        # estimating each area's prevalence. An estimator that shrinks toward a
        # smoother target compresses the ranking, which is what eb_stack's
        # correlation drop was showing. Scoring lift here makes that visible in
        # the tournament rather than in an argument about it afterwards.
        #
        # lift = mean TRUE prevalence in the areas this estimator ranks in the
        # top fifth, over the mean across all areas. 1.0 is random allocation.
        k <- max(1L, round(0.20 * sum(f2)))
        sel <- order(p[f2], decreasing = TRUE)[seq_len(k)]
        tv <- truth[f2]
        ov <- mean(tv)
        lift <- if (is.finite(ov) && ov > 0) mean(tv[sel]) / ov else NA_real_
        # The ceiling: what a perfectly informed allocation would reach in this
        # cell. Lift is scale-free but not spread-free, so a bare 1.2 means
        # something different where the ceiling is 1.3 than where it is 5.
        osel <- order(tv, decreasing = TRUE)[seq_len(k)]
        lift_oracle <- if (is.finite(ov) && ov > 0) mean(tv[osel]) / ov else NA_real_
        acc[[nm]] <- rbind(acc[[nm]], data.frame(
          r_truth = suppressWarnings(stats::cor(truth[f2], p[f2])),
          mae_truth = 100 * mean(abs(p[f2] - truth[f2])),
          bias_truth = 100 * mean(p[f2] - truth[f2]),
          lift_truth = lift,
          lift_oracle = lift_oracle,
          lift_share = if (is.finite(lift_oracle) && lift_oracle > 1)
            (lift - 1) / (lift_oracle - 1) else NA_real_,
          r_obs = suppressWarnings(stats::cor(p_obs[f2], p[f2])),
          mae_obs = 100 * mean(abs(p[f2] - p_obs[f2]))))
      }
    }
    for (nm in names(acc)) {
      A <- acc[[nm]]
      out[[length(out) + 1L]] <- data.frame(
        estimator = nm, rho = rh, sd_logit = sd_logit, replicates = nrow(A),
        r_truth = round(mean(A$r_truth, na.rm = TRUE), 4),
        r_truth_lo = round(stats::quantile(A$r_truth, .1, na.rm = TRUE), 4),
        r_truth_hi = round(stats::quantile(A$r_truth, .9, na.rm = TRUE), 4),
        mae_truth = round(mean(A$mae_truth, na.rm = TRUE), 3),
        bias_truth = round(mean(A$bias_truth, na.rm = TRUE), 3),
        # Targeting lift, carried through the per-replicate accumulator into
        # the returned table. Adding it to acc[] alone was not enough: this
        # aggregation names its columns explicitly, so a new metric that is not
        # listed here is computed and then silently discarded.
        lift_truth = round(mean(A$lift_truth, na.rm = TRUE), 4),
        lift_oracle = round(mean(A$lift_oracle, na.rm = TRUE), 4),
        lift_share = round(mean(A$lift_share, na.rm = TRUE), 4),
        r_obs = round(mean(A$r_obs, na.rm = TRUE), 4),
        mae_obs = round(mean(A$mae_obs, na.rm = TRUE), 3),
        median_lambda = NA_real_, icc = round(icc, 4),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}
