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
#   eb_blend      empirical Bayes between the two, lambda_d = tau2/(tau2 + v_d)
#   ridge         the covariate ridge under leave-one-region-out
#   region_tilt   the jackknifed region mean plus a covariate fit to
#                 within-region residuals
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
                           R = 200L, seed = 20260922L, k_screen = 20L) {
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
      # empirical Bayes toward the jackknifed region mean
      v_d <- area_sampling_var_shrunk(p_obs, pmax(n_obs, 1), deff = 1.5)
      tau2 <- max(stats::var(p_obs[fin]) - mean(v_d[fin]), 1e-8)
      lam <- tau2 / (tau2 + v_d)
      lam[!is.finite(lam)] <- 0
      est$eb_blend <- ifelse(is.finite(fr), lam * p_obs + (1 - lam) * fr, p_obs)
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

      for (nm in names(est)) {
        p <- est[[nm]]; f2 <- fin & is.finite(p)
        if (sum(f2) < 10 || stats::sd(truth[f2]) == 0) next
        acc[[nm]] <- rbind(acc[[nm]], data.frame(
          r_truth = suppressWarnings(stats::cor(truth[f2], p[f2])),
          mae_truth = 100 * mean(abs(p[f2] - truth[f2])),
          bias_truth = 100 * mean(p[f2] - truth[f2]),
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
        r_obs = round(mean(A$r_obs, na.rm = TRUE), 4),
        mae_obs = round(mean(A$mae_obs, na.rm = TRUE), 3),
        median_lambda = NA_real_, icc = round(icc, 4),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}
