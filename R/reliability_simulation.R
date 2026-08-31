# =============================================================================
# R/reliability_simulation.R
#
# IS THE CEILING ESTIMATOR BIASED? SIMULATE A WORLD WHERE THE ANSWER IS KNOWN.
#
# WHY THIS EXISTS
# ---------------
# WS1a measures reliability by splitting the real sample in half. WS1b compares
# that measurement against the analytic ceiling. Neither can say which of the two
# is closer to the truth, because on real data the truth is not observed.
#
# Here it is. District prevalences are chosen by the simulation, so the quantity
# every estimator is trying to recover is known exactly. The survey structure is
# taken from the real data (actual districts, actual clusters, actual number of
# respondents per cluster) and only the outcome is simulated, so the sampling
# geometry that drives the whole problem is preserved.
#
# WHAT IS BEING ESTIMATED
# -----------------------
# r_max is meant to bound the correlation any predictor can achieve against a
# NOISY district estimate. The quantity it should equal is therefore
#
#   r_oracle = cor(pi_d, p_d)
#
# the correlation between the TRUE district prevalence and the OBSERVED one. A
# predictor that recovered the truth perfectly would score exactly this and no
# more. The simulation computes r_oracle directly, then asks what the analytic
# and empirical estimators say about the same simulated data.
#
# A NOTE ON WHAT THIS DOES NOT TEST
# ---------------------------------
# The simulation generates outcomes from a model with a district effect and a
# cluster effect. If the real world has structure this model omits, the answer
# transfers only as far as that model does. It is a test of the ESTIMATOR, not a
# validation of the data-generating assumptions.
# =============================================================================

#' Intracluster correlation, estimated from districts that hold more than one
#' cluster.
#'
#' Uses the ANOVA estimator for a binary outcome, pooled over districts, which
#' needs no model fit and degrades gracefully when most districts hold a single
#' cluster. Returns NA when no district has two clusters.
#'
#' @param y binary outcome
#' @param cl cluster id
#' @param a2 district id
#' @return list(icc, n_districts_used, n_clusters_used)
estimate_icc <- function(y, cl, a2) {
  ok <- is.finite(y) & !is.na(cl) & !is.na(a2)
  y <- y[ok]; cl <- cl[ok]; a2 <- a2[ok]
  ncl <- tapply(cl, a2, function(z) length(unique(z)))
  use <- names(which(ncl >= 2L))
  if (!length(use)) return(list(icc = NA_real_, n_districts_used = 0L,
                                n_clusters_used = 0L))
  sel <- a2 %in% use
  y <- y[sel]; cl <- cl[sel]; a2 <- a2[sel]
  # One-way ANOVA of y on cluster, nested inside district: remove the district
  # mean first so between-DISTRICT variation is not counted as between-cluster.
  dm <- tapply(y, a2, mean)
  yc <- y - dm[a2]
  g  <- paste(a2, cl, sep = "\r")
  nj <- table(g)
  if (length(nj) < 2) return(list(icc = NA_real_, n_districts_used = length(use),
                                  n_clusters_used = length(nj)))
  mj <- tapply(yc, g, mean)
  msb <- sum(nj * (mj - mean(yc))^2) / (length(nj) - 1)
  wss <- sum((yc - mj[g])^2)
  dfw <- length(yc) - length(nj)
  msw <- if (dfw > 0) wss / dfw else NA_real_
  n0  <- (sum(nj) - sum(nj^2) / sum(nj)) / (length(nj) - 1)
  if (!is.finite(msw) || msw <= 0 || !is.finite(n0) || n0 <= 0)
    return(list(icc = NA_real_, n_districts_used = length(use),
                n_clusters_used = length(nj)))
  s2b <- max((msb - msw) / n0, 0)
  list(icc = s2b / (s2b + msw), n_districts_used = length(use),
       n_clusters_used = length(nj))
}

#' One simulation replicate: known truth, real survey structure.
#' @keywords internal
.sim_once <- function(a2, cl, w, sd_logit, icc, base_p) {
  du <- sort(unique(a2)); cu <- unique(paste(a2, cl, sep = "\r"))
  # True district prevalence: logit-normal around the cell's own base rate, so
  # the simulated prevalences occupy the same part of the scale as the real ones.
  mu <- stats::qlogis(min(max(base_p, 0.01), 0.99))
  pi_d <- stats::plogis(mu + stats::rnorm(length(du), 0, sd_logit))
  names(pi_d) <- du
  # Cluster effect scaled so the within-district intracluster correlation is icc.
  # On the probability scale, Var_between_cluster = icc * p(1-p).
  key <- paste(a2, cl, sep = "\r")
  pbar <- pi_d[a2]
  sd_c <- sqrt(pmax(icc, 0) * pbar * (1 - pbar))
  eps  <- stats::rnorm(length(cu), 0, 1); names(eps) <- cu
  p_i  <- pmin(pmax(pbar + sd_c * eps[key], 0.001), 0.999)
  y    <- stats::rbinom(length(p_i), 1L, p_i)
  list(y = y, pi_d = pi_d)
}

#' Bias of the ceiling estimators against a known attainable correlation.
#'
#' @param d individual-level data supplying the real survey structure
#' @param a2_col,cl_col,w_col,y_col column names
#' @param sd_logit between-district standard deviation of the true prevalence on
#'   the logit scale; sweeping this sweeps the true reliability
#' @param R replicates
#' @param deff design effect handed to the analytic estimator, matching what
#'   admin2_reliability() assumes in production
#' @return one row per (sd_logit) with the mean of each quantity over replicates
simulate_ceiling_bias <- function(d, a2_col, cl_col, w_col, y_col,
                                  sd_logit = c(0.2, 0.5, 1.0),
                                  R = 100L, deff = 1.5, seed = 20260902L,
                                  split_B = 40L) {
  y0 <- suppressWarnings(as.numeric(haven::zap_labels(d[[y_col]])))
  a2 <- as.character(d[[a2_col]])
  cl <- if (cl_col %in% names(d)) as.character(d[[cl_col]]) else a2
  w  <- if (w_col %in% names(d))
    suppressWarnings(as.numeric(haven::zap_labels(d[[w_col]]))) else rep(1, nrow(d))
  w[!is.finite(w) | w <= 0] <- NA_real_
  ok <- is.finite(y0) & !is.na(a2) & !is.na(cl) & is.finite(w)
  y0 <- y0[ok]; a2 <- a2[ok]; cl <- cl[ok]; w <- w[ok]
  if (length(y0) < 100) return(NULL)

  ic <- estimate_icc(y0, cl, a2)
  icc <- if (is.finite(ic$icc)) ic$icc else 0
  base_p <- mean(y0)

  set.seed(seed)
  out <- list()
  for (s in sd_logit) {
    r_or <- r_an <- r_em <- lam_true <- rep(NA_real_, R)
    n_zero <- 0L
    for (r in seq_len(R)) {
      sm <- .sim_once(a2, cl, w, s, icc, base_p)
      ys <- sm$y; pi_d <- sm$pi_d
      lev <- names(pi_d)
      g <- factor(a2, levels = lev)
      p_obs <- .wprev(ys, w, g, lev)
      n_obs <- as.integer(table(g))
      fin <- is.finite(p_obs) & n_obs >= 2
      if (sum(fin) < 5 || stats::sd(p_obs[fin]) == 0) next

      # The quantity r_max is supposed to bound.
      r_or[r] <- suppressWarnings(stats::cor(pi_d[fin], p_obs[fin]))
      # The true reliability of the observed estimate, by construction.
      lam_true[r] <- stats::var(pi_d[fin]) / stats::var(p_obs[fin])

      # The analytic estimator, exactly as production computes it.
      sv <- data.frame(svy_prev = as.numeric(p_obs[fin]),
                       n_svy = n_obs[fin])
      an <- admin2_reliability(sv, deff = deff, boot = 0)
      if (!is.null(an)) {
        r_an[r] <- an$r_max
        if (is.finite(an$r_max) && an$r_max == 0) n_zero <- n_zero + 1L
      }
      # The empirical estimator, on the same simulated data.
      dd <- data.frame(a2 = a2, cl = cl, w = w, y = ys, stringsAsFactors = FALSE)
      em <- tryCatch(split_half_reliability(dd, "a2", "cl", "w", "y",
                                            scheme = "within", B = split_B,
                                            seed = seed + r),
                     error = function(e) NULL)
      if (!is.null(em)) r_em[r] <- em$r_max_emp
    }
    m <- function(x) round(mean(x, na.rm = TRUE), 4)
    out[[length(out) + 1L]] <- data.frame(
      sd_logit = s, replicates = sum(is.finite(r_or)),
      icc = round(icc, 4), icc_districts = ic$n_districts_used,
      lambda_true = m(lam_true),
      r_oracle = m(r_or),
      r_max_analytic = m(r_an),
      r_max_empirical = m(r_em),
      bias_analytic = round(mean(r_an - r_or, na.rm = TRUE), 4),
      bias_empirical = round(mean(r_em - r_or, na.rm = TRUE), 4),
      pct_analytic_zero = round(100 * n_zero / max(R, 1), 1),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
