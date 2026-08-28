# =============================================================================
# R/area_reliability.R
#
# How much of the between-Admin-2 variation in survey prevalence is real signal,
# and what ceiling does that put on any model's correlation?
#
# WHY THIS EXISTS
# ---------------
# A median Admin-2 district in these surveys contributes 6-36 biomarker
# measurements. Much of the spread in `svy_prev` across districts is therefore
# sampling noise, not geography. Reporting a Pearson r of 0.30 without saying
# whether the achievable maximum was 0.9 or 0.0 makes a saturated model and a
# useless one look identical.
#
#   lambda = (Var(p) - E[sampling variance]) / Var(p)     reliability
#   r_max  = sqrt(lambda)                                 NOISE CEILING on Pearson r
#
# NAMING. r_max is referred to as the NOISE CEILING in the docs, which is the
# term this quantity carries in the encoding-model literature. The statistics
# are Spearman's 1904 correction for attenuation: a model is scored against a
# noisy estimate of district prevalence, so even a model that predicted the
# truth exactly could only correlate sqrt(lambda) with the yardstick. The
# function and column names below keep `reliability` and `r_max` because five
# committed result tables depend on them.
#
# NOTE ON DESIGN EFFECTS. `svy_prev_se` cannot be used to estimate deff at this
# resolution: most districts hold a single PSU, so survey::svymean cannot form a
# between-cluster variance and returns ~1e-17. deff is therefore an assumption,
# defaulting to 1.5, and r_max should be read as a band -- see
# admin2_reliability(deff = ...) and sandbox_parsimony/R/01_noise_audit.R.
# =============================================================================

#' Sampling variance of an area prevalence, for AVERAGING over areas.
#'
#' E[p(1-p)] = pi(1-pi)(n-1)/n for a binomial proportion, so p(1-p)/(n-1) is
#' unbiased for pi(1-pi)/n. Individual areas observed at 0% still return 0, but
#' the mean over areas -- all the variance decomposition needs -- is unbiased.
#'
#' @param p vector of area prevalences
#' @param n vector of area sample sizes
#' @param deff assumed design effect
area_sampling_var <- function(p, n, deff = 1.5) {
  n <- pmax(n, 2)
  deff * p * (1 - p) / (n - 1)
}

#' Per-area sampling variance, for PRECISION WEIGHTS.
#'
#' The plug-in p(1-p)/n collapses to exactly zero for any area observed at 0% or
#' 100%, and inverting that hands those areas infinite weight. Use the
#' Agresti-Coull shrunk proportion so a district observed at 0/6 is treated as
#' about as informative as six observations allow.
area_sampling_var_shrunk <- function(p, n, deff = 1.5) {
  n <- pmax(n, 1)
  p_s <- (p * n + 0.5) / (n + 1)
  deff * p_s * (1 - p_s) / n
}

#' Reliability of a set of Admin-2 survey prevalences.
#'
#' @param svy_admin2 data.frame with svy_prev and n_svy
#' @param deff assumed design effect (default 1.5)
#' @param boot bootstrap replicates for the interval (0 disables)
#' @return list(lambda, r_max, lambda_lo, lambda_hi, r_max_hi, n_areas, deff)
admin2_reliability <- function(svy_admin2, deff = 1.5, boot = 500L) {
  if (is.null(svy_admin2) || !all(c("svy_prev", "n_svy") %in% names(svy_admin2)))
    return(NULL)
  p <- as.numeric(svy_admin2$svy_prev); n <- as.numeric(svy_admin2$n_svy)
  ok <- is.finite(p) & is.finite(n)
  p <- p[ok]; n <- n[ok]
  if (length(p) < 5) return(NULL)

  lam1 <- function(pp, nn) {
    v_obs <- stats::var(pp)
    if (!is.finite(v_obs) || v_obs <= 0) return(NA_real_)
    max(0, v_obs - mean(area_sampling_var(pp, nn, deff))) / v_obs
  }
  lam <- lam1(p, n)

  lo <- hi <- NA_real_
  if (boot > 0 && length(p) >= 8) {
    set.seed(4242L)
    bs <- replicate(boot, { i <- sample.int(length(p), replace = TRUE)
                            lam1(p[i], n[i]) })
    bs <- bs[is.finite(bs)]
    if (length(bs) > 20) {
      lo <- stats::quantile(bs, .025, names = FALSE)
      hi <- stats::quantile(bs, .975, names = FALSE)
    }
  }
  list(lambda = lam, r_max = sqrt(lam),
       lambda_lo = lo, lambda_hi = hi,
       r_max_hi = if (is.finite(hi)) sqrt(hi) else NA_real_,
       n_areas = length(p), median_n = stats::median(n), deff = deff)
}

#' Attach the ceiling to a metrics table produced elsewhere.
#'
#' Adds r_max, r_max_hi and r_share = pearson_r / r_max. A cell whose r_max is
#' indistinguishable from zero has no predictable between-district variation at
#' all; `signal` flags that so a table can grey it out rather than presenting a
#' meaningless correlation as a model failure.
#'
#' @param metrics data.frame containing pearson_r
#' @param svy_admin2 the survey table the metrics were scored against
add_reliability_columns <- function(metrics, svy_admin2, deff = 1.5) {
  if (is.null(metrics) || !nrow(metrics)) return(metrics)
  rel <- admin2_reliability(svy_admin2, deff = deff)
  if (is.null(rel)) {
    metrics$r_max <- NA_real_; metrics$r_share <- NA_real_
    metrics$signal <- NA; return(metrics)
  }
  metrics$r_max    <- round(rel$r_max, 3)
  metrics$r_max_hi <- round(rel$r_max_hi, 3)
  metrics$r_share  <- if (is.finite(rel$r_max) && rel$r_max > 0.05)
    round(metrics$pearson_r / rel$r_max, 2) else NA_real_
  # "is there anything here to predict?" -- uses the OPTIMISTIC bound so a cell
  # is only called signal-free when even the upper bound says so.
  metrics$signal <- is.finite(rel$r_max_hi) && rel$r_max_hi > 0.15
  metrics
}
