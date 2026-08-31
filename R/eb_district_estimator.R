# =============================================================================
# R/eb_district_estimator.R
#
# WS-B3. The shipping district estimator: an empirical Bayes blend between the
# survey's own district estimate and its region's estimate.
#
# WHY THIS ONE
# ------------
# The tournament in R/estimator_tournament.R scores every candidate against a
# SIMULATED TRUTH rather than against the survey's noisy district estimates.
# That distinction is the whole reason this estimator was not chosen earlier:
# scored against the observed estimate, the direct estimate has a correlation of
# exactly 1 by construction and was therefore never treated as a candidate, and
# any estimator that shrinks was penalised for the yardstick's own sampling
# noise.
#
# Against truth the ordering is different. The empirical Bayes blend is first on
# mean absolute error at every level of covariate signal tested, and second on
# correlation, close behind the direct estimate. It beats the flat regional mean
# on both metrics at every setting, and it beats the covariate ridge everywhere
# in the range of covariate signal this project has actually measured.
#
# WHAT IT IS NOT
# --------------
# It is not a covariate model. It uses no predictor at all: only the survey's own
# district estimates, their sample sizes, and the region each district sits in.
# WS-A found two predictors surviving a family-wise permutation correction out of
# 294, and WS-C3 found the one mechanistically motivated covariate set failing.
# Shipping a covariate model on that evidence would not be supportable.
#
# THE JACKKNIFE IS NOT OPTIONAL
# -----------------------------
# The regional target is computed from the region's OTHER districts. Without
# that, a district contributes to the quantity it is shrunk toward, which is the
# defect that withdrew the anchoring claim (register 4.6) and then the
# covariate-free regional claim (register 4.12). The un-jackknifed variant
# scored mean r 0.516 and the jackknifed one 0.076 on the same cells.
# =============================================================================

#' Empirical Bayes district estimates, shrunk toward a jackknifed region mean.
#'
#' @param svy data.frame with Admin1, Admin2, svy_prev and n_svy
#' @param deff design effect used for the per-district sampling variance
#' @param min_region_districts a region needs at least this many OTHER districts
#'   for a jackknifed regional target to be formed; below it the district falls
#'   back to the national mean of the other districts
#' @return `svy` with added columns:
#'   region_target   the jackknifed regional (or national) shrink target
#'   v_d             per-district sampling variance, Agresti-Coull stabilised
#'   tau2            between-district variance net of sampling noise
#'   lambda          weight on the district's own estimate
#'   eb_prev         the blended estimate
#'   eb_source       "region" or "national", recording which target was used
#' @param lambda_emp the cell's empirical split-half reliability from
#'   results/tables/reliability_empirical.csv. Supplying it is strongly
#'   preferred; see the tau2 note in the body for what happens without it.
eb_district_estimate <- function(svy, deff = 1.5, min_region_districts = 2L,
                                 lambda_emp = NA_real_) {
  need <- c("Admin1", "Admin2", "svy_prev", "n_svy")
  if (is.null(svy) || !all(need %in% names(svy))) return(NULL)
  d <- svy
  p <- as.numeric(d$svy_prev); n <- pmax(as.numeric(d$n_svy), 1)
  reg <- trimws(as.character(d$Admin1))
  ok <- is.finite(p) & is.finite(n)
  if (sum(ok) < 5) return(NULL)

  # Per-district sampling variance. The Agresti-Coull stabilisation matters:
  # the plug-in p(1-p)/n collapses to exactly zero for a district observed at
  # 0 or 100 percent, and inverting that hands it infinite weight, so such a
  # district would never be shrunk at all.
  v_d <- area_sampling_var_shrunk(p, n, deff = deff)
  v_d[!is.finite(v_d) | v_d <= 0] <- NA_real_

  # ── tau2, and the trap that the moment estimator walks into ───────────────
  #
  # The obvious estimator is var(p) - mean(v_d), the observed between-district
  # variance net of sampling noise. It is the same quantity the analytic
  # reliability ceiling uses, and it inherits the same defect: WS1 measured that
  # estimator returning exactly zero in 10 of 24 cells where the empirical
  # split-half reliability finds real signal in 21 of 24.
  #
  # Floored at zero, tau2 = 0 forces lambda = 0 for every district, and the
  # blend collapses to the flat regional mean. Measured on this data before the
  # fix: lambda was exactly 0.000 in FOURTEEN of the 24 cells, so the shipping
  # estimator would have silently become the very arm that WS-B1 withdrew at a
  # jackknifed mean r of 0.076.
  #
  # The fix uses this session's own result. reliability_empirical.csv holds a
  # split-half lambda per cell that needs no design-effect assumption and no
  # truncation, and reliability IS tau2 / (tau2 + vbar), so
  #
  #     tau2 = lambda_emp / (1 - lambda_emp) * vbar
  #
  # The moment estimator remains the fallback where no empirical reliability is
  # available, and which route was taken is recorded on every row.
  vbar <- mean(v_d[ok], na.rm = TRUE)
  tau2_moment <- max(stats::var(p[ok]) - vbar, 0)
  tau2 <- tau2_moment; tau2_source <- "moment"
  if (is.finite(lambda_emp) && lambda_emp > 0 && lambda_emp < 0.999) {
    tau2 <- lambda_emp / (1 - lambda_emp) * vbar
    tau2_source <- "split_half_reliability"
  }
  if (!is.finite(tau2) || tau2 <= 0) { tau2 <- 1e-8; tau2_source <- "degenerate" }

  target <- rep(NA_real_, nrow(d)); src <- rep(NA_character_, nrow(d))
  for (i in which(ok)) {
    j <- which(ok & reg == reg[i] & seq_len(nrow(d)) != i)
    if (length(j) >= min_region_districts) {
      target[i] <- stats::weighted.mean(p[j], n[j]); src[i] <- "region"
    } else {
      j2 <- which(ok & seq_len(nrow(d)) != i)
      if (length(j2) >= min_region_districts) {
        target[i] <- stats::weighted.mean(p[j2], n[j2]); src[i] <- "national"
      }
    }
  }
  lam <- tau2 / (tau2 + v_d)
  lam[!is.finite(lam)] <- 1                     # no variance estimate: do not shrink
  eb <- ifelse(is.finite(target), lam * p + (1 - lam) * target, p)

  d$region_target <- round(target, 6)
  d$v_d      <- round(v_d, 8)
  d$tau2     <- round(tau2, 8)
  d$tau2_source <- tau2_source
  d$lambda_emp  <- round(lambda_emp, 4)
  d$lambda   <- round(lam, 4)
  d$eb_prev  <- round(eb, 6)
  d$eb_source <- src
  d
}

#' District rank intervals from the blended estimate.
#'
#' A point rank from an estimator at this precision is not distinguishable from
#' many neighbouring ranks, so `docs/FITNESS_FOR_USE.md` requires a rank interval
#' wherever a district ranking is shown. The interval is a parametric bootstrap
#' on the blended estimate's own posterior variance, which for this blend is
#' lambda * v_d.
#'
#' @param eb output of eb_district_estimate()
#' @param B bootstrap replicates
#' @param seed random seed
#' @return `eb` with rank, rank_lo and rank_hi added (1 = highest prevalence)
eb_rank_intervals <- function(eb, B = 2000L, seed = 20260924L) {
  if (is.null(eb) || !all(c("eb_prev", "lambda", "v_d") %in% names(eb))) return(eb)
  p <- eb$eb_prev; se <- sqrt(pmax(eb$lambda * eb$v_d, 0))
  ok <- is.finite(p) & is.finite(se)
  if (sum(ok) < 3) return(eb)
  set.seed(seed)
  R <- matrix(NA_real_, sum(ok), B)
  pp <- p[ok]; ss <- se[ok]
  for (b in seq_len(B))
    R[, b] <- rank(-(pp + stats::rnorm(length(pp), 0, ss)), ties.method = "average")
  eb$rank <- NA_real_; eb$rank_lo <- NA_real_; eb$rank_hi <- NA_real_
  eb$rank[ok]    <- rank(-pp, ties.method = "average")
  eb$rank_lo[ok] <- apply(R, 1, stats::quantile, 0.05, names = FALSE)
  eb$rank_hi[ok] <- apply(R, 1, stats::quantile, 0.95, names = FALSE)
  eb
}
