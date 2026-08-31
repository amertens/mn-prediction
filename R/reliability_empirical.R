# =============================================================================
# R/reliability_empirical.R
#
# HOW RELIABLE IS A DISTRICT'S SURVEY PREVALENCE, MEASURED RATHER THAN ASSUMED?
#
# WHY THIS EXISTS
# ---------------
# admin2_reliability() estimates the noise ceiling analytically: it subtracts an
# assumed binomial sampling variance, inflated by a design effect fixed at 1.5,
# from the observed between-district variance. Section 9 of
# docs/SESSION_FINDINGS_FOR_REVIEW.md records that the resulting ceiling is
# exceeded by measured skill in about half of all cells, by a factor of two for
# the best anchoring arm, and states the question as unresolved: is the ceiling
# biased low, or is the evaluation optimistic?
#
# An analytic estimator cannot answer a question about itself. This file
# supplies the independent measurement: split the respondents in each district
# into halves, estimate the prevalence twice, and see how well the two estimates
# agree across districts. That correlation is the reliability of a half-sized
# survey, and the Spearman-Brown formula corrects it to full size. Nothing about
# it depends on a design-effect assumption or on a binomial model.
#
# WHAT THE SPLIT DOES AND DOES NOT REPRODUCE
# ------------------------------------------
# Two schemes are computed, and the difference between them is informative.
#
#   "within"   Respondents are split inside each cluster, so both halves draw on
#              every cluster the district contains. Between-cluster variation is
#              therefore NOT reproduced in the split, and the halves agree more
#              than two independent surveys of the district would. This scheme
#              yields an UPPER bound on reliability.
#
#   "cluster"  Whole clusters are assigned to halves, which does reproduce
#              between-cluster variation. It is the design-faithful scheme and
#              yields the more honest number, but it needs at least two clusters
#              in a district. Measured 2026-08-31: the median district holds ONE
#              cluster in three of the four countries (Gambia 17 of 30 districts,
#              Ghana 62 of 75, Malawi 74 of 87 hold exactly one), so this scheme
#              is only broadly computable for Sierra Leone.
#
# The prompt for this work specifies the "within" scheme. It is the primary
# output. The bound direction is stated here because a ceiling that is an upper
# bound cannot be used to argue that a measured skill is impossible.
#
# THE TRUNCATION THAT MATTERS
# ---------------------------
# lambda can come out negative, because a correlation between two noisy halves
# can be negative by chance. admin2_reliability() applies max(0, .) and then
# takes a square root, so its r_max is truncated at zero and its sampling
# distribution is not symmetric. Measured on this data: 10 of 24 cells return
# r_max EXACTLY 0.000 (source: results/tables/frozen_2026-09/admin1_arms.csv,
# column r_max). Anything that divides by r_max inherits that. Both the
# truncated and the untruncated lambda are therefore reported here, so the
# effect of the truncation can be seen rather than inferred.
# =============================================================================

#' Weighted prevalence by group, fast enough for a few hundred replicates.
#' @keywords internal
.wprev <- function(y, w, g, lev) {
  sw  <- vapply(split(w,     g), sum, numeric(1))
  swy <- vapply(split(w * y, g), sum, numeric(1))
  out <- rep(NA_real_, length(lev)); names(out) <- lev
  nm <- names(sw)
  out[nm] <- ifelse(sw > 0, swy / sw, NA_real_)
  out
}

#' Empirical split-half reliability of district prevalence for one cell.
#'
#' @param d individual-level data frame
#' @param a2_col,cl_col,w_col,y_col column names for district, cluster, survey
#'   weight and the binary outcome
#' @param scheme "within" splits respondents inside each cluster;
#'   "cluster" assigns whole clusters to halves
#' @param B number of random splits
#' @param min_half minimum respondents a district must retain in EACH half for
#'   that district to enter the correlation
#' @param seed random seed
#' @return one-row data.frame, or NULL when the cell cannot support the scheme
#' @param min_units minimum number of areas that must survive a split for that
#'   split to contribute. Defaults to 5, which is what WS1a used; the WS4a
#'   resolution sweep lowers it to 4 so that Sierra Leone's four Admin-1 regions
#'   can be estimated at all. A reliability from four points is noisy and the
#'   sweep labels it as such.
split_half_reliability <- function(d, a2_col, cl_col, w_col, y_col,
                                   scheme = c("within", "cluster"),
                                   B = 200L, min_half = 3L, seed = 20260901L,
                                   min_units = 5L) {
  scheme <- match.arg(scheme)
  y  <- suppressWarnings(as.numeric(haven::zap_labels(d[[y_col]])))
  a2 <- as.character(d[[a2_col]])
  cl <- if (!is.null(cl_col) && cl_col %in% names(d)) as.character(d[[cl_col]]) else a2
  w  <- if (!is.null(w_col) && w_col %in% names(d))
    suppressWarnings(as.numeric(haven::zap_labels(d[[w_col]]))) else rep(1, nrow(d))
  w[!is.finite(w) | w <= 0] <- NA_real_

  ok <- is.finite(y) & !is.na(a2) & !is.na(cl) & is.finite(w)
  y <- y[ok]; a2 <- a2[ok]; cl <- cl[ok]; w <- w[ok]
  if (length(y) < 30 || length(unique(y)) < 2) return(NULL)

  # Districts must be able to give both halves something to estimate from.
  keep_d <- names(which(table(a2) >= 2L * min_half))
  if (length(keep_d) < min_units) return(NULL)
  sel <- a2 %in% keep_d
  y <- y[sel]; a2 <- a2[sel]; cl <- cl[sel]; w <- w[sel]

  # The cluster scheme needs at least two clusters in a district to be able to
  # put one on each side; districts that cannot are dropped, and the count is
  # reported so a reader can see how much of the cell survived.
  n_cl_per_d <- tapply(cl, a2, function(z) length(unique(z)))
  if (scheme == "cluster") {
    keep_d2 <- names(which(n_cl_per_d >= 2L))
    if (length(keep_d2) < min_units) return(NULL)
    sel <- a2 %in% keep_d2
    y <- y[sel]; a2 <- a2[sel]; cl <- cl[sel]; w <- w[sel]
  }

  lev <- sort(unique(a2))
  # Strata for the random assignment: cluster within district for "within",
  # district for "cluster" (whole clusters are the unit assigned).
  strat <- if (scheme == "within") paste(a2, cl, sep = "\r") else a2
  idx_by_strat <- split(seq_along(y), strat)

  set.seed(seed)
  rs <- rep(NA_real_, B); n_used <- rep(NA_integer_, B)
  for (b in seq_len(B)) {
    half <- integer(length(y))
    if (scheme == "within") {
      # Balanced assignment inside each cluster: a permutation split rather than
      # independent coin flips, so a small cluster cannot land entirely on one
      # side and leave the other half of that district unrepresented.
      for (ii in idx_by_strat) {
        n <- length(ii)
        half[ii] <- sample(rep_len(c(1L, 2L), n), n)
      }
    } else {
      for (dd in lev) {
        ic <- unique(cl[a2 == dd])
        asg <- sample(rep_len(c(1L, 2L), length(ic)), length(ic))
        names(asg) <- ic
        j <- which(a2 == dd)
        half[j] <- asg[cl[j]]
      }
    }
    i1 <- which(half == 1L); i2 <- which(half == 2L)
    n1 <- table(factor(a2[i1], levels = lev))
    n2 <- table(factor(a2[i2], levels = lev))
    good <- lev[n1 >= min_half & n2 >= min_half]
    if (length(good) < min_units) next
    p1 <- .wprev(y[i1], w[i1], factor(a2[i1], levels = lev), lev)[good]
    p2 <- .wprev(y[i2], w[i2], factor(a2[i2], levels = lev), lev)[good]
    fin <- is.finite(p1) & is.finite(p2)
    if (sum(fin) < min_units || stats::sd(p1[fin]) == 0 || stats::sd(p2[fin]) == 0) next
    rs[b] <- suppressWarnings(stats::cor(p1[fin], p2[fin]))
    n_used[b] <- sum(fin)
  }

  rs <- rs[is.finite(rs)]
  if (length(rs) < 20) return(NULL)

  # Spearman-Brown steps a half-length correlation up to full length.
  lam_raw  <- 2 * rs / (1 + rs)
  lam_raw[!is.finite(lam_raw)] <- NA_real_
  lam_raw  <- lam_raw[is.finite(lam_raw)]
  lam_trunc <- pmin(pmax(lam_raw, 0), 1)

  q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
  data.frame(
    scheme          = scheme,
    n_splits        = length(rs),
    n_areas_median  = as.integer(stats::median(n_used, na.rm = TRUE)),
    r_halfhalf      = round(stats::median(rs), 4),
    lambda_emp_raw  = round(stats::median(lam_raw), 4),
    lambda_emp      = round(stats::median(lam_trunc), 4),
    lambda_lo       = round(q(lam_trunc, 0.025), 4),
    lambda_hi       = round(q(lam_trunc, 0.975), 4),
    r_max_emp       = round(sqrt(stats::median(lam_trunc)), 4),
    r_max_emp_lo    = round(sqrt(q(lam_trunc, 0.025)), 4),
    r_max_emp_hi    = round(sqrt(q(lam_trunc, 0.975)), 4),
    frac_neg_lambda = round(mean(lam_raw < 0), 3),
    stringsAsFactors = FALSE)
}

#' Run the empirical reliability over every cell in the store.
#'
#' @return data.frame with one row per cell and scheme, carrying the analytic
#'   ceiling alongside so 1a and 1b can be read off one table
build_reliability_empirical <- function(store = here::here("_targets_full"),
                                        B = 200L, seed = 20260901L,
                                        countries = NULL, outcomes = NULL,
                                        schemes = c("within", "cluster")) {
  cfgs <- get_country_configs()
  if (!is.null(countries)) cfgs <- cfgs[intersect(names(cfgs), countries)]
  rows <- list()
  for (cn in names(cfgs)) {
    cc <- cfgs[[cn]]
    ocs <- names(cc$outcomes)
    if (!is.null(outcomes)) ocs <- intersect(ocs, outcomes)
    for (ocn in ocs) {
      oc <- cc$outcomes[[ocn]]
      sfx <- paste0(tolower(cn), "_", ocn)
      od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", sfx), store = store),
                     error = function(e) NULL)
      sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", sfx), store = store),
                     error = function(e) NULL)
      if (is.null(od) || is.null(oc$binary)) next
      d <- od$data
      if (!all(c(cc$admin2_col, oc$binary) %in% names(d))) next

      # The analytic ceiling, recomputed here on the same survey table the
      # published tables used, so 1b compares like with like.
      ana <- if (!is.null(sv)) admin2_reliability(sv, deff = 1.5, boot = 500L) else NULL

      for (s in schemes) {
        r <- tryCatch(split_half_reliability(
          d, cc$admin2_col, cc$cluster_id, cc$weight_col, oc$binary,
          scheme = s, B = B, seed = seed), error = function(e) NULL)
        if (is.null(r)) next
        r$country <- cc$country; r$outcome <- oc$tag
        r$lambda_analytic <- if (!is.null(ana)) round(ana$lambda, 4) else NA_real_
        r$r_max_analytic  <- if (!is.null(ana)) round(ana$r_max, 4) else NA_real_
        r$n_areas_analytic <- if (!is.null(ana)) ana$n_areas else NA_integer_
        r$median_n_analytic <- if (!is.null(ana)) ana$median_n else NA_real_
        r$deff_analytic <- if (!is.null(ana)) ana$deff else NA_real_
        rows[[length(rows) + 1L]] <- r
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  front <- c("country", "outcome", "scheme")
  out[, c(front, setdiff(names(out), front))]
}

#' What design effect would reconcile the analytic ceiling with the empirical one?
#'
#' Solving lambda_emp = 1 - deff * E[p(1-p)/(n-1)] / Var(p) for deff turns the
#' disagreement into a single interpretable number: the design effect the
#' analytic formula would need in order to agree with the measurement. A value
#' below 1 means the assumed 1.5 was too large and the ceiling was pushed down.
#'
#' @param sv survey Admin-2 table with svy_prev and n_svy
#' @param lambda_emp the empirical reliability for that cell
#' @return implied design effect, or NA
implied_deff <- function(sv, lambda_emp) {
  if (is.null(sv) || !is.finite(lambda_emp)) return(NA_real_)
  p <- as.numeric(sv$svy_prev); n <- as.numeric(sv$n_svy)
  ok <- is.finite(p) & is.finite(n)
  p <- p[ok]; n <- n[ok]
  if (length(p) < 5) return(NA_real_)
  v_obs <- stats::var(p)
  unit  <- mean(area_sampling_var(p, n, deff = 1))
  if (!is.finite(v_obs) || v_obs <= 0 || unit <= 0) return(NA_real_)
  (1 - lambda_emp) * v_obs / unit
}
