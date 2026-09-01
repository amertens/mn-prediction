# =============================================================================
# WS1d. Noise injection against a well-measured district quantity.
#
# Section 3 of docs/SESSION_FINDINGS_FOR_REVIEW.md carries a positive control:
# earth observation predicts district education, which shows the geospatial
# linkage works. This test uses the same kind of target for a different purpose.
# District education is well measured, so it can stand in for a known truth.
# Adding binomial noise of a known size then produces a target of known
# reliability, and the ceiling estimators can be checked against it.
#
# The property under test is the one Section 9 leaves open: does the ceiling
# track the correlation an oracle predictor actually achieves? WS1c answers it
# by simulation over the real survey structure. This test pins the answer so a
# future change to admin2_reliability() or to split_half_reliability() cannot
# silently reintroduce the bias.
#
# Runs on the smoke profile: one country, one indicator, a few hundred
# districts' worth of arithmetic and no model fitting.
# =============================================================================

.mnp_edu_truth <- function() {
  p <- here::here("data", "DHS", "clean", "Ghana_2014_dhs_custom_admin2_wide.rds")
  if (!file.exists(p)) return(NULL)
  w <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(w) || !"dhs2014_w_no_education_adm2" %in% names(w)) return(NULL)
  v <- suppressWarnings(as.numeric(w[["dhs2014_w_no_education_adm2"]]))
  v <- v[is.finite(v) & v > 0.01 & v < 0.99]
  if (length(v) < 30) return(NULL)
  v
}

test_that("the education target is available and is genuinely variable", {
  truth <- .mnp_edu_truth()
  skip_if(is.null(truth), "Ghana 2014 DHS custom Admin-2 wide file absent")
  expect_gte(length(truth), 30L)
  # Without real between-district spread the whole exercise is vacuous.
  expect_gt(stats::sd(truth), 0.05)
})

test_that("measured oracle correlation falls as injected noise rises", {
  truth <- .mnp_edu_truth()
  skip_if(is.null(truth), "education target absent")
  set.seed(20260904L)
  ns <- c(500L, 100L, 25L)
  r_oracle <- vapply(ns, function(n) {
    obs <- stats::rbinom(length(truth), n, truth) / n
    suppressWarnings(stats::cor(truth, obs))
  }, numeric(1))
  # Strictly decreasing in the amount of noise. This is the behaviour a ceiling
  # is supposed to describe; if it does not hold the rest cannot be interpreted.
  expect_true(all(diff(r_oracle) < 0),
              info = paste("r by n:", paste(round(r_oracle, 3), collapse = ", ")))
  expect_gt(r_oracle[1], 0.9)
})

test_that("the analytic ceiling under-states the attainable correlation", {
  skip_if_not(exists("admin2_reliability"), "admin2_reliability not loaded")
  truth <- .mnp_edu_truth()
  skip_if(is.null(truth), "education target absent")
  set.seed(20260905L)
  # A single moderate noise level, repeated, so the comparison is not one draw.
  n <- 40L; R <- 40L
  d_an <- numeric(R)
  for (i in seq_len(R)) {
    obs <- stats::rbinom(length(truth), n, truth) / n
    r_or <- suppressWarnings(stats::cor(truth, obs))
    an <- admin2_reliability(
      data.frame(svy_prev = obs, n_svy = rep(n, length(obs))),
      deff = 1.5, boot = 0)
    d_an[i] <- an$r_max - r_or
  }
  # Measured in WS1c over the real survey structure: mean bias -0.161 with
  # deff fixed at 1.5. The direction is the pinned property; the tolerance is
  # loose because this fixture is simple random sampling, where an assumed
  # design effect of 1.5 is by construction too large.
  expect_lt(mean(d_an), 0,
            label = sprintf("mean analytic bias %.3f", mean(d_an)))
})

test_that("the split-half estimator recovers the attainable correlation", {
  skip_if_not(exists("split_half_reliability"),
              "split_half_reliability not loaded")
  truth <- .mnp_edu_truth()
  skip_if(is.null(truth), "education target absent")
  set.seed(20260906L)
  n <- 40L
  nd <- length(truth)
  # Expand to respondent level so the split-half estimator has rows to split.
  d <- data.frame(
    a2 = rep(seq_len(nd), each = n),
    cl = rep(seq_len(nd), each = n),          # one cluster per district
    w  = 1,
    y  = stats::rbinom(nd * n, 1L, rep(truth, each = n)))
  obs <- vapply(split(d$y, d$a2), mean, numeric(1))
  r_or <- suppressWarnings(stats::cor(truth, obs))

  em <- split_half_reliability(d, "a2", "cl", "w", "y", scheme = "within",
                               B = 60L, seed = 20260907L)
  expect_false(is.null(em))
  # WS1c measures the empirical estimator's mean bias at +0.007 over the real
  # survey structure. A tolerance of 0.15 here allows for this fixture being a
  # single draw at one noise level.
  expect_lt(abs(em$r_max_emp - r_or), 0.15,
            label = sprintf("empirical %.3f against oracle %.3f",
                            em$r_max_emp, r_or))
})

test_that("the analytic ceiling collapses to zero where signal remains", {
  skip_if_not(exists("admin2_reliability"), "admin2_reliability not loaded")
  truth <- .mnp_edu_truth()
  skip_if(is.null(truth), "education target absent")
  # Compress the between-district spread toward the mean so true reliability is
  # low but not zero, then check how often the estimator reports exactly zero.
  # WS1c measured this at 56 to 97 percent of replicates in that regime, while
  # the oracle correlation was still 0.16 to 0.37.
  set.seed(20260908L)
  flat <- mean(truth) + (truth - mean(truth)) * 0.15
  n <- 25L; R <- 40L
  zero <- 0L; r_or <- numeric(R)
  for (i in seq_len(R)) {
    obs <- stats::rbinom(length(flat), n, flat) / n
    r_or[i] <- suppressWarnings(stats::cor(flat, obs))
    an <- admin2_reliability(data.frame(svy_prev = obs, n_svy = rep(n, length(obs))),
                             deff = 1.5, boot = 0)
    if (is.finite(an$r_max) && an$r_max == 0) zero <- zero + 1L
  }
  # The oracle still achieves a real correlation in this regime.
  expect_gt(mean(r_or), 0.1)
  # And the analytic estimator still calls the cell signal-free a large part of
  # the time. This is the mechanism behind "4 of 24 cells have no detectable
  # signal above noise"; the pin is that it is a property of the estimator.
  expect_gt(zero, 0L,
            label = sprintf("%d of %d replicates returned r_max exactly 0", zero, R))
})
