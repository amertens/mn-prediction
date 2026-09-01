# =============================================================================
# The LOCO extrapolation diagnostic.
#
# The failure this must not have is the one that makes the measure useless
# rather than wrong: reporting everything as outlying (which halfspace depth
# does at this p and n), or reporting nothing as outlying because a direction
# with zero spread divided by zero.
# =============================================================================

.mnp_lo_cloud <- function(n = 150, p = 40, shift = 0, seed = 5L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  # give it correlated structure, so the joint question differs from the
  # marginal one
  X[, 2] <- X[, 1] * 0.9 + stats::rnorm(n, 0, 0.4)
  X[, 3] <- X[, 1] * 0.8 + stats::rnorm(n, 0, 0.5)
  X + shift
}

test_that("a point drawn from the training distribution is not flagged", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  X <- .mnp_lo_cloud(n = 180)
  lo <- loco_outlyingness(X[1:150, ], X[151:180, ])
  expect_false(is.null(lo))
  # the held-out rows come from the SAME distribution, so almost none of them
  # should exceed the training maximum
  expect_lt(mean(lo$out_test > lo$thresh_max), 0.25)
})

test_that("a shifted cloud is flagged as outlying", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  X <- .mnp_lo_cloud(n = 150)
  Y <- .mnp_lo_cloud(n = 30, shift = 6, seed = 9L)
  lo <- loco_outlyingness(X, Y)
  expect_gt(mean(lo$flag_test), 0.8)
  expect_gt(stats::median(lo$out_test), stats::median(lo$out_train))
})

test_that("it catches a JOINT outlier that every marginal range accepts", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  # Columns 1 and 2 are strongly correlated in training. A test point at
  # (+2, -2) sits inside both marginal ranges and far off the correlation
  # ridge. This is the case the marginal check provably cannot see, and the
  # reason the joint measure was ported at all.
  set.seed(3)
  n <- 200
  a <- stats::rnorm(n)
  X <- cbind(a, a + stats::rnorm(n, 0, 0.15), matrix(stats::rnorm(n * 8), n, 8))
  te <- matrix(0, nrow = 1, ncol = ncol(X))
  te[1, 1] <- 2; te[1, 2] <- -2
  expect_lt(loco_marginal_out(X, te, tol = 0.10), 0.5)   # marginal: not flagged
  lo <- loco_outlyingness(X, te)
  expect_true(lo$flag_test[1])                           # joint: flagged
  # and it fires on the ORTHOGONAL distance, which is the diagnostic-specific
  # claim: the point breaks the correlation structure rather than being an
  # extreme version of it.
  expect_gt(lo$od_test[1], lo$cut_od)
})

test_that("the measure does not saturate the way halfspace depth does", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  X <- .mnp_lo_cloud(n = 160, p = 60)
  lo <- loco_outlyingness(X[1:130, ], X[131:160, ])
  # a useful measure separates points; a saturated one returns near-constant
  # values and cannot rank anything
  expect_gt(stats::sd(lo$out_test), 0)
  expect_gt(length(unique(round(lo$out_test, 3))), 5)
})

test_that("a zero-spread direction cannot divide by zero", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  X <- .mnp_lo_cloud(n = 120, p = 20)
  X[, 5] <- 1                       # constant column
  lo <- loco_outlyingness(X[1:100, ], X[101:120, ])
  expect_false(is.null(lo))
  expect_true(all(is.finite(lo$out_test)))
})

test_that("the reduction is fitted on training only", {
  skip_if_not(exists("loco_outlyingness"), "diagnostic not loaded")
  X <- .mnp_lo_cloud(n = 150)
  te1 <- .mnp_lo_cloud(n = 20, seed = 11L)
  te2 <- rbind(te1, .mnp_lo_cloud(n = 20, shift = 9, seed = 12L))
  a <- loco_outlyingness(X, te1)
  b <- loco_outlyingness(X, te2)
  # adding wildly outlying rows to the TEST set must not move the scores of the
  # rows already there; if it does, the test set is defining its own reference
  expect_equal(a$out_test, b$out_test[seq_along(a$out_test)], tolerance = 1e-8)
  expect_equal(a$cut_od, b$cut_od, tolerance = 1e-8)
})
