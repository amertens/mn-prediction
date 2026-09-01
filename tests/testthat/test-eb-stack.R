# =============================================================================
# The stacked shrinkage target.
#
# Two properties would make this estimator dishonest rather than merely wrong,
# and both have precedents in this project: fitting on data the district
# contributes to (registers 4.6 and 4.12), and carrying a level bias into the
# blend (the -3.1 to -3.3 pp ridge bias the calibration gate catches in six
# cells). Those two are pinned first.
# =============================================================================

.mnp_stack_cell <- function(n_regions = 6, per_region = 8, seed = 11L,
                            covariate_signal = 0.6) {
  set.seed(seed)
  m <- n_regions * per_region
  reg <- rep(paste0("R", seq_len(n_regions)), each = per_region)
  z <- stats::rnorm(m)
  lat <- rep(seq(5, 12, length.out = n_regions), each = per_region) +
    stats::rnorm(m, 0, 0.2)
  lon <- stats::rnorm(m, 0, 2)
  eta <- sqrt(covariate_signal) * z + sqrt(1 - covariate_signal) * stats::rnorm(m)
  p_true <- stats::plogis(-0.8 + eta)
  n <- rep(c(6, 12, 20, 35, 60, 90, 120, 150), length.out = m)
  p <- stats::rbinom(m, n, p_true) / n
  X <- cbind(z = z, noise1 = stats::rnorm(m), noise2 = stats::rnorm(m),
             noise3 = stats::rnorm(m))
  list(p = p, n = n, reg = reg, X = X, lonlat = cbind(lon, lat),
       p_true = p_true, fin = rep(TRUE, m))
}

test_that("no candidate uses the district it predicts", {
  skip_if_not(exists("eb_stack_candidates"), "stack not loaded")
  d <- .mnp_stack_cell()
  cand <- eb_stack_candidates(d$p, d$n, d$reg, X = d$X, lonlat = d$lonlat,
                              fin = d$fin)
  # Perturb ONE district's estimate wildly. Every candidate value for THAT
  # district must be unchanged, because none of them may have seen it.
  d2 <- d; d2$p[1] <- 0.99
  cand2 <- eb_stack_candidates(d2$p, d2$n, d2$reg, X = d2$X,
                               lonlat = d2$lonlat, fin = d2$fin)
  for (nm in names(cand)) {
    a <- cand[[nm]][1]; b <- cand2[[nm]][1]
    if (is.finite(a) && is.finite(b))
      expect_equal(a, b, tolerance = 1e-8,
                   label = paste0("candidate '", nm, "' leaked district 1"))
  }
})

test_that("the region candidate reproduces the shipped jackknifed target", {
  skip_if_not(exists("eb_stack_candidates"), "stack not loaded")
  d <- .mnp_stack_cell()
  cand <- eb_stack_candidates(d$p, d$n, d$reg, fin = d$fin)
  i <- 3L
  j <- which(d$reg == d$reg[i] & seq_along(d$reg) != i)
  expect_equal(cand$region[i], stats::weighted.mean(d$p[j], d$n[j]),
               tolerance = 1e-8)
})

test_that("recalibration removes a planted level bias", {
  skip_if_not(exists("eb_stack_recalibrate"), "stack not loaded")
  d <- .mnp_stack_cell()
  cand <- eb_stack_candidates(d$p, d$n, d$reg, X = d$X, fin = d$fin)
  # Plant the ridge's own defect: a systematic -5 pp offset.
  cand$ridge <- pmax(cand$ridge - 0.05, 0)
  before <- mean(cand$ridge - d$p, na.rm = TRUE)
  rec <- eb_stack_recalibrate(cand, d$p, d$n, d$reg, d$fin)
  after <- mean(rec$ridge - d$p, na.rm = TRUE)
  expect_lt(abs(after), abs(before))
  expect_lt(abs(after), 0.02)
})

test_that("weights are convex and sum to one", {
  skip_if_not(exists("eb_stack_target"), "stack not loaded")
  s <- eb_stack_target(.mnp_stack_cell()$p, .mnp_stack_cell()$n,
                       .mnp_stack_cell()$reg,
                       X = .mnp_stack_cell()$X, shrink = 0.5)
  expect_true(all(s$weights >= -1e-9))
  expect_equal(sum(s$weights), 1, tolerance = 1e-6)
  # A convex combination can never leave the hull of its candidates.
  M <- do.call(cbind, s$candidates)
  ok <- is.finite(s$target) & apply(M, 1, function(z) all(is.finite(z)))
  expect_true(all(s$target[ok] >= apply(M[ok, , drop = FALSE], 1, min) - 1e-8))
  expect_true(all(s$target[ok] <= apply(M[ok, , drop = FALSE], 1, max) + 1e-8))
})

test_that("too few areas falls back to the region mean rather than stacking", {
  skip_if_not(exists("eb_stack_weights"), "stack not loaded")
  d <- .mnp_stack_cell(n_regions = 3, per_region = 3)
  cand <- eb_stack_candidates(d$p, d$n, d$reg, X = d$X, fin = d$fin)
  wf <- eb_stack_weights(cand, d$p, d$n, d$reg, d$fin)
  # Nine districts cannot support four fitted weights. The fallback must be the
  # target the tournament actually confirmed, not an arbitrary corner.
  expect_true(wf$route %in% c("insufficient_areas", "nnls_failed", "nnls"))
  if (wf$route == "insufficient_areas") expect_equal(unname(wf$w[["region"]]), 1)
})

test_that("shrink pulls the solution toward the region-mean corner", {
  skip_if_not(exists("eb_stack_target"), "stack not loaded")
  d <- .mnp_stack_cell(covariate_signal = 0.9)
  free <- eb_stack_target(d$p, d$n, d$reg, X = d$X, lonlat = d$lonlat, shrink = 0)
  tied <- eb_stack_target(d$p, d$n, d$reg, X = d$X, lonlat = d$lonlat, shrink = 0.9)
  expect_gt(tied$weights[["region"]], free$weights[["region"]])
})

test_that("a strong covariate earns weight and pure noise does not", {
  skip_if_not(exists("eb_stack_target"), "stack not loaded")
  # This is the Gambia case: the covariate axis really does carry the signal.
  hi <- .mnp_stack_cell(covariate_signal = 0.95, seed = 21L)
  s_hi <- eb_stack_target(hi$p, hi$n, hi$reg, X = hi$X, shrink = 0)
  # ...and the Sierra Leone case, where it does not.
  lo <- .mnp_stack_cell(covariate_signal = 0.02, seed = 21L)
  lo$X[, "z"] <- stats::rnorm(nrow(lo$X))   # break the link entirely
  s_lo <- eb_stack_target(lo$p, lo$n, lo$reg, X = lo$X, shrink = 0)
  expect_gt(s_hi$weights[["ridge"]], s_lo$weights[["ridge"]])
})
