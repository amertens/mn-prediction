# =============================================================================
# WS-B3. The shipping district estimator.
#
# The properties pinned here are the ones that would silently break it. The
# tau2 test in particular exists because the first implementation used the
# moment estimator and degenerated to the flat regional mean in 14 of 24 cells,
# which is the arm WS-B1 withdrew at a jackknifed mean r of 0.076.
# =============================================================================

.mnp_fake_cell <- function(n_regions = 5, per_region = 6, seed = 1L) {
  set.seed(seed)
  reg <- rep(paste0("R", seq_len(n_regions)), each = per_region)
  d <- data.frame(
    Admin1 = reg,
    Admin2 = paste0("D", seq_along(reg)),
    svy_prev = pmin(pmax(rep(seq(0.1, 0.5, length.out = n_regions),
                             each = per_region) + stats::rnorm(length(reg), 0, 0.05),
                         0.01), 0.99),
    n_svy = rep(c(4, 8, 15, 25, 40, 60), length.out = length(reg)),
    stringsAsFactors = FALSE)
  d
}

test_that("lambda rises with sample size and the blend sits between its inputs", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  eb <- eb_district_estimate(.mnp_fake_cell(), lambda_emp = 0.6)
  expect_false(is.null(eb))
  ok <- is.finite(eb$lambda) & is.finite(eb$n_svy)
  expect_gt(stats::cor(eb$n_svy[ok], eb$lambda[ok]), 0.8)
  # every blended value lies between the district estimate and its target
  b <- is.finite(eb$eb_prev) & is.finite(eb$region_target)
  lo <- pmin(eb$svy_prev[b], eb$region_target[b])
  hi <- pmax(eb$svy_prev[b], eb$region_target[b])
  expect_true(all(eb$eb_prev[b] >= lo - 1e-9 & eb$eb_prev[b] <= hi + 1e-9))
})

test_that("the regional target is jackknifed", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  d <- .mnp_fake_cell()
  eb <- eb_district_estimate(d, lambda_emp = 0.6)
  # Recompute one district's target by hand from its region's OTHER districts.
  i <- 1L
  j <- which(d$Admin1 == d$Admin1[i] & seq_len(nrow(d)) != i)
  manual <- stats::weighted.mean(d$svy_prev[j], d$n_svy[j])
  expect_equal(eb$region_target[i], round(manual, 6), tolerance = 1e-6)
  # and it must NOT equal the target that includes the district itself
  incl <- stats::weighted.mean(d$svy_prev[c(i, j)], d$n_svy[c(i, j)])
  expect_false(isTRUE(all.equal(eb$region_target[i], round(incl, 6))))
})

test_that("tau2 comes from the empirical reliability when it is supplied", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  d <- .mnp_fake_cell()
  with_rel <- eb_district_estimate(d, lambda_emp = 0.6)
  expect_equal(with_rel$tau2_source[1], "split_half_reliability")
  # Without a reliability the estimator falls back to the moment route, and
  # where that route yields a non-positive tau2 it says "degenerate" rather than
  # pretending. Both are the fallback; neither may claim the reliability route.
  without  <- eb_district_estimate(d, lambda_emp = NA_real_)
  expect_true(without$tau2_source[1] %in% c("moment", "degenerate"))
  expect_false(without$tau2_source[1] == "split_half_reliability")
  # A cell with wide between-district spread and large samples exercises the
  # moment branch proper.
  wide <- .mnp_fake_cell(n_regions = 6, per_region = 8, seed = 3L)
  wide$svy_prev <- pmin(pmax(rep(seq(0.05, 0.75, length.out = 6), each = 8), 0.01), 0.99)
  wide$n_svy <- 200L
  expect_equal(eb_district_estimate(wide, lambda_emp = NA_real_)$tau2_source[1], "moment")
})

test_that("the moment estimator can collapse where the reliability route does not", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  # A cell with no real between-district variation: every district shares one
  # true rate and differs only by sampling. var(p) - mean(v_d) goes to zero or
  # below, so the moment route floors tau2 and lambda collapses to zero.
  set.seed(7)
  d <- data.frame(Admin1 = rep(paste0("R", 1:5), each = 6),
                  Admin2 = paste0("D", 1:30),
                  n_svy = rep(8L, 30), stringsAsFactors = FALSE)
  d$svy_prev <- stats::rbinom(30, 8, 0.3) / 8
  m <- eb_district_estimate(d, lambda_emp = NA_real_)
  r <- eb_district_estimate(d, lambda_emp = 0.5)
  expect_lt(stats::median(m$lambda, na.rm = TRUE),
            stats::median(r$lambda, na.rm = TRUE))
  # This is the failure mode that mattered: with the moment route the blend
  # moves materially toward the flat regional mean, which is the withdrawn arm.
  # The property is the GAP between routes, not an absolute level, because the
  # moment estimate is itself noisy at thirty districts.
  expect_gt(stats::median(r$lambda, na.rm = TRUE) -
              stats::median(m$lambda, na.rm = TRUE), 0.1)
})

test_that("a zero-reliability cell shrinks fully, which is correct", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  eb <- eb_district_estimate(.mnp_fake_cell(), lambda_emp = 0)
  # lambda_emp of exactly 0 means the survey's two halves do not agree about
  # which district is worse, so the district's own estimate carries nothing and
  # full shrinkage is the right answer. Those cells are suppressed downstream.
  expect_lt(stats::median(eb$lambda, na.rm = TRUE), 0.01)
})

test_that("a district observed at zero percent is not given infinite weight", {
  skip_if_not(exists("eb_district_estimate"), "estimator not loaded")
  d <- .mnp_fake_cell()
  d$svy_prev[1] <- 0; d$n_svy[1] <- 1L
  eb <- eb_district_estimate(d, lambda_emp = 0.6)
  # The plug-in p(1-p)/n is exactly zero here; without stabilisation lambda
  # would be 1 and a one-respondent district would be published unshrunk.
  expect_lt(eb$lambda[1], 0.5)
  expect_gt(eb$eb_prev[1], 0.02)
})

test_that("rank intervals are ordered and cover the point rank", {
  skip_if_not(exists("eb_rank_intervals"), "rank intervals not loaded")
  eb <- eb_rank_intervals(eb_district_estimate(.mnp_fake_cell(), lambda_emp = 0.6),
                          B = 300L)
  ok <- is.finite(eb$rank)
  expect_true(all(eb$rank_lo[ok] <= eb$rank[ok] + 1e-6))
  expect_true(all(eb$rank_hi[ok] >= eb$rank[ok] - 1e-6))
  # A ranking at this precision must not claim to separate every district.
  expect_gt(stats::median(eb$rank_hi[ok] - eb$rank_lo[ok]), 1)
})
