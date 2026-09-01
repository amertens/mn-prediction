# =============================================================================
# The Admin-2 join lint.
#
# Section 8 of docs/SESSION_FINDINGS_FOR_REVIEW.md records the tenth instance of
# one defect: joining Admin-2 tables on the district name, which fans rows in
# Malawi because six district names occur in more than one region. The document
# recommends a test that fails on `by = "Admin2"`. This is that test, as a
# ratchet rather than a ban; R/lint_admin2_joins.R explains why.
# =============================================================================

test_that("the lint finds the sites it is supposed to find", {
  # A self-test, so a broken regex fails loudly rather than reporting a clean
  # codebase. Without this the whole suite passes when the scanner matches
  # nothing at all.
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "R"))
  writeLines(c(
    'a <- dplyr::left_join(x, y, by = "Admin2")',      # by_name_only
    'b <- merge(x, y, by = c("Admin2"))',              # by_name_only
    'cc <- dplyr::inner_join(x, y, "Admin2")',         # positional_join
    'd <- merge(x, y, "Admin2")',                      # positional_merge
    '# e <- left_join(x, y, by = "Admin2")   commented, must NOT count',
    'g <- dplyr::left_join(x, y, by = c("Admin1", "Admin2"))  # pair key, fine',
    'i <- merge(x, y, c("Admin1", "Admin2"))  # positional pair key, also fine',
    'h <- dplyr::select(x, Admin2, svy_prev)  # not a join'
  ), file.path(tmp, "R", "fake.R"))

  s <- scan_admin2_joins(root = tmp, dirs = "R")
  expect_equal(nrow(s), 4L)
  expect_setequal(s$kind,
                  c("by_name_only", "by_name_only", "positional_join", "positional_merge"))
  expect_false(any(grepl("^#", s$code)))
})

test_that("comment stripping does not break on a hash inside a string", {
  expect_equal(.strip_r_comment('x <- "a#b"  # trailing'), 'x <- "a#b"  ')
  expect_equal(.strip_r_comment("# whole line"), "")
  expect_equal(.strip_r_comment('y <- 1'), 'y <- 1')
})

test_that("no name-only Admin-2 join exists outside the recorded baseline", {
  base_path <- testthat::test_path("admin2_join_baseline.csv")
  skip_if_not(file.exists(base_path), "baseline absent")
  baseline <- utils::read.csv(base_path, stringsAsFactors = FALSE)
  scan <- scan_admin2_joins()
  d <- diff_admin2_baseline(scan, baseline)

  # A NEW name-only join is a failure. Fix it with admin2_join_by(), or, if the
  # join is provably safe, add it to the baseline deliberately by running
  # scripts/accuracy_impact/ws7a_baseline.R and saying why in `audit`.
  if (nrow(d$new)) {
    msg <- paste(sprintf("  %s:%d  [%s]  %s", d$new$file, d$new$line,
                         d$new$kind, d$new$code), collapse = "\n")
    fail(paste0(nrow(d$new), " name-only Admin-2 join(s) not in the baseline:\n", msg))
  }
  expect_equal(nrow(d$new), 0L)
})

test_that("the baseline does not rot", {
  base_path <- testthat::test_path("admin2_join_baseline.csv")
  skip_if_not(file.exists(base_path), "baseline absent")
  baseline <- utils::read.csv(base_path, stringsAsFactors = FALSE)
  d <- diff_admin2_baseline(scan_admin2_joins(), baseline)
  # Fixed sites are good news, and the baseline should shrink to match so it
  # never grandfathers a join that no longer exists. This warns rather than
  # fails: a fix should not break the build.
  if (nrow(d$fixed))
    warning(sprintf("%d baseline site(s) no longer present; re-run ws7a_baseline.R",
                    nrow(d$fixed)))
  succeed()
})
