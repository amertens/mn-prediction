# =============================================================================
# Malawi selenium and iodine outcomes (MW-SE / MW-IO, 2026-09-15).
#
# The MNS 2015-16 measured plasma selenium (all groups) and urinary iodine
# (women, school-age children). Three Malawi-only outcomes were added to
# R/config.R; their binary indicators are derived in load_merged_data() from
# the configured cut-offs, the way folate_def and b12_def already are.
# =============================================================================

if (file.exists(here::here("R", "malawi_outcomes.R"))) source(here::here("R", "malawi_outcomes.R"))

test_that("the Malawi config carries the selenium and iodine outcomes with the documented cut-offs", {
  src <- readLines(here::here("R", "config.R"), warn = FALSE)
  expect_true(any(grepl('tag *= *"child_selenium"', src)))
  expect_true(any(grepl('tag *= *"women_selenium"', src)))
  expect_true(any(grepl('tag *= *"women_iodine"', src)))
  expect_true(any(grepl("cutoff *= *84\\.6", src)))
})

test_that("derive_malawi_binary() applies a less-than cut-off and keeps NA", {
  skip_if_not(exists("derive_malawi_binary"))
  d <- data.frame(sel = c(50, 84.6, 120, NA), iod = c(50, 100, 400, NA))
  d <- derive_malawi_binary(d, "sel", "sel_def", 84.6)
  d <- derive_malawi_binary(d, "iod", "iod_def", 100)
  expect_equal(d$sel_def, c(1L, 0L, 0L, NA_integer_))
  expect_equal(d$iod_def, c(1L, 0L, 0L, NA_integer_))
  # an existing column is left alone
  d2 <- data.frame(sel = c(50, 120), sel_def = c(9L, 9L))
  expect_equal(derive_malawi_binary(d2, "sel", "sel_def", 84.6)$sel_def, c(9L, 9L))
})

test_that("the Malawi merged dataset yields the survey's own prevalences", {
  skip_if_not(exists("derive_malawi_binary"))
  p <- here::here("data", "IPD", "Malawi", "Malawi_merged_dataset.rds")
  skip_if_not(file.exists(p), "Malawi merged dataset absent")
  m <- readRDS(p)
  m <- derive_malawi_binary(m, "sel", "sel_def", 84.6)
  m <- derive_malawi_binary(m, "iod", "iod_def", 100)
  w <- m[m$population == "women", ]; c <- m[m$population == "preschool children", ]
  # Phiri et al. 2019 (same survey): about 62% of women and >80% of preschool
  # children below 84.6 ug/L; MNS report: iodine insufficiency in women ~10-15%
  expect_gt(mean(w$sel_def, na.rm = TRUE), 0.55); expect_lt(mean(w$sel_def, na.rm = TRUE), 0.70)
  expect_gt(mean(c$sel_def, na.rm = TRUE), 0.80)
  expect_gt(mean(w$iod_def, na.rm = TRUE), 0.08); expect_lt(mean(w$iod_def, na.rm = TRUE), 0.20)
  expect_gt(sum(is.finite(w$sel)), 750); expect_gt(sum(is.finite(w$iod)), 750)
})
