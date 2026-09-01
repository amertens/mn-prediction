# =============================================================================
# Per-country Admin-2 unit counts.
#
# The pair-key fan is silent in code and loud in the row count: the tenth
# instance surfaced as 90 districts where the prediction table held 87. The
# assertion is one-sided, because a table may legitimately hold fewer units than
# the survey supports but never more. R/unit_counts.R explains the asymmetry.
#
# These tests read the artefacts written by scripts/accuracy_impact/ws7a_guards.R
# rather than the targets store, so the suite runs in seconds and does not need
# a built store present.
# =============================================================================

test_that("the reference counts were derived, and are self-consistent", {
  ref <- mnp_artefact("admin2_unit_reference.csv")
  expect_true(all(c("country", "analytic_units", "n_outcomes") %in% names(ref)))
  expect_gte(nrow(ref), 4L)
  expect_true(all(ref$analytic_units > 0))
  # Every country must contribute at least the four shared outcomes, otherwise
  # the reference was built from a partially populated store and the maxima
  # below would be too low to catch anything.
  expect_true(all(ref$n_outcomes >= 4L))
})

test_that("no consumed table exceeds its country's analytic unit count", {
  chk <- mnp_artefact("admin2_unit_check.csv")
  over <- chk[chk$status == "OVER", , drop = FALSE]
  if (nrow(over)) {
    msg <- paste(sprintf("  %s / %s: %d rows against a reference of %d",
                         over$label, over$country, over$observed_max, over$reference),
                 collapse = "\n")
    fail(paste0("Admin-2 unit over-count, the pair-key fan signature:\n", msg))
  }
  expect_equal(nrow(over), 0L)
})

test_that("every consumed table resolves against the reference", {
  chk <- mnp_artefact("admin2_unit_check.csv")
  # A country present in a result table but absent from the reference means the
  # check silently passed that country. That is the failure mode this catches.
  expect_equal(sum(chk$status == "no reference"), 0L)
  # All four countries appear for the anchoring table specifically, which is the
  # one Section 4 rests on.
  arms <- chk[chk$label == "admin1_arms.csv", ]
  expect_gte(nrow(arms), 4L)
})

test_that("Malawi is the country where a name-only join can fan", {
  # The premise of the whole guard. If GADM ever stops having duplicate Malawi
  # district names, the baseline's grandfathered joins become safe and this test
  # should be revisited rather than silently continuing to protect nothing.
  ref <- mnp_artefact("admin2_unit_reference.csv")
  mw <- ref[ref$country == "Malawi", ]
  skip_if(nrow(mw) == 0, "Malawi absent from the reference")
  expect_gt(mw$analytic_units[1], 80)
})
