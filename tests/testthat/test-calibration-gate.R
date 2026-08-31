# =============================================================================
# WS-D. The calibration gate.
#
# WS2c measured the un-anchored national aggregate against the survey's own
# national estimate across 24 cells: mean absolute gap 9.60 pp, worst case 77.57
# pp where the ridge predicts 89.6 percent against a measured 12.0 percent. That
# model reached committed result tables with nothing to stop it.
# =============================================================================

test_that("the gate passes a well-calibrated cell and fails a wild one", {
  skip_if_not(exists("calibration_gate"), "calibration_gate not loaded")
  good <- calibration_gate(pred = rep(0.20, 30), national_prev = 0.21,
                           national_se = 0.01)
  expect_equal(good$status, "passed")
  expect_lt(good$abs_gap_pp, good$threshold_pp)

  bad <- calibration_gate(pred = rep(0.896, 14), national_prev = 0.120,
                          national_se = 0.01)
  expect_equal(bad$status, "calibration_failed")
  expect_gt(bad$abs_gap_pp, 70)
})

test_that("the threshold widens with an imprecise survey estimate", {
  skip_if_not(exists("calibration_gate"), "calibration_gate not loaded")
  # A 12 pp gap against a precise survey estimate fails on the absolute floor.
  tight <- calibration_gate(pred = rep(0.32, 20), national_prev = 0.20,
                            national_se = 0.005)
  expect_equal(tight$status, "calibration_failed")
  # The same gap against an estimate whose own CI half-width is 8 pp does not,
  # because the model is being asked to reproduce a number the survey does not
  # itself pin down. Punishing it there would charge the model for survey noise.
  loose <- calibration_gate(pred = rep(0.32, 20), national_prev = 0.20,
                            national_se = 0.08 / 1.96)
  expect_gt(loose$threshold_pp, tight$threshold_pp)
  expect_equal(loose$status, "passed")
})

test_that("a missing standard error falls back to the absolute floor", {
  skip_if_not(exists("calibration_gate"), "calibration_gate not loaded")
  g <- calibration_gate(pred = rep(0.35, 10), national_prev = 0.20,
                        national_se = NA_real_)
  expect_equal(g$threshold_pp, 10)
  expect_equal(g$status, "calibration_failed")
})

test_that("population weights are honoured", {
  skip_if_not(exists("calibration_gate"), "calibration_gate not loaded")
  # Two districts, one tiny and extreme. Unweighted the aggregate is 0.50 and
  # fails; weighted by population it is close to 0.10 and passes. A gate that
  # ignored the weights would reject a correct model.
  p <- c(0.10, 0.90); w <- c(1e6, 1e3)
  expect_equal(calibration_gate(p, pop = w, national_prev = 0.10)$status, "passed")
  expect_equal(calibration_gate(p, pop = NULL, national_prev = 0.10)$status,
               "calibration_failed")
})

test_that("apply_calibration_gate removes exactly the failed cells", {
  skip_if_not(exists("apply_calibration_gate"), "apply_calibration_gate not loaded")
  tbl <- data.frame(country = c("Ghana", "SierraLeone", "Ghana"),
                    outcome = c("child_iron", "child_vitA", "women_iron"),
                    value = 1:3, stringsAsFactors = FALSE)
  gate <- data.frame(country = c("Ghana", "Sierra Leone"),
                     outcome = c("child_iron", "child_vitA"),
                     status = c("passed", "calibration_failed"),
                     stringsAsFactors = FALSE)
  out <- apply_calibration_gate(tbl, gate)
  expect_equal(nrow(out), 2L)
  # Country labels differ across tables ("SierraLeone" against "Sierra Leone"),
  # so the match must be on a squashed key or the gate silently passes it.
  expect_false(any(out$outcome == "child_vitA"))
  expect_equal(attr(out, "gate_excluded"), "Sierra Leone child_vitA")
})

test_that("the retroactive report exists and fails Sierra Leone child vitamin A", {
  p <- here::here("results", "tables", "calibration_gate_report.csv")
  skip_if_not(file.exists(p),
              "calibration_gate_report.csv absent; run wsd_calibration_gate.R")
  r <- utils::read.csv(p, stringsAsFactors = FALSE)
  expect_true(all(c("country", "outcome", "status", "abs_gap_pp") %in% names(r)))
  sl <- r[grepl("[Ss]ierra", r$country) & r$outcome == "child_vitA", ]
  expect_equal(nrow(sl), 1L)
  expect_equal(sl$status[1], "calibration_failed")
  expect_gt(sl$abs_gap_pp[1], 70)
})
