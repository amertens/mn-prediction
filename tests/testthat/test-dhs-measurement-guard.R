# =============================================================================
# WS-C4. The DHS standard-recode measurement bands.
#
# Every leak guard written before this one matches on ANALYTE NAMES. DHS
# recodes are named by questionnaire position, so none of them fires. These
# tests exist because admitting Malawi to the individual-level arms puts
# several hundred DHS-coded columns under the `gw_` prefix for the first time,
# and `gw_hw53` is a child's haemoglobin reading that no prior pattern catches.
# =============================================================================

test_that("DHS haemoglobin and anaemia codes are not questionnaire columns", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  hb <- c("gw_hw53",   # child haemoglobin, KR recode
          "gw_hw57",   # child anaemia level
          "gw_v453",   # women's haemoglobin, IR recode
          "gw_v457", "gw_v458",   # women's anaemia level, measurement result
          "gw_ha53", "gw_ha57",   # women in the PR person recode
          "gw_hc53", "gw_hc57",   # children in the PR person recode
          "gw_hb53")              # men
  cl <- biomarker_column_class(hb)
  expect_true(all(cl != "questionnaire"))
  expect_true(all(is_biomarker_column(hb)))
  # and they must be excluded from the questionnaire arm specifically
  expect_true(all(!allowed_under_arm(hb, "questionnaire")))
})

test_that("the malaria blood test is classed as a draw, not a field measure", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  expect_equal(biomarker_column_class("gw_hml35"), "blood_draw")
  expect_equal(biomarker_column_class("gw_hml32"), "blood_draw")
  # A draw is excluded from BOTH arms; a field measure only from the first.
  expect_false(allowed_under_arm("gw_hml35", "questionnaire_hb"))
  expect_true(allowed_under_arm("gw_hw53", "questionnaire_hb"))
})

test_that("bednet counts are questionnaire items, not malaria tests", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  # hml1 and hml2 are numbers of mosquito nets. A guard written as ^hml[0-9]+$
  # would take them, and they are exactly the kind of household item the
  # questionnaire arm is meant to contain.
  expect_equal(biomarker_column_class(c("gw_hml1", "gw_hml2")),
               c("questionnaire", "questionnaire"))
})

test_that("anthropometry is a field measurement, not a questionnaire response", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  anthro <- c("gw_hw1", "gw_hw2", "gw_hw3", "gw_hw70", "gw_hw71", "gw_hw72",
              "gw_ha2", "gw_ha3", "gw_ha40",
              # the IR block, which the first version of this test wrongly
              # listed as questionnaire items below
              "gw_v437", "gw_v438", "gw_v445")
  expect_true(all(biomarker_column_class(anthro) == "hb_field"))
  # Height and weight come from the fieldworker's visit, so they belong with
  # field haemoglobin rather than with the questionnaire.
  expect_true(all(!allowed_under_arm(anthro, "questionnaire")))
  expect_true(all(allowed_under_arm(anthro, "questionnaire_hb")))
})

test_that("ordinary DHS questionnaire codes survive the guard", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  # Education, water, toilet, floor, assets, vaccination, breastfeeding.
  ok <- c("gw_v106", "gw_v113", "gw_v116", "gw_v127", "gw_v161",
          "gw_hv201", "gw_hv205", "gw_hv206", "gw_hv213", "gw_hv270",
          "gw_h3", "gw_h33", "gw_h43", "gw_m45_1", "gw_v404", "gw_v459")
  expect_true(all(biomarker_column_class(ok) == "questionnaire"))
  expect_true(all(allowed_under_arm(ok, "questionnaire")))
})

test_that("the guard is scoped to the survey prefix", {
  skip_if_not(.mnp_load_biomarker_guard(), "classifier not loaded")
  # Admin-2 aggregates derived from OTHER DHS rounds are proxy covariates:
  # available to a country that never ran a micronutrient survey, and they must
  # not be swept up by a pattern aimed at the concurrent survey.
  expect_equal(dhs_measurement_class(c("dhs_hw53", "hw53", "map_hb_c")),
               c("", "", ""))
  expect_equal(dhs_measurement_class("gw_hw53"), "hb_field")
  # v459 sits inside the numeric band but asks whether the household has a
  # bednet. A band guard must not swallow the question next to the measurement.
  expect_equal(dhs_measurement_class("gw_v459"), "")
})
