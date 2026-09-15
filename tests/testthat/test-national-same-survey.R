# =============================================================================
# National track (VMNIS prevalence model): the same-survey rule for anaemia.
#
# LK-02. The country x year panel carries DHS-measured anaemia / haemoglobin
# columns. They are legitimate external predictors of a VMNIS prevalence EXCEPT
# where the VMNIS survey is the DHS round itself (a micronutrient module carried
# by a DHS: Malawi 2015-16, Tanzania 2010, Uganda 2006/2011/2016, Cambodia
# 2014 ...), because then the anaemia and the deficiency were measured on the
# same people. VMNIS records that linkage in its Surveymethodology text
# ("Survey carried out jointly with the 2015-16 Malawi Demographic and Health
# Survey"). The rule: keep the columns, blank the cells of the linked
# country-years, let the nearest-year carry fill them from another round.
# =============================================================================

source(here::here("R", "national_covariates.R"))

test_that("DHS-linked VMNIS surveys are recognised from the methodology text", {
  skip_if_not(exists("dhs_linked_vmnis_surveys"))
  nat <- data.frame(
    iso3c = c("MWI", "MWI", "SLE", "GHA", "UGA"),
    Beginyear = c(2015, 2015, 2013, 2017, 2011),
    Endyear   = c(2016, 2016, 2013, 2017, 2011),
    Surveymethodology = c(
      "Survey carried out jointly with the 2015-16 Malawi Demographic and Health Survey (Survey ID 10901)",
      "Survey carried out jointly with the 2015-16 Malawi Demographic and Health Survey (Survey ID 10901)",
      "Two-stage cluster survey, 30 clusters per district, independent of the DHS sampling frame",
      "Survey conducted based on probability sample to produce estimates nationally and at three belts",
      "Micronutrient module of the 2011 UDHS"),
    stringsAsFactors = FALSE)
  d <- dhs_linked_vmnis_surveys(nat, curated = NULL)   # text rule alone
  expect_setequal(d$iso3c, c("MWI", "UGA"))
  expect_equal(d$year_from[d$iso3c == "MWI"], 2015)
  expect_equal(d$year_to[d$iso3c == "MWI"], 2016)
  expect_equal(nrow(d), 2L)          # duplicates collapsed
})

test_that("the curated list adds DHS-carried modules the text does not name", {
  skip_if_not(exists("dhs_linked_vmnis_surveys"))
  # Cambodia 2014's VMNIS text describes the CDHS design without naming it
  nat <- data.frame(iso3c = "KHM", Beginyear = 2014, Endyear = 2014,
                    Surveymethodology = "Survey representative at the national level and for each of the 19 sampling domains",
                    stringsAsFactors = FALSE)
  d <- dhs_linked_vmnis_surveys(nat)
  expect_true(any(d$iso3c == "KHM" & d$year_from == 2014))
  # and the curated table is the project's single source for that knowledge
  expect_true(file.exists(here::here("metadata", "vmnis_dhs_linked_surveys.csv")))
  expect_equal(nrow(dhs_linked_vmnis_surveys(nat, curated = NULL)), 0L)
})

test_that("a survey that merely mentions being independent of the DHS is not linked", {
  skip_if_not(exists("dhs_linked_vmnis_surveys"))
  nat <- data.frame(iso3c = "SLE", Beginyear = 2013, Endyear = 2013,
                    Surveymethodology = "Sample drawn independently of the DHS frame",
                    stringsAsFactors = FALSE)
  expect_equal(nrow(dhs_linked_vmnis_surveys(nat, curated = NULL)), 0L)
  # and Sierra Leone 2013 is not on the curated list either: SLMS 2013 drew its own sample
  d <- dhs_linked_vmnis_surveys(nat)
  expect_false(any(d$iso3c == "SLE"))
})

test_that("null_same_survey_cells blanks only the linked country-years of the anaemia columns", {
  skip_if_not(exists("null_same_survey_cells"))
  M <- matrix(1:12, nrow = 6, dimnames = list(NULL, c("anaemia_w", "wealth")))
  iso <- c("MWI", "MWI", "MWI", "SLE", "SLE", "GHA")
  yr  <- c(2010, 2015, 2016, 2013, 2019, 2017)
  linked <- data.frame(iso3c = "MWI", year_from = 2015, year_to = 2016)
  out <- null_same_survey_cells(M, iso, yr, cols = "anaemia_w", linked = linked)
  expect_true(all(is.na(out[2:3, "anaemia_w"])))
  expect_equal(out[c(1, 4, 5, 6), "anaemia_w"], M[c(1, 4, 5, 6), "anaemia_w"])
  expect_equal(out[, "wealth"], M[, "wealth"])           # other columns untouched
  expect_equal(attr(out, "n_cells_blanked"), 2L)
})

test_that("the anaemia column detector reads the Stata labels, not the codes", {
  skip_if_not(exists("anaemia_status_columns"))
  labs <- c(AN_ANEM_W_ANY = "Women with any anemia", ML_HEMO_C_HL8 = "Children with hemoglobin < 8 g/dL",
            AN_MNIM_W_IRN = "Women who took iron tablets during pregnancy",
            CN_NUTS_C_HA2 = "Children stunted", HB_KNOW = "Heard of anaemia")
  hit <- anaemia_status_columns(labs)
  expect_setequal(hit, c("AN_ANEM_W_ANY", "ML_HEMO_C_HL8", "HB_KNOW"))
})
