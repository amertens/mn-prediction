# =============================================================================
# The assay guard and the leakage report.
#
# Section 7 of docs/SESSION_FINDINGS_FOR_REVIEW.md records three attempts at the
# guard, two of which produced publishable-looking numbers that were leakage.
# WS7a found an eleventh instance, gw_wm_whbc and gw_gchb, two Ghana haemoglobin
# columns that the case-sensitive `Hb` pattern let through. These tests pin both
# the specific columns and the general property.
# =============================================================================

test_that("the guard predicate could be loaded", {
  expect_true(.mnp_load_biomarker_guard())
})

test_that("the guard blocks the columns Section 7 identified", {
  skip_if_not(.mnp_load_biomarker_guard())
  leaks <- c(
    # Section 7.2, the measured top of the correlation ranking
    "gw_cAnemiaYN", "gw_cAnemiaCat", "gw_bs2", "gw_bis", "gw_cID_NoAdj",
    # Section 7.1, the raw assay panel that produced r = 0.973
    "gw_cFER", "gw_cSTFR", "gw_cCRP", "gw_cAGP", "gw_wFER",
    # derived indicators
    "gw_cVitADef", "gw_wFerAdj", "gw_cVAI", "gw_rpb1"
  )
  expect_true(all(is_biomarker_column(leaks)),
              info = paste("unblocked:", paste(leaks[!is_biomarker_column(leaks)],
                                               collapse = ", ")))
})

test_that("the guard blocks the lower-case haemoglobin columns found in WS7a", {
  skip_if_not(.mnp_load_biomarker_guard())
  # Regression pin. gw_wm_whbc is women's haemoglobin in g/dL (measured: n 981,
  # range 7.0-17.0, median 12.9) and gw_gchb is child haemoglobin (n 1159,
  # range 6.2-15.0, median 11.4). Both sat inside the no-blood-draw arm.
  expect_true(is_biomarker_column("gw_wm_whbc"))
  expect_true(is_biomarker_column("gw_gchb"))
  expect_true(is_biomarker_column("gw_bs2_hb"))
  expect_true(is_biomarker_column("gw_bs3_hba1c"))
})

test_that("the guard spares the columns Section 7.4 verified individually", {
  skip_if_not(.mnp_load_biomarker_guard())
  keep <- c(
    # identifiers
    "gw_HHID", "gw_MotherID", "gw_SAMPLE_ID", "gw_WomanID", "gw_hID",
    "gw_indivID", "gw_momID", "gw_old_indivID", "gw_hHRostCaretakerID02",
    # knowledge and belief items
    "gw_wHeardAnemia", "gw_wFFAnemia", "gw_wFFGivesBlood", "gw_wFFImpHealth",
    # exposures
    "gw_wIodSalt", "gw_wFeSuppl", "gw_wFeSupplTime", "gw_wSupplastPregFeFA",
    "gw_wVitASuppl", "gw_wFeRichFood", "gw_wFvoPreferBuy",
    # the lower-case household block the case-sensitivity was chosen to spare
    "gw_hBuyOilTimes", "gw_hBreadType", "gw_hBuyBreadTimes", "gw_hBirdsNum"
  )
  bad <- keep[is_biomarker_column(keep)]
  expect_equal(length(bad), 0L,
               info = paste("over-blocked:", paste(bad, collapse = ", ")))
})

test_that("no cell carries a flagged predictor", {
  sm <- mnp_artefact("leakage_report_summary.csv")
  flagged <- sm[sm$n_flagged > 0, , drop = FALSE]
  if (nrow(flagged)) {
    msg <- paste(sprintf("  %s / %s / %s: %d column(s), max |r| %.3f, top %s",
                         flagged$country, flagged$outcome, flagged$predictor_set,
                         flagged$n_flagged, flagged$max_abs_r, flagged$top_column),
                 collapse = "\n")
    fail(paste0("Predictors above the leakage threshold:\n", msg))
  }
  expect_equal(nrow(flagged), 0L)
})

test_that("the report covers every cell in both predictor sets", {
  sm <- mnp_artefact("leakage_report_summary.csv")
  # A report that silently covered three cells would pass the test above while
  # protecting nothing. Both predictor sets must be present for every cell.
  expect_setequal(unique(sm$predictor_set), c("proxy", "questionnaire"))
  cells <- unique(paste(sm$country, sm$outcome))
  expect_gte(length(cells), 24L)
  expect_equal(nrow(sm), 2L * length(cells))
  # Every flagged decision must rest on a real sample.
  expect_true(all(sm$max_abs_r <= 1))
})
