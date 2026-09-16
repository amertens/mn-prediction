# =============================================================================
# The fit-time predictor policy (TP-01 / TA-01, 2026-09-15).
#
# Every protocol script takes its column list through drop_near_outcome_v2().
# Three things are enforced there and stamped in the metadata:
#   tiers          open / survey_public / survey_dhs from predictor_tiers.csv;
#                  V2_PREDICTOR_TIERS selects the tiers a run may use
#   national       subnational == FALSE columns are dropped unless
#                  V2_KEEP_NATIONAL=1 (they cannot rank districts)
#   alignment      every column has a temporal_alignment.csv rule; the stamp
#                  writes year_used and year_offset_max_abs
# =============================================================================

source(here::here("R", "protocol_v2.R"))

.fake_meta <- function() data.frame(
  column = c("dhs_c_stunting", "hces_food_share", "mics_heat_vbcg_sy", "aef_A00", "gfdx_wheat_mandatory_sy", "soil_zinc_mean_0_20"),
  domain = c("Adult nutrition", "Household diet and consumption (HCES)", "Immunisation", "Satellite embedding", "Food fortification and supplementation", "Soil characteristics"),
  source = c("DHS", "HCES microdata (IHS4, IHS 2015/16, SLIHS 2018, GLSS7)",
             "WHO Health Inequality Data Repository (HEAT) - UNICEF MICS by subnational region", "AlphaEarth (GEE)",
             "GFDx (Global Fortification Data Exchange) programme fields", "SoilGrids / iSDA"),
  subnational = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE)

test_that("tiers are assigned from the source by the rules file, first match wins", {
  m <- .fake_meta()
  expect_equal(assign_tier_v2(m), c("survey_dhs", "survey_public", "survey_public", "open", "open", "open"))
  expect_error(assign_tier_v2(data.frame(column = "x", source = NA_character_)), NA)   # NA source falls to the catch-all
})

test_that("the default policy keeps every tier and drops the national constants", {
  m <- .fake_meta()
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "", V2_KEEP_NATIONAL = ""), {
    expect_equal(sort(suppressMessages(drop_near_outcome_v2(m$column, m))),
                 sort(setdiff(m$column, "gfdx_wheat_mandatory_sy")))
  })
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "", V2_KEEP_NATIONAL = "1"), {
    expect_setequal(suppressMessages(drop_near_outcome_v2(m$column, m)), m$column)
  })
})

test_that("V2_PREDICTOR_TIERS selects tiers; an unknown tier is an error", {
  m <- .fake_meta()
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "open"), {
    expect_setequal(suppressMessages(drop_near_outcome_v2(m$column, m)), c("aef_A00", "soil_zinc_mean_0_20"))
  })
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "open,survey_public"), {
    expect_setequal(suppressMessages(drop_near_outcome_v2(m$column, m)),
                    c("aef_A00", "soil_zinc_mean_0_20", "hces_food_share", "mics_heat_vbcg_sy"))
  })
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "open,dhs"), expect_error(predictor_tiers_v2(), "unknown tier"))
  # a stamped `tier` column is used as-is
  m$tier <- c("open", "open", "open", "open", "open", "open")
  withr::with_envvar(c(V2_PREDICTOR_TIERS = "open"), {
    expect_setequal(suppressMessages(drop_near_outcome_v2(m$column, m)), setdiff(m$column, "gfdx_wheat_mandatory_sy"))
  })
})

test_that("domain labels in the live metadata are unique on their first 12 characters (PC column keys)", {
  p <- here::here("data", "covariates", "harmonized", "predictors_admin2_shared_metadata.csv")
  skip_if_not(file.exists(p), "shared metadata absent")
  M <- read.csv(p, stringsAsFactors = FALSE)
  dm <- sort(unique(stats::na.omit(M$domain))); pref <- make.names(substr(dm, 1, 12))
  expect_false(anyDuplicated(pref) > 0, info = paste(dm[pref %in% pref[duplicated(pref)]], collapse = " | "))
})

test_that("the live metadata is stamped and every column has a tier and an alignment rule", {
  p <- here::here("data", "covariates", "harmonized", "predictors_admin2_shared_metadata.csv")
  skip_if_not(file.exists(p), "shared metadata absent")
  M <- read.csv(p, stringsAsFactors = FALSE)
  skip_if_not(all(c("tier", "year_used", "year_offset_max_abs", "alignment_rule") %in% names(M)), "metadata not stamped (run scripts/covariates/stamp_predictor_metadata.R)")
  expect_true(all(M$tier %in% V2_TIERS_ALL))
  expect_equal(assign_tier_v2(M), M$tier)
  expect_true(all(nzchar(M$year_used)))
  expect_true(all(nzchar(M$alignment_rule)))
  expect_true(all(is.na(M$year_offset_max_abs) | M$year_offset_max_abs >= 0))
  # every column in the live set matches a rule in the rules file
  R <- read.csv(here::here("metadata", "covariates", "temporal_alignment.csv"), stringsAsFactors = FALSE)
  hit <- Reduce(`|`, lapply(R$column_regex, function(rx) grepl(rx, M$column, perl = TRUE)))
  expect_true(all(hit), info = paste("unmatched:", paste(M$column[!hit], collapse = ", ")))
  # DHS is survey_dhs, HCES and HEAT are survey_public, nothing else is survey-tier
  expect_true(all(M$tier[grepl("^dhs_", M$column)] == "survey_dhs"))
  expect_true(all(M$tier[grepl("^hces_|^mics_", M$column)] == "survey_public"))   # HCES, MICS microdata, MICS via HEAT
  expect_true(all(M$tier[!grepl("^dhs_|^hces_|^mics_", M$column)] == "open"))
})
