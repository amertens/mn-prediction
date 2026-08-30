# =============================================================================
# scripts/covariates/11_vocabulary_compatibility.R
#
# Prove that every consumer of the admin-2 covariate table works under BOTH
# vocabularies, and in particular that none of them silently selects ZERO
# covariates.
#
# WHY
# ---
# The harmonised switch was launched three times before this test existed. Each
# failure was the same shape: a consumer selected covariates with a hardcoded
# `grep("^gee_", ...)`, which returns nothing when names carry no prefix. Some
# of those raised an error; the dangerous ones did not -- the individual-level
# merge in _targets.R guards its merge with `if (length(gee_cols) > 0)`, so it
# would have dropped every Earth-observation covariate out of the person-level
# model without a word.
#
# A covariate selector returning zero is never correct here. This test asserts
# that directly, for each selector, in each vocabulary.
#
#   Rscript scripts/covariates/11_vocabulary_compatibility.R
# Exit status is non-zero on failure, so it can gate a launch.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
targets::tar_source("R/")

STORE <- here("_targets_full")
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)
fails <- character()
check <- function(label, ok, detail = "") {
  cat(sprintf("[%-4s] %-52s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) fails <<- c(fails, label)
}

# Two tables for the same country, one in each vocabulary.
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
harm <- H[H$country == "Ghana", setdiff(names(H), "country")]

# The legacy table is SYNTHESISED from the harmonised one by re-prefixing,
# rather than read from the targets store. Reading the store made this test
# unreliable: the store is _targets_full (see _targets.yaml), so an aborted
# harmonised run had already overwritten gee_admin2_ghana, and the "legacy"
# fixture was silently harmonised -- the test then only ever exercised one
# vocabulary. What is under test here is the SELECTOR's handling of two name
# shapes, and a synthetic fixture tests that deterministically.
legacy <- harm
covs <- setdiff(names(legacy), c("Admin1", "Admin2"))
names(legacy)[match(covs, names(legacy))] <- paste0("gee_", covs)
# A few real legacy name shapes, so the soil/crop patterns are exercised as they
# actually appear rather than only as re-prefixed canonical names.
legacy$gee_soilzinc_mean_0_20 <- harm[[grep("^soil_zinc_mean", names(harm), value = TRUE)[1]]]
legacy$gee_esa_worldcereal_32121_tc_maize_main_annual_mean <- runif(nrow(harm))

cat("\n=== vocabulary detection ===\n")
check("legacy table detected as legacy", area_vocab_of(legacy) == "legacy", area_vocab_of(legacy))
check("harmonised table detected as harmonised", area_vocab_of(harm) == "harmonized", area_vocab_of(harm))

cat("\n=== selector returns a non-empty set in BOTH vocabularies ===\n")
for (nm in c("legacy", "harmonized")) {
  d <- if (nm == "legacy") legacy else harm
  a <- area_covariate_cols(d)
  u <- area_covariate_cols(d, regime = "universal")
  s <- area_soil_cols(d)
  cr <- area_crop_cols(d)
  check(sprintf("%s: area_covariate_cols(all)", nm), length(a) > 20, paste(length(a), "cols"))
  check(sprintf("%s: area_covariate_cols(universal)", nm), length(u) > 20, paste(length(u), "cols"))
  check(sprintf("%s: area_soil_cols", nm), length(s) > 5, paste(length(s), "cols"))
  check(sprintf("%s: universal excludes survey sources", nm),
        !any(grepl("^(dhs_|espen_|MAP_|ihme_)", u)), "")
  check(sprintf("%s: no key column selected as a covariate", nm),
        !any(c("Admin1", "Admin2", "svy_prev", "n_svy") %in% a), "")
  if (nm == "harmonized")
    check("harmonized: crop columns found", length(cr) > 0, paste(length(cr), "cols"))
}

cat("\n=== name-level classifier (pooled individual design matrix) ===\n")
check("legacy names classified as area covariates",
      all(is_area_covariate_name(head(grep("^gee_", names(legacy), value = TRUE), 20))), "")
check("canonical names classified as area covariates",
      all(is_area_covariate_name(head(area_covariate_cols(harm), 20))), "")
check("survey items NOT classified as area covariates",
      !any(is_area_covariate_name(c("gw_cSex", "gw_hWealthquintile", "Admin2"))), "")

cat("\n=== downstream consumers accept both vocabularies ===\n")
cfgs <- get_country_configs(); cc <- cfgs[["Ghana"]]; oc <- cc$outcomes$child_iron
svy <- rd("svy_admin2_ghana_child_iron")
od  <- rd("outcome_data_ghana_child_iron")

for (nm in c("legacy", "harmonized")) {
  d <- if (nm == "legacy") legacy else harm
  # area recipe frame: the target that failed on the third launch
  fr <- tryCatch(ar_build_frame(od, svy, d, cc, oc, mode = "universal"),
                 error = function(e) structure(list(), err = conditionMessage(e)))
  ok <- is.data.frame(fr) && nrow(fr) > 0 &&
    length(area_covariate_cols(fr, regime = "universal")) > 10
  check(sprintf("%s: ar_build_frame(universal)", nm), ok,
        if (is.data.frame(fr)) sprintf("%d rows", nrow(fr)) else attr(fr, "err"))

  # feature engineering
  fe <- tryCatch(engineer_admin2_features(d), error = function(e) NULL)
  check(sprintf("%s: engineer_admin2_features", nm),
        !is.null(fe) && ncol(fe) > 1, if (is.null(fe)) "NULL" else paste(ncol(fe) - 1, "features"))

  # weighted area ensemble
  w <- tryCatch(fit_area_weighted_sl(svy, d, weighted = TRUE), error = function(e) NULL)
  check(sprintf("%s: fit_area_weighted_sl", nm),
        !is.null(w) && nrow(w$area_preds) > 0,
        if (is.null(w)) "NULL" else sprintf("%d vars", w$n_vars))
}

cat(sprintf("\n%d checks failed\n", length(fails)))
if (length(fails)) { print(fails); quit(status = 1) }
cat("PASS: both vocabularies work everywhere tested\n")
