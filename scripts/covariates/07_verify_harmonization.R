# =============================================================================
# scripts/covariates/07_verify_harmonization.R
#
# Adversarial verification of the covariate harmonisation. Every check below
# exists because it is a way the harmoniser could be WRONG while still looking
# right -- a renamed column that quietly took its values from the wrong source,
# a unit conversion applied twice, a key join that matched the wrong district.
#
# Run after any change to the rules, the engine, or a source:
#   Rscript scripts/covariates/07_verify_harmonization.R
#
# Exit status is non-zero if any check FAILS, so this can gate a rebuild.
# Output: results/tables/harmonization_verification.csv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
source(here("R", "covariates", "canonicalize.R"))

DIR <- here("data", "covariates", "harmonized")
RAW <- here("data", "covariates", "country")
checks <- list()
note <- function(id, what, status, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(check = id, what = what,
                                               status = status, detail = detail,
                                               stringsAsFactors = FALSE)
  cat(sprintf("[%-4s] %-46s %s\n", status, id, detail))
}

pooled <- suppressMessages(readr::read_csv(file.path(DIR, "predictors_admin2_harmonized.csv"),
                                           show_col_types = FALSE)) |> as.data.frame()
map    <- suppressMessages(readr::read_csv(file.path(DIR, "column_map.csv"), show_col_types = FALSE))
covg   <- suppressMessages(readr::read_csv(file.path(DIR, "coverage.csv"), show_col_types = FALSE))
dict   <- suppressMessages(readr::read_csv(file.path(DIR, "data_dictionary.csv"), show_col_types = FALSE))
excl   <- cov_load_exclusions(); units <- cov_load_units()
vars   <- setdiff(names(pooled), c("country", "Admin1", "Admin2"))
countries <- sort(unique(pooled$country))

cat("\n===== STRUCTURAL =====\n")

# 1. No raw column left unclassified. An unmatched name is a covariate that
#    silently disappeared from the model.
n_unmatched <- sum(map$action == "unmatched")
note("1-no-unmatched", "every raw column is classified",
     if (n_unmatched == 0) "PASS" else "FAIL", paste(n_unmatched, "unmatched"))

# 2. One canonical name must never be fed by two DIFFERENT source families in
#    the same country -- that would silently blend two measurements.
fam_of <- function(x) sub("_.*$", "", cov_base_name(x))
collide <- map %>% filter(action == "keep", !is.na(canonical)) %>%
  mutate(fam = fam_of(variable)) %>%
  group_by(country, canonical) %>% summarise(n_fam = dplyr::n_distinct(fam), .groups = "drop") %>%
  filter(n_fam > 1)
note("2-no-family-collision", "canonical fed by one family per country",
     if (!nrow(collide)) "PASS" else "FAIL",
     if (nrow(collide)) paste(utils::head(collide$canonical, 3), collapse = ", ") else "")

# 3. collapse='latest' must be unambiguous: if several raw columns tie on the
#    highest year, the pick is arbitrary and therefore country-dependent.
ties <- map %>% filter(action == "keep", collapse == "latest") %>%
  group_by(country, canonical) %>%
  summarise(n_at_max = if (all(is.na(year))) 1L else
              sum(year == suppressWarnings(max(year, na.rm = TRUE)), na.rm = TRUE),
            n_tot = dplyr::n(), .groups = "drop") %>%
  filter(n_tot > 1, n_at_max > 1)
note("3-latest-unambiguous", "collapse='latest' has a unique winner",
     if (!nrow(ties)) "PASS" else "WARN",
     if (nrow(ties)) paste(nrow(ties), "ambiguous, e.g.", ties$canonical[1]) else "")

# 4. collapse='none' must be genuinely 1:1.
dupnone <- map %>% filter(action == "keep", collapse == "none") %>%
  count(country, canonical) %>% filter(n > 1)
note("4-none-is-1to1", "collapse='none' maps 1:1",
     if (!nrow(dupnone)) "PASS" else "FAIL",
     if (nrow(dupnone)) paste(utils::head(dupnone$canonical, 3), collapse = ", ") else "")

# 5. Excluded variables must not appear in the shipped data.
in_out <- vars[Reduce(`|`, lapply(excl$canonical_regex,
                                  function(p) grepl(p, vars, perl = TRUE)),
                      init = rep(FALSE, length(vars)))]
note("5-exclusions-honoured", "excluded variables absent from output",
     if (!length(in_out)) "PASS" else "FAIL", paste(in_out, collapse = ", "))

# 6. Every shipped column must be documented.
undoc <- setdiff(vars, dict$canonical)
nofam <- dict$canonical[is.na(dict$family)]
nodom <- dict$canonical[is.na(dict$domain)]
note("6-documented", "every predictor has a dictionary row",
     if (!length(undoc)) "PASS" else "FAIL", paste(length(undoc), "undocumented"))
note("6b-domain-assigned", "every predictor has family + domain",
     if (!length(nofam) && !length(nodom)) "PASS" else "FAIL",
     sprintf("%d without family, %d without domain", length(nofam), length(nodom)))

# 7. The pooled table must be one row per (country, Admin1, Admin2).
key <- paste(pooled$country, pooled$Admin1, pooled$Admin2, sep = "")
note("7-unique-keys", "one row per country x admin-2",
     if (!anyDuplicated(key)) "PASS" else "FAIL",
     paste(sum(duplicated(key)), "duplicate keys"))

# 8. No all-NA and no zero-variance columns.
allna <- vars[vapply(vars, function(v) all(is.na(pooled[[v]])), logical(1))]
novar <- vars[vapply(vars, function(v) {
  x <- pooled[[v]]; sum(is.finite(x)) > 0 && stats::sd(x, na.rm = TRUE) == 0 }, logical(1))]
note("8-no-dead-columns", "no all-NA or zero-variance predictors",
     if (!length(allna) && !length(novar)) "PASS" else "WARN",
     sprintf("%d all-NA, %d constant", length(allna), length(novar)))

# 9. Exact-duplicate columns: two canonical names carrying identical values are
#    the redundancy this exercise was meant to remove.
sig <- vapply(vars, function(v) paste(utils::head(round(pooled[[v]], 10), 200), collapse = "|"),
              character(1))
dups <- names(sig)[duplicated(sig)]
note("9-no-duplicate-columns", "no two predictors are value-identical",
     if (!length(dups)) "PASS" else "WARN",
     if (length(dups)) paste(utils::head(dups, 4), collapse = ", ") else "")

cat("\n===== VALUES: DID THE HARMONISER MOVE THE RIGHT NUMBERS? =====\n")

# 10. Re-derive a sample of canonical values straight from the per-country raw
#     table and compare. This is the check that catches a rule pointing at the
#     wrong source column, or a unit conversion applied twice.
recheck <- list()
for (ctry in countries) {
  f <- file.path(RAW, sprintf("%s_predictors_admin2_raw.csv", ctry))
  if (!file.exists(f)) next
  raw <- suppressMessages(readr::read_csv(f, show_col_types = FALSE)) |> as.data.frame()
  mp  <- map[map$country == ctry & map$action == "keep" & !is.na(map$canonical), , drop = FALSE]
  po  <- pooled[pooled$country == ctry, , drop = FALSE]
  idx <- match(paste(po$Admin1, po$Admin2, sep = ""),
               paste(trimws(raw$Admin1), trimws(raw$Admin2), sep = ""))
  uf  <- cov_unit_factors(intersect(vars, mp$canonical), ctry, units)
  for (cn in intersect(vars, unique(mp$canonical))) {
    grp <- mp[mp$canonical == cn, , drop = FALSE]
    src <- if (nrow(grp) == 1L) grp$variable[1] else
      grp$variable[order(grp$year, decreasing = TRUE, na.last = TRUE)][1]
    if (!src %in% names(raw)) next
    x <- suppressWarnings(as.numeric(raw[[src]]))[idx]
    u <- uf[uf$canonical == cn, , drop = FALSE]
    if (identical(grp$collapse[1], "mean") && nrow(grp) > 1L) {
      M <- vapply(grp$variable, function(v) suppressWarnings(as.numeric(raw[[v]]))[idx],
                  numeric(length(idx)))
      x <- rowMeans(M, na.rm = TRUE); x[is.nan(x)] <- NA_real_
    }
    expect <- x * (u$multiply %||% 1) + (u$add %||% 0)
    # The harmoniser clamps physically impossible negatives to zero for the
    # four count/area variables below (canonicalize.R). GHSL population comes
    # back from raster resampling with small negatives -- Ghana min -7.29 --
    # and left alone they reach models as real data. The verifier has to apply
    # the same clamp or it reports a 7.29 discrepancy on a deliberate fix, and
    # a permanently red check trains everyone to ignore it.
    if (cn %in% c("ghs_pop", "built_surface", "built_surface_nres",
                  "grassland_frac")) {
      expect[is.finite(expect) & expect < 0] <- 0
    }
    got <- po[[cn]]
    d <- suppressWarnings(max(abs(expect - got), na.rm = TRUE))
    recheck[[length(recheck) + 1L]] <- data.frame(
      country = ctry, canonical = cn, max_abs_diff = if (is.finite(d)) d else NA_real_)
  }
}
rc <- dplyr::bind_rows(recheck)
bad <- rc %>% filter(is.finite(max_abs_diff), max_abs_diff > 1e-9)
note("10-values-reproduce", "canonical values re-derive from raw + declared units",
     if (!nrow(bad)) "PASS" else "FAIL",
     sprintf("%d of %d (country, variable) pairs re-derive exactly",
             nrow(rc) - nrow(bad), nrow(rc)))
if (nrow(bad)) print(utils::head(as.data.frame(bad), 8))

# 11. Unit conversions must have actually landed: post-conversion, the countries
#     that received a conversion must agree with those that did not.
lst <- grep("^lst_night_(m[0-9]{2}|annual_(mean|min|max))_", vars, value = TRUE)
lst_ok <- if (!length(lst)) NA else {
  r <- range(unlist(pooled[, lst]), na.rm = TRUE); r[1] > -10 && r[2] < 45 }
note("11a-lst-celsius", "night LST is in a plausible Celsius range",
     if (isTRUE(lst_ok)) "PASS" else "FAIL",
     if (length(lst)) sprintf("range %.1f to %.1f C over %d columns",
                              min(unlist(pooled[, lst]), na.rm = TRUE),
                              max(unlist(pooled[, lst]), na.rm = TRUE), length(lst)) else "absent")

ndvi <- grep("^ndvi_y", vars, value = TRUE)
ndvi_ok <- length(ndvi) > 0 && all(unlist(pooled[, ndvi]) >= -1, na.rm = TRUE) &&
           all(unlist(pooled[, ndvi]) <= 1, na.rm = TRUE)
note("11b-ndvi-scaled", "NDVI lies in [-1, 1]",
     if (ndvi_ok) "PASS" else "FAIL",
     sprintf("range %.3f to %.3f over %d years",
             min(unlist(pooled[, ndvi]), na.rm = TRUE),
             max(unlist(pooled[, ndvi]), na.rm = TRUE), length(ndvi)))

frac <- grep("^lcover_.*_frac_", vars, value = TRUE)
frac_ok <- !length(frac) || (all(unlist(pooled[, frac]) >= 0, na.rm = TRUE) &&
                             all(unlist(pooled[, frac]) <= 100, na.rm = TRUE))
note("11c-cover-fractions", "cover fractions lie in [0, 100]",
     if (frac_ok) "PASS" else "FAIL",
     sprintf("range %.2f to %.2f", min(unlist(pooled[, frac]), na.rm = TRUE),
             max(unlist(pooled[, frac]), na.rm = TRUE)))

aef <- grep("^aef_A", vars, value = TRUE)
nrm <- sqrt(rowSums(as.matrix(pooled[, aef])^2, na.rm = TRUE))
note("11d-aef-unit-norm", "AlphaEarth vectors are unit length",
     if (length(aef) == 64 && max(abs(nrm - 1)) < 1e-8) "PASS" else "FAIL",
     sprintf("%d dims, max |norm-1| = %.2e", length(aef), max(abs(nrm - 1))))

# 12. Same measurement, same scale across countries. A conversion that was
#     declared but did not work shows up as a spread ratio that is still huge.
spread <- vapply(vars, function(v) {
  iq <- tapply(pooled[[v]], pooled$country, function(x) stats::IQR(x, na.rm = TRUE))
  iq <- iq[is.finite(iq) & iq > 0]
  if (length(iq) < 2) NA_real_ else max(iq) / min(iq)
}, numeric(1))
declared <- vars[vapply(vars, function(v)
  any(vapply(seq_len(nrow(units)), function(i)
    grepl(units$canonical_regex[i], v, perl = TRUE), logical(1))), logical(1))]
unresolved <- declared[is.finite(spread[declared]) & spread[declared] > 20]
note("12-conversions-worked", "no declared conversion leaves a >20x spread gap",
     if (!length(unresolved)) "PASS" else "FAIL", paste(unresolved, collapse = ", "))

cat("\n===== COVERAGE =====\n")

# 13. Every country contributes every shipped predictor (that is what "shared"
#     means) and the coverage table agrees with the data.
miss <- vapply(vars, function(v)
  sum(tapply(pooled[[v]], pooled$country, function(x) all(is.na(x)))), numeric(1))
note("13-shared-really-shared", "no shipped predictor is all-NA in a country",
     if (all(miss == 0)) "PASS" else "FAIL",
     paste(sum(miss > 0), "predictors empty in >=1 country"))

agree <- setdiff(covg$canonical[covg$in_all], vars)
note("13b-coverage-agrees", "coverage.csv matches the shipped columns",
     if (!length(agree)) "PASS" else "FAIL", paste(length(agree), "listed but not shipped"))

res <- dplyr::bind_rows(checks)
dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(res, here("results", "tables", "harmonization_verification.csv"))

cat(sprintf("\n%d checks: %d PASS, %d WARN, %d FAIL\n", nrow(res),
            sum(res$status == "PASS"), sum(res$status == "WARN"), sum(res$status == "FAIL")))
if (any(res$status == "FAIL")) {
  cat("\nFAILING:\n"); print(as.data.frame(res[res$status == "FAIL", ]), row.names = FALSE)
  quit(status = 1)
}
