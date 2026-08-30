# =============================================================================
# scripts/covariates/04_document_and_qc.R                           [STAGE 4]
#
# Make the harmonised dataset understandable and shareable, and prove that the
# harmonisation actually worked.
#
# Produces (data/covariates/harmonized/):
#   data_dictionary.csv   one row per canonical variable: what it is, where it
#                         came from, its unit, its causal pathway, its licence,
#                         and the raw column it was built from in each country
#   qc_report.csv         per-variable cross-country distribution check
#   qc_summary.md         the human-readable version, including anything the
#                         harmonisation did NOT fix
#
# THE QC TEST THAT MATTERS
# ------------------------
# Two countries can carry the same column name and still be on different scales
# -- that is how the pre-harmonisation data ended up with Tanzanian night-time
# temperature in Kelvin (292) sitting in the same column as Celsius (20)
# elsewhere. A model reads that as a country effect and calls it signal. So
# every canonical variable is re-checked AFTER harmonisation, and a variable
# whose countries still disagree is reported rather than quietly shipped.
#
# The detector is the ratio of country SPREADS, not of country medians: a unit
# error multiplies the whole distribution, whereas geography moves the level and
# leaves spreads comparable. A variable with a declared unit conversion that
# still disagrees is UNRESOLVED -- the conversion did not do its job.
#
#   Rscript scripts/covariates/04_document_and_qc.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
source(here("R", "covariates", "canonicalize.R"))

DIR <- here("data", "covariates", "harmonized")

pooled <- suppressMessages(readr::read_csv(file.path(DIR, "predictors_admin2_harmonized.csv"),
                                           show_col_types = FALSE))
map    <- suppressMessages(readr::read_csv(file.path(DIR, "column_map.csv"), show_col_types = FALSE))
covg   <- suppressMessages(readr::read_csv(file.path(DIR, "coverage.csv"), show_col_types = FALSE))
reg    <- cov_load_registry()

vars <- setdiff(names(pooled), c("country", "Admin1", "Admin2"))

# ── Data dictionary ─────────────────────────────────────────────────────────
# Family is recovered from the canonical stem, then matched to the registry.
canon_family <- function(v) {
  dplyr::case_when(
    grepl("^aef_", v)           ~ "alphaearth",
    grepl("^ndvi_anomaly", v)   ~ "ndvi_mean_anomaly",
    grepl("^ndvi_y", v)         ~ "ndvi",
    grepl("^popdens_", v)       ~ "popdensity",
    grepl("^precip_", v)        ~ "trmm",
    grepl("^tclim_", v)         ~ "terraclimate",
    grepl("^lst_night_", v)     ~ "lst_night",
    grepl("^lcover_", v)        ~ "landcoverlayers",
    grepl("^npp_", v)           ~ "productivity",
    grepl("^lai_", v)           ~ "lai8days",
    grepl("^evi_", v)           ~ "dailyevi",
    grepl("^wapor_", v)         ~ "wapor",
    grepl("^aod_", v)           ~ "aerosoloptical",
    grepl("^soilgrids_", v)     ~ "soilgrids",
    grepl("^soil_", v)          ~ "soil",
    grepl("^spam_", v)          ~ "mapspam",
    grepl("^espen_", v)         ~ "espen",
    grepl("^dhs_", v)           ~ "dhs",
    grepl("^built_surface", v)  ~ "ghsbuilts",
    grepl("^ghs_pop", v)        ~ "ghspop",
    grepl("^ntl_ccnl", v)       ~ "ccnl",
    grepl("^wsf", v)            ~ "wsf",
    grepl("^elevation", v)      ~ "elevation",
    grepl("^access_", v)        ~ "accessibility",
    grepl("^human_modification", v) ~ "globalhumanmodification",
    grepl("^grassland_frac", v) ~ "gpw_grasslands",
    TRUE ~ NA_character_)
}

src_names <- map %>% filter(action == "keep", !is.na(canonical)) %>%
  group_by(canonical) %>%
  summarise(source_columns = paste(sort(unique(variable)), collapse = " | "),
            n_source_columns = dplyr::n(),
            collapse_policy = paste(sort(unique(collapse)), collapse = ","),
            unit_conversion = paste(sort(unique(stats::na.omit(unit_reason))), collapse = " | "),
            .groups = "drop")

dict <- tibble::tibble(canonical = vars, family = canon_family(vars)) %>%
  left_join(reg, by = "family") %>%
  left_join(src_names, by = "canonical") %>%
  left_join(covg %>% select(canonical, n_countries), by = "canonical") %>%
  mutate(across(where(is.character), ~ ifelse(.x == "", NA_character_, .x))) %>%
  select(canonical, family, provider, domain, pathway, canonical_unit, temporal_kind,
         licence, n_countries, n_source_columns, collapse_policy, unit_conversion,
         source_columns, source_note = notes)
readr::write_csv(dict, file.path(DIR, "data_dictionary.csv"))
message("-> data_dictionary.csv (", nrow(dict), " canonical variables)")

# ── QC: cross-country distributions after harmonisation ─────────────────────
#
# The detector is the ratio of country SPREADS (IQR), not of country medians.
# A unit or scale error multiplies the whole distribution, so it inflates the
# spread ratio; genuine geography (Gambia is flat, Tanzania is not) moves the
# LEVEL while spreads stay comparable. Using medians alone also misfires badly
# on zero-centred variables -- every AlphaEarth dimension and the PDSI index sit
# near zero, so their median ratio is enormous and meaningless.
#
# Level separation is reported separately and labelled as possibly real, so a
# reader can tell "these are on different scales" from "these are different
# places".
SPREAD_FLAG <- 20
LEVEL_FLAG  <- 10

qc <- lapply(vars, function(v) {
  s <- pooled %>% group_by(country) %>%
    summarise(med = stats::median(.data[[v]], na.rm = TRUE),
              iqr = stats::IQR(.data[[v]], na.rm = TRUE),
              n_na = sum(is.na(.data[[v]])), n = dplyr::n(), .groups = "drop")
  meds <- s$med[is.finite(s$med)]; iqrs <- s$iqr[is.finite(s$iqr) & s$iqr > 0]
  spread_ratio <- if (length(iqrs) >= 2) max(iqrs) / min(iqrs) else NA_real_
  typ_iqr <- if (length(iqrs)) stats::median(iqrs) else NA_real_
  level_sep <- if (length(meds) >= 2 && is.finite(typ_iqr) && typ_iqr > 0)
    (max(meds) - min(meds)) / typ_iqr else NA_real_
  tibble::tibble(canonical = v, n_countries_with_data = sum(is.finite(s$med)),
                 med_min = suppressWarnings(min(meds)), med_max = suppressWarnings(max(meds)),
                 spread_ratio = spread_ratio, level_separation_iqr = level_sep,
                 pct_missing = round(100 * sum(s$n_na) / sum(s$n), 1))
}) %>% dplyr::bind_rows() %>%
  left_join(dict %>% select(canonical, family, canonical_unit, unit_conversion), by = "canonical") %>%
  mutate(has_declared_conversion = !is.na(unit_conversion),
         flag = dplyr::case_when(
           !is.na(spread_ratio) & spread_ratio > SPREAD_FLAG & has_declared_conversion ~
             "UNRESOLVED - conversion declared but spreads still differ by >20x",
           !is.na(spread_ratio) & spread_ratio > SPREAD_FLAG ~
             "SCALE - country spreads differ by >20x; check for a unit/export defect or a genuinely different range",
           pct_missing > 20 ~ "SPARSE - >20% missing",
           !is.na(level_separation_iqr) & level_separation_iqr > LEVEL_FLAG ~
             "LEVEL - countries differ in level by >10 within-country IQRs (may be real geography)",
           TRUE ~ "ok")) %>%
  arrange(desc(spread_ratio))
readr::write_csv(qc, file.path(DIR, "qc_report.csv"))

unmatched <- map %>% filter(action == "unmatched") %>% count(base, sort = TRUE)
dropped   <- map %>% filter(action == "drop") %>% count(reason, sort = TRUE)

lines <- c(
  "# Harmonised admin-2 predictor dataset - QC summary",
  "",
  sprintf("Generated: %s", Sys.Date()),
  sprintf("- Countries: %s", paste(sort(unique(pooled$country)), collapse = ", ")),
  sprintf("- Admin-2 areas: %d", nrow(pooled)),
  sprintf("- Canonical predictors shared by every country: %d", length(vars)),
  "",
  "## Cross-country scale check",
  "",
  paste("Detector: ratio of country IQRs (a unit error multiplies the spread;",
        "geography moves the level). Level separation is reported separately."),
  "",
  paste0("- ", qc %>% filter(grepl("^(UNRESOLVED|SCALE|SPARSE)", flag)) %>%
           mutate(l = sprintf("`%s` (%s): spread ratio %.1fx, medians %.4g to %.4g -- %s",
                              canonical, family, spread_ratio, med_min, med_max, flag)) %>%
           pull(l)),
  "",
  "Level differences (likely real geography, listed for transparency):",
  "",
  paste0("- ", qc %>% filter(grepl("^LEVEL", flag)) %>%
           mutate(l = sprintf("`%s` (%s): medians %.4g to %.4g, %.0f IQRs apart",
                              canonical, family, med_min, med_max,
                              level_separation_iqr)) %>% pull(l)),
  "",
  sprintf("%d of %d variables pass unflagged.", sum(qc$flag == "ok"), nrow(qc)),
  "",
  "## Columns dropped, by reason",
  "",
  paste0("- ", dropped$n, "x ", dropped$reason),
  "",
  "## Raw names matched by no rule",
  "",
  if (nrow(unmatched)) paste0("- ", unmatched$n, "x `", unmatched$base, "`") else
    "- none: every raw column was classified.")
writeLines(unlist(lines), file.path(DIR, "qc_summary.md"))
message("-> qc_report.csv and qc_summary.md")
message(sprintf("   %d/%d variables pass; %d flagged for review",
                sum(qc$flag == "ok"), nrow(qc), sum(qc$flag != "ok")))
message("\nStage 4 complete.")
