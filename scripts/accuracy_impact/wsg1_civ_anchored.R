# =============================================================================
# scripts/accuracy_impact/wsg1_civ_anchored.R
#
# WS-G. Cote d'Ivoire: a region-level anchored export, with NO covariate tilt.
#
# WHY NO COVARIATE TILT
# ---------------------
# Two independent reasons, and either alone would be sufficient.
#
# 1. IT IS NOT FEASIBLE. The pooled model matches covariates by EXACT NAME, and
#    a name that is absent is silently zero rather than an error. The harmonized
#    predictor set is 294 columns dominated by blocks Cote d'Ivoire has none of:
#    dhs 97, aef 64, soil 38, tclim 13, spam 12, precip 10. Its 83 cached GEE
#    rasters would populate roughly 20 to 30 names and zero the rest, which
#    produces a confident number with no content behind it.
#
# 2. IT WOULD NOT HELP. Under leave-one-country-out, the flat anchor beats every
#    covariate arm on error: 9.85 pp against 13.50 for the un-anchored
#    transported model and 10.81 for the best anchored covariate arm. The
#    covariates buy about 0.20 of correlation and cost accuracy.
#
# WHY NOT ANCHOR ON GHANA
# -----------------------
# It was considered and it is wrong. VMNIS holds both countries:
#
#              non-pregnant women     Cote d'Ivoire      Ghana      gap
#   folate deficiency                        86.1%       53.8%    32 pp
#   vitamin B12 deficiency                   18.1%        6.9%    11 pp
#
# Borrowing Ghana's level would inject a 32-point error, which is the exact
# failure mode that makes cross-country transport fail. Cote d'Ivoire's own
# VMNIS entry is used instead, and it is available at NINE SUB-NATIONAL UNITS,
# so the anchor sits at region rather than national resolution. The anchoring
# budget curve shows Admin-1 error running well below Admin-2 at every survey
# fraction, so a region-level anchor is the better product regardless.
#
# WHAT THIS IS NOT
# ----------------
# It is not a validated prediction. There are no Cote d'Ivoire biomarker labels
# in this project -- the rdhs cache holds two GPS files and nothing else -- so
# nothing here is scored. Every row carries `scored = FALSE`.
#
#   Rscript scripts/accuracy_impact/wsg1_civ_anchored.R
# -> results/deliverables/civ_anchored_regions.csv
# -> docs/requests/civ_labels_request.md
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})

OUT <- here("results", "deliverables", "civ_anchored_regions.csv")
REQ <- here("docs", "requests", "civ_labels_request.md")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(REQ), showWarnings = FALSE, recursive = TRUE)

v <- as.data.frame(readRDS(here("data", "national", "vmnis_national_slim.rds")))
civ <- v[v$CountryCode == "CIV", , drop = FALSE]
if (!nrow(civ)) stop("[wsg1] no Cote d'Ivoire rows in VMNIS")

# The national rows have a blank Representativenessname; the sub-national rows
# name one of the nine eco-regions. Both are kept and labelled, because a
# consumer needs to know which resolution a number came from.
civ$unit <- trimws(as.character(civ$Representativenessname))
civ$level <- ifelse(nchar(civ$unit) == 0, "national", "eco_region")

keep <- c("mn_group", "Population", "year", "unit", "level", "Samplesize",
          "prev", "Indicator", "Indicatorunit", "Deficiencycutoff",
          "Dataadjustedfor", "Surveymethodology")
out <- civ[, intersect(keep, names(civ)), drop = FALSE]
out <- out[order(out$mn_group, out$level, out$unit), ]

# ---------------------------------------------------------------------------
# National reference for each nutrient, and each region's deviation from it.
# The DEVIATION is the targeting-relevant quantity: it says which eco-regions
# are worse than the country average, which is what a programme allocates on.
# ---------------------------------------------------------------------------
nat <- out |> filter(level == "national") |> group_by(mn_group) |>
  summarise(national_prev = stats::weighted.mean(prev, pmax(Samplesize, 1), na.rm = TRUE),
            national_n = sum(Samplesize, na.rm = TRUE), .groups = "drop")
out <- out |> left_join(nat, by = "mn_group") |>
  mutate(deviation_pp = round(100 * (prev - national_prev), 2),
         prev_pp = round(100 * prev, 2),
         national_prev_pp = round(100 * national_prev, 2))

# Rank the eco-regions within each nutrient, worst first.
# Ranked among the ECO-REGIONS ONLY. Ranking over the whole nutrient group
# includes the national rows in the denominator and leaves gaps in the
# sequence (1,2,3,5,6,8,...), which reads as though regions were dropped.
out$rank_within_nutrient <- NA_integer_
for (g in unique(out$mn_group)) {
  i <- which(out$mn_group == g & out$level == "eco_region")
  if (!length(i)) next
  out$rank_within_nutrient[i] <- rank(-out$prev[i], ties.method = "min")
}
out$n_eco_regions <- ave(as.integer(out$level == "eco_region"), out$mn_group, FUN = sum)

# ---------------------------------------------------------------------------
# The provenance columns. Every one of these is a reason not to over-read the
# numbers, and they travel with the rows rather than living in a README.
# ---------------------------------------------------------------------------
out$scored <- FALSE
out$scoring_blocked_reason <- "no Cote d'Ivoire biomarker labels in this project"
out$covariate_tilt_applied <- FALSE
out$covariate_tilt_reason <- paste(
  "harmonized predictors unavailable for CIV (dhs/aef/soil/tclim/spam/precip",
  "blocks absent); and under LOCO the flat anchor outperforms every covariate",
  "arm on error, 9.85 pp against 13.50")
out$anchor_source <- "VMNIS, Cote d'Ivoire's own survey"
out$anchor_not_ghana_reason <- "Ghana folate 53.8% against CIV 86.1%, a 32 pp level error"
out$cutoff_verified <- is.finite(suppressWarnings(as.numeric(out$Deficiencycutoff)))
out$survey_year <- out$year
out$years_stale <- 2026L - as.integer(out$year)
out$unit_is_administrative <- FALSE
out$unit_note <- "survey eco-regions, not GADM admin units; crosswalk outstanding"

readr::write_csv(out, OUT)

cat(sprintf("[wsg1] %d rows: %d national, %d eco-region\n", nrow(out),
            sum(out$level == "national"), sum(out$level == "eco_region")))
cat(sprintf("[wsg1] nutrients: %s\n", paste(sort(unique(out$mn_group)), collapse = ", ")))
cat(sprintf("[wsg1] survey year %s, %d years stale\n",
            paste(unique(out$year), collapse = "/"), max(out$years_stale, na.rm = TRUE)))
cat(sprintf("[wsg1] cutoff verified on %d of %d rows\n",
            sum(out$cutoff_verified), nrow(out)))

cat("\n=== eco-regions, worst first ===\n")
show <- out |> filter(level == "eco_region") |>
  select(mn_group, unit, prev_pp, national_prev_pp, deviation_pp, rank_within_nutrient) |>
  arrange(mn_group, rank_within_nutrient)
print(as.data.frame(show), row.names = FALSE)

writeLines(c(
  "# Request: Cote d'Ivoire micronutrient biomarker labels",
  "",
  "## What is already in hand",
  "",
  "- VMNIS national and nine eco-region prevalences for **folate** and",
  "  **vitamin B12** in non-pregnant women, 2007.",
  "- 83 cached GEE rasters covering NDVI, EVI, LST, population, land cover,",
  "  elevation and accessibility.",
  "- Two DHS GPS files (CIGE61FL, CIGE81FL) in the rdhs cache.",
  "",
  "## What is missing, and what each would unlock",
  "",
  "| Needed | Unlocks |",
  "|---|---|",
  "| Individual biomarker records with cluster IDs | Any *scored* estimate. Nothing in this export is validated. |",
  "| Iron and vitamin A national or sub-national estimates | Those outcomes are absent from VMNIS for CIV; only folate and B12 are covered. |",
  "| The deficiency cutoffs used in the 2007 survey | `Deficiencycutoff` is NA on every row, so the 86.1% folate figure cannot be checked against Ghana's 53.8% for definitional comparability. |",
  "| A crosswalk from the nine eco-regions to GADM Admin-1/Admin-2 | Joining this anchor to any spatial layer. The eco-regions are a survey stratification. |",
  "| A survey more recent than 2007 | The current anchor is 19 years old. |",
  "",
  "## What is deliberately not attempted",
  "",
  "A covariate rank-tilt. The harmonized predictor set is 294 columns and is",
  "dominated by blocks CIV has none of (dhs 97, aef 64, soil 38, tclim 13,",
  "spam 12, precip 10). Because the pooled model matches by exact column name",
  "and treats an absent name as zero rather than as an error, a tilt built on",
  "the ~20-30 available names would be silently mostly-zeros.",
  "",
  "Independently, the leave-one-country-out evidence says the tilt would not",
  "help: a flat anchor scores 9.85 pp mean absolute error against 13.50 for the",
  "un-anchored transported covariate model.",
  "",
  "Generated by `scripts/accuracy_impact/wsg1_civ_anchored.R`."
), REQ)
cat(sprintf("\n-> %s\n-> %s\n", file.path("results","deliverables", basename(OUT)),
            file.path("docs","requests", basename(REQ))))
