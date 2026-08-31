# =============================================================================
# scripts/accuracy_impact/ws0_freeze.R
#
# WS0. Freeze the headline tables this work will revise, before any refit.
#
# WHY THIS EXISTS
# ---------------
# Every workstream that follows either adds rows to an existing table or writes
# a new one. Guardrail 4 forbids overwriting a result under results/tables/, and
# WS9 has to classify every regenerated number as unchanged, a new scheme row,
# or a changed baseline. Both require a byte-exact record of the starting state.
#
# The manifest schema is the one results/tables/frozen_2026-08/MANIFEST.csv
# already uses, so the two freezes can be read by the same code.
#
#   Rscript scripts/accuracy_impact/ws0_freeze.R
# -> results/tables/frozen_2026-09/            (copies)
# -> results/tables/frozen_2026-09/MANIFEST.csv
# =============================================================================
suppressPackageStartupMessages({library(here); library(tools)})

DEST <- here("results", "tables", "frozen_2026-09")
STAMP <- as.character(Sys.Date())

# The Stage 0 table set: everything behind the claims in Sections 2.3, 3, 4, 5,
# 6 and 9 of docs/SESSION_FINDINGS_FOR_REVIEW.md, plus the corrected-layer
# tables that WS3, WS5 and WS7a build on. `why` records which claim it backs, so
# the manifest doubles as the claim-to-file map the register needs.
SPEC <- list(
  list("results/tables/area_comparison_all.csv", "section3",
       "individual-level aggregate, area-level arms, spatial smoother, r_share"),
  list("results/tables/resolution_comparison.csv", "section3",
       "admin-1 vs admin-2 comparison and the WHO-category accuracy claim"),
  list("results/tables/national_estimates_all.csv", "section3",
       "national recovery: 24/24 inside the survey 95% CI, mean abs error 0.96 pp"),
  list("results/tables/bivariate_fdr.csv", "section3",
       "predictors surviving FDR control: 0 of 294 in all 24 cells"),
  list("results/tables/penalized_retained.csv", "section3",
       "penalised regression retained predictors, median 0"),
  list("results/tables/benchmarks_all.csv", "section3",
       "SuperLearner against the best of 21 comparators, 0 of 16 LOCO holdouts"),
  list("results/tables/covariate_validation.csv", "section3",
       "covariate bounds and known-relationship checks"),
  list("results/tables/admin1_arms.csv", "section4",
       "the five anchoring arms; WS2 rescores these cells"),
  list("results/tables/individual_anchor.csv", "section5",
       "proxy against questionnaire, district and cluster; WS3 rescores these"),
  list("results/tables/national_composition.csv", "section6",
       "composition arms including the oracle national level"),
  list("results/tables/national_composition_levels.csv", "section6",
       "predicted against true national levels"),
  list("results/tables/national_vmnis_loco.csv", "section6",
       "the VMNIS LOCO model panels"),
  list("results/tables/national_vmnis_ceiling.csv", "section6",
       "VMNIS variance components and r_max; WS6a recomputes the sampling term"),
  list("results/tables/corrected/protocol_reconciliation.csv", "corrected",
       "per-cell r under each fold protocol; the WS3a starting point"),
  list("results/tables/corrected/protocol_reconciliation_medians.csv", "corrected",
       "median r by protocol, including indiv_region_wt 0.031"),
  list("results/tables/corrected/cv_honesty_compare.csv", "corrected",
       "section 2.3 selection optimism, honest against optimistic schemes"),
  list("results/tables/corrected/admin2_unit_counts.csv", "corrected",
       "per-country polygon, unique-name and analytic unit counts for WS7a"),
  list("results/tables/corrected/subsample_cost_of_accuracy.csv", "corrected",
       "subsample replicate structure WS5 reuses"),
  list("results/tables/corrected/subsample_summary.csv", "corrected",
       "subsample summary by coverage stratum")
)

rows <- list()
missing <- character(0)
for (s in SPEC) {
  src <- s[[1]]; grp <- s[[2]]; why <- s[[3]]
  abs_src <- here(src)
  if (!file.exists(abs_src)) { missing <- c(missing, src); next }
  # Mirror the source path under DEST so frozen_path is unambiguous when two
  # tables share a basename across results/tables and results/tables/corrected.
  rel <- sub("^results/tables/", "", src)
  dst <- file.path(DEST, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(abs_src, dst, overwrite = TRUE)
  fi <- file.info(abs_src)
  rows[[length(rows) + 1L]] <- data.frame(
    frozen_path  = file.path("results/tables/frozen_2026-09", rel),
    source_path  = src,
    group        = grp,
    why          = why,
    md5          = unname(tools::md5sum(abs_src)),
    bytes        = as.integer(fi$size),
    source_mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S"),
    frozen_date  = STAMP,
    stringsAsFactors = FALSE)
}

man <- do.call(rbind, rows)
write.csv(man, file.path(DEST, "MANIFEST.csv"), row.names = FALSE)

cat(sprintf("frozen %d tables -> %s\n", nrow(man), "results/tables/frozen_2026-09"))
if (length(missing)) {
  cat("MISSING (not frozen, recorded here and in the WS0 findings):\n")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
}
print(man[, c("source_path", "group", "bytes")], row.names = FALSE)
