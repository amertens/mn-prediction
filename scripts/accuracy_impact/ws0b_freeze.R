# =============================================================================
# scripts/accuracy_impact/ws0b_freeze.R
#
# WS0 of the signal-and-shipping session. Snapshot the tables this session will
# extend or supersede, before any refit, following the frozen_2026-09 schema.
#
#   Rscript scripts/accuracy_impact/ws0b_freeze.R
# -> results/tables/frozen_2026-09b/ + MANIFEST.csv
# =============================================================================
suppressPackageStartupMessages({library(here); library(tools)})
DEST <- here("results", "tables", "frozen_2026-09b"); STAMP <- as.character(Sys.Date())
SPEC <- list(
  c("results/tables/anchor_controls.csv", "ws2", "anchoring arms; WS-B1 adds a jackknifed flat-region arm"),
  c("results/tables/reliability_empirical.csv", "ws1", "empirical ceilings; WS-A disattenuates with these"),
  c("results/tables/reliability_simulation.csv", "ws1", "estimator bias; WS-B2 extends the simulator"),
  c("results/tables/reliability_analytic_vs_empirical.csv", "ws1", "analytic against empirical ceiling"),
  c("results/tables/reliability_skill_curve.csv", "ws4b", "skill curve; WS-E2 adds decomposition columns"),
  c("results/tables/positive_control_targets.csv", "ws4b", "education and the other pseudo-targets"),
  c("results/tables/individual_arms_2026-09_PARITY.csv", "ws3", "four-cell arms; WS-C1 extends to sixteen"),
  c("results/tables/anchor_implied_shifts.csv", "ws2", "implied shifts; WS-D reuses for the calibration gate"),
  c("results/tables/anchoring_design_curve.csv", "ws5", "anchoring budget replicates"),
  c("results/tables/resolution_sweep.csv", "ws4a", "resolution sweep"),
  c("results/tables/r_share_revised.csv", "ws1", "r_share under both ceilings"),
  c("results/tables/leakage_report_summary.csv", "ws7a", "leakage report; regenerated when predictors change"),
  c("results/tables/assay_guard_label_comparison.csv", "ws7b", "label-derived guard comparison"))
rows <- list(); missing <- character(0)
for (s in SPEC) {
  abs_src <- here(s[1])
  if (!file.exists(abs_src)) { missing <- c(missing, s[1]); next }
  rel <- sub("^results/tables/", "", s[1]); dst <- file.path(DEST, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(abs_src, dst, overwrite = TRUE)
  fi <- file.info(abs_src)
  rows[[length(rows) + 1L]] <- data.frame(
    frozen_path = file.path("results/tables/frozen_2026-09b", rel), source_path = s[1],
    group = s[2], why = s[3], md5 = unname(tools::md5sum(abs_src)),
    bytes = as.integer(fi$size), source_mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S"),
    frozen_date = STAMP, stringsAsFactors = FALSE)
}
man <- do.call(rbind, rows)
utils::write.csv(man, file.path(DEST, "MANIFEST.csv"), row.names = FALSE)
cat(sprintf("frozen %d tables -> results/tables/frozen_2026-09b\n", nrow(man)))
if (length(missing)) cat("MISSING:", paste(missing, collapse = ", "), "\n")
