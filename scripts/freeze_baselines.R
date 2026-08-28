# =============================================================================
# scripts/freeze_baselines.R
#
# WS0. Snapshot the headline result tables to results/tables/frozen_2026-08/
# BEFORE any refit, so every downstream diff has a fixed baseline.
#
# This script COPIES. It never moves, edits or deletes a source file, and it
# refuses to overwrite an existing frozen snapshot unless FREEZE_FORCE=true.
#
# The frozen tree is grouped by source location rather than mirroring the full
# repo-relative path, which would nest results/tables/ inside itself:
#
#   results/tables/*.csv            -> frozen_2026-08/tables/
#   results/tables/corrected/*      -> frozen_2026-08/tables_corrected/
#   results/sensitivity/*.csv       -> frozen_2026-08/sensitivity/
#   sandbox_parsimony/out/*.csv     -> frozen_2026-08/sandbox_parsimony_out/
#
# MANIFEST.csv records the source path for each frozen file, so
# scripts/regression_gate.R (WS9) can map a regenerated table back to its
# baseline without guessing.
#
# Usage:
#   Rscript scripts/freeze_baselines.R
# =============================================================================

suppressWarnings(suppressMessages(library(here)))

FROZEN_DIR <- here::here("results", "tables", "frozen_2026-08")
FORCE      <- tolower(Sys.getenv("FREEZE_FORCE", "false")) %in% c("1", "true", "yes")

# ── The baselines, as identified in the Stage 0 audit ────────────────────────
# `group` selects the subdirectory under FROZEN_DIR; `why` is carried into the
# manifest so the reason a file was frozen survives without a separate doc.
.spec <- function(group, why, paths) {
  data.frame(group = group, why = why, source_path = paths, stringsAsFactors = FALSE)
}

BASELINES <- do.call(rbind, list(

  .spec("tables", "within-country Admin-2 metrics and reliability-ceiling share",
        c("results/tables/area_comparison_all.csv",
          "results/tables/area_comparison_all_PREFULL.csv",
          "results/tables/area_loco_comparison.csv")),

  .spec("tables", "individual-level SuperLearner CV performance",
        c("results/tables/cv_performance_all.csv",
          "results/tables/admin2_error_all.csv",
          "results/tables/all_methods_child_iron.csv")),

  .spec("tables", "national recovery",
        "results/tables/national_estimates_all.csv"),

  .spec("tables", "LOCO benchmark table, includes spatial_plus_soil",
        "results/tables/benchmarks_all.csv"),

  .spec("tables", "parsimonious / transportability model outputs",
        c("results/tables/transportability_loco_results.csv",
          "results/tables/transportability_area_loco_metrics.csv",
          "results/tables/transportability_area_loco_predictions.csv",
          "results/tables/transportability_area_selected_vars.csv",
          "results/tables/transportability_experiments_folds.csv",
          "results/tables/transportability_experiments_summary.csv",
          "results/tables/transportability_recipe_bakeoff.csv",
          "results/tables/sl_prescreened_main.csv",
          "results/tables/single_var_importance.csv")),

  .spec("tables", "within-country reliability ceiling",
        "results/tables/transportability_within_country_ceiling.csv"),

  .spec("tables", "calibration and discrimination diagnostics",
        c("results/tables/calibration_tables.csv",
          "results/tables/diagnostics_binary.csv",
          "results/tables/diagnostics_binary_calibrated.csv",
          "results/tables/diagnostics_continuous.csv",
          "results/tables/domain_ablation_all.csv")),

  .spec("sensitivity", "covariate-hygiene v1 to v2 comparison",
        "results/sensitivity/covariate_hygiene_comparison.csv"),

  .spec("sensitivity", "distributional-estimator prototype (WS4 precedent)",
        c("results/sensitivity/distributional_prototype_metrics.csv",
          "results/sensitivity/distributional_transport_childiron.csv",
          "results/sensitivity/distributional_within_heterosked.csv")),

  .spec("sensitivity", "anchored calibration (WS7 precedent)",
        "results/sensitivity/gradient_anchored_calibration.csv"),

  .spec("sensitivity", "GEE legacy-parity validation reference procedure",
        "results/sensitivity/gee_legacy_parity_validation_sierraleone.csv"),

  .spec("sandbox_parsimony_out", "national recovery, sandbox critical review",
        "sandbox_parsimony/out/national_recovery.csv"),

  .spec("sandbox_parsimony_out",
        "spatial_plus_soil under original vs leakage-corrected selection",
        "sandbox_parsimony/out/spatial_plus_soil_rescored.csv"),

  .spec("sandbox_parsimony_out", "parsimonious LOCO headline",
        c("sandbox_parsimony/out/loco_headline.csv",
          "sandbox_parsimony/out/loco_bakeoff.csv",
          "sandbox_parsimony/out/loco_spatial.csv")),

  .spec("sandbox_parsimony_out", "covariate-hygiene effect, sandbox",
        "sandbox_parsimony/out/hygiene_effect.csv"),

  .spec("sandbox_parsimony_out", "decision accuracy and level reliability",
        c("sandbox_parsimony/out/decision_accuracy.csv",
          "sandbox_parsimony/out/level_reliability.csv",
          "sandbox_parsimony/out/level_reliability_ext.csv"))
))

# Everything the corrected-methods layer emits, taken as a directory so a new
# corrected table cannot silently escape the freeze.
corrected_files <- list.files(here::here("results", "tables", "corrected"),
                              full.names = FALSE)
if (length(corrected_files)) {
  BASELINES <- rbind(BASELINES, .spec(
    "tables_corrected", "corrected-methods (P1-P7, P9) bundle",
    file.path("results", "tables", "corrected", corrected_files)
  ))
}

BASELINES$source_path <- gsub("\\\\", "/", BASELINES$source_path)

# ── Guard against clobbering an existing freeze ──────────────────────────────
if (dir.exists(FROZEN_DIR) && !FORCE) {
  stop("Frozen baseline already exists at results/tables/frozen_2026-08/. ",
       "A freeze is meant to happen once, before any refit. ",
       "Set FREEZE_FORCE=true only if you are deliberately re-freezing.",
       call. = FALSE)
}

# ── Copy ─────────────────────────────────────────────────────────────────────
rows    <- list()
missing <- character(0)

for (i in seq_len(nrow(BASELINES))) {
  src_rel <- BASELINES$source_path[i]
  src_abs <- here::here(src_rel)

  if (!file.exists(src_abs)) {
    missing <- c(missing, src_rel)
    next
  }

  dest_dir <- file.path(FROZEN_DIR, BASELINES$group[i])
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_abs <- file.path(dest_dir, basename(src_rel))

  ok <- file.copy(src_abs, dest_abs, overwrite = TRUE, copy.date = TRUE)
  if (!ok) stop("Failed to copy ", src_rel, call. = FALSE)

  # md5 of the COPY, so the manifest verifies what is actually stored.
  rows[[length(rows) + 1L]] <- data.frame(
    frozen_path = file.path("results", "tables", "frozen_2026-08",
                            BASELINES$group[i], basename(src_rel)),
    source_path = src_rel,
    group       = BASELINES$group[i],
    why         = BASELINES$why[i],
    md5         = unname(tools::md5sum(dest_abs)),
    bytes       = as.integer(file.info(dest_abs)$size),
    source_mtime = format(file.info(src_abs)$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, rows)
manifest$frozen_date <- format(Sys.Date(), "%Y-%m-%d")
manifest <- manifest[order(manifest$group, manifest$source_path), ]

write.csv(manifest, file.path(FROZEN_DIR, "MANIFEST.csv"), row.names = FALSE)

# ── Report ───────────────────────────────────────────────────────────────────
cat(sprintf("Frozen %d files to results/tables/frozen_2026-08/\n", nrow(manifest)))
for (g in sort(unique(manifest$group))) {
  cat(sprintf("  %-22s %d\n", g, sum(manifest$group == g)))
}
if (length(missing)) {
  cat(sprintf("\nNOT FOUND (%d) -- recorded as absent, not fabricated:\n", length(missing)))
  cat(paste0("  ", missing, collapse = "\n"), "\n")
}
cat("\nManifest: results/tables/frozen_2026-08/MANIFEST.csv\n")
