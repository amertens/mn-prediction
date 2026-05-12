# =============================================================================
# scripts/run_full_pipeline.R
#
# Master runner for the Gambia micronutrient deficiency prediction pipeline.
# Executes all steps in dependency order, with timing and skip-if-done logic.
#
# Pipeline steps:
#   Step 1:  Individual-level SL pipeline (01-05)
#            - Load data, fit SL, predict Admin1, bootstrap, domain ablation
#   Step 2:  Admin2 predictions & area-level model (06)
#            - Individual SL -> Admin2 prevalence, survey-weighted estimates,
#              area-level GEE model for unsurveyed areas, error analysis
#   Step 3:  DHS indicator aggregation at Admin2 (surveyPrev)
#            - Direct, Fay-Herriot, and cluster-level model estimates
#   Step 4:  Prepare stakeholder walkthrough data
#            - Bridge from pipeline outputs to QMD-ready .rds files
#   Step 5:  Render stakeholder walkthrough QMD (optional)
#
# Usage:
#   # RECOMMENDED: use targets (handles caching automatically):
#   targets::tar_make()                    # fast mode (default)
#   Sys.setenv(PIPELINE_MODE = "full")     # switch to publication mode
#   targets::tar_make()
#
#   # ALTERNATIVE: this legacy sequential runner:
#   source("scripts/run_full_pipeline.R")          # fast mode (default)
#   PIPELINE_MODE <- "full"
#   source("scripts/run_full_pipeline.R")           # full mode
#
#   # Or set flags to skip steps you don't need:
#   SKIP_SL_PIPELINE    <- TRUE   # skip if models already saved
#   SKIP_ADMIN2         <- FALSE
#   SKIP_DHS_SURVEYPREV <- TRUE   # skip if DHS admin2 tables already saved
#   SKIP_WALKTHROUGH    <- FALSE
#   RENDER_QMD          <- TRUE   # set TRUE to also render the QMD
#   source("scripts/run_full_pipeline.R")
# =============================================================================

# ── Default flags (set these BEFORE sourcing to override) ────────────────────
if (!exists("SKIP_SL_PIPELINE"))    SKIP_SL_PIPELINE    <- FALSE
if (!exists("SKIP_ADMIN2"))         SKIP_ADMIN2         <- FALSE
if (!exists("SKIP_DHS_SURVEYPREV")) SKIP_DHS_SURVEYPREV <- FALSE
if (!exists("SKIP_WALKTHROUGH"))    SKIP_WALKTHROUGH    <- FALSE
if (!exists("RENDER_QMD"))          RENDER_QMD          <- FALSE
if (!exists("PIPELINE_MODE"))       PIPELINE_MODE       <- "fast"

library(here)

cat("\n")
cat("=============================================================\n")
cat("  Gambia Micronutrient Prediction — Full Pipeline\n")
cat("=============================================================\n")
cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  Working directory: %s\n\n", here::here()))

# Helper: run a step with timing
run_step <- function(step_num, step_name, expr, skip = FALSE) {
  cat(sprintf("\n--- Step %s: %s ", step_num, step_name))
  if (skip) {
    cat("[SKIPPED] ---\n")
    return(invisible(NULL))
  }
  cat("---\n")
  t0 <- proc.time()
  tryCatch(
    expr,
    error = function(e) {
      cat(sprintf("\n  *** FAILED: %s ***\n", conditionMessage(e)))
      cat(sprintf("  You can skip this step next time with SKIP_... <- TRUE\n"))
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]
  mins <- floor(elapsed / 60)
  secs <- round(elapsed %% 60)
  cat(sprintf("  [Step %s completed in %dm %ds]\n", step_num, mins, secs))
}


# =============================================================================
# STEP 1: Individual-level SuperLearner pipeline (scripts 01-05)
# =============================================================================
# Fits continuous + binary SL models for all outcomes, predicts to Admin1,
# runs cluster bootstrap uncertainty, and domain ablation.
#
# Inputs:  data/IPD/Gambia/Gambia_merged_dataset.rds
# Outputs: results/models/*.rds, results/tables/*.csv, results/figures/*.png
# =============================================================================
run_step("1", "Individual-level SL pipeline (01-05)", {
  source(here::here("scripts", "run_gambia_pipeline.R"))
}, skip = SKIP_SL_PIPELINE)


# =============================================================================
# STEP 2: Admin2 predictions, mapping, & area-level model (script 06)
# =============================================================================
# Uses the SL model from Step 1 (or fits a minimal fallback) to produce
# Admin2-level predictions.  Then fits an area-level elastic net model
# using GEE raster covariates to predict unsurveyed Admin2 areas.
#
# Inputs:  data/IPD/Gambia/Gambia_merged_dataset.rds
#          results/models/res_bin_GW_Gambia_SL_child_vitA.rds (optional)
#          data/Gambia_GEE_rasters/*.tif
# Outputs: results/admin2/tables/*.csv, results/admin2/figures/*.png
# =============================================================================
run_step("2", "Admin2 predictions & area-level model (06)", {
  source(here::here("scripts", "06_admin2_predictions_map.R"))
}, skip = SKIP_ADMIN2)


# =============================================================================
# STEP 3: DHS indicator aggregation at Admin2 via surveyPrev
# =============================================================================
# Downloads DHS microdata via rdhs, computes Admin2-level prevalence for
# 22 standard DHS indicators using direct, Fay-Herriot, and cluster-level
# model estimators from the surveyPrev package.
#
# Requirements: rdhs credentials configured, INLA installed, internet access
# Inputs:  DHS API (downloads microdata)
# Outputs: data/DHS/clean/gambia_admin2_*.csv
#
# NOTE: This step requires DHS API credentials and INLA.  Skip if those
#       are not set up, or if the output files already exist.
# =============================================================================
run_step("3", "DHS Admin2 aggregation (surveyPrev)", {
  source(here::here("src", "Gambia", "gambia_DHS_aggregation_admin2.R"))
}, skip = SKIP_DHS_SURVEYPREV)


# =============================================================================
# STEP 4: Prepare stakeholder walkthrough data
# =============================================================================
# Bridges pipeline outputs (Steps 1-2) into QMD-ready .rds files.
# If full pipeline outputs are not yet available, creates synthetic
# predictions so the walkthrough still renders with realistic visuals.
#
# Inputs:  data/IPD/Gambia/Gambia_merged_dataset.rds
#          results/admin2/tables/* (from Step 2, if available)
# Outputs: docs/helpers/walkthrough_data/*.rds
# =============================================================================
run_step("4", "Prepare walkthrough data", {
  source(here::here("docs", "helpers", "prepare_gambia_admin2_walkthrough.R"))
}, skip = SKIP_WALKTHROUGH)


# =============================================================================
# STEP 5: Render stakeholder walkthrough QMD (optional)
# =============================================================================
# Renders the Quarto document to HTML.  Only runs if RENDER_QMD = TRUE.
#
# Inputs:  docs/stakeholder_walkthrough_proxy_vad.qmd
#          docs/helpers/walkthrough_data/*.rds (from Step 4)
# Outputs: docs/stakeholder_walkthrough_proxy_vad.html
# =============================================================================
run_step("5", "Render stakeholder walkthrough QMD", {
  qmd_path <- here::here("docs", "stakeholder_walkthrough_proxy_vad.qmd")
  if (!file.exists(qmd_path)) {
    cat("  QMD file not found — skipping render.\n")
  } else {
    quarto::quarto_render(qmd_path)
    cat(sprintf("  Rendered: %s\n",
                sub("\\.qmd$", ".html", basename(qmd_path))))
  }
}, skip = !RENDER_QMD)


# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("=============================================================\n")
cat("  Pipeline complete\n")
cat(sprintf("  Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n\n")

# List all output files
output_dirs <- c(
  here::here("results", "models"),
  here::here("results", "tables"),
  here::here("results", "figures"),
  here::here("results", "admin2", "tables"),
  here::here("results", "admin2", "figures"),
  here::here("docs", "helpers", "walkthrough_data")
)

cat("  Output files:\n")
for (d in output_dirs) {
  if (!dir.exists(d)) next
  files <- list.files(d, full.names = TRUE)
  if (length(files) == 0) next
  cat(sprintf("\n  %s/\n", basename(d)))
  for (f in sort(files)) {
    sz <- file.size(f)
    cat(sprintf("    %-55s %s\n", basename(f),
                if (sz < 1e6) sprintf("%.1f KB", sz / 1e3)
                else          sprintf("%.1f MB", sz / 1e6)))
  }
}

cat("\n  Quick-start for next time (skip completed steps):\n")
cat("    SKIP_SL_PIPELINE <- TRUE\n")
cat("    SKIP_DHS_SURVEYPREV <- TRUE\n")
cat("    source(\"scripts/run_full_pipeline.R\")\n\n")
