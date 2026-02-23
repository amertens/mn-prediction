# =============================================================================
# scripts/run_gambia_pipeline.R
#
# End-to-end runner for the Gambia micronutrient deficiency prediction pipeline.
#
# Steps executed:
#   00  Library loading + parallel plan
#   01  Load data and construct per-outcome datasets
#   02  Fit SuperLearner models (continuous + binary)
#   03  Predict and aggregate to Admin1 prevalence + CV performance
#   04  Bootstrap uncertainty intervals (Admin1 + national)
#   05  Domain ablation study
#
# Outputs (see src/analysis/*.R for full details):
#   results/models/res_GW_Gambia_SL_child_vitA.rds
#   results/models/res_GW_Gambia_SL_women_vitA.rds
#   results/models/res_GW_Gambia_SL_child_iron.rds
#   results/models/res_GW_Gambia_SL_mom_iron.rds
#   results/models/res_bin_GW_Gambia_SL_*.rds
#   results/tables/gambia_admin1_prevalence_*.csv
#   results/tables/gambia_admin1_ci_*.csv
#   results/tables/gambia_national_ci_*.csv
#   results/tables/gambia_cv_performance.csv
#   results/tables/gambia_domain_ablation_summary.csv
#   results/figures/admin1_scatter_*.png
#
# Usage:
#   Rscript scripts/run_gambia_pipeline.R
#   # or interactively:
#   source(here::here("scripts/run_gambia_pipeline.R"))
# =============================================================================

set.seed(12345)

# ============================================================================
# 00. Libraries and parallel plan
# ============================================================================
cat("=============================================================\n")
cat("  Gambia Micronutrient Prediction Pipeline\n")
cat("=============================================================\n\n")

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyverse)
  library(haven)
  library(purrr)
  library(labelled)
  library(sl3)
  library(origami)
  library(tlverse)
  library(caret)
  library(data.table)
  library(ck37r)
  library(SuperLearner)
  library(future)
  library(future.apply)
  library(washb)
  library(recipes)
  library(readr)
  library(ggplot2)
  library(pROC)
})

# Optional packages (used for scatter plots; graceful degradation if absent)
has_ggrepel <- requireNamespace("ggrepel",  quietly = TRUE)
has_scales  <- requireNamespace("scales",   quietly = TRUE)
if (!has_ggrepel) message("[info] ggrepel not installed; labels will overlap in scatter plots")
if (!has_scales)  message("[info] scales not installed; axis labels will be unformatted")

# Parallel processing
plan(multicore, workers = max(1L, floor(availableCores() / 2)))
options(future.globals.maxSize = 5 * 1024^3)   # 5 GB

cat(sprintf("[00] Workers: %d\n\n", nbrOfWorkers()))

# ============================================================================
# Source helper functions and SL setup
# ============================================================================
source(here::here("src/0-functions.R"))
source(here::here("src/0-SL-setup.R"))          # defines DHS_SL, slmod, slmod2_bin, …
source(here::here("src/DHS/DHS_functions.R"))    # stub; DHS_SL is in 0-SL-setup.R
source(here::here("src/DHS/DHS_variable_recode.R"))

# ============================================================================
# Load pipeline configuration
# ============================================================================
source(here::here("src/analysis/config.R"))      # defines `cfg`

# ============================================================================
# Step 01 – Load and construct datasets
# ============================================================================
source(here::here("src/analysis/01_load_and_construct.R"))
# Produces: gw_data_list, Xvars_full, Xvars, domain_vars

# ============================================================================
# Step 02 – Fit SuperLearner models
# ============================================================================
source(here::here("src/analysis/02_fit_sl_models.R"))
# Produces: sl_results_cont, sl_results_bin

# ============================================================================
# Step 03 – Predict, aggregate to Admin1, compute CV performance
# ============================================================================
source(here::here("src/analysis/03_predict_and_aggregate_admin1.R"))
# Produces: admin1_prev_list, cv_perf_table

# ============================================================================
# Step 04 – Bootstrap uncertainty intervals
# ============================================================================
source(here::here("src/analysis/04_bootstrap_uncertainty.R"))

# ============================================================================
# Step 05 – Domain ablation
# ============================================================================
source(here::here("src/analysis/05_domain_ablation.R"))

# ============================================================================
# Summary: list output files
# ============================================================================
cat("\n=============================================================\n")
cat("  Pipeline complete.  Output files:\n")
cat("=============================================================\n\n")

out_files <- c(
  list.files(cfg$out_models,  full.names = TRUE, pattern = "\\.rds$"),
  list.files(cfg$out_tables,  full.names = TRUE, pattern = "\\.csv$"),
  list.files(cfg$out_figures, full.names = TRUE, pattern = "\\.png$")
)

for (f in sort(out_files)) {
  sz <- file.size(f)
  cat(sprintf("  %-60s  %s\n", basename(f),
              if (sz < 1e6) sprintf("%.1f KB", sz / 1e3)
              else          sprintf("%.1f MB", sz / 1e6)))
}

cat("\nDone.\n")
