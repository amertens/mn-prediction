# =============================================================================
# scripts/run_nested_loco.R
#
# WS1. Run the nested-versus-original LOCO selection comparison
# (R/corrected/p10_nested_loco.R) off the cached {targets} store.
#
# Reads merged_<country> from the store rather than rebuilding the pipeline,
# the same pattern scripts/compare_covariate_hygiene.R uses. The selection
# procedure is glm-scored and cheap; the expensive upstream merges are already
# built.
#
# Profiles (ANALYSIS_PROFILE, also recorded in pipeline_params):
#   smoke  child_iron only, two held-out countries. For development.
#   full   all four shared outcomes, every held-out country.
#
# Run:
#   Rscript scripts/run_nested_loco.R
#   ANALYSIS_PROFILE=smoke Rscript scripts/run_nested_loco.R
#
# Writes results/tables/corrected/loco_nested_selection*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

params  <- get_pipeline_params()
PROFILE <- params$analysis_profile
SEED    <- params$seed

# Which estimator scores the two schemes. "glm" is the workhorse: cheap enough
# for the full grid and identical to the metric the selection search optimises.
# "sl" reproduces the estimator behind the published transport numbers and is
# the one a claim about those numbers must be made on. NESTED_SCORERS=glm,sl
# runs both into one table, distinguished by the `scorer` column.
SCORERS <- trimws(strsplit(Sys.getenv("NESTED_SCORERS", "glm"), ",")[[1]])
sl_learners <- if ("sl" %in% SCORERS) setup_mlr3_learners(params) else NULL

configs   <- get_country_configs()
countries <- names(configs)

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)
all_merged <- setNames(lapply(countries, function(c) read_target(paste0("merged_", tolower(c)))),
                       countries)
missing <- countries[vapply(all_merged, is.null, logical(1))]
if (length(missing))
  stop("merged_* not in the {targets} store for: ", paste(missing, collapse = ", "),
       ". Build those targets first.", call. = FALSE)

if (PROFILE == "smoke") {
  outcomes <- "child_iron"
  holdouts <- countries[seq_len(min(2L, length(countries)))]
} else {
  outcomes <- NESTED_LOCO_OUTCOMES
  holdouts <- countries
}

cat(sprintf("\n[ws1] profile=%s  pipeline_mode=%s  scorers=%s  seed=%d\n",
            PROFILE, params$mode, paste(SCORERS, collapse = "+"), SEED))
cat(sprintf("[ws1] countries: %s\n", paste(countries, collapse = ", ")))
cat(sprintf("[ws1] outcomes:  %s\n", paste(outcomes, collapse = ", ")))
cat(sprintf("[ws1] held out:  %s\n", paste(holdouts, collapse = ", ")))
cat(sprintf("[ws1] merged rows: %s\n\n",
            paste(sprintf("%s=%d", countries, vapply(all_merged, nrow, integer(1))),
                  collapse = ", ")))

t0 <- Sys.time()
res <- run_nested_loco(all_merged, configs,
                       outcomes = outcomes,
                       held_out_countries = holdouts,
                       seed = SEED,
                       country_weighted_also = TRUE,
                       scorers = SCORERS,
                       sl_learners = sl_learners,
                       params = params)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (is.null(res$metrics)) stop("No nested-LOCO metrics produced.", call. = FALSE)

# Recorded per row so a table holding both scorers still says which SL stack
# produced the "sl" rows. The glm scorer does not depend on it.
res$metrics$pipeline_mode <- params$mode
res$metrics$analysis_profile <- PROFILE

# The smoke profile is for development and must not be mistaken for the
# reportable grid, so it writes to its own filenames.
out_dir <- here("results", "tables", "corrected")
paths <- if (PROFILE == "smoke") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  s <- summarize_nested_loco(res$metrics)
  p <- file.path(out_dir, c("loco_nested_selection_SMOKE.csv",
                            "loco_nested_selection_paired_SMOKE.csv",
                            "loco_nested_selection_summary_SMOKE.csv",
                            "loco_nested_selected_vars_SMOKE.csv"))
  readr::write_csv(res$metrics, p[1]); readr::write_csv(s$paired, p[2])
  readr::write_csv(s$summary,  p[3]); readr::write_csv(res$selections, p[4])
  p
} else {
  write_nested_loco_tables(res, out_dir)
}

s <- summarize_nested_loco(res$metrics)
cat("\n================ selected sets ================\n")
print(as.data.frame(res$selections[, c("scheme", "held_out", "n_selected",
                                       "n_candidates", "final_glm_auc")]),
      row.names = FALSE)
cat("\n========= paired deltas (nested - original) =========\n")
print(as.data.frame(s$summary), row.names = FALSE)
cat(sprintf("\n[ws1] elapsed %.1f min\n", elapsed))
cat("[ws1] wrote:\n"); cat(paste0("  ", paths, collapse = "\n"), "\n")
