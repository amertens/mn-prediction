# =============================================================================
# sensitivity/12_aggregation_level_national.R
#
# Ad-hoc runner for the aggregation-level national-prevalence sweep. The core
# logic lives in R/aggregation_level_sweep.R (shared with the
# `aggregation_level_national` _targets target); this script just reads cached
# outcome_data_* targets for a chosen panel and prints/writes the comparison.
#
# QUESTION: does the level at which we fit (cluster / admin-1 / admin-2) change
# NATIONAL-PREVALENCE accuracy? For each country x outcome and each level we
# aggregate the individuals, fit a spatial-block OOF ridge on level-aggregated
# proxies, roll OOF predictions up to a national estimate, and compare to the
# observed national prevalence (|error| in pp).
#
# Reads cached `outcome_data_<lower(key)>_<outcome>` targets (no cache mutation).
# Writes results/sensitivity/aggregation_level_national_{<slug>_<outcome>,summary}.csv
#
# Usage (country arg is the CONFIG KEY, e.g. "SierraLeone" not "Sierra Leone"):
#   Rscript sensitivity/12_aggregation_level_national.R                 # default panel
#   Rscript sensitivity/12_aggregation_level_national.R Ghana women_iron
#
# CAVEATS: cluster rows share ~admin-2-resolution predictors (pseudo-replication
# — this tests the roll-up, not within-district signal; swap in
# R/cluster_aggregation.R::build_cluster_dataset() for that). Ridge is used, not
# the full SuperLearner — consistent across levels, which is what makes the
# comparison fair; absolute numbers won't match national_estimates_all.csv.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(targets); library(dplyr); library(tidyr); library(readr)
})
source(here::here("R", "config.R"))
source(here::here("R", "aggregation_level_sweep.R"))

cfgs <- get_country_configs()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  panel <- list(list(key = args[[1]], outcome = args[[2]]))
} else {
  panel <- list(
    list(key = "Gambia",      outcome = "child_vitA"),
    list(key = "Ghana",       outcome = "women_iron"),
    list(key = "Ghana",       outcome = "child_iron"),
    list(key = "SierraLeone", outcome = "child_iron"),
    list(key = "Malawi",      outcome = "women_iron")
  )
}

out_dir <- here::here("results", "sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

all_res <- list()
for (p in panel) {
  cc <- cfgs[[p$key]]
  if (is.null(cc)) { message("unknown config key: ", p$key,
                             " (use e.g. 'SierraLeone', not 'Sierra Leone')"); next }
  oc <- cc$outcomes[[p$outcome]]
  if (is.null(oc)) { message("no outcome ", p$outcome, " for ", p$key); next }
  od <- tryCatch(tar_read_raw(paste0("outcome_data_", tolower(p$key), "_", p$outcome)),
                 error = function(e) NULL)
  if (is.null(od)) { message("no cached outcome_data for ", p$key, " ", p$outcome,
                             " — build the targets pipeline first"); next }

  res <- aggregation_sweep_one(od, cc, oc)
  if (is.null(res) || !nrow(res)) next
  readr::write_csv(res, file.path(out_dir,
    sprintf("aggregation_level_national_%s_%s.csv", tolower(p$key), p$outcome)))
  cat(sprintf("\n[%s | %s] national-prevalence error (pp) by level:\n", cc$country, p$outcome))
  print(res[, c("level", "n_units", "obs_national_pct", "pred_national_pct",
                "abs_error_pp", "rank_spearman")], row.names = FALSE)
  all_res[[length(all_res) + 1L]] <- res
}

all_res <- dplyr::bind_rows(all_res)
if (nrow(all_res)) {
  readr::write_csv(all_res, file.path(out_dir, "aggregation_level_national_summary.csv"))
  cat("\n============ |national error| (pp) by aggregation level ============\n")
  print(as.data.frame(tidyr::pivot_wider(
    all_res[, c("country", "outcome", "level", "abs_error_pp")],
    names_from = level, values_from = abs_error_pp)), row.names = FALSE)
  cat("\nMean |error| and mean rank-r by level:\n")
  print(as.data.frame(all_res |> dplyr::group_by(level) |>
    dplyr::summarise(mean_abs_error_pp = mean(abs_error_pp, na.rm = TRUE),
                     mean_rank_r = mean(rank_spearman, na.rm = TRUE), .groups = "drop")),
    row.names = FALSE)
  cat(sprintf("\nWrote %s\n", file.path(out_dir, "aggregation_level_national_summary.csv")))
} else {
  cat("\nNo results — is the targets cache built? (outcome_data_* targets needed)\n")
}
