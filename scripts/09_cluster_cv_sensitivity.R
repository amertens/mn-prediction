# =============================================================================
# scripts/09_cluster_cv_sensitivity.R
#
# Refit one or more country x outcome SLs under the patched
# mlr3_SL_clustered (cluster_id passed via group=) and diff against the
# prior cv_perf snapshot (results/tables/cv_performance_all.csv).
#
# Does NOT touch the targets cache. Reads cached outcome_data targets,
# re-runs the fits in-process with the patched code, writes:
#   results/sensitivity/cluster_cv_<country_slug>_<outcome>.csv  (one per fit)
#   results/sensitivity/cluster_cv_summary.csv                   (combined)
#
# Usage:
#   # default: run a multi-outcome sensitivity panel
#   Rscript scripts/09_cluster_cv_sensitivity.R
#
#   # single fit
#   Rscript scripts/09_cluster_cv_sensitivity.R <Country> <outcome_tag>
#
# Country strings must match get_country_configs() names exactly
# (e.g. "Ghana", "Malawi", "Sierra Leone"). Outcome tags match
# cc$outcomes (e.g. "women_iron", "child_iron", "women_b12").
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(targets); library(dplyr); library(readr)
  library(mlr3superlearner); library(mlr3); library(mlr3learners)
  library(mlr3extralearners); library(caret); library(ck37r)
  library(recipes); library(data.table); library(origami); library(pROC)
})

source(here::here("R", "config.R"))
source(here::here("R", "mlr3_fitting.R"))
source(here::here("R", "admin1_analysis.R"))   # extract_cv_performance
source(here::here("R", "diagnostics.R"))

# ── Args ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  panel <- list(list(country = args[[1]], outcome = args[[2]]))
} else {
  # Default multi-outcome panel (chosen to span sample sizes and prevalences)
  panel <- list(
    list(country = "Ghana",        outcome = "women_iron"),
    list(country = "Malawi",       outcome = "women_iron"),
    list(country = "Ghana",        outcome = "women_b12"),
    list(country = "Sierra Leone", outcome = "child_iron")
  )
}

cfgs     <- get_country_configs()
sl       <- tar_read(sl_learners)
params   <- tar_read(pipeline_params)
baseline <- readr::read_csv(
  here::here("results", "tables", "cv_performance_all.csv"),
  show_col_types = FALSE)

country_slug <- function(country) tolower(gsub("[^a-z]+", "", tolower(country)))

run_one <- function(country, outcome) {
  slug <- country_slug(country)
  cat(sprintf("\n=========================================================\n"))
  cat(sprintf("  [%s] %s\n", country, outcome))
  cat(sprintf("=========================================================\n"))

  cc <- cfgs[[country]]
  oc <- cc$outcomes[[outcome]]

  od_target <- paste0("outcome_data_", slug, "_", outcome)
  od <- tar_read_raw(od_target)

  cat(sprintf("Outcome data: n=%d, %d X (proxy-only)\n",
              nrow(od$d), length(od$Xvars)))

  set.seed(params$seed %||% 12345)
  t0  <- Sys.time()
  fit <- fit_mlr3_models(od, cc, oc, sl, params)
  t1  <- Sys.time()
  cat(sprintf("Fit %.1f min\n", as.numeric(t1 - t0, units = "mins")))

  cont <- extract_cv_performance(fit$cont_fit, oc, "continuous")
  bin  <- extract_cv_performance(fit$bin_fit,  oc, "binary")
  new_perf <- dplyr::bind_rows(cont, bin) %>%
    dplyr::mutate(country = cc$country)

  base <- baseline %>%
    dplyr::filter(country == cc$country, outcome == oc$tag)

  metrics <- c("auc", "brier", "rmse", "mae", "r2", "r")
  jcols   <- c("country", "outcome", "model_type")
  diff_df <- new_perf %>%
    dplyr::select(dplyr::all_of(c(jcols, metrics))) %>%
    dplyr::rename_with(~ paste0(.x, "_patched"),
                       .cols = dplyr::all_of(metrics)) %>%
    dplyr::left_join(
      base %>%
        dplyr::select(dplyr::all_of(c(jcols, metrics))) %>%
        dplyr::rename_with(~ paste0(.x, "_baseline"),
                           .cols = dplyr::all_of(metrics)),
      by = jcols)
  for (m in metrics) {
    diff_df[[paste0(m, "_delta")]] <-
      diff_df[[paste0(m, "_patched")]] - diff_df[[paste0(m, "_baseline")]]
  }

  out_dir <- here::here("results", "sensitivity")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir,
                        sprintf("cluster_cv_%s_%s.csv", slug, outcome))
  readr::write_csv(diff_df, out_path)
  cat(sprintf("Wrote %s\n", out_path))
  print(diff_df %>% dplyr::select(model_type,
                                  dplyr::ends_with("_delta")))
  diff_df
}

# ── Run panel ─────────────────────────────────────────────────────────────────
all_diffs <- lapply(panel, function(p) {
  tryCatch(run_one(p$country, p$outcome),
           error = function(e) {
             cat(sprintf("\n[ERROR %s|%s] %s\n",
                         p$country, p$outcome, conditionMessage(e)))
             NULL
           })
})
all_diffs <- dplyr::bind_rows(all_diffs)

if (nrow(all_diffs) > 0 && length(panel) > 1) {
  summary_path <- here::here("results", "sensitivity",
                             "cluster_cv_summary.csv")
  readr::write_csv(all_diffs, summary_path)
  cat(sprintf("\nWrote combined summary: %s\n", summary_path))

  cat("\n=== Panel summary ===\n")
  cat("AUC deltas (binary):\n")
  print(all_diffs %>%
          dplyr::filter(model_type == "binary") %>%
          dplyr::select(country, outcome, auc_baseline, auc_patched, auc_delta))
  cat("\nRMSE deltas (continuous):\n")
  print(all_diffs %>%
          dplyr::filter(model_type == "continuous") %>%
          dplyr::select(country, outcome, rmse_baseline, rmse_patched, rmse_delta))

  auc_deltas <- all_diffs %>%
    dplyr::filter(model_type == "binary", !is.na(auc_delta)) %>%
    dplyr::pull(auc_delta)
  if (length(auc_deltas) > 0) {
    cat(sprintf("\nAUC delta: mean=%+.3f, min=%+.3f, max=%+.3f (n=%d outcomes)\n",
                mean(auc_deltas), min(auc_deltas), max(auc_deltas),
                length(auc_deltas)))
  }
}

cat("\nDone.\n")
