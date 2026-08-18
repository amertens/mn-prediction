# =============================================================================
# _targets.R
#
# Gambia (+ future countries) Micronutrient Deficiency Prediction Pipeline
#
# Uses the {targets} package for dependency-aware, cached execution.
# Only steps whose inputs changed will rerun.
#
# TWO MODES:
#
#   FAST (default) — ~5-10 min (parallel). Minimal SL stack (mean +
#     lasso + ranger), conformal CIs. Good for development, debugging,
#     checking relative performance across outcomes.
#
#   FULL — ~30-60 min (parallel). Evidence-based 6-learner stack
#     (mean + lasso + 2×ranger + lasso→xgb + xgb), conformal CIs.
#     Publication-quality results.
#
# UNCERTAINTY:
#   Conformal prediction intervals (not bootstrap) — uses CV residuals
#   to construct distribution-free intervals. No model refitting needed.
#
# PARALLELISM:
#   The pipeline parallelises at the targets level:
#   1. Targets-level: independent outcomes run concurrently (4 workers default)
#   3. Domain ablation: domain fits run in parallel via future.apply
#   Set TARGETS_WORKERS env var to control targets-level workers.
#   Within-target workers auto-scale based on available cores.
#   On Windows, uses multisession (socket clusters); needs ~2-4 GB RAM/worker.
#
# Usage:
#   # Fast mode with parallelism (default):
#   targets::tar_make_future(workers = 4)
#
#   # Full mode with parallelism:
#   Sys.setenv(PIPELINE_MODE = "full")
#   targets::tar_make_future(workers = 4)
#
#   # Or from the command line:
#   PIPELINE_MODE=full Rscript -e 'targets::tar_make_future(workers = 4)'
#
#   # Sequential (no parallelism, lower memory):
#   targets::tar_make()
#
#   # Switch modes (invalidates params + everything downstream):
#   Sys.setenv(PIPELINE_MODE = "full")
#   targets::tar_invalidate(names = "pipeline_params")
#   targets::tar_make()
#
#   # Other useful commands:
#   targets::tar_visnetwork()    # visualize the dependency graph
#   targets::tar_read(cv_perf)   # read a cached result
#
# Adding a new country:
#   1. Add an entry to get_country_configs() in R/config.R
#   2. Run tar_make() — new country targets are auto-generated
#
# Adding a new outcome:
#   1. Add to the outcomes list in the country config
#   2. Run tar_make() — new outcome targets are auto-generated
# =============================================================================

library(targets)
library(future)

# ── Parallel execution setup ─────────────────────────────────────────────────
# Targets-level parallelism: run independent outcome targets concurrently.
# On Windows, must use multisession (not multicore).
#
# N_WORKERS controls how many targets run in parallel. With 4 outcomes per
# country, 4 workers keeps all outcomes busy simultaneously. Increase if
# running multiple countries. Set via env var to override:
#   TARGETS_WORKERS=8 Rscript -e 'targets::tar_make_future(workers = 8)'
n_workers <- as.integer(Sys.getenv("TARGETS_WORKERS", "4"))

# Set the future plan for tar_make_future(). On Windows, multisession spawns
# socket-based R processes. Each needs ~2-4 GB for SL objects.
future::plan(future::multisession, workers = n_workers)
options(future.globals.maxSize = 3 * 1024^3)  # 3 GB — SL objects are large

# ── Note: Bootstrap replaced by conformal prediction intervals ────────────
# Conformal intervals use existing CV residuals — no refitting needed.
# The old bootstrap code is preserved in R/bootstrap.R for reference.

# ── Packages available to all targets ────────────────────────────────────────
tar_option_set(
  packages = c(
    "here", "dplyr", "tidyr", "readr", "tibble", "rlang",
    "ggplot2", "sf", "geodata", "terra", "exactextractr",
    "srvyr", "survey", "scales", "viridis", "patchwork", "ggrepel",
    "mlr3", "mlr3learners", "mlr3extralearners", "mlr3superlearner",
    "origami", "caret", "data.table", "ck37r", "labelled",
    "recipes", "future.apply", "glmnet", "pROC", "ROCR", "haven", "readxl",
    # used by the folded-in corrected-methods learners (R/corrected/)
    "ranger", "rpart",
    # unit-level SAE (nested-error GLMM) in the protocol-reconciliation target
    "lme4"
  ),
  # Increase memory limit for SL fitting
  memory = "transient",
  garbage_collection = TRUE,
  # Error handling: workspace for debugging
  workspace_on_error = TRUE
)

# ── Source all function files in R/ ──────────────────────────────────────────
tar_source("R/")

# Also source legacy helper files that define learner objects, etc.
# (These are needed by sl_helpers.R / DHS_SL_clustered)
source(here::here("src", "0-functions.R"))

# ── Display mode ────────────────────────────────────────────────────────────
pipeline_mode <- Sys.getenv("PIPELINE_MODE", "fast")
message(sprintf(
  paste0("\n=== Pipeline mode: %s | %d parallel workers ===\n",
         "  %s\n  Set PIPELINE_MODE=%s to switch.\n",
         "  Run with: targets::tar_make_future(workers = %d)\n"),
  toupper(pipeline_mode),
  n_workers,
  if (pipeline_mode == "fast")
    "Minimal SL stack (3 learners), conformal CIs. ~5-10 min with parallelism."
  else
    "Evidence-based 6-learner stack, conformal CIs. ~30-60 min with parallelism.",
  if (pipeline_mode == "fast") "full" else "fast",
  n_workers
))


# =============================================================================
# DYNAMIC TARGET FACTORY  —  SENSITIVITY ANALYSIS (individual-level SuperLearner)
#
# For each country x outcome, generates the full set of individual-level
# SuperLearner targets (sl_fit_*, conformal CIs, ablations, diagnostics, SHAP,
# national estimates, and the Admin-2 aggregation of person-level predictions).
#
# NOTE ON ROLE (2026-06-15): the individual-level SL pipeline is now a
# SENSITIVITY analysis, not the primary one. The PRIMARY analysis is the
# Admin-2 area-level small-area estimation (SAE) — see `benchmark_targets`
# (area-level SuperLearner, Fay-Herriot, BYM2) further below. The individual-SL
# model code lives in R/sensitivity/{sl_fitting,mlr3_fitting}.R; this factory is
# retained so its outputs remain available for comparison, but headline claims
# and the dashboard default to the area-level SAE predictions.
# =============================================================================

#' Generate targets for one country x one outcome combination
make_outcome_targets <- function(country_name, outcome_name, cc, oc, params) {

  # Create unique target name suffixes: e.g., "gambia_child_vitA"
  suffix <- paste0(tolower(country_name), "_", outcome_name)

  list(
    # ── 1. Build outcome-specific dataset ──────────────────────────────────
    tar_target_raw(
      paste0("outcome_data_", suffix),
      substitute(
        build_outcome_dataset(merged_data, cc_val, oc_val),
        list(
          merged_data = as.symbol(paste0("merged_", tolower(country_name))),
          cc_val = cc,
          oc_val = oc
        )
      )
    ),

    # ── 2. Fit SuperLearner models (EXPENSIVE — cached) ───────────────────
    tar_target_raw(
      paste0("sl_fit_", suffix),
      substitute(
        fit_mlr3_models(outcome_data, cc_val, oc_val, sl_learners, params_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          cc_val       = cc,
          oc_val       = oc,
          sl_learners  = as.symbol("sl_learners"),
          params_val   = params
        )
      )
    ),

    # ── 3. CV performance metrics ─────────────────────────────────────────
    tar_target_raw(
      paste0("cv_perf_", suffix),
      substitute({
        cont <- extract_cv_performance(sl_fit$cont_fit, oc_val, "continuous")
        bin  <- extract_cv_performance(sl_fit$bin_fit, oc_val, "binary")
        dplyr::bind_rows(cont, bin) %>%
          dplyr::mutate(country = cc_val$country)
      },
      list(
        sl_fit  = as.symbol(paste0("sl_fit_", suffix)),
        oc_val  = oc,
        cc_val  = cc
      ))
    ),

    # ── 4. Admin1 aggregation ─────────────────────────────────────────────
    tar_target_raw(
      paste0("admin1_prev_", suffix),
      substitute({
        fit <- sl_fit$bin_fit %||% sl_fit$cont_fit
        model_type <- if (!is.null(sl_fit$bin_fit)) "binary_prob" else "cont_threshold"
        aggregate_admin1(fit, outcome_data, cc_val, oc_val, model_type)
      },
      list(
        sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
        outcome_data = as.symbol(paste0("outcome_data_", suffix)),
        cc_val       = cc,
        oc_val       = oc
      ))
    ),

    # ── 5. Conformal prediction intervals (FAST — uses existing CV preds) ──
    # Replaces the expensive bootstrap (B × SL refits) with distribution-free
    # conformal intervals derived from the cross-validated residuals.
    # Runs for ALL country×outcome combos since it's essentially free.
    tar_target_raw(
      paste0("conformal_ci_", suffix),
      substitute(
        compute_conformal_ci(outcome_data, sl_fit, cc_val, oc_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    ),

    # ── 6. Domain ablation (EXPENSIVE — cached) ───────────────────────────
    tar_target_raw(
      paste0("ablation_", suffix),
      substitute(
        run_domain_ablation(outcome_data, sl_fit, cc_val, oc_val, sl_learners, params_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          cc_val       = cc,
          oc_val       = oc,
          sl_learners  = as.symbol("sl_learners"),
          params_val   = params
        )
      )
    ),

    # ── 7. Admin2 individual-level SL aggregation ─────────────────────────
    tar_target_raw(
      paste0("admin2_sl_", suffix),
      substitute(
        aggregate_admin2_sl(sl_fit, outcome_data, cc_val, oc_val),
        list(
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    ),

    # ── 8. Admin2 survey-weighted prevalence ──────────────────────────────
    tar_target_raw(
      paste0("svy_admin2_", suffix),
      substitute(
        compute_svy_admin2(outcome_data, cc_val, oc_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    ),

    # ── 9. Admin2 error analysis ──────────────────────────────────────────
    tar_target_raw(
      paste0("admin2_error_", suffix),
      substitute(
        compute_admin2_error(admin2_sl, svy_admin2, cc_val, oc_val),
        list(
          admin2_sl  = as.symbol(paste0("admin2_sl_", suffix)),
          svy_admin2 = as.symbol(paste0("svy_admin2_", suffix)),
          cc_val     = cc,
          oc_val     = oc
        )
      )
    ),

    # ── 10. Plots (Admin-1 scatter — runs for all combos) ────────────────
    tar_target_raw(
      paste0("plot_admin1_scatter_", suffix),
      substitute(
        plot_admin1_scatter(admin1_prev, oc_val, cc_val),
        list(
          admin1_prev = as.symbol(paste0("admin1_prev_", suffix)),
          oc_val      = oc,
          cc_val      = cc
        )
      )
    ),

    # ── 12. Post-hoc diagnostics (cheap — uses existing predictions) ────
    tar_target_raw(
      paste0("diagnostics_", suffix),
      substitute(
        run_diagnostics(sl_fit, cc_val, oc_val),
        list(
          sl_fit = as.symbol(paste0("sl_fit_", suffix)),
          cc_val = cc,
          oc_val = oc
        )
      )
    ),

    # ── 12b. Calibrated diagnostics: Platt recalibration of OOF binary
    # predictions before computing Brier-skill / calibration-slope metrics.
    # Repairs the ~25% of binary models with negative Brier skill identified
    # in diagnostics_binary.csv. Cheap — reuses the same sl_fit object.
    tar_target_raw(
      paste0("diagnostics_calibrated_", suffix),
      substitute(
        run_diagnostics_calibrated(sl_fit, cc_val, oc_val),
        list(
          sl_fit = as.symbol(paste0("sl_fit_", suffix)),
          cc_val = cc,
          oc_val = oc
        )
      )
    ),

    # ── 13. National estimates: survey-weighted observed vs predicted ─────
    tar_target_raw(
      paste0("national_est_", suffix),
      substitute(
        compute_national_estimates(outcome_data, sl_fit, cc_val, oc_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    ),

    # ── 14. SHAP feature attributions (district top factors + global imp) ─
    tar_target_raw(
      paste0("shap_", suffix),
      substitute(
        compute_shap_explanations(outcome_data, sl_fit, cc_val, oc_val),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    ),

    # ── 15. Single-variable permutation importance (top-30 vars) ──────────
    tar_target_raw(
      paste0("varimp_", suffix),
      substitute(
        run_single_var_ablation(outcome_data, sl_fit, cc_val, oc_val,
                                 top_n = 30L, n_perm = 3L),
        list(
          outcome_data = as.symbol(paste0("outcome_data_", suffix)),
          sl_fit       = as.symbol(paste0("sl_fit_", suffix)),
          cc_val       = cc,
          oc_val       = oc
        )
      )
    )
  )
}


# =============================================================================
# BUILD THE TARGET LIST
# =============================================================================

# Get configs
all_country_configs <- get_country_configs()
params              <- get_pipeline_params(pipeline_mode)

# ── Static targets (shared across all outcomes) ─────────────────────────────
static_targets <- list(

  # Pipeline parameters — depends on mode string so switching modes
  # invalidates all downstream targets automatically
  tar_target(pipeline_params, get_pipeline_params(pipeline_mode)),

  # SL learner stacks (shared across all outcomes)
  tar_target(sl_learners, setup_mlr3_learners(pipeline_params))
)

# ── Per-country targets ─────────────────────────────────────────────────────
country_targets <- list()

# Food security data paths (shared across countries)
hfid_path <- here::here("data", "raw", "HFID", "hfid_hv1.csv")
ch_path   <- here::here("data", "raw", "CadreHarmonise",
                         "cadre_harmonise_caf_ipc_dec25.xlsx")

for (country_name in names(all_country_configs)) {
  cc <- all_country_configs[[country_name]]
  lc <- tolower(country_name)

  # ── Track the merged dataset FILE so on-disk changes invalidate downstream ──
  # A bare path string is NOT tracked by {targets} — only a format = "file"
  # target hashes the file's contents. Without this, edits to the source
  # *_merged_dataset.rds (e.g. the 2026-05-13 admin-1 IHME addition) are
  # silently ignored and the pipeline keeps serving a stale cached snapshot.
  path_target_name <- paste0("path_merged_", lc)
  country_targets <- c(country_targets, list(
    tar_target_raw(
      path_target_name,
      substitute(p, list(p = cc$data_path)),
      format = "file"
    )
  ))

  # Load raw merged data. Depends on the file target above, so it reloads
  # whenever the underlying .rds changes on disk.
  raw_target_name <- paste0("merged_raw_", lc)
  raw_data_target <- tar_target_raw(
    raw_target_name,
    substitute(
      load_merged_data(data_path_val),
      list(data_path_val = as.symbol(path_target_name))
    ),
    format = "rds"
  )
  country_targets <- c(country_targets, list(raw_data_target))

  # Merge food security data (HFID + Cadre Harmonisé)
  fsec_target_name <- paste0("merged_fsec_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      fsec_target_name,
      substitute(
        merge_food_security(raw_data, cc_val, hfid_val, ch_val),
        list(
          raw_data  = as.symbol(raw_target_name),
          cc_val    = cc,
          hfid_val  = hfid_path,
          ch_val    = ch_path
        )
      )
    )
  ))

  # ── Track the legacy-parity GEE CSV so edits to it invalidate downstream ────
  # A country with no data/<Country>_GEE_rasters/ gets its Admin-2 covariates
  # from data/GEE/<ISO3>_legacy_parity_admin2_gee.csv instead (see
  # .append_legacy_parity_cols in R/admin2_analysis.R). That file is read inside
  # extract_gee_admin2(), so {targets} cannot see it — without this stamp, a
  # freshly-written or corrected parity CSV would be ignored and the pipeline
  # would keep serving a covariate-free cached snapshot. Same hazard the
  # path_merged_* file targets above guard against.
  #
  # format = "file" is not usable here: the four raster countries have no such
  # CSV and a missing file target is an error. Instead the stamp re-computes
  # every run (cheap md5 of a ~200 KB CSV) and only invalidates downstream when
  # the file's contents actually change.
  parity_stamp_name <- paste0("gee_parity_stamp_", lc)
  country_targets <- c(country_targets, list(
    tar_target_raw(
      parity_stamp_name,
      substitute({
        p <- legacy_parity_csv_path(cc_val)
        if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
      }, list(cc_val = cc)),
      cue = tar_cue(mode = "always")
    )
  ))

  # Extract GEE raster zonal means for Admin-2 polygons (once per country).
  # Depends on the parity stamp so a new/updated parity CSV re-triggers it.
  gee_admin2_target_name <- paste0("gee_admin2_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      gee_admin2_target_name,
      substitute({
        force(parity_stamp_val)
        extract_gee_admin2(cc_val)
      }, list(cc_val = cc, parity_stamp_val = as.symbol(parity_stamp_name)))
    )
  ))

  # ── Track the DHS clean inputs the external merge reads internally ─────────
  # merge_external_predictors() calls load_dhs_admin2(), which reads
  # data/DHS/clean/<country>_<year>_dhs_(custom_)admin2_wide.rds directly. Those
  # reads are invisible to {targets}, so regenerating them (re-running
  # src/DHS/DHS_admin2_aggregation.R) would otherwise leave merged_ext_* serving
  # a stale cached result. Same md5 + always-cue pattern as gee_parity_stamp_*.
  dhs_stamp_name <- paste0("dhs_stamp_", lc)
  country_targets <- c(country_targets, list(
    tar_target_raw(
      dhs_stamp_name,
      substitute(dhs_clean_stamp(cc_val), list(cc_val = cc)),
      cue = tar_cue(mode = "always")
    )
  ))

  # Merge external predictors (CHIRPS, WorldPop, MAP, HarvestStat, etc.)
  # AND GEE Admin-2 zonal means into individual-level data.
  merged_ext_target_name <- paste0("merged_ext_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      merged_ext_target_name,
      substitute({
        force(dhs_stamp_val)
        merge_external_predictors(fsec_data, cc_val, cache_dir_val)
      },
        list(
          fsec_data     = as.symbol(fsec_target_name),
          cc_val        = cc,
          cache_dir_val = here::here("data", "external_cache"),
          dhs_stamp_val = as.symbol(dhs_stamp_name)
        )
      )
    )
  ))

  # Merge GEE Admin-2 raster variables into individual records
  merged_target_name <- paste0("merged_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      merged_target_name,
      substitute({
        d <- ext_data
        gee <- gee_admin2_data
        if (!is.null(gee) && is.data.frame(gee) && nrow(gee) > 0) {
          # Only merge gee_ columns (not Admin1/Admin2 which already exist)
          gee_cols <- grep("^gee_", colnames(gee), value = TRUE)
          if (length(gee_cols) > 0) {
            # Remove any gee_ columns already in the data to avoid .x/.y suffixes
            existing_gee <- grep("^gee_", colnames(d), value = TRUE)
            new_gee <- setdiff(gee_cols, existing_gee)
            if (length(new_gee) > 0) {
              admin2_col <- cc_val$admin2_col
              n_before <- nrow(d)
              d <- merge(d, gee[, c("Admin2", new_gee), drop = FALSE],
                         by.x = admin2_col, by.y = "Admin2",
                         all.x = TRUE, sort = FALSE)
              if (nrow(d) != n_before) {
                warning(sprintf("[merge_gee] Row count changed: %d -> %d", n_before, nrow(d)))
              }
              cat(sprintf("[merge_gee] Added %d GEE Admin-2 raster columns to %s\n",
                          length(new_gee), cc_val$country))
            }
          }
        }
        d
      },
      list(
        ext_data        = as.symbol(merged_ext_target_name),
        gee_admin2_data = as.symbol(gee_admin2_target_name),
        cc_val          = cc
      ))
    )
  ))

  # Per-outcome targets
  for (outcome_name in names(cc$outcomes)) {
    oc <- cc$outcomes[[outcome_name]]
    outcome_targets <- make_outcome_targets(country_name, outcome_name, cc, oc, params)
    country_targets <- c(country_targets, outcome_targets)
  }
}

# ── Area-level model: ALL country × outcome combinations ──────────────────
# Expensive GEE raster extraction + area-level elastic net. Produces
# Admin-2 prevalence predictions for *all* polygons in the country
# (surveyed + unsurveyed), enabling full choropleth coverage in the
# dashboard. Previously this ran only for Ghana × women_iron as an exemplar.
for (cc_name_area in names(all_country_configs)) {
  cc_area_local <- all_country_configs[[cc_name_area]]

  # Per-country GADM download + GEE raster zonal extraction. Runs ONCE per
  # country and is reused by every area_model_<outcome> target below — the
  # extraction depends only on the country, not the outcome.
  area_cov_target <- paste0("area_covariates_", tolower(cc_name_area))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      area_cov_target,
      # Same untracked-input hazard as gee_admin2_*: extract_area_covariates()
      # reads the legacy-parity CSV internally for countries with no rasters, so
      # depend on that file's stamp explicitly.
      substitute({
        force(parity_stamp_val)
        extract_area_covariates(cc_val)
      }, list(cc_val = cc_area_local,
              parity_stamp_val = as.symbol(paste0("gee_parity_stamp_",
                                                  tolower(cc_name_area)))))
    )
  ))

  for (oc_name_area in names(cc_area_local$outcomes)) {
    oc_area_local <- cc_area_local$outcomes[[oc_name_area]]
    area_suffix <- paste0(tolower(cc_name_area), "_", oc_name_area)

    country_targets <- c(country_targets, list(
      tar_target_raw(
        paste0("area_model_", area_suffix),
        substitute(
          fit_area_level_model(svy_admin2, area_cov, cc_val, oc_val, params_val),
          list(
            svy_admin2 = as.symbol(paste0("svy_admin2_", area_suffix)),
            area_cov   = as.symbol(area_cov_target),
            cc_val     = cc_area_local,
            oc_val     = oc_area_local,
            params_val = params
          )
        )
      )
    ))

    # Plot targets for the original Ghana × women_iron exemplar only —
    # the dashboard renders all other combinations interactively.
    if (area_suffix == "ghana_women_iron") {
      country_targets <- c(country_targets, list(
        tar_target_raw(
          paste0("plot_admin2_coverage_", area_suffix),
          substitute(
            plot_admin2_coverage(area_model, oc_val, cc_val),
            list(
              area_model = as.symbol(paste0("area_model_", area_suffix)),
              oc_val     = oc_area_local,
              cc_val     = cc_area_local
            )
          )
        ),
        tar_target_raw(
          paste0("plot_admin2_forest_", area_suffix),
          substitute(
            plot_admin2_forest(area_model, oc_val, cc_val),
            list(
              area_model = as.symbol(paste0("area_model_", area_suffix)),
              oc_val     = oc_area_local,
              cc_val     = cc_area_local
            )
          )
        )
      ))
    }
  }
}

# ── Level-2 conceptual ablation: Ghana × women_iron exemplar only ──────────
# Finer-grained domain importance at Level-2 (e.g., Vaccinations, Precipitation,
# Cooking fuel). Single exemplar for the results document.
{
  cc_l2 <- all_country_configs[["Ghana"]]
  oc_l2 <- cc_l2$outcomes[["women_iron"]]

  country_targets <- c(country_targets, list(
    tar_target_raw(
      "ablation_l2_ghana_women_iron",
      substitute(
        run_conceptual_permutation(sl_fit, outcome_data, oc_val,
                                    level = "level2", n_perm = 5),
        list(
          sl_fit       = as.symbol("sl_fit_ghana_women_iron"),
          outcome_data = as.symbol("outcome_data_ghana_women_iron"),
          oc_val       = oc_l2
        )
      )
    )
  ))
}


# ── Area-level comparison targets ──────────────────────────────────────────
# Compare individual-level SL (aggregated to Admin-2) vs area-level SL
# (directly predicting Admin-2 prevalence from GEE covariates) using
# both MSE and NLL (beta-inspired) loss functions.

area_comparison_targets <- list()
for (cc_name in names(all_country_configs)) {
  cc_local <- all_country_configs[[cc_name]]
  ctry_lower <- tolower(cc_name)

  for (oc_name in names(cc_local$outcomes)) {
    oc_local <- cc_local$outcomes[[oc_name]]
    suffix <- paste0(ctry_lower, "_", oc_name)

    area_comparison_targets <- c(area_comparison_targets, list(
      tar_target_raw(
        paste0("area_comparison_", suffix),
        substitute(
          run_area_comparison(svy_admin2_val, sl_admin2_val, gee_admin2_val,
                              cc_val, oc_val),
          list(
            svy_admin2_val = as.symbol(paste0("svy_admin2_", suffix)),
            sl_admin2_val  = as.symbol(paste0("admin2_sl_", suffix)),
            gee_admin2_val = as.symbol(paste0("gee_admin2_", ctry_lower)),
            cc_val = cc_local,
            oc_val = oc_local
          )
        )
      )
    ))
  }
}

# Combined area comparison summary
area_comp_names <- paste0("area_comparison_",
  unlist(lapply(names(all_country_configs), function(cc_name) {
    paste0(tolower(cc_name), "_", names(all_country_configs[[cc_name]]$outcomes))
  })))
area_comp_syms <- lapply(area_comp_names, as.symbol)
area_comp_list_expr <- as.call(c(list(as.symbol("list")), setNames(area_comp_syms, area_comp_names)))

area_comparison_targets <- c(area_comparison_targets, list(
  tar_target_raw(
    "area_comparison_all",
    substitute({
      comps <- all_comps_val
      rows <- lapply(comps, function(x) x$comparison)
      rows <- Filter(Negate(is.null), rows)
      result <- dplyr::bind_rows(rows)
      write.csv(result, here::here("results", "tables", "area_comparison_all.csv"),
                row.names = FALSE)
      result
    }, list(all_comps_val = area_comp_list_expr))
  )
))

# ── Area-level LOCO targets ──────────────────────────────────────────────
# Area-level LOCO for the 4 shared outcomes
area_loco_targets <- list()
if (length(all_country_configs) >= 2) {
  shared_outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

  for (otag in shared_outcomes) {
    # Build area LOCO dataset: pool svy_admin2 and gee_admin2 across countries
    loco_svy_syms <- lapply(names(all_country_configs), function(cc_name) {
      as.symbol(paste0("svy_admin2_", tolower(cc_name), "_", otag))
    })
    names(loco_svy_syms) <- names(all_country_configs)
    loco_svy_expr <- as.call(c(list(as.symbol("list")), loco_svy_syms))

    loco_gee_syms <- lapply(names(all_country_configs), function(cc_name) {
      as.symbol(paste0("gee_admin2_", tolower(cc_name)))
    })
    names(loco_gee_syms) <- names(all_country_configs)
    loco_gee_expr <- as.call(c(list(as.symbol("list")), loco_gee_syms))

    area_loco_targets <- c(area_loco_targets, list(
      tar_target_raw(
        paste0("area_loco_", otag),
        substitute({
          pooled <- build_area_loco_dataset(svy_list_val, gee_list_val)
          run_area_loco(pooled)
        }, list(
          svy_list_val = loco_svy_expr,
          gee_list_val = loco_gee_expr
        ))
      )
    ))
  }

  # Combined area LOCO summary
  area_loco_names <- paste0("area_loco_", shared_outcomes)
  area_loco_syms <- setNames(lapply(area_loco_names, as.symbol), shared_outcomes)
  area_loco_list_expr <- as.call(c(list(as.symbol("list")), area_loco_syms))

  area_loco_targets <- c(area_loco_targets, list(
    tar_target_raw(
      "area_loco_comparison",
      substitute({
        all_loco <- loco_list_val
        rows <- list()
        for (oc in names(all_loco)) {
          if (is.data.frame(all_loco[[oc]]) && nrow(all_loco[[oc]]) > 0) {
            all_loco[[oc]]$outcome <- oc
            rows[[oc]] <- all_loco[[oc]]
          }
        }
        result <- dplyr::bind_rows(rows)
        write.csv(result, here::here("results", "tables", "area_loco_comparison.csv"),
                  row.names = FALSE)
        result
      }, list(loco_list_val = area_loco_list_expr))
    )
  ))
}


# ── PRIMARY ANALYSIS: Admin-2 area-level small-area estimation (SAE) ───────
# Benchmark suite — the PRIMARY analysis of the project. Compares the
# area-level SuperLearner against (1) country-mean baseline, (2) OLS GLM,
# (3) elastic-net penalised GLM, (4) Fay-Herriot, (5) BYM2 (INLA). Runs both
# within-country k-fold CV and leave-one-country-out (LOCO) transportability.
# Implementation in R/benchmark_models.R. These area-level SAE predictions —
# not the individual-level SuperLearner factory above — are the headline
# results and the dashboard default. The individual-level SL pipeline is a
# sensitivity analysis (see the DYNAMIC TARGET FACTORY comment above).
#
# Outputs:
#   - area_adjacency_list : country -> spdep::nb (built once, reused).
#   - benchmarks_<outcome>: per-outcome data.frame with one row per
#                           (method, holdout, model_type) eval combo.
#   - benchmarks_all      : stacked rollup with CSV save for the dashboard.
benchmark_targets <- list()
if (length(all_country_configs) >= 2) {
  bench_outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
  cc_lower_b <- tolower(names(all_country_configs))
  cc_label_b <- names(all_country_configs)  # display names (Gambia, ...)

  gee_syms_b <- setNames(lapply(paste0("gee_admin2_", cc_lower_b), as.symbol),
                          cc_label_b)
  gee_expr_b <- as.call(c(list(as.symbol("list")), gee_syms_b))

  # Adjacency built once across all training countries (used by BYM2 LOCO
  # and within-country runs). Depends on configs + at least one svy_admin2
  # frame per country so the GADM polygons align with the admin-2 ordering.
  any_svy_syms_b <- setNames(
    lapply(cc_lower_b, function(c) as.symbol(paste0("svy_admin2_", c, "_child_vitA"))),
    cc_label_b)
  any_svy_expr_b <- as.call(c(list(as.symbol("list")), any_svy_syms_b))

  benchmark_targets <- c(benchmark_targets, list(
    tar_target_raw(
      "area_adjacency_list",
      substitute({
        build_adjacency_list(get_country_configs(), any_svy_val)
      }, list(any_svy_val = any_svy_expr_b))
    )
  ))

  # Per-outcome benchmark target.
  for (otag in bench_outcomes) {
    svy_syms_b <- setNames(
      lapply(cc_lower_b, function(c) as.symbol(paste0("svy_admin2_", c, "_", otag))),
      cc_label_b)
    svy_expr_b <- as.call(c(list(as.symbol("list")), svy_syms_b))

    benchmark_targets <- c(benchmark_targets, list(
      tar_target_raw(
        paste0("benchmarks_", otag),
        substitute({
          pooled <- build_area_loco_dataset(svy_list_val, gee_list_val)
          # Run LOCO only (within-country requires a per-country adjacency
          # slice, which is already available via area_adjacency_list).
          out <- run_area_benchmarks_loco(
            pooled_data    = pooled$pooled_data,
            gee_vars       = pooled$common_gee_vars,
            country_names  = pooled$country_names,
            adjacency_list = adj_val,
            outcome_label  = otag_val)
          if (nrow(out) > 0) out$eval_type <- "loco"
          out
        }, list(
          svy_list_val = svy_expr_b,
          gee_list_val = gee_expr_b,
          adj_val      = as.symbol("area_adjacency_list"),
          otag_val     = otag
        ))
      )
    ))
  }

  # Combined benchmarks rollup with CSV save.
  bench_names <- paste0("benchmarks_", bench_outcomes)
  bench_syms <- setNames(lapply(bench_names, as.symbol), bench_outcomes)
  bench_list_expr <- as.call(c(list(as.symbol("list")), bench_syms))

  benchmark_targets <- c(benchmark_targets, list(
    tar_target_raw(
      "benchmarks_all",
      substitute({
        rows <- list()
        all_b <- bench_list_val
        for (oc in names(all_b)) {
          if (is.data.frame(all_b[[oc]]) && nrow(all_b[[oc]]) > 0) {
            all_b[[oc]]$outcome <- oc
            rows[[oc]] <- all_b[[oc]]
          }
        }
        result <- dplyr::bind_rows(rows)
        out_dir <- here::here("results", "tables")
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(result, file.path(out_dir, "benchmarks_all.csv"),
                  row.names = FALSE)
        result
      }, list(bench_list_val = bench_list_expr))
    )
  ))

  # ── Fast primary SL run: ONE method (sl_prescreened) per outcome ──────────
  # Optimized main-analysis workflow. fit_predict_sl_prescreened uses
  # five-stage prescreening + spatial coords + a fast learner library
  # (no BART, no deep ranger). Per-outcome wall time ~30s per LOCO fold
  # vs. the historical 17-hour SL.
  #
  # The 4 main outcomes are in every training country; women_b12 and
  # women_folate are in 3 of 4 (Gambia missing), so we wrap their
  # svy_admin2 list in tryCatch and let build_area_loco_dataset drop
  # missing countries.
  sl_outcomes <- c(bench_outcomes, "women_b12", "women_folate")
  # Hardcode known country-availability: women_b12 / women_folate are absent
  # from Gambia. Other outcomes are 4-country. Determined from the targets
  # manifest at session start (archive/sandbox/00_setup.R load_pooled() also handles
  # this case).
  outcome_countries <- list(
    child_vitA   = cc_lower_b,
    women_vitA   = cc_lower_b,
    child_iron   = cc_lower_b,
    women_iron   = cc_lower_b,
    women_b12    = setdiff(cc_lower_b, "gambia"),
    women_folate = setdiff(cc_lower_b, "gambia")
  )
  for (otag in sl_outcomes) {
    avail <- outcome_countries[[otag]]
    if (length(avail) < 2) next  # Need >=2 countries for LOCO.
    avail_labels <- cc_label_b[match(avail, cc_lower_b)]
    sl_svy_syms <- setNames(
      lapply(avail, function(c) as.symbol(paste0("svy_admin2_", c, "_", otag))),
      avail_labels)
    sl_svy_expr <- as.call(c(list(as.symbol("list")), sl_svy_syms))
    sl_gee_syms <- setNames(
      lapply(avail, function(c) as.symbol(paste0("gee_admin2_", c))),
      avail_labels)
    sl_gee_expr <- as.call(c(list(as.symbol("list")), sl_gee_syms))

    benchmark_targets <- c(benchmark_targets, list(
      tar_target_raw(
        paste0("sl_prescreened_", otag),
        substitute({
          pooled <- build_area_loco_dataset(svy_list_val, gee_list_val)
          # Add admin-2 centroids so the prescreened SL can use (lon, lat).
          pooled$pooled_data <- add_admin2_centroids(
            pooled$pooled_data, get_country_configs(), svy_list_val)
          out <- run_area_benchmarks_loco(
            pooled_data    = pooled$pooled_data,
            gee_vars       = pooled$common_gee_vars,
            country_names  = pooled$country_names,
            adjacency_list = NULL,            # not needed for sl_prescreened
            outcome_label  = otag_val,
            methods        = "sl_prescreened",
            augment_features = FALSE)
          if (nrow(out) > 0) {
            out$eval_type <- "loco_sl_only"
            out$outcome   <- otag_val
          }
          out
        }, list(
          svy_list_val = sl_svy_expr,
          gee_list_val = sl_gee_expr,
          otag_val     = otag
        ))
      )
    ))
  }

  # Rollup: combined CSV consumed by dashboard prep.
  sl_pre_names <- paste0("sl_prescreened_", sl_outcomes)
  sl_pre_syms  <- setNames(lapply(sl_pre_names, as.symbol), sl_outcomes)
  sl_pre_expr  <- as.call(c(list(as.symbol("list")), sl_pre_syms))
  benchmark_targets <- c(benchmark_targets, list(
    tar_target_raw(
      "sl_prescreened_all",
      substitute({
        rows <- list()
        all_pre <- pre_list_val
        for (oc in names(all_pre)) {
          if (is.data.frame(all_pre[[oc]]) && nrow(all_pre[[oc]]) > 0) {
            rows[[oc]] <- all_pre[[oc]]
          }
        }
        result <- dplyr::bind_rows(rows)
        out_dir <- here::here("results", "tables")
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(result,
                   file.path(out_dir, "sl_prescreened_main.csv"),
                   row.names = FALSE)
        result
      }, list(pre_list_val = sl_pre_expr))
    )
  ))
}


# ── Enriched area-level transportability targets ──────────────────────────
# Universal / parsimonious within-country-centered elastic-net LOCO using
# harmonized GEE + IHME + Malaria-Atlas + food-security proxies aggregated to
# Admin-2. Functions live in R/transportability_area.R. Produces both LOCO
# metrics and per-Admin-2 out-of-sample predictions (the latter power the
# dashboard "transportability error" difference map).
area_transport_targets <- list()
if (length(all_country_configs) >= 2) {
  shared_outcomes_t <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
  cc_lower_t <- tolower(names(all_country_configs))

  gee_syms_t <- setNames(lapply(paste0("gee_admin2_", cc_lower_t), as.symbol), cc_lower_t)
  gee_expr_t <- as.call(c(list(as.symbol("list")), gee_syms_t))
  merged_syms_t <- setNames(lapply(paste0("merged_", cc_lower_t), as.symbol), cc_lower_t)
  merged_expr_t <- as.call(c(list(as.symbol("list")), merged_syms_t))

  # Enriched Admin-2 covariate tables (one per country)
  area_transport_targets <- c(area_transport_targets, list(
    tar_target_raw(
      "area_transport_covariates",
      substitute({
        gl <- gee_list_val; ml <- merged_list_val
        cov <- lapply(names(gl), function(cn) build_admin2_covariates(gl[[cn]], ml[[cn]]))
        names(cov) <- names(gl)
        cov
      }, list(gee_list_val = gee_expr_t, merged_list_val = merged_expr_t))
    )
  ))

  for (otag in shared_outcomes_t) {
    svy_syms_t <- setNames(
      lapply(cc_lower_t, function(cn) as.symbol(paste0("svy_admin2_", cn, "_", otag))),
      cc_lower_t)
    svy_expr_t <- as.call(c(list(as.symbol("list")), svy_syms_t))
    area_transport_targets <- c(area_transport_targets, list(
      tar_target_raw(
        paste0("area_transport_", otag),
        substitute({
          pooled <- assemble_area_transport(svy_list_val, area_transport_covariates, otag_val)
          if (is.null(pooled)) NULL else run_area_transport_loco(pooled, AREA_TRANSPORT_RECIPE)
        }, list(svy_list_val = svy_expr_t, otag_val = otag))
      )
    ))
  }

  at_names <- paste0("area_transport_", shared_outcomes_t)
  at_syms  <- setNames(lapply(at_names, as.symbol), shared_outcomes_t)
  at_list_expr <- as.call(c(list(as.symbol("list")), at_syms))
  area_transport_targets <- c(area_transport_targets, list(
    tar_target_raw(
      "area_transport_summary",
      substitute({
        res <- res_list_val
        mets  <- dplyr::bind_rows(lapply(res, function(r) if (!is.null(r)) r$metrics))
        preds <- dplyr::bind_rows(lapply(res, function(r) if (!is.null(r)) r$predictions))
        dir.create(here::here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
        dir.create(here::here("results", "transportability"), showWarnings = FALSE, recursive = TRUE)
        write.csv(mets,  here::here("results", "tables", "transportability_area_loco_metrics.csv"), row.names = FALSE)
        write.csv(preds, here::here("results", "tables", "transportability_area_loco_predictions.csv"), row.names = FALSE)
        saveRDS(preds, here::here("results", "transportability", "area_loco_predictions.rds"))
        list(metrics = mets, predictions = preds)
      }, list(res_list_val = at_list_expr))
    )
  ))
}


# ── Transportability targets (cross-country LOCO CV) ──────────────────────
# Train on N-1 countries, predict on the held-out country. Uses only proxy
# predictors shared across ALL countries (IHME, MAP, GEE).
# Restricted to women_iron only for speed (most informative outcome for
# cross-country comparison). To run all outcomes, change the line below.

# Only build transportability targets if we have 2+ countries
transport_targets <- list()
if (length(all_country_configs) >= 2) {

  # Outcomes measured in ALL 4 countries for cross-country LOCO
  first_cc <- all_country_configs[[1]]
  outcome_tags <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

  # Build symbol list for merged data targets (e.g., merged_gambia, merged_ghana, ...)
  country_names_lower <- tolower(names(all_country_configs))
  merged_syms <- lapply(paste0("merged_", country_names_lower), as.symbol)
  names(merged_syms) <- names(all_country_configs)

  for (otag in outcome_tags) {
    # Pooled dataset for this outcome.
    # We pass each merged_* target as a direct dependency by building a call
    # that constructs the list from the actual target values (not tar_read).
    pooled_name <- paste0("pooled_", otag)

    # Build the expression: list(Gambia = merged_gambia, Ghana = merged_ghana, ...)
    merged_list_expr <- as.call(c(list(as.symbol("list")), merged_syms))

    transport_targets <- c(transport_targets, list(
      tar_target_raw(
        pooled_name,
        substitute({
          all_merged <- merged_list_val
          build_pooled_dataset(all_merged, get_country_configs(), otag_val)
        }, list(
          merged_list_val = merged_list_expr,
          otag_val        = otag
        ))
      )
    ))

    # LOCO cross-validation for this outcome
    loco_name <- paste0("loco_", otag)
    transport_targets <- c(transport_targets, list(
      tar_target_raw(
        loco_name,
        substitute(
          run_loco_cv(pooled_data, sl_learners, pipeline_params),
          list(pooled_data = as.symbol(pooled_name))
        )
      )
    ))
  }

  # Combined transportability summary
  loco_target_names <- paste0("loco_", outcome_tags)
  loco_syms <- setNames(lapply(loco_target_names, as.symbol), outcome_tags)
  loco_list_expr <- as.call(c(list(as.symbol("list")), loco_syms))

  transport_targets <- c(transport_targets, list(
    tar_target_raw(
      "transportability_summary",
      substitute(
        summarize_transportability(loco_list_val),
        list(loco_list_val = loco_list_expr)
      )
    )
  ))

  # ── GEE-only LOCO: uses only remote sensing predictors ──────────────────
  for (otag in outcome_tags) {
    pooled_gee_name <- paste0("pooled_gee_", otag)

    transport_targets <- c(transport_targets, list(
      tar_target_raw(
        pooled_gee_name,
        substitute({
          all_merged <- merged_list_val
          build_pooled_gee_only(all_merged, get_country_configs(), otag_val)
        }, list(
          merged_list_val = merged_list_expr,
          otag_val        = otag
        ))
      )
    ))

    loco_gee_name <- paste0("loco_gee_", otag)
    transport_targets <- c(transport_targets, list(
      tar_target_raw(
        loco_gee_name,
        substitute(
          run_loco_gee_only(pooled_gee_data, sl_learners, pipeline_params),
          list(pooled_gee_data = as.symbol(pooled_gee_name))
        )
      )
    ))
  }

  # Combined GEE-only transportability summary
  loco_gee_target_names <- paste0("loco_gee_", outcome_tags)
  loco_gee_syms <- setNames(lapply(loco_gee_target_names, as.symbol), outcome_tags)
  loco_gee_list_expr <- as.call(c(list(as.symbol("list")), loco_gee_syms))

  transport_targets <- c(transport_targets, list(
    tar_target_raw(
      "transportability_gee_summary",
      substitute(
        summarize_transportability(loco_list_val),
        list(loco_list_val = loco_gee_list_expr)
      )
    )
  ))
}


# ── Out-of-sample prediction: Côte d'Ivoire ─────────────────────────────────
# Uses pooled area-level models from all surveyed countries to predict
# micronutrient deficiency prevalence in an unsurveyed country.
oos_targets <- list()

civ_raster_dir <- here::here("data", "Cote_dIvoire_GEE_rasters")
if (dir.exists(civ_raster_dir) && length(all_country_configs) >= 2) {

  # Extract GEE covariates for Côte d'Ivoire Admin-2 polygons (once)
  oos_targets <- c(oos_targets, list(
    tar_target(
      oos_gee_civ,
      extract_gee_for_country("CIV", civ_raster_dir)
    )
  ))

  # Build svy_admin2 symbol list for all countries.
  # OOS prediction requires svy_admin2 targets from ALL countries, so we
  # restrict to outcomes that exist in EVERY country config (intersection).
  all_outcome_sets <- lapply(all_country_configs, function(cc) names(cc$outcomes))
  outcome_tags_oos <- Reduce(intersect, all_outcome_sets)

  for (otag in outcome_tags_oos) {
    # Collect svy_admin2 targets for all countries into a named list
    svy_syms <- list()
    for (cn in names(all_country_configs)) {
      svy_name <- paste0("svy_admin2_", tolower(cn), "_", otag)
      svy_syms[[cn]] <- as.symbol(svy_name)
    }
    svy_list_expr <- as.call(c(list(as.symbol("list")), svy_syms))

    oos_name <- paste0("oos_civ_", otag)
    oos_targets <- c(oos_targets, list(
      tar_target_raw(
        oos_name,
        substitute({
          svy_list <- svy_list_val
          oc <- get_country_configs()[[1]]$outcomes[[otag_val]]
          predict_oos_pooled(
            svy_admin2_list = svy_list,
            country_configs = get_country_configs(),
            oos_gadm_code   = "CIV",
            oos_raster_dir  = raster_dir_val,
            oc              = oc,
            params          = pipeline_params,
            ext_cache_dir   = here::here("data", "external_cache"),
            oos_country_name = "Cote d'Ivoire"
          )
        }, list(
          svy_list_val   = svy_list_expr,
          otag_val       = otag,
          raster_dir_val = civ_raster_dir
        ))
      )
    ))
  }
}


# ── Summary targets (combine across outcomes) ───────────────────────────────
# These collect per-outcome results into combined tables

# Dynamically build the list of CV performance target names
cv_perf_names <- character()
admin2_error_names <- character()
ablation_names <- character()
diagnostics_names <- character()
diagnostics_calibrated_names <- character()
national_est_names <- character()
shap_names <- character()
varimp_names <- character()
area_model_names <- character()

for (country_name in names(all_country_configs)) {
  cc <- all_country_configs[[country_name]]
  for (outcome_name in names(cc$outcomes)) {
    suffix <- paste0(tolower(country_name), "_", outcome_name)
    cv_perf_names <- c(cv_perf_names, paste0("cv_perf_", suffix))
    admin2_error_names <- c(admin2_error_names, paste0("admin2_error_", suffix))
    ablation_names <- c(ablation_names, paste0("ablation_", suffix))
    diagnostics_names <- c(diagnostics_names, paste0("diagnostics_", suffix))
    diagnostics_calibrated_names <- c(diagnostics_calibrated_names,
                                       paste0("diagnostics_calibrated_", suffix))
    national_est_names <- c(national_est_names, paste0("national_est_", suffix))
    shap_names <- c(shap_names, paste0("shap_", suffix))
    varimp_names <- c(varimp_names, paste0("varimp_", suffix))
    area_model_names <- c(area_model_names, paste0("area_model_", suffix))
  }
}

# Helper: build a list(...) expression from target name strings
make_list_expr <- function(target_names) {
  syms <- lapply(target_names, as.symbol)
  as.call(c(list(as.symbol("list")), syms))
}

summary_targets <- list(

  # Combined CV performance table
  tar_target_raw(
    "cv_perf",
    substitute(
      dplyr::bind_rows(target_vals),
      list(target_vals = make_list_expr(cv_perf_names))
    )
  ),

  # Combined Admin2 error table
  tar_target_raw(
    "admin2_error_all",
    substitute(
      dplyr::bind_rows(target_vals),
      list(target_vals = make_list_expr(admin2_error_names))
    )
  ),

  # Combined domain ablation table
  tar_target_raw(
    "ablation_all",
    substitute(
      dplyr::bind_rows(target_vals),
      list(target_vals = make_list_expr(ablation_names))
    )
  ),

  # Combined diagnostics (metrics, calibration tables, PR curves)
  tar_target_raw(
    "diagnostics_all",
    substitute({
      diag_list <- diag_vals
      # Extract binary metrics
      bin_metrics <- dplyr::bind_rows(lapply(diag_list, function(x) x$binary_metrics))
      cont_metrics <- dplyr::bind_rows(lapply(diag_list, function(x) x$continuous_metrics))
      cal_tables <- dplyr::bind_rows(lapply(diag_list, function(x) x$calibration_table))
      pr_curves <- dplyr::bind_rows(lapply(diag_list, function(x) x$pr_curve))
      roc_curves <- dplyr::bind_rows(lapply(diag_list, function(x) x$roc_curve))
      list(
        binary_metrics = bin_metrics,
        continuous_metrics = cont_metrics,
        calibration_tables = cal_tables,
        pr_curves = pr_curves,
        roc_curves = roc_curves
      )
    },
    list(diag_vals = make_list_expr(diagnostics_names)))
  ),

  # Combined calibrated diagnostics — same shape as diagnostics_all but the
  # binary OOF probabilities have been re-passed through a Platt logistic
  # recalibrator before metrics are recomputed. Writes a CSV side-by-side
  # with the raw diagnostics so the user can see Brier-skill recovery.
  tar_target_raw(
    "diagnostics_calibrated_all",
    substitute({
      diag_list <- diag_vals
      bin_metrics <- dplyr::bind_rows(lapply(diag_list, function(x) x$binary_metrics))
      cal_tables <- dplyr::bind_rows(lapply(diag_list, function(x) x$calibration_table))
      out_dir <- here::here("results", "tables")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      if (!is.null(bin_metrics) && nrow(bin_metrics) > 0)
        readr::write_csv(bin_metrics,
                          file.path(out_dir, "diagnostics_binary_calibrated.csv"))
      list(binary_metrics = bin_metrics, calibration_tables = cal_tables)
    },
    list(diag_vals = make_list_expr(diagnostics_calibrated_names)))
  ),

  # Combined national estimates
  tar_target_raw(
    "national_estimates_all",
    substitute(
      dplyr::bind_rows(target_vals),
      list(target_vals = make_list_expr(national_est_names))
    )
  ),

  # Combined SHAP outputs (district top factors + global importance)
  tar_target_raw(
    "shap_all",
    substitute({
      shap_list <- target_vals
      district_factors <- dplyr::bind_rows(lapply(shap_list,
                                                    function(x) x$district_factors))
      global_importance <- dplyr::bind_rows(lapply(shap_list,
                                                     function(x) x$global_importance))
      list(district_factors = district_factors,
           global_importance = global_importance)
    }, list(target_vals = make_list_expr(shap_names)))
  ),

  # Combined per-variable importance
  tar_target_raw(
    "varimp_all",
    substitute(
      dplyr::bind_rows(target_vals),
      list(target_vals = make_list_expr(varimp_names))
    )
  ),

  # Save combined CSV tables
  tar_target(
    save_tables,
    {
      out_dir <- here::here("results", "tables")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      if (!is.null(cv_perf) && nrow(cv_perf) > 0)
        readr::write_csv(cv_perf, file.path(out_dir, "cv_performance_all.csv"))

      if (!is.null(admin2_error_all) && nrow(admin2_error_all) > 0)
        readr::write_csv(admin2_error_all, file.path(out_dir, "admin2_error_all.csv"))

      if (!is.null(ablation_all) && nrow(ablation_all) > 0)
        readr::write_csv(ablation_all, file.path(out_dir, "domain_ablation_all.csv"))

      if (!is.null(diagnostics_all$binary_metrics) && nrow(diagnostics_all$binary_metrics) > 0)
        readr::write_csv(diagnostics_all$binary_metrics,
                         file.path(out_dir, "diagnostics_binary.csv"))
      if (!is.null(diagnostics_all$continuous_metrics) && nrow(diagnostics_all$continuous_metrics) > 0)
        readr::write_csv(diagnostics_all$continuous_metrics,
                         file.path(out_dir, "diagnostics_continuous.csv"))
      # Calibration tables and PR/ROC curve data (for dashboard reuse)
      if (!is.null(diagnostics_all$calibration_tables) &&
          nrow(diagnostics_all$calibration_tables) > 0)
        readr::write_csv(diagnostics_all$calibration_tables,
                         file.path(out_dir, "calibration_tables.csv"))
      if (!is.null(diagnostics_all$pr_curves) &&
          nrow(diagnostics_all$pr_curves) > 0)
        readr::write_csv(diagnostics_all$pr_curves,
                         file.path(out_dir, "pr_curves.csv"))
      if (!is.null(diagnostics_all$roc_curves) &&
          nrow(diagnostics_all$roc_curves) > 0)
        readr::write_csv(diagnostics_all$roc_curves,
                         file.path(out_dir, "roc_curves.csv"))

      # Save national estimates
      if (!is.null(national_estimates_all) && nrow(national_estimates_all) > 0)
        readr::write_csv(national_estimates_all,
                         file.path(out_dir, "national_estimates_all.csv"))

      # SHAP and per-variable importance — for dashboard
      if (!is.null(shap_all$district_factors) &&
          nrow(shap_all$district_factors) > 0)
        readr::write_csv(shap_all$district_factors,
                         file.path(out_dir, "shap_district_factors.csv"))
      if (!is.null(shap_all$global_importance) &&
          nrow(shap_all$global_importance) > 0)
        readr::write_csv(shap_all$global_importance,
                         file.path(out_dir, "shap_global_importance.csv"))
      if (!is.null(varimp_all) && nrow(varimp_all) > 0)
        readr::write_csv(varimp_all,
                         file.path(out_dir, "single_var_importance.csv"))

      list(cv_perf = cv_perf, admin2_error = admin2_error_all,
           ablation = ablation_all, diagnostics = diagnostics_all,
           national = national_estimates_all,
           shap = shap_all, varimp = varimp_all)
    }
  )
)

# =============================================================================
# CORRECTED METHODS (P1–P8)  —  folded in from the former _targets_corrected.R
# =============================================================================
# The honest "corrected-methods" layer, previously a separate parallel pipeline
# (_targets_corrected.R, now archived), is a section of this single pipeline.
# It consumes the SAME per-slice data targets built above (outcome_data_*,
# svy_admin2_*, gee_admin2_*) and applies:
#   P1 leakage-free in-fold preprocessing + cluster/spatial-block CV
#   P2 out-of-fold Platt calibration
#   P3 decision-value targeting
#   P4 split-conformal + design-based intervals
#   P5 sampling-error-aware admin-2 error
#   P6 out-of-support trust flags
#   P7 partial-pooling small-area estimator
#   P8 join/centering bug fixes (enforced by construction in R/corrected/)
# It emits the canonical `corrected_methods_comparison` bundle, which
# build_methods_comparison() writes to results/tables/corrected/*.csv and
# dashboard/data/methods_comparison.rds. Functions live in R/corrected/.
# Config objects (cc, oc) are embedded as literals exactly as make_outcome_targets
# does, which also avoids the get_country_configs()[[label]] keying pitfall
# (the registry is keyed "SierraLeone", not the display label "Sierra Leone").
CORRECTED_LIB  <- c("mean", "glmnet", "ranger", "rpart")
CORRECTED_MAXP <- 800L
CORRECTED_V    <- 5L

corrected_targets     <- list()
corrected_result_syms <- list()
corrected_recon_syms  <- list()   # protocol-reconciliation (P9) per-slice rows
area_recipe_syms      <- list()   # area-level recipe per-slice within-country rows
ar_frame_univ_syms    <- list()   # universal-mode area frames, grouped by outcome (transport)
for (country_name in names(all_country_configs)) {
  cc  <- all_country_configs[[country_name]]
  low <- tolower(country_name)
  for (outcome_name in names(cc$outcomes)) {
    oc      <- cc$outcomes[[outcome_name]]
    suffix  <- paste0(low, "_", outcome_name)
    od_sym  <- as.symbol(paste0("outcome_data_", suffix))
    svy_sym <- as.symbol(paste0("svy_admin2_",  suffix))
    gee_sym <- as.symbol(paste0("gee_admin2_",  low))
    sl_nm   <- paste0("corrected_sl_",  suffix)
    int_nm  <- paste0("corrected_int_", suffix)
    res_nm  <- paste0("corrected_result_", suffix)
    corrected_result_syms[[suffix]] <- as.symbol(res_nm)
    corrected_recon_syms[[suffix]]  <- as.symbol(paste0("corrected_recon_", suffix))
    area_recipe_syms[[suffix]]      <- as.symbol(paste0("area_recipe_", suffix))
    ar_frame_univ_syms[[outcome_name]] <- c(
      ar_frame_univ_syms[[outcome_name]],
      setNames(list(as.symbol(paste0("ar_frame_univ_", suffix))), country_name))

    corrected_targets <- c(corrected_targets, list(
      # P1: leakage-free SL fit (honest cluster + spatial-block + optimistic CV)
      tar_target_raw(sl_nm,
        substitute(fit_corrected_sl(OD, CCV, OCV, LAB, library = LIBR,
                                    max_pred = MP, V = VV),
          list(OD = od_sym, CCV = cc, OCV = oc, LAB = cc$country,
               LIBR = CORRECTED_LIB, MP = CORRECTED_MAXP, VV = CORRECTED_V))),
      tar_target_raw(paste0("corrected_cvperf_", suffix),
        substitute(extract_cv_perf_corrected(SL), list(SL = as.symbol(sl_nm)))),
      tar_target_raw(paste0("corrected_prodcv_", suffix),
        substitute(prod_cv_for(LAB, OCN),
                   list(LAB = cc$country, OCN = outcome_name))),
      # P2 out-of-fold calibration
      tar_target_raw(paste0("corrected_calib_", suffix),
        substitute(diagnostics_calibrated_oof(SL), list(SL = as.symbol(sl_nm)))),
      # P5 sampling-error-aware admin-2 error
      tar_target_raw(paste0("corrected_err_", suffix),
        substitute(admin2_error_corrected(SL, SV),
                   list(SL = as.symbol(sl_nm), SV = svy_sym))),
      # P3 decision value
      tar_target_raw(paste0("corrected_dec_", suffix),
        substitute(decision_value_corrected(SL, SV),
                   list(SL = as.symbol(sl_nm), SV = svy_sym))),
      # P4 split-conformal + design-based intervals
      tar_target_raw(int_nm,
        substitute(intervals_corrected(SL, SV),
                   list(SL = as.symbol(sl_nm), SV = svy_sym))),
      # P6 out-of-support trust flags
      tar_target_raw(paste0("corrected_trust_", suffix),
        substitute(trust_flags_corrected(SL, SV, GEE, INT),
                   list(SL = as.symbol(sl_nm), SV = svy_sym, GEE = gee_sym,
                        INT = as.symbol(int_nm)))),
      # P7 partial-pooling area estimator
      tar_target_raw(paste0("corrected_area_", suffix),
        substitute(area_partial_pooling_corrected(SL, SV, GEE),
                   list(SL = as.symbol(sl_nm), SV = svy_sym, GEE = gee_sym))),
      # P9 individual-vs-area district-ranking on ONE matched protocol (reuses SL fit)
      tar_target_raw(paste0("corrected_recon_", suffix),
        substitute(reconcile_protocols(SL, OD, SV, GEE, CCV, OCV),
                   list(SL = as.symbol(sl_nm), OD = od_sym, SV = svy_sym,
                        GEE = gee_sym, CCV = cc, OCV = oc))),
      # Area-level recipe (docs/AREA_LEVEL_RECIPE_SPEC.md): within-country
      # district+region honest CV, enriched features, CI/perm-null metrics.
      tar_target_raw(paste0("area_recipe_", suffix),
        substitute(ar_within_country(OD, SV, GEE, CCV, OCV, mode = "enriched"),
                   list(OD = od_sym, SV = svy_sym, GEE = gee_sym, CCV = cc, OCV = oc))),
      # universal-mode area frame (gee only) for the cross-country transport rollups
      tar_target_raw(paste0("ar_frame_univ_", suffix),
        substitute(ar_build_frame(OD, SV, GEE, CCV, OCV, mode = "universal"),
                   list(OD = od_sym, SV = svy_sym, GEE = gee_sym, CCV = cc, OCV = oc))),
      # per-slice result bundle
      tar_target_raw(res_nm,
        substitute(list(cv_perf = CV, prod_cv = PCV, calibration = CAL,
                        admin2_error = ERR, decision = DEC, trust = TR,
                        area_pp = AR, interval_summary = INT$summary),
          list(CV  = as.symbol(paste0("corrected_cvperf_", suffix)),
               PCV = as.symbol(paste0("corrected_prodcv_", suffix)),
               CAL = as.symbol(paste0("corrected_calib_", suffix)),
               ERR = as.symbol(paste0("corrected_err_", suffix)),
               DEC = as.symbol(paste0("corrected_dec_", suffix)),
               TR  = as.symbol(paste0("corrected_trust_", suffix)),
               AR  = as.symbol(paste0("corrected_area_", suffix)),
               INT = as.symbol(int_nm))))
    ))
  }
}
# Final canonical bundle (corrected metrics + production reference, side by side)
corrected_result_list_expr <- as.call(c(list(as.symbol("list")),
  setNames(corrected_result_syms, names(corrected_result_syms))))
corrected_targets <- c(corrected_targets, list(
  tar_target_raw("corrected_methods_comparison",
    substitute(build_methods_comparison(RL),
               list(RL = corrected_result_list_expr)))
))

# P9 roll-up: matched-protocol reconciliation of individual-vs-area district
# ranking (writes results/tables/corrected/protocol_reconciliation*.csv).
corrected_recon_list_expr <- as.call(c(list(as.symbol("list")),
  setNames(corrected_recon_syms, names(corrected_recon_syms))))
corrected_targets <- c(corrected_targets, list(
  tar_target_raw("protocol_reconciliation",
    substitute(build_protocol_reconciliation(RL),
               list(RL = corrected_recon_list_expr)))
))

# Area-level recipe roll-ups (docs/AREA_LEVEL_RECIPE_SPEC.md). Cross-country LOCO
# transport per outcome (universal features) for outcomes in >=3 countries, then
# a table writer for results/tables/corrected/area_recipe_{within,transport}.csv.
ar_transport_syms <- list()
for (oc_nm in names(ar_frame_univ_syms)) {
  fr <- ar_frame_univ_syms[[oc_nm]]
  if (length(fr) >= 3) {
    frame_list_expr <- as.call(c(list(as.symbol("list")), fr))   # named by country
    corrected_targets <- c(corrected_targets, list(
      tar_target_raw(paste0("area_recipe_transport_", oc_nm),
        substitute(ar_transport_loco(FR), list(FR = frame_list_expr)))))
    ar_transport_syms[[oc_nm]] <- as.symbol(paste0("area_recipe_transport_", oc_nm))
  }
}
area_recipe_within_expr <- as.call(c(list(as.symbol("list")),
  setNames(area_recipe_syms, names(area_recipe_syms))))
area_recipe_transport_expr <- as.call(c(list(as.symbol("list")),
  setNames(ar_transport_syms, names(ar_transport_syms))))
corrected_targets <- c(corrected_targets, list(
  tar_target_raw("area_recipe_tables",
    substitute(build_area_recipe_tables(W, TR),
               list(W = area_recipe_within_expr, TR = area_recipe_transport_expr)))))

# ── Combine everything ──────────────────────────────────────────────────────
# ── Aggregation-level national-prevalence sweep ─────────────────────────────
# Refreshes with the pipeline: for each country x outcome, compares national-
# prevalence accuracy when the model is fit at cluster / admin-1 / admin-2 level
# (see R/aggregation_level_sweep.R). Depends on every outcome_data_* target, so
# it re-runs whenever those change. Writes the summary CSV as a side effect.
aggregation_targets <- list()
if (length(all_country_configs) >= 1) {
  od_syms <- list()
  for (ck in names(all_country_configs)) {
    # Skip countries whose merged dataset isn't built yet (e.g. Tanzania WIP) so
    # this target doesn't depend on outcome_data_* targets that can't build.
    if (!file.exists(all_country_configs[[ck]]$data_path)) next
    for (tag in names(all_country_configs[[ck]]$outcomes)) {
      sfx <- paste0(tolower(ck), "_", tag)
      od_syms[[sfx]] <- as.symbol(paste0("outcome_data_", sfx))
    }
  }
  od_list_expr <- as.call(c(list(as.symbol("list")), od_syms))
  aggregation_targets <- list(
    tar_target_raw(
      "aggregation_level_national",
      substitute({
        res <- run_aggregation_level_sweep(od_list_val, get_country_configs())
        dir.create(here::here("results", "sensitivity"),
                   showWarnings = FALSE, recursive = TRUE)
        if (nrow(res)) readr::write_csv(res, here::here(
          "results", "sensitivity", "aggregation_level_national_summary.csv"))
        res
      }, list(od_list_val = od_list_expr))
    )
  )
}

c(static_targets, country_targets, area_comparison_targets, area_loco_targets,
  benchmark_targets, area_transport_targets, transport_targets, oos_targets,
  summary_targets, corrected_targets, aggregation_targets)
