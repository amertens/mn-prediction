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
    "recipes", "future.apply", "glmnet", "pROC", "ROCR", "haven", "readxl"
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
# DYNAMIC TARGET FACTORY
#
# For each country x outcome, generates the full set of targets.
# This makes it trivial to add countries/outcomes without editing this file.
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

  # Load raw merged data (one per country — cached, only reloads if file changes)
  raw_target_name <- paste0("merged_raw_", tolower(country_name))
  raw_data_target <- tar_target_raw(
    raw_target_name,
    substitute(
      load_merged_data(data_path_val),
      list(data_path_val = cc$data_path)
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

  # Extract GEE raster zonal means for Admin-2 polygons (once per country)
  gee_admin2_target_name <- paste0("gee_admin2_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      gee_admin2_target_name,
      substitute(
        extract_gee_admin2(cc_val),
        list(cc_val = cc)
      )
    )
  ))

  # Merge external predictors (CHIRPS, WorldPop, MAP, HarvestStat, etc.)
  # AND GEE Admin-2 zonal means into individual-level data.
  merged_ext_target_name <- paste0("merged_ext_", tolower(country_name))
  country_targets <- c(country_targets, list(
    tar_target_raw(
      merged_ext_target_name,
      substitute(
        merge_external_predictors(fsec_data, cc_val, cache_dir_val),
        list(
          fsec_data     = as.symbol(fsec_target_name),
          cc_val        = cc,
          cache_dir_val = here::here("data", "external_cache")
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
  for (oc_name_area in names(cc_area_local$outcomes)) {
    oc_area_local <- cc_area_local$outcomes[[oc_name_area]]
    area_suffix <- paste0(tolower(cc_name_area), "_", oc_name_area)

    country_targets <- c(country_targets, list(
      tar_target_raw(
        paste0("area_model_", area_suffix),
        substitute(
          fit_area_level_model(svy_admin2, cc_val, oc_val, params_val),
          list(
            svy_admin2 = as.symbol(paste0("svy_admin2_", area_suffix)),
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

# ── Combine everything ──────────────────────────────────────────────────────
c(static_targets, country_targets, area_comparison_targets, area_loco_targets,
  transport_targets, oos_targets, summary_targets)
