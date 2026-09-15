# =============================================================================
# National Prediction Pipeline — targets definition
# =============================================================================
# Standalone pipeline for BRINDA and VMNIS national-level prevalence prediction.
# Run with: targets::tar_make(script = "national_prediction/_targets.R",
#                              store  = "national_prediction/_targets")
# =============================================================================

library(targets)
library(tarchetypes)

# Source all R function files
tar_source("national_prediction/R/")

# Package dependencies
tar_option_set(
  packages = c(
    "here", "dplyr", "tidyr", "readr", "tibble", "rlang", "stringr", "purrr",
    "ggplot2", "patchwork", "scales",
    "haven", "labelled", "data.table",
    "sl3", "origami", "caret", "recipes", "ck37r", "washb",
    "glmnet", "ranger", "xgboost"
  ),
  memory = "transient",
  garbage_collection = TRUE
)


# =============================================================================
# Pipeline targets
# =============================================================================
list(

  # =========================================================================
  # Stage 1: Configuration
  # =========================================================================
  tar_target(data_paths, get_data_paths()),
  tar_target(sl_params, get_sl_params()),
  tar_target(brinda_outcomes, get_brinda_outcomes()),
  tar_target(vmnis_outcomes, get_vmnis_outcomes()),
  tar_target(lmic_countries, get_lmic_countries()),

  # =========================================================================
  # Stage 2: Data loading & cleaning
  # =========================================================================
  tar_target(
    predictor_data,
    load_predictor_data(data_paths$predictor_data)
  ),
  tar_target(
    brinda_clean,
    clean_brinda_data(data_paths$brinda_data, predictor_data)
  ),
  tar_target(
    vmnis_clean,
    clean_vmnis_data(data_paths$vmnis_data, predictor_data)
  ),
  tar_target(
    vmnis_heatmap_data,
    clean_vmnis_for_heatmap(data_paths$vmnis_data, lmic_countries)
  ),

  # Predictor variable lists
  tar_target(brinda_xvars, get_brinda_xvars(brinda_clean)),
  tar_target(vmnis_xvars, get_vmnis_xvars(vmnis_clean)),

  # SuperLearner stacks (built once, reused across outcomes)
  tar_target(sl_full, build_sl_full()),
  tar_target(sl_simple, build_sl_simple()),

  # =========================================================================
  # Stage 3: BRINDA model fitting
  # =========================================================================
  tar_target(
    brinda_fit_prev_ferr,
    run_national_sl(brinda_clean, "prev_ferr", brinda_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    brinda_fit_prev_stfr,
    run_national_sl(brinda_clean, "prev_stfr", brinda_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    brinda_fit_prev_rbp,
    run_national_sl(brinda_clean, "prev_rbp", brinda_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    brinda_fit_prev_zinc,
    run_national_sl(brinda_clean, "prev_zinc", brinda_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),

  # =========================================================================
  # Stage 4: VMNIS model fitting
  # =========================================================================

  # --- Zinc ---
  tar_target(
    vmnis_d_zinc,
    filter_vmnis_indicator(vmnis_clean, "Zinc (plasma or serum)")
  ),
  tar_target(
    vmnis_fit_mean_zinc,
    run_national_sl(vmnis_d_zinc, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    vmnis_fit_prev_zinc,
    run_national_sl(vmnis_d_zinc, "Prevalenceofdeficiency", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),

  # --- Vitamin B12 (sl_simple) ---
  tar_target(
    vmnis_d_b12,
    filter_vmnis_indicator(vmnis_clean, "Vitamin B12")
  ),
  tar_target(
    vmnis_fit_mean_b12,
    run_national_sl(vmnis_d_b12, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_simple)
  ),
  tar_target(
    vmnis_fit_prev_b12,
    run_national_sl(vmnis_d_b12, "Prevalenceofdeficiency", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_simple)
  ),

  # --- Retinol ---
  tar_target(
    vmnis_d_retinol,
    filter_vmnis_indicator(vmnis_clean, "Retinol (plasma or serum)")
  ),
  tar_target(
    vmnis_fit_mean_retinol,
    run_national_sl(vmnis_d_retinol, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    vmnis_fit_prev_retinol,
    run_national_sl(vmnis_d_retinol, "Prevalenceofdeficiency", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),

  # --- Vitamin D (sl_simple) ---
  tar_target(
    vmnis_d_vitD,
    filter_vmnis_indicator(vmnis_clean, "25-Hydroxyvitamin D")
  ),
  tar_target(
    vmnis_fit_mean_vitD,
    run_national_sl(vmnis_d_vitD, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_simple)
  ),
  tar_target(
    vmnis_fit_prev_vitD,
    run_national_sl(vmnis_d_vitD, "Prevalenceofdeficiency", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_simple)
  ),

  # --- Folate ---
  tar_target(
    vmnis_d_folate,
    filter_vmnis_indicator(vmnis_clean, "Folate (plasma or serum)")
  ),
  tar_target(
    vmnis_fit_mean_folate,
    run_national_sl(vmnis_d_folate, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),
  tar_target(
    vmnis_fit_prev_folate,
    run_national_sl(vmnis_d_folate, "Prevalenceofdeficiency", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_full)
  ),

  # --- Haemoglobin (mean only; sl_simple due to N=27k rows) ---
  tar_target(
    vmnis_d_haemoglobin,
    filter_vmnis_indicator(vmnis_clean, "Haemoglobin")
  ),
  tar_target(
    vmnis_fit_mean_haemoglobin,
    run_national_sl(vmnis_d_haemoglobin, "Mean", vmnis_xvars,
                    folds = sl_params$folds, CV = sl_params$CV, sl = sl_simple)
  ),

  # =========================================================================
  # Stage 5: Performance evaluation
  # =========================================================================
  tar_target(
    brinda_perf,
    dplyr::bind_rows(
      safe_evaluate_performance(brinda_fit_prev_ferr, "prev_ferr"),
      safe_evaluate_performance(brinda_fit_prev_stfr, "prev_stfr"),
      safe_evaluate_performance(brinda_fit_prev_rbp, "prev_rbp"),
      safe_evaluate_performance(brinda_fit_prev_zinc, "prev_zinc")
    ) %>% dplyr::mutate(source = "BRINDA")
  ),

  tar_target(
    vmnis_perf,
    dplyr::bind_rows(
      safe_evaluate_performance(vmnis_fit_mean_zinc, "mean_zinc"),
      safe_evaluate_performance(vmnis_fit_prev_zinc, "prev_zinc"),
      safe_evaluate_performance(vmnis_fit_mean_b12, "mean_b12"),
      safe_evaluate_performance(vmnis_fit_prev_b12, "prev_b12"),
      safe_evaluate_performance(vmnis_fit_mean_retinol, "mean_retinol"),
      safe_evaluate_performance(vmnis_fit_prev_retinol, "prev_retinol"),
      safe_evaluate_performance(vmnis_fit_mean_vitD, "mean_vitD"),
      safe_evaluate_performance(vmnis_fit_prev_vitD, "prev_vitD"),
      safe_evaluate_performance(vmnis_fit_mean_folate, "mean_folate"),
      safe_evaluate_performance(vmnis_fit_prev_folate, "prev_folate"),
      safe_evaluate_performance(vmnis_fit_mean_haemoglobin, "mean_haemoglobin")
    ) %>% dplyr::mutate(source = "VMNIS")
  ),

  tar_target(
    all_perf,
    dplyr::bind_rows(brinda_perf, vmnis_perf)
  ),

  # =========================================================================
  # Stage 6: Variable importance
  # =========================================================================

  # BRINDA variable importance (one per outcome)
  tar_target(varimp_brinda_ferr, safe_calc_importance(brinda_fit_prev_ferr, "prev_ferr")),
  tar_target(varimp_brinda_stfr, safe_calc_importance(brinda_fit_prev_stfr, "prev_stfr")),
  tar_target(varimp_brinda_rbp,  safe_calc_importance(brinda_fit_prev_rbp, "prev_rbp")),
  tar_target(varimp_brinda_zinc, safe_calc_importance(brinda_fit_prev_zinc, "prev_zinc")),

  # VMNIS variable importance (selected outcomes)
  tar_target(varimp_vmnis_zinc,    safe_calc_importance(vmnis_fit_mean_zinc, "mean_zinc")),
  tar_target(varimp_vmnis_retinol, safe_calc_importance(vmnis_fit_mean_retinol, "mean_retinol")),
  tar_target(varimp_vmnis_folate,  safe_calc_importance(vmnis_fit_mean_folate, "mean_folate")),
  tar_target(varimp_vmnis_haemoglobin, safe_calc_importance(vmnis_fit_mean_haemoglobin, "mean_haemoglobin")),

  # =========================================================================
  # Stage 7: Visualization
  # =========================================================================
  tar_target(
    fig_heatmap_pw,
    plot_vmnis_heatmap(vmnis_heatmap_data, "Pregnant women")
  ),
  tar_target(
    fig_heatmap_pw_alt,
    plot_vmnis_heatmap(vmnis_heatmap_data, "Pregnant women",
                       exclude_haemoglobin = TRUE)
  ),
  tar_target(
    fig_heatmap_all,
    plot_vmnis_heatmap_all(vmnis_heatmap_data)
  ),
  tar_target(
    fig_performance,
    plot_performance_comparison(all_perf)
  ),
  tar_target(
    fig_obs_vs_pred_brinda,
    plot_obs_vs_pred(
      list(prev_ferr = brinda_fit_prev_ferr,
           prev_stfr = brinda_fit_prev_stfr,
           prev_rbp  = brinda_fit_prev_rbp,
           prev_zinc = brinda_fit_prev_zinc),
      outcome_labels = c(prev_ferr = "Iron (Ferritin)",
                         prev_stfr = "Iron (sTfR)",
                         prev_rbp  = "Vitamin A (RBP)",
                         prev_zinc = "Zinc")
    )
  ),
  tar_target(
    fig_obs_vs_pred_vmnis,
    plot_obs_vs_pred(
      list(mean_zinc     = vmnis_fit_mean_zinc,
           mean_retinol  = vmnis_fit_mean_retinol,
           mean_folate   = vmnis_fit_mean_folate,
           mean_b12      = vmnis_fit_mean_b12,
           mean_vitD     = vmnis_fit_mean_vitD,
           mean_haemoglobin = vmnis_fit_mean_haemoglobin)
    )
  ),

  # Variable importance plots (BRINDA)
  tar_target(
    fig_varimp_brinda_ferr,
    if (!is.null(varimp_brinda_ferr))
      plot_varimp_bar(varimp_brinda_ferr$varimp_top, "Iron (Ferritin)")
  ),
  tar_target(
    fig_varimp_brinda_zinc,
    if (!is.null(varimp_brinda_zinc))
      plot_varimp_bar(varimp_brinda_zinc$varimp_top, "Zinc (BRINDA)")
  ),

  # Variable importance plots (VMNIS)
  tar_target(
    fig_varimp_vmnis_zinc,
    if (!is.null(varimp_vmnis_zinc))
      plot_varimp_bar(varimp_vmnis_zinc$varimp_top, "Zinc (VMNIS)")
  ),
  tar_target(
    fig_varimp_vmnis_retinol,
    if (!is.null(varimp_vmnis_retinol))
      plot_varimp_bar(varimp_vmnis_retinol$varimp_top, "Retinol (VMNIS)")
  ),

  # =========================================================================
  # Stage 8: Save outputs
  # =========================================================================
  tar_target(
    save_tables,
    {
      readr::write_csv(all_perf,
        here::here("national_prediction", "results", "tables",
                   "national_sl_performance.csv"))
      readr::write_csv(brinda_perf,
        here::here("national_prediction", "results", "tables",
                   "brinda_performance.csv"))
      readr::write_csv(vmnis_perf,
        here::here("national_prediction", "results", "tables",
                   "vmnis_performance.csv"))
      "tables_saved"
    }
  ),

  tar_target(
    save_figures,
    {
      fig_dir <- here::here("national_prediction", "results", "figures")

      ggplot2::ggsave(file.path(fig_dir, "vmnis_heatmap_pw.png"),
                      fig_heatmap_pw, width = 10, height = 8)
      ggplot2::ggsave(file.path(fig_dir, "vmnis_heatmap_pw_alt.png"),
                      fig_heatmap_pw_alt, width = 10, height = 8)
      ggplot2::ggsave(file.path(fig_dir, "vmnis_heatmap_all.png"),
                      fig_heatmap_all, width = 14, height = 10)
      ggplot2::ggsave(file.path(fig_dir, "performance_comparison.png"),
                      fig_performance, width = 12, height = 8)
      ggplot2::ggsave(file.path(fig_dir, "obs_vs_pred_brinda.png"),
                      fig_obs_vs_pred_brinda, width = 10, height = 8)
      ggplot2::ggsave(file.path(fig_dir, "obs_vs_pred_vmnis.png"),
                      fig_obs_vs_pred_vmnis, width = 12, height = 8)

      # Variable importance (save if non-null)
      if (!is.null(fig_varimp_brinda_ferr))
        ggplot2::ggsave(file.path(fig_dir, "varimp_brinda_ferr.png"),
                        fig_varimp_brinda_ferr, width = 8, height = 6)
      if (!is.null(fig_varimp_brinda_zinc))
        ggplot2::ggsave(file.path(fig_dir, "varimp_brinda_zinc.png"),
                        fig_varimp_brinda_zinc, width = 8, height = 6)
      if (!is.null(fig_varimp_vmnis_zinc))
        ggplot2::ggsave(file.path(fig_dir, "varimp_vmnis_zinc.png"),
                        fig_varimp_vmnis_zinc, width = 8, height = 6)
      if (!is.null(fig_varimp_vmnis_retinol))
        ggplot2::ggsave(file.path(fig_dir, "varimp_vmnis_retinol.png"),
                        fig_varimp_vmnis_retinol, width = 8, height = 6)

      "figures_saved"
    }
  ),

  # =========================================================================
  # Stage 9: Results report
  # =========================================================================
  tar_target(
    report,
    {
      # Ensure all outputs are available before rendering
      force(save_tables)
      force(save_figures)
      quarto::quarto_render(
        here::here("national_prediction", "docs", "results_report.qmd"),
        output_file = "results_report.html"
      )
      here::here("national_prediction", "docs", "results_report.html")
    }
  )
)
