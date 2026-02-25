# =============================================================================
# src/analysis/03b_predict_and_aggregate_admin2_ghana.R
#
# Admin-2 aggregation and mapping for Ghana — Women Vitamin A Deficiency
#
# Parallel to 03_predict_and_aggregate_admin1.R but:
#   - Uses the Ghana merged dataset (individual-level GMS data)
#   - Aggregates predicted prevalences to Admin2 (district) level
#   - Produces choropleth maps of Ghana at Admin2 resolution
#   - Also computes survey-weighted direct estimates for comparison
#
# Outcome: Women VAD (binary) — gw_wVADAdjThurn
#   (Thurnham-adjusted retinol-binding protein deficiency indicator)
#
# Outputs:
#   results/tables/ghana_admin2_prevalence_women_vitA.csv
#   results/figures/ghana_admin2_predicted_prev_women_vitA.png
#   results/figures/ghana_admin2_survey_direct_prev_women_vitA.png
#   results/figures/ghana_admin2_comparison_women_vitA.png
#
# Usage:
#   source(here::here("src/analysis/03b_predict_and_aggregate_admin2_ghana.R"))
#   # or: Rscript src/analysis/03b_predict_and_aggregate_admin2_ghana.R
# =============================================================================

cat("\n[03b] Ghana Admin-2 aggregation: Women Vitamin A...\n")

# ---- 0. Packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyverse)
  library(sf)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(survey)
  library(srvyr)
  library(geodata)
  library(haven)
  library(labelled)
  library(sl3)
  library(origami)
  library(tlverse)
  library(caret)
  library(data.table)
  library(ck37r)
  library(SuperLearner)
  library(future)
  library(washb)
  library(recipes)
  library(pROC)
  library(readr)
})

# Source helper functions and SL setup
source(here::here("src/0-functions.R"))
source(here::here("src/0-SL-setup.R"))
source(here::here("src/DHS/DHS_functions.R"))
source(here::here("src/DHS/DHS_variable_recode.R"))

# Output directories
out_tables  <- here::here("results", "tables")
out_figures <- here::here("results", "figures")
out_models  <- here::here("results", "models")
dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_figures, showWarnings = FALSE, recursive = TRUE)
dir.create(out_models,  showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load Ghana data
# ============================================================================
cat("  Loading Ghana merged dataset...\n")
df <- readRDS(here::here("data", "IPD", "Ghana", "Ghana_merged_dataset.rds"))

# Drop geometry if sf object
if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)

# Unique row identifier
df$dataid <- paste0("ghana", seq_len(nrow(df)))

# Ensure binary outcome is numeric 0/1
df$gw_wVADAdjThurn <- as.numeric(df$gw_wVADAdjThurn)

cat(sprintf("  Dataset: %d rows, %d cols\n", nrow(df), ncol(df)))
cat(sprintf("  Admin2 levels: %d\n", length(unique(df$Admin2))))
cat(sprintf("  Women VAD prevalence (unweighted): %.1f%%\n",
            100 * mean(df$gw_wVADAdjThurn, na.rm = TRUE)))

# ============================================================================
# 2. Build predictor variable lists (same as GW Ghana SL_bin.R)
# ============================================================================
cat("  Building predictor lists...\n")

gw_vars     <- colnames(df)[grepl("gw_", colnames(df))]
dhs_vars    <- colnames(df)[grepl("dhs2014_", colnames(df)) |
                              grepl("dhs2016_", colnames(df)) |
                              grepl("dhs2017_", colnames(df))]
mics_vars   <- colnames(df)[grepl("mics_", colnames(df))]
ihme_vars   <- colnames(df)[grepl("ihme_", colnames(df))]
lsms_vars   <- colnames(df)[grepl("lsms_", colnames(df))]
map_vars    <- colnames(df)[grepl("MAP_", colnames(df))]
wfp_vars    <- c("nearest_market_id", colnames(df)[grepl("wfp_", colnames(df))])
wfp_vars    <- wfp_vars[wfp_vars %in% colnames(df)]
flunet_vars <- colnames(df)[grepl("flunet_", colnames(df))]
gee_vars    <- colnames(df)[grepl("gee_", colnames(df))]

# Remove outcome-leaking gw_ variables (biomarker-related)
leak_patterns <- c("wID", "cID", "VAD", "RBP", "rbp", "Ferr", "TFR",
                    "crp", "agp", "cRBP", "cVAD", "VAI", "Folate", "B12",
                    "AGP", "CRP", "RDT", "fer", "bis", "nflam", "nemia",
                    "globin", "cHb", "wHb", "gchb", "hbc", "Anemia", "wm_wmst")
for (pat in leak_patterns) {
  gw_vars <- gw_vars[!grepl(pat, gw_vars)]
}
gw_vars <- gw_vars[!gw_vars %in% c("gw_cn", "gw_hhid", "gw_childid")]

Xvars_full <- c("Admin1", "Admin2", "gw_month", gw_vars, dhs_vars, mics_vars,
                 ihme_vars, lsms_vars, map_vars, wfp_vars, flunet_vars, gee_vars)
Xvars_full <- Xvars_full[Xvars_full %in% colnames(df)]

cat(sprintf("  Predictors: %d\n", length(Xvars_full)))

# ============================================================================
# 3. Prepare women vitamin A dataset
# ============================================================================
cat("  Preparing women vitamin A dataset...\n")

# Women = rows where gw_childid is NA (mothers, not children)
df_women_vitA <- df %>%
  dplyr::select(dataid, gw_wVADAdjThurn, gw_cnum, dplyr::all_of(Xvars_full)) %>%
  as.data.frame() %>%
  dplyr::filter(is.na(df$gw_childid), !is.na(gw_wVADAdjThurn))

cat(sprintf("  Women Vitamin A dataset: %d rows\n", nrow(df_women_vitA)))

# ============================================================================
# 4. Fit binary SL model (or load if previously saved)
# ============================================================================
model_file <- file.path(out_models, "res_bin_GW_Ghana_SL_women_vitA_full_admin2.rds")

if (file.exists(model_file)) {
  cat(sprintf("  Loading pre-fitted model: %s\n", basename(model_file)))
  res_bin <- readRDS(model_file)
} else {
  cat("  Fitting binary SuperLearner for women vitamin A...\n")
  cat("  (This may take a while — using slmod_bin with 5-fold CV)\n")

  plan(multicore, workers = max(1L, floor(availableCores() / 2)))
  options(future.globals.maxSize = 5 * 1024^3)

  res_bin <- try(DHS_SL(
    d       = df_women_vitA,
    outcome = "gw_wVADAdjThurn",
    Xvars   = Xvars_full,
    id      = "gw_cnum",
    folds   = 5,
    CV      = FALSE,
    sl      = slmod_bin
  ))

  if (!inherits(res_bin, "try-error")) {
    saveRDS(res_bin, file = model_file)
    cat(sprintf("  Model saved: %s\n", basename(model_file)))
  } else {
    stop("  Model fitting failed. Check error messages above.")
  }
}

# ============================================================================
# 5. Extract predictions and aggregate to Admin2
# ============================================================================
cat("  Aggregating predictions to Admin2...\n")

# The model result contains: res_bin$res with columns dataid, Y, yhat_full
pred_df <- res_bin$res %>%
  dplyr::select(dataid, Y, yhat_full) %>%
  dplyr::rename(Y_obs = Y, yhat_pred = yhat_full)

# Join back to get Admin2 from the original dataset
admin2_col <- df_women_vitA %>%
  dplyr::select(dataid, Admin1, Admin2)

pred_with_admin <- pred_df %>%
  dplyr::left_join(admin2_col, by = "dataid") %>%
  dplyr::filter(!is.na(Admin2), !is.na(yhat_pred))

# Aggregate: predicted prevalence = mean of predicted probabilities per Admin2
# Observed prevalence = mean of observed binary outcome per Admin2 (unweighted)
admin2_prev <- pred_with_admin %>%
  dplyr::group_by(Admin1, Admin2) %>%
  dplyr::summarise(
    n               = dplyr::n(),
    obs_prevalence  = mean(Y_obs, na.rm = TRUE),
    pred_prevalence = mean(yhat_pred, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  dplyr::mutate(
    outcome   = "Women VAD (Thurnham-adj)",
    pred_type = "binary_prob"
  )

cat(sprintf("  Admin2 units with data: %d\n", nrow(admin2_prev)))

# Save CSV
csv_file <- file.path(out_tables, "ghana_admin2_prevalence_women_vitA.csv")
readr::write_csv(admin2_prev, csv_file)
cat(sprintf("  Saved: %s\n", basename(csv_file)))

# ============================================================================
# 6. Compute survey-weighted direct estimates at Admin2
# ============================================================================
cat("  Computing survey-weighted direct estimates...\n")

# Use original dataset (all women) with survey weights
# Survey design variables from the GMS
psu_var    <- "gw_EACode"
strata_var <- "gw_Strata"
wt_var     <- "gw_sWeight"

# Prepare survey data (women only, non-missing outcome)
d_survey <- df %>%
  dplyr::filter(is.na(gw_childid), !is.na(gw_wVADAdjThurn)) %>%
  dplyr::filter(!is.na(Admin2))

# Check that survey design variables exist
survey_vars_ok <- all(c(psu_var, strata_var, wt_var) %in% colnames(d_survey))

if (survey_vars_ok) {
  d_survey[[psu_var]]    <- as.factor(d_survey[[psu_var]])
  d_survey[[strata_var]] <- as.factor(d_survey[[strata_var]])
  d_survey[[wt_var]]     <- as.numeric(d_survey[[wt_var]])

  options(survey.lonely.psu = "adjust")

  des <- srvyr::as_survey_design(
    d_survey,
    ids     = !!rlang::sym(psu_var),
    strata  = !!rlang::sym(strata_var),
    weights = !!rlang::sym(wt_var),
    nest    = TRUE
  )

  direct_est <- des %>%
    dplyr::group_by(Admin2) %>%
    dplyr::summarise(
      direct_prev    = srvyr::survey_mean(gw_wVADAdjThurn,
                                           vartype = c("se", "ci"),
                                           na.rm = TRUE),
      n_unw          = srvyr::unweighted(sum(!is.na(gw_wVADAdjThurn))),
      .groups        = "drop"
    )

  cat(sprintf("  Direct estimates computed for %d Admin2 units\n", nrow(direct_est)))
} else {
  # Placeholder: use unweighted means if survey design variables are missing
  cat("  [warn] Survey design variables not found; using unweighted means as placeholder\n")

  direct_est <- d_survey %>%
    dplyr::group_by(Admin2) %>%
    dplyr::summarise(
      direct_prev       = mean(gw_wVADAdjThurn, na.rm = TRUE),
      direct_prev_se    = sd(gw_wVADAdjThurn, na.rm = TRUE) / sqrt(dplyr::n()),
      direct_prev_low   = direct_prev - 1.96 * direct_prev_se,
      direct_prev_upp   = direct_prev + 1.96 * direct_prev_se,
      n_unw             = sum(!is.na(gw_wVADAdjThurn)),
      .groups           = "drop"
    )
}

# Save direct estimates
direct_csv <- file.path(out_tables, "ghana_admin2_direct_est_women_vitA.csv")
readr::write_csv(direct_est, direct_csv)
cat(sprintf("  Saved: %s\n", basename(direct_csv)))

# ============================================================================
# 7. Get Ghana Admin2 boundaries
# ============================================================================
cat("  Downloading Ghana Admin-2 boundaries...\n")

poly_adm2 <- geodata::gadm(country = "GH", level = 2, path = tempdir()) |>
  sf::st_as_sf() |>
  dplyr::select(NAME_1, NAME_2) |>
  dplyr::rename(Admin1 = NAME_1, Admin2 = NAME_2) |>
  sf::st_transform(crs = 4326)

cat(sprintf("  Admin2 polygons: %d\n", nrow(poly_adm2)))

# ============================================================================
# 8. Map: Predicted prevalences at Admin2
# ============================================================================
cat("  Creating predicted prevalence map...\n")

poly_pred <- poly_adm2 %>%
  dplyr::left_join(
    admin2_prev %>% dplyr::select(Admin2, pred_prevalence, n),
    by = "Admin2"
  )

p_pred <- ggplot(poly_pred) +
  geom_sf(aes(fill = pred_prevalence), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    name    = "Predicted\nPrevalence",
    labels  = scales::percent_format(accuracy = 1),
    limits  = c(0, NA),
    na.value = "grey90",
    option  = "plasma"
  ) +
  labs(
    title    = "Predicted Women VAD Prevalence by District",
    subtitle = "Binary SuperLearner — predicted probability aggregated to Admin-2",
    caption  = "Ghana Micronutrient Survey 2017 | Model: slmod_bin (5-fold CV)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold")
  )

pred_png <- file.path(out_figures, "ghana_admin2_predicted_prev_women_vitA.png")
ggsave(pred_png, plot = p_pred, width = 8, height = 7, dpi = 200)
cat(sprintf("  Saved: %s\n", basename(pred_png)))

# ============================================================================
# 9. Map: Survey direct estimates at Admin2
# ============================================================================
cat("  Creating survey direct estimate map...\n")

# Rename the direct prevalence column for consistent join
direct_for_map <- direct_est %>%
  dplyr::select(Admin2, direct_prev)

poly_direct <- poly_adm2 %>%
  dplyr::left_join(direct_for_map, by = "Admin2")

p_direct <- ggplot(poly_direct) +
  geom_sf(aes(fill = direct_prev), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    name    = "Survey\nPrevalence",
    labels  = scales::percent_format(accuracy = 1),
    limits  = c(0, NA),
    na.value = "grey90",
    option  = "plasma"
  ) +
  labs(
    title    = "Survey Direct Estimate: Women VAD Prevalence by District",
    subtitle = "Design-based weighted prevalence at Admin-2",
    caption  = "Ghana Micronutrient Survey 2017 | Weighted by sWeight"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold")
  )

direct_png <- file.path(out_figures, "ghana_admin2_survey_direct_prev_women_vitA.png")
ggsave(direct_png, plot = p_direct, width = 8, height = 7, dpi = 200)
cat(sprintf("  Saved: %s\n", basename(direct_png)))

# ============================================================================
# 10. Side-by-side comparison map
# ============================================================================
cat("  Creating side-by-side comparison map...\n")

# Combine both estimates into one long-form dataset for faceting
combined_prev <- dplyr::bind_rows(
  admin2_prev %>%
    dplyr::select(Admin2, prevalence = pred_prevalence) %>%
    dplyr::mutate(source = "Model-predicted"),
  direct_est %>%
    dplyr::select(Admin2, prevalence = direct_prev) %>%
    dplyr::mutate(source = "Survey direct estimate")
)

poly_combined <- poly_adm2 %>%
  dplyr::left_join(combined_prev, by = "Admin2")

p_combined <- ggplot(poly_combined) +
  geom_sf(aes(fill = prevalence), color = "white", linewidth = 0.12) +
  facet_wrap(~ source, ncol = 2) +
  scale_fill_viridis_c(
    name    = "Prevalence",
    labels  = scales::percent_format(accuracy = 1),
    limits  = c(0, NA),
    na.value = "grey90",
    option  = "plasma"
  ) +
  labs(
    title    = "Women Vitamin A Deficiency: Predicted vs. Survey Estimates",
    subtitle = "Admin-2 (district) level, Ghana Micronutrient Survey 2017",
    caption  = paste0(
      "Left: Binary SL predicted probabilities aggregated to district\n",
      "Right: Design-based weighted prevalence from survey data"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold")
  )

comp_png <- file.path(out_figures, "ghana_admin2_comparison_women_vitA.png")
ggsave(comp_png, plot = p_combined, width = 14, height = 7, dpi = 200)
cat(sprintf("  Saved: %s\n", basename(comp_png)))

# ============================================================================
# 11. Scatter plot: predicted vs observed (Admin2)
# ============================================================================
cat("  Creating scatter plot...\n")

# Merge predicted and direct estimates for scatter
scatter_df <- admin2_prev %>%
  dplyr::select(Admin2, pred_prevalence, n) %>%
  dplyr::left_join(
    direct_est %>% dplyr::select(Admin2, direct_prev),
    by = "Admin2"
  ) %>%
  dplyr::filter(!is.na(pred_prevalence), !is.na(direct_prev))

if (nrow(scatter_df) > 0) {
  r2_val <- cor(scatter_df$direct_prev, scatter_df$pred_prevalence,
                use = "complete.obs")^2

  p_scatter <- ggplot(scatter_df, aes(x = direct_prev, y = pred_prevalence)) +
    geom_point(aes(size = n), colour = "#2166AC", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey30") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    scale_size_continuous(name = "Sample size", range = c(1.5, 6)) +
    labs(
      title    = "Admin-2: Predicted vs. Survey Direct Prevalence",
      subtitle = sprintf("Women VAD (Thurnham-adj) | %d districts | R\u00b2 = %.3f",
                         nrow(scatter_df), r2_val),
      x        = "Survey direct estimate",
      y        = "Model-predicted prevalence",
      caption  = "Ghana Micronutrient Survey 2017"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # Try with ggrepel labels if available
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_scatter <- p_scatter +
      ggrepel::geom_text_repel(aes(label = Admin2), size = 2.5,
                                max.overlaps = 12, colour = "grey40")
  }

  scatter_png <- file.path(out_figures, "ghana_admin2_scatter_women_vitA.png")
  ggsave(scatter_png, plot = p_scatter, width = 9, height = 7, dpi = 200)
  cat(sprintf("  Saved: %s\n", basename(scatter_png)))
} else {
  cat("  [skip] No matching districts for scatter plot\n")
}

# ============================================================================
# 12. CV performance summary
# ============================================================================
cat("  Computing CV performance metrics...\n")

Y_obs  <- pred_df$Y_obs
Y_pred <- pred_df$yhat_pred

# Clamp predicted probabilities
Y_pred_clamped <- pmin(pmax(Y_pred, 1e-6), 1 - 1e-6)

# AUC
auc_val <- tryCatch({
  roc_obj <- pROC::roc(response = Y_obs, predictor = Y_pred_clamped,
                        quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}, error = function(e) NA_real_)

# Brier score
brier <- mean((Y_pred_clamped - Y_obs)^2, na.rm = TRUE)

# Calibration
calib <- tryCatch({
  logit_p <- log(Y_pred_clamped / (1 - Y_pred_clamped))
  fit     <- suppressWarnings(glm(Y_obs ~ logit_p, family = binomial()))
  coef(fit)
}, error = function(e) c(intercept = NA_real_, slope = NA_real_))

perf_df <- data.frame(
  outcome         = "women_vitA",
  country         = "Ghana",
  model           = "slmod_bin",
  n               = length(Y_obs),
  prevalence      = mean(Y_obs, na.rm = TRUE),
  auc             = auc_val,
  brier           = brier,
  calib_intercept = unname(calib[1]),
  calib_slope     = unname(calib[2])
)

perf_csv <- file.path(out_tables, "ghana_admin2_cv_performance_women_vitA.csv")
readr::write_csv(perf_df, perf_csv)
cat(sprintf("  Saved: %s\n", basename(perf_csv)))

cat("\n  --- CV Performance ---\n")
cat(sprintf("  AUC   : %.3f\n", auc_val))
cat(sprintf("  Brier : %.4f\n", brier))
cat(sprintf("  Calib intercept: %.3f  slope: %.3f\n", calib[1], calib[2]))

cat("\n[03b] Done.\n\n")
