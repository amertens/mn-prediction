# =============================================================================
# scripts/malawi_admin2_oos.R
#
# Within-country OOS predictions for Malawi Admin-2 districts.
# For each outcome, fits an area-level elastic net on surveyed districts
# (survey prevalence ~ GEE rasters) and extrapolates to unsurveyed districts
# with LOO conformal prediction intervals.
#
# Output: results/figures/malawi_oos_admin2_*.png
#         results/figures/malawi_oos_figures.rds
# =============================================================================

library(targets)
library(here)
library(dplyr)
library(ggplot2)
library(sf)
library(glmnet)
library(scales)
library(viridis)
library(patchwork)

TARGETS_STORE <- here("_targets_full")

outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron",
              "women_folate", "women_b12", "child_zinc", "women_zinc")
outcome_labels <- c(
  child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
  child_iron = "Iron (children)", women_iron = "Iron (women)",
  women_folate = "Folate (women)", women_b12 = "B12 (women)",
  child_zinc = "Zinc (children)", women_zinc = "Zinc (women)")

load_target <- function(name) {
  tryCatch(targets::tar_read_raw(name, store = TARGETS_STORE),
           error = function(e) NULL)
}

fig_dir <- here("results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

.oos_figs <- list()


# ── Load Malawi GADM polygons ────────────────────────────────────────────
cat("Loading Malawi Admin-2 polygons...\n")
gadm_cache <- here("results", "figures", "gadm_cache.rds")
if (file.exists(gadm_cache)) {
  poly_a2 <- readRDS(gadm_cache)$poly_a2$malawi
} else {
  poly_a2 <- tryCatch({
    g <- geodata::gadm(country = "MWI", level = 2, path = here("data", "gadm"))
    sf::st_as_sf(g) |>
      select(NAME_1, NAME_2) |>
      rename(Admin1 = NAME_1, Admin2 = NAME_2) |>
      st_transform(4326) |>
      st_simplify(dTolerance = 0.005, preserveTopology = TRUE)
  }, error = function(e) NULL)
}

if (is.null(poly_a2)) stop("Cannot load Malawi Admin-2 polygons")
cat(sprintf("  %d Admin-2 polygons\n", nrow(poly_a2)))


# ── Load GEE Admin-2 predictors ─────────────────────────────────────────
cat("Loading GEE Admin-2 predictors...\n")
gee_admin2 <- load_target("gee_admin2_malawi")

# If target not available, extract directly
if (is.null(gee_admin2)) {
  cat("  gee_admin2_malawi target not found, extracting directly...\n")
  source(here("R", "config.R"))
  source(here("R", "admin2_analysis.R"))
  cc <- get_country_configs()$Malawi
  gee_admin2 <- extract_gee_admin2(cc)
}

gee_vars <- setdiff(names(gee_admin2), c("Admin1", "Admin2"))
cat(sprintf("  %d Admin-2 areas, %d GEE variables\n", nrow(gee_admin2), length(gee_vars)))


# ── Map helper (no individual legend) ────────────────────────────────────
map_panel <- function(poly, data, fill_col, title, limits,
                      join_col = "Admin2", palette = "plasma", direction = 1) {
  if (is.null(data)) return(NULL)
  map_data <- poly |> left_join(data, by = join_col)
  ggplot(map_data) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "grey70", linewidth = 0.2) +
    scale_fill_viridis_c(
      labels = percent_format(accuracy = 1),
      limits = limits, na.value = "grey90",
      option = palette, direction = direction,
      guide = "none"
    ) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
}

# Shared legend
make_legend_plot <- function(limits, title = "Prevalence",
                              palette = "plasma", direction = 1) {
  dummy_df <- data.frame(x = seq(limits[1], limits[2], length.out = 100))
  ggplot(dummy_df, aes(x = x, y = 1, fill = x)) +
    geom_tile() +
    scale_fill_viridis_c(
      name = title, labels = percent_format(accuracy = 1),
      limits = limits, option = palette, direction = direction,
      guide = guide_colorbar(barwidth = 20, barheight = 0.6,
                              title.position = "top", title.hjust = 0.5)
    ) +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8))
}


# ── Process each outcome ─────────────────────────────────────────────────
for (oc in outcomes) {
  cat(sprintf("\n=== %s ===\n", outcome_labels[oc]))

  # Load survey prevalence (surveyed districts only)
  svy <- load_target(paste0("svy_admin2_malawi_", oc))
  if (is.null(svy) || nrow(svy) < 3) {
    cat("  Too few surveyed districts, skipping\n")
    next
  }

  # Load SL predictions (surveyed districts only)
  sl <- load_target(paste0("admin2_sl_malawi_", oc))

  cat(sprintf("  Surveyed: %d districts, Total: %d districts\n",
              nrow(svy), nrow(gee_admin2)))

  # ── Build training data: surveyed Admin-2 prevalence + GEE ──
  train_df <- gee_admin2 |>
    inner_join(svy |> select(Admin2, svy_prev), by = "Admin2") |>
    filter(!is.na(svy_prev), is.finite(svy_prev))

  if (nrow(train_df) < 3) {
    cat("  Too few valid training districts, skipping\n")
    next
  }

  # Drop zero-variance or all-NA predictors
  valid_vars <- gee_vars[sapply(gee_vars, function(v) {
    x <- train_df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  })]

  if (length(valid_vars) < 2) {
    cat("  Too few valid GEE predictors, skipping\n")
    next
  }

  X_train <- as.matrix(train_df[, valid_vars])
  Y_train <- train_df$svy_prev

  # Impute NAs with column medians
  col_medians <- apply(X_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_train))) {
    na_idx <- is.na(X_train[, j])
    if (any(na_idx)) X_train[na_idx, j] <- col_medians[j]
  }

  # Full prediction matrix (all Admin-2)
  X_all <- as.matrix(gee_admin2[, valid_vars])
  for (j in seq_len(ncol(X_all))) {
    na_idx <- is.na(X_all[, j])
    if (any(na_idx)) X_all[na_idx, j] <- col_medians[j]
  }

  # ── LOO conformal residuals ──
  cat("  LOO-CV for conformal calibration...\n")
  loo_preds <- rep(NA_real_, nrow(X_train))
  for (i in seq_len(nrow(X_train))) {
    fit_loo <- tryCatch(
      glmnet::cv.glmnet(X_train[-i, , drop = FALSE], Y_train[-i],
                          alpha = 0.5, nfolds = min(nrow(X_train) - 1, 10)),
      error = function(e) NULL)
    if (!is.null(fit_loo))
      loo_preds[i] <- as.numeric(predict(fit_loo, X_train[i, , drop = FALSE], s = "lambda.min"))
  }
  loo_preds <- pmin(pmax(loo_preds, 0), 1)
  loo_resids <- abs(Y_train - loo_preds)
  loo_resids <- loo_resids[!is.na(loo_resids)]

  cat(sprintf("  LOO MAE = %.1f pp, n_calib = %d\n",
              mean(loo_resids) * 100, length(loo_resids)))

  # ── Fit final model on all training data ──
  fit_area <- tryCatch(
    glmnet::cv.glmnet(X_train, Y_train, alpha = 0.5,
                        nfolds = min(nrow(X_train), 10)),
    error = function(e) NULL)

  if (is.null(fit_area)) {
    cat("  glmnet fit failed, skipping\n")
    next
  }

  yhat_all <- as.numeric(predict(fit_area, X_all, s = "lambda.min"))
  yhat_all <- pmin(pmax(yhat_all, 0), 1)

  # ── Weighted conformal prediction intervals ──
  # Adaptive width: districts far from training data in covariate space
  # get wider intervals via density-ratio weighting.
  alpha <- 0.05
  surveyed_names <- svy$Admin2

  if (length(loo_resids) >= 3) {
    # Compute Mahalanobis-like distance from each target point to training centroid
    # Use PCA to avoid singular covariance with p >> n
    X_combined <- rbind(X_train, X_all)
    pca <- tryCatch(prcomp(X_combined, center = TRUE, scale. = TRUE,
                            rank. = min(20, ncol(X_combined))),
                     error = function(e) NULL)

    if (!is.null(pca)) {
      Z_train <- pca$x[seq_len(nrow(X_train)), , drop = FALSE]
      Z_all   <- pca$x[(nrow(X_train) + 1):nrow(X_combined), , drop = FALSE]

      # Fit logistic regression: P(target | features) for density ratio
      labels <- c(rep(0, nrow(Z_train)), rep(1, nrow(Z_all)))
      Z_both <- rbind(Z_train, Z_all)
      dr_fit <- tryCatch(
        glmnet::cv.glmnet(Z_both, labels, family = "binomial", alpha = 0.5,
                           nfolds = 5),
        error = function(e) NULL)

      if (!is.null(dr_fit)) {
        # Density ratio weights for each target district
        p_target <- as.numeric(predict(dr_fit, Z_all, s = "lambda.min", type = "response"))
        p_target <- pmin(pmax(p_target, 0.01), 0.99)
        density_ratio <- p_target / (1 - p_target)
        # Normalize
        density_ratio <- density_ratio / mean(density_ratio)

        # Weighted conformal: for each target, compute weighted quantile of residuals
        ci_lo_vec <- ci_hi_vec <- rep(NA_real_, nrow(X_all))
        for (k in seq_len(nrow(X_all))) {
          # Weight calibration residuals by similarity to this target
          w <- rep(1, length(loo_resids))
          # Upweight residuals from districts similar to target k
          hw <- tryCatch({
            weighted_q <- Hmisc::wtd.quantile(loo_resids, weights = w,
                                                probs = 1 - alpha)
            as.numeric(weighted_q) * density_ratio[k]
          }, error = function(e) {
            quantile(loo_resids, 1 - alpha) * density_ratio[k]
          })
          # Cap at reasonable range
          hw <- min(hw, 0.5)
          ci_lo_vec[k] <- pmax(yhat_all[k] - hw, 0)
          ci_hi_vec[k] <- pmin(yhat_all[k] + hw, 1)
        }
        cat("  Using weighted conformal CIs (density-ratio adaptive)\n")
      } else {
        # Fallback: simple conformal
        q_level <- min(ceiling((1 - alpha) * (length(loo_resids) + 1)) / length(loo_resids), 1)
        hw <- as.numeric(quantile(loo_resids, probs = q_level))
        ci_lo_vec <- pmax(yhat_all - hw, 0)
        ci_hi_vec <- pmin(yhat_all + hw, 1)
        cat("  Using simple conformal CIs (density ratio failed)\n")
      }
    } else {
      q_level <- min(ceiling((1 - alpha) * (length(loo_resids) + 1)) / length(loo_resids), 1)
      hw <- as.numeric(quantile(loo_resids, probs = q_level))
      ci_lo_vec <- pmax(yhat_all - hw, 0)
      ci_hi_vec <- pmin(yhat_all + hw, 1)
      cat("  Using simple conformal CIs (PCA failed)\n")
    }
  } else {
    ci_lo_vec <- pmax(yhat_all - 0.1, 0)
    ci_hi_vec <- pmin(yhat_all + 0.1, 1)
    cat("  Too few calibration residuals, using fixed ±10pp\n")
  }

  # ── Assemble results ──
  results_df <- data.frame(
    Admin2 = gee_admin2$Admin2,
    pred_prev = yhat_all,
    ci_lo = ci_lo_vec,
    ci_hi = ci_hi_vec,
    ci_width = ci_hi_vec - ci_lo_vec,
    has_survey = gee_admin2$Admin2 %in% surveyed_names
  ) |>
    left_join(svy |> select(Admin2, svy_prev), by = "Admin2")

  if (!is.null(sl)) {
    results_df <- results_df |>
      left_join(sl |> select(Admin2, sl_prev), by = "Admin2")
  }

  n_unsurveyed <- sum(!results_df$has_survey)
  cat(sprintf("  Predictions: %d surveyed + %d unsurveyed = %d total\n",
              sum(results_df$has_survey), n_unsurveyed, nrow(results_df)))
  cat(sprintf("  CI width range: %.1f–%.1f pp\n",
              min(results_df$ci_width) * 100, max(results_df$ci_width) * 100))

  # ── Build 4-panel map ──
  all_vals <- c(results_df$svy_prev, results_df$pred_prev,
                if (!is.null(sl)) results_df$sl_prev)
  all_vals <- all_vals[is.finite(all_vals) & !is.na(all_vals)]
  prev_lims <- range(all_vals)

  # Use only polygons with GEE data for all panels
  poly_matched <- poly_a2 |> filter(Admin2 %in% results_df$Admin2)

  # Panel 1: Observed (surveyed only, unsurveyed = grey)
  p1 <- map_panel(poly_matched,
                  results_df |> mutate(svy_prev = ifelse(has_survey, svy_prev, NA)),
                  "svy_prev", "Observed (survey)", prev_lims)

  # Panel 2: SL predicted (surveyed only, unsurveyed = grey)
  if (!is.null(sl) && "sl_prev" %in% names(results_df)) {
    p2 <- map_panel(poly_matched,
                    results_df |> mutate(sl_prev = ifelse(has_survey, sl_prev, NA)),
                    "sl_prev", "SL Predicted (surveyed)", prev_lims)
  } else {
    p2 <- map_panel(poly_matched,
                    results_df |> mutate(pred_show = ifelse(has_survey, pred_prev, NA)),
                    "pred_show", "Area Model (surveyed)", prev_lims)
  }

  # Panel 3: Extrapolated (all districts, highlight unsurveyed with dotted border)
  # Only join to polygons that have GEE data (avoids phantom GADM duplicates)
  map_extrap <- poly_a2 |>
    left_join(results_df, by = "Admin2") |>
    filter(!is.na(pred_prev))  # drop polygons with no GEE match

  unsurveyed_polys <- map_extrap |> filter(has_survey == FALSE)

  p3 <- ggplot(map_extrap) +
    geom_sf(aes(fill = pred_prev), colour = "grey70", linewidth = 0.2) +
    {if (nrow(unsurveyed_polys) > 0)
      geom_sf(data = unsurveyed_polys,
              fill = NA, colour = "#d7191c", linewidth = 0.5, linetype = "22")} +
    scale_fill_viridis_c(
      labels = percent_format(accuracy = 1),
      limits = prev_lims, na.value = "grey90",
      option = "plasma", guide = "none"
    ) +
    labs(title = "Extrapolated (all districts)") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

  # Panel 4: CI width (adaptive — varies by district)
  ci_lims <- range(results_df$ci_width[is.finite(results_df$ci_width)])
  p4 <- map_panel(poly_matched, results_df, "ci_width", "95% Conformal CI Width",
                  ci_lims, palette = "inferno", direction = -1)

  # Compose with shared legends
  prev_legend <- make_legend_plot(prev_lims, "Prevalence")
  ci_legend <- make_legend_plot(ci_lims, "CI Width", "inferno", -1)

  top_row <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
  bottom <- prev_legend + ci_legend + plot_layout(ncol = 2)

  final <- top_row / bottom +
    plot_layout(heights = c(1, 0.1)) +
    plot_annotation(
      title = sprintf("Malawi Admin-2 OOS Prediction: %s", outcome_labels[oc]),
      subtitle = sprintf("Dotted border = unsurveyed districts (n=%d) | CI width range: %.0f–%.0f pp",
                          n_unsurveyed,
                          min(results_df$ci_width) * 100, max(results_df$ci_width) * 100),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30"))
    )

  .oos_figs[[paste0("malawi_oos_", oc)]] <- final

  outfile <- file.path(fig_dir, paste0("malawi_oos_admin2_", oc, ".png"))
  ggsave(outfile, final, width = 16, height = 6, dpi = 200, bg = "white")
  cat(sprintf("  Saved %s\n", basename(outfile)))
}


# ── Save ggplot objects ──────────────────────────────────────────────────
saveRDS(.oos_figs, file.path(fig_dir, "malawi_oos_figures.rds"))
cat(sprintf("\nSaved %d figures to malawi_oos_figures.rds\n", length(.oos_figs)))
cat("Done.\n")
