# =============================================================================
# scripts/national_prevalence_comparison.R
#
# Three-way comparison of national prevalence:
#   1. Observed (survey-weighted)
#   2. Within-country SL predicted (from conformal CI targets)
#   3. LOCO predicted (from transportability — trained on OTHER 3 countries)
#
# If LOCO pred_prev is not yet stored (older pipeline run), we approximate
# it from the admin1 predictions attribute or mark as NA.
# =============================================================================

library(targets)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

TARGETS_STORE <- here("_targets_full")
fig_dir <- here("results", "figures")

outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
outcome_labels <- c(
  child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
  child_iron = "Iron (children)", women_iron = "Iron (women)")

countries_lc <- c(Gambia = "gambia", Ghana = "ghana",
                   SierraLeone = "sierraleone", Malawi = "malawi")
country_display <- c(Gambia = "Gambia", Ghana = "Ghana",
                      SierraLeone = "Sierra Leone", Malawi = "Malawi")

rows <- list()

for (oc in outcomes) {
  # ── LOCO data ──
  loco <- tryCatch(
    targets::tar_read_raw(paste0("loco_", oc), store = TARGETS_STORE),
    error = function(e) NULL)
  if (is.null(loco)) next

  for (i in seq_len(nrow(loco))) {
    ho <- loco$held_out[i]
    ctry_lc <- countries_lc[ho]

    # 1. Observed national prevalence
    obs_prev <- loco$prev_test[i]

    # 2. Within-country SL predicted national prevalence
    within_pred <- NA_real_
    ci <- tryCatch(
      targets::tar_read_raw(paste0("conformal_ci_", ctry_lc, "_", oc),
                             store = TARGETS_STORE),
      error = function(e) NULL)
    if (!is.null(ci) && !is.null(ci$national_ci)) {
      nc <- ci$national_ci
      if ("boot_mean" %in% names(nc) && nc$boot_mean[1] >= 0 && nc$boot_mean[1] <= 1) {
        within_pred <- nc$boot_mean[1]
      }
    }

    # 3. LOCO predicted national prevalence
    loco_pred <- NA_real_
    # Try the new pred_prev column (available after transportability.R update)
    if ("pred_prev" %in% names(loco)) {
      loco_pred <- loco$pred_prev[i]
    }
    # Fallback: try admin1 attribute (only works for last held-out country)
    if (is.na(loco_pred)) {
      a1 <- attr(loco, "admin1_preds")
      if (!is.null(a1) && is.data.frame(a1) && "pred_prev" %in% names(a1) && "n" %in% names(a1)) {
        # This is only valid if this row corresponds to the last held-out country
        if (i == nrow(loco)) {
          loco_pred <- weighted.mean(a1$pred_prev, a1$n, na.rm = TRUE)
        }
      }
    }

    rows[[paste0(oc, "_", ho)]] <- data.frame(
      outcome = oc,
      outcome_label = outcome_labels[oc],
      held_out = ho,
      country = country_display[ho],
      obs_prev = obs_prev,
      within_pred = within_pred,
      loco_pred = loco_pred,
      loco_auc = loco$auc[i],
      stringsAsFactors = FALSE
    )
  }
}

df <- bind_rows(rows)

cat("=== National Prevalence Comparison ===\n")
print(df |> select(country, outcome_label, obs_prev, within_pred, loco_pred, loco_auc) |>
        mutate(across(where(is.numeric), ~round(., 4))))

# ── Pivot to long format for plotting ────────────────────────────────────
df_long <- df |>
  pivot_longer(
    cols = c(obs_prev, within_pred, loco_pred),
    names_to = "estimate_type",
    values_to = "prevalence"
  ) |>
  mutate(
    estimate_label = case_when(
      estimate_type == "obs_prev"    ~ "Observed (survey)",
      estimate_type == "within_pred" ~ "Within-country SL",
      estimate_type == "loco_pred"   ~ "LOCO (cross-country)"
    ),
    estimate_label = factor(estimate_label,
                             levels = c("Observed (survey)",
                                        "Within-country SL",
                                        "LOCO (cross-country)"))
  ) |>
  filter(!is.na(prevalence))

# ── Figure 1: Grouped bar chart ──────────────────────────────────────────
p1 <- ggplot(df_long, aes(x = country, y = prevalence, fill = estimate_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    name = NULL,
    values = c("Observed (survey)" = "#2c7bb6",
               "Within-country SL" = "#abd9e9",
               "LOCO (cross-country)" = "#d7191c")
  ) +
  facet_wrap(~ outcome_label, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "National prevalence",
    title = "National Prevalence: Survey vs. Within-Country vs. Cross-Country Prediction",
    subtitle = "Blue = survey | Light blue = within-country SL | Red = LOCO (trained on other 3 countries)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p1)
ggsave(file.path(fig_dir, "national_prev_3way_bars.png"), p1,
       width = 11, height = 8, dpi = 200, bg = "white")
cat("Saved national_prev_3way_bars.png\n")


# ── Figure 2: Scatter — observed vs predicted (within + LOCO) ────────────
df_scatter <- df |>
  filter(!is.na(within_pred)) |>
  pivot_longer(
    cols = c(within_pred, loco_pred),
    names_to = "model_type",
    values_to = "pred_prev"
  ) |>
  mutate(
    model_label = ifelse(model_type == "within_pred",
                          "Within-country SL", "LOCO (cross-country)"),
    model_label = factor(model_label,
                          levels = c("Within-country SL", "LOCO (cross-country)"))
  ) |>
  filter(!is.na(pred_prev))

p2 <- ggplot(df_scatter, aes(x = obs_prev, y = pred_prev)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = country, shape = model_label), size = 3.5, alpha = 0.8) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_colour_brewer(name = "Country", palette = "Set2") +
  scale_shape_manual(name = "Model", values = c("Within-country SL" = 16,
                                                   "LOCO (cross-country)" = 1)) +
  facet_wrap(~ outcome_label, scales = "free") +
  labs(
    x = "Observed national prevalence (survey)",
    y = "Predicted national prevalence",
    title = "National Prevalence Recovery: Within-Country vs. LOCO",
    subtitle = "Filled = within-country SL | Open = LOCO (trained on other 3 countries) | Dashed = perfect"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold")
  )

print(p2)
ggsave(file.path(fig_dir, "national_prev_3way_scatter.png"), p2,
       width = 10, height = 8, dpi = 200, bg = "white")
cat("Saved national_prev_3way_scatter.png\n")


# ── Figure 3: Forest plot — all three side by side ───────────────────────
p3 <- ggplot(df_long, aes(x = prevalence,
                            y = reorder(paste0(country, " — ", outcome_label), prevalence),
                            colour = estimate_label, shape = estimate_label)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_colour_manual(
    name = NULL,
    values = c("Observed (survey)" = "#2c7bb6",
               "Within-country SL" = "#abd9e9",
               "LOCO (cross-country)" = "#d7191c")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Observed (survey)" = 15,
               "Within-country SL" = 16,
               "LOCO (cross-country)" = 4)
  ) +
  labs(
    x = "National prevalence",
    y = NULL,
    title = "National Prevalence: Three-Way Comparison",
    subtitle = "Square = survey | Circle = within-country SL | X = LOCO (cross-country)"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

print(p3)
ggsave(file.path(fig_dir, "national_prev_3way_forest.png"), p3,
       width = 10, height = 9, dpi = 200, bg = "white")
cat("Saved national_prev_3way_forest.png\n")

# ── Save ggplot objects ──────────────────────────────────────────────────
saveRDS(list(
  bars = p1,
  scatter = p2,
  forest = p3,
  data = df
), file.path(fig_dir, "national_prev_3way_figures.rds"))
cat("Saved national_prev_3way_figures.rds\n")
cat("Done.\n")
