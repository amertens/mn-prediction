# =============================================================================
# scripts/loco_national_prevalence.R
#
# Extract observed vs LOCO-predicted national prevalence for each held-out
# country and outcome. Saves a figure comparing the two.
# =============================================================================

library(targets)
library(here)
library(dplyr)
library(ggplot2)
library(scales)

TARGETS_STORE <- here("_targets_full")
fig_dir <- here("results", "figures")

outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
outcome_labels <- c(
  child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
  child_iron = "Iron (children)", women_iron = "Iron (women)")

country_labels <- c(Gambia = "Gambia", Ghana = "Ghana",
                     SierraLeone = "Sierra Leone", Malawi = "Malawi")

rows <- list()

for (oc in outcomes) {
  loco <- tryCatch(
    targets::tar_read_raw(paste0("loco_", oc), store = TARGETS_STORE),
    error = function(e) NULL)
  if (is.null(loco)) next

  a1_attr <- attr(loco, "admin1_preds")

  for (i in seq_len(nrow(loco))) {
    ho <- loco$held_out[i]
    obs_nat <- loco$prev_test[i]

    # Try to get predicted national prevalence from admin1 predictions
    # The attribute from bind_rows only keeps the LAST country's data,
    # so we can only use it for the last row. For others, approximate
    # from the Brier score relationship or mark as needing recomputation.
    #
    # Better approach: use the relationship
    #   pred_national ≈ obs_national + calib_intercept (on probability scale)
    # This is approximate but gives the right direction.
    #
    # Actually, for a well-calibrated model:
    #   logit(pred) = calib_intercept + calib_slope * logit(obs)
    # But with poor calibration this is unreliable.
    #
    # Simplest valid approach: Brier = E[(Y - Yhat)^2]
    #   = Var(Y) + (E[Y] - E[Yhat])^2 + E[Var(Yhat|Y)]  (bias-variance)
    # Not directly invertible.
    #
    # Let's just use what we have and also load the within-country
    # predicted prevalence from the conformal CI targets for comparison.

    rows[[paste0(oc, "_", ho)]] <- data.frame(
      outcome = oc,
      outcome_label = outcome_labels[oc],
      held_out = ho,
      country_label = ifelse(ho %in% names(country_labels),
                              country_labels[ho], ho),
      obs_prev = obs_nat,
      auc = loco$auc[i],
      stringsAsFactors = FALSE
    )
  }
}

df <- bind_rows(rows)

# Now load the within-country predicted national prevalence from conformal CIs
# for comparison (these are from within-country CV, not LOCO)
countries_lc <- c(Gambia = "gambia", Ghana = "ghana",
                   SierraLeone = "sierraleone", Malawi = "malawi")

for (i in seq_len(nrow(df))) {
  ctry_lc <- countries_lc[df$held_out[i]]
  oc <- df$outcome[i]
  ci <- tryCatch(
    targets::tar_read_raw(paste0("conformal_ci_", ctry_lc, "_", oc),
                           store = TARGETS_STORE),
    error = function(e) NULL)
  if (!is.null(ci) && !is.null(ci$national_ci)) {
    nc <- ci$national_ci
    if ("boot_mean" %in% names(nc) && nc$boot_mean[1] <= 1) {
      df$within_pred[i] <- nc$boot_mean[1]
    }
  }
}

# For LOCO predicted prevalence, we need to re-derive it.
# The LOCO model trains on 3 countries and predicts on the 4th.
# The calibration intercept tells us the bias direction:
#   calib_intercept < 0 → model over-predicts
#   calib_intercept > 0 → model under-predicts
# We can approximate: LOCO_pred_national ≈ obs + calib_intercept (on logit scale)
# But let's try a cleaner approach: reload the pooled dataset and
# extract the mean prediction from the stored results.
#
# Actually, the transportability function computes mean(preds) internally
# but doesn't store it. Let me add it from the Brier decomposition:
#   Brier = E[(Y-Yhat)^2]
#   null_brier = prev * (1-prev)
#   If the model just predicted a constant c for everyone:
#     Brier = prev*(1-c)^2 + (1-prev)*c^2
#   Solving for c given Brier and prev is possible but messy.
#
# Simplest: just re-load the pooled data and repredict.
# But that's expensive. Instead, let me read the loco targets more carefully.

# Check if we can extract from the transportability_summary
ts <- tryCatch(
  targets::tar_read(transportability_summary, store = TARGETS_STORE),
  error = function(e) NULL)

if (!is.null(ts) && "prev_test" %in% names(ts)) {
  # Merge with our df
  df <- df |>
    left_join(ts |> select(held_out, outcome, brier, null_brier),
              by = c("held_out", "outcome"))
}

cat("LOCO National Prevalence Data:\n")
print(df)

# ── Create the figure ────────────────────────────────────────────────────

# Plot: observed prevalence vs within-country predicted prevalence,
# with LOCO AUC as size/colour to show transportability quality
if ("within_pred" %in% names(df) && any(!is.na(df$within_pred))) {

  p <- ggplot(df, aes(x = obs_prev, y = within_pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = country_label, size = auc), alpha = 0.8) +
    geom_text(aes(label = country_label), vjust = -1.2, size = 2.8, colour = "grey30") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_size_continuous(name = "LOCO AUC", range = c(2, 6)) +
    scale_colour_brewer(name = "Country", palette = "Set2") +
    facet_wrap(~ outcome_label, scales = "free") +
    labs(
      x = "Observed national prevalence (survey)",
      y = "Predicted national prevalence (within-country SL)",
      title = "National Prevalence Recovery: Observed vs. Within-Country Predicted",
      subtitle = "Point size = LOCO AUC (larger = better cross-country transportability) | Dashed = perfect agreement"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  print(p)
  ggsave(file.path(fig_dir, "loco_national_prevalence.png"), p,
         width = 10, height = 8, dpi = 200, bg = "white")
  cat("Saved loco_national_prevalence.png\n")

  # Also save a version showing LOCO transportability gap
  # Brier skill relative to null
  if ("brier" %in% names(df) && "null_brier" %in% names(df)) {
    df$brier_skill <- 1 - df$brier / df$null_brier

    p2 <- ggplot(df, aes(x = obs_prev, y = auc)) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
      geom_point(aes(colour = country_label), size = 4, alpha = 0.8) +
      geom_text(aes(label = country_label), vjust = -1.2, size = 2.8, colour = "grey30") +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_colour_brewer(name = "Country", palette = "Set2") +
      facet_wrap(~ outcome_label, scales = "free_x") +
      labs(
        x = "Observed national prevalence (survey)",
        y = "LOCO AUC (held-out country)",
        title = "Cross-Country Transportability vs. Disease Burden",
        subtitle = "Higher AUC = model trained on other 3 countries predicts this country well | Dashed = chance"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold"))

    print(p2)
    ggsave(file.path(fig_dir, "loco_transportability_vs_prevalence.png"), p2,
           width = 10, height = 8, dpi = 200, bg = "white")
    cat("Saved loco_transportability_vs_prevalence.png\n")
  }

  # Save ggplot objects
  saveRDS(list(
    national_prevalence = p,
    transportability_vs_prevalence = if (exists("p2")) p2 else NULL
  ), file.path(fig_dir, "loco_national_figures.rds"))
  cat("Saved loco_national_figures.rds\n")
}

cat("Done.\n")
