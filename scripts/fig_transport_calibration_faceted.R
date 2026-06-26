# =============================================================================
# scripts/fig_transport_calibration_faceted.R
#
# Faceted transport-calibration figure: one panel per outcome, free scales so
# each outcome's prevalence range is visible. Each panel shows predicted vs.
# true national prevalence across LOCO folds, the 45-degree perfect-calibration
# line, residual drop-lines, and a per-outcome fitted calibration slope (in the
# strip label). Outcome-agnostic: facets whatever outcomes are in the CSV, so it
# picks up folate/B12 automatically once those folds are added.
#
# Input : results/tables/transportability_loco_results.csv
# Output: results/figures/transport_calibration_loco_faceted.png
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(dplyr) })

d <- read.csv("results/tables/transportability_loco_results.csv", stringsAsFactors = FALSE)
d <- d[is.finite(d$prev_test) & is.finite(d$pred_prev), , drop = FALSE]
d$true_pct <- d$prev_test * 100
d$pred_pct <- d$pred_prev * 100
d$country_lab <- gsub("([a-z])([A-Z])", "\\1 \\2", d$held_out)

oc_lab <- c(child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
            child_iron = "Iron (children)",      women_iron = "Iron (women)",
            women_folate = "Folate (women)",     women_b12 = "Vitamin B12 (women)",
            child_zinc = "Zinc (children)",      women_zinc = "Zinc (women)")
d$outcome_pretty <- ifelse(d$outcome %in% names(oc_lab), unname(oc_lab[d$outcome]), d$outcome)

# Per-outcome OLS calibration slope (perfect = 1); shown in the strip label.
slopes <- d |>
  group_by(outcome, outcome_pretty) |>
  summarise(slope = if (n() >= 2 && sd(true_pct) > 0)
              unname(coef(lm(pred_pct ~ true_pct))[2]) else NA_real_,
            n = n(), .groups = "drop")
slopes$strip <- ifelse(is.na(slopes$slope),
                       sprintf("%s  (n=%d)", slopes$outcome_pretty, slopes$n),
                       sprintf("%s  —  calib. slope = %.2f", slopes$outcome_pretty, slopes$slope))
d <- left_join(d, slopes[, c("outcome", "strip")], by = "outcome")
# Order panels by ascending slope (worst calibration first) for narrative flow.
ord <- slopes$strip[order(slopes$slope, na.last = TRUE)]
d$strip <- factor(d$strip, levels = ord)

ctry_colors <- c("Gambia" = "#1b9e77", "Ghana" = "#d95f02",
                 "Sierra Leone" = "#7570b3", "Malawi" = "#e7298a")

p <- ggplot(d, aes(true_pct, pred_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey45", linewidth = 0.5) +
  geom_segment(aes(xend = true_pct, yend = true_pct, color = country_lab),
               linewidth = 0.4, alpha = 0.55) +
  geom_point(aes(color = country_lab), size = 3) +
  geom_text_repel(aes(label = country_lab, color = country_lab),
                  size = 2.7, show.legend = FALSE, max.overlaps = 20,
                  seed = 1, min.segment.length = 0) +
  facet_wrap(~ strip, scales = "free", ncol = 2) +
  scale_color_manual(values = ctry_colors, name = "Held-out country") +
  labs(
    title = "Cross-country transport calibration by outcome",
    subtitle = "Predicted vs. true national prevalence, leave-one-country-out folds. Dashed = perfect calibration (y = x).",
    x = "True national prevalence (survey, %)",
    y = "Predicted national prevalence (%)",
    caption = "Each panel: one outcome, one point per held-out country (trained on the others). Slope < 1 = regression toward the training-set mean."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey40", hjust = 0),
    strip.text = element_text(face = "bold", size = 9.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

n_oc <- length(unique(d$strip))
out <- "results/figures/transport_calibration_loco_faceted.png"
ggsave(out, p, width = 9, height = 3.6 * ceiling(n_oc / 2) + 0.6, dpi = 200, bg = "white")
cat(sprintf("Wrote %s  (%d outcome panels)\n", out, n_oc))
print(slopes[, c("outcome_pretty", "slope", "n")])
