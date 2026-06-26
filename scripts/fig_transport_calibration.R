# =============================================================================
# scripts/fig_transport_calibration.R
#
# Transport-calibration figure: predicted vs. true NATIONAL prevalence across
# leave-one-country-out (LOCO) folds. Each point is one held-out country ×
# outcome; the model is trained on the other three countries. Distance from the
# 45° line is the calibration error. The shallow fitted slope (< 1) is the
# headline finding: transported models regress toward the training-set mean
# prevalence — over-predicting rare outcomes, under-predicting common ones.
#
# Input : results/tables/transportability_loco_results.csv
# Output: results/figures/transport_calibration_loco.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel)
})

inp <- "results/tables/transportability_loco_results.csv"
d <- read.csv(inp, stringsAsFactors = FALSE)

# True (survey) vs predicted national prevalence, in percent.
d$true_pct <- d$prev_test * 100
d$pred_pct <- d$pred_prev * 100

# Pretty labels (covers every outcome that may appear in the LOCO table).
oc_lab <- c(child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
            child_iron = "Iron (children)",      women_iron = "Iron (women)",
            women_folate = "Folate (women)",     women_b12 = "Vitamin B12 (women)",
            child_zinc = "Zinc (children)",      women_zinc = "Zinc (women)")
present <- oc_lab[names(oc_lab) %in% d$outcome]
d$outcome_lab <- factor(unname(oc_lab[d$outcome]), levels = unname(present))
d$country_lab <- gsub("([a-z])([A-Z])", "\\1 \\2", d$held_out)

# Overall calibration fit across all folds (the regression-to-mean line).
fit   <- lm(pred_pct ~ true_pct, data = d)
slope <- unname(coef(fit)[2]); intc <- unname(coef(fit)[1])

lim <- c(0, max(d$true_pct, d$pred_pct) * 1.08)

oc_colors <- c("Vitamin A (children)" = "#E69F00", "Vitamin A (women)" = "#D55E00",
               "Iron (children)" = "#0072B2", "Iron (women)" = "#56B4E9",
               "Folate (women)" = "#009E73", "Vitamin B12 (women)" = "#CC79A7",
               "Zinc (children)" = "#999999", "Zinc (women)" = "#666666")

p <- ggplot(d, aes(true_pct, pred_pct)) +
  # perfect-calibration reference
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey45", linewidth = 0.5) +
  # residual drop-lines from each point to the 45° line (the calibration gap)
  geom_segment(aes(xend = true_pct, yend = true_pct, color = outcome_lab),
               linewidth = 0.4, alpha = 0.55) +
  # fitted calibration line across all folds
  geom_abline(slope = slope, intercept = intc, color = "black", linewidth = 0.7) +
  geom_point(aes(color = outcome_lab), size = 3.2) +
  geom_text_repel(aes(label = country_lab, color = outcome_lab),
                  size = 3, show.legend = FALSE, max.overlaps = 20,
                  seed = 1, min.segment.length = 0) +
  annotate("text", x = lim[2] * 0.97, y = lim[2] * 0.985,
           label = "perfect calibration (y = x)", hjust = 1, vjust = 1,
           color = "grey45", size = 3, fontface = "italic") +
  annotate("text", x = lim[2] * 0.55, y = intc + slope * lim[2] * 0.55 - lim[2] * 0.05,
           label = sprintf("fitted slope = %.2f\n(perfect = 1.00)", slope),
           hjust = 0, vjust = 1, size = 3.3, fontface = "bold") +
  scale_color_manual(values = oc_colors, name = "Outcome") +
  scale_x_continuous(limits = lim, expand = c(0, 0)) +
  scale_y_continuous(limits = lim, expand = c(0, 0)) +
  coord_equal() +
  labs(
    title = "Transported models regress toward the training-set mean prevalence",
    subtitle = "Predicted vs. true national prevalence, leave-one-country-out folds (trained on the other three)",
    x = "True national prevalence (survey, %)",
    y = "Predicted national prevalence (%)",
    caption = "Points above the dashed line over-predict; below under-predict. Shallow fit slope = regression to the mean."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey40", hjust = 0),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
out <- "results/figures/transport_calibration_loco.png"
ggsave(out, p, width = 7.2, height = 7.6, dpi = 200, bg = "white")
cat(sprintf("Wrote %s\n", out))
cat(sprintf("Calibration fit: pred = %.2f + %.2f * true  (slope; perfect = 1)\n",
            intc, slope))
