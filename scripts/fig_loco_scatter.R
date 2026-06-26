# =============================================================================
# scripts/fig_loco_scatter.R
#
# Observed vs. predicted admin-2 deficiency prevalence under leave-one-country-
# out (LOCO) cross-validation, faceted by outcome (rows) and held-out country
# (columns). Each point is one admin-2 polygon in the held-out country; the
# prediction comes from the area-level model trained on the other three
# countries. Dashed line is y = x (perfect prediction). A per-panel Pearson r
# (survey vs. modelled prevalence, n_svy-weighted) is annotated in each panel.
#
# Input : results/tables/transportability_area_loco_predictions.csv
# Output: results/figures/loco_scatter.png
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

inp <- "results/tables/transportability_area_loco_predictions.csv"
d <- read.csv(inp, stringsAsFactors = FALSE)
d <- d[is.finite(d$survey_prev) & is.finite(d$modeled_prev), , drop = FALSE]
d$obs_pct  <- d$survey_prev  * 100
d$pred_pct <- d$modeled_prev * 100

# Pretty labels ------------------------------------------------------------
oc_lab <- c(child_vitA = "Vitamin A (children)", women_vitA = "Vitamin A (women)",
            child_iron = "Iron (children)",      women_iron = "Iron (women)",
            women_folate = "Folate (women)",     women_b12 = "Vitamin B12 (women)",
            child_zinc = "Zinc (children)",      women_zinc = "Zinc (women)")
d$outcome_pretty <- ifelse(d$outcome %in% names(oc_lab),
                           unname(oc_lab[d$outcome]), d$outcome)
# CamelCase / lowercase country codes -> display names
ctry_lab <- c(gambia = "Gambia", ghana = "Ghana",
              sierraleone = "Sierra Leone", malawi = "Malawi")
key <- tolower(gsub("[^a-z]", "", tolower(d$country)))
d$country_lab <- ifelse(key %in% names(ctry_lab), unname(ctry_lab[key]), d$country)
d$country_lab <- factor(d$country_lab,
                        levels = c("Gambia", "Ghana", "Sierra Leone", "Malawi"))

# Per-panel Pearson r (unweighted, NA when prediction is constant). Matches the
# definition used in results/tables/transportability_area_loco_metrics.csv so the
# annotated r equals the value reported in the manuscript's @tbl-loco.
upearson <- function(x, y) {
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y))
}
ann <- d |>
  group_by(outcome_pretty, country_lab) |>
  summarise(r = upearson(obs_pct, pred_pct), n = n(),
            x = min(obs_pct), y = max(pred_pct), .groups = "drop")
ann$lab <- ifelse(is.na(ann$r), sprintf("n=%d", ann$n),
                  sprintf("r=%.2f", ann$r))

ctry_colors <- c("Gambia" = "#1b9e77", "Ghana" = "#d95f02",
                 "Sierra Leone" = "#7570b3", "Malawi" = "#e7298a")

p <- ggplot(d, aes(obs_pct, pred_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey45", linewidth = 0.4) +
  geom_point(aes(color = country_lab, size = n_svy), alpha = 0.6,
             show.legend = c(color = FALSE, size = TRUE)) +
  geom_text(data = ann, aes(x = x, y = y, label = lab),
            hjust = 0, vjust = 1, size = 2.8, color = "grey20") +
  facet_grid(outcome_pretty ~ country_lab, scales = "free") +
  scale_color_manual(values = ctry_colors, guide = "none") +
  scale_size_area(max_size = 4, name = "Survey n") +
  labs(
    title = "Cross-country transport of admin-2 deficiency prevalence (LOCO)",
    subtitle = paste("Observed vs. predicted admin-2 prevalence; each point one district in the held-out",
                     "country.\nDashed = perfect prediction (y = x). Panel r = Pearson correlation (matches the LOCO transportability table)."),
    x = "Observed admin-2 prevalence (survey, %)",
    y = "Predicted admin-2 prevalence (LOCO, %)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.6, "lines")
  )

out <- "results/figures/loco_scatter.png"
ggsave(out, p, width = 9.5, height = 9, dpi = 200, bg = "white")
cat(sprintf("Wrote %s  (%d points, %d panels)\n",
            out, nrow(d), nrow(ann)))
print(ann[, c("outcome_pretty", "country_lab", "r", "n")])
