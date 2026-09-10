# =============================================================================
# scripts/policy_deck/03_figure_geostat.R
#
# Figure 9: the proxy index against the DHS-style geostatistical model.
# Reads results/tables/cluster_level/mbg_comparison_cells.csv only.
#
#   Rscript scripts/policy_deck/03_figure_geostat.R
# -> results/figures/policy_deck/fig9_geostatistical.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(grid)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

OUT   <- "results/figures/policy_deck"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROXY <- "#0F7B8A"; GEO <- "#6A51A3"; INK <- "#1A1A1A"

M <- read.csv("results/tables/cluster_level/mbg_comparison_cells.csv", stringsAsFactors = FALSE)

agg <- function(est, tgt, arm, col)
  mean(M[[col]][M$estimand == est & M$target == tgt & M$arm == arm], na.rm = TRUE)

D <- bind_rows(
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "Districts inside\na surveyed country",
             model = "Proxy index",            v = agg("infill", "level", "domain_index", "rho_agg")),
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "Districts inside\na surveyed country",
             model = "Geostatistical model",   v = agg("infill", "level", "mbg", "rho_agg")),
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "A whole region\nheld out",
             model = "Proxy index",            v = agg("region", "level", "domain_index", "rho_agg")),
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "A whole region\nheld out",
             model = "Geostatistical model",   v = agg("region", "level", "mbg", "rho_agg")),
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "Ranking on\nprevalence",
             model = "Proxy index",            v = agg("infill", "prev", "domain_index", "rho_agg")),
  data.frame(panel = "Ranking accuracy (higher is better)",
             measure = "Ranking on\nprevalence",
             model = "Geostatistical model",   v = agg("infill", "prev", "mbg", "rho_agg")),
  data.frame(panel = "Error in the prevalence number (lower is better)",
             measure = "Percentage points\nout, per district",
             model = "Proxy index",            v = agg("infill", "prev", "domain_index", "wmae_agg")),
  data.frame(panel = "Error in the prevalence number (lower is better)",
             measure = "Percentage points\nout, per district",
             model = "Geostatistical model",   v = agg("infill", "prev", "mbg", "wmae_agg")))

D$measure <- factor(D$measure, levels = rev(c("Districts inside
a surveyed country",
                                              "A whole region
held out",
                                              "Ranking on
prevalence",
                                              "Percentage points
out, per district")))
D$model <- factor(D$model, levels = c("Proxy index", "Geostatistical model"))
FILL <- c("Proxy index" = PROXY, "Geostatistical model" = GEO)

base_theme <- function(base = 18) theme_minimal(base_size = base) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 15.5, margin = margin(b = 6)),
        legend.position = "none",
        axis.text.y = element_text(size = 14, lineheight = 0.95),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 13.5, colour = "grey35"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

RK <- D[D$panel == "Ranking accuracy (higher is better)", ]
ER <- D[D$panel != "Ranking accuracy (higher is better)", ]

pL <- ggplot(RK, aes(v, measure, fill = model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", v)), position = position_dodge(width = 0.72),
            hjust = -0.18, size = 5.2, colour = INK) +
  scale_fill_manual(values = FILL) +
  scale_x_continuous(limits = c(0, 0.47), expand = c(0, 0)) +
  labs(title = "Getting the order right", x = "ranking accuracy - higher is better", y = NULL) +
  base_theme()

pR <- ggplot(ER, aes(v, measure, fill = model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", v)), position = position_dodge(width = 0.72),
            hjust = -0.18, size = 5.2, colour = INK) +
  scale_fill_manual(values = FILL) +
  scale_x_continuous(limits = c(0, 14.5), expand = c(0, 0)) +
  labs(title = "Getting the number right", x = "points out - lower is better", y = NULL) +
  base_theme()

leg <- ggplot(RK, aes(v, measure, fill = model)) + geom_col() +
  scale_fill_manual(values = FILL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 16))

p <- patchwork::wrap_plots(pL, pR, widths = c(1.75, 1)) +
  patchwork::plot_annotation(
    subtitle = "It also needs survey clusters, so it cannot run where there is no survey at all.",
    caption = "Same 24 country-outcome combinations, same district folds. The proxy index ranks better in 20 of the 24; on prevalence the two tie.",
    theme = theme(plot.subtitle = element_text(size = 17, colour = "grey20", margin = margin(b = 10)),
                  plot.caption = element_text(size = 12, colour = "grey45", hjust = 0),
                  plot.margin = margin(12, 18, 8, 12))) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 16))

ggsave(file.path(OUT, "fig9_geostatistical.png"), p, width = 12.2, height = 6.3, dpi = 200, bg = "white")
cat("wrote fig9_geostatistical.png
")
