# =============================================================================
# scripts/policy_deck/05_civ_map.R
#
# Cote d'Ivoire: the deliverable for a country with no biomarker survey.
#
# LEFT   the 33 districts ranked by the climate + soil index, trained on the
#        four surveyed countries and never on Cote d'Ivoire.
# RIGHT  how tightly each district is pinned. Width of the 90% rank interval
#        from 400 bootstrap refits of the training set (script 06), in places
#        out of 33.
#
# The right panel is estimation uncertainty: how much the ranking depends on
# which training districts we happened to learn from. It is NOT a test of
# whether the model transports to Cote d'Ivoire, which no CIV-internal quantity
# can be; the external bound on that is the held-out transport accuracy of 0.37.
#
#   Rscript scripts/policy_deck/06_civ_rank_uncertainty.R   # first
#   Rscript scripts/policy_deck/05_civ_map.R
# -> results/figures/policy_deck/fig10_civ_ranking.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

OUT   <- "results/figures/policy_deck"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROXY <- "#0F7B8A"; INK <- "#1A1A1A"

U <- read.csv("results/tables/policy_deck/civ_rank_uncertainty.csv",
              stringsAsFactors = FALSE, fileEncoding = "UTF-8")
B <- readRDS("dashboard/data/oos_cote_divoire.rds")$boundaries
nD <- nrow(U)

U$pct <- 100 * (U$rank_med - 0.5) / nD
g <- dplyr::left_join(B, U[, c("Admin2", "pct", "rank_med", "rank_width")], by = "Admin2")
cat("districts joined:", sum(is.finite(g$pct)), "of", nrow(B), "\n")

cent <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(g))))
g$X <- cent[, 1]; g$Y <- cent[, 2]
lab <- g[order(g$rank_med), ][1:3, ]

base_map <- function(fill_var, cols, name, breaks, labels, title) {
  ggplot(g) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = cols, na.value = "grey92", name = NULL,
                         breaks = breaks, labels = labels,
                         guide = guide_colourbar(barwidth = 10, barheight = 0.8,
                                                 ticks = FALSE)) +
    labs(title = title) +
    theme_void(base_size = 17) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5,
                                    margin = margin(b = 6)),
          legend.position = "bottom",
          legend.text = element_text(size = 12))
}

pL <- base_map("pct", c("#0B4F5A", PROXY, "#8FC7CF", "#E7EFF0"),
               NULL, c(6, 94), c("worst", "best"),
               "Which districts to reach first") +
  ggrepel::geom_text_repel(data = lab, aes(X, Y, label = Admin2), size = 4.4,
                           fontface = "bold", colour = INK, seed = 1,
                           bg.colour = "white", bg.r = 0.16,
                           box.padding = 0.5, min.segment.length = 0.2,
                           segment.colour = "grey40")

wmax <- round(max(g$rank_width, na.rm = TRUE))
pR <- base_map("rank_width", c("#F2F2F2", "#BFC6CC", "#7C8B95", "#3D4A54"),
               NULL, c(1, wmax), c("pinned down", paste0("could move ", wmax, " places")),
               "How firmly each district is placed")

p <- patchwork::wrap_plots(pL, pR, widths = c(1, 1)) +
  patchwork::plot_annotation(
    subtitle = "C\u00f4te d'Ivoire has no biomarker survey. The index is trained on the other four countries and has never seen it.",
    caption = paste0(
      "Children's iron, from climate and soil layers only. Right: width of the 90% rank range over 400 refits
",
      "on resampled training districts, median 5 places out of 33. The five worst-ranked districts stay in the worst third every time."),
    theme = theme(plot.subtitle = element_text(size = 17, colour = "grey20", margin = margin(b = 12)),
                  plot.caption  = element_text(size = 11.5, colour = "grey45", hjust = 0),
                  plot.margin   = margin(12, 18, 8, 12)))

ggsave(file.path(OUT, "fig10_civ_ranking.png"), p, width = 12.2, height = 6.4, dpi = 200, bg = "white")
cat("wrote fig10_civ_ranking.png\n")
