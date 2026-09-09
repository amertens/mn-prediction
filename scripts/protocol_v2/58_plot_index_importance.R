# =============================================================================
# scripts/protocol_v2/58_plot_index_importance.R   [WS-02]
#
# Faceted bar chart of what the zero-tuning index weights: for each outcome the
# ten predictors with the largest |beta| in the pooled four-country fit (script
# 57), bar length = beta x training SD (signed; + = more deficiency), bar
# colour = data source of the predictor. One panel per outcome; one figure per
# target.
#   Rscript scripts/protocol_v2/58_plot_index_importance.R
# -> results/figures/protocol_v2/index_importance_top10_{level,prev}.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
P2 <- "results/tables/protocol_v2"; FIG <- "results/figures/protocol_v2"; dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
IT <- read.csv(file.path(P2, "index_importance_top.csv"), stringsAsFactors = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
src_short <- function(s) dplyr::case_when(
  s == "DHS" ~ "DHS survey aggregates",
  grepl("AlphaEarth", s) ~ "AlphaEarth embedding",
  s == "GEE" ~ "Earth Engine rasters",
  grepl("SoilGrids", s) ~ "SoilGrids / iSDA",
  grepl("IHME", s) ~ "IHME modelled surfaces",
  grepl("Livestock", s) ~ "Livestock (GLW4)",
  grepl("MapSPAM", s) ~ "MapSPAM crops",
  grepl("Malaria Atlas", s) ~ "Malaria Atlas",
  grepl("market prices", s) ~ "WFP market prices",
  grepl("ESPEN", s) ~ "ESPEN helminths",
  grepl("WorldPop", s) | grepl("RWI", s) ~ "WorldPop / GHSL / RWI",
  grepl("Surface Water", s) ~ "Water and coast distance",
  grepl("FAOSTAT", s) | grepl("GFDx", s) | grepl("HFID", s) ~ "National food supply",
  TRUE ~ "Other")
outcome_label <- c(child_vitA = "Child vitamin A", women_vitA = "Women's vitamin A", child_iron = "Child iron",
                   women_iron = "Women's iron", women_folate = "Women's folate", women_b12 = "Women's B12")
pal <- c("DHS survey aggregates" = "#6a3d9a", "Earth Engine rasters" = "#1f78b4", "AlphaEarth embedding" = "#a6cee3",
         "SoilGrids / iSDA" = "#b15928", "IHME modelled surfaces" = "#e31a1c", "Livestock (GLW4)" = "#ff7f00",
         "MapSPAM crops" = "#33a02c", "Malaria Atlas" = "#fb9a99", "WFP market prices" = "#fdbf6f", "ESPEN helminths" = "#cab2d6",
         "WorldPop / GHSL / RWI" = "#b2df8a", "Water and coast distance" = "#08519c", "National food supply" = "#ffff99", "Other" = "grey60")
for (tg in c("level", "prev")) {
  d <- IT |> filter(target == tg, rank <= 10) |> left_join(MD[, c("column", "source")], by = "column") |>
    mutate(source = src_short(source), Outcome = factor(outcome_label[outcome], levels = outcome_label),
           replicated = incountry_sign_agree == incountry_fits,
           short = ifelse(nchar(column) > 34, paste0(substr(column, 1, 32), ".."), column),
           lab = ifelse(replicated, short, paste0(short, " *"))) |>
    group_by(Outcome) |> arrange(Outcome, beta_std) |> mutate(key = factor(paste(lab, Outcome, sep = " | "), levels = paste(lab, Outcome, sep = " | "))) |> ungroup()
  g <- ggplot(d, aes(x = beta_std, y = key, fill = source)) +
    geom_col(width = 0.75) + geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    facet_wrap(~ Outcome, scales = "free_y", ncol = 3) +
    scale_y_discrete(labels = function(k) sub(" [|] .*$", "", k)) +
    scale_fill_manual(values = pal, drop = TRUE, name = "Data source") +
    labs(x = "Index weight per predictor (beta x training SD; + = more deficiency)", y = NULL,
         title = sprintf("What the zero-tuning index weights: ten largest predictors per outcome, %s target",
                         if (tg == "level") "biomarker level" else "prevalence"),
         subtitle = "Pooled four-country fit projected back onto the rank-normalised predictors; * = sign not reproduced in every country's own fit") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"),
                                          axis.text.y = element_text(size = 8), plot.title.position = "plot") +
    guides(fill = guide_legend(nrow = 2))
  ggsave(file.path(FIG, sprintf("index_importance_top10_%s.png", tg)), g, width = 13, height = 8.5, dpi = 170, bg = "white")
  cat("wrote", file.path(FIG, sprintf("index_importance_top10_%s.png", tg)), "\n")
}
