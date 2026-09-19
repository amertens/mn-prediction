# =============================================================================
# scripts/policy_deck/24_civ_2007_survey_check.R
#
# An external check on the transported Cote d'Ivoire ranking with REAL biomarker
# data. The 2007 national nutrition survey reported vitamin B12 deficiency in
# non-pregnant women by nine survey eco-regions (WHO VMNIS; the rows are
# documented in docs/findings/WSG_COTE_DIVOIRE.md): North 48.6% down to
# Abidjan 0%. Our climate-and-soil index, trained on the four surveyed countries
# and never shown a Cote d'Ivoire biomarker, ranks the 33 GADM regions. This
# script assigns each region to a survey zone by centroid, averages the model
# rank within zone, and compares.
#
# Suggested by S. Hess at the 18 Sep 2026 check-in ("isn't Cote d'Ivoire on
# VMNIS? even regional").
#
# CAVEATS, which travel with the figure: the survey is 2007 (19 years old); the
# nine eco-regions are not administrative units, so the zone assignment below
# is an approximate compass crosswalk from region centroids; the survey cut-off
# is unverified; nine points.
#
#   Rscript scripts/policy_deck/24_civ_2007_survey_check.R
# -> results/figures/mnf15/figQ_civ_2007_b12_check.png
#    results/figures/mnf15/civ_b12_vs_2007_survey.csv
# =============================================================================
suppressPackageStartupMessages({library(sf); library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/data_prep.R")
OUT <- "results/figures/mnf15"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROXY <- "#0F7B8A"; WARM <- "#B45309"

# 2007 survey, B12 deficiency in NPW by eco-region (WSG_COTE_DIVOIRE.md)
VM <- c(North = 48.6, `North East` = 40.4, `North West` = 36.4, West = 29.8,
        Central = 20.0, `Central West` = 11.5, `South East` = 7.0, South = 6.2, Abidjan = 0.0)

g <- sf::st_as_sf(load_gadm_cached("CIV", level = 2))
g <- sf::st_make_valid(sf::st_transform(g, 4326))
ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(g))))
c2 <- data.frame(Admin1 = g$NAME_1, Admin2 = g$NAME_2, lon = ctr[, 1], lat = ctr[, 2])

R <- read.csv("results/tables/policy_deck/civ_climate_soil_ranking.csv", stringsAsFactors = FALSE)
b <- R[R$outcome == "women_b12", c("Admin1", "Admin2", "index", "rank")]
m <- dplyr::inner_join(c2, b, by = c("Admin1", "Admin2"))
stopifnot(nrow(m) == 33)

zone <- function(a1, a2, lat, lon) {
  if (grepl("Abidjan", a1) || grepl("Abidjan", a2)) return("Abidjan")
  if (lat >= 8.7) return(if (lon <= -6.6) "North West" else if (lon >= -4.4) "North East" else "North")
  if (lat >= 6.7) return(if (lon <= -6.6) "West" else if (lon <= -5.2) "Central West" else "Central")
  if (lon >= -4.4) "South East" else "South"
}
m$zone <- mapply(zone, m$Admin1, m$Admin2, m$lat, m$lon)

z <- m |> group_by(zone) |>
  summarise(model_rank = mean(rank), n_regions = dplyr::n(), .groups = "drop") |>
  mutate(survey_b12_pct = VM[zone]) |>
  arrange(desc(survey_b12_pct))
rho <- cor(z$survey_b12_pct, z$model_rank, method = "spearman")
cat(sprintf("Spearman(2007 survey B12 %%, model rank) = %.2f over %d zones (negative = agreement)\n", rho, nrow(z)))
print(z)
write.csv(z, file.path(OUT, "civ_b12_vs_2007_survey.csv"), row.names = FALSE)
writeLines(sprintf("%.2f", rho), file.path(OUT, "civ_b12_vs_2007_rho.txt"))

# Figure: survey prevalence on x, model's mean rank on y (1 = worst). Agreement is a
# falling line.
p <- ggplot(z, aes(survey_b12_pct, model_rank)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey75", linewidth = 0.8, formula = y ~ x) +
  geom_point(aes(size = n_regions), colour = PROXY) +
  geom_text(aes(label = zone), vjust = -1.1, size = 4.6) +
  scale_y_reverse(breaks = c(1, 5, 10, 15, 20, 25, 30, 33)) +
  scale_size_continuous(range = c(4, 9), guide = "none") +
  scale_x_continuous(limits = c(-3, 55)) +
  labs(subtitle = sprintf("Nine survey zones. Rank agreement %.2f. The model never saw a Cote d'Ivoire biomarker.", abs(rho)),
       x = "Women's B12 deficiency measured by the 2007 national survey (%)",
       y = "Model's mean district rank in the zone (1 = worst)",
       caption = "2007 survey zones are not administrative units; regions are assigned to zones by centroid. Point size = regions in the zone.") +
  theme_minimal(base_size = 17) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(colour = "grey25", size = 15),
        plot.caption = element_text(colour = "grey45", size = 11, hjust = 0))
ggsave(file.path(OUT, "figQ_civ_2007_b12_check.png"), p, width = 10.5, height = 5.6, dpi = 200, bg = "white")
cat("wrote figQ_civ_2007_b12_check.png\n")
