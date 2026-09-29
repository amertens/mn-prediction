# =============================================================================
# scripts/policy_deck/33_mnf15_v4_vmnis_coverage.R
#
# How few micronutrient biomarker surveys sub-Saharan Africa has: for the v4
# MNF15 introduction slide (Andrew, 27 Sep: the Ghana gap maps do not fit there).
# Source: WHO Vitamin and Mineral Nutrition Information System (VMNIS)
# Micronutrients Database, export of 25 February 2025 (data/RA_2026-09/VMNIS.zip).
# Nationally representative rows only; blood biomarkers of the five nutrients in
# the talk: ferritin (iron), retinol or RBP (vitamin A), zinc, folate (serum or
# red cell), vitamin B12. Survey year = the row's end year.
#
#   Rscript scripts/policy_deck/33_mnf15_v4_vmnis_coverage.R
# -> results/figures/mnf15_v4/v4_vmnis_coverage.png, results/tables/policy_deck/v4_vmnis_coverage.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(readxl); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15_v4"; PD <- "results/tables/policy_deck"
PROXY <- "#0F7B8A"; INK <- "#1A1A1A"; GREY <- "grey40"
tmp <- file.path(tempdir(), "vmnis"); unzip("data/RA_2026-09/VMNIS.zip", exdir = tmp)
FILES <- c(Iron = "Ferritin", `Vitamin A` = "Retinol (plasma or serum)", `Vitamin A ` = "Retinol binding protein",
           Zinc = "Zinc (plasma or serum)", Folate = "Folate (plasma or serum)", `Folate ` = "Folate (red blood cell)", `Vitamin B12` = "Vitamin B12")
V <- bind_rows(lapply(names(FILES), function(k) {
  f <- file.path(tmp, "VMNIS", sprintf("VMNISIndicator_%s_25022025.xlsx", FILES[[k]]))
  x <- read_excel(f, guess_max = 10000)
  data.frame(nutrient = trimws(k), iso3 = x$CountryCode, country = x$Country, rep = x$`Representativeness.`,
             year = suppressWarnings(as.integer(x$`End year`)), survey = x$SurveyId, stringsAsFactors = FALSE)
})) |> filter(rep == "national", is.finite(year))
w <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
ssa <- w[w$region_wb == "Sub-Saharan Africa" & w$continent == "Africa", ]
ssa$iso3 <- ifelse(ssa$iso_a3 == "-99", ssa$adm0_a3, ssa$iso_a3)
S <- V |> filter(iso3 %in% ssa$iso3)
# the four surveys the model learns from must be in the database, or the map undercounts
for (k in c("GHA", "GMB", "SLE", "MWI")) cat(k, "latest national survey in VMNIS:", max(c(S$year[S$iso3 == k], NA), na.rm = TRUE), "\n")
last <- S |> group_by(iso3) |> summarise(last = max(year), n_surveys = n_distinct(survey), .groups = "drop")
ssa <- ssa |> left_join(last, by = "iso3") |>
  mutate(band = factor(case_when(is.na(last) ~ "None on record", last >= 2015 ~ "2015 or later", last >= 2005 ~ "2005 to 2014", TRUE ~ "Before 2005"),
                       levels = c("2015 or later", "2005 to 2014", "Before 2005", "None on record")))
n_ssa <- nrow(ssa)
tab <- table(ssa$band); print(tab)
by_nut <- S |> group_by(nutrient) |> summarise(ever = n_distinct(iso3), recent = n_distinct(iso3[year >= 2015]), .groups = "drop") |>
  mutate(nutrient = factor(nutrient, levels = rev(c("Iron", "Vitamin A", "Folate", "Vitamin B12", "Zinc"))))
print(by_nut)
write.csv(sf::st_drop_geometry(ssa)[, c("iso3", "admin", "last", "n_surveys", "band")], file.path(PD, "v4_vmnis_coverage.csv"), row.names = FALSE)
pal <- c("2015 or later" = "#0B4F5A", "2005 to 2014" = "#5FA8B3", "Before 2005" = "#CFE5E8", "None on record" = "#F2F2F2")
# simple version for the slide (Andrew, 27 Sep): no caption (the source is a footnote on the slide), no subtitle,
# legends at the bottom
pm <- ggplot(ssa) + geom_sf(aes(fill = band), colour = "grey55", linewidth = 0.15) +
  scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
  coord_sf(xlim = c(-18, 52), ylim = c(-35, 25), expand = FALSE) +
  labs(title = "Latest national survey with\nblood biomarkers") +
  theme_void(base_size = 15) +
  theme(plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5, margin = margin(b = 4)),
        legend.position = "bottom", legend.text = element_text(size = 12.5), legend.key.size = unit(0.5, "cm")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
bl <- rbind(data.frame(nutrient = by_nut$nutrient, v = by_nut$ever, what = "Ever"),
            data.frame(nutrient = by_nut$nutrient, v = by_nut$recent, what = "Since 2015"))
bl$what <- factor(bl$what, levels = c("Since 2015", "Ever"))
pb <- ggplot(by_nut, aes(y = nutrient)) +
  geom_col(data = bl[bl$what == "Ever", ], aes(x = v, fill = what), width = 0.66) +
  geom_col(data = bl[bl$what == "Since 2015", ], aes(x = v, fill = what), width = 0.66) +
  geom_text(aes(x = ever, label = ever), hjust = -0.3, size = 5, colour = GREY) +
  geom_text(aes(x = recent, label = recent), hjust = 1.3, size = 5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Since 2015" = "#0B4F5A", "Ever" = "#CFE5E8"), name = NULL, breaks = c("Since 2015", "Ever")) +
  scale_x_continuous(limits = c(0, n_ssa), breaks = c(0, 10, 20, 30, 40, n_ssa)) +
  labs(title = sprintf("Countries with a national survey,\nof %d in sub-Saharan Africa", n_ssa), x = NULL, y = NULL) +
  theme_minimal(base_size = 15) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                        axis.text.y = element_text(size = 14.5, colour = INK), plot.title = element_text(size = 14.5, face = "bold"),
                                        plot.title.position = "plot", legend.position = "bottom", legend.text = element_text(size = 12.5))
ggsave(file.path(OUT, "v4_vmnis_coverage.png"), pm + pb + plot_layout(widths = c(1.2, 1)),
       width = 9.6, height = 5.8, dpi = 220, bg = "white")
cat("wrote v4_vmnis_coverage.png\n")
