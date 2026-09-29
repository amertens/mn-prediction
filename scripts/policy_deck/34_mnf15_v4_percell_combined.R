# =============================================================================
# scripts/policy_deck/34_mnf15_v4_percell_combined.R
#
# "How well does it rank districts, nutrient by nutrient?" as ONE panel: for
# every nutrient-country combination, the in-country estimate and the
# whole-country-held-out estimate side by side on the same row, with the
# survey's regional average on the in-country line. Same numbers as the
# two-panel v3_percell_pairs.png of script 26 (read from the same table, no
# refitting); only the layout changes (Andrew, 28 September).
#
#   Rscript scripts/policy_deck/34_mnf15_v4_percell_combined.R
# <- results/tables/policy_deck/v3_percell_pairs.csv   (script 28)
#    results/figures/mnf15/cell_master.csv             (measurability screen)
# -> results/figures/mnf15_v4/v4_percell_pairs_combined.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- Sys.getenv("FIG_OUT", "results/figures/mnf15_v4"); PD <- "results/tables/policy_deck"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROXY <- "#0F7B8A"; WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)

CM <- read.csv("results/figures/mnf15/cell_master.csv")
LAB <- c(in_index = "Inside the country (each district hidden in turn)",
         ho_index = "Whole country held out",
         in_region = "Survey's regional average")
NUDGE <- 0.19   # in-country above the row line, held out below
PP <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |>
  left_join(CM[, c("country", "outcome", "keep", "nutrient", "pop", "prev")], by = c("country", "outcome")) |>
  filter(keep, country != "SierraLeone") |>   # 14 districts: no in-country test (fewer than 12 left to train on), so dropped

  mutate(ctry = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"),
         row = sprintf("%s, %s (%s)", ctry, tolower(pop), ifelse(100 * prev < 10, sprintf("%.1f%%", 100 * prev), sprintf("%.0f%%", 100 * prev))),   # national prevalence in brackets (Andrew, 28 Sep)
         nutrient = factor(nutrient, levels = c("Vitamin B12", "Vitamin A", "Iron", "Folate", "Zinc")),
         who = factor(LAB[ifelse(arm == "region_mean_jk", "in_region", ifelse(estimand == "infill", "in_index", "ho_index"))], levels = LAB),
         dy = ifelse(estimand == "infill", NUDGE, -NUDGE))
# zinc (Malawi only) has an in-country test and no held-out one
chk(sum(PP$who == LAB[["in_index"]]) == 14 && sum(PP$who == LAB[["ho_index"]]) == 14 - sum(PP$nutrient == "Zinc" & PP$estimand == "infill" & PP$arm == "domain_index"),
    "14 in-country combinations, each held out unless measured in one country only")
# FAIR_REGION=1 (v6): the regional average scored as a planner would use it, every district
# in a region given the figure of the region's other surveyed districts, same-region pairs a
# coin toss (script 39). The leave-one-out figure plotted before reverses two surveyed
# districts of the same region by construction and under-rates the survey's regional figures.
if (Sys.getenv("FAIR_REGION") == "1") {
  FR <- read.csv(file.path(PD, "v6_fair_regional_pairs.csv"))
  PP <- PP |> left_join(FR[, c("country", "outcome", "regional_fair")], by = c("country", "outcome")) |>
    mutate(pairs = ifelse(arm == "region_mean_jk", regional_fair, pairs)) |> select(-regional_fair)
  chk(!anyNA(PP$pairs[PP$arm == "region_mean_jk"]), "a fair regional figure for every in-country combination")
}
# rows ordered within nutrient by the mean of the model's two estimates; a row is
# keyed by nutrient AND label ("Ghana, women" recurs under every nutrient)
ordr <- PP |> filter(arm == "domain_index") |> group_by(nutrient, row) |> summarise(m = mean(pairs), .groups = "drop") |>
  arrange(nutrient, m) |> mutate(pos = seq_len(n()))
PP <- PP |> left_join(ordr[, c("nutrient", "row", "pos")], by = c("nutrient", "row"))
PP$y <- PP$pos + PP$dy   # numeric positions so the two estimates sit either side of the row line

p <- ggplot(PP, aes(pairs * 100, y)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = GREY) +
  geom_segment(data = filter(PP, arm == "domain_index"), aes(x = lo * 100, xend = hi * 100, yend = y, colour = who), linewidth = 1.1) +
  geom_point(aes(shape = who, colour = who, fill = who), size = 4, stroke = 1.2) +
  scale_shape_manual(values = setNames(c(21, 24, 23), LAB), name = NULL) +
  scale_colour_manual(values = setNames(c(PROXY, NAVY, WARM), LAB), name = NULL) +
  scale_fill_manual(values = setNames(c(PROXY, NAVY, "white"), LAB), name = NULL) +
  scale_y_continuous(breaks = ordr$pos, labels = ordr$row, expand = expansion(add = 0.45)) +
  facet_grid(nutrient ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(28, 92), breaks = c(30, 50, 70, 90), labels = c("30%", "50%\ncoin toss", "70%", "90%")) +
  labs(x = "Share of district pairs put in the survey's order (95% interval)", y = NULL) +
  theme_minimal(base_size = 15) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 13),
        strip.text.y.left = element_text(face = "bold", size = 13, angle = 0, hjust = 1), strip.placement = "outside",
        axis.text.y = element_text(size = 12.5, colour = INK), panel.spacing.y = unit(0.35, "lines"))
# free_y keeps each nutrient's own rows: breaks outside a facet's range are dropped
ggsave(file.path(OUT, "v4_percell_pairs_combined.png"), p, width = 13, height = 6.4, dpi = 220, bg = "white")
cat("wrote v4_percell_pairs_combined.png\n")
s <- PP |> group_by(who) |> summarise(mean_pairs = round(100 * mean(pairs)), n = n(), .groups = "drop"); print(as.data.frame(s))
