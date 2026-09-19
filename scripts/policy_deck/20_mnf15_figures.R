# =============================================================================
# scripts/policy_deck/20_mnf15_figures.R
#
# Figures for the MNF 2026 15-20 minute talk (docs/slides/MNF15-talk-2026-09.qmd).
# Every figure reads a committed result table; nothing is re-estimated here
# except where a figure needs per-district predictions, which are produced by
# 02_figure_ghana_map.R and 21_mnf15_ghana_loco.R.
#
#   Rscript scripts/policy_deck/20_mnf15_figures.R
# -> results/figures/mnf15/*.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2 <- "results/tables/protocol_v2"
PROXY <- "#0F7B8A"; WARM <- "#B45309"; INK <- "#1A1A1A"; GREY <- "grey45"
base <- function(sz = 18) theme_minimal(base_size = sz) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = sz + 3),
        plot.subtitle = element_text(colour = "grey25", size = sz - 1),
        plot.caption = element_text(colour = GREY, size = sz - 6, hjust = 0),
        legend.position = "bottom")
sv <- function(p, f, w = 12.2, h = 5.6) {
  ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 200, bg = "white")
  cat("wrote", f, "\n")
}

CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"))
NCE  <- read.csv(file.path(P2, "nce_targeting_metrics.csv"))
VC   <- read.csv(file.path(P2, "variance_components_ceiling.csv"))
TG   <- read.csv(file.path(P2, "targets_v2.csv"))
DOM  <- read.csv(file.path(P2, "index_importance_domains.csv"))
TC   <- read.csv(file.path(P2, "training_curve_climate_soil.csv"))

arm_sp <- function(est, tgt, a) CELL |> filter(estimand == est, target == tgt, arm == a) |>
  select(country, outcome, sp = spearman)

# The measurability screen: survey-side only (prevalence >= 2%, ceiling >= 0.30)
prev <- TG |> group_by(country, outcome) |>
  summarise(prev = weighted.mean(y_prev, n_eff + 1e-9, na.rm = TRUE), .groups = "drop")
ceil <- VC |> filter(rung == "admin2", target == "prev") |>
  select(country, outcome, ceiling = ceiling_vc, clusters = mean_clusters_per_unit,
         single = share_single_cluster)
M <- prev |> left_join(ceil, by = c("country", "outcome")) |>
  left_join(arm_sp("infill",  "level", "domain_index")   |> rename(infill = sp), by = c("country","outcome")) |>
  left_join(arm_sp("infill",  "level", "region_mean_jk") |> rename(jk = sp),     by = c("country","outcome")) |>
  left_join(arm_sp("country", "level", "domain_index")   |> rename(tr = sp),     by = c("country","outcome")) |>
  mutate(nutrient = recode(sub("^(child|women)_", "", outcome),
                           iron = "Iron", vitA = "Vitamin A", folate = "Folate",
                           b12 = "Vitamin B12", zinc = "Zinc"),
         pop = ifelse(grepl("^child", outcome), "Children", "Women"),
         keep = prev >= 0.02 & !is.na(ceiling) & ceiling >= 0.30)
write.csv(M, file.path(OUT, "cell_master.csv"), row.names = FALSE)

# ---------------------------------------------------------------- FIG A
# Which deficiencies can be mapped at district level at all?
scr <- M |> mutate(reason = case_when(
    prev < 0.02 & (is.na(ceiling) | ceiling < 0.30) ~ "Too rare AND no district signal",
    prev < 0.02                                      ~ "Too rare to rank",
    TRUE                                             ~ "Survey finds no district signal"),
    lab = paste0(country, ", ", tolower(pop), "'s ", tolower(nutrient)))
dropped <- scr |> filter(!keep) |> arrange(reason, prev)
p <- ggplot(dropped, aes(x = reorder(lab, prev), y = pmax(prev, 0.001) * 100, fill = reason)) +
  geom_col(width = .68) + coord_flip() +
  scale_fill_manual(values = c("Too rare to rank" = WARM,
                               "Survey finds no district signal" = PROXY,
                               "Too rare AND no district signal" = "#7A4E8C"), name = NULL) +
  labs(title = "Eight of 24 nutrient-country combinations cannot be mapped by anyone",
       subtitle = "Decided before any model is fitted, from the survey alone",
       y = "National prevalence (%)", x = NULL,
       caption = "Kept if national prevalence is at least 2% (WHO public-health threshold) and the survey's own district estimates carry district-level signal.") +
  base()
sv(p, "figA_measurability_screen.png", 12.2, 5.2)

# ---------------------------------------------------------------- FIG B
# Accuracy by nutrient, inside a country and across a border, against the ceiling
nut <- M |> filter(keep) |>
  mutate(Nutrient = paste0(nutrient, ", ", tolower(pop))) |>
  group_by(Nutrient) |>
  summarise(`Inside a surveyed country` = mean(infill, na.rm = TRUE),
            `Country held out of training` = mean(tr, na.rm = TRUE),
            ceiling = mean(ceiling, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
nl <- nut |> pivot_longer(c(`Inside a surveyed country`, `Country held out of training`),
                          names_to = "what", values_to = "rho") |>
  mutate(what = factor(what, levels = c("Inside a surveyed country", "Country held out of training")))
p <- ggplot(nl, aes(rho, reorder(Nutrient, rho))) +
  geom_col(aes(x = ceiling), fill = "grey90", width = .62) +
  geom_point(aes(colour = what), size = 5) +
  geom_text(aes(label = sprintf("%.2f", rho), colour = what,
                vjust = ifelse(what == "Inside a surveyed country", -1.4, 2.2)),
            size = 4.4, show.legend = FALSE) +
  scale_colour_manual(values = c(PROXY, WARM), name = NULL) +
  scale_x_continuous(limits = c(0, .8)) +
  labs(subtitle = "Grey bar: the best score any model could reach against the survey's own district numbers",
       x = "Ranking accuracy (0 = chance, 1 = perfect)", y = NULL,
       caption = "Measurable combinations only. Folate's collapse across borders is an assay difference, not a model failure.") +
  base()
sv(p, "figB_accuracy_by_nutrient.png", 12.2, 6.0)

# ---------------------------------------------------------------- FIG C
# Accuracy by country, with what the survey could see
cty <- M |> filter(keep) |> group_by(Country = country) |>
  summarise(infill = mean(infill, na.rm = TRUE), tr = mean(tr, na.rm = TRUE),
            ceiling = mean(ceiling, na.rm = TRUE), clusters = mean(clusters, na.rm = TRUE),
            .groups = "drop") |>
  mutate(Country = recode(Country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"),
         lab = sprintf("%s\n%.1f clusters per district", Country, clusters))
cl <- cty |> pivot_longer(c(infill, tr), names_to = "what", values_to = "rho") |>
  mutate(what = recode(what, infill = "Inside the country", tr = "Country held out entirely"),
         what = factor(what, levels = c("Inside the country", "Country held out entirely")))
p <- ggplot(cl, aes(reorder(lab, rho), rho, fill = what)) +
  geom_col(position = position_dodge(.72), width = .66) +
  geom_hline(yintercept = 0.079, linetype = "dashed", colour = GREY) +
  annotate("text", x = 0.7, y = 0.10, label = "chance", colour = GREY, size = 4.2, hjust = 0) +
  scale_fill_manual(values = c(PROXY, WARM), name = NULL) +
  labs(title = "How well the model does follows what each survey could see",
       subtitle = "The Gambia put two clusters in a district; Ghana and Malawi put one",
       x = NULL, y = "Ranking accuracy", caption = "Measurable combinations only. Sierra Leone has too few districts to run the within-country test.") +
  base()
sv(p, "figC_accuracy_by_country.png", 12.2, 5.4)

# ---------------------------------------------------------------- FIG D
# Most important domains of proxy predictors: the top five domains by share of
# the index's weight, per outcome, INCLUDING the three environmental blocks so
# their size relative to everything else is visible. Colour = domain, shared
# across facets. Share is the domain's share of the index weight in the pooled
# four-country fit (biomarker level).
pool <- DOM |> filter(scope == "pooled", target == "level")
top <- pool |> group_by(outcome) |> slice_max(share, n = 5) |> ungroup() |>
  mutate(olab = recode(outcome, child_iron = "Children's iron", child_vitA = "Children's vitamin A",
                       women_b12 = "Women's B12", women_folate = "Women's folate",
                       women_iron = "Women's iron", women_vitA = "Women's vitamin A"),
         dlab = recode(domain,
            "Satellite embedding" = "Satellite imagery",
            "Climate and weather" = "Climate",
            "Soil characteristics" = "Soil",
            "Ecosystem productivity/greenness" = "Plant productivity",
            "Infection and inflammation burden" = "Infectious disease",
            "Education, employment, SES" = "Schooling and employment",
            "Child anthropometry" = "Child wasting and underweight",
            "Infant and young child feeding" = "Infant feeding (meat and fish)",
            "Anaemia and haemoglobin" = "Anaemia (modelled)",
            "Built environment" = "Built environment"))
DPAL <- c("Climate" = "#0F7B8A", "Soil" = "#7C4B24", "Satellite imagery" = "#274C77",
          "Plant productivity" = "#4C9A2A", "Infectious disease" = "#B23A48",
          "Schooling and employment" = "#B45309", "Child wasting and underweight" = "#8A5FA8",
          "Infant feeding (meat and fish)" = "#C2185B", "Anaemia (modelled)" = "#D98E04",
          "Built environment" = "#5E6472")
top$dlab <- factor(top$dlab, levels = names(DPAL))
top <- top |> arrange(olab, share) |> mutate(row = factor(seq_len(dplyr::n())))
p <- ggplot(top, aes(share * 100, row, fill = dlab)) +
  geom_col(width = .72) +
  scale_y_discrete(breaks = top$row, labels = top$dlab) +
  facet_wrap(~ olab, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = DPAL, name = NULL) +
  labs(subtitle = "Share of the model's weight carried by each domain, top five per outcome",
       x = "Share of the model's weight (%)", y = NULL,
       caption = "Pooled four-country fit. These say where deficiency is. They are not causes and not things to change.") +
  base(15) +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10.5)) +
  guides(fill = guide_legend(nrow = 2))
sv(p, "figD_what_each_model_uses.png", 12.6, 6.4)

# Main-body version: women's outcomes only (all four nutrients are measured in women),
# so the slide carries four panels rather than six. Six-panel version stays for the appendix.
topw <- top |> filter(grepl("^Women", olab)) |> arrange(olab, share) |> mutate(row = factor(seq_len(dplyr::n())))
pw <- ggplot(topw, aes(share * 100, row, fill = dlab)) +
  geom_col(width = .72) +
  scale_y_discrete(breaks = topw$row, labels = topw$dlab) +
  facet_wrap(~ olab, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = DPAL, name = NULL) +
  labs(subtitle = "Share of the model's weight carried by each domain, top five per outcome, women",
       x = "Share of the model's weight (%)", y = NULL,
       caption = "Pooled four-country fit. These say where deficiency is. They are not causes and not things to change.") +
  base(16) +
  theme(strip.text = element_text(face = "bold", size = 15), axis.text.y = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 2))
sv(pw, "figD_women_domains.png", 12.6, 6.2)

# ---------------------------------------------------------------- FIG E
# The odds ruler (slide 6, option A)
# The last two ticks are the SAME estimator: the mean over all 18 scored
# nutrient-country combinations, and the mean over the six of those in which it
# scores best. They are not an average over different models.
od <- data.frame(
  what = c("A coin toss", "The survey's own regional averages", "A map of neighbouring districts",
           "Our model (all nutrients, all countries)",
           "Our model, best-predicted nutrients (B12, vitamin A)"),
  v = c(50, 55.5, 59.1, 59.8, 69.9))
od$what <- factor(od$what, levels = rev(od$what))
p <- ggplot(od, aes(v, what)) +
  geom_segment(aes(x = 50, xend = v, yend = what), colour = "grey85", linewidth = 2.6) +
  geom_point(aes(colour = grepl("^Our model", what)), size = 6.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.0f%%", v)), hjust = -0.45, size = 5.4) +
  scale_colour_manual(values = c("grey55", PROXY)) +
  scale_x_continuous(limits = c(48, 76)) +
  labs(title = "How often is the worse district identified correctly?",
       x = "Share of district pairs put in the survey's order", y = NULL,
       caption = "Each district predicted with itself held out. 18 nutrient-country combinations.") +
  base()
sv(p, "figE_odds_ruler.png", 12.2, 5.0)

# ---------------------------------------------------------------- FIG F
# Four paired bars (slide 6, option B)
pb <- data.frame(
  metric = factor(rep(c("Which of two districts is worse",
                        "Worst-third list is correct",
                        "Deficiency in the flagged fifth",
                        "District percentage, points off"), each = 2),
                  levels = c("Which of two districts is worse", "Worst-third list is correct",
                             "Deficiency in the flagged fifth", "District percentage, points off")),
  who = rep(c("The model", "What you use today"), 4),
  v = c(59.8, 55.5, 52.4, 33.0, 32.2, 30.1, 10.7, 11.8))
p <- ggplot(pb, aes(v, who, fill = who)) + geom_col(width = .6) +
  geom_text(aes(label = sprintf("%.1f", v)), hjust = -0.2, size = 4.6) +
  facet_wrap(~ metric, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("The model" = PROXY, "What you use today" = "grey60"), guide = "none") +
  labs(title = "Four ways of asking how good it is", x = NULL, y = NULL,
       caption = "First three: higher is better, per cent. Last: lower is better, percentage points off the survey's district figure (calibrated model).") +
  base(15)
sv(p, "figF_four_metrics.png", 12.2, 5.6)

# ---------------------------------------------------------------- FIG G
# Two-districts pictogram (slide 6, option C)
grid <- expand.grid(x = 1:10, y = c(2, 1))
grid$who <- ifelse(grid$y == 2, "The model", "The survey's regional averages")
grid$lit <- c(rep(TRUE, 6), rep(FALSE, 4), rep(TRUE, 6), rep(FALSE, 4))
p <- ggplot(grid, aes(x, y)) +
  geom_point(aes(colour = lit), size = 11, show.legend = FALSE) +
  geom_text(data = data.frame(x = 0.2, y = c(2, 1), lab = c("The model: 6 of 10", "Regional averages: 5.5 of 10")),
            aes(label = lab), hjust = 1, size = 5, inherit.aes = TRUE) +
  scale_colour_manual(values = c(`TRUE` = PROXY, `FALSE` = "grey88")) +
  scale_x_continuous(limits = c(-4.5, 10.6)) + scale_y_continuous(limits = c(0.5, 2.5)) +
  labs(title = "Out of every ten pairs of districts, how many are put in the right order?",
       caption = "A coin toss gets five. Each district predicted with itself held out.") +
  theme_void(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 21, margin = margin(b = 14)),
        plot.caption = element_text(colour = GREY, size = 12, hjust = 0))
sv(p, "figG_pictogram.png", 12.2, 4.2)

# ---------------------------------------------------------------- FIG H
# Per-country grade strip (slide 6 option B/C bottom)
gs <- M |> filter(keep, !is.na(tr)) |> group_by(country) |>
  summarise(tr = mean(tr), .groups = "drop") |>
  mutate(country = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"))
p <- ggplot(gs, aes(reorder(country, tr), tr)) +
  geom_col(fill = PROXY, width = .6) +
  geom_hline(yintercept = 0.079, linetype = "dashed", colour = GREY) +
  geom_text(aes(label = sprintf("%.2f", tr)), vjust = -0.5, size = 5.4) +
  annotate("text", x = 0.6, y = 0.105, label = "chance", colour = GREY, size = 4, hjust = 0) +
  labs(title = "Each country predicted with its own survey deleted",
       x = NULL, y = "Ranking accuracy", caption = "Trained on the other three countries only.") +
  base()
sv(p, "figH_country_strip.png", 10.5, 4.6)

# ---------------------------------------------------------------- FIG I
# Cost-to-acquire ladder
tiers <- data.frame(
  what = c("Download today, no account", "Add public survey microdata (free account)",
           "Add DHS microdata (data-use agreement)"),
  v = c(0.316, 0.302, 0.281), cells = c("19 of 22", "18 of 22", "17 of 22"))
tiers$what <- factor(tiers$what, levels = rev(tiers$what))
p <- ggplot(tiers, aes(v, what)) + geom_col(fill = PROXY, width = .58) +
  geom_hline(yintercept = 0.079, linetype = "dashed", colour = GREY) +
  geom_text(aes(label = sprintf("%.2f  (%s beat chance)", v, cells)), hjust = -0.06, size = 5) +
  scale_x_continuous(limits = c(0, 0.46)) +
  labs(subtitle = "The data that travel across a border are the free ones",
       x = "Ranking accuracy in a country never surveyed", y = NULL,
       caption = "Dashed line: chance. Survey microdata costs weeks and a registration and makes cross-border prediction slightly worse.") +
  base()
sv(p, "figI_cost_to_acquire.png", 12.2, 4.4)

# ---------------------------------------------------------------- FIG J
# Ceiling staircase
st <- data.frame(clusters = c(1, 2, 3, 5), ceiling = c(0.539, 0.663, 0.731, 0.806))
pins <- data.frame(clusters = c(1.2, 2.33, 1.18, 4.29),
                   ceiling = c(0.538, 0.699, 0.470, 0.622),
                   lab = c("Ghana", "The Gambia", "Malawi", "Sierra Leone"))
p <- ggplot(st, aes(clusters, ceiling)) +
  geom_step(linewidth = 1.5, colour = PROXY) +
  geom_point(data = pins, aes(clusters, ceiling), size = 5, colour = WARM) +
  geom_text(data = pins, aes(label = lab), vjust = -1.1, size = 4.8, colour = WARM) +
  scale_y_continuous(limits = c(0.4, 0.9)) +
  labs(subtitle = "The best score any model could reach against the survey's own district numbers",
       x = "Survey clusters per district", y = "Best attainable ranking accuracy",
       caption = "Projected from these surveys' own variance structure. Two clusters per district raises the bar for the survey's estimates as much as for any model's.") +
  base()
sv(p, "figJ_ceiling_staircase.png", 11.5, 5.4)

# ---------------------------------------------------------------- FIG K
# Pooling curve with "your survey here"
# training_curve_climate_soil.csv carries `set` (full / climate-soil), not `arm`
pc <- TC |> filter(target == "level", set == "full") |>
  group_by(n = n_train_countries) |> summarise(rho = mean(spearman, na.rm = TRUE), .groups = "drop")
nxt <- data.frame(n = max(pc$n) + 1, rho = NA_real_)
p <- ggplot(pc, aes(n, rho)) +
  geom_line(linewidth = 1.4, colour = PROXY) + geom_point(size = 6, colour = PROXY) +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -1.3, size = 5.2) +
  geom_point(data = nxt, aes(n, 0.35), shape = 21, size = 6, colour = GREY, fill = "white", stroke = 1.3) +
  annotate("segment", x = max(pc$n), xend = max(pc$n) + 1, y = max(pc$rho), yend = 0.35,
           linetype = "dashed", colour = GREY) +
  annotate("text", x = max(pc$n) + 1, y = 0.30, label = "your survey here", colour = GREY, size = 5) +
  scale_x_continuous(breaks = 1:4, limits = c(0.85, 4.45)) + scale_y_continuous(limits = c(0.1, 0.42)) +
  labs(x = "Number of biomarker surveys the model learned from",
       y = "Ranking accuracy in a new country",
       caption = "The curve has not flattened at three surveys. Every additional biomarker survey improves the map for every country in the pool.") +
  base()
sv(p, "figK_pooling_curve.png", 11.5, 5.2)

cat("\nall MNF15 figures written to", OUT, "\n")

# ---------------------------------------------------------------- FIG L
# What a country has today: the national number and the regional means, no model.
# Ghana children's iron, from the survey alone. Sonja's opening slide.
suppressPackageStartupMessages({library(sf)})
B <- readRDS("dashboard/data/admin2_boundaries.rds")[["ghana"]]
g <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_prev), is.finite(n_eff))
nat <- weighted.mean(g$y_prev, g$n_eff)
reg <- g |> group_by(Admin1) |> summarise(reg = weighted.mean(y_prev, n_eff), .groups = "drop")
gg <- B |> left_join(reg, by = "Admin1") |>
  left_join(g |> select(Admin1, Admin2, y_prev), by = c("Admin1", "Admin2"))
mk <- function(v, panel) data.frame(sf::st_drop_geometry(gg)[, c("Admin1","Admin2")],
                                    value = v, panel = panel, geometry = sf::st_geometry(gg))
p1 <- "One national number"; p2 <- "Regional averages"; p3 <- "Limited clusters
at the district level"
long <- sf::st_as_sf(rbind(mk(rep(nat, nrow(gg)), p1), mk(gg$reg, p2), mk(gg$y_prev, p3)))
long$panel <- factor(long$panel, levels = c(p1, p2, p3))
p <- ggplot(long) + geom_sf(aes(fill = value * 100), colour = "white", linewidth = 0.15) +
  facet_wrap(~ panel) +
  scale_fill_gradientn(colours = c("#E7EFF0", "#8FC7CF", PROXY, "#0B4F5A"),
                       na.value = "grey92", name = "% of children iron deficient",
                       guide = guide_colourbar(barwidth = 15, barheight = 0.8, ticks = FALSE)) +
  labs(subtitle = "Ghana, children's iron, 2017 survey. Grey: no survey cluster reached the district.") +
  theme_void(base_size = 17) +
  theme(strip.text = element_text(face = "bold", size = 15, margin = margin(b = 7), vjust = 0),
        plot.subtitle = element_text(colour = "grey25", size = 14, margin = margin(b = 8)),
        plot.caption = element_text(colour = GREY, size = 11, hjust = 0),
        legend.position = "bottom")
sv(p, "figL_what_you_have_today.png", 12.6, 5.6)

# ---------------------------------------------------------------- FIG N
# The four training countries on a map of Africa (slide 3).
suppressPackageStartupMessages({library(rnaturalearth)})
af <- rnaturalearth::ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")
TRAIN <- c("Gambia", "Ghana", "Sierra Leone", "Malawi")
af$grp <- ifelse(af$admin %in% TRAIN, "Model training surveys", "Rest of Africa")
lab <- af[af$admin %in% TRAIN, ]
lab$nm <- ifelse(lab$admin == "Gambia", "The Gambia", lab$admin)
ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(lab))))
lab$x <- ctr[, 1]; lab$y <- ctr[, 2]
# nudge the labels off the small West African countries so they do not sit on top of each other
lab$xn <- lab$x + c(Gambia = -13, Ghana = -1, `Sierra Leone` = -12, Malawi = 11)[lab$admin]
lab$yn <- lab$y + c(Gambia = 7, Ghana = -7, `Sierra Leone` = -6, Malawi = 2)[lab$admin]
p <- ggplot(af) +
  geom_sf(aes(fill = grp), colour = "white", linewidth = 0.18) +
  geom_segment(data = lab, aes(x = xn, y = yn, xend = x, yend = y),
               colour = "grey45", linewidth = 0.4) +
  geom_label(data = lab, aes(xn, yn, label = nm), size = 4.6, label.size = 0,
             fill = "white", colour = PROXY, fontface = "bold", label.padding = unit(0.12, "lines")) +
  scale_fill_manual(values = c("Model training surveys" = PROXY, "Rest of Africa" = "grey88"),
                    guide = "none") +
  coord_sf(xlim = c(-40, 54), ylim = c(-36, 38), expand = FALSE) +
  theme_void()
sv(p, "figN_africa_training.png", 7.4, 6.2)

# ---------------------------------------------------------------- FIG O
# Where the predictors come from: every source and every domain, in one figure.
# Lifted from the data-sources deck's "The sources" chart (MN-proxy-data-sources
# -2026-09.qmd) and paired with the same count by domain, which that deck shows
# separately. Colour = whether a column exists in all four countries, i.e.
# whether it can serve a model for a country we have never surveyed.
MDp <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
source_label <- function(s) dplyr::case_when(
  grepl("^DHS", s) ~ "DHS indicators", grepl("^MICS microdata", s) ~ "MICS microdata",
  grepl("HEAT", s) ~ "MICS via WHO HEAT",
  grepl("^HCES", s) ~ "Household budget surveys", grepl("AlphaEarth", s) ~ "AlphaEarth satellite embedding",
  grepl("SoilGrids", s) ~ "SoilGrids / iSDA soil", grepl("^GEE", s) ~ "Earth Engine: climate, land, built",
  grepl("IHME", s) ~ "IHME modelled surfaces", grepl("Malaria Atlas", s) ~ "Malaria Atlas Project",
  grepl("MapSPAM", s) ~ "MapSPAM crops", grepl("Koppen", s) ~ "Koppen / agro-ecological zones",
  grepl("HFID", s) ~ "HFID food security", grepl("Real-Time Food Prices", s) ~ "World Bank RTFP prices",
  grepl("Tang et al", s) ~ "MIMI dietary inadequacy (Ghana)", grepl("WFP", s) ~ "WFP market prices",
  grepl("FAOSTAT", s) ~ "FAOSTAT food balance", grepl("Livestock", s) ~ "Gridded Livestock of the World",
  grepl("ESPEN", s) ~ "WHO ESPEN helminths", grepl("Surface Water", s) ~ "Earth Engine: water and coast",
  grepl("GFDx", s) ~ "GFDx fortification", grepl("Global Anaemia", s) ~ "WHO anaemia estimates",
  grepl("Global Data Lab", s) ~ "Global Data Lab HDI", grepl("VAS", s) ~ "UNICEF vitamin A supplementation",
  grepl("FluNet", s) ~ "WHO FluNet", grepl("ACLED", s) ~ "ACLED conflict events",
  TRUE ~ "WorldPop / GHS / RWI")
cov_lab <- function(n) ifelse(n == 4, "In all four countries", "In one to three countries")
mk_cnt <- function(lab, panel) {
  d <- data.frame(lab = lab, cov = cov_lab(MDp$n_countries), panel = panel) |>
    count(lab, cov, panel, name = "n")
  tot <- tapply(d$n, d$lab, sum)
  d$lab <- factor(d$lab, levels = names(sort(tot)))
  list(d = d, tot = data.frame(lab = factor(names(tot), levels = names(sort(tot))),
                               n = as.vector(tot), panel = panel))
}
a <- mk_cnt(source_label(MDp$source), "By source (28)")
b <- mk_cnt(MDp$domain, "By domain (29)")
dd <- rbind(a$d, b$d); tt <- rbind(a$tot, b$tot)
dd$cov <- factor(dd$cov, levels = c("In all four countries", "In one to three countries"))
PANELS <- c("By source (28)", "By domain (29)")
dd$panel <- factor(dd$panel, levels = PANELS); tt$panel <- factor(tt$panel, levels = PANELS)
p <- ggplot(dd, aes(n, lab, fill = cov)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  geom_text(data = tt, aes(n, lab, label = n), inherit.aes = FALSE, hjust = -0.25, size = 3.1) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(`In all four countries` = PROXY,
                               `In one to three countries` = "#9FCBD3"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Number of district-level predictors", y = NULL,
       caption = sprintf("%d predictors in total. Colour: whether the layer exists in all four countries, and so can serve a model for a country we have never surveyed.", nrow(MDp))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold", size = 13),
        plot.caption = element_text(colour = GREY, size = 9.5, hjust = 0),
        axis.text.y = element_text(size = 9.5))
sv(p, "figO_sources_and_domains.png", 12.6, 6.4)

# ---------------------------------------------------------------- FIG P
# For Sonja's 3-minute landscape talk (asked for at the 18 Sep check-in):
# every data source we drew on, how many district-level predictors it gave us,
# stacked by what the predictors describe. Green / orange / turquoise to match
# her deck. Reads in 20 seconds: most of the volume is environmental and free;
# the survey-derived sources are many but thin.
GROUP <- function(d) dplyr::case_when(
  d %in% c("Satellite embedding", "Climate and weather", "Soil characteristics",
           "Ecosystem productivity/greenness", "Water and coast proximity") ~ "Environment (satellite, climate, soil)",
  d %in% c("Agricultural production, land use", "Livestock density", "Food prices and supply",
           "Market prices (RTFP)", "Household diet and consumption (HCES)",
           "Dietary inadequacy (MODELLED SURFACE, HCES)", "Food fortification and supplementation") ~ "Food system and diet",
  d %in% c("Infection and inflammation burden", "Malaria incidence and treatment", "Helminth burden and control",
           "Infant and child morbidity/mortality", "Child mortality", "Adult and maternal mortality",
           "Immunisation", "Anaemia and haemoglobin", "Child anthropometry", "Adult nutrition",
           "Nutrition status (MODELLED SURFACE)", "Infant and young child feeding") ~ "Health and nutrition status",
  TRUE ~ "Household and socio-economic status")
SONJA <- c("Environment (satellite, climate, soil)" = "#2A9D8F",   # turquoise
           "Food system and diet"                   = "#F4A261",   # orange
           "Health and nutrition status"            = "#7FB069",   # green
           "Household and socio-economic status"     = "#264653")   # dark teal
sp <- data.frame(src = source_label(MDp$source), grp = GROUP(MDp$domain)) |>
  count(src, grp, name = "n")
stot <- tapply(sp$n, sp$src, sum)
sp$src <- factor(sp$src, levels = names(sort(stot)))
sp$grp <- factor(sp$grp, levels = names(SONJA))
p <- ggplot(sp, aes(n, src, fill = grp)) +
  geom_col(width = 0.74) +
  geom_text(data = data.frame(src = factor(names(stot), levels = levels(sp$src)), n = as.vector(stot)),
            aes(n, src, label = n), inherit.aes = FALSE, hjust = -0.25, size = 3.6) +
  scale_fill_manual(values = SONJA, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "District-level predictors drawn from the source", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.caption = element_text(colour = GREY, size = 10, hjust = 0)) +
  guides(fill = guide_legend(nrow = 2))
sv(p, "figP_sources_for_sonja.png", 11.0, 6.6)

# ---------------------------------------------------------------- FIG R
# Simplified "where the predictors come from" for the main talk (Sonja: the
# 28 x 29 version cannot be read quickly). Eight domain groups, with counts and
# two example sources each.
G8 <- function(d) dplyr::case_when(
  d %in% c("Satellite embedding", "Climate and weather", "Ecosystem productivity/greenness", "Water and coast proximity") ~ "Climate, satellite imagery, greenness",
  d == "Soil characteristics" ~ "Soil chemistry",
  d %in% c("Agricultural production, land use", "Livestock density") ~ "Crops and livestock",
  d %in% c("Food prices and supply", "Market prices (RTFP)", "Household diet and consumption (HCES)",
           "Dietary inadequacy (MODELLED SURFACE, HCES)", "Food fortification and supplementation") ~ "Food prices, diet, fortification",
  d %in% c("Infection and inflammation burden", "Malaria incidence and treatment", "Helminth burden and control",
           "Infant and child morbidity/mortality", "Child mortality", "Adult and maternal mortality", "Immunisation") ~ "Infectious disease and immunisation",
  d %in% c("Anaemia and haemoglobin", "Child anthropometry", "Adult nutrition",
           "Nutrition status (MODELLED SURFACE)", "Infant and young child feeding") ~ "Nutrition status and infant feeding",
  d %in% c("Household assets and characteristics", "Education, employment, SES", "Water and sanitation",
           "Fertility, reproductive health", "Healthcare access") ~ "Households, schooling, water, health access",
  TRUE ~ "Settlement, population, conflict")
EX <- c("Climate, satellite imagery, greenness" = "TerraClimate, MODIS, AlphaEarth",
        "Soil chemistry" = "SoilGrids, iSDA",
        "Crops and livestock" = "MapSPAM, Gridded Livestock of the World",
        "Food prices, diet, fortification" = "WFP, World Bank, budget surveys, GFDx",
        "Infectious disease and immunisation" = "Malaria Atlas, WHO ESPEN, IHME",
        "Nutrition status and infant feeding" = "IHME, MICS, DHS",
        "Households, schooling, water, health access" = "MICS, DHS, WorldPop",
        "Settlement, population, conflict" = "WorldPop, GHSL, ACLED")
g8 <- data.frame(grp = G8(MDp$domain), cov = cov_lab(MDp$n_countries)) |> count(grp, cov, name = "n")
gtot <- tapply(g8$n, g8$grp, sum)
g8$grp <- factor(g8$grp, levels = names(sort(gtot)))
g8$cov <- factor(g8$cov, levels = c("In all four countries", "In one to three countries"))
lab8 <- data.frame(grp = factor(names(gtot), levels = levels(g8$grp)), n = as.vector(gtot),
                   ex = EX[names(gtot)])
p <- ggplot(g8, aes(n, grp, fill = cov)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  geom_text(data = lab8, aes(n, grp, label = sprintf("%d   %s", n, ex)), inherit.aes = FALSE,
            hjust = -0.05, size = 4.2, colour = "grey30") +
  scale_fill_manual(values = c(`In all four countries` = PROXY, `In one to three countries` = "#9FCBD3"), name = NULL) +
  scale_x_continuous(limits = c(0, 330), expand = expansion(mult = c(0, 0))) +
  labs(x = "Number of district-level predictors", y = NULL,
       caption = sprintf("%d predictors from 28 public sources, linked to 554 districts. Examples of sources in grey.", nrow(MDp))) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.caption = element_text(colour = GREY, size = 11, hjust = 0),
        axis.text.y = element_text(size = 14))
sv(p, "figR_predictor_groups.png", 12.6, 5.6)
