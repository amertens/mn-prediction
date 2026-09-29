# =============================================================================
# scripts/policy_deck/43_mnf15_v6_targeting_and_transport.R
#
# Two figures for the v6 MNF15 talk, from tables already on disk (no fitting):
#
#   v6_targeting_rules.png   (review item B2, with TC-01 / script 74)
#     Share of a country's deficient people reached under two budget rules, by
#     ranking rule, inside a surveyed country and in a country with no survey.
#     Budget of districts (a fifth of them): random, model by rate, the most
#     populous districts (no model), model by expected cases, perfect knowledge.
#     Budget of people (a fifth of the target population): random, the survey's
#     regional figures (in-country only), model by rate, perfect knowledge.
#     Grey dots: each measurable nutrient-country combination; large dot: mean.
#   v6_incountry_vs_heldout.png   (review item B1)
#     District pairs in the survey's order, per measurable combination: inside
#     the country (each district hidden) and with the whole country held out,
#     with the survey's regional figures scored fairly (script 39) as a tick.
#
#   Rscript scripts/policy_deck/43_mnf15_v6_targeting_and_transport.R
# <- results/tables/protocol_v2/tc01_targeting_by_cases_cells.csv (script 74)
#    results/tables/policy_deck/v3_percell_pairs.csv (28), v6_fair_regional_pairs.csv (39)
#    results/figures/mnf15/cell_master.csv
# -> results/figures/mnf15_v6/v6_targeting_rules.png, v6_incountry_vs_heldout.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROXY <- "#0F7B8A"; NAVY <- "#274C77"; WARM <- "#B45309"; INK <- "#1A1A1A"; GREY <- "grey40"
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
th <- theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                                             strip.text = element_text(face = "bold", size = 14), axis.text.y = element_text(colour = INK, size = 13))

# ---------------------------------------------------------------- B2: targeting rules
TC <- read.csv("results/tables/protocol_v2/tc01_targeting_by_cases_cells.csv") |> filter(measurable, is.finite(capture))
# in-country: only the combinations with an in-country test (Sierra Leone has none), for every rule alike
ok_in <- TC |> filter(estimand == "infill", arm == "model_rate") |> transmute(key = paste(country, outcome))
TC <- TC |> filter(estimand == "country" | paste(country, outcome) %in% ok_in$key)
RULE <- c(random = "Random choice", model_rate = "Model, highest rates", region_rate = "Survey's regional figures, highest rates",
          pop_only = "Most populous districts (no model)", model_cases = "Model, most expected cases",
          oracle_cases = "Perfect knowledge", oracle_rate = "Perfect knowledge")
d <- TC |> filter((framing == "districts" & arm %in% c("random", "model_rate", "pop_only", "model_cases", "oracle_cases")) |
                  (framing == "people" & arm %in% c("random", "region_rate", "model_rate", "oracle_rate"))) |>
  mutate(rule = RULE[arm],
         col = case_when(arm %in% c("model_rate", "model_cases") ~ "model", arm == "pop_only" ~ "pop",
                         arm == "region_rate" ~ "region", TRUE ~ "ref"))
n_in <- n_distinct(paste(d$country, d$outcome)[d$estimand == "infill"]); n_lc <- n_distinct(paste(d$country, d$outcome)[d$estimand == "country"])
cat(sprintf("targeting: %d in-country, %d held-out measurable combinations\n", n_in, n_lc))
chk(n_in == 14 && n_lc == 15, "14 in-country and 15 held-out measurable combinations (zinc: one country)")
d <- d |> mutate(est = factor(ifelse(estimand == "infill", sprintf("Inside a surveyed country (%d)", n_in), sprintf("Country with no survey (%d)", n_lc)),
                              levels = c(sprintf("Inside a surveyed country (%d)", n_in), sprintf("Country with no survey (%d)", n_lc))),
                 fr = factor(ifelse(framing == "districts", "Budget: a fifth of the districts", "Budget: a fifth of the people"),
                             levels = c("Budget: a fifth of the districts", "Budget: a fifth of the people")))
ord <- c("Perfect knowledge", "Model, most expected cases", "Most populous districts (no model)", "Model, highest rates",
         "Survey's regional figures, highest rates", "Random choice")
d$rule <- factor(d$rule, levels = rev(ord))
m <- d |> group_by(est, fr, rule, col) |> summarise(v = mean(capture), .groups = "drop")
print(as.data.frame(m |> mutate(v = round(100 * v, 1))), row.names = FALSE)
COL <- c(model = PROXY, pop = NAVY, region = WARM, ref = "grey55")
p2 <- ggplot(d, aes(100 * capture, rule)) +
  geom_vline(xintercept = 20, linetype = "dashed", colour = GREY) +
  geom_point(colour = "grey70", size = 1.8, alpha = 0.7, position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_point(data = m, aes(100 * v, rule, colour = col), size = 5.2) +
  geom_text(data = m, aes(100 * v, rule, label = sprintf("%.0f%%", 100 * v), colour = col), vjust = -0.95, size = 4.6, fontface = "bold") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100), labels = function(x) paste0(x, "%")) +
  facet_grid(fr ~ est, scales = "free_y", space = "free_y") +
  labs(x = "Share of the country's deficient people in the districts chosen", y = NULL) + th +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 13), panel.spacing = unit(1.1, "lines"))
ggsave(file.path(OUT, "v6_targeting_rules.png"), p2, width = 13.5, height = 6.6, dpi = 220, bg = "white")
cat("wrote v6_targeting_rules.png\n")

# ---------------------------------------------------------------- B1: in-country vs held out
CM <- read.csv("results/figures/mnf15/cell_master.csv")
PP <- read.csv("results/tables/policy_deck/v3_percell_pairs.csv") |> filter(arm == "domain_index") |>
  transmute(country, outcome, estimand = recode(estimand, infill = "inside", country = "held"), pairs) |>
  tidyr::pivot_wider(names_from = estimand, values_from = pairs)
FR <- read.csv("results/tables/policy_deck/v6_fair_regional_pairs.csv") |> select(country, outcome, regional = regional_fair)
NUT <- c("Vitamin B12" = "B12", "Vitamin A" = "vitamin A", "Iron" = "iron", "Folate" = "folate", "Zinc" = "zinc")
b <- PP |> inner_join(CM |> filter(keep) |> select(country, outcome, nutrient, pop), by = c("country", "outcome")) |>
  left_join(FR, by = c("country", "outcome")) |>
  mutate(ctry = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"),
         lab = sprintf("%s, %s's %s", ctry, tolower(pop), NUT[nutrient]),
         gap = held - inside)
chk(nrow(b) == 16 && sum(is.finite(b$inside)) == 14 && sum(is.finite(b$held)) == 15 && !anyNA(b$lab), "16 measurable: 14 in-country, 15 held out")
b <- b |> arrange(is.finite(gap), gap) |> mutate(lab = factor(lab, levels = lab))
n_ok <- sum(b$gap >= 0, na.rm = TRUE)
both <- is.finite(b$inside) & is.finite(b$held)
s_in <- mean(b$inside, na.rm = TRUE); s_ho <- mean(b$held[both]); s_rg <- mean(b$regional, na.rm = TRUE)
cat(sprintf("B1: over the 14: inside %.1f, held out %.1f, regional (fair) %.1f; held out as good or better in %d of %d with both; held out over all %d %.1f\n",
            100 * s_in, 100 * s_ho, 100 * s_rg, n_ok, sum(both), sum(is.finite(b$held)), 100 * mean(b$held, na.rm = TRUE)))
LEG <- c(inside = "Inside the country (each district hidden in turn)", held = "Whole country held out", reg = "Survey's regional figures")
p1 <- ggplot(b, aes(y = lab)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = GREY) +
  geom_segment(data = filter(b, is.finite(inside)), aes(x = 100 * inside, xend = 100 * held, yend = lab), colour = "grey70", linewidth = 1.3) +
  geom_point(data = filter(b, is.finite(regional)), aes(x = 100 * regional, shape = LEG[["reg"]], colour = LEG[["reg"]]), size = 5, stroke = 1.6) +
  geom_point(data = filter(b, is.finite(inside)), aes(x = 100 * inside, shape = LEG[["inside"]], colour = LEG[["inside"]]), size = 4.4) +
  geom_point(aes(x = 100 * held, shape = LEG[["held"]], colour = LEG[["held"]]), size = 4.4) +
  scale_shape_manual(values = setNames(c(16, 17, 124), LEG), name = NULL) +
  scale_colour_manual(values = setNames(c(PROXY, NAVY, WARM), LEG), name = NULL) +
  scale_x_continuous(limits = c(40, 85), breaks = c(40, 50, 60, 70, 80), labels = c("40%", "50%\ncoin toss", "60%", "70%", "80%")) +
  labs(x = "Share of district pairs put in the survey's order", y = NULL) + th +
  theme(legend.position = "bottom", legend.text = element_text(size = 13))
ggsave(file.path(OUT, "v6_incountry_vs_heldout.png"), p1, width = 12.5, height = 6.6, dpi = 220, bg = "white")
cat("wrote v6_incountry_vs_heldout.png\n")

# ---------------------------------------------------------------- anchor designs, with the flat national figure (LV-02)
# "Could a small national survey plus the model replace a district survey?" redrawn on the
# post-fix AR-01 table, adding the design AR-01 lacked: every district given the national figure.
AR <- read.csv("results/tables/protocol_v2/anchor_and_rank_summary.csv") |> filter(rank_from == "level", set == "climate_soil")
L0 <- read.csv("results/tables/protocol_v2/lv02_internal_summary.csv") |> filter(rank_from == "level", set == "climate_soil")
chk(isTRUE(all.equal(sort(round(L0$median_mae_A1, 2)), sort(AR$mae[AR$design == "A1_anchor_rank"]), tolerance = 0.011)), "LV-02 A1 matches AR-01 A1")
DES <- c(B_district_survey = "Survey measuring every district directly", C_regional_survey = "Regional survey (one figure per region)",
         A2_region_anchor_rank = "Regional figures + model ranking", A1_anchor_rank = "National figure + model ranking",
         A0 = "National figure alone (no model)")
ad <- bind_rows(AR |> filter(design %in% names(DES)) |> transmute(fraction, design, mae),
                L0 |> transmute(fraction, design = "A0", mae = median_mae_A0)) |>
  mutate(lab = factor(DES[design], levels = DES))
print(as.data.frame(ad |> filter(fraction %in% c(0.05, 0.25, 1)) |> arrange(fraction, design)), row.names = FALSE)
COLD <- setNames(c(WARM, "grey60", "#7FB8C0", PROXY, NAVY), DES)
p3 <- ggplot(ad, aes(100 * fraction, mae, colour = lab, linetype = lab)) +
  geom_line(linewidth = 1.3) + geom_point(size = 2.6) +
  scale_colour_manual(values = COLD, name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "solid", "22"), DES), name = NULL) +
  scale_x_continuous(breaks = c(5, 25, 40, 60, 80, 100), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Survey size, share of a full national survey", y = "Median error in district\nprevalence (points)") +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 11.5),
                                        legend.key.width = unit(1.6, "lines"), legend.margin = margin(0, 0, 0, 0))
# sized for the picture box of the copied slide (5.67 x 4.31 in)
ggsave(file.path(OUT, "v6_anchor_designs.png"), p3, width = 6.6, height = 5.4, dpi = 260, bg = "white")
cat("wrote v6_anchor_designs.png\n")

# ---------------------------------------------------------------- pooled pairs ruler (post-fix; replaces figE_odds_ruler)
SP <- read.csv("results/tables/protocol_v2/sep01_summary.csv") |> filter(cellset == "14 measurable", target == "level")
b12 <- read.csv("results/tables/policy_deck/v3_percell_pairs.csv") |> filter(arm == "domain_index", estimand == "infill", outcome == "women_b12")
chk(nrow(b12) == 2, "two in-country B12 combinations (Ghana, Malawi)")
PR <- data.frame(lab = c("A coin toss", "The survey's own regional figures", "A map of neighbouring districts (from the survey)",
                         "Our model, all nutrients", "Our model, vitamin B12"),
                 v = 100 * c(0.5, SP$pairs_all[SP$arm == "regional_fair"], SP$pairs_all[SP$arm == "spatial"],
                             SP$pairs_all[SP$arm == "domain_index"], mean(b12$pairs)),
                 col = c("ref", "region", "region", "model", "model"))
PR$lab <- factor(PR$lab, levels = rev(PR$lab)); print(PR)
p4 <- ggplot(PR, aes(v, lab, colour = col)) +
  geom_segment(aes(x = 50, xend = v, yend = lab), colour = "grey80", linewidth = 1.4) +
  geom_point(size = 6) + geom_text(aes(label = sprintf("%.0f%%", v)), hjust = -0.55, size = 5.4, fontface = "bold") +
  scale_colour_manual(values = c(ref = "grey55", region = WARM, model = PROXY), guide = "none") +
  scale_x_continuous(limits = c(48, 82), breaks = c(50, 60, 70, 80), labels = function(x) paste0(x, "%")) +
  labs(x = "Share of district pairs put in the survey's order, inside a surveyed country", y = NULL) + th
ggsave(file.path(OUT, "v6_pooled_pairs.png"), p4, width = 12, height = 4.6, dpi = 220, bg = "white")
cat("wrote v6_pooled_pairs.png\n")
