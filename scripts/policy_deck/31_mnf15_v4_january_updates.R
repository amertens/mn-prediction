# =============================================================================
# scripts/policy_deck/31_mnf15_v4_january_updates.R
#
# Updated versions of the best slides of the January 2026 Ghana presentation
# (docs/slides/MN-proxy-january-Ghana-presentation-best slides.pptx, Andrew's
# file with his comments), rebuilt on the current protocol-v2 results for the v4
# MNF15 talk (docs/slides/MNF15-talk-2026-09-v4.qmd). Reads committed tables only.
#
#   v4_cv_schematic.png     January slide 6 comment: show the cross-validation process,
#                           with held-out districts, a region and a country on a map
#   v4_skill_by_outcome.png January slide 7: percent improvement over the national average,
#                           now district prevalence error, by nutrient and information source
#   v4_national.png         January slide 9: can proxy models recover NATIONAL prevalence?
#                           now without the country's own survey (national indicators, VMNIS)
#   v4_sources.png          January slide 11: which data sources matter, now as a whole
#                           country held out (drop one source / that source alone)
#   v4_surveys.png          January slide 4: the four surveys and years, with the survey-team credit
#   printed numbers         January slide 8: the policy metrics for the best model (women's B12)
#
#   Rscript scripts/policy_deck/31_mnf15_v4_january_updates.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R")
OUT <- Sys.getenv("FIG_OUT", "results/figures/mnf15_v4"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
V5 <- grepl("mnf15_v5|mnf15_v6", OUT)   # v5 (28 Sep, Andrew's comments on v4-ANM): plainer labels, no captions baked into images
P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"
PROXY <- "#0F7B8A"; WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"; LGREY <- "grey62"
sv <- function(p, f, w, h) { ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 220, bg = "white"); cat("wrote", f, "\n") }
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
th <- function(base = 16) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(size = 12, colour = GREY, hjust = 0),
        plot.subtitle = element_text(size = 14, colour = GREY), plot.title.position = "plot", plot.caption.position = "plot")
CM <- read.csv("results/figures/mnf15/cell_master.csv")
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"))

# ============================================================================ cross-validation schematic
sf::sf_use_s2(FALSE)
BND <- readRDS("dashboard/data/admin2_boundaries.rds"); B <- BND[["ghana"]]
g <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_level), is.finite(n_eff_cont))
chk(nrow(g) == 75, "75 surveyed Ghana districts")
folds <- make_folds_v2("kfold_district", nrow(g), k = 5, rep_id = 1)   # the benchmark's first fold draw
g$fold <- folds
REG <- "Northern"
chk(REG %in% g$Admin1, "held-out region present")
COLS <- c(fit = PROXY, hidden = WARM, none = "grey82")
LAB <- c(fit = "Used to fit the model", hidden = if (V5) "Held out test data" else "Hidden, then predicted and scored", none = "No survey data")
mapfill <- function(status_of) {
  b <- B |> left_join(g |> select(Admin1, Admin2, fold), by = c("Admin1", "Admin2"))
  b$st <- status_of(b); b$st <- factor(b$st, levels = names(COLS)); b
}
m1 <- mapfill(function(b) ifelse(is.na(b$fold), "none", ifelse(b$fold == 1, "hidden", "fit")))
m2 <- mapfill(function(b) ifelse(is.na(b$fold), "none", ifelse(b$Admin1 == REG, "hidden", "fit")))
mt <- theme_void(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 14.5, hjust = 0.5, margin = margin(b = 4)),
                                        plot.subtitle = element_text(size = 12, hjust = 0.5, colour = GREY), legend.position = "none")
pm1 <- ggplot(m1) + geom_sf(aes(fill = st), colour = "white", linewidth = 0.1) + scale_fill_manual(values = COLS) +
  labs(title = "One district in five hidden", subtitle = if (V5) NULL else "Ghana, one of five parts") + mt
pm2 <- ggplot(m2) + geom_sf(aes(fill = st), colour = "white", linewidth = 0.1) + scale_fill_manual(values = COLS) +
  labs(title = "A whole region hidden", subtitle = if (V5) NULL else sprintf("Ghana, %s Region", REG)) + mt
w <- rnaturalearth::ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
w$st <- factor(dplyr::case_when(w$admin == "Ghana" ~ "hidden", w$admin %in% c("Gambia", "Sierra Leone", "Malawi") ~ "fit", TRUE ~ "none"), levels = names(COLS))
pm3 <- ggplot(w) + geom_sf(aes(fill = st), colour = "white", linewidth = 0.15) + scale_fill_manual(values = COLS) +
  coord_sf(xlim = c(-18, 37), ylim = c(-18, 17), expand = FALSE) +
  labs(title = "A whole country hidden", subtitle = if (V5) NULL else "Ghana, learned from the other three") + mt
fg <- expand.grid(part = 1:5, round = 1:5) |> mutate(st = factor(ifelse(part == round, "hidden", "fit"), levels = names(COLS)))
fg <- rbind(fg, data.frame(part = 99, round = 1, st = factor("none", levels = names(COLS))))   # off-canvas, so the legend draws a grey "No survey data" key
pfg <- ggplot(fg, aes(part, -round, fill = st)) + geom_tile(colour = "white", linewidth = 1.4, width = 0.95, height = 0.8) +
  scale_fill_manual(values = COLS, labels = LAB, drop = FALSE, name = NULL) +
  annotate("text", x = 0.3, y = -(1:5), label = paste("Round", 1:5), hjust = 1, size = 3.9, colour = GREY) +
  coord_cartesian(xlim = c(-0.6, 5.5), clip = "off") +
  labs(title = "Cross-validation", subtitle = if (V5) NULL else "Surveyed districts split in five parts;\nrepeated with ten random splits") + mt +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 12.5)) +
  guides(fill = guide_legend(override.aes = list(fill = unname(COLS), colour = "grey55", linewidth = 0.4)))   # every key drawn, the grey one outlined
sv(pfg + pm1 + pm2 + pm3 + plot_layout(widths = c(1.05, 1, 1, 1.25)), "v4_cv_schematic.png", 13, 5.0)

# ============================================================================ skill over the national average, by outcome
# District PREVALENCE error, each district hidden in turn (in-fill, prevalence target), as the January slide's
# "percent improvement over the null model" but at district level: 1 - MAE(method) / MAE(national average).
ARM <- c(domain_index_cal = "Public data only (Domain-PC index)", spatial = "Neighbour map (needs the survey)",
         region_mean_jk = "Regional averages (needs the survey)")
sk <- CELL |> filter(estimand == "infill", target == "prev", arm %in% c(names(ARM), "null_train_mean")) |>
  inner_join(CM |> filter(keep) |> select(country, outcome, nutrient, pop), by = c("country", "outcome")) |>
  select(country, outcome, nutrient, pop, arm, mae) |> pivot_wider(names_from = arm, values_from = mae) |>
  filter(if_all(all_of(c(names(ARM), "null_train_mean")), is.finite)) |>   # Sierra Leone is not scored in-country
  pivot_longer(all_of(names(ARM)), names_to = "arm", values_to = "mae") |>
  mutate(skill = 100 * (1 - mae / null_train_mean), who = factor(ARM[arm], levels = rev(ARM)),
         nutrient = factor(nutrient, levels = c("Vitamin A", "Iron", "Folate", "Vitamin B12", "Zinc")),
         panel = factor(paste0(nutrient, ", ", tolower(pop)), levels = c("Vitamin A, children", "Vitamin A, women", "Iron, children", "Iron, women", "Folate, women", "Vitamin B12, women", "Zinc, women")),
         row = factor(c(domain_index_cal = "Public data only", spatial = "Neighbour map", region_mean_jk = "Regional averages")[arm],
                      levels = rev(c("Public data only", "Neighbour map", "Regional averages"))))
cat("  skill cells:", n_distinct(paste(sk$country, sk$outcome)), "
")
sm <- sk |> group_by(panel, row, who, arm) |> summarise(skill = mean(skill), n = dplyr::n(), .groups = "drop")
cat("  skill over the national average (mean over measurable cells):\n"); print(as.data.frame(sk |> group_by(arm) |> summarise(skill = mean(skill))), digits = 3)
psk <- ggplot(sk, aes(skill, row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#C0392B") +
  geom_point(aes(colour = who), size = 2.4, alpha = 0.45) +
  geom_point(data = sm, aes(colour = who), size = 5.2) +
  facet_wrap(~ panel, ncol = 3) +
  scale_colour_manual(values = setNames(c(PROXY, NAVY, WARM), unname(ARM)), name = NULL, breaks = unname(ARM)) +
  labs(x = "Percent improvement in district prevalence error over the national average", y = NULL,
       caption = if (V5) NULL else "Large dot: average over countries; small dots: each country. Each district hidden in turn; measurable combinations only.") +
  th(15) + theme(strip.text = element_text(face = "bold", size = 14.5), legend.position = "bottom", legend.text = element_text(size = 13.5),
                 axis.text.y = element_text(size = 13, colour = INK), panel.spacing = unit(1.2, "lines"))
sv(psk, "v4_skill_by_outcome.png", 13, 5.8)

# ============================================================================ national prevalence without the country's survey
NL <- read.csv("results/tables/national_composition_levels.csv") |>
  mutate(ctry = recode(country, gambia = "The Gambia", ghana = "Ghana", sierraleone = "Sierra Leone", malawi = "Malawi"),
         pop = ifelse(grepl("^child", outcome), "children", "women"), lab = paste0(ctry, ", ", pop))
chk(nrow(NL) == 8, "eight vitamin A national cells")
nl <- NL |> select(lab, pop, Survey = true_national_pp, `Predicted from national indicators` = vmnis_level_pp) |>
  pivot_longer(-c(lab, pop), names_to = "what", values_to = "v")
ordr <- NL |> arrange(pop, true_national_pp) |> pull(lab)
nl$lab <- factor(nl$lab, levels = ordr)
pn1 <- ggplot(nl, aes(v, lab)) +
  geom_line(aes(group = lab), colour = "grey70", linewidth = 1) +
  geom_point(aes(colour = what), size = 5) +
  geom_text(data = NL |> mutate(lab = factor(lab, levels = ordr)), aes(x = pmax(true_national_pp, vmnis_level_pp), y = lab, label = sprintf("%+.0f", vmnis_level_pp - true_national_pp)),
            hjust = -0.6, size = 4.6, colour = GREY) +
  scale_colour_manual(values = c(Survey = WARM, `Predicted from national indicators` = NAVY), name = NULL) +
  scale_x_continuous(limits = c(0, 46)) +
  labs(title = "Vitamin A deficiency, each country predicted\nwithout its own survey", x = "National prevalence (%)", y = NULL) +
  th(15) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 15),
                 axis.text.y = element_text(size = 13.5, colour = INK))
VL <- read.csv("results/tables/national_vmnis_loco.csv") |> mutate(model = ifelse(is.na(model) | model == "", "null", model)) |>
  filter(model %in% c("null", "ridge", "rf")) |>
  mutate(panel = sub(" \\| Preschool-age children", ", children", sub(" \\| Non-pregnant women \\(NPW\\)", ", women", panel)),
         panel = sprintf("%s\n(%d countries)", panel, n_countries),
         who = factor(c(null = "Average of the other countries", ridge = "Model on national indicators (ridge)", rf = "Model on national indicators (forest)")[model],
                      levels = c("Average of the other countries", "Model on national indicators (ridge)", "Model on national indicators (forest)")))
pn2 <- ggplot(VL, aes(mae_pp, panel, fill = who)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f", mae_pp)), position = position_dodge(width = 0.8), hjust = -0.3, size = 4.2, colour = GREY) +
  scale_fill_manual(values = c("grey70", "#7FA7C9", NAVY), name = NULL) + scale_x_continuous(limits = c(0, 22)) +
  labs(title = "WHO national survey data, each country\nleft out in turn", x = "Typical error in national prevalence (percentage points)", y = NULL) +
  th(15) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 15),
                 axis.text.y = element_text(size = 12.5, colour = INK), panel.grid.major.y = element_blank())
sv(pn1 + pn2 + plot_layout(widths = c(1, 1.15)) +
     plot_annotation(caption = "Levels do not transport: that is why the model's ranking is anchored to a measured national figure (a small national sample, about 40 blood draws per nutrient).",
                     theme = theme(plot.caption = element_text(size = 13, colour = INK, hjust = 0))),
   "v4_national.png", 13, 5.8)
cat("  national vitamin A misses (pp):\n"); print(NL[, c("lab", "true_national_pp", "vmnis_level_pp", "vmnis_err_pp")])

# ============================================================================ data sources, a whole country held out
SRC <- c("SoilGrids/iSDA" = "Soil maps (SoilGrids, iSDA)", GEE = "Climate and greenness (Earth Engine)", MapSPAM = "Crop mix (MapSPAM)",
         IHME = "Modelled health surfaces (IHME)", "Malaria Atlas" = "Malaria Atlas Project", WorldPop.GHS = "Population and built-up land (WorldPop, GHSL)",
         "Koppen/AEZ" = "Climate zones (Koppen, agro-ecological)", AlphaEarth = "Satellite imagery summary (AlphaEarth)", WFP.prices = "Food prices (WFP)",
         SHDI = "Subnational human development (Global Data Lab)", ESPEN = "Worm surveys (WHO ESPEN)", ACLED = "Conflict events (ACLED)",
         HCES = "Household budget surveys (HCES)", MICS = "MICS survey microdata", HEAT = "MICS regional indicators (WHO HEAT)",
         GLW = "Livestock density (FAO GLW4)", JRC = "Surface water and coast (JRC)")
key <- function(s) dplyr::case_when(grepl("^Gridded Livestock", s) ~ "GLW", grepl("^JRC", s) ~ "JRC", grepl("HEAT", s) ~ "HEAT",
                                    grepl("SHDI", s) ~ "SHDI", grepl("ESPEN", s) ~ "ESPEN", grepl("ACLED", s) ~ "ACLED",
                                    grepl("^HCES", s) ~ "HCES", grepl("^MICS", s) ~ "MICS", s == "WorldPop/GHS" ~ "WorldPop.GHS",
                                    s == "WFP prices" ~ "WFP.prices", TRUE ~ s)
SA <- read.csv(file.path(P2, "source_ablation_loco_summary.csv")) |> filter(target == "level") |>
  mutate(k = key(source), lab = SRC[k])
chk(!anyNA(SA$lab), "every source has a plain name")
full <- unique(SA$full); chk(length(full) == 1, "one full-model score")
SA <- SA |> arrange(only) |> mutate(lab = factor(lab, levels = lab))
ps1 <- ggplot(SA, aes(only, lab)) + geom_col(aes(fill = only > full), width = 0.7) +
  geom_vline(xintercept = full, linetype = "dashed", colour = INK) +
  annotate("text", x = full - 0.006, y = 2.6, label = sprintf("all sources\ntogether %.2f", full), hjust = 1, size = 3.9, colour = INK, lineheight = 0.9) +
  scale_fill_manual(values = c(`TRUE` = PROXY, `FALSE` = "#9DB8BD"), guide = "none") +
  labs(title = "That source alone", x = "Agreement in a country held out (correlation)", y = NULL) +
  th(13.5) + theme(plot.title = element_text(face = "bold", size = 14), axis.text.y = element_text(size = 12, colour = INK))
ps2 <- ggplot(SA, aes(delta_drop, lab)) + geom_col(aes(fill = delta_drop > 0), width = 0.7) + geom_vline(xintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(`TRUE` = WARM, `FALSE` = "grey75"), guide = "none") +
  labs(title = "Lost when it is removed", x = "Drop in agreement", y = NULL) +
  th(13.5) + theme(plot.title = element_text(face = "bold", size = 14), axis.text.y = element_blank())
sv(ps1 + ps2 + plot_layout(widths = c(1.25, 1)) +
     plot_annotation(caption = "Average status, 22 nutrient-country combinations, each country predicted from the other three (the Domain-PC index). Survey microdata and household budget surveys\nmake cross-border prediction slightly worse; soil and climate carry it.",
                     theme = theme(plot.caption = element_text(size = 11.5, colour = GREY, hjust = 0))),
   "v4_sources.png", 13, 5.6)

# ============================================================================ policy metrics for the best model (women's B12)
PP <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(outcome == "women_b12")
cat("\n  women's B12, pairs:\n"); print(PP[, c("country", "estimand", "arm", "n", "spearman", "pairs", "lo", "hi")], digits = 3)
NT <- read.csv(file.path(P2, "nce_targeting_metrics.csv")) |> filter(estimand == "infill", outcome == "women_b12") |>
  group_by(country, arm) |> summarise(capture = mean(capture_top20, na.rm = TRUE), lift = mean(lift, na.rm = TRUE), .groups = "drop")
cat("  women's B12, worst-fifth capture:\n"); print(as.data.frame(NT), digits = 3)
CP <- read.csv(file.path(P2, "conformal_prev_cells.csv")) |> filter(outcome == "women_b12")
cat("  women's B12, conformal bands:\n"); print(CP, digits = 3)
WF <- read.csv(file.path(PD, "worst_fifth_probability.csv")) |> filter(outcome == "women_b12")
top <- WF |> group_by(country) |> slice_max(p_worst_fifth, n = 5, with_ties = FALSE) |> summarise(hit_top5 = mean(obs_worst_fifth), .groups = "drop")
cat("  women's B12, share of the model's five most likely worst-fifth districts that are in the survey's worst fifth:\n"); print(as.data.frame(top))
XT <- read.csv("results/tables/external_validation/xv_transport.csv") |> filter(arm == "domain_index", grepl("b12", outcome), is.finite(spearman))
cat("  women's B12, external:\n"); print(XT[, intersect(c("country", "outcome", "target", "soil", "arm_group", "n_units", "spearman"), names(XT))], digits = 3)
# ============================================================================ the four surveys (January slide 4: list the surveys and years)
SY <- read.csv("metadata/survey_years.csv", stringsAsFactors = FALSE) |> filter(in_protocol)
nd <- TG |> filter(outcome %in% c("child_iron", "women_iron"), is.finite(y_level)) |> group_by(country) |> summarise(n = n_distinct(paste(Admin1, Admin2)), .groups = "drop")
chk(all(c(30, 75, 87, 14) %in% nd$n), "district counts 30 / 75 / 87 / 14")
ST <- data.frame(country = c("Gambia", "Ghana", "SierraLeone", "Malawi"),
                 name = c("The Gambia", "Ghana", "Sierra Leone", "Malawi"),
                 districts = c("30 districts", "75 of 260 districts", "14 districts", "87 areas
(Traditional Authorities)"),
                 nutrients = c("Iron, vitamin A", "Iron, vitamin A,
folate, B12", "Iron, vitamin A,
folate, B12", "Iron, vitamin A,
folate, B12, zinc")) |>
  left_join(SY |> select(country, survey, fieldwork_start), by = "country") |>
  mutate(survey = sub(" (", "\n(", survey, fixed = TRUE))
chk(nrow(ST) == 4 && !anyNA(ST$survey), "four surveys named")
if (V5) {   # "n of N" for every country, from the boundary files the maps use
  n_all <- c(Gambia = nrow(BND[["gambia"]]), Ghana = nrow(BND[["ghana"]]), SierraLeone = nrow(BND[["sierraleone"]]), Malawi = nrow(BND[["malawi"]]))
  n_srv <- setNames(nd$n, nd$country)[ST$country]
  ST$districts <- sprintf("%d of %d %s", n_srv, n_all[ST$country], ifelse(ST$country == "Malawi", "areas", "districts"))
  cat("  surveys, n of N:", paste(ST$name, ST$districts, collapse = "; "), "\n")
}
w2 <- rnaturalearth::ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
w2$st <- ifelse(w2$admin %in% c("Gambia", "Ghana", "Sierra Leone", "Malawi"), "yes", "no")
lb <- w2[w2$st == "yes", ]; lc <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(lb)))); lb$X <- lc[, 1]; lb$Y <- lc[, 2]
lb$nm <- ifelse(lb$admin == "Gambia", "The Gambia", lb$admin)
pmap <- ggplot(w2) + geom_sf(aes(fill = st), colour = "white", linewidth = 0.15) +
  ggrepel::geom_text_repel(data = sf::st_drop_geometry(lb), aes(X, Y, label = nm), size = 5, fontface = "bold", colour = PROXY,
                           seed = 4, bg.colour = "white", bg.r = 0.15, box.padding = 0.6, min.segment.length = 0) +
  scale_fill_manual(values = c(yes = PROXY, no = "grey88"), guide = "none") + coord_sf(xlim = c(-19, 42), ylim = c(-30, 24), expand = FALSE) + theme_void()
tb <- ST |> mutate(row = 4:1) |> select(row, name, survey, districts, nutrients) |> pivot_longer(-row, names_to = "col", values_to = "txt") |>
  mutate(x = c(name = 0, survey = 1.5, districts = 5.0, nutrients = 7.4)[col], face = ifelse(col == "name", "bold", "plain"))
hd <- data.frame(x = c(0, 1.5, 5.0, 7.4), txt = c("Country", "Survey", "Where blood was drawn", "Nutrients measured"))
ptb <- ggplot(tb) + geom_text(aes(x, row, label = txt, fontface = face), hjust = 0, size = 4.6, lineheight = 0.9, colour = INK) +
  geom_text(data = hd, aes(x, 4.75, label = txt), hjust = 0, size = 5, fontface = "bold", colour = GREY) +
  geom_hline(yintercept = c(1.5, 2.5, 3.5), colour = "grey88") +
  coord_cartesian(xlim = c(0, 9.6), ylim = c(0.5, 5), expand = FALSE) + theme_void()
sv(pmap + ptb + plot_layout(widths = c(1, 2.6)) +
     plot_annotation(caption = if (V5) NULL else "206 districts with biomarker data in all. With thanks to the national survey teams and their partners, including, for Ghana's 2017 survey, the University of Ghana,\nGroundWork, the University of Wisconsin-Madison and KEMRI-Wellcome Trust, with UNICEF and Global Affairs Canada.",
                     theme = theme(plot.caption = element_text(size = 12, colour = GREY, hjust = 0))),
   "v4_surveys.png", 13, 4.6)
# ============================================================================ full-talk updates
# (docs/slides/MN-proxy-full-talk-best slides.pptx, Andrew's comments of 27 September)
# FT slide 18: prevalence error by outcome, with a legend and points an audience can read
OL <- c(child_vitA = "children's vitamin A", women_vitA = "women's vitamin A", child_iron = "children's iron", women_iron = "women's iron",
        women_folate = "women's folate", women_b12 = "women's B12", child_zinc = "children's zinc", women_zinc = "women's zinc")
cl <- function(cn, on) paste0(recode(cn, SierraLeone = "Sierra Leone", Gambia = "The Gambia"), ", ", OL[on])
ordc <- CELL |> filter(estimand == "infill", arm == "domain_index", target == "level", is.finite(spearman)) |> arrange(spearman) |> mutate(cell = cl(country, outcome))
ER <- c(domain_index_cal = "Model (calibrated for level)", region_mean_jk = "Survey's regional average",
        null_train_mean = if (V5) "National estimate" else "One national number for every district")
de <- CELL |> filter(target == "prev", estimand == "infill", arm %in% names(ER), is.finite(wmae)) |>
  mutate(cell = factor(cl(country, outcome), levels = ordc$cell), who = factor(ER[arm], levels = ER)) |> filter(!is.na(cell))
rng <- de |> group_by(cell) |> summarise(lo = min(wmae), hi = max(wmae), .groups = "drop")
strong <- ordc$cell[ordc$spearman >= 0.5]
pe <- ggplot(de, aes(wmae, cell)) +
  geom_segment(data = rng, aes(x = lo, xend = hi, y = cell, yend = cell), colour = "grey85", linewidth = 2.2, inherit.aes = FALSE) +
  geom_point(aes(shape = who, colour = who, size = who), stroke = 1.3) +
  scale_shape_manual(values = setNames(c(16, 18, if (V5) 4 else 124), ER), name = NULL) +
  scale_colour_manual(values = setNames(c(PROXY, WARM, "grey35"), ER), name = NULL) +
  scale_size_manual(values = setNames(c(4.8, 5.6, if (V5) 4.2 else 6.5), ER), name = NULL) +
  labs(x = "Typical error in a district's predicted prevalence (percentage points; each district hidden in turn)", y = NULL,
       caption = if (V5) NULL else "Rows ordered by how well the model ranks each combination; bold: ranking of 0.5 or more. Population-weighted mean absolute error.") +
  th(15) + theme(legend.position = "top", legend.text = element_text(size = 14), panel.grid.major.y = element_blank(),
                 axis.text.y = element_text(size = 12.5, colour = if (V5) INK else ifelse(levels(de$cell) %in% strong, PROXY, "grey30"),
                                            face = if (V5) "plain" else ifelse(levels(de$cell) %in% strong, "bold", "plain")))
sv(pe, "v4_ft_prev_error.png", 13, 6.0)
cat(sprintf("  prevalence error means: %s\n", paste(sprintf("%s %.1f", levels(de$who), tapply(de$wmae, de$who, mean)), collapse = "; ")))

# FT slide 22: the four-panel map on a stronger combination, Malawi women's B12 (tables: 10_viz_tables.R A,B2 with overrides)
VZ <- "results/tables/policy_deck/viz"
oof <- read.csv(file.path(VZ, "oof_women_b12.csv")) |> filter(country == "Malawi")
dep <- read.csv(file.path(VZ, "deploy_malawi_women_b12_cal.csv"))
Bm <- sf::st_as_sf(BND[["malawi"]])
gm <- Bm |> left_join(oof[, c("Admin1", "Admin2", "y_prev", "pred")], by = c("Admin1", "Admin2")) |>
  left_join(dep[, c("Admin1", "Admin2", "prev_anchored", "p_modplus_cal", "th_modplus")], by = c("Admin1", "Admin2"))
chk(sum(is.finite(gm$y_prev)) >= 80 && sum(is.finite(gm$prev_anchored)) >= 200, "Malawi B12 tables joined")
rho_m <- cor(gm$y_prev, gm$pred, method = "spearman", use = "complete.obs")
pr_m <- { o <- gm$y_prev; p <- gm$pred; k <- is.finite(o) & is.finite(p); o <- o[k]; p <- p[k]
  so <- sign(outer(o, o, "-")); sp <- sign(outer(p, p, "-")); m <- so * sp; u <- upper.tri(m); sum(m[u] > 0) / sum(m[u] != 0) }
lim <- c(0, max(c(gm$y_prev, gm$prev_anchored), na.rm = TRUE))
scm <- scale_fill_gradientn(colours = c("#F3EEE6", "#E9B77F", WARM, "#6E2F05"), limits = lim, labels = scales::percent, na.value = "grey88", name = NULL,
                            guide = guide_colourbar(barwidth = 7, barheight = 0.6, ticks = FALSE))
mth <- theme_void(base_size = 13) + theme(plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, lineheight = 0.95), legend.position = "bottom")
q1 <- ggplot(gm) + geom_sf(aes(fill = y_prev), colour = "white", linewidth = 0.05) + scm + labs(title = "Survey 2015-16\n(areas with\nblood samples)") + mth
q2 <- ggplot(gm) + geom_sf(aes(fill = pred), colour = "white", linewidth = 0.05) + scm + labs(title = sprintf("Each area hidden\nin turn (%.0f of 100\npairs in order)", 100 * pr_m)) + mth
gm$q_all <- 100 * (rank(-gm$prev_anchored, na.last = "keep") - 0.5) / sum(is.finite(gm$prev_anchored))   # the deployed order, 0 = worst
q3 <- ggplot(gm) + geom_sf(aes(fill = q_all), colour = "white", linewidth = 0.05) +
  scale_fill_gradientn(colours = c("#0B4F5A", PROXY, "#8FC7CF", "#E7EFF0"), limits = c(0, 100), breaks = c(8, 92), labels = c("worst", "best"), na.value = "grey88", name = NULL,
                       guide = guide_colourbar(barwidth = 7, barheight = 0.6, ticks = FALSE)) +
  labs(title = "Every area ranked,\nincluding those\nwith no survey") + mth
q4 <- ggplot(gm) + geom_sf(aes(fill = p_modplus_cal), colour = "white", linewidth = 0.05) +
  scale_fill_gradient2(low = "#FFF7BC", mid = "#F0A050", high = "#7A0177", midpoint = 0.5, limits = c(0, 1), labels = scales::percent, na.value = "grey88", name = NULL,
                       guide = guide_colourbar(barwidth = 7, barheight = 0.6, ticks = FALSE)) +
  labs(title = sprintf("Chance of %.0f%%\nor more deficient\n(calibrated)", 100 * dep$th_modplus[1])) + mth
sv(q1 + q2 + q3 + q4 + plot_layout(nrow = 1) +
     plot_annotation(caption = if (V5) "Malawi, women's vitamin B12 deficiency. Grey: no blood samples. The chance map uses the anchored, calibrated\nprevalence, shrunk toward the national figure of about 11 per cent, so few areas are confidently above 20 per cent." else "Malawi, women's vitamin B12 deficiency, the best-ranked combination. Grey: no blood samples. The chance map uses the anchored, calibrated prevalence\n(shrunk toward the national figure of about 11 per cent, so few areas are confidently above 20 per cent).",
                     theme = theme(plot.caption = element_text(size = 12, colour = GREY, hjust = 0), plot.caption.position = "plot")),
   "v4_ft_malawi_b12_four.png", 11, 6.2)
cat(sprintf("  Malawi B12 four-panel: in-fill Spearman %.2f, pairs %.0f%%; threshold %.0f%%\n", rho_m, 100 * pr_m, 100 * dep$th_modplus[1]))

# FT slide 21: interpretable measures for the strongest combinations (printed for the slide)
for (cc in list(c("Malawi", "women_b12"), c("Gambia", "women_vitA"), c("Gambia", "child_vitA"), c("Gambia", "women_iron"))) {
  pp <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(country == cc[1], outcome == cc[2], arm == "domain_index", estimand == "infill")
  nt <- read.csv(file.path(P2, "nce_targeting_metrics.csv")) |> filter(country == cc[1], outcome == cc[2], estimand == "infill", arm %in% c("domain_index", "region_mean_jk")) |>
    group_by(arm) |> summarise(c = mean(capture_top20, na.rm = TRUE), .groups = "drop")
  cat(sprintf("  %s %s: pairs %.0f%%; worst fifth holds %.0f%% of deficient (regional %.0f%%)\n", cc[1], cc[2], 100 * pp$pairs,
              100 * nt$c[nt$arm == "domain_index"], 100 * nt$c[nt$arm == "region_mean_jk"]))
}
cat("\nall January-update figures written to", OUT, "\n")
