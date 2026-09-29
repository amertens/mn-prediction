# =============================================================================
# scripts/policy_deck/29_mnf15_v4_figures.R
#
# New figures for the v4 MNF15 talk (docs/slides/MNF15-talk-2026-09-v4.qmd),
# answering Andrew's comments on the reordered deck (27 September). Reads
# committed result tables only (no refitting) and writes results/figures/mnf15_v4/.
#
#   v4_ml_comparison.png      why a simple model: the SL-06 run, every learner on
#                             identical folds, three tests (district, region, country)
#   v4_region_heldout.png     a whole region hidden, wide, both targets named, with
#                             the SuperLearner (same SL-06 run for every arm)
#   v4_accuracy_by_nutrient.png  share of district pairs in the survey's order, with
#                             the survey's own ceiling converted to the same scale
#   v4_domains_child_iron.png top-10 domains for children's iron in Ghana, Malawi and
#                             the four-country model, plus that model's top-10 layers
#   v4_measurability.png      figA relabelled ("little real difference between districts")
#   v4_targeting.png          fig6 on the 16 mappable combinations only
#   v4_civ_b12.png            Cote d'Ivoire B12: anchored predicted prevalence against
#                             the 2007 survey zones, the prevalence map, and how firmly
#                             each district is placed
#
#   Rscript scripts/policy_deck/29_mnf15_v4_figures.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/predictor_plain_names.R")
source("R/admin2_key_hygiene.R")
OUT <- Sys.getenv("FIG_OUT", "results/figures/mnf15_v4"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
V5 <- grepl("mnf15_v5|mnf15_v6", OUT)   # v5 (28 Sep): no text under the accuracy chart
P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"
PROXY <- "#0F7B8A"; WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"; LGREY <- "grey62"
sv <- function(p, f, w, h) { ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 220, bg = "white"); cat("wrote", f, "\n") }
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
th <- function(base = 16) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(size = 12, colour = GREY, hjust = 0),
        plot.subtitle = element_text(size = 14, colour = GREY), plot.title.position = "plot", plot.caption.position = "plot")
CM <- read.csv("results/figures/mnf15/cell_master.csv")   # measurability screen (keep) per cell, from 20_mnf15_figures.R

# mean over cells with a 95% interval (normal approximation across cells), as the full talk's cell_ci
cell_ci <- function(d, by) d |> group_by(across(all_of(by))) |>
  summarise(n = dplyr::n(), est = mean(v), se = sd(v) / sqrt(n), lo = est - 1.96 * se, hi = est + 1.96 * se, .groups = "drop")

# ============================================================================ SL-06 run
# ML_SRC=v6: the post-fix SuperLearner runs (NS-01 rank shards for in-fill and region, the v6 country run)
ML_V6 <- Sys.getenv("ML_SRC") == "v6"
fs <- if (ML_V6) c(list.files(P2, pattern = "^weight_sources_raw_ns01_sl_rank_[a-d]\\.csv$", full.names = TRUE),
                   file.path(P2, "weight_sources_raw_v6_sl_country.csv")) else
  list.files(P2, pattern = "^weight_sources_raw_sl_(gambia|ghana|malawi[0-9]|sierraleone|country)\\.csv$", full.names = TRUE)
if (ML_V6) chk(length(fs) == 5 && all(file.exists(fs)), "four NS-01 rank shards and the v6 country run")
SL <- bind_rows(lapply(fs, read.csv)) |> filter(is.finite(spearman)) |>
  group_by(country, outcome, target, estimand, arm) |> summarise(v = mean(spearman), .groups = "drop")   # in-fill: mean over its 3 fold draws
chk(abs(mean(SL$v[SL$estimand == "infill" & SL$target == "level" & SL$arm == "domain_index"]) - (if (ML_V6) 0.389 else 0.392)) < 0.005, "in-fill index 0.39")
EST <- c(infill = "One district in five hidden", region = "A whole region hidden", country = "A whole country hidden")

# --- why a simple model: every learner on the same folds, average status
ML <- c(domain_index = "Domain-PC index (nothing tuned)", sl_nnls = "SuperLearner ensemble (12 methods)",
        sl_lrn_ranger = "Random forest", sl_lrn_xgb = "Gradient boosting", sl_lrn_lasso = "Lasso", sl_lrn_enet = "Elastic net")
dm <- SL |> filter(target == "level", arm %in% names(ML)) |> cell_ci(c("estimand", "arm")) |>
  mutate(Model = factor(ML[arm], levels = rev(ML)), est_lab = factor(EST[estimand], levels = EST),
         col = ifelse(arm == "domain_index", PROXY, ifelse(arm == "sl_nnls", WARM, LGREY)))
pml <- ggplot(dm, aes(est, Model)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = col), height = 0, linewidth = 1.1) +
  geom_point(aes(colour = col), size = 4.6) +
  scale_colour_identity() + facet_wrap(~ est_lab, nrow = 1) +
  scale_x_continuous(limits = c(-0.02, 0.6), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Agreement with the survey's order of districts (correlation; 0 = chance)", y = NULL,
       caption = "Average status (biomarker level); mean over nutrient-country combinations with a 95% interval. Every method scored on the same hidden districts.") +
  th(16) + theme(strip.text = element_text(face = "bold", size = 14.5), axis.text.y = element_text(size = 14.5, colour = INK),
                 panel.spacing.x = unit(1.4, "lines"))
sv(pml, "v4_ml_comparison.png", 13, 4.6)
cat("  ML comparison (level):\n"); print(as.data.frame(dm |> select(estimand, arm, est) |> tidyr::pivot_wider(names_from = estimand, values_from = est)), digits = 3)

# --- a whole region hidden, both targets, with the SuperLearner (one run for every arm)
RG <- c(domain_index = "Domain-PC index", sl_nnls = "SuperLearner ensemble (12 methods)",
        sl_lrn_spatial_plus_domain = "Neighbour smoother + Domain-PC index", sl_lrn_spatial_gam = "Neighbour smoother (map of nearby districts)",
        domain_enet = "Elastic net on the domain components", sl_lrn_enet = "Elastic net on all layers")
if (ML_V6) RG <- RG[names(RG) != "domain_enet"]   # not an arm of the post-fix SuperLearner runs; dropped rather than mixed in from another run
TGT <- c(level = "Average status (mean biomarker level)", prev = "Share deficient (prevalence)")
common <- SL |> filter(estimand == "region", arm %in% names(RG)) |> group_by(target, country, outcome) |>
  filter(n_distinct(arm) == length(RG)) |> ungroup()   # cells every arm scored (one prevalence cell lacks the elastic net)
dr <- common |> cell_ci(c("target", "arm")) |>
  mutate(Model = factor(RG[arm], levels = rev(RG)), Target = factor(TGT[target], levels = TGT),
         col = ifelse(arm == "domain_index", PROXY, ifelse(arm == "sl_nnls", WARM, LGREY)))
chk(all(dr$n[dr$target == "level"] == 18) && all(dr$n >= 17), "18 (level) and 17+ (prevalence) region cells per arm")
prg <- ggplot(dr, aes(est, Model)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = col), height = 0, linewidth = 1.2) +
  geom_point(aes(colour = col), size = 5) +
  geom_text(aes(x = hi, label = sprintf("%.2f", est), colour = col), hjust = -0.3, size = 5, fontface = "bold") +
  scale_colour_identity() + facet_wrap(~ Target, nrow = 1) +
  scale_x_continuous(limits = c(-0.2, 0.62), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Agreement with the survey's order of districts (correlation; 0 = chance)", y = NULL,
       caption = "All six nutrient outcomes: mean over the 18 nutrient-country combinations (17 for prevalence), 95% interval.
Each region of each country hidden in turn; every method on the same folds.") +
  th(17) + theme(strip.text = element_text(face = "bold", size = 16), axis.text.y = element_text(size = 15.5, colour = INK),
                 panel.spacing.x = unit(2, "lines"))
sv(prg, "v4_region_heldout.png", 13.2, 5.0)
cat("  region hold-out:\n"); print(as.data.frame(dr |> select(target, arm, est, lo, hi)), digits = 3)

# ============================================================================ accuracy by nutrient, pairs
# The ceiling is the variance-components correlation between the survey's district value and
# the true district value (average status). On the pairs scale a correlation r corresponds,
# for bivariate-normal values, to 1/2 + asin(r)/pi of pairs in the same order (Greiner's relation).
VC <- read.csv(file.path(P2, "variance_components_ceiling.csv")) |> filter(rung == "admin2", target == "level") |>
  select(country, outcome, ceil = ceiling_vc)
PP <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(arm == "domain_index") |>
  inner_join(CM |> filter(keep) |> select(country, outcome, nutrient, pop), by = c("country", "outcome")) |>
  left_join(VC, by = c("country", "outcome"))
nb <- PP |> mutate(Nutrient = paste0(nutrient, ", ", tolower(pop))) |> group_by(Nutrient) |>
  summarise(inc = mean(pairs[estimand == "infill"]), ho = mean(pairs[estimand == "country"]),
            ceil = mean(0.5 + asin(ceil[estimand == "infill"]) / pi), .groups = "drop") |>
  filter(is.finite(inc)) |> arrange(inc) |> mutate(Nutrient = factor(Nutrient, levels = Nutrient))
nl <- nb |> pivot_longer(c(inc, ho), names_to = "what", values_to = "v") |>
  mutate(what = factor(c(inc = "Inside a surveyed country", ho = "Whole country held out")[what], levels = c("Inside a surveyed country", "Whole country held out")))
pan <- ggplot(nl, aes(y = Nutrient)) +
  geom_col(data = nb, aes(x = ceil * 100), fill = "grey88", width = 0.66) +
  geom_text(data = nb, aes(x = ceil * 100, label = sprintf("best possible %.0f%%", 100 * ceil)), hjust = -0.08, size = 4.5, colour = GREY) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = GREY) +
  geom_point(aes(x = v * 100, colour = what), size = 5.4, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(PROXY, WARM), name = NULL) +
  coord_cartesian(xlim = c(40, 100)) +
  scale_x_continuous(breaks = c(50, 60, 70, 80, 90, 100), labels = c("50%\ncoin toss", "60%", "70%", "80%", "90%", "100%")) +
  labs(x = "Share of district pairs put in the survey's order", y = NULL,
       subtitle = "Grey bar: the best any model could do. Most districts' survey figure rests on one community of about a dozen people,
so no map, however good, can match the survey perfectly.",
       caption = if (V5) NULL else "Measurable combinations only. Ceiling from the surveys' own variance components (average status), converted to pairs assuming normal errors.") +
  th(17) + theme(legend.position = "bottom", legend.text = element_text(size = 15), axis.text.y = element_text(size = 15.5, colour = INK),
                 plot.subtitle = element_text(size = 14.5, colour = INK))
sv(pan, "v4_accuracy_by_nutrient.png", 13, 5.6)
cat("  pairs by nutrient:\n"); print(as.data.frame(nb), digits = 3)

# ============================================================================ domains, children's iron
DOM <- read.csv(file.path(P2, "index_importance_domains.csv")) |> filter(outcome == "child_iron", target == "level")
COL <- read.csv(file.path(P2, "index_importance_columns.csv")) |> filter(outcome == "child_iron", target == "level", scope == "pooled")
DL <- c("Satellite embedding" = "Satellite imagery summary", "Climate and weather" = "Climate", "Soil characteristics" = "Soil",
        "Ecosystem productivity/greenness" = "Vegetation greenness", "Infection and inflammation burden" = "Infection (malaria, worms)",
        "Water and sanitation" = "Water and sanitation", "Infant and young child feeding" = "Infant and young child feeding",
        "Child anthropometry" = "Child growth and overweight", "Education, employment, SES" = "Schooling and employment",
        "Household diet and consumption (HCES)" = "Household diet (budget surveys)", "Livestock density" = "Livestock",
        "Agricultural production, land use" = "Crops and land use", "Built environment" = "Built-up land",
        "Household assets and characteristics" = "Household assets", "Immunisation" = "Immunisation", "Anaemia and haemoglobin" = "Anaemia",
        "Food prices and supply" = "Food prices", "Water and coast proximity" = "Distance to water and coast",
        "Fertility, reproductive health" = "Fertility", "Healthcare access" = "Access to care")
FITS <- c(Ghana = "Ghana
only", Malawi = "Malawi
only", all = "All four
countries")
dd <- DOM |> filter((scope == "country" & fit %in% c("Ghana", "Malawi")) | scope == "pooled") |>
  mutate(key = ifelse(scope == "pooled", "all", fit)) |> group_by(key) |> arrange(desc(share)) |> mutate(rk = row_number()) |>
  filter(rk <= 10) |> ungroup() |>
  mutate(dlab = ifelse(domain %in% names(DL), DL[domain], domain), col = factor(FITS[key], levels = FITS),
         env = domain %in% c("Satellite embedding", "Climate and weather", "Soil characteristics", "Ecosystem productivity/greenness"))
chk(!any(grepl("modelled|MODELLED", dd$dlab)), "no 'modelled' in the domain labels")
ord <- DOM |> filter(scope == "pooled") |> mutate(dlab = ifelse(domain %in% names(DL), DL[domain], domain)) |> arrange(share)
dd$dlab <- factor(dd$dlab, levels = intersect(ord$dlab, unique(dd$dlab)))
pdm <- ggplot(dd, aes(col, dlab)) +
  geom_tile(aes(fill = env, alpha = share), colour = "white", linewidth = 1.2) +
  geom_text(aes(label = rk, colour = share > 0.09), size = 5.2, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = PROXY, `FALSE` = NAVY), labels = c(`TRUE` = "Measured everywhere, every year
(satellite, climate, soil)", `FALSE` = "Other domains"), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Top ten domains in each model (1 = most weight)") +
  th(14) + theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold", size = 13.5, colour = INK),
                 axis.text.y = element_text(size = 13, colour = INK), legend.position = "bottom", legend.text = element_text(size = 12),
                 plot.title = element_text(face = "bold", size = 14.5)) + guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 0.9)))
SHORT <- c(mics_c_fg_vita_fv = "Young children eating vitamin A-rich fruit and vegetables",
           ihme_overweightprevalence = "Child overweight (IHME)", lcover_crops_frac_t0 = "Cropland cover (satellite)")
tc <- COL |> arrange(rank) |> slice_head(n = 10) |>
  mutate(lab = ifelse(column %in% names(SHORT), SHORT[column], vapply(column, plain_of, "")),
         lab = stringr::str_wrap(lab, 30), w = beta_std / max(abs(beta_std)),
         dir = ifelse(beta_std > 0, "Higher value, more deficiency", "Higher value, less deficiency"))
chk(nrow(tc) == 10 && all(tc$loco_sign_agree == 4), "top-10 pooled layers, direction kept in every leave-one-country-out fit")
tc$lab <- factor(tc$lab, levels = rev(tc$lab))
ptc <- ggplot(tc, aes(w, lab, fill = dir)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c("Higher value, more deficiency" = WARM, "Higher value, less deficiency" = NAVY), name = NULL) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, 0, 1), labels = c("less\ndeficiency", "0", "more\ndeficiency")) +
  labs(x = "Relative weight in the four-country model", y = NULL, title = "Top ten layers, four-country model",
       caption = "Children's iron, average status. Markers of where deficiency is, not causes;
every one keeps its direction whichever country is left out.") +
  th(13) + theme(axis.text.y = element_text(size = 11.5, colour = INK, lineheight = 0.9), legend.position = "bottom",
                 legend.text = element_text(size = 11.5), plot.title = element_text(face = "bold", size = 14)) + guides(fill = guide_legend(ncol = 1))
sv(pdm + ptc + plot_layout(widths = c(0.75, 1)), "v4_domains_child_iron.png", 13, 6.2)
cat("  pooled top-10 layers:\n"); print(as.data.frame(tc[, c("column", "beta_std", "lab")]))

# ============================================================================ measurability, relabelled
scr <- CM |> mutate(reason = case_when(
    prev < 0.02 & (is.na(ceiling) | ceiling < 0.30) ~ "Too rare, and little real difference between districts",
    prev < 0.02                                      ~ "Too rare to rank (under 2%)",
    TRUE                                             ~ "Survey shows little real difference between districts"),
    ctry = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"),
    lab = paste0(ctry, ", ", tolower(pop), "'s ", sub("^vitamin b12$", "vitamin B12", sub("^vitamin a$", "vitamin A", tolower(nutrient)))))
dropped <- scr |> filter(!keep) |> arrange(reason, prev)
chk(nrow(dropped) == 8, "eight combinations fail the screen")
pA <- ggplot(dropped, aes(x = reorder(lab, prev), y = pmax(prev, 0.001) * 100, fill = reason)) +
  geom_col(width = .68) + coord_flip() +
  scale_fill_manual(values = c("Too rare to rank (under 2%)" = WARM, "Survey shows little real difference between districts" = PROXY,
                               "Too rare, and little real difference between districts" = "#7A4E8C"), name = NULL) +
  labs(y = "National prevalence (%)", x = NULL, subtitle = "Decided from the survey alone, before any model is fitted",
       caption = "\"Little real difference\": the survey's district figures differ, but no more than sampling noise from one or two communities per district would produce.") +
  th(16) + theme(legend.position = "bottom", legend.text = element_text(size = 13), axis.text.y = element_text(size = 14, colour = INK)) +
  guides(fill = guide_legend(ncol = 1))
sv(pA, "v4_measurability.png", 12.2, 5.4)

# ============================================================================ targeting, mappable cells only
NT <- read.csv(file.path(P2, "nce_targeting_metrics.csv")) |> filter(estimand == "infill") |>
  inner_join(CM |> filter(keep) |> select(country, outcome), by = c("country", "outcome"))
chk(n_distinct(paste(NT$country, NT$outcome)) == 16, "16 mappable combinations")
cap <- NT |> group_by(country, outcome, arm) |> summarise(c = mean(capture_top20, na.rm = TRUE), .groups = "drop") |> group_by(arm) |> summarise(c = mean(c, na.rm = TRUE), .groups = "drop")
gc6 <- function(a) cap$c[cap$arm == a]
F6 <- data.frame(who = c("Perfect knowledge (unreachable)", "Domain-PC index", "Survey's regional averages"),
                 v = c(gc6("oracle_ceiling"), gc6("domain_index"), gc6("region_mean_jk")), col = c("grey70", PROXY, WARM))
F6$who <- factor(F6$who, levels = rev(F6$who))
p6 <- ggplot(F6, aes(v, who, fill = col)) + geom_col(width = 0.62) +
  geom_vline(xintercept = 0.2, linetype = "dashed", colour = GREY) +
  annotate("text", x = 0.205, y = 0.45, label = "a random fifth of districts: 20%", hjust = 0, size = 5, colour = GREY) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * v)), hjust = -0.2, size = 6.4, colour = INK) +
  scale_fill_identity() + scale_x_continuous(labels = scales::percent, limits = c(0, 0.5), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.5))) +
  labs(x = "Share of the country's deficient people living in the fifth of districts each method picks", y = NULL,
       caption = "Average over the 16 mappable nutrient-country combinations in The Gambia, Ghana, Sierra Leone and Malawi
(the 8 that fail the screen are left out); each district hidden in turn.") +
  th(17) + theme(axis.text.y = element_text(size = 16, colour = INK), panel.grid.major.y = element_blank())
sv(p6, "v4_targeting.png", 13, 4.8)
cat(sprintf("  targeting, 16 mappable cells: index %.3f, regional %.3f, perfect %.3f, null %.3f\n",
            gc6("domain_index"), gc6("region_mean_jk"), gc6("oracle_ceiling"), gc6("null_train_mean")))

# ============================================================================ Cote d'Ivoire B12
# Anchored prevalence, design A1 of scripts/protocol_v2/35_anchor_and_rank.R:
#   p_d = expit(logit(anchor) + rho * sd * z_d)
# z_d = the transported climate-and-soil B12 index standardised within Cote d'Ivoire (higher = worse);
# sd = between-district SD of logit prevalence in the training surveys (Malawi at district level, as AR-01);
# rho = how well that ranking held up for women's B12 when each training country was held out (AR-01
# rho_obs, climate-and-soil set, level ranking). anchor = the 2007 survey's national 18.1%, so the level is
# the survey's and the district pattern is the model's: the design the talk recommends.
sf::sf_use_s2(FALSE)
ANCHOR <- 0.181
VM <- c(North = 48.6, `North East` = 40.4, `North West` = 36.4, West = 29.8, Central = 20.0,
        `Central West` = 11.5, `South East` = 7.0, South = 6.2, Abidjan = 0.0)
R <- read.csv(file.path(PD, "civ_climate_soil_ranking.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8") |> filter(outcome == "women_b12")
chk(nrow(R) == 33 && cor(R$index, R$rank) < -0.9, "33 CIV districts; higher index = worse (rank 1)")
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE) |> filter(outcome == "women_b12", is.finite(y_prev), is.finite(n_eff), n_eff > 0)
lg <- function(p) { p <- pmin(pmax(p, 0.005), 0.995); log(p / (1 - p)) }
sds <- TG |> mutate(unit = ifelse(country == "Malawi", Admin1, paste(Admin1, Admin2))) |> group_by(country, unit) |>
  summarise(p = weighted.mean(y_prev, n_eff), .groups = "drop") |> group_by(country) |> summarise(sd = sd(lg(p)), n = dplyr::n(), .groups = "drop")
AR <- read.csv(file.path(P2, "anchor_and_rank.csv")) |> filter(outcome == "women_b12", set == "climate_soil", rank_from == "level", design == "A1_anchor_rank") |>
  distinct(country, rho_obs)
rho <- max(0, mean(AR$rho_obs)); sdv <- mean(sds$sd)
R$z <- as.numeric(scale(R$index)); R$p <- plogis(qlogis(ANCHOR) + rho * sdv * R$z)
cat(sprintf("  CIV B12: rho %.3f (%s), sd %.3f (%s); district range %.1f%% to %.1f%%\n", rho, paste(sprintf("%s %.2f", AR$country, AR$rho_obs), collapse = ", "),
            sdv, paste(sprintf("%s %.2f", sds$country, sds$sd), collapse = ", "), 100 * min(R$p), 100 * max(R$p)))
Bc <-readRDS("dashboard/data/oos_cote_divoire.rds")$boundaries
gc <- dplyr::left_join(Bc, R[, c("Admin1", "Admin2", "p", "rank")], by = admin2_join_by(Bc, R))
chk(sum(is.finite(gc$p)) == 33, "33 CIV districts joined to boundaries")
ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(gc))))
gc$lon <- ctr[, 1]; gc$lat <- ctr[, 2]
zone <- function(a1, a2, lat, lon) {   # the compass crosswalk of 24_civ_2007_survey_check.R
  if (grepl("Abidjan", a1) || grepl("Abidjan", a2)) return("Abidjan")
  if (lat >= 8.7) return(if (lon <= -6.6) "North West" else if (lon >= -4.4) "North East" else "North")
  if (lat >= 6.7) return(if (lon <= -6.6) "West" else if (lon <= -5.2) "Central West" else "Central")
  if (lon >= -4.4) "South East" else "South"
}
gc$zone <- mapply(zone, gc$Admin1, gc$Admin2, gc$lat, gc$lon)
z <- sf::st_drop_geometry(gc) |> group_by(zone) |> summarise(pred = 100 * mean(p), mrank = mean(rank), n = dplyr::n(), .groups = "drop") |>
  mutate(meas = VM[zone])
rho_z <- cor(z$meas, z$mrank, method = "spearman")
chk(abs(abs(rho_z) - 0.95) < 0.02, "zone crosswalk reproduces figQ (0.95)")
cat(sprintf("  CIV B12 zones: Spearman %.2f; predicted zone range %.1f to %.1f%%, measured 0 to 48.6%%; MAE %.1f points\n",
            rho_z, min(z$pred), max(z$pred), mean(abs(z$pred - z$meas))))
psc <- ggplot(z, aes(meas, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(size = n), colour = PROXY) +
  ggrepel::geom_text_repel(aes(label = zone), size = 4.4, seed = 2, box.padding = 0.4, min.segment.length = 0.3) +
  scale_size_continuous(range = c(3, 7), guide = "none") +
  coord_equal(xlim = c(0, 52), ylim = c(0, 52)) +
  labs(x = "Measured, 2007 survey (%)", y = "Predicted by the model (%)") + th(15)
U <- read.csv(file.path(PD, "civ_rank_uncertainty_all.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8") |> filter(outcome == "women_b12", domain_set == "cs")
gc <- dplyr::left_join(gc, U[, c("Admin1", "Admin2", "rank_width")], by = admin2_join_by(gc, U))
chk(sum(is.finite(gc$rank_width)) == 33, "33 CIV rank intervals joined")
mt <- theme_void(base_size = 15) + theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5), legend.position = "bottom",
                                        legend.text = element_text(size = 12))
pmap <- ggplot(gc) + geom_sf(aes(fill = 100 * p), colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(colours = c("#F3EEE6", "#E9B77F", WARM, "#6E2F05"), limits = c(0, NA), name = "% deficient",
                       guide = guide_colourbar(barwidth = 11, barheight = 0.7, ticks = FALSE, title.position = "top")) +
  labs(title = "Predicted B12 deficiency, women") + mt
wmax <- round(max(gc$rank_width, na.rm = TRUE))
punc <- ggplot(gc) + geom_sf(aes(fill = rank_width), colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(colours = rev(c("#F2F2F2", "#BFC6CC", "#7C8B95", "#3D4A54")), limits = c(0, wmax), breaks = c(0.5, wmax - 0.5),
                       labels = c("firm", sprintf("could move\n%d places", wmax)), name = NULL,
                       guide = guide_colourbar(barwidth = 11, barheight = 0.7, ticks = FALSE)) +
  labs(title = "How firmly each district is placed") + mt
sv(psc + pmap + punc + plot_layout(widths = c(1, 1, 1)), "v4_civ_b12.png", 13, 5.0)
write.csv(sf::st_drop_geometry(gc)[, c("Admin1", "Admin2", "zone", "p", "rank", "rank_width")], file.path(PD, "v4_civ_b12_anchored.csv"), row.names = FALSE)
write.csv(z, file.path(PD, "v4_civ_b12_zones.csv"), row.names = FALSE)
# ============================================================================ appendix: women's domains, "Anaemia" not "Anaemia (modelled)"
# As figD_women_domains.png of 20_mnf15_figures.R (pooled fit, top five domains per women's outcome), relabelled.
DW <- read.csv(file.path(P2, "index_importance_domains.csv")) |> filter(scope == "pooled", target == "level", grepl("^women_", outcome)) |>
  group_by(outcome) |> slice_max(share, n = 5) |> ungroup() |>
  mutate(olab = recode(outcome, women_b12 = "Women's B12", women_folate = "Women's folate", women_iron = "Women's iron", women_vitA = "Women's vitamin A"),
         dlab = recode(domain, "Satellite embedding" = "Satellite imagery", "Climate and weather" = "Climate", "Soil characteristics" = "Soil",
                       "Ecosystem productivity/greenness" = "Plant productivity", "Infection and inflammation burden" = "Infectious disease",
                       "Education, employment, SES" = "Schooling and employment", "Child anthropometry" = "Child wasting and underweight",
                       "Infant and young child feeding" = "Infant feeding (meat and fish)", "Anaemia and haemoglobin" = "Anaemia",
                       "Built environment" = "Built environment"))
DPAL <- c("Climate" = "#0F7B8A", "Soil" = "#7C4B24", "Satellite imagery" = "#274C77", "Plant productivity" = "#4C9A2A",
          "Infectious disease" = "#B23A48", "Schooling and employment" = "#B45309", "Child wasting and underweight" = "#8A5FA8",
          "Infant feeding (meat and fish)" = "#C2185B", "Anaemia" = "#D98E04", "Built environment" = "#5E6472")
chk(all(DW$dlab %in% names(DPAL)), "every women's domain has a colour")
DW <- DW |> mutate(dlab = factor(dlab, levels = names(DPAL))) |> arrange(olab, share) |> mutate(row = factor(seq_len(dplyr::n())))
pw <- ggplot(DW, aes(share * 100, row, fill = dlab)) + geom_col(width = .72) +
  scale_y_discrete(breaks = DW$row, labels = DW$dlab) + facet_wrap(~ olab, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = DPAL, name = NULL) +
  labs(subtitle = "Share of the model's weight carried by each domain, top five per outcome, women", x = "Share of the model's weight (%)", y = NULL,
       caption = "Pooled four-country fit. These say where deficiency is. They are not causes and not things to change.") +
  th(16) + theme(strip.text = element_text(face = "bold", size = 15), axis.text.y = element_text(size = 12), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
sv(pw, "v4_women_domains.png", 12.6, 6.2)

cat("\nall v4 figures written to", OUT, "\n")
