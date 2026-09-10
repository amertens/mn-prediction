# =============================================================================
# scripts/policy_deck/01_figures_main.R
#
# Figures 1-6 and 8 for the 15-minute policy deck
# (docs/slides/MN-proxy-policy-deck-2026-09.qmd).
#
# READS ONLY result tables. Fits nothing. Every number traces to a CSV named in
# docs/slides/POLICY_DECK_BRIEF_2026-09.md.
#
#   Rscript scripts/policy_deck/01_figures_main.R
# -> results/figures/policy_deck/*.png   (16:9, base font 18pt)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(grid)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

P2  <- "results/tables/protocol_v2"
CL  <- "results/tables/cluster_level"
OUT <- "results/figures/policy_deck"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rd <- function(...) read.csv(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)

# ── house style ──────────────────────────────────────────────────────────────
PROXY  <- "#0F7B8A"   # the model of record
SURVEY <- "#D2691E"   # survey-based baseline
GEO    <- "#6A51A3"   # DHS-style geostatistical model
OTHER  <- "#9AA0A6"   # everything else
INK    <- "#1A1A1A"
CHANCE_D <- 0.08      # 95th pct of the no-signal null, districts
CHANCE_R <- 0.16      # ... regions

theme_deck <- function(base = 18) {
  theme_minimal(base_size = base) +
    theme(text            = element_text(colour = INK),
          plot.title      = element_text(face = "bold", size = base * 1.15, margin = margin(b = 4)),
          plot.subtitle   = element_text(size = base * 0.98, colour = "grey20", margin = margin(b = 12)),
          plot.caption    = element_text(size = base * 0.62, colour = "grey45", hjust = 0),
          axis.title      = element_text(size = base * 0.85),
          strip.text      = element_text(face = "bold", size = base * 0.9),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          plot.margin     = margin(12, 18, 8, 12))
}
save16x9 <- function(p, file, w = 12.2, h = 6.3, dpi = 200) {
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = dpi, bg = "white")
  cat("wrote", file, "\n")
}

# =============================================================================
# FIGURE 1 - model comparison lollipop, three panels
# =============================================================================
BM <- rd(P2, "benchmarks_v2_summary.csv")
WS <- rd(P2, "weight_sources_summary.csv")

lab_map <- c(domain_index        = "Proxy index (what we propose)",
             spatial_plus_domain = "Neighbour smoother + proxies",
             spatial             = "Neighbour smoother alone",
             domain_enet         = "Penalised regression",
             region_mean_jk      = "Survey's own regional average")

get_bm <- function(est, arm, tgt = "level") {
  r <- BM[BM$estimand == est & BM$arm == arm & BM$target == tgt, ]
  if (!nrow(r)) NA_real_ else r$mean_spearman[1]
}
get_ws <- function(est, arm, tgt = "level") {
  r <- WS[WS$estimand == est & WS$arm == arm & WS$target == tgt, ]
  if (!nrow(r)) NA_real_ else r$mean_spearman[1]
}

panels <- c(infill = "Inside a surveyed country",
            region = "A region held out",
            country = "A country with no survey")

rows <- list()
for (est in names(panels)) {
  for (a in names(lab_map)) rows[[length(rows) + 1L]] <-
    data.frame(estimand = est, method = lab_map[[a]], value = get_bm(est, a))
  rows[[length(rows) + 1L]] <-
    data.frame(estimand = est, method = "Twenty public layers", value = get_ws(est, "sparse20"))
}
F1 <- bind_rows(rows)

# region_mean_jk and the smoothers are undefined outside a surveyed country;
# the benchmark table simply has no row, so keep them as an explicit label.
F1$cannot <- !is.finite(F1$value)
F1$value_plot <- ifelse(F1$cannot, NA_real_, F1$value)

ord <- F1 |> filter(estimand == "infill") |> arrange(value_plot) |> pull(method)
F1$method   <- factor(F1$method, levels = ord)
F1$estimand <- factor(panels[F1$estimand], levels = unname(panels))
F1$col <- ifelse(F1$method == "Proxy index (what we propose)", PROXY,
          ifelse(F1$method == "Survey's own regional average", SURVEY, OTHER))

chance <- data.frame(estimand = factor(unname(panels), levels = unname(panels)),
                     x = c(CHANCE_D, CHANCE_R, CHANCE_D))

p1 <- ggplot(F1, aes(x = value_plot, y = method)) +
  geom_rect(data = chance, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = x, ymin = -Inf, ymax = Inf),
            fill = "grey88", alpha = 0.9) +
  geom_segment(aes(x = 0, xend = value_plot, yend = method, colour = col), linewidth = 1.5) +
  geom_point(aes(colour = col), size = 6) +
  geom_text(aes(label = ifelse(is.na(value_plot), "", sprintf("%.2f", value_plot))),
            hjust = -0.45, size = 5.1, colour = INK) +
  geom_text(data = subset(F1, cannot), aes(x = 0.02, label = "needs a survey here"),
            hjust = 0, size = 4.6, colour = "grey45", fontface = "italic") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 0.56), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  facet_wrap(~ estimand) +
  labs(
       subtitle = "Does the predicted order of districts match the survey's? Grey band = chance.",
       x = "Ranking accuracy", y = NULL) +
  theme_deck()
save16x9(p1, "fig1_model_comparison.png", h = 6.6)

# =============================================================================
# FIGURE 2 - "simpler wins" strip
# =============================================================================
sl <- rd(P2, "sl_rank_loss_scores.csv")
sl_mean <- sl |> group_by(arm) |> summarise(v = mean(spearman, na.rm = TRUE), .groups = "drop")
gv <- function(a) { x <- sl_mean$v[sl_mean$arm == a]; if (!length(x)) NA_real_ else x }

# All five arms come from ONE table and ONE set of folds, so they are directly
# comparable. The eighth arm in that file is the no-model national average
# (-0.19); it is the null, not a method, so it is described in the caption
# rather than plotted on the same axis.
F2 <- data.frame(
  method = c("Proxy index (no tuning)",
             "Machine-learning ensemble,\ntuned for ranking",
             "Random forest",
             "Machine-learning ensemble,\ntuned for error",
             "Penalised regression"),
  value  = c(gv("domain_index"),
             max(c(gv("rank_nnls"), gv("rank_discrete")), na.rm = TRUE),
             gv("rf"),
             max(c(gv("mse_nnls"), gv("mse_discrete")), na.rm = TRUE),
             gv("enet")))
F2 <- F2[is.finite(F2$value), ]
F2$method <- factor(F2$method, levels = F2$method[order(F2$value)])
F2$col <- ifelse(grepl("^Proxy index", F2$method), PROXY, OTHER)

p2 <- ggplot(F2, aes(x = value, y = method)) +
  annotate("rect", xmin = -Inf, xmax = CHANCE_D, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_segment(aes(x = 0, xend = value, yend = method, colour = col), linewidth = 1.6) +
  geom_point(aes(colour = col), size = 7) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.4, size = 5.4, colour = INK) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 0.36), breaks = seq(0, 0.3, 0.1), expand = c(0, 0)) +
  labs(
       subtitle = "With 14 to 87 districts per country, methods that tune themselves overfit.",
       x = "Ranking accuracy", y = NULL,
       caption = "Same folds for all five. Person-level prediction is not shown: it is no better than a coin toss.") +
  theme_deck()
save16x9(p2, "fig2_simpler_wins.png", h = 6.2)

# =============================================================================
# FIGURE 3 - top five predictors per outcome
# =============================================================================
TOP <- rd(P2, "index_importance_top.csv") |> filter(target == "level")

PLAIN <- c(
  # added 2026-09-09 after RR-10 (cluster-model DHS columns moved several survey aggregates into the leading lists)
  dhs_CN_NUTS_C_HA2 = "Stunted children (survey)", dhs_CN_NUTS_C_WH2 = "Wasted children (survey)",
  dhs_w_height_low = "Short-stature women", dhs_w_birth_interval_short = "Short birth intervals",
  dhs_w_health_insurance = "Women with health insurance", dhs_FP_CUSA_W_MOD = "Modern contraceptive use",
  dhs_w_modern_fp = "Modern family planning use", dhs_hh_soap_available = "Households with soap",
  dhs_w_barrier_permission = "Need permission to seek care", dhs_w_bmi_low = "Thin women (BMI)",
  dhs_AN_NUTS_W_THN = "Thin women (survey)", dhs_c_fg_roots = "Children eating roots and tubers",
  dhs_c_fg_legumes = "Children eating legumes", dhs_c_fg_grains = "Children eating grains",
  dhs_c_fg_dairy = "Children eating dairy", dhs_c_fg_flesh = "Children eating meat or fish",
  dhs_c_fg_eggs = "Children eating eggs", dhs_c_fg_other_fruitveg = "Children eating other fruit and vegetables",
  dhs_CN_BRFS_C_EXB = "Exclusive breastfeeding", dhs_hh_cows = "Household cattle ownership",
  dhs_hh_cows_any = "Any household cattle", dhs_hh_cattle = "Household cattle", dhs_hh_cattle_any = "Any household cattle",
  dhs_hh_goats = "Household goats", dhs_hh_goats_any = "Any household goats", dhs_hh_sheep = "Household sheep",
  dhs_w_decides_earnings = "Women deciding on earnings", dhs_w_occ_agric = "Women in farm work",
  dhs_hh_improved_water = "Improved water source", dhs_hh_crowding = "Household crowding",
  dhs_w_deworm_pregnancy = "Deworming in pregnancy", dhs_hh_itn_any = "Bed net in household",
  dhs_ML_NETP_H_IT2 = "Bed net per two people", dhs_c_deworm_6mo = "Child deworming",
  dhs_w_sib_maternal_any = "Sibling maternal death", ihme_severeanemia = "Modelled severe anaemia",
  ihme_allanemia = "Modelled anaemia (any)", ihme_moderateanemia = "Modelled moderate anaemia",
  ihme_mildanemia = "Modelled mild anaemia", ihme_wastingprevalence = "Modelled child wasting",
  espen_sth_cov_mean = "Deworming coverage", tclim_pdsi_t0 = "Drought index", lcover_crops_frac_t0 = "Cropland cover",
  wapor_sd_t0 = "Vegetation seasonality", grassland_frac = "Grassland share", wdist_coast_km_mean = "Distance to coast",
  wdist_coast_km_min = "Distance to coast (nearest)", wdist_perm_km_mean = "Distance to permanent water",
  wdist_any_km_mean = "Distance to any water", elevation = "Elevation", glw_ruminant_share = "Ruminant share of livestock",
  glw_cattle_km2 = "Cattle density", glw_pigs_km2 = "Pig density", spam_share_cereals = "Cereal share of cropland",
  spam_share_oilcrops = "Oil-crop share of cropland", soil_phosphorus_stdev_0_20 = "Soil phosphorus variability",
  soil_zinc_stdev_0_20 = "Soil zinc variability", map_blooddisorders201201africahbcallelefrequency = "Haemoglobin C gene frequency",
  map_blooddisorders201201globalsicklehaemoglobinhbsallelefrequency = "Sickle-cell gene frequency",
  map_sy_pf_mortality_rate = "Malaria mortality", map_sy_pf_incidence_rate = "Malaria incidence", map_sy_pf_parasite_rate = "Malaria parasite rate",
  rwi_sd = "Wealth index spread", lcover_water_seasonal_frac_t0 = "Seasonal water cover",
  lcover_grass_frac_t0 = "Grassland cover", npp_gpp_t0 = "Vegetation growth (gross)",
  npp_npp_t0 = "Vegetation growth (net)", ihme_stuntingprevalence = "Modelled child stunting",
  ihme_underweightprevalence = "Modelled child underweight", wpop_share_under5 = "Share of people under five",
  wpop_dependency_ratio = "Dependency ratio", dhs_w_primary_edu = "Women with primary schooling",
  dhs_w_no_education = "Women with no schooling", fprice_staple_rel = "Relative staple food price",
  glw_tlu_per_capita = "Livestock per person", glw_sheep_km2 = "Sheep density",
  glw_cattle_km2 = "Cattle density", glw_pigs_km2 = "Pig density",
  glw_ruminant_share = "Ruminant share of livestock", dhs_w_owns_house = "Women owning their home",
  dhs_w_working = "Women in paid work", dhs_hh_cows = "Household cattle ownership",
  dhs_hh_cows_any = "Household owns any cattle", spam_share_cereals = "Cereal share of cropland",
  spam_share_roots = "Root-crop share of cropland", spam_share_oilcrops = "Oil-crop share of cropland",
  map_blooddisorders201201africahbcallelefrequency = "Haemoglobin C gene frequency",
  soil_phosphorus_stdev_0_20 = "Soil phosphorus variability",
  ihme_severeanemia = "Modelled severe anaemia", ihme_anemia = "Modelled anaemia",
  dhs_FP_CUSA_W_MOD = "Modern contraceptive use", espen_sth_cov_mean = "Deworming coverage",
  dhs_w_health_insurance = "Health insurance", dhs_AN_NUTS_W_THN = "Thin women",
  dhs_CN_NUTS_C_HA2 = "Stunted children (survey)", dhs_CN_NUTS_C_WH2 = "Wasted children (survey)",
  tclim_pdsi_t0 = "Drought index", wdist_perm_km_mean = "Distance to permanent water",
  wdist_coast_km_mean = "Distance to the coast", wdist_coast_km_min = "Distance to the coast (nearest)",
  lcover_crops_frac_t0 = "Cropland cover", wapor_sd_t0 = "Vegetation seasonality",
  grassland_frac = "Grassland share", elevation = "Elevation",
  map_sy_pf_parasite_rate = "Malaria parasite rate", map_sy_pf_incidence_rate = "Malaria incidence",
  dhs_hh_improved_water = "Households with improved water")

plain_of <- function(x) {
  out <- unname(PLAIN[x])
  if (any(is.na(out))) {   # fall back to a cleaned code rather than failing the whole figure set (2026-09-09)
    warning("no plain-language name for: ", paste(x[is.na(out)], collapse = ", "), " (cleaned code used)")
    out[is.na(out)] <- gsub("_", " ", sub("_t0$", "", sub("^(dhs|glw|ihme|map|spam|lcover|wdist|tclim|soil|wpop)_", "", x[is.na(out)])))
  }
  out
}
group_of <- function(x) ifelse(grepl("^dhs_", x), "Household survey (DHS)",
                        ifelse(grepl("^(ihme_|glw_|map_|fprice_|espen_|wpop_)", x),
                               "Modelled or administrative", "Remotely sensed environment"))
GRPCOL <- c("Remotely sensed environment" = PROXY,
            "Household survey (DHS)"      = SURVEY,
            "Modelled or administrative"  = GEO)

OUTLAB <- c(child_vitA = "Children: vitamin A", child_iron = "Children: iron",
            women_vitA = "Women: vitamin A",   women_iron = "Women: iron",
            women_folate = "Women: folate",    women_b12 = "Women: B12")

F3 <- TOP |> group_by(outcome) |> slice_min(rank, n = 5) |> ungroup() |>
  mutate(name = plain_of(column), grp = group_of(column),
         panel = factor(OUTLAB[outcome], levels = unname(OUTLAB)),
         key = paste0(panel, "§", name))
F3 <- F3 |> arrange(panel, beta_std) |> mutate(key = factor(key, levels = key))

p3 <- ggplot(F3, aes(x = beta_std, y = key, fill = grp)) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.6) +
  scale_fill_manual(values = GRPCOL) +
  scale_y_discrete(labels = function(k) sub("^.*?§", "", k)) +
  scale_x_continuous(limits = c(-5.6, 5.6), breaks = seq(-4, 4, 2)) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  labs(
       subtitle = "Bars to the right mean the district ranks worse. Colour = the kind of data.",
       x = "Direction and weight in the model", y = NULL,
       caption = "These weights say where deficiency is. They are not causes and not things to change.") +
  theme_deck(base = 16) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        panel.grid.major.x = element_line(colour = "grey92"),
        panel.spacing.x = unit(1.4, "lines"))
save16x9(p3, "fig3_top_predictors.png", h = 7.0)

# =============================================================================
# FIGURE 4 - domain scatter: share of the model vs cost of dropping
# =============================================================================
DOM <- rd(P2, "index_importance_domains.csv") |>
  filter(scope == "pooled", target == "level") |>
  group_by(domain) |> summarise(share = mean(share, na.rm = TRUE), .groups = "drop")
ABL <- rd(P2, "domain_ablation_loco_summary.csv") |> filter(target == "level") |>
  select(domain, delta_drop)
F4 <- inner_join(DOM, ABL, by = "domain") |> filter(is.finite(share), is.finite(delta_drop))

SHORT <- c("Agricultural production, land use" = "Crops & land use",
           "Ruralness, population density, built environment" = "Population density",
           "Nutrition status (MODELLED SURFACE)" = "Modelled child nutrition",
           "Infant and child morbidity/mortality" = "Child illness",
           "Household assets and characteristics" = "Household assets",
           "Education, employment, SES" = "Education & work",
           "Fertility, reproductive health" = "Reproductive health",
           "Ecosystem productivity/greenness" = "Vegetation",
           "Water and coast proximity" = "Water & coast",
           "Malaria incidence and treatment" = "Malaria",
           "Food fortification and supplementation" = "Fortification",
           "Helminth burden and control" = "Worms",
           "Adult and maternal mortality" = "Adult mortality",
           "Food prices and supply" = "Food prices",
           "Climate and weather" = "Climate", "Soil characteristics" = "Soil",
           "Satellite embedding" = "Satellite imagery", "Livestock density" = "Livestock",
           "Water and sanitation" = "Water & sanitation", "Built environment" = "Land cover",
           "Dietary diversity" = "Diet diversity", "Healthcare access" = "Health access",
           "Adult nutrition" = "Adult nutrition")
F4$lab <- ifelse(F4$domain %in% names(SHORT), SHORT[F4$domain], F4$domain)
# two domains carry near-identical names in the metadata; make them distinct
F4$lab[trimws(F4$domain) == "Malaria"] <- "Malaria (single layer)"
F4$lab[F4$domain == "Malaria incidence and treatment"] <- "Malaria"
F4$key <- ifelse(F4$delta_drop >= 0.008, "carries a new country",
          ifelse(F4$delta_drop <= -0.008, "holds a new country back", "neither"))
F4$col <- ifelse(F4$key == "carries a new country", PROXY,
          ifelse(F4$key == "holds a new country back", SURVEY, OTHER))
xm <- median(F4$share, na.rm = TRUE)

p4 <- ggplot(F4, aes(share, delta_drop)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.6) +
  geom_vline(xintercept = xm, colour = "grey85", linewidth = 0.5, linetype = "22") +
  geom_point(aes(colour = col), size = 5.4, alpha = 0.95) +
  ggrepel::geom_text_repel(aes(label = lab, colour = col), size = 4.5, seed = 1,
                           max.overlaps = 30, box.padding = 0.42, min.segment.length = 0.25,
                           segment.colour = "grey70") +
  scale_colour_identity() +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.14))) +
  labs(
       subtitle = "Right = a big share of the model. Above the line = removing it hurts a new country.",
       x = "Share of the model, inside a surveyed country",
       y = "Accuracy lost in a new country if removed",
       caption = "Soil and climate matter most where there is no survey. The household-survey domains do the opposite.") +
  annotate("text", x = quantile(F4$share, 0.62), y = max(F4$delta_drop) * 0.99,
           label = "carries a new country", hjust = 0.5, size = 5.0,
           colour = PROXY, fontface = "bold") +
  annotate("text", x = quantile(F4$share, 0.90), y = min(F4$delta_drop) * 0.88,
           label = "holds a new country back", hjust = 0.5, size = 5.0,
           colour = SURVEY, fontface = "bold") +
  theme_deck() + theme(panel.grid.major.y = element_line(colour = "grey93"))
save16x9(p4, "fig4_domain_scatter.png", h = 6.5)

# =============================================================================
# FIGURE 5 - learning curve over training countries
# =============================================================================
TC <- rd(P2, "training_country_curve.csv") |>
  filter(arm == "domain_index", target == "level", is.finite(spearman)) |>
  group_by(n_train_countries) |>
  summarise(m = mean(spearman), lo = m - sd(spearman)/sqrt(n()), hi = m + sd(spearman)/sqrt(n()),
            n = n(), .groups = "drop")

p5 <- ggplot(TC, aes(n_train_countries, m)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = CHANCE_D, fill = "grey88") +
  annotate("text", x = 1, y = CHANCE_D - 0.018, label = "chance", hjust = 0, size = 4.6, colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = PROXY, alpha = 0.16) +
  geom_line(colour = PROXY, linewidth = 1.8) +
  geom_point(colour = PROXY, size = 6) +
  geom_text(aes(label = sprintf("%.2f", m)), vjust = -1.25, size = 5.2, colour = INK) +
  scale_x_continuous(breaks = TC$n_train_countries) +
  scale_y_continuous(limits = c(0, max(TC$hi) * 1.22)) +
  labs(
       subtitle = "Accuracy in a country held out of training, by number of training countries.",
       x = "Countries used for training", y = "Ranking accuracy in the new country",
       caption = "Shaded band = standard error across held-out country-outcome combinations.") +
  theme_deck() + theme(panel.grid.major.y = element_line(colour = "grey93"))
save16x9(p5, "fig5_learning_curve.png", h = 6.2)

# =============================================================================
# FIGURE 6 - targeting bars
# =============================================================================
NT <- rd(P2, "nce_targeting_summary.csv") |> filter(estimand == "infill")
gcap <- function(a) { x <- NT$mean_capture[NT$arm == a]; if (!length(x)) NA_real_ else x }
F6 <- data.frame(
  who = c("Perfect knowledge", "Proxy model", "Survey's regional averages", "No information"),
  v   = c(gcap("oracle_ceiling"), gcap("domain_index"), gcap("region_mean_jk"), gcap("null_train_mean")))
F6 <- F6[is.finite(F6$v), ]
F6$who <- factor(F6$who, levels = rev(F6$who))
F6$col <- c("Perfect knowledge" = OTHER, "Proxy model" = PROXY,
            "Survey's regional averages" = SURVEY, "No information" = OTHER)[as.character(F6$who)]

p6 <- ggplot(F6, aes(v, who, fill = col)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * v)), hjust = -0.16, size = 6.4, colour = INK) +
  scale_fill_identity() +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.56), expand = c(0, 0)) +
  labs(
       subtitle = "Share of a country's deficient people in the fifth of districts each method picks.",
       x = "Share of deficient people reached", y = NULL,
       caption = "'Perfect knowledge' is the unreachable best case; 'no information' means using the national average.") +
  theme_deck() + theme(panel.grid.major.x = element_line(colour = "grey92"))
save16x9(p6, "fig6_targeting.png", h = 6.0)

# =============================================================================
# FIGURE 8 - the twenty public layers behind child iron
# =============================================================================
F8 <- TOP |> filter(outcome == "child_iron", rank <= 20) |>
  mutate(name = plain_of(column), grp = group_of(column),
         dir = ifelse(beta_std > 0, "more deficiency", "less deficiency")) |>
  arrange(desc(rank))
F8$name <- factor(F8$name, levels = F8$name)

p8 <- ggplot(F8, aes(x = 1, y = name, fill = grp)) +
  geom_tile(width = 0.9, height = 0.82, colour = "white", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%s  %s", ifelse(beta_std > 0, "\u25b2", "\u25bc"), name)),
            x = 1, size = 4.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = GRPCOL) +
  facet_wrap(~ grp, scales = "free_y", ncol = 3) +
  labs(
       subtitle = "The child-iron composite. Up = the district ranks worse, down = it ranks better.",
       x = NULL, y = NULL,
       caption = sprintf("Direction, not size, is what holds. %s of the twenty also keep their direction in each country's own separate fit.", c("Nine","Ten","Eleven","Twelve","Thirteen","Fourteen","Fifteen","Sixteen","Seventeen","Eighteen","Nineteen","Twenty")[max(1, min(12, sum(TOP$incountry_sign_agree[TOP$outcome == "child_iron" & TOP$rank <= 20] == TOP$incountry_fits[TOP$outcome == "child_iron" & TOP$rank <= 20]) - 8))])) +
  theme_void(base_size = 17) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 17, colour = "grey20", margin = margin(b = 12)),
        plot.caption = element_text(size = 11, colour = "grey45", hjust = 0),
        strip.text = element_text(face = "bold", size = 15, margin = margin(b = 6)),
        legend.position = "none", plot.margin = margin(12, 18, 8, 12))
save16x9(p8, "fig8_twenty_layers.png", h = 6.6)

cat("\nDONE - figures in", OUT, "\n")
