# =============================================================================
# scripts/policy_deck/08_civ_candidates.R   [CV-01, 2026-09-17]
#
# THE TWO PRE-REGISTERED TRANSPORT CANDIDATES ON COTE D'IVOIRE, SIDE BY SIDE
#
# Scripts 04 and 06 now run for two domain sets (CIV_SET=cs, cs_top5) and six
# outcomes on the 212-column CIV database (script 03b). This script summarises
# them: per outcome, how far the two rankings agree, how firmly each places its
# districts (median 90% rank-interval width, districts with P(worst third)
# >= 0.8), and whether the worst-ranked fifth is the same; and draws the child
# iron comparison map for the deck.
#
#   Rscript scripts/policy_deck/08_civ_candidates.R
# -> results/tables/policy_deck/civ_candidates_summary.csv
#    results/tables/policy_deck/civ_candidates_agreement.csv
#    results/figures/policy_deck/fig13_civ_candidates.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTT <- "results/tables/policy_deck"; FDIR <- "results/figures/policy_deck"
R1 <- read.csv(file.path(OUTT, "civ_climate_soil_ranking.csv"), stringsAsFactors = FALSE) |> mutate(domain_set = "cs")
R2 <- read.csv(file.path(OUTT, "civ_cs_top5_ranking.csv"), stringsAsFactors = FALSE) |> mutate(domain_set = "cs_top5")
U  <- read.csv(file.path(OUTT, "civ_rank_uncertainty_all.csv"), stringsAsFactors = FALSE)
G1 <- read.csv(file.path(OUTT, "civ_transport_guards.csv"), stringsAsFactors = FALSE) |> mutate(domain_set = "cs")
G2 <- read.csv(file.path(OUTT, "civ_transport_guards_cs_top5.csv"), stringsAsFactors = FALSE) |> mutate(domain_set = "cs_top5")
R <- bind_rows(R1, R2); n <- dplyr::n_distinct(R$Admin2); k5 <- ceiling(n / 5); k3 <- ceiling(n / 3)

# ── agreement between the two candidates, per outcome ────────────────────────
W <- R |> select(outcome, Admin1, Admin2, domain_set, rank) |> pivot_wider(names_from = domain_set, values_from = rank)
AG <- W |> group_by(outcome) |> summarise(
  spearman = cor(cs, cs_top5, method = "spearman"),
  worst_fifth_overlap = length(intersect(Admin2[cs <= k5], Admin2[cs_top5 <= k5])) / k5,
  worst_third_overlap = length(intersect(Admin2[cs <= k3], Admin2[cs_top5 <= k3])) / k3, .groups = "drop")
# agreement of each candidate ACROSS outcomes (is it one map or six?)
across <- function(set) { M <- R |> filter(domain_set == set) |> select(outcome, Admin2, rank) |> pivot_wider(names_from = outcome, values_from = rank)
  C <- cor(as.matrix(M[, -1]), method = "spearman"); mean(C[upper.tri(C)]) }
AG$across_outcomes_cs <- across("cs"); AG$across_outcomes_cs_top5 <- across("cs_top5")
write.csv(AG, file.path(OUTT, "civ_candidates_agreement.csv"), row.names = FALSE)

# ── firmness and guards, per outcome x set ───────────────────────────────────
SU <- U |> group_by(outcome, domain_set) |> summarise(n = n(), median_width = median(rank_width), p80 = sum(p_worst3rd >= 0.8), p20 = sum(p_worst3rd <= 0.2),
  worst5_stay = sum(rank_med <= 5 & rank_hi <= k3), .groups = "drop")
GU <- bind_rows(G1, G2) |> group_by(domain_set, arm) |> summarise(guard_mean = mean(spearman, na.rm = TRUE), guard_pos = sum(spearman > 0, na.rm = TRUE), n_common = median(n_common), .groups = "drop")
SU <- SU |> left_join(GU |> filter(grepl("with CIV", arm)) |> select(domain_set, guard_mean, guard_pos, n_common), by = "domain_set")
write.csv(SU, file.path(OUTT, "civ_candidates_summary.csv"), row.names = FALSE)
cat("agreement between candidates:\n"); print(as.data.frame(AG), digits = 3)
cat("\nfirmness (median 90% width, P>=0.8 count) and the with-CIV guard:\n"); print(as.data.frame(SU), digits = 3)
cat("\nguards:\n"); print(as.data.frame(GU), digits = 3)

# ── figure: child iron, both candidates, ranking and P(worst third) ──────────
B <- sf::st_read("data/external_cache/gee_geoms/civ_admin2.gpkg", quiet = TRUE)
ci <- U |> filter(outcome == "child_iron")
g <- B |> left_join(ci |> select(Admin1, Admin2, domain_set, rank_med, p_worst3rd, rank_width), by = c("Admin1", "Admin2"), relationship = "many-to-many")
lab <- c(cs = "Climate + soil (pre-registered)", cs_top5 = "+ anaemia, agriculture, infection (DA-03)")
g$set <- factor(lab[g$domain_set], levels = lab)
top <- ci |> filter(domain_set == "cs", rank_med <= 5) |> pull(Admin2)
p1 <- ggplot(g) + geom_sf(aes(fill = rank_med), colour = "white", linewidth = 0.2) + facet_wrap(~ set) +
  scale_fill_distiller(palette = "PuBuGn", direction = 1, name = "Rank (1 = worst)", trans = "reverse") +
  geom_sf_text(data = g |> filter(Admin2 %in% top, domain_set == "cs"), aes(label = Admin2), size = 3, fontface = "bold") +
  labs(title = "Which districts to reach first (median rank over 400 refits)") + theme_void(base_size = 12) + theme(legend.position = "right", strip.text = element_text(face = "bold"))
p2 <- ggplot(g) + geom_sf(aes(fill = p_worst3rd), colour = "white", linewidth = 0.2) + facet_wrap(~ set) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1), name = "P(worst third)") +
  labs(title = "How sure: probability of being in the worst third") + theme_void(base_size = 12) + theme(legend.position = "right", strip.text = element_text(face = "bold"))
ag <- AG[AG$outcome == "child_iron", ]
cap <- sprintf("Children's iron. The two candidates agree at Spearman %.2f across the 33 districts; %.0f%% of the worst fifth is the same.
Uncertainty: 400 stratified resamples of the training districts, index refitted and Cote d'Ivoire re-ranked each time.",
               ag$spearman, 100 * ag$worst_fifth_overlap)
p <- (p1 / p2) + plot_annotation(caption = cap, theme = theme(plot.caption = element_text(hjust = 0, size = 10, colour = "grey30")))
ggsave(file.path(FDIR, "fig13_civ_candidates.png"), p, width = 11, height = 8.6, dpi = 160, bg = "white")
cat("\nwrote", file.path(FDIR, "fig13_civ_candidates.png"), "\nDONE\n")
