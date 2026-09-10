# =============================================================================
# scripts/policy_deck/02_figure_ghana_map.R
#
# Figure 7: Ghana, survey rank versus held-out predicted rank.
#
# NOTE ON SCOPE. Every other figure in the deck only reads result tables. This
# one cannot: no committed table holds a per-district predicted rank for Ghana.
# It reuses the held-out map code from
# docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd verbatim in substance -
# a 5-fold, 10-repeat cross-validated fit of the index within Ghana, so every
# district is predicted with itself held out. Nothing is re-estimated for any
# other figure and the pipeline is not run.
#
#   Rscript scripts/policy_deck/02_figure_ghana_map.R
# -> results/figures/policy_deck/fig7_ghana_map.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUT <- "results/figures/policy_deck"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2  <- "results/tables/protocol_v2"
PROXY <- "#0F7B8A"; INK <- "#1A1A1A"
set.seed(20260909L)

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
B <- readRDS("dashboard/data/admin2_boundaries.rds")[["ghana"]]

ON  <- "child_vitA"
LAB <- "Children's vitamin A deficiency, Ghana"

t <- TG[TG$country == "Ghana" & TG$outcome == ON & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
m <- dplyr::inner_join(t, S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)],
                       by = c("Admin1", "Admin2"))
cat("districts scored:", nrow(m), "\n")

Xr  <- prep_predictors_v2(as.matrix(m[, PREDS]))
Y   <- .v2_logit(m$y_prev)
aux <- list(Admin1 = m$Admin1, y_nat = Y)

pred <- matrix(NA_real_, nrow(m), 10)
for (r in 1:10) {
  folds <- make_folds_v2("kfold_district", nrow(m), k = 5, rep_id = r)
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f)
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
    pred[te, r] <- ARMS_V2[["domain_index"]](tr, te, Y, NULL, D, aux)
  }
}
m$pred <- rowMeans(pred, na.rm = TRUE)
rho <- cor(m$y_prev, m$pred, method = "spearman", use = "complete.obs")
cat(sprintf("held-out ranking accuracy: %.2f\n", rho))

# Ranks, not levels: the ranking is the product, and the two panels are then on
# one comparable scale. 1 = worst district.
m$r_survey <- rank(-m$y_prev, ties.method = "average")
m$r_model  <- rank(-m$pred,   ties.method = "average")
n <- nrow(m)
m$q_survey <- 100 * (m$r_survey - 0.5) / n
m$q_model  <- 100 * (m$r_model  - 0.5) / n

g <- dplyr::left_join(B, m[, c("Admin1", "Admin2", "q_survey", "q_model")],
                      by = c("Admin1", "Admin2"))
p1 <- "What the survey measured"
p2 <- "What the model predicted"
long <- rbind(
  data.frame(sf::st_drop_geometry(g)[, c("Admin1", "Admin2")], value = g$q_survey,
             panel = p1, geometry = sf::st_geometry(g)),
  data.frame(sf::st_drop_geometry(g)[, c("Admin1", "Admin2")], value = g$q_model,
             panel = p2, geometry = sf::st_geometry(g)))
long <- sf::st_as_sf(long)
long$panel <- factor(long$panel, levels = c(p1, p2))

p <- ggplot(long) +
  geom_sf(aes(fill = value), colour = "white", linewidth = 0.18) +
  facet_wrap(~ panel) +
  scale_fill_gradientn(
    colours = c("#0B4F5A", PROXY, "#8FC7CF", "#E7EFF0"),
    limits = c(0, 100), na.value = "grey92",
    breaks = c(5, 95), labels = c("worst districts", "best districts"),
    name = NULL, guide = guide_colourbar(barwidth = 16, barheight = 0.9,
                                         ticks = FALSE, title.position = "top")) +
  labs(
       subtitle = sprintf("Ghana, children's vitamin A, each district held out. Accuracy %.2f against %.2f for chance.", rho, 0.08),
       caption = "Darker = worse. Grey districts had no survey clusters. The north stands out in both maps; district-by-district within it, the model and the survey disagree.") +
  theme_void(base_size = 18) +
  theme(plot.title    = element_text(face = "bold", size = 22, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 17, colour = "grey20", margin = margin(b = 10)),
        plot.caption  = element_text(size = 12, colour = "grey45", hjust = 0),
        strip.text    = element_text(face = "bold", size = 15, margin = margin(b = 8)),
        legend.position = "bottom",
        plot.margin   = margin(12, 18, 8, 12))

ggsave(file.path(OUT, "fig7_ghana_map.png"), p, width = 12.2, height = 6.6, dpi = 200, bg = "white")
cat("wrote fig7_ghana_map.png\n")
writeLines(sprintf("%.4f", rho), file.path(OUT, "fig7_ghana_rho.txt"))
