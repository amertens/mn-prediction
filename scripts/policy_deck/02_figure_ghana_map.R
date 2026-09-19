# =============================================================================
# scripts/policy_deck/02_figure_ghana_map.R
#
# Figure 7: Ghana, survey rank, held-out predicted rank, and the model's prediction
# for every district including the ones the survey never reached.
#
# NOTE ON SCOPE. Every other figure in the deck only reads result tables. This
# one cannot: no committed table holds a per-district predicted rank for Ghana.
# It reuses the held-out map code from
# docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd verbatim in substance -
# a 5-fold, 10-repeat cross-validated fit of the index within Ghana, so every
# district is predicted with itself held out. Nothing is re-estimated for any
# other figure and the pipeline is not run.
#
# PARAMETERISED 2026-09-18. The outcome and target are read from the environment
# so the same code can build the figure for any Ghana cell, the way script 09
# does. Defaults are children's iron on the prevalence target: child vitamin A,
# the previous default, is the one Ghana cell where the survey's own regional
# average beats the index on both targets (level 0.277 vs 0.410, prevalence
# 0.265 vs 0.314), so it is the wrong worked example to put in front of a Ghana
# audience.
#
#   Rscript scripts/policy_deck/02_figure_ghana_map.R
#   FIG7_OUTCOME=child_vitA FIG7_TARGET=level Rscript scripts/policy_deck/02_figure_ghana_map.R
#
# -> results/figures/policy_deck/fig7_ghana_map_<outcome>_<target>.png
#    plus a copy at fig7_ghana_map.png (the name the decks reference)
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

ON  <- Sys.getenv("FIG7_OUTCOME", "child_iron")
TGT <- Sys.getenv("FIG7_TARGET",  "prev")
stopifnot(TGT %in% c("prev", "level"))

OLAB <- c(child_iron = "children's iron", child_vitA = "children's vitamin A",
          child_zinc = "children's zinc", women_iron = "women's iron",
          women_vitA = "women's vitamin A", women_b12 = "women's B12",
          women_folate = "women's folate", women_zinc = "women's zinc")
olab <- if (ON %in% names(OLAB)) OLAB[[ON]] else ON
tlab <- if (TGT == "prev") "districts ranked by the share deficient" else "districts ranked by average status"

ycol <- if (TGT == "prev") "y_prev" else "y_level"
ncol_eff <- if (TGT == "prev") "n_eff" else "n_eff_cont"

t <- TG[TG$country == "Ghana" & TG$outcome == ON &
        is.finite(TG[[ycol]]) & is.finite(TG[[ncol_eff]]), ]
m <- dplyr::inner_join(t, S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)],
                       by = c("Admin1", "Admin2"))
cat("outcome:", ON, "| target:", TGT, "| districts scored:", nrow(m), "\n")

Xr  <- prep_predictors_v2(as.matrix(m[, PREDS]))
# The modelling target: logit for a prevalence, the level as supplied otherwise
# (02_run_benchmarks_v2.R build_cell()). Higher y_level = worse status already.
m$y_obs <- m[[ycol]]
Y   <- if (TGT == "prev") .v2_logit(m$y_obs) else m$y_obs
aux <- list(Admin1 = m$Admin1, y_nat = Y)

# The comparator that matters for this figure is the survey's own regional mean
# with the district held out. It is scored on the SAME folds as the model, using
# the pipeline's own arm (arm_region_mean_jk_v2), so the two numbers in the
# subtitle are exactly the comparison benchmarks_v2_cells.csv makes.
# The hard-coded "0.08 for chance" the subtitle used to carry is the permutation
# null for CROSS-COUNTRY transport and is the wrong baseline for an in-fill map.

pred <- matrix(NA_real_, nrow(m), 10)
predjk <- matrix(NA_real_, nrow(m), 10)
rho_rep <- rho_jk_rep <- numeric(10)
for (r in 1:10) {
  folds <- make_folds_v2("kfold_district", nrow(m), k = 5, rep_id = r)
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f)
    D <- domain_representation_v2(Xr, domain_of, sign_rows = tr)
    pred[te, r]   <- ARMS_V2[["domain_index"]](tr, te, Y, NULL, D, aux)
    predjk[te, r] <- ARMS_V2[["region_mean_jk"]](tr, te, Y, NULL, D, aux)
  }
  # scored per draw and averaged, as benchmarks_v2_cells.csv does
  rho_rep[r]    <- cor(m$y_obs, pred[, r],   method = "spearman", use = "complete.obs")
  rho_jk_rep[r] <- cor(m$y_obs, predjk[, r], method = "spearman", use = "complete.obs")
}
m$pred <- rowMeans(pred, na.rm = TRUE)
rho    <- mean(rho_rep)
rho_jk <- mean(rho_jk_rep)
cat(sprintf("held-out ranking accuracy: %.3f (sd over draws %.3f)  |  survey's regional average: %.3f\n",
            rho, stats::sd(rho_rep), rho_jk))

# Ranks, not levels: the ranking is the product, and the two panels are then on
# one comparable scale. 1 = worst district.
m$r_survey <- rank(-m$y_obs, ties.method = "average")
m$r_model  <- rank(-m$pred,  ties.method = "average")
n <- nrow(m)
m$q_survey <- 100 * (m$r_survey - 0.5) / n
m$q_model  <- 100 * (m$r_model  - 0.5) / n

# Per-district export, so any rank quoted on a slide can be checked. r_jk is the
# rank under the survey's own regional average with the district left out, which
# is the number a programme has for an unsurveyed district today.
m$jk <- rowMeans(predjk, na.rm = TRUE)
m$r_jk <- rank(-m$jk, ties.method = "average")
utils::write.csv(m[, c("Admin1", "Admin2", "n_psu", ycol, "pred", "jk",
                       "r_survey", "r_model", "r_jk")],
                 file.path(OUT, sprintf("fig7_ghana_district_ranks_%s_%s.csv", ON, TGT)),
                 row.names = FALSE)

# Third panel: the deployment case. Rank-normalise over ALL of Ghana's districts,
# orient the domain components on the surveyed rows, fit the index on the surveyed
# rows and predict every district. Surveyed districts are in-sample here (the
# honest check is the middle panel); the unsurveyed ones are what a programme
# would receive.
all_s <- S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)]
Xr_all <- prep_predictors_v2(as.matrix(all_s[, PREDS]))
key_all <- paste(all_s$Admin1, all_s$Admin2); key_m <- paste(m$Admin1, m$Admin2)
tr_all <- match(key_m, key_all); tr_all <- tr_all[is.finite(tr_all)]
Y_all <- rep(NA_real_, nrow(all_s)); Y_all[tr_all] <- Y[match(key_all[tr_all], key_m)]
stopifnot(sum(is.finite(Y_all)) == nrow(m))
D_all <- domain_representation_v2(Xr_all, domain_of, sign_rows = tr_all)
pred_all <- ARMS_V2[["domain_index"]](tr_all, seq_len(nrow(all_s)), Y_all, NULL, D_all, list(Admin1 = all_s$Admin1, y_nat = Y_all))
all_s$q_all <- 100 * (rank(-pred_all, ties.method = "average") - 0.5) / nrow(all_s)
cat(sprintf("all districts predicted: %d (surveyed %d, unsurveyed %d)
", nrow(all_s), length(tr_all), nrow(all_s) - length(tr_all)))
g <- dplyr::left_join(B, m[, c("Admin1", "Admin2", "q_survey", "q_model")], by = c("Admin1", "Admin2")) |>
  dplyr::left_join(all_s[, c("Admin1", "Admin2", "q_all")], by = c("Admin1", "Admin2"))
p1 <- "Survey"
p2 <- "Model, held out"
p3 <- "Model, every district"
long <- rbind(
  data.frame(sf::st_drop_geometry(g)[, c("Admin1", "Admin2")], value = g$q_survey,
             panel = p1, geometry = sf::st_geometry(g)),
  data.frame(sf::st_drop_geometry(g)[, c("Admin1", "Admin2")], value = g$q_model,
             panel = p2, geometry = sf::st_geometry(g)),
  data.frame(sf::st_drop_geometry(g)[, c("Admin1", "Admin2")], value = g$q_all,
             panel = p3, geometry = sf::st_geometry(g)))
long <- sf::st_as_sf(long)
long$panel <- factor(long$panel, levels = c(p1, p2, p3))

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
       subtitle = sprintf(
         "Ghana, %s: %s, each district held out.\nRanking accuracy %.2f — the survey's own regional average %.2f.",
         olab, tlab, rho, rho_jk),
       caption = paste("Darker = worse. Grey: no survey clusters. Right: fitted on all surveyed districts and applied to every district, which is what a programme would receive.",
                       "Neighbouring districts resemble each other, so geography alone fills part of the gap. Adding environmental, dietary and health data predicts district",
                       "deficiency better, and it reaches the districts and countries no survey has visited.", sep = "
")) +
  theme_void(base_size = 18) +
  theme(plot.title    = element_text(face = "bold", size = 22, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 17, colour = "grey20", margin = margin(b = 10)),
        plot.caption  = element_text(size = 12, colour = "grey45", hjust = 0),
        strip.text    = element_text(face = "bold", size = 16, margin = margin(b = 8)),
        legend.position = "bottom",
        plot.margin   = margin(12, 18, 8, 12))

tag  <- sprintf("%s_%s", ON, TGT)
f_tag <- file.path(OUT, sprintf("fig7_ghana_map_%s.png", tag))
ggsave(f_tag, p, width = 14.5, height = 6.4, dpi = 200, bg = "white")
file.copy(f_tag, file.path(OUT, "fig7_ghana_map.png"), overwrite = TRUE)
writeLines(sprintf("%.4f", rho), file.path(OUT, sprintf("fig7_ghana_rho_%s.txt", tag)))
writeLines(sprintf("%.4f", rho), file.path(OUT, "fig7_ghana_rho.txt"))
cat(sprintf("wrote %s and copied to fig7_ghana_map.png\n", basename(f_tag)))
