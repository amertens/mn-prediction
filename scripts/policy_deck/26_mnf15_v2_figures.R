# =============================================================================
# scripts/policy_deck/26_mnf15_v2_figures.R
#
# Figures for the v2 MNF15 talk (docs/slides/MNF15-talk-2026-09-v2.qmd), the
# restructure proposed in docs/slides/MNF15-talk-outline-2026-09-27.md. Writes
# to results/figures/mnf15_v2/ only, so the current deck's figures
# (results/figures/mnf15/, policy_deck/) are untouched.
#
# One scale on every accuracy figure: the match score (Spearman between the
# model's district order and the survey's, 0 = guessing, 1 = perfect match),
# the scale the dashboard's start page uses. One rank palette on every ranking
# map (darker = worse). Every number is read from a committed result table or
# recomputed and checked against one; the script stops if a check fails.
#
# Needs results/tables/policy_deck/ghana_heldout_child_iron_level*.csv from
# scripts/policy_deck/25_mnf15_v2_ghana_heldout.R.
#
#   Rscript scripts/policy_deck/26_mnf15_v2_figures.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")   # the benchmark's headline set
source("R/protocol_v2.R")
source("R/admin2_key_hygiene.R")

VER <- Sys.getenv("FIG_VER", "v2")   # v3 (27 Sep, later): plain titles, pairs, named models, no ask line; writes mnf15_v3/
V3 <- VER %in% c("v3", "v4", "v5", "v6")
V4 <- VER %in% c("v4", "v5", "v6")
V5 <- VER %in% c("v5", "v6")   # v6 (28 Sep evening): v5 figures on the post-fix measurability screen (script 20 rerun 05:34); v5 (28 Sep, comments on v4-ANM; post-fix tables): plain map titles, no text baked under the maps   # v4 (27 Sep, comments on the reordered deck): better Ghana example, pairs only, no neighbour-map row
OUT <- file.path("results/figures", paste0("mnf15_", VER)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"; XV <- "results/tables/external_validation"
PROXY <- "#0F7B8A"; WARM <- "#B45309"; NAVY <- "#274C77"; INK <- "#1A1A1A"; GREY <- "grey40"; LGREY <- "grey60"
RANKPAL <- c("#0B4F5A", PROXY, "#8FC7CF", "#E7EFF0")   # 1 = worst district = darkest
sv <- function(p, f, w, h) { ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 220, bg = "white"); cat("wrote", f, "\n") }
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)

CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"))
TG   <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
BND  <- readRDS("dashboard/data/admin2_boundaries.rds")
B_gh <- BND[["ghana"]]
cellv <- function(cn, on, tgt, est, arm) CELL$spearman[CELL$country == cn & CELL$outcome == on &
                                                        CELL$target == tgt & CELL$estimand == est & CELL$arm == arm]

map_theme <- theme_void(base_size = 18) +
  theme(strip.text = element_text(face = "bold", size = 17, margin = margin(b = 6)),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 6)),
        legend.position = "bottom", legend.text = element_text(size = 15))
rank_fill <- function(lab_lo = "worst districts", lab_hi = "best districts")
  scale_fill_gradientn(colours = RANKPAL, limits = c(0, 100), na.value = "grey90", breaks = c(6, 94),
                       labels = c(lab_lo, lab_hi), name = NULL,
                       guide = guide_colourbar(barwidth = 15, barheight = 0.8, ticks = FALSE))

#' The match-score scale as a small dot plot, one row per method, on 0 (guessing) to 1 (perfect)
ruler <- function(pts, title) {
  pts$row <- factor(pts$lab, levels = rev(pts$lab))
  ggplot(pts, aes(v, row)) +
    geom_segment(aes(x = 0, xend = 1, yend = row), colour = "grey88", linewidth = 2.4, lineend = "round") +
    geom_point(aes(colour = col), size = 6.5) +
    geom_text(aes(label = if ("txt" %in% names(pts)) txt else sprintf("%.2f", v), colour = col), hjust = -0.12, nudge_x = 0.02, size = 5.8, fontface = "bold") +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, if ("txt" %in% names(pts)) 1.25 else 1.02), breaks = c(0, 0.5, 1),
                       labels = if (V3) c("0 = chance", "0.5", "1 = same order") else c("0 = guessing", "0.5", "1 = perfect match"),
                       expand = expansion(add = c(0.01, 0.01))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 15.5, colour = INK, hjust = 1),
          axis.text.x = element_text(size = 13, colour = GREY),
          plot.title = element_text(size = 15, colour = GREY, hjust = 0, margin = margin(b = 4)),
          plot.title.position = "plot")
}

#' v4: the same ruler on the share of district pairs in the survey's order (50 = coin toss)
ruler4 <- function(pts, title) {
  pts$row <- factor(pts$lab, levels = rev(pts$lab))
  ggplot(pts, aes(v, row)) +
    geom_segment(aes(x = 0.5, xend = 1, yend = row), colour = "grey88", linewidth = 2.4, lineend = "round") +
    geom_point(aes(colour = col), size = 6.5) +
    geom_text(aes(label = sprintf("%.0f of 100 pairs", 100 * v), colour = col), hjust = -0.12, nudge_x = 0.01, size = 5.8, fontface = "bold") +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0.5, 1.12), breaks = c(0.5, 0.75, 1), labels = c("50 = coin toss", "75", "100 = same order"),
                       expand = expansion(add = c(0.005, 0.01))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 15.5, colour = INK, hjust = 1),
          axis.text.x = element_text(size = 13, colour = GREY),
          plot.title = element_text(size = 15, colour = GREY, hjust = 0, margin = margin(b = 4)), plot.title.position = "plot")
}

# ============================================================================ S2
# A national survey measures regions. Programmes act on districts.
g <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_prev), is.finite(n_eff))
n_surv <- nrow(g); n_all <- nrow(B_gh); n_single <- sum(g$n_psu == 1); kids <- median(g$n_raw[g$n_psu == 1])
chk(n_surv == 75 && n_all == 260 && n_single == 62, "Ghana counts (75 / 260 / 62)")
nat <- weighted.mean(g$y_prev, g$n_eff)
reg <- g |> group_by(Admin1) |> summarise(reg = weighted.mean(y_prev, n_eff), .groups = "drop")
gg <- B_gh |> left_join(reg, by = "Admin1") |> left_join(g |> select(Admin1, Admin2, y_prev), by = c("Admin1", "Admin2"))
mk <- function(v, panel) data.frame(sf::st_drop_geometry(gg)[, c("Admin1", "Admin2")], value = v, panel = panel,
                                    geometry = sf::st_geometry(gg))
P <- c("One national\nnumber", "Regional\naverages", "What the survey\nmeasured")
long <- sf::st_as_sf(rbind(mk(rep(nat, nrow(gg)), P[1]), mk(gg$reg, P[2]), mk(gg$y_prev, P[3])))
long$panel <- factor(long$panel, levels = P)
pm <- ggplot(long) + geom_sf(aes(fill = value * 100), colour = "white", linewidth = 0.12) + facet_wrap(~ panel) +
  scale_fill_gradientn(colours = c("#F3EEE6", "#E9B77F", WARM, "#6E2F05"), na.value = "grey88",
                       name = "% of children iron deficient   (grey: never visited)",
                       guide = guide_colourbar(barwidth = 16, barheight = 0.8, ticks = FALSE, title.position = "top")) +
  map_theme + theme(legend.title = element_text(size = 14, colour = GREY), strip.text = element_text(size = 16, lineheight = 0.95))
callout <- function(big, small) ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  annotate("text", x = 0, y = 0.68, label = big, hjust = 0, size = 13, fontface = "bold", colour = WARM) +
  annotate("text", x = 0, y = 0.22, label = small, hjust = 0, vjust = 0.5, size = 6.1, colour = INK, lineheight = 0.92)
side <- callout(sprintf("%d of %d", n_all - n_surv, n_all), "districts have no survey\ncluster at all") /
  callout(sprintf("%d of %d", n_single, n_surv), "visited districts rest on\none cluster, about a\ndozen children")
ask <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Can free public data rank every district, even in a country with no survey?",
           hjust = 0.5, size = 7.8, fontface = "bold", colour = PROXY)
if (V3) sv((pm | side) + plot_layout(widths = c(3, 1.05)), "v2_gap.png", 13, 5.2) else
  sv(((pm | side) + plot_layout(widths = c(3, 1.05))) / ask + plot_layout(heights = c(5, 0.62)), "v2_gap.png", 13, 5.7)
if (V4) {   # v4: a compact version for the introduction slide, beside Sonja's text
  c1 <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +   # one callout only: less about clusters (Andrew, 27 Sep)
    annotate("text", x = 0.05, y = 0.5, label = sprintf("%d of %d", n_all - n_surv, n_all), hjust = 0, size = 12, fontface = "bold", colour = WARM) +
    annotate("text", x = 0.47, y = 0.5, label = "districts had no survey
measurement of their own", hjust = 0, size = 6.2, colour = INK, lineheight = 0.92)
  sv(pm / c1 + plot_layout(heights = c(3.6, 0.7)), "v4_gap_compact.png", 7.8, 5.2)
}
cat(sprintf("  S2: %d of %d unvisited; %d of %d single-cluster; median %.1f children in a single-cluster district\n",
            n_all - n_surv, n_all, n_single, n_surv, kids))

# ============================================================================ S5
# With a survey: the in-fill benchmark reproduced exactly (02_run_benchmarks_v2.R:
# headline tiers, domain components built once on the surveyed rows, 5-fold x 10
# draws), keeping the per-district predictions.
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
t5 <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_level), is.finite(n_eff_cont))
m5 <- inner_join(t5, S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
X5 <- prep_predictors_v2(as.matrix(m5[, PREDS])); D5 <- domain_representation_v2(X5, domain_of)
cx5 <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(B_gh))))
c5 <- data.frame(Admin1 = B_gh$Admin1, Admin2 = B_gh$Admin2, lon = cx5[, 1], lat = cx5[, 2])
ll5 <- left_join(m5[, c("Admin1", "Admin2")], c5, by = c("Admin1", "Admin2"))
Y5 <- m5$y_level; aux5 <- list(Admin1 = m5$Admin1, y_nat = Y5, lon = ll5$lon, lat = ll5$lat)
pr <- pj <- ps <- matrix(NA_real_, nrow(m5), 10); rr <- rj <- rs <- numeric(10)
for (r in 1:10) {
  folds <- make_folds_v2("kfold_district", nrow(m5), k = 5, rep_id = r)
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f)
    pr[te, r] <- ARMS_V2[["domain_index"]](tr, te, Y5, X5, D5, aux5)
    pj[te, r] <- ARMS_V2[["region_mean_jk"]](tr, te, Y5, X5, D5, aux5)
    ps[te, r] <- ARMS_V2[["spatial"]](tr, te, Y5, X5, D5, aux5)
  }
  rr[r] <- score_v2(Y5, pr[, r], m5$n_eff_cont, scale = "level")$spearman
  rj[r] <- score_v2(Y5, pj[, r], m5$n_eff_cont, scale = "level")$spearman
  rs[r] <- score_v2(Y5, ps[, r], m5$n_eff_cont, scale = "level")$spearman
}
rho_model <- cellv("Ghana", "child_iron", "level", "infill", "domain_index")
rho_jk    <- cellv("Ghana", "child_iron", "level", "infill", "region_mean_jk")
rho_sp    <- cellv("Ghana", "child_iron", "level", "infill", "spatial")
cat(sprintf("  S5 in-fill reproduced: model %.4f (table %.4f), regional average %.4f (table %.4f); neighbour map (table) %.4f\n",
            mean(rr), rho_model, mean(rj), rho_jk, rho_sp))
chk(abs(mean(rr) - rho_model) < 0.002 && abs(mean(rj) - rho_jk) < 0.002, "in-fill does not reproduce benchmarks_v2_cells.csv")
m5$pred <- rowMeans(pr); m5$jk <- rowMeans(pj)
m5$r_survey <- rank(-m5$y_level, ties.method = "average"); m5$r_model <- rank(-m5$pred, ties.method = "average")
m5$r_jk <- rank(-m5$jk, ties.method = "average")
write.csv(m5[, c("Admin1", "Admin2", "n_psu", "y_level", "pred", "jk", "r_survey", "r_model", "r_jk")],
          file.path(PD, "v2_ghana_infill_child_iron_level.csv"), row.names = FALSE)
# every district: fit on the 75 surveyed rows, apply to all 260 (deployment)
all5 <- S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)]
Xa5 <- prep_predictors_v2(as.matrix(all5[, PREDS]))
tr_all <- match(paste(m5$Admin1, m5$Admin2), paste(all5$Admin1, all5$Admin2)); chk(!anyNA(tr_all), "surveyed rows in the all-district table")
Ya5 <- rep(NA_real_, nrow(all5)); Ya5[tr_all] <- Y5
Da5 <- domain_representation_v2(Xa5, domain_of, sign_rows = tr_all)
all5$q_all <- 100 * (rank(-ARMS_V2[["domain_index"]](tr_all, seq_len(nrow(all5)), Ya5, Xa5, Da5, list(y_nat = Ya5)),
                          ties.method = "average") - 0.5) / nrow(all5)
n75 <- nrow(m5)
m5$q_survey <- 100 * (m5$r_survey - 0.5) / n75; m5$q_model <- 100 * (m5$r_model - 0.5) / n75
g5 <- B_gh |> left_join(m5[, c("Admin1", "Admin2", "q_survey", "q_model")], by = c("Admin1", "Admin2")) |>
  left_join(all5[, c("Admin1", "Admin2", "q_all")], by = c("Admin1", "Admin2"))
P5 <- c("What the survey\nmeasured", "Model, each\ndistrict hidden", "Model, all\n260 districts")
L5 <- sf::st_as_sf(rbind(
  data.frame(sf::st_drop_geometry(g5)[, c("Admin1", "Admin2")], value = g5$q_survey, panel = P5[1], geometry = sf::st_geometry(g5)),
  data.frame(sf::st_drop_geometry(g5)[, c("Admin1", "Admin2")], value = g5$q_model,  panel = P5[2], geometry = sf::st_geometry(g5)),
  data.frame(sf::st_drop_geometry(g5)[, c("Admin1", "Admin2")], value = g5$q_all,    panel = P5[3], geometry = sf::st_geometry(g5))))
L5$panel <- factor(L5$panel, levels = P5)
ctr <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(g5))))
g5$X <- ctr[, 1]; g5$Y <- ctr[, 2]
nth <- function(k) { k <- round(k); paste0(k, ifelse(k %% 100 %in% 11:13, "th", c("th", "st", "nd", "rd", rep("th", 6))[k %% 10 + 1])) }
EX <- if (V4) c("Central Gonja", "Prestea-Huni Valley") else c("Binduri", "Atiwa East")   # v4: two districts the regional average misplaces
ex <- sf::st_drop_geometry(g5)[g5$Admin2 %in% EX, c("Admin1", "Admin2", "X", "Y")] |>
  left_join(m5[, c("Admin1", "Admin2", "r_survey", "r_model", "r_jk")], by = c("Admin1", "Admin2"))   # pair key, never the name alone
chk(nrow(ex) == 2, "worked-example districts found")
lab5 <- rbind(data.frame(panel = P5[1], ex[, c("Admin2", "X", "Y")], lab = paste(ex$Admin2, nth(ex$r_survey))),
              data.frame(panel = P5[2], ex[, c("Admin2", "X", "Y")], lab = paste(ex$Admin2, nth(ex$r_model))))
lab5$panel <- factor(lab5$panel, levels = P5)
if (V4) lab5$lab <- sub(" ", "\n", lab5$lab)   # two lines, so the labels stay inside their panel
lab5$dx <- ifelse(lab5$Admin2 == EX[1], -0.30, if (V4) 0.30 else -0.30); lab5$dy <- ifelse(lab5$Admin2 == EX[1], 0.45, if (V4) -0.15 else -0.40)
lab5$hj <- ifelse(lab5$dx < 0, 1, 0)
pm5 <- ggplot(L5) + geom_sf(aes(fill = value), colour = "white", linewidth = 0.12) + facet_wrap(~ panel) +
  geom_point(data = lab5, aes(X, Y), shape = 21, size = 8, stroke = 1.6, colour = WARM, fill = NA) +
  geom_label(data = lab5, aes(X + dx, Y + dy, label = lab, hjust = hj), size = 5.0, colour = WARM, fontface = "bold",
             label.size = 0, fill = scales::alpha("white", 0.85), label.padding = unit(0.1, "lines")) +
  rank_fill() + map_theme + coord_sf(clip = "off") +
  theme(legend.position = "right", strip.text = element_text(size = 16, lineheight = 0.95)) +
  guides(fill = guide_colourbar(barwidth = 0.9, barheight = 9, ticks = FALSE))
concord <- function(o, p) { so <- sign(outer(o, o, "-")); sp <- sign(outer(p, p, "-")); m <- so * sp; u <- upper.tri(m); sum(m[u] > 0) / sum(m[u] != 0) }
pair5 <- c(jk = concord(Y5, rowMeans(pj)), sp = concord(Y5, rowMeans(ps)), model = concord(Y5, rowMeans(pr)))
cat(sprintf("  S5 neighbour map reproduced: %.4f (table %.4f); pairs in the survey's order: regional %.0f%%, neighbour %.0f%%, model %.0f%%
",
            mean(rs), rho_sp, 100 * pair5["jk"], 100 * pair5["sp"], 100 * pair5["model"]))
chk(abs(mean(rs) - rho_sp) < 0.01, "neighbour map does not reproduce the table")
ptxt <- function(v, p) sprintf("%.2f  (%.0f%% of pairs)", v, 100 * p)
rl5 <- ruler(data.frame(lab = c("Regional average", if (V3) "Neighbour map (survey only)" else "Neighbour map (no public data)", if (V3) "Domain-PC index" else "Model"),
                        v = c(rho_jk, rho_sp, rho_model), col = c(LGREY, LGREY, PROXY),
                        txt = if (V3) c(ptxt(rho_jk, pair5["jk"]), ptxt(rho_sp, pair5["sp"]), ptxt(rho_model, pair5["model"])) else sprintf("%.2f", c(rho_jk, rho_sp, rho_model))),
             if (V3) "Agreement with the survey's order of districts, each district hidden in turn (children's iron)" else "Match with the survey's order, each district hidden in turn (children's iron)")
if (V4) rl5 <- ruler4(data.frame(lab = c("Regional average", "Domain-PC index"), v = c(pair5["jk"], pair5["model"]), col = c(LGREY, PROXY)),
                      "District pairs put in the survey's order, each district hidden in turn (children's iron)")
sv(pm5 / ((plot_spacer() | rl5 | plot_spacer()) + plot_layout(widths = c(0.12, 1, 0.12))) + plot_layout(heights = c(4.4, if (V4) 1.0 else 1.35)),
   "v2_ghana_infill.png", 13, 6.0)
cat("  S5 ranks (survey / regional average / model):\n"); print(ex[, c("Admin2", "r_survey", "r_jk", "r_model")])

# ============================================================================ S6
# Without a single Ghanaian blood sample
H6 <- read.csv(file.path(PD, "ghana_heldout_child_iron_level.csv"))
A6 <- read.csv(file.path(PD, "ghana_heldout_child_iron_level_all.csv"))
rho_ho <- cellv("Ghana", "child_iron", "level", "country", "domain_index")
chk(abs(cor(H6$y_obs, H6$pred, method = "spearman") - rho_ho) < 0.002, "held-out predictions reproduce the table")
H6$q_survey <- 100 * (H6$r_survey - 0.5) / nrow(H6)
g6 <- B_gh |> left_join(H6[, c("Admin1", "Admin2", "q_survey")], by = c("Admin1", "Admin2")) |>
  left_join(A6[, c("Admin1", "Admin2", "q_all")], by = c("Admin1", "Admin2"))
P6 <- if (V5) c("Survey Child Iron\nRankings", "Modeled Child Iron\nRankings") else c("What the 2017\nsurvey measured", "Model with no\nGhanaian data")
L6 <- sf::st_as_sf(rbind(
  data.frame(sf::st_drop_geometry(g6)[, c("Admin1", "Admin2")], value = g6$q_survey, panel = P6[1], geometry = sf::st_geometry(g6)),
  data.frame(sf::st_drop_geometry(g6)[, c("Admin1", "Admin2")], value = g6$q_all,    panel = P6[2], geometry = sf::st_geometry(g6))))
L6$panel <- factor(L6$panel, levels = P6)
pm6 <- ggplot(L6) + geom_sf(aes(fill = value), colour = "white", linewidth = 0.12) + facet_wrap(~ panel) +
  rank_fill() + map_theme + theme(strip.text = element_text(size = 15.5, lineheight = 0.95))
# each country held out in turn, averaged over all its transported nutrients (the 22 cells of the 0.30 headline)
LO <- CELL |> filter(estimand == "country", target == "level", arm == "domain_index", is.finite(spearman))
chk(nrow(LO) == 22 && (V5 || abs(mean(LO$spearman) - 0.302) < 0.002), "22 transport cells, mean 0.30")
cat(sprintf("  transport mean over 22 cells: %.3f\n", mean(LO$spearman)))
null95 <- read.csv(file.path(P2, "transport_null_calibration.csv"))
pc6 <- LO |> group_by(country) |> summarise(v = mean(spearman), n = n(), .groups = "drop") |>
  mutate(country = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"))
if (V3) {   # the same pairwise scale as the rest of v3: share of district pairs in the survey's order, all 22 transported cells
  PPc <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(estimand == "country", arm == "domain_index")
  chk(nrow(PPc) == 22, "22 transported cells in v3_percell_pairs.csv")
  pc6 <- PPc |> group_by(country) |> summarise(v = mean(pairs), n = n(), .groups = "drop") |>
    mutate(country = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"))
  cat(sprintf("  S6 v3 country bars (pairs): %s\n", paste(sprintf("%s %.0f%%", pc6$country, 100 * pc6$v), collapse = ", ")))
}
bar6 <- ggplot(pc6, aes(v, reorder(country, v))) +
  geom_col(fill = "#8FC7CF", width = 0.62) +
  geom_vline(xintercept = if (V3) 0.5 else 0.079, linetype = "dashed", colour = GREY) +
  geom_text(aes(label = if (V3) sprintf("%.0f%%", 100 * v) else sprintf("%.2f", v)), hjust = -0.25, size = 5.8) +
  annotate("text", x = if (V3) 0.51 else 0.09, y = 4.55, label = if (V3) "coin toss" else "guessing", hjust = 0, size = 4.6, colour = GREY) +
  scale_x_continuous(limits = if (V3) c(0, 0.95) else c(0, 0.9), breaks = NULL) + scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(title = "Each country held out\nin turn, all nutrients", x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 15.5, colour = INK),
        plot.title = element_text(size = 15, face = "bold", colour = INK))
rl6 <- ruler(data.frame(lab = c("Regional average (needs Ghana's survey)", "Neighbour map (needs Ghana's survey)",
                                "Model trained on The Gambia, Sierra Leone, Malawi"),
                        v = c(rho_jk, rho_sp, rho_ho), col = c(LGREY, LGREY, PROXY),
                        txt = if (V3) c(ptxt(rho_jk, pair5["jk"]), ptxt(rho_sp, pair5["sp"]), ptxt(rho_ho, concord(H6$y_obs, H6$pred))) else sprintf("%.2f", c(rho_jk, rho_sp, rho_ho))),
             if (V3) "Agreement with the survey's order of districts, children's iron in Ghana" else "Match with the survey's order, children's iron in Ghana")
top6 <- ((pm6 + theme(legend.position = "right") + guides(fill = guide_colourbar(barwidth = 0.9, barheight = 8, ticks = FALSE))) |
           bar6) + plot_layout(widths = c(2.3, 1))
if (V5) {   # maps and country bars only: no ruler, no caption (Andrew, 28 Sep)
  sv(top6, "v2_ghana_heldout.png", 13, 5.2)
  cat(sprintf("  S6 v5: pairs, model held out %.3f; regional average %.3f\n", concord(H6$y_obs, H6$pred), pair5["jk"]))
} else if (V4) {   # pairs only, no neighbour-map row, and why the held-out model is no worse here than the in-country one
  gh <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(country == "Ghana", arm == "domain_index")
  gh_in <- mean(gh$pairs[gh$estimand == "infill"]); gh_ho <- mean(gh$pairs[gh$estimand == "country" & gh$outcome %in% gh$outcome[gh$estimand == "infill"]])
  n_tr <- sum(TG$outcome == "child_iron" & TG$country != "Ghana" & is.finite(TG$y_level) & is.finite(TG$n_eff_cont))
  cat(sprintf("  S6 v4: Ghana, all six outcomes: in-country %.1f%%, held out %.1f%%; held-out training districts %d\n", 100 * gh_in, 100 * gh_ho, n_tr))
  rl6 <- ruler4(data.frame(lab = c("Regional average (needs Ghana's survey)", "Model trained on The Gambia, Sierra Leone, Malawi"),
                           v = c(pair5["jk"], concord(H6$y_obs, H6$pred)), col = c(LGREY, PROXY)),
                "District pairs put in the survey's order, children's iron in Ghana")
  why <- sprintf("Why about as good as the model that saw Ghana's survey (68 of 100 pairs)? It learned from %d districts in three countries, not about 60 in Ghana. Over all six of Ghana's outcomes it does a little worse: %.0f pairs in 100 against %.0f.",
                 n_tr, 100 * gh_ho, 100 * gh_in)
  sv((top6 / ((plot_spacer() | rl6 | plot_spacer()) + plot_layout(widths = c(0.05, 1, 0.12))) + plot_layout(heights = c(4.3, 1.0))) +
       plot_annotation(caption = stringr::str_wrap(why, 150), theme = theme(plot.caption = element_text(size = 13.5, colour = INK, hjust = 0))),
     "v2_ghana_heldout.png", 13, 6.2)
} else
sv(top6 / ((plot_spacer() | rl6 | plot_spacer()) + plot_layout(widths = c(0.05, 1, 0.12))) + plot_layout(heights = c(4.3, 1.35)),
   "v2_ghana_heldout.png", 13, 6.0)
cat(sprintf("  S6: held out %.4f; per country: %s\n", rho_ho, paste(sprintf("%s %.2f", pc6$country, pc6$v), collapse = ", ")))

# ============================================================================ S7
# No survey at all: Cote d'Ivoire's ranking, and where the model has been checked
sf::sf_use_s2(FALSE)
U <- read.csv(file.path(PD, "civ_rank_uncertainty.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
Bc <- readRDS("dashboard/data/oos_cote_divoire.rds")$boundaries
gc <- dplyr::left_join(Bc, U[, c("Admin1", "Admin2", "rank_med", "rank_width")], by = admin2_join_by(Bc, U))
chk(sum(is.finite(gc$rank_med)) == 33, "33 CIV districts joined")
gc$pct <- 100 * (gc$rank_med - 0.5) / nrow(U)
cc <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(gc)))); gc$X <- cc[, 1]; gc$Y <- cc[, 2]
top3 <- gc[order(gc$rank_med), ][1:3, ]
pciv <- ggplot(gc) + geom_sf(aes(fill = pct), colour = "white", linewidth = 0.25) +
  ggrepel::geom_label_repel(data = sf::st_drop_geometry(top3), aes(X, Y, label = Admin2), size = 5.2, fontface = "bold",
                            colour = INK, fill = scales::alpha("white", 0.85), label.size = 0, seed = 3,
                            box.padding = 0.6, min.segment.length = 0, segment.colour = "grey30") +
  rank_fill() + labs(title = if (V3) "C\u00f4te d'Ivoire, children's iron (climate-and-soil\nmodel): no Ivorian data used" else "C\u00f4te d'Ivoire, children's iron:\nno Ivorian data used") + map_theme
XT <- read.csv(file.path(XV, "xv_transport.csv"))
XP <- read.csv(file.path(XV, "xv_transport_pooled.csv"))
xv <- XT |> filter(arm == "domain_index", is.finite(spearman),
                   (soil == "isda" & arm_group == "africa") | (soil == "sgrid" & country %in% c("India", "Pakistan")))
n_tests <- nrow(xv); n_pos <- sum(xv$spearman > 0)
chk(n_tests == 41 && (V5 || n_pos == 36), "36 of 41 external tests")
cat(sprintf("  external tests: %d of %d positive\n", n_pos, n_tests))
w <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
TRAIN <- c("Gambia", "Ghana", "Sierra Leone", "Malawi"); CHECK <- c("Zambia", "Ethiopia", "Sudan", "Nigeria", "Pakistan", "India")
w$grp <- dplyr::case_when(w$admin %in% TRAIN ~ "Learned from (4 national surveys)",
                          w$admin %in% CHECK ~ "Checked against WHO survey data (6)",
                          w$admin == "Ivory Coast" & !V5 ~ "Mapped with no survey",
                          TRUE ~ "other")
GRP <- c("Learned from (4 national surveys)" = PROXY, "Checked against WHO survey data (6)" = WARM,
         "Mapped with no survey" = NAVY, "other" = "grey90")
lab <- w[w$admin %in% c(CHECK, if (!V5) "Ivory Coast"), ]
lc <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(lab))))
lab$X <- lc[, 1]; lab$Y <- lc[, 2]; lab$nm <- ifelse(lab$admin == "Ivory Coast", "C\u00f4te d'Ivoire", lab$admin)
pw <- ggplot(w) + geom_sf(aes(fill = grp), colour = "white", linewidth = 0.15) +
  ggrepel::geom_text_repel(data = sf::st_drop_geometry(lab), aes(X, Y, label = nm, colour = grp), size = 4.8, fontface = "bold",
                           seed = 5, bg.colour = "white", bg.r = 0.15, box.padding = 0.35, min.segment.length = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = GRP, breaks = if (V5) names(GRP)[1:2] else names(GRP)[1:3], name = NULL) +
  scale_colour_manual(values = GRP, guide = "none") +
  coord_sf(xlim = c(-19, 97), ylim = c(-35, 38), expand = FALSE) +
  labs(title = if (V3) "Checked against WHO survey results in six countries it never saw" else sprintf("Checked in six countries it never saw: %d of %d tests pointed the right way", n_pos, n_tests),
       caption = if (V5) "Survey results by region, not district." else if (V3) "Regions, not districts. Africa: 12 of 12 tests right way round on average status, 17 of 21 on prevalence. South Asia: 7 of 8." else NULL) +
  map_theme + theme(legend.position = "bottom", legend.text = element_text(size = 13.5),
                    plot.title = element_text(size = 15.5)) + guides(fill = guide_legend(nrow = 2))
sv(pciv + pw + plot_layout(widths = c(1, 1.55)), "v2_civ_checked.png", 13, 5.6)
cat(sprintf("  S7: %d of %d external tests positive; pooled: %s\n", n_pos, n_tests,
            paste(sprintf("%s %.3f", XP$set, XP$mean_rho)[1:5], collapse = "; ")))

# ============================================================================ S8
# Some deficiencies can be mapped from public data; others cannot yet
CM <- read.csv("results/figures/mnf15/cell_master.csv")   # written by 20_mnf15_figures.R: in-fill / transport, level target
nut <- CM |> filter(keep) |> group_by(nutrient, pop) |>
  summarise(infill = mean(infill, na.rm = TRUE), tr = mean(tr, na.rm = TRUE), .groups = "drop")
gv <- function(nu, po, col) nut[[col]][nut$nutrient == nu & nut$pop == po]
xa <- XT |> filter(arm == "domain_index", soil == "isda", arm_group == "africa", is.finite(spearman)) |>
  mutate(nut = sub("^(child|women)_", "", outcome)) |> group_by(nut) |>
  summarise(v = mean(spearman), n = n(), pos = sum(spearman > 0), .groups = "drop")
xg <- function(k) { r <- xa[xa$nut == k, ]; sprintf("%.2f, %d of %d", r$v, r$pos, r$n) }
xn <- function(k) xa$v[xa$nut == k]
band_col <- function(v) { v <- round(v, 2); ifelse(is.na(v), "grey", ifelse(v >= 0.40, "g", ifelse(v >= 0.20, "a", "r"))) }   # band on the printed value
f2 <- function(v) sprintf("%.2f", v)
rows <- list(
  list("Vitamin B12, women", f2(gv("Vitamin B12", "Women", "infill")), gv("Vitamin B12", "Women", "infill"),
       f2(gv("Vitamin B12", "Women", "tr")), gv("Vitamin B12", "Women", "tr"), xg("b12"), xn("b12"), "Yes", "g"),
  list("Iron, women / children", paste(f2(gv("Iron", "Women", "infill")), "/", f2(gv("Iron", "Children", "infill"))),
       mean(c(gv("Iron", "Women", "infill"), gv("Iron", "Children", "infill"))),
       paste(f2(gv("Iron", "Women", "tr")), "/", f2(gv("Iron", "Children", "tr"))),
       mean(c(gv("Iron", "Women", "tr"), gv("Iron", "Children", "tr"))), xg("iron"), xn("iron"), "Yes", "g"),
  list("Vitamin A, women / children", paste(f2(gv("Vitamin A", "Women", "infill")), "/", f2(gv("Vitamin A", "Children", "infill"))),
       mean(c(gv("Vitamin A", "Women", "infill"), gv("Vitamin A", "Children", "infill"))),
       paste(f2(gv("Vitamin A", "Women", "tr")), "/", f2(gv("Vitamin A", "Children", "tr"))),
       mean(c(gv("Vitamin A", "Women", "tr"), gv("Vitamin A", "Children", "tr"))),
       paste0(xg("vitA"), "\nboth clear misses"), xn("vitA"), if (V4) "Not yet" else "Not yet, for\nsupplementation", "a"),
  list("Folate, women", f2(gv("Folate", "Women", "infill")), gv("Folate", "Women", "infill"),
       f2(gv("Folate", "Women", "tr")), gv("Folate", "Women", "tr"), xg("folate"), xn("folate"), "Within a\ncountry only", "a"),
  if (V4) list("Zinc (Malawi only)", "model not better\nthan chance", NA, "not testable (measured\nin 1 country only)", NA, "not testable (measured\nin 1 country only)", NA, "No", "r") else
  list("Zinc (Malawi only)", "no district signal\nin the survey", NA, "not testable", NA, "not testable", NA, "No", "r"))
tb <- do.call(rbind, lapply(seq_along(rows), function(i) { r <- rows[[i]]
  data.frame(row = i, col = 1:5, txt = c(r[[1]], r[[2]], r[[4]], r[[6]], r[[8]]),
             band = c("label", band_col(r[[3]]), band_col(r[[5]]), band_col(r[[7]]), r[[9]])) }))
hdr <- data.frame(row = 0, col = 1:5, txt = c("", "Inside a surveyed\ncountry", "Whole country\nheld out", if (V3) "Four African\ncountries (WHO data)" else "Six external\ncountries", "Use it?"), band = "hdr")
tb <- rbind(hdr, tb)
FILL <- c(label = "white", hdr = "white", g = "#D8EDD5", a = "#FBEBC2", r = "#F4D3D0", grey = "grey93")
W8 <- c(3.3, 2.2, 2.2, 2.6, 2.6); x8 <- cumsum(W8) - W8 / 2
tb$x <- x8[tb$col]; tb$w <- W8[tb$col]
p8 <- ggplot(tb) +
  geom_tile(aes(x, -row, width = w - 0.08, height = ifelse(row == 0, 0.9, 0.9), fill = band), colour = NA) +
  geom_text(aes(ifelse(col == 1, x - w / 2 + 0.1, x), -row, label = txt, hjust = ifelse(col == 1, 0, 0.5),
                fontface = ifelse(row == 0 | col %in% c(1, 5), "bold", "plain")),
            size = ifelse(tb$row == 0, 5.4, 5.6), lineheight = 0.9, colour = INK) +
  scale_fill_manual(values = FILL, guide = "none") +
  annotate("text", x = 0.1, y = -5.95, hjust = 0, vjust = 0.5, size = 4.7, colour = GREY, lineheight = 0.95,
           label = if (V3) "Agreement with the survey's order (0 = chance, 1 = same order): green 0.40 or more, amber 0.20 to 0.40, red below 0.20. Measurable\ncombinations only. External: regions of Zambia, Ethiopia, Sudan and Nigeria in WHO survey data; tests pointing the right way." else
             "Match score with the survey (0 = guessing, 1 = perfect): green 0.40 or more, amber 0.20 to 0.40, red below 0.20.\nExternal: regions of Zambia, Ethiopia, Sudan and Nigeria in WHO survey deposits; tests pointing the right way.") +
  coord_cartesian(xlim = c(0, sum(W8)), ylim = c(-6.35, 0.5), expand = FALSE) + theme_void()
sv(p8, "v2_nutrients.png", 13, 5.4)
cat("  S8 external by nutrient (Africa, iSDA):\n"); print(as.data.frame(xa))

# ============================================================================ S9
# Clusters per district: the ceiling (projected, as figJ)
st <- data.frame(k = factor(c("1 cluster", "2 clusters", "3 clusters"), levels = c("1 cluster", "2 clusters", "3 clusters")),
                 v = c(0.539, 0.663, 0.731))
p9 <- ggplot(st, aes(k, v)) + geom_col(fill = c(WARM, PROXY, PROXY), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", v)), vjust = -0.45, size = 7, fontface = "bold") +
  annotate("text", x = 1, y = 0.13, label = "most\ndistricts\ntoday", colour = "white", size = 5.4, fontface = "bold", lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 0.85), expand = c(0, 0)) +
  labs(title = "Best score any model could reach\nagainst the survey's district numbers", x = "survey clusters per district", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 17, colour = INK),
        plot.title = element_text(size = 17, face = "bold", colour = INK), axis.title.x = element_text(size = 15, colour = GREY))
sv(p9, "v2_ceiling.png", 5.6, 4.9)

# ============================================================================ appendix, corrected versions
# pooling curve: no projected point, and the per-country qualifier on the figure
TC <- read.csv(file.path(P2, "training_curve_climate_soil.csv"))
pc <- TC |> filter(target == "level", set == "full") |> group_by(n = n_train_countries) |> summarise(rho = mean(spearman, na.rm = TRUE), .groups = "drop")
pk <- ggplot(pc, aes(n, rho)) + geom_line(linewidth = 1.4, colour = PROXY) + geom_point(size = 6, colour = PROXY) +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -1.3, size = 5.4) +
  scale_x_continuous(breaks = 1:3, limits = c(0.8, 3.2)) + scale_y_continuous(limits = c(0.1, 0.4)) +
  labs(x = "Number of biomarker surveys the model learned from", y = "Match score in a country held out",
       caption = "Average over held-out countries and nutrients. The gain so far is in The Gambia and Ghana as held-out countries; Malawi and Sierra Leone have not improved.") +
  theme_minimal(base_size = 18) + theme(panel.grid.minor = element_blank(), plot.caption = element_text(size = 13, colour = GREY, hjust = 0))
sv(pk, "v2_pooling_curve.png", 12, 5.2)

# Cote d'Ivoire: the same ranking, and how firmly each district is placed (dark = firm)
wmax <- round(max(gc$rank_width, na.rm = TRUE))
pf <- ggplot(gc) + geom_sf(aes(fill = rank_width), colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(colours = rev(c("#F2F2F2", "#BFC6CC", "#7C8B95", "#3D4A54")), limits = c(1, wmax), breaks = c(1.4, wmax - 0.4),
                       labels = c("firm", sprintf("could move %d places", wmax)), name = NULL,
                       guide = guide_colourbar(barwidth = 15, barheight = 0.8, ticks = FALSE)) +
  labs(title = "How firmly each district is placed\n(stability under retraining, not a confidence interval)") + map_theme
sv((if (V4) pciv + labs(title = "Côte d'Ivoire, children's iron
(climate-and-soil model)") else pciv) + pf, "v2_civ_firmness.png", 12.6, 5.6)

# accuracy by nutrient with the folate caption corrected
nl <- CM |> filter(keep) |> mutate(Nutrient = paste0(nutrient, ", ", tolower(pop))) |> group_by(Nutrient) |>
  summarise(`Inside a surveyed country` = mean(infill, na.rm = TRUE), `Country held out of training` = mean(tr, na.rm = TRUE),
            ceiling = mean(ceiling, na.rm = TRUE), .groups = "drop") |>
  pivot_longer(c(`Inside a surveyed country`, `Country held out of training`), names_to = "what", values_to = "rho") |>
  mutate(what = factor(what, levels = c("Inside a surveyed country", "Country held out of training")))
pb <- ggplot(nl, aes(rho, reorder(Nutrient, rho))) + geom_col(aes(x = ceiling), fill = "grey90", width = .62) +
  geom_point(aes(colour = what), size = 5) +
  geom_text(aes(label = sprintf("%.2f", rho), colour = what, vjust = ifelse(what == "Inside a surveyed country", -1.4, 2.2)), size = 4.6, show.legend = FALSE) +
  scale_colour_manual(values = c(PROXY, WARM), name = NULL) + scale_x_continuous(limits = c(0, .8)) +
  labs(subtitle = "Grey bar: the best score any model could reach against the survey's own district numbers",
       x = "Match score (0 = guessing, 1 = perfect)", y = NULL,
       caption = "Measurable combinations only. Folate was measured by immunoassay in Ghana and Sierra Leone and by microbiologic assay in Malawi: not the same measurement.") +
  theme_minimal(base_size = 18) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                        plot.caption = element_text(size = 12.5, colour = GREY, hjust = 0))
sv(pb, "v2_accuracy_by_nutrient.png", 12.2, 6.0)

# ============================================================================ v3 only
if (V3) {
  # How close could any model get? Ceiling and score on the SAME target (average status). The v2 and
  # MNF15 figures set prevalence-target ceilings against average-status scores, which flatters the model.
  VCl <- read.csv(file.path(P2, "variance_components_ceiling.csv")) |> filter(rung == "admin2", target == "level") |>
    select(country, outcome, ceil_level = ceiling_vc)
  CML <- CM |> filter(keep, is.finite(infill)) |> left_join(VCl, by = c("country", "outcome"))
  cn <- CML |> group_by(nutrient) |> summarise(ceiling = mean(ceil_level), model = mean(infill), .groups = "drop") |>
    mutate(nutrient = factor(nutrient, levels = rev(c("Vitamin B12", "Vitamin A", "Iron", "Folate", "Zinc"))))
  cat(sprintf("  v3 ceiling (average status): %s\n", paste(sprintf("%s %.2f vs model %.2f", cn$nutrient, cn$ceiling, cn$model), collapse = "; ")))
  pc9 <- ggplot(cn, aes(y = nutrient)) +
    geom_col(aes(x = ceiling), fill = "grey85", width = 0.62) +
    geom_point(aes(x = model), colour = PROXY, size = 7) +
    geom_text(aes(x = ceiling, label = sprintf("%.2f", ceiling)), hjust = -0.25, size = 5.6, colour = GREY) +
    geom_text(aes(x = model, label = sprintf("%.2f", model)), vjust = -1.25, size = 5.6, colour = PROXY, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0 = chance", "0.5", "1 = same order"), expand = c(0, 0)) +
    labs(title = "Grey: the best any model could score against\nthe survey's own district numbers. Teal: ours.", x = NULL, y = NULL) +
    theme_minimal(base_size = 18) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 17, colour = INK), axis.text.x = element_text(size = 13, colour = GREY),
          plot.title = element_text(size = 15.5, face = "bold", colour = INK), plot.title.position = "plot")
  sv(pc9, "v3_ceiling_by_nutrient.png", 5.8, 4.9)

  # Accuracy for every nutrient and country: share of district pairs in the survey's order, 95% interval
  PP <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> left_join(CM[, c("country", "outcome", "keep", "nutrient", "pop")], by = c("country", "outcome")) |>
    filter(keep) |>
    mutate(ctry = recode(country, SierraLeone = "Sierra Leone", Gambia = "The Gambia"),
           row = sprintf("%s, %s", ctry, tolower(pop)),
           nutrient = factor(nutrient, levels = c("Vitamin B12", "Vitamin A", "Iron", "Folate")),
           panel = factor(ifelse(estimand == "infill", "Inside a surveyed country\n(each district hidden in turn)", "Whole country held out\n(no survey from that country)"),
                          levels = c("Inside a surveyed country\n(each district hidden in turn)", "Whole country held out\n(no survey from that country)")),
           who = ifelse(arm == "domain_index", "Domain-PC index (95% interval)", "Survey's regional average"))
  ordr <- PP |> filter(arm == "domain_index") |> group_by(nutrient, row) |> summarise(m = mean(pairs), .groups = "drop") |> arrange(nutrient, m)
  PP$row <- factor(PP$row, levels = unique(ordr$row))
  chk(nrow(filter(PP, arm == "domain_index")) == (if (VER == "v6") 29 else 30), "measurable combinations: 14 in-country + 16 held out (v6: 15, zinc has no held-out test)")
  pcell <- ggplot(PP, aes(pairs * 100, row)) +
    geom_vline(xintercept = 50, linetype = "dashed", colour = GREY) +
    geom_errorbarh(data = filter(PP, arm == "domain_index"), aes(xmin = lo * 100, xmax = hi * 100), height = 0, colour = PROXY, linewidth = 1.1) +
    geom_point(aes(shape = who, colour = who, fill = who), size = 4.2, stroke = 1.2) +
    scale_shape_manual(values = c("Domain-PC index (95% interval)" = 21, "Survey's regional average" = 23), name = NULL) +
    scale_colour_manual(values = c("Domain-PC index (95% interval)" = PROXY, "Survey's regional average" = WARM), name = NULL) +
    scale_fill_manual(values = c("Domain-PC index (95% interval)" = PROXY, "Survey's regional average" = "white"), name = NULL) +
    facet_grid(nutrient ~ panel, scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(28, 92), breaks = c(30, 50, 70, 90), labels = c("30%", "50%\ncoin toss", "70%", "90%")) +
    labs(x = "Share of district pairs put in the survey's order", y = NULL) +
    theme_minimal(base_size = 15) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 14), strip.text.x = element_text(face = "bold", size = 14),
          strip.text.y.left = element_text(face = "bold", size = 13, angle = 0, hjust = 1), strip.placement = "outside",
          axis.text.y = element_text(size = 12.5, colour = INK), panel.spacing.y = unit(0.35, "lines"))
  sv(pcell, "v3_percell_pairs.png", 13, 6.4)
  s_in <- PP |> filter(arm == "domain_index", estimand == "infill"); s_rg <- PP |> filter(arm == "region_mean_jk"); s_lc <- PP |> filter(arm == "domain_index", estimand == "country")
  cat(sprintf("  v3 per-cell pairs: in-country %.0f%% (regional average %.0f%%), held out %.0f%%; in-country intervals clear 50%% in %d of %d, held out in %d of %d\n",
              100 * mean(s_in$pairs), 100 * mean(s_rg$pairs), 100 * mean(s_lc$pairs), sum(s_in$lo > 0.5), nrow(s_in), sum(s_lc$lo > 0.5), nrow(s_lc)))

  # accuracy by nutrient against the matching (average-status) ceiling, for the appendix
  nl3 <- CML |> mutate(Nutrient = paste0(nutrient, ", ", tolower(pop))) |> group_by(Nutrient) |>
    summarise(`Inside a surveyed country` = mean(infill), ceiling = mean(ceil_level), .groups = "drop") |>
    left_join(CM |> filter(keep) |> mutate(Nutrient = paste0(nutrient, ", ", tolower(pop))) |> group_by(Nutrient) |>
                summarise(`Country held out of training` = mean(tr, na.rm = TRUE), .groups = "drop"), by = "Nutrient") |>
    pivot_longer(c(`Inside a surveyed country`, `Country held out of training`), names_to = "what", values_to = "rho")
  pb3 <- ggplot(nl3, aes(rho, reorder(Nutrient, rho))) + geom_col(aes(x = ceiling), fill = "grey90", width = .62) +
    geom_point(aes(colour = what), size = 5) +
    geom_text(aes(label = sprintf("%.2f", rho), colour = what, vjust = ifelse(what == "Inside a surveyed country", -1.4, 2.2)), size = 4.6, show.legend = FALSE) +
    scale_colour_manual(values = c(PROXY, WARM), name = NULL) + scale_x_continuous(limits = c(0, 1)) +
    labs(subtitle = "Grey bar: the best any model could score against the survey's own district numbers (average status)",
         x = "Agreement with the survey's order (0 = chance, 1 = same order)", y = NULL,
         caption = "Measurable combinations only. Folate was measured by immunoassay in Ghana and Sierra Leone and by microbiologic assay in Malawi: not the same measurement.") +
    theme_minimal(base_size = 18) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                          plot.caption = element_text(size = 12.5, colour = GREY, hjust = 0))
  sv(pb3, "v2_accuracy_by_nutrient.png", 12.2, 6.0)
}
cat("\nall", VER, "figures written to", OUT, "\n")
