# =============================================================================
# scripts/policy_deck/45_mnf15_v6_accuracy_appendix.R
#
# Five appendix figures of the v6 MNF15 talk, redrawn on the post-fix tables
# (survey outcome definitions fixed, every result table rebuilt 27-28 Sep).
# Nothing is fitted here: every number is read from a result table on disk, and
# each source is checked to be newer than the post-fix targets_v2.csv.
#
#   a6_training_curve.png            <- protocol_v2/training_curve_climate_soil.csv (script 30)
#                                       full vocabulary checked identical to training_country_curve.csv (script 15)
#   a6_pred_vs_observed.png          <- policy_deck/viz/oof_child_iron.csv (script 10, block A), checked
#                                       against targets_v2.csv (script 73's own outputs cover other cells)
#   a6_where_rankings_work.png       <- protocol_v2/benchmarks_v2_raw.csv (rep level; mean = benchmarks_v2_cells.csv)
#                                       + results/figures/mnf15/cell_master.csv (prevalence, measurability)
#   a6_mean_concentration_error.png  <- protocol_v2/benchmarks_v2_raw.csv (rmse_sd is not carried in _cells.csv)
#   a6_domain_ablation.png           <- protocol_v2/domain_ablation_loco.csv (script 23)
#
# Old figures: full_talk_extracts/ft_training_curve.png, ft_pred_vs_observed.png and v4-ANM slides
# 21, 26, 42 (drawn by docs/slides/MN-proxy-full-talk-2026-09.qmd chunks curve, pred-obs,
# cells-scatter, cells-level, ablation). Same estimands, targets and summaries; deck colours and labels.
#
#   Rscript scripts/policy_deck/45_mnf15_v6_accuracy_appendix.R
#   A6_LEVEL_ARM=domain_index  draws the concentration figure with the uncalibrated index, which is
#   what the v4 slide showed; the default, domain_index_cal, matches the prevalence sibling
#   (v4_ft_prev_error.png, script 31).
# -> results/figures/mnf15_v6/a6_*.png (220 dpi), numbers for the slide text printed to the console
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2 <- "results/tables/protocol_v2"; VZ <- "results/tables/policy_deck/viz"
PROXY <- "#0F7B8A"; NAVY <- "#274C77"; WARM <- "#B45309"; INK <- "#1A1A1A"; GREY <- "grey40"
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
sv <- function(p, f, w, h) { ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 220, bg = "white"); cat(sprintf("wrote %s (%.2f x %.2f in)\n", f, w, h)) }
options(width = 200)

# ---- post-fix guard: every source must be newer than the rebuilt targets -------------------------
FIX_TIME <- file.mtime(file.path(P2, "targets_v2.csv"))
cat("post-fix targets_v2.csv:", format(FIX_TIME, "%d %b %H:%M"), "\n")
postfix <- function(f) {
  t <- file.mtime(f); chk(!is.na(t) && t > FIX_TIME, paste(f, "is not newer than the post-fix targets"))
  cat(sprintf("  source %-58s %s\n", f, format(t, "%d %b %H:%M")))
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
OL <- c(child_vitA = "children's vitamin A", women_vitA = "women's vitamin A", child_iron = "children's iron", women_iron = "women's iron",
        women_folate = "women's folate", women_b12 = "women's B12", child_zinc = "children's zinc", women_zinc = "women's zinc")
CNL <- c(Gambia = "The Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
cl <- function(cn, on) paste0(CNL[cn], ", ", OL[on])
TLAB <- c(level = "Biomarker level", prev = "Prevalence")
CM <- postfix("results/figures/mnf15/cell_master.csv")
chk(sum(CM$keep) == 16, "16 measurable combinations in cell_master.csv")
KEEP <- paste(CM$country, CM$outcome)[CM$keep]

# =================================================================== 1. training-country curve
cat("\n===== 1. Every combination of training countries =====\n")
TCS <- postfix(file.path(P2, "training_curve_climate_soil.csv"))
TC  <- postfix(file.path(P2, "training_country_curve.csv"))
j <- inner_join(TC |> filter(arm == "domain_index") |> select(target, outcome, heldout, train_set, a = spearman),
                TCS |> filter(set == "full") |> select(target, outcome, heldout, train_set, b = spearman),
                by = c("target", "outcome", "heldout", "train_set"))
chk(nrow(j) == sum(TCS$set == "full") && max(abs(j$a - j$b), na.rm = TRUE) < 1e-12,
    "full vocabulary: training_curve_climate_soil.csv (set full) equals training_country_curve.csv (domain_index)")
SETLAB <- c(climate_soil = "Climate + soil", full = "Full vocabulary")
curve <- function(D) D |> group_by(target, set, n_train_countries) |>
  summarise(fits = n(), m = mean(spearman, na.rm = TRUE), lo = quantile(spearman, 0.25, na.rm = TRUE),
            hi = quantile(spearman, 0.75, na.rm = TRUE), positive = mean(spearman > 0, na.rm = TRUE), .groups = "drop")
d1 <- curve(TCS)
print(as.data.frame(d1 |> mutate(across(c(m, lo, hi, positive), ~ round(.x, 3)))), row.names = FALSE)
sl <- function(D, tg, s) { v <- D$m[D$target == tg & D$set == s]; (v[3] - v[1]) / 2 }
cat(sprintf("gain per added training country (k=1 to 3): level full %+.3f, climate+soil %+.3f; prevalence full %+.3f, climate+soil %+.3f\n",
            sl(d1, "level", "full"), sl(d1, "level", "climate_soil"), sl(d1, "prev", "full"), sl(d1, "prev", "climate_soil")))
cat("outcome-country combinations in the curve:", nrow(distinct(TCS, outcome, heldout)), "(",
    sum(paste(distinct(TCS, outcome, heldout)$heldout, distinct(TCS, outcome, heldout)$outcome) %in% KEEP), "of them measurable )\n")
cat("per held-out country, level, mean by number of training countries:\n")
print(as.data.frame(TCS |> filter(target == "level") |> group_by(set, heldout, n_train_countries) |>
  summarise(v = round(mean(spearman, na.rm = TRUE), 3), .groups = "drop") |> pivot_wider(names_from = n_train_countries, values_from = v)), row.names = FALSE)
cat("measurable combinations only (not drawn):\n")
print(as.data.frame(curve(TCS |> filter(paste(heldout, outcome) %in% KEEP)) |> mutate(across(c(m, lo, hi, positive), ~ round(.x, 3)))), row.names = FALSE)
d1p <- d1 |> mutate(target = factor(TLAB[target], levels = TLAB), set = factor(SETLAB[set], levels = SETLAB))
pd <- position_dodge(width = 0.14)
p1 <- ggplot(d1p, aes(n_train_countries, m, colour = set, group = set)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.9, alpha = 0.45, position = pd) +
  geom_line(linewidth = 1.1, position = pd) + geom_point(size = 3, position = pd) +
  facet_wrap(~ target) + scale_x_continuous(breaks = 1:3) +
  scale_colour_manual(values = c("Climate + soil" = WARM, "Full vocabulary" = NAVY)) +
  labs(x = "Training countries", y = "Transported Spearman (mean, IQR)", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(size = 13, colour = INK))
sv(p1, "a6_training_curve.png", 5.8, 4.4)

# =================================================================== 2. predicted against observed
cat("\n===== 2. Does the model track the survey, district by district? =====\n")
OOF <- postfix(file.path(VZ, "oof_child_iron.csv"))
TG  <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
jt <- OOF |> inner_join(TG |> filter(outcome == "child_iron") |> select(country, Admin1, Admin2, y_tg = y_prev, psu_tg = n_psu, raw_tg = n_raw),
                        by = c("country", "Admin1", "Admin2"))
chk(nrow(jt) == nrow(OOF) && max(abs(jt$y_prev - jt$y_tg)) < 1e-12 && all(jt$n_psu == jt$psu_tg) && all(jt$n_raw == jt$raw_tg),
    "oof_child_iron.csv carries the post-fix survey values, clusters and respondents")
pairs_in_order <- function(o, p) { k <- is.finite(o) & is.finite(p); o <- o[k]; p <- p[k]
  m <- sign(outer(o, o, "-")) * sign(outer(p, p, "-")); u <- upper.tri(m); sum(m[u] > 0) / sum(m[u] != 0) }
st <- OOF |> group_by(country) |>
  summarise(n = n(), rho = cor(y_prev, pred, method = "spearman"), pairs = pairs_in_order(y_prev, pred),
            mae = 100 * mean(abs(y_prev - pred)), several = sum(n_psu > 1), sd_ratio = sd(pred) / sd(y_prev),
            size_vs_error = cor(n_raw, abs(y_prev - pred), method = "spearman"), draw_range = 100 * median(pred_hi - pred_lo), .groups = "drop")
print(as.data.frame(st |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
BC <- postfix(file.path(P2, "benchmarks_v2_cells.csv"))
cat("benchmark (mean over draws of each draw's Spearman), child iron, prevalence, in-fill, index:\n")
print(BC |> filter(outcome == "child_iron", target == "prev", estimand == "infill", arm == "domain_index", country %in% OOF$country) |>
        transmute(country, spearman = round(spearman, 3), mae = round(mae, 1)), row.names = FALSE)
PANELS <- c("Ghana", "Malawi", "The Gambia")
d2 <- OOF |> mutate(Country = factor(CNL[country], levels = PANELS), grade = ifelse(n_psu > 1, "several clusters", "one cluster"))
lab2 <- st |> left_join(OOF |> group_by(country) |> summarise(x0 = min(y_prev), .groups = "drop"), by = "country") |>
  mutate(Country = factor(CNL[country], levels = PANELS),
         lab = sprintf("Spearman %.2f\n%.0f of 100 pairs in order\nAverage error %.1f points\n%d districts", rho, 100 * pairs, mae, n))
p2 <- ggplot(d2, aes(y_prev, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi), width = 0, colour = "grey75") +
  geom_point(aes(shape = grade, size = n_raw), colour = PROXY, alpha = 0.85, stroke = 0.7) +
  geom_text(data = lab2, aes(x = x0, y = Inf, label = lab), hjust = 0, vjust = 1.15, size = 3.2, colour = "grey20", lineheight = 0.95) +
  scale_shape_manual(values = c(`several clusters` = 16, `one cluster` = 1), name = NULL) +
  scale_size_continuous(range = c(1.5, 5), name = "Respondents") +
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.04, 0.3))) +
  facet_wrap(~ Country, scales = "free") +
  labs(x = "Survey estimate of child iron deficiency", y = "Held-out prediction (bar: range over draws)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(size = 13, colour = INK))
sv(p2, "a6_pred_vs_observed.png", 11.5, 4.6)

# =================================================================== 3 and 4: per-combination accuracy (rep level)
BR <- postfix(file.path(P2, "benchmarks_v2_raw.csv"))
chk(identical(unique(BR$predictor_tiers), "open,survey_public"), "headline predictor set (no DHS)")
CELL <- BR |> filter(is.finite(spearman) | is.finite(rmse_sd)) |> group_by(country, outcome, target, estimand, arm) |>
  summarise(reps = n(), across(c(spearman, rmse_sd), list(m = ~ mean(.x, na.rm = TRUE), se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(is.finite(.x))))),
            .groups = "drop")
cc <- CELL |> inner_join(BC, by = c("country", "outcome", "target", "estimand", "arm")) |> filter(is.finite(spearman))
chk(max(abs(cc$spearman_m - cc$spearman)) < 1e-9, "rep-level means reproduce benchmarks_v2_cells.csv")
STRONG <- 0.5

# =================================================================== 3. where the rankings work
cat("\n===== 3. Where do the rankings work? =====\n")
NUT <- c(iron = "Iron", vitA = "Vitamin A", folate = "Folate", b12 = "B12", zinc = "Zinc")
ix <- CELL |> filter(estimand == "infill", arm == "domain_index", target == "level", is.finite(spearman_m)) |>
  left_join(BC |> filter(estimand == "infill", arm == "domain_index", target == "level") |> select(country, outcome, n_areas), by = c("country", "outcome")) |>
  left_join(CM |> select(country, outcome, prev, ceiling, clusters, single, keep), by = c("country", "outcome")) |>
  mutate(cell = cl(country, outcome), nutrient = factor(NUT[sub("^(child|women)_", "", outcome)], levels = NUT),
         measurable = factor(ifelse(keep, "Measurable", "Not measurable"), levels = c("Measurable", "Not measurable")))
chk(nrow(ix) == 18 && all(is.finite(ix$prev)), "18 in-fill combinations with a prevalence")
print(as.data.frame(ix |> arrange(desc(spearman_m)) |>
  transmute(cell, spearman = round(spearman_m, 3), se = round(spearman_se, 3), prev = round(100 * prev, 1), ceiling = round(ceiling, 2),
            clusters = clusters, single = round(single, 2), keep)), row.names = FALSE)
s <- ix |> filter(spearman_m >= STRONG); g <- ix |> filter(country == "Gambia"); o <- ix |> filter(country != "Gambia")
cat(sprintf("at or above %.1f: %d of %d (mean %.2f): %s\n", STRONG, nrow(s), nrow(ix), mean(s$spearman_m), paste(s$cell, collapse = "; ")))
cat(sprintf("The Gambia: %.2f to %.2f, %d of %d at or above the line | Ghana and Malawi: %d of %d\n",
            min(g$spearman_m), max(g$spearman_m), sum(g$spearman_m >= STRONG), nrow(g), sum(o$spearman_m >= STRONG), nrow(o)))
cat(sprintf("measurable: %d of %d in-fill combinations; at or above the line among them: %d\n", sum(ix$keep), nrow(ix), sum(ix$keep & ix$spearman_m >= STRONG)))
cat(sprintf("lowest prevalence among in-fill combinations: %s %.1f%%; combinations under 1%%: %d\n",
            ix$cell[which.min(ix$prev)], 100 * min(ix$prev), sum(ix$prev < 0.01)))
cat(sprintf("zinc: %s (mean %.2f)\n", paste(sprintf("%s %.2f", ix$cell[grepl("zinc", ix$outcome)], ix$spearman_m[grepl("zinc", ix$outcome)]), collapse = ", "),
            mean(ix$spearman_m[grepl("zinc", ix$outcome)])))
mv <- ix |> filter(country == "Malawi", grepl("vitA", outcome))
cat(sprintf("Malawi vitamin A: %s\n", paste(sprintf("%s prev %.1f%% ceiling %.2f Spearman %.2f", mv$cell, 100 * mv$prev, mv$ceiling, mv$spearman_m), collapse = "; ")))
wv <- ix |> filter(outcome == "women_vitA")
cat(sprintf("women's vitamin A: %s\n", paste(sprintf("%s %.1f%% -> %.2f", wv$cell, 100 * wv$prev, wv$spearman_m), collapse = "; ")))
print(as.data.frame(ix |> group_by(country) |> summarise(combinations = n(), mean_spearman = round(mean(spearman_m), 3), clusters_per_district = first(clusters),
                                                         single_cluster_share = round(first(single), 3), .groups = "drop")), row.names = FALSE)
NUTCOL <- c(Iron = NAVY, `Vitamin A` = WARM, B12 = PROXY, Folate = "#7A4E8C", Zinc = "grey50")
NY <- c("Ghana, children's iron" = 0.045, "The Gambia, children's iron" = -0.05, "Malawi, children's iron" = -0.035)   # label placement only
ix$ny <- ifelse(ix$cell %in% names(NY), NY[ix$cell], 0)
p3 <- ggplot(ix, aes(prev, spearman_m)) +
  geom_hline(yintercept = STRONG, linetype = 3, colour = "grey55") + geom_hline(yintercept = 0, colour = "grey55") +
  geom_errorbar(aes(ymin = spearman_m - 1.96 * spearman_se, ymax = spearman_m + 1.96 * spearman_se), width = 0, colour = "grey70") +
  geom_point(aes(size = n_areas, colour = nutrient, shape = measurable), alpha = 0.9, stroke = 1.1) +
  ggrepel::geom_text_repel(aes(label = cell), size = 2.9, colour = "grey20", max.overlaps = 40, box.padding = 0.4, point.padding = 0.6,
                           force = 2, min.segment.length = 0.3, segment.colour = "grey55", seed = 3, nudge_y = ix$ny) +
  scale_x_log10(labels = function(x) paste0(100 * x, "%"), breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7), limits = c(0.012, 1.1)) +
  scale_size_continuous(range = c(2, 6), breaks = c(30, 75, 87), name = "Surveyed districts") +
  scale_colour_manual(values = NUTCOL, name = NULL) +
  scale_shape_manual(values = c(Measurable = 16, `Not measurable` = 1), name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)), shape = guide_legend(order = 2, override.aes = list(size = 3)),
         size = guide_legend(order = 3)) +
  labs(x = "National prevalence of the deficiency (log scale)", y = "Held-out Spearman, biomarker level (in-fill)") +
  theme_minimal(base_size = 12) + theme(legend.position = "right", panel.grid.minor = element_blank())
sv(p3, "a6_where_rankings_work.png", 7.4, 5.2)

# =================================================================== 4. mean concentration error
cat("\n===== 4. How far off is the predicted mean concentration? =====\n")
LEVEL_ARM <- Sys.getenv("A6_LEVEL_ARM", "domain_index_cal")
chk(LEVEL_ARM %in% c("domain_index_cal", "domain_index"), "A6_LEVEL_ARM is domain_index_cal or domain_index")
ER <- stats::setNames(c(if (LEVEL_ARM == "domain_index_cal") "Model (calibrated for level)" else "Model (the index)", "Survey's regional average"),
                      c(LEVEL_ARM, "region_mean_jk"))
ordc <- ix |> arrange(spearman_m)
lv <- CELL |> filter(target == "level", estimand == "infill", arm %in% c("domain_index", "domain_index_cal", "region_mean_jk", "null_train_mean"), is.finite(rmse_sd_m)) |>
  mutate(cell = cl(country, outcome))
w4 <- lv |> select(cell, arm, rmse_sd_m) |> pivot_wider(names_from = arm, values_from = rmse_sd_m) |> mutate(cell = factor(cell, levels = rev(ordc$cell))) |> arrange(cell)
print(as.data.frame(w4 |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
for (a in c("domain_index_cal", "domain_index", "region_mean_jk"))
  cat(sprintf("  %-17s mean %.3f; below 1 in %d of %d; below the regional average in %s\n", a, mean(w4[[a]]), sum(w4[[a]] < 1), nrow(w4),
              if (a == "region_mean_jk") "-" else sprintf("%d of %d", sum(w4[[a]] < w4$region_mean_jk), nrow(w4))))
d4 <- lv |> filter(arm %in% names(ER)) |> mutate(cell = factor(cell, levels = ordc$cell), who = factor(ER[arm], levels = ER))
chk(nrow(d4) == 36, "18 combinations x 2 arms")
p4 <- ggplot(d4, aes(rmse_sd_m, cell)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(data = d4 |> filter(arm == LEVEL_ARM), aes(xmin = rmse_sd_m - 1.96 * rmse_sd_se, xmax = rmse_sd_m + 1.96 * rmse_sd_se),
                width = 0, colour = "grey60", orientation = "y") +
  geom_point(aes(shape = who, colour = who, size = who), stroke = 1.2) +
  scale_shape_manual(values = stats::setNames(c(16, 18), ER), name = NULL) +
  scale_colour_manual(values = stats::setNames(c(PROXY, WARM), ER), name = NULL) +
  scale_size_manual(values = stats::setNames(c(3.2, 3.8), ER), name = NULL) +
  labs(x = "Error in the district mean log concentration, relative to the spread between districts (1 = the national mean does as well)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.text = element_text(size = 13), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey93"),
        axis.title.x = element_text(size = 11.5), axis.text.y = element_text(colour = INK, size = 12))
sv(p4, "a6_mean_concentration_error.png", 11, 5.4)

# =================================================================== 5. domain ablation
cat("\n===== 5. What travels across borders, and what does not? =====\n")
DA <- postfix(file.path(P2, "domain_ablation_loco.csv"))
full <- DA |> filter(variant == "full", domain == "ALL") |> select(target, outcome, heldout, full = spearman)
d5 <- DA |> filter(variant == "drop", domain != "ALL") |> inner_join(full, by = c("target", "outcome", "heldout")) |> mutate(delta = full - spearman) |>
  filter(target == "level") |> group_by(domain) |>
  summarise(n = sum(is.finite(delta)), est = mean(delta, na.rm = TRUE), se = sd(delta, na.rm = TRUE) / sqrt(n), .groups = "drop") |>
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se) |> arrange(desc(est))
chk(all(d5$n == 22) && nrow(d5) == 21, "21 domains, each over 22 outcome-country combinations")
print(as.data.frame(d5 |> mutate(across(c(est, se, lo, hi), ~ round(.x, 4)))), row.names = FALSE)
cat(sprintf("full index transported Spearman (level): %.3f over %d combinations\n", mean(full$full[full$target == "level"]), sum(full$target == "level")))
cat(sprintf("top three by mean loss: %s | at or below zero: %d of %d | interval clear of zero: %s | interval below zero: %s\n",
            paste(head(d5$domain, 3), collapse = ", "), sum(d5$est <= 0), nrow(d5), paste(d5$domain[d5$lo > 0], collapse = ", "),
            paste(d5$domain[d5$hi < 0], collapse = ", ")))
d5p <- d5 |> arrange(est) |> mutate(label = sub(" \\(HCES\\)$", "", domain), label = factor(label, levels = label))
p5 <- ggplot(d5p, aes(est, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, colour = "grey35", orientation = "y") +
  geom_point(size = 3, colour = NAVY) +
  scale_x_continuous(breaks = seq(-0.02, 0.04, 0.02), labels = function(x) formatC(x, format = "f", digits = 2)) +
  labs(x = "Loss in transported Spearman when dropped\n(mean, 95% interval)", y = NULL) +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), axis.text.y = element_text(colour = INK),
                                        axis.title.x = element_text(hjust = 1), plot.margin = margin(6, 14, 6, 6))
sv(p5, "a6_domain_ablation.png", 5.8, 5.6)

# ---- numbers in the old notes of the transport slides (post-fix sources) ---------------------------
BM  <- postfix(file.path(P2, "benchmarks_v2_summary.csv")); BMD <- postfix(file.path(P2, "benchmarks_v2_summary_withdhs.csv"))
ND  <- postfix(file.path(P2, "nested_domain_selection.csv")); CS1 <- postfix(file.path(P2, "climate_soil_admin1.csv")); A1 <- postfix(file.path(P2, "admin1_transport.csv"))
g1 <- function(X) X$mean_spearman[X$estimand == "country" & X$target == "level" & X$arm == "domain_index"]
cat(sprintf("transport, level: index %.3f without the DHS domains, %.3f with them; climate+soil %.3f district, %.3f regional; full vocabulary regional %.3f\n",
            g1(BM), g1(BMD), mean(ND$spearman[ND$arm == "fixed_cs" & ND$target == "level"], na.rm = TRUE),
            mean(CS1$spearman[CS1$set == "climate_soil" & CS1$target == "level"], na.rm = TRUE),
            mean(A1$spearman[A1$arm == "domain_index" & A1$target == "level"], na.rm = TRUE)))
cat("\nDONE\n")
