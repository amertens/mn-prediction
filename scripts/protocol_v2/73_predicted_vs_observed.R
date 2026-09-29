# =============================================================================
# scripts/protocol_v2/73_predicted_vs_observed.R   [IL-02f]
#
# PREDICTED AGAINST OBSERVED, FOR THE TWO CLEAREST WITHIN-COUNTRY LEVEL CELLS
#   Gambia women's iron deficiency (prevalence), the calibrated PCA index
#   Malawi women's vitamin B12 (district mean log concentration), elastic net
#   on the proxy columns
# Both under the in-fill estimand (a district nobody surveyed, inside a surveyed
# country), rebuilt exactly as 02_run_benchmarks_v2.R builds them: the same
# targets (targets_v2.csv), predictors, build_cell(), the ten
# make_folds_v2("kfold_district") draws, and every in-fill arm run in the same
# order after each draw, so the elastic-net arms consume the same random stream
# and the scores reproduce the benchmark table (checked and printed below). Each
# district's prediction is the mean of its ten out-of-fold predictions.
#
#   Rscript scripts/protocol_v2/73_predicted_vs_observed.R
# -> results/tables/protocol_v2/il02_pred_vs_obs_districts.csv   district-level (as targets_v2.csv)
# -> results/tables/protocol_v2/il02_pred_vs_obs_scores.csv      per-arm scores, mean over draws
# -> results/figures/il02_pred_vs_obs_examples.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
# The headline benchmark (benchmarks_v2_cells.csv) is the no-DHS predictor set (rerun_benchmarks.sh);
# match it unless V2_PREDICTOR_TIERS says otherwise, and check against the table for that set.
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
BENCH <- switch(Sys.getenv("V2_PREDICTOR_TIERS"), "open,survey_public" = "benchmarks_v2_cells.csv",
                "open" = "benchmarks_v2_cells_open.csv", "benchmarks_v2_cells_withdhs.csv")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; FIGDIR <- "results/figures"; REPS <- 10L
# Subtitles and the bottom note are off for the slides (PI, 28 Sep); IL02_FIG_NOTES=1 restores them.
NOTES <- identical(Sys.getenv("IL02_FIG_NOTES", "0"), "1")

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))
# verbatim from 02_run_benchmarks_v2.R, plus the district names in the return value
build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  ncol_eff <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[ncol_eff]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  D  <- domain_representation_v2(Xr, domain_of)
  y_nat <- m[[ycol]]
  y_mod <- if (target == "prev") .v2_logit(y_nat) else y_nat
  list(country = cn, outcome = on, target = target, y_nat = y_nat, y_mod = y_mod, X = Xr, D = D,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = y_nat, target = target),
       w = m[[ncol_eff]], Admin1 = m$Admin1, Admin2 = m$Admin2, n = nrow(m), n_raw = m$n_raw)
}

CELLS <- list(
  list(country = "Gambia", outcome = "women_iron", target = "prev", arm = "domain_index_cal", model = "Calibrated PCA index",
       title = "The Gambia: iron deficiency in women"),
  list(country = "Malawi", outcome = "women_b12", target = "level", arm = "raw_enet", model = "Elastic net on the proxy columns",
       title = "Malawi: vitamin B12 in women"))
KEEP <- c("null_train_mean", "region_mean_jk")

dist_rows <- list(); score_rows <- list()
for (cc in CELLS) {
  cell <- build_cell(cc$country, cc$outcome, cc$target); stopifnot(!is.null(cell))
  scale <- if (cell$target == "prev") "prev" else "level"
  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cell$n, k = 5, rep_id = r)
    for (a in arms_for_estimand_v2("infill")) {          # 02's order, so the random stream matches
      fn <- ARMS_V2[[a]]; pred_mod <- rep(NA_real_, cell$n)
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12 || !length(te)) next
        p <- tryCatch(fn(tr, te, cell$y_mod, cell$X, cell$D, cell$aux), error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred_mod[te] <- p
      }
      pred <- if (cell$target == "prev") .v2_expit(pred_mod) else pred_mod
      if (a == "null_train_mean") {                       # the null on the natural scale, as in 02
        pred <- rep(NA_real_, cell$n)
        for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12 || !length(te)) next; pred[te] <- mean(cell$y_nat[tr]) }
      }
      score_rows[[length(score_rows) + 1]] <- data.frame(country = cc$country, outcome = cc$outcome, target = cc$target, arm = a, rep = r,
                                                         score_v2(cell$y_nat, pred, cell$w, scale = scale))
      if (a %in% c(cc$arm, KEEP))
        dist_rows[[length(dist_rows) + 1]] <- data.frame(country = cc$country, outcome = cc$outcome, target = cc$target, arm = a, rep = r,
                                                         Admin1 = cell$Admin1, Admin2 = cell$Admin2, observed = cell$y_nat, predicted = pred,
                                                         n_eff = cell$w, n_raw = cell$n_raw)
    }
  }
}
SC <- bind_rows(score_rows) |> group_by(country, outcome, target, arm) |>
  summarise(across(c(spearman, pearson, mae, wmae, bias), ~ mean(.x, na.rm = TRUE)), reps = n(), .groups = "drop")
write.csv(SC, file.path(OUTDIR, "il02_pred_vs_obs_scores.csv"), row.names = FALSE)
DR <- bind_rows(dist_rows) |> group_by(country, outcome, target, arm, Admin1, Admin2) |>
  summarise(observed = first(observed), predicted = mean(predicted, na.rm = TRUE), n_eff = first(n_eff), n_raw = first(n_raw), .groups = "drop")
write.csv(DR, file.path(OUTDIR, "il02_pred_vs_obs_districts.csv"), row.names = FALSE)

# the check: these scores must reproduce the benchmark table the other scripts read
cat("predictor tiers:", Sys.getenv("V2_PREDICTOR_TIERS"), "| checked against", BENCH, "\n")
B <- read.csv(file.path(OUTDIR, BENCH)) |> filter(estimand == "infill")
chk <- SC |> inner_join(B, by = c("country", "outcome", "target", "arm"), suffix = c("", "_bench")) |>
  filter(arm %in% c(sapply(CELLS, `[[`, "arm"), KEEP)) |>
  transmute(country, outcome, arm, spearman = round(spearman, 3), spearman_bench = round(spearman_bench, 3),
            wmae = round(wmae, 3), wmae_bench = round(wmae_bench, 3), mae = round(mae, 3), mae_bench = round(mae_bench, 3))
print(as.data.frame(chk), row.names = FALSE)

# ---- figure --------------------------------------------------------------------------------------
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8983"; GRID <- "#e9e8e4"; SURF <- "#fcfcfb"; BLUE <- "#2a78d6"
theme_v <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.background = element_rect(fill = SURF, colour = NA),
        text = element_text(colour = INK), axis.text = element_text(colour = INK2), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.4), plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = INK2, size = 10), legend.position = "bottom", legend.text = element_text(colour = INK2),
        legend.title = element_text(colour = INK2, size = 10))
panel <- function(cc) {
  s <- SC |> filter(country == cc$country, outcome == cc$outcome)
  g <- function(a, col) s[[col]][s$arm == a]
  d <- DR |> filter(country == cc$country, outcome == cc$outcome, arm == cc$arm)
  nul <- DR |> filter(country == cc$country, outcome == cc$outcome, arm == "null_train_mean")
  prev <- cc$target == "prev"
  tr <- if (prev) function(v) 100 * v else function(v) exp(-v)          # level = mean of -log(pmol/L)
  d$x <- tr(d$observed); d$y <- tr(d$predicted); national <- tr(mean(nul$predicted))
  err <- if (prev) function(a) sprintf("%.1f pts (%.1f)", g(a, "mae"), g(a, "wmae")) else
                   function(a) sprintf("%.0f%% (%.0f%%)", 100 * (exp(g(a, "mae")) - 1), 100 * (exp(g(a, "wmae")) - 1))
  lab <- paste0("Average error, district unweighted (weighted by sample)\n",
                cc$model, ": ", err(cc$arm), "\n",
                "Average of surveyed districts in the region: ", err("region_mean_jk"), "\n",
                "National average: ", err("null_train_mean"))
  lim <- range(c(d$x, d$y)); lim <- if (prev) c(max(0, lim[1] - 3), min(100, lim[2] + 3)) else lim * c(0.9, 1.1)
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.5) +
    geom_hline(yintercept = national, colour = MUTED, linetype = "22", linewidth = 0.5) +
    geom_point(aes(size = n_eff), colour = BLUE, alpha = 0.75, stroke = 0) +
    annotate("text", x = lim[1], y = lim[2], label = lab, hjust = 0, vjust = 1, size = 3.1, colour = INK2, lineheight = 1.05) +
    annotate("text", x = lim[2], y = national, label = "national average", hjust = 1, vjust = -0.5, size = 3, colour = MUTED) +
    scale_size_area(max_size = 6, name = "District effective sample size") +
    labs(title = cc$title,
         subtitle = if (NOTES) sprintf("%d districts, each predicted without its own survey data. Rank correlation %.2f.", nrow(d), g(cc$arm, "spearman")),
         x = if (prev) "Survey prevalence in the district (%)" else "Survey geometric mean in the district (pmol/L)",
         y = if (prev) "Predicted prevalence (%)" else "Predicted geometric mean (pmol/L)") +
    theme_v
  if (prev) p + coord_equal(xlim = lim, ylim = lim) + scale_x_continuous(labels = function(v) paste0(v, "%")) + scale_y_continuous(labels = function(v) paste0(v, "%"))
  else p + scale_x_log10() + scale_y_log10() + coord_equal(xlim = lim, ylim = lim)
}
p1 <- panel(CELLS[[1]]); p2 <- panel(CELLS[[2]])
cap <- paste0("Each dot is a district; its prediction is the mean of ten out-of-fold predictions (5 folds blocked by district, 10 draws). Diagonal: prediction equals the survey value.\n",
              "Survey values are themselves noisy (median 24 women per district in The Gambia, 8 in Malawi), so part of every error is survey noise; the large dots are the best-measured districts.")
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(p1, p2, ncol = 2) + patchwork::plot_annotation(caption = if (NOTES) cap, theme = theme(plot.caption = element_text(colour = INK2, size = 9, hjust = 0), plot.background = element_rect(fill = SURF, colour = NA)))
  ggsave(file.path(FIGDIR, "il02_pred_vs_obs_examples.png"), fig, width = 13, height = 7, dpi = 200, bg = SURF)
} else {
  png(file.path(FIGDIR, "il02_pred_vs_obs_examples.png"), width = 13, height = 7, units = "in", res = 200, bg = SURF)
  gridExtra::grid.arrange(p1, p2, ncol = 2, bottom = if (NOTES) grid::textGrob(cap, x = 0.01, hjust = 0, gp = grid::gpar(col = INK2, fontsize = 9)))
  dev.off()
}
cat("figure written\n")
