# =============================================================================
# scripts/policy_deck/46_mnf15_v6_maps_methods_appendix.R
#
# Post-fix redraws of three appendix items of the v6 MNF15 talk, and a
# post-fix recompute of one appendix text row. Reads post-fix result tables
# (rebuilt 27-28 September); the only fitting is block "coverage" (the index,
# climate + soil, refit 100 times per held-out country, a few minutes).
# Never writes to an existing file.
#
#   models    a6_model_comparison.png: the method-comparison lollipop of
#             01_figures_main.R figure 1 (three hold-outs, average status),
#             redrawn on the post-fix SuperLearner runs (NS-01 rank shards
#             a-d for in-fill and region, the v6 country run). Every arm that
#             exists in those runs is taken from them (same folds as the SL):
#             the index, the SuperLearner, both neighbour smoothers and the
#             elastic net on the domain components. The survey's regional
#             average comes from benchmarks_v2_cells.csv and "twenty public
#             layers" from weight_sources_cells.csv (both post-fix), over the
#             same country-outcome combinations.
#             -> results/tables/policy_deck/a6_model_comparison_values.csv
#   maps      a6_regional_level_maps.png: child iron, the held-out district
#             predictions (viz/oof_child_iron.csv, rebuilt 28 Sep 05:35 on the
#             post-fix targets) and the survey, both aggregated to Admin-1 with
#             respondent weights, as in the full talk's admin1-pairs chunk.
#             Malawi's Admin-1 units are its 27 districts with data (the
#             project's 87 Malawi units are sub-district areas).
#             -> results/tables/policy_deck/a6_regional_level_values.csv
#   coverage  the "Say how sure each map is?" row: block C of 10_viz_tables.R
#             (coverage of the 90% resampling rank interval under
#             leave-one-country-out, climate + soil, average status) on the
#             post-fix targets_v2.csv.
#             -> results/tables/policy_deck/viz/rank_interval_coverage_postfix.csv
#             -> results/tables/policy_deck/viz/rank_interval_districts_postfix.csv
#             A6_TG=<path> A6_COV_TAG=<tag> reruns it on another targets file
#             (e.g. the pre-fix snapshot, to check the code reproduces the
#             slide's 38% / 55% / 52%).
#
#   Rscript scripts/policy_deck/46_mnf15_v6_maps_methods_appendix.R [models,maps,coverage]
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"; VZ <- file.path(PD, "viz")
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
BLOCKS <- if (length(commandArgs(TRUE))) strsplit(commandArgs(TRUE)[1], ",")[[1]] else c("models", "maps", "coverage")
rd  <- function(...) read.csv(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
no_clobber <- function(f) chk(!file.exists(f) || identical(Sys.getenv("A6_OVERWRITE"), "1"), paste("would overwrite", f))
wcsv <- function(d, f) { no_clobber(f); write.csv(d, f, row.names = FALSE); cat("wrote", f, "\n") }
sv <- function(p, f, w, h) { ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 220, bg = "white"); cat("wrote", file.path(OUT, f), "\n") }

PROXY <- "#0F7B8A"; NAVY <- "#274C77"; WARM <- "#B45309"; OTHER <- "#9AA0A6"; INK <- "#1A1A1A"

# =============================================================================
# models: the method comparison, post-fix
# =============================================================================
if ("models" %in% BLOCKS) {
  cat("\n== models\n")
  fs <- c(list.files(P2, pattern = "^weight_sources_raw_ns01_sl_rank_[a-d][.]csv$", full.names = TRUE),
          file.path(P2, "weight_sources_raw_v6_sl_country.csv"))
  chk(length(fs) == 5 && all(file.exists(fs)), "four NS-01 rank shards and the v6 country run")
  SLc <- bind_rows(lapply(fs, read.csv, stringsAsFactors = FALSE)) |> filter(target == "level", is.finite(spearman)) |>
    group_by(estimand, country, outcome, arm) |> summarise(v = mean(spearman), draws = dplyr::n(), .groups = "drop")   # in-fill: mean over its 3 fold draws
  SL_ARMS <- c(domain_index = "idx", sl_nnls = "sl", sl_lrn_spatial_plus_domain = "sp_dom", sl_lrn_spatial_gam = "sp",
               sl_lrn_domain_pc_enet = "enet")
  cells <- SLc |> filter(arm == "domain_index") |> distinct(estimand, country, outcome)
  chk(all(table(cells$estimand)[c("infill", "region", "country")] == c(18, 18, 22)), "18 / 18 / 22 scored combinations")
  chk(abs(mean(SLc$v[SLc$estimand == "infill" & SLc$arm == "domain_index"]) - 0.389) < 0.005, "in-fill index 0.389 as the main talk's slide 9")

  BC <- rd(P2, "benchmarks_v2_cells.csv") |> filter(target == "level", is.finite(spearman))
  WS <- rd(P2, "weight_sources_cells.csv") |> filter(target == "level", is.finite(spearman))
  other <- bind_rows(
    BC |> filter(arm %in% c("region_mean_jk", "domain_index", "spatial", "spatial_plus_domain", "domain_enet")) |>
      transmute(estimand, country, outcome, arm = paste0("bench_", arm), v = spearman),
    WS |> filter(arm %in% c("sparse20", "domain_index")) |> transmute(estimand, country, outcome, arm = paste0("ws_", arm), v = spearman)) |>
    semi_join(cells, by = c("estimand", "country", "outcome"))
  ALL <- bind_rows(SLc |> select(estimand, country, outcome, arm, v), other)

  CM <- read.csv("results/figures/mnf15/cell_master.csv", stringsAsFactors = FALSE)
  keep <- CM |> filter(keep) |> select(country, outcome)
  summ <- function(d) d |> group_by(estimand, arm) |> summarise(cells = dplyr::n(), mean_spearman = mean(v), cells_positive = sum(v > 0), .groups = "drop")
  VAL <- bind_rows(summ(ALL) |> mutate(set = "all scored"),
                   summ(ALL |> semi_join(keep, by = c("country", "outcome"))) |> mutate(set = "measurable (cell_master keep)"))
  # paired against the index of the same run
  pair <- function(a, ref) { x <- ALL |> filter(arm %in% c(a, ref)) |> select(estimand, country, outcome, arm, v) |>
      tidyr::pivot_wider(names_from = arm, values_from = v) |> filter(is.finite(.data[[a]]), is.finite(.data[[ref]]))
    x |> group_by(estimand) |> summarise(arm = a, ref = ref, cells = dplyr::n(), diff = mean(.data[[a]] - .data[[ref]]), better = sum(.data[[a]] > .data[[ref]]), .groups = "drop") }
  PR <- bind_rows(pair("sl_nnls", "domain_index"), pair("sl_lrn_spatial_plus_domain", "domain_index"), pair("sl_lrn_spatial_gam", "domain_index"),
                  pair("sl_lrn_domain_pc_enet", "domain_index"), pair("bench_region_mean_jk", "bench_domain_index"),
                  pair("bench_spatial_plus_domain", "bench_domain_index"), pair("bench_spatial", "bench_domain_index"),
                  pair("bench_domain_enet", "bench_domain_index"), pair("ws_sparse20", "ws_domain_index"))
  wcsv(bind_rows(VAL |> mutate(kind = "mean"), PR |> transmute(estimand, arm, cells, mean_spearman = diff, cells_positive = better, set = paste("paired minus", ref), kind = "paired")),
       file.path(PD, "a6_model_comparison_values.csv"))
  cat("\n  means over all scored combinations:\n")
  print(as.data.frame(VAL |> filter(set == "all scored") |> select(estimand, arm, mean_spearman) |> tidyr::pivot_wider(names_from = estimand, values_from = mean_spearman)), digits = 3)
  cat("\n  paired:\n"); print(as.data.frame(PR), digits = 3)
  cat("\n  measurable set only:\n")
  print(as.data.frame(VAL |> filter(set != "all scored") |> select(estimand, arm, cells, mean_spearman) |> tidyr::pivot_wider(names_from = estimand, values_from = c(cells, mean_spearman))), digits = 3)

  LBL <- "Domain-PC index (what we propose)"
  rows_def <- c(domain_index = LBL, sl_lrn_spatial_plus_domain = "Neighbour smoother + proxies", sl_lrn_spatial_gam = "Neighbour smoother alone",
                sl_nnls = "Full SuperLearner (12 to 14 methods)", ws_sparse20 = "Twenty public layers",
                sl_lrn_domain_pc_enet = "Penalised regression", bench_region_mean_jk = "Survey's own regional average")
  panels <- c(infill = "Inside a surveyed country", region = "A region held out", country = "A country with no survey")
  S1 <- VAL |> filter(set == "all scored", arm %in% names(rows_def))
  F1 <- tidyr::expand_grid(estimand = names(panels), arm = names(rows_def)) |> left_join(S1, by = c("estimand", "arm")) |>
    mutate(method = rows_def[arm], value = mean_spearman)
  F1$cannot <- !is.finite(F1$value)
  chk(sum(F1$cannot) == 4 && all(F1$arm[F1$cannot] %in% c("sl_lrn_spatial_plus_domain", "sl_lrn_spatial_gam", "bench_region_mean_jk")), "only the survey-based arms are missing, and only outside a surveyed country")
  ord <- F1 |> filter(estimand == "infill") |> arrange(value) |> pull(method)
  ord <- c(setdiff(ord, LBL), LBL)                      # the proposal stays on top
  F1$method   <- factor(F1$method, levels = ord)
  F1$estimand <- factor(panels[F1$estimand], levels = unname(panels))
  F1$col <- ifelse(F1$arm == "domain_index", PROXY, ifelse(F1$arm == "bench_region_mean_jk", WARM, OTHER))
  NUL <- rd(P2, "transport_null_calibration.csv")
  chance <- data.frame(estimand = factor(unname(panels), levels = unname(panels)),
                       x = c(NUL$null_mean_q95[NUL$tier == "admin2"], NUL$null_mean_q95[NUL$tier == "admin1"], NUL$null_mean_q95[NUL$tier == "admin2"]))
  cat(sprintf("  chance bands (95th pct of the permutation null): districts %.3f, regions %.3f\n", chance$x[1], chance$x[2]))
  p1 <- ggplot(F1, aes(x = value, y = method)) +
    geom_rect(data = chance, inherit.aes = FALSE, aes(xmin = -Inf, xmax = x, ymin = -Inf, ymax = Inf), fill = "grey88", alpha = 0.9) +
    geom_segment(data = subset(F1, !cannot), aes(x = 0, xend = value, yend = method, colour = col), linewidth = 1.5) +
    geom_point(data = subset(F1, !cannot), aes(colour = col), size = 6) +
    geom_text(data = subset(F1, !cannot), aes(label = sprintf("%.2f", value)), hjust = -0.45, size = 5.1, colour = INK) +
    geom_text(data = subset(F1, cannot), aes(x = 0.02, label = "needs a survey here"), hjust = 0, size = 4.6, colour = "grey45", fontface = "italic") +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 0.56), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
    facet_wrap(~ estimand) +
    labs(subtitle = "Does the predicted order of districts match the survey's? Grey band = chance.", x = "Ranking accuracy", y = NULL) +
    theme_minimal(base_size = 18) +
    theme(text = element_text(colour = INK), plot.subtitle = element_text(size = 17.6, colour = "grey20", margin = margin(b = 12)),
          axis.title = element_text(size = 15.3), strip.text = element_text(face = "bold", size = 15.5), strip.clip = "off",
          panel.spacing.x = unit(1.6, "lines"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "none",
          plot.margin = margin(12, 18, 8, 12))
  sv(p1, "a6_model_comparison.png", 12.2, 6.4)
}

# =============================================================================
# maps: survey against held-out prediction at the regional level, child iron
# =============================================================================
if ("maps" %in% BLOCKS) {
  cat("\n== maps\n")
  suppressPackageStartupMessages(library(sf)); sf::sf_use_s2(FALSE)
  f_oof <- file.path(VZ, "oof_child_iron.csv")
  chk(file.info(f_oof)$mtime > as.POSIXct("2026-09-27 16:41:00"), "oof_child_iron.csv is post-fix (rebuilt after 27 Sep 16:41)")
  OOF <- rd(f_oof)
  a1 <- OOF |> group_by(country, Admin1) |>
    summarise(Survey = weighted.mean(y_prev, n_raw), Predicted = weighted.mean(pred, n_raw), districts = dplyr::n(), respondents = sum(n_raw), .groups = "drop")
  RHO <- a1 |> group_by(country) |> summarise(regions = dplyr::n(), spearman = cor(Survey, Predicted, method = "spearman"),
                                              pearson = cor(Survey, Predicted), mae_pp = 100 * mean(abs(Survey - Predicted)), .groups = "drop")
  dis <- OOF |> group_by(country) |> summarise(districts = dplyr::n(), spearman_district = cor(y_prev, pred, method = "spearman"), .groups = "drop")
  RHO <- left_join(RHO, dis, by = "country"); print(as.data.frame(RHO), digits = 3)
  wcsv(left_join(a1, RHO |> select(country, spearman_region = spearman), by = "country"), file.path(PD, "a6_regional_level_values.csv"))

  BND <- readRDS("dashboard/data/admin2_boundaries.rds"); LC <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi")
  CL <- c(Gambia = "The Gambia", Ghana = "Ghana", Malawi = "Malawi")
  mt <- theme_void(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, margin = margin(b = 4)),
          strip.text = element_text(size = 12, colour = "grey25", margin = margin(b = 3)),
          legend.position = "bottom", legend.text = element_text(size = 10.5), legend.title = element_text(size = 11, colour = "grey25"),
          plot.margin = margin(4, 6, 4, 6))
  maps <- list()
  for (cn in names(LC)) {
    B <- sf::st_as_sf(BND[[LC[[cn]]]]) |> sf::st_make_valid() |> group_by(Admin1) |> summarise(geometry = sf::st_union(geometry), .groups = "drop")
    chk(nrow(B) == RHO$regions[RHO$country == cn], paste(cn, "every Admin-1 unit has survey data"))
    g <- left_join(B, a1[a1$country == cn, ], by = "Admin1") |>
      tidyr::pivot_longer(c(Survey, Predicted), names_to = "panel", values_to = "value") |>
      mutate(panel = factor(panel, levels = c("Survey", "Predicted")))
    lim <- range(100 * g$value, na.rm = TRUE); lim <- c(floor(lim[1] / 5) * 5, ceiling(lim[2] / 5) * 5)
    maps[[cn]] <- ggplot(sf::st_as_sf(g)) + geom_sf(aes(fill = 100 * value), colour = "white", linewidth = 0.25) +
      facet_wrap(~ panel, nrow = if (cn == "Gambia") 2 else 1) +
      scale_fill_gradientn(colours = c("#F3EEE6", "#E9B77F", WARM, "#6E2F05"), limits = lim, name = "% deficient",
                           labels = function(x) paste0(x, "%"),
                           guide = guide_colourbar(barwidth = 8.5, barheight = 0.55, ticks = FALSE, title.position = "top", title.hjust = 0.5)) +
      labs(title = sprintf("%s: %d regions\nSpearman %.2f", CL[[cn]], nrow(B), RHO$spearman[RHO$country == cn])) + mt
  }
  p3 <- patchwork::wrap_plots(maps, ncol = 3, widths = c(1.2, 1.1, 0.8))
  sv(p3, "a6_regional_level_maps.png", 10.17, 4.6)
}

# =============================================================================
# coverage: the stability row, post-fix (block C of 10_viz_tables.R, unchanged logic)
# =============================================================================
if ("coverage" %in% BLOCKS) {
  cat("\n== coverage\n")
  if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")   # as 10_viz_tables.R
  source("R/protocol_v2.R"); source("R/protocol_v2_weights.R"); source("R/protocol_v2_importance.R")
  set.seed(20260917L)                                                                                     # as 10_viz_tables.R
  f_tg <- Sys.getenv("A6_TG", file.path(P2, "targets_v2.csv")); TAG <- Sys.getenv("A6_COV_TAG", "postfix")
  OUTV <- Sys.getenv("A6_COV_DIR", VZ)
  if (TAG == "postfix") chk(file.info(f_tg)$mtime > as.POSIXct("2026-09-27 16:41:00"), "targets_v2.csv is post-fix")
  TG <- read.csv(f_tg, stringsAsFactors = FALSE)
  S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
  MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
  domain_of <- stats::setNames(MD$domain, MD$column)
  NB <- 100L; t0 <- Sys.time()
  CS <- drop_near_outcome_v2(intersect(MD$column[MD$domain %in% c("Climate and weather", "Soil characteristics")], names(S)), MD)
  ALL4 <- c("Gambia", "Ghana", "Malawi", "SierraLeone"); rows <- list(); drows <- list()
  for (on in c("child_iron", "child_vitA", "women_iron", "women_vitA", "women_folate", "women_b12")) {
    cl <- list()
    for (cn in ALL4) { t <- TG[TG$country == cn & TG$outcome == on & is.finite(TG$y_level), ]
      m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", CS)], by = c("Admin1", "Admin2")); if (nrow(m) < 12) next
      cl[[cn]] <- list(n = nrow(m), y = m$y_level, X = prep_predictors_v2(as.matrix(m[, CS]))) }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE])); ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y)))); yobs <- unlist(lapply(cl, function(z) z$y))
    for (h in names(cl)) { te <- which(ctry == h); pool <- setdiff(names(cl), h); truth <- rank(-yobs[te])
      R <- matrix(NA_real_, length(te), NB)
      for (b in seq_len(NB)) { tr <- unlist(lapply(pool, function(g) { i <- which(ctry == g); i[sample.int(length(i), length(i), replace = TRUE)] }))
        D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        p <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, D, NULL), error = function(e) rep(NA_real_, length(te)))
        if (all(is.finite(p))) R[, b] <- rank(-p) }
      R <- R[, colSums(is.finite(R)) == length(te), drop = FALSE]; if (!ncol(R)) next
      lo <- apply(R, 1, quantile, 0.05); hi <- apply(R, 1, quantile, 0.95); med <- apply(R, 1, median)
      drows[[paste(on, h)]] <- data.frame(outcome = on, heldout = h, n_districts = length(te), truth_rank = truth, rank_med = med, rank_lo = lo, rank_hi = hi,
        score = abs(truth - med) / length(te), stringsAsFactors = FALSE)
      rows[[paste(on, h)]] <- data.frame(outcome = on, heldout = h, n_districts = length(te), refits = ncol(R), coverage_90 = mean(truth >= lo & truth <= hi),
        median_width = median(hi - lo), width_share = median(hi - lo) / length(te), spearman_med = cor(truth, apply(R, 1, median), method = "spearman"), stringsAsFactors = FALSE)
      cat(sprintf("  %-13s %-12s coverage %.2f width %.0f of %d\n", on, h, mean(truth >= lo & truth <= hi), median(hi - lo), length(te))) } }
  COV <- bind_rows(rows); DR <- bind_rows(drows)
  wcsv(COV, file.path(OUTV, sprintf("rank_interval_coverage_%s.csv", TAG)))
  wcsv(DR,  file.path(OUTV, sprintf("rank_interval_districts_%s.csv", TAG)))
  # the slide's numbers, computed as the full talk's qmd does (lines 286-289)
  p3 <- DR$rank_med <= DR$n_districts / 3; t3 <- DR$truth_rank <= DR$n_districts / 3
  cat(sprintf("\n  SLIDE ROW (%s): coverage %.0f%% (range %.0f-%.0f%%, %d combinations; median width %.0f%% of the list); calibrated 90%% half-width %.0f%% of the list (resampling half-width %.0f%%); worst-third call right %.0f%% (chance %.0f%%); top-five call %.0f%%; %d held-out districts; %.1f min\n",
              TAG, 100 * mean(COV$coverage_90), 100 * min(COV$coverage_90), 100 * max(COV$coverage_90), nrow(COV), 100 * median(COV$width_share),
              100 * quantile(DR$score, 0.9), 100 * median((DR$rank_hi - DR$rank_lo) / 2 / DR$n_districts),
              100 * mean(t3[p3]), 100 * mean(t3), 100 * mean(t3[DR$rank_med <= 5]), nrow(DR), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  CM <- read.csv("results/figures/mnf15/cell_master.csv", stringsAsFactors = FALSE) |> filter(keep)
  k <- paste(COV$heldout, COV$outcome) %in% paste(CM$country, CM$outcome); kd <- paste(DR$heldout, DR$outcome) %in% paste(CM$country, CM$outcome)
  p3k <- p3 & kd
  cat(sprintf("  measurable set only (%d combinations): coverage %.0f%%, calibrated half-width %.0f%%, worst-third right %.0f%% (chance %.0f%%)\n",
              sum(k), 100 * mean(COV$coverage_90[k]), 100 * quantile(DR$score[kd], 0.9), 100 * mean(t3[p3k]), 100 * mean(t3[kd])))
}
