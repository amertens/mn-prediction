# =============================================================================
# scripts/protocol_v2/28_climate_soil_admin1.R   [DA-03]
#
# THE CLIMATE + SOIL INDEX AT THE FIRST SUB-NATIONAL TIER
#
# The headline transport result (12 of 12 combinations, mean Spearman ~0.50 at
# Admin-1, script 16) used all 18 domains. DA-01/DA-02 found at the district
# rung that a two-domain remotely-sensed index (climate + soil) transports
# better than the full index (0.368 vs 0.252 on the biomarker level, positive
# in 22/22). Does the same hold one tier up, where the units are larger and
# the survey target less noisy?
#
# Same harness as script 16 (Admin-1 aggregation, LOCO, within-country
# z-scored outcome, PCs oriented from training countries), three domain sets:
#   full        all domains (the published 12/12)
#   climate_soil  Climate and weather + Soil characteristics
#   csa           + Agricultural production, land use
# and the nested-selected set from DA-02 is NOT re-run here: with 4-27 units
# per country the inner LOCO would be selecting on noise.
#
# CAVEAT stated up front: the two-domain set was identified on these four
# countries (DA-01). This is a consistency check across tiers, not an
# independent validation. The validation is the next country added.
#
#   Rscript scripts/protocol_v2/28_climate_soil_admin1.R
# -> results/tables/protocol_v2/climate_soil_admin1.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; TOPFRAC <- 0.20; set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
POP <- readRDS("dashboard/data/admin2_population.rds"); BND <- readRDS("dashboard/data/admin2_boundaries.rds")
POP$country <- gsub(" ", "", POP$country)  # FIX 2026-09-04: the file spells "Sierra Leone" with a space; without this the join silently dropped the country
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
SETS <- list(full = domains, climate_soil = c("Climate and weather", "Soil characteristics"),
             csa = c("Climate and weather", "Soil characteristics", "Agricultural production, land use"))

build_a1 <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]
  if (!nrow(t)) return(NULL)
  pp <- POP[POP$country == cn, c("Admin2", pop_for(on))]; names(pp)[2] <- "pop"
  t <- left_join(t, pp, by = "Admin2"); t <- t[is.finite(t$y_prev) & is.finite(t$pop) & t$pop > 0, ]; if (!nrow(t)) return(NULL)
  a1 <- t |> group_by(Admin1) |> summarise(y = stats::weighted.mean(.data[[ycol]], .data[[wcol]]), prev = stats::weighted.mean(y_prev, n_eff),
                                           w = sum(.data[[wcol]]), pop = sum(pop), .groups = "drop") |> filter(is.finite(y), is.finite(prev))
  x1 <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
  m <- a1 |> inner_join(x1, by = "Admin1") |> inner_join(c1, by = "Admin1"); m <- m[is.finite(m$lon), ]; if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), y_nat = m$y, prev = m$prev, pop = m$pop, y_mod = if (target == "prev") .v2_logit(m$y) else m$y,
       X = Xr, w = m$w, Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}
capture <- function(prev, pop, s) { ok <- is.finite(prev) & is.finite(pop) & is.finite(s); if (sum(ok) < 3) return(NA_real_)
  b <- prev[ok] * pop[ok]; k <- max(1L, round(TOPFRAC * sum(ok))); sel <- order(s[ok], decreasing = TRUE)[seq_len(k)]; sum(b[sel]) / sum(b) }

rows <- list()
for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_a1(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat))
  prv <- unlist(lapply(cl, function(z) z$prev)); pp <- unlist(lapply(cl, function(z) z$pop)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
              Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  folds <- as.integer(factor(ctry))
  for (sname in names(SETS)) {
    pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
      D <- domain_representation_v2(Xm, domain_of, sign_rows = tr); keep <- which(col_domain(colnames(D)) %in% SETS[[sname]])
      if (!length(keep)) next
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, D[, keep, drop = FALSE], aux), error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) pred[te] <- p }
    for (cn in unique(ctry)) { k <- which(ctry == cn)
      s <- score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")
      rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, country = cn, set = sname, n_units = length(k),
        spearman = s$spearman, capture = capture(prv[k], pp[k], pred[k]), thin = length(k) < 8, stringsAsFactors = FALSE) }
  }
  cat("a1 sets done", target, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "climate_soil_admin1.csv"), row.names = FALSE)
cat("\n===== DA-03: first sub-national tier, LOCO, by domain set =====\n")
print(as.data.frame(R |> group_by(target, set) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
  median_rho = round(median(spearman, na.rm = TRUE), 3), cells_positive = sum(spearman > 0, na.rm = TRUE),
  capture = round(mean(capture, na.rm = TRUE), 3), .groups = "drop") |> arrange(target, desc(mean_rho))), row.names = FALSE)
cat("\n(excluding thin countries, n_units < 8)\n")
print(as.data.frame(R[!R$thin, ] |> group_by(target, set) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
  cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |> arrange(target, desc(mean_rho))), row.names = FALSE)
for (tg in unique(R$target)) { W <- pivot_wider(R[R$target == tg, c("outcome", "country", "set", "spearman")], names_from = set, values_from = spearman)
  for (a in c("climate_soil", "csa")) { d <- W[[a]] - W$full
    cat(sprintf("  %-5s %-12s vs full: better in %2d of %2d | median %+.3f\n", tg, a, sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE))) } }
cat("\nDONE\n")
