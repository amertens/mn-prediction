# =============================================================================
# scripts/protocol_v2/16_admin1_transport.R   [GAP 2 + GAP 3a]
#
# DOES TRANSPORT WORK BETTER AT ADMIN-1, AND DOES IT MOVE BURDEN?
#
# Two questions answered from one set of predictions, because computing
# leave-one-country-out at admin-1 twice would be waste.
#
#   GAP 2  ranking. Transport is established at ADMIN-2 (21 of 22 combinations,
#          mean Spearman 0.281). Regional estimates are what a ministry acts on
#          first, and aggregating to admin-1 averages away target noise, so the
#          admin-1 number should read BETTER. The P3 probe hinted at 0.309.
#          If it holds, the NCE's regional claim gets direct evidence.
#
#   GAP 3a burden. Ranking is not reach. Under transport at admin-2 the burden
#          lift was 1.007 across 12 cells -- no better than random -- which is
#          why transport is currently a ranking claim only. Whether that also
#          holds at admin-1, where the units are larger and fewer, is untested.
#
# PROTOCOL. Identical to the country estimand in 02b_merge_and_loco.R: outcome
# z-scored WITHIN country (a ranking claim), predictors rank-normalised within
# country and pooled on common columns, domain PCs oriented from the TRAINING
# countries only.
#
# BURDEN is population-weighted: burden = prevalence x population for the
# outcome's target group, so a small very-deficient region does not outrank a
# large moderately-deficient one. Capture is the share of national burden inside
# the worst-ranked 20% of regions; lift is capture / 0.20.
#
# HONEST LIMIT, STATED UP FRONT. Admin-1 has 4-27 units per country. The worst-
# ranked 20% of 4 regions is ONE region, so Sierra Leone's burden figure is a
# single draw and is reported but flagged. Ranking with 4 units is likewise
# thin. n_units is carried on every row so no reader has to guess.
#
#   Rscript scripts/protocol_v2/16_admin1_transport.R
# -> results/tables/protocol_v2/admin1_transport.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")

OUTDIR  <- "results/tables/protocol_v2"
TOPFRAC <- 0.20
set.seed(20260903L)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
POP <- readRDS("dashboard/data/admin2_population.rds")
POP$country <- gsub(" ", "", POP$country)  # FIX 2026-09-04: the file spells "Sierra Leone" with a space; without this the join silently dropped the country
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"

#' Aggregate outcome, predictors, centroids and population from admin-2 to
#' admin-1. Outcome is precision-weighted by effective n; population is summed.
build_a1 <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]
  if (!nrow(t)) return(NULL)
  pc <- pop_for(on)
  pp <- POP[POP$country == cn, c("Admin2", pc)]
  names(pp)[2] <- "pop"
  t <- left_join(t, pp, by = "Admin2")
  # prevalence is needed for burden even when the modelled target is the level
  t <- t[is.finite(t$y_prev) & is.finite(t$pop) & t$pop > 0, ]
  if (!nrow(t)) return(NULL)

  a1 <- t |> group_by(Admin1) |>
    summarise(y   = stats::weighted.mean(.data[[ycol]], .data[[wcol]]),
              prev = stats::weighted.mean(y_prev, n_eff),
              w   = sum(.data[[wcol]]), pop = sum(pop), .groups = "drop") |>
    filter(is.finite(y), is.finite(prev))
  x1 <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |>
    summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |>
    summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE),
              .groups = "drop")
  m <- a1 |> inner_join(x1, by = "Admin1") |> inner_join(c1, by = "Admin1")
  m <- m[is.finite(m$lon), ]
  if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), y_nat = m$y, prev = m$prev, pop = m$pop,
       y_mod = if (target == "prev") .v2_logit(m$y) else m$y,
       X = Xr, w = m$w, Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}

capture <- function(prev, pop, score) {
  ok <- is.finite(prev) & is.finite(pop) & is.finite(score)
  if (sum(ok) < 3) return(c(NA, NA))
  prev <- prev[ok]; pop <- pop[ok]; score <- score[ok]
  b <- prev * pop
  k <- max(1, round(TOPFRAC * length(prev)))
  sel <- order(score, decreasing = TRUE)[seq_len(k)]
  cap <- sum(b[sel]) / sum(b)
  c(cap, cap / TOPFRAC)
}

rows <- list()
for (target in c("level", "prev")) {
  for (on in unique(TG$outcome)) {
    cl <- list()
    for (cn in COUNTRIES) {
      z <- tryCatch(build_a1(cn, on, target), error = function(e) NULL)
      if (!is.null(z)) cl[[cn]] <- z
    }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
    if (length(common) < 20) next

    Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
    Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
    ynat <- unlist(lapply(cl, function(z) z$y_nat))
    prv  <- unlist(lapply(cl, function(z) z$prev))
    pp   <- unlist(lapply(cl, function(z) z$pop))
    wv   <- unlist(lapply(cl, function(z) z$w))
    aux  <- list(lon = unlist(lapply(cl, function(z) z$lon)),
                 lat = unlist(lapply(cl, function(z) z$lat)),
                 Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))),
                 y_nat = Y)
    folds <- as.integer(factor(ctry))

    for (a in c("null_train_mean", "domain_index", "domain_enet")) {
      pred <- rep(NA_real_, length(Y))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 12) next
        Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
        p <- tryCatch(ARMS_V2[[a]](tr, te, Y, Xm, Dm, aux),
                      error = function(e) rep(NA_real_, length(te)))
        if (length(p) == length(te)) pred[te] <- p
      }
      for (cn in unique(ctry)) {
        k <- which(ctry == cn)
        s <- score_v2(ynat[k], pred[k], wv[k],
                      scale = if (target == "prev") "prev" else "level")
        cp <- capture(prv[k], pp[k], pred[k])
        rows[[length(rows) + 1L]] <- data.frame(
          country = cn, outcome = on, target = target, arm = a,
          n_units = length(k), spearman = s$spearman, topk = s$topk,
          capture_top20 = cp[1], lift = cp[2],
          thin = length(k) < 8, stringsAsFactors = FALSE)
      }
    }
    cat("a1 transport done", target, on, "\n")
  }
}

R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "admin1_transport.csv"), row.names = FALSE)

cat("\n===== GAP 2: ADMIN-1 TRANSPORT, RANKING =====\n")
for (tg in unique(R$target)) {
  d <- R[R$target == tg & R$arm != "null_train_mean", ]
  s <- d |> group_by(arm) |>
    summarise(cells = dplyr::n(),
              mean_rho = round(mean(spearman, na.rm = TRUE), 3),
              median_rho = round(median(spearman, na.rm = TRUE), 3),
              cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop")
  cat("\n--", tg, "--\n"); print(as.data.frame(s), row.names = FALSE)
  cat("   (excluding thin countries, n_units < 8)\n")
  d2 <- d[!d$thin, ]
  s2 <- d2 |> group_by(arm) |>
    summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
              cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop")
  print(as.data.frame(s2), row.names = FALSE)
}

cat("\n===== GAP 3a: ADMIN-1 TRANSPORT, BURDEN CAPTURED =====\n")
cat("(random selection captures 0.20; lift = capture / 0.20)\n")
b <- R[R$target == "prev", ] |> group_by(arm) |>
  summarise(cells = dplyr::n(),
            mean_capture = round(mean(capture_top20, na.rm = TRUE), 3),
            mean_lift = round(mean(lift, na.rm = TRUE), 3),
            cells_lift_gt1 = sum(lift > 1, na.rm = TRUE), .groups = "drop")
print(as.data.frame(b), row.names = FALSE)
cat("\nper-cell (prev), model arms only:\n")
print(as.data.frame(R[R$target == "prev" & R$arm != "null_train_mean",
                      c("country","outcome","arm","n_units","spearman","lift","thin")]),
      row.names = FALSE)
cat("\nDONE\n")
