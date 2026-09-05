# =============================================================================
# scripts/protocol_v2/30_training_curve_climate_soil.R   [TC-02]
#
# DOES EACH ADDED COUNTRY STILL BUY ~0.05 WITH A TWO-DOMAIN INDEX?
#
# Script 15 measured the training-country curve for the full index: transport
# rises ~0.05 per added training country, monotone in all four specifications.
# DA-01/02 found a climate+soil index transports at least as well at the
# district rung. The NCE's "add countries" argument and the "collect the
# remotely sensed layers" argument are made separately; this asks whether they
# hold together -- whether a parsimonious index also improves with countries,
# or whether it is already saturated at three.
#
# Identical to script 15 except the domain set; both sets run on the same
# subsets so the curves are directly comparable.
#
#   Rscript scripts/protocol_v2/30_training_curve_climate_soil.R
# -> results/tables/protocol_v2/training_curve_climate_soil.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; MIN_TRAIN <- 20L; set.seed(20260903L)
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
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
SETS <- list(full = domains, climate_soil = c("Climate and weather", "Soil characteristics"))
build_cell <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; if (!all(c(ycol, wcol) %in% names(t))) return(NULL)
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]; if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  yn <- m[[ycol]]
  list(country = cn, n = nrow(m), y_nat = yn, y_mod = if (target == "prev") .v2_logit(yn) else yn,
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}
rows <- list()
for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
              Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  for (h in names(cl)) { pool <- setdiff(names(cl), h); te <- which(ctry == h)
    for (k in seq_along(pool)) for (sb in utils::combn(pool, k, simplify = FALSE)) {
      tr <- which(ctry %in% sb); if (length(tr) < MIN_TRAIN) next
      Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr); dom <- col_domain(colnames(Dm))
      for (sname in names(SETS)) { keep <- which(dom %in% SETS[[sname]]); if (!length(keep)) next
        p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm[, keep, drop = FALSE], aux), error = function(e) rep(NA_real_, length(te)))
        if (length(p) != length(te)) next
        s <- score_v2(ynat[te], p, wv[te], scale = if (target == "prev") "prev" else "level")
        rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = h, set = sname, n_train_countries = k,
          train_set = paste(sb, collapse = "+"), spearman = s$spearman, stringsAsFactors = FALSE) } } }
  cat("curve done", target, on, "\n")
}
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "training_curve_climate_soil.csv"), row.names = FALSE)
cat("\n===== TC-02: transport vs number of training countries, full index vs climate+soil =====\n")
for (tg in unique(R$target)) for (sn in names(SETS)) { d <- R[R$target == tg & R$set == sn, ]
  s <- d |> group_by(n_train_countries) |> summarise(fits = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
        pct_positive = round(100 * mean(spearman > 0, na.rm = TRUE)), .groups = "drop")
  fit <- lm(spearman ~ n_train_countries, data = d)
  cat(sprintf("\n-- %s / %-12s  slope %+.4f per added country --\n", tg, sn, coef(fit)[2])); print(as.data.frame(s), row.names = FALSE) }
cat("\nDONE\n")
