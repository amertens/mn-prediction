# =============================================================================
# scripts/protocol_v2/33_transport_null_calibration.R   [NC-01]
#
# HOW SURPRISING IS "12 OF 12" UNDER NO TRANSPORT AT ALL?
#
# The headline regional result is 12 of 12 country-outcome cells with positive
# leave-one-country-out Spearman. Read as twelve independent coin flips that is
# p = 0.0002; but the cells are not independent -- several outcomes share a
# country, and outcomes within a country are correlated (iron with iron,
# vitamin A with wasting). The honest question is: if the covariates carried NO
# cross-country information, how often would the SAME procedure return 12/12?
#
# Null construction. Keep everything -- predictors, folds, PC orientation,
# aggregation -- and permute the held-out country's OUTCOME across its own
# units before scoring, independently per outcome but with the SAME permutation
# of units applied to every outcome of that country. That preserves within-
# country correlation between outcomes (the reason cells are not independent)
# while destroying any relation between covariates and outcome. Repeat 500
# times; record the distribution of (a) cells positive and (b) mean Spearman,
# at both the regional tier and the district rung.
#
# The predictions do not change under the null (the model never sees the
# held-out outcome), so this is cheap: fit once per fold, permute the truth.
#
#   Rscript scripts/protocol_v2/33_transport_null_calibration.R
# -> results/tables/protocol_v2/transport_null_calibration.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; NPERM <- 500L; set.seed(20260903L)
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

build <- function(cn, on, tier) {
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont) & t$n_eff_cont > 0, ]
  if (nrow(t) < 12) return(NULL)
  if (tier == "admin1") {
    a <- t |> group_by(Admin1) |> summarise(y = stats::weighted.mean(y_level, n_eff_cont), w = sum(n_eff_cont), .groups = "drop")
    x <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1") |> inner_join(c1, by = "Admin1"); unit <- m$Admin1
  } else {
    m <- t[, c("Admin1", "Admin2", "y_level", "n_eff_cont")]; names(m)[3:4] <- c("y", "w")
    m <- m |> inner_join(S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
      inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2")); unit <- m$Admin2
  }
  m <- m[is.finite(m$lon) & is.finite(m$y), ]; if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), unit = unit[is.finite(m$lon) & is.finite(m$y)], y_nat = m$y, y_mod = m$y, X = Xr, w = m$w, Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}

out <- list()
for (tier in c("admin1", "admin2")) {
  # one LOCO fit per outcome; keep predictions per held-out country
  preds <- list()   # preds[[outcome]][[country]] = list(pred, y, w, unit)
  for (on in unique(TG$outcome)) {
    cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on, tier), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
    if (length(cl) < 3) next
    common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
    Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
    ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
    unit <- unlist(lapply(cl, function(z) z$unit))
    aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
                Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
    folds <- as.integer(factor(ctry)); pred <- rep(NA_real_, length(Y))
    for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
      Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
      p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, aux), error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) pred[te] <- p }
    for (cn in unique(ctry)) { k <- which(ctry == cn); preds[[on]][[cn]] <- list(pred = pred[k], y = ynat[k], w = wv[k], unit = unit[k]) }
    cat(tier, "fitted", on, "\n")
  }
  cells <- do.call(rbind, lapply(names(preds), function(on) data.frame(outcome = on, country = names(preds[[on]]))))
  sp <- function(y, p) { ok <- is.finite(y) & is.finite(p); if (sum(ok) < 3 || stats::sd(p[ok]) == 0) NA_real_ else suppressWarnings(stats::cor(y[ok], p[ok], method = "spearman")) }
  obs <- vapply(seq_len(nrow(cells)), function(i) { z <- preds[[cells$outcome[i]]][[cells$country[i]]]; sp(z$y, z$pred) }, numeric(1))
  obs_pos <- sum(obs > 0, na.rm = TRUE); obs_mean <- mean(obs, na.rm = TRUE)
  # null: permute each country's units once per replicate, apply to all its outcomes
  null_pos <- integer(NPERM); null_mean <- numeric(NPERM)
  for (b in seq_len(NPERM)) {
    perm <- list()
    for (cn in unique(cells$country)) { units_cn <- unique(unlist(lapply(preds, function(po) if (!is.null(po[[cn]])) po[[cn]]$unit else NULL)))
      perm[[cn]] <- stats::setNames(sample(units_cn), units_cn) }
    r <- vapply(seq_len(nrow(cells)), function(i) { z <- preds[[cells$outcome[i]]][[cells$country[i]]]
      idx <- match(perm[[cells$country[i]]][z$unit], z$unit); sp(z$y[idx], z$pred) }, numeric(1))
    null_pos[b] <- sum(r > 0, na.rm = TRUE); null_mean[b] <- mean(r, na.rm = TRUE)
  }
  out[[tier]] <- data.frame(tier = tier, cells = nrow(cells), obs_positive = obs_pos, obs_mean_rho = round(obs_mean, 3),
    p_positive = mean(null_pos >= obs_pos), p_mean_rho = mean(null_mean >= obs_mean),
    null_positive_q50 = median(null_pos), null_positive_q95 = unname(stats::quantile(null_pos, 0.95)),
    null_mean_q95 = round(unname(stats::quantile(null_mean, 0.95)), 3), stringsAsFactors = FALSE)
  cat(sprintf("\n===== NC-01 · %s · %d cells =====\nobserved: %d positive, mean rho %.3f\nnull (country-block permutation, %d reps): positive median %d, 95th pct %d; mean rho 95th pct %.3f\np(null >= observed): cells positive %.4f | mean rho %.4f\n",
              tier, nrow(cells), obs_pos, obs_mean, NPERM, median(null_pos), unname(stats::quantile(null_pos, 0.95)), unname(stats::quantile(null_mean, 0.95)), mean(null_pos >= obs_pos), mean(null_mean >= obs_mean)))
}
R <- bind_rows(out); write.csv(R, file.path(OUTDIR, "transport_null_calibration.csv"), row.names = FALSE); print(R, row.names = FALSE); cat("\nDONE\n")
