# =============================================================================
# scripts/protocol_v2/25_nested_domain_selection.R   [DA-02]
#
# DOES A FEW-DOMAIN INDEX TRANSPORT BETTER, CHOSEN HONESTLY?
#
# DA-01 found that climate alone (0.351) and soil alone (0.318) each transport
# better than the full 18-domain index (0.252), and that ten domains are net
# dead weight. But those "only" sets were identified on the same LOCO folds
# they were then scored on, so quoting them is selection on the test.
#
# This script selects domains NESTED: for each held-out country, an inner
# leave-one-country-out over the three TRAINING countries scores each domain
# and greedily adds domains while inner transport improves (cap 5). The chosen
# set is then fitted on all three training countries (PCs oriented from them)
# and scored once on the held-out country. The held-out country never
# influences which domains are chosen.
#
# Comparators on the same folds:
#   full        all domains (the production index)
#   nested      domains chosen by inner LOCO   <- the honest few-domain index
#   fixed_cs    climate + soil                  } pre-specified from DA-01, so
#   fixed_csa   climate + soil + agriculture    } optimistic; reported labelled
#
#   Rscript scripts/protocol_v2/25_nested_domain_selection.R
# -> results/tables/protocol_v2/nested_domain_selection.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; MAXDOM <- 5L; set.seed(20260903L)

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
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2],
             stringsAsFactors = FALSE) }))
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
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
CS  <- c("Climate and weather", "Soil characteristics")
CSA <- c(CS, "Agricultural production, land use")

# transport score of a domain set: fit index on rows `tr`, predict rows `te`,
# with D oriented from `tr`; returns Spearman on te (natural scale)
score_set <- function(doms, tr, te, Y, Xm, ynat, wv, aux, target) {
  D <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  keep <- which(col_domain(colnames(D)) %in% doms); if (!length(keep)) return(NA_real_)
  p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, D[, keep, drop = FALSE], aux), error = function(e) NULL)
  if (is.null(p) || length(p) != length(te)) return(NA_real_)
  score_v2(ynat[te], p, wv[te], scale = if (target == "prev") "prev" else "level")$spearman
}

rows <- list(); chosen_log <- list()
for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
  Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
              Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  avail <- intersect(domains, unique(stats::na.omit(col_domain(colnames(
    domain_representation_v2(Xm, domain_of, sign_rows = seq_along(Y)))))))

  for (h in names(cl)) {
    te <- which(ctry == h); pool <- setdiff(names(cl), h); trall <- which(ctry %in% pool)
    if (length(trall) < 20 || length(pool) < 2) next
    # inner LOCO over the training countries only
    inner <- function(doms) mean(vapply(pool, function(g) {
      itr <- which(ctry %in% setdiff(pool, g)); ite <- which(ctry == g)
      if (length(itr) < 20) return(NA_real_)
      score_set(doms, itr, ite, Y, Xm, ynat, wv, aux, target) }, numeric(1)), na.rm = TRUE)
    single <- vapply(avail, function(d) inner(d), numeric(1))
    sel <- names(which.max(single)); best <- max(single, na.rm = TRUE)
    repeat {
      if (length(sel) >= MAXDOM) break
      cand <- setdiff(avail, sel); if (!length(cand)) break
      gain <- vapply(cand, function(d) inner(c(sel, d)), numeric(1))
      if (all(!is.finite(gain)) || max(gain, na.rm = TRUE) <= best + 1e-4) break
      sel <- c(sel, names(which.max(gain))); best <- max(gain, na.rm = TRUE)
    }
    chosen_log[[length(chosen_log) + 1L]] <- data.frame(target = target, outcome = on, heldout = h,
      n_chosen = length(sel), chosen = paste(sel, collapse = " | "), inner_rho = round(best, 3), stringsAsFactors = FALSE)
    for (a in c("full", "nested", "fixed_cs", "fixed_csa")) {
      doms <- switch(a, full = avail, nested = sel, fixed_cs = intersect(CS, avail), fixed_csa = intersect(CSA, avail))
      rows[[length(rows) + 1L]] <- data.frame(target = target, outcome = on, heldout = h, arm = a,
        n_domains = length(doms), spearman = score_set(doms, trall, te, Y, Xm, ynat, wv, aux, target), stringsAsFactors = FALSE)
    }
  }
  cat("nested done", target, on, "\n")
}
R <- bind_rows(rows); CH <- bind_rows(chosen_log)
write.csv(R, file.path(OUTDIR, "nested_domain_selection.csv"), row.names = FALSE)
write.csv(CH, file.path(OUTDIR, "nested_domain_selection_chosen.csv"), row.names = FALSE)

cat("\n===== DA-02: transport (LOCO Spearman) by domain set =====\n")
print(as.data.frame(R |> group_by(target, arm) |> summarise(cells = dplyr::n(),
  mean_rho = round(mean(spearman, na.rm = TRUE), 3), median_rho = round(median(spearman, na.rm = TRUE), 3),
  cells_positive = sum(spearman > 0, na.rm = TRUE), mean_n_domains = round(mean(n_domains), 1), .groups = "drop") |>
  arrange(target, desc(mean_rho))), row.names = FALSE)
for (tg in unique(R$target)) {
  W <- tidyr::pivot_wider(R[R$target == tg, c("outcome", "heldout", "arm", "spearman")], names_from = arm, values_from = spearman)
  hh <- function(a, b) { d <- W[[a]] - W[[b]]; sprintf("  %-5s %-10s vs %-10s better in %2d of %2d | median %+.3f", tg, a, b,
    sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE)) }
  cat(hh("nested", "full"), "\n", hh("fixed_cs", "full"), "\n", hh("fixed_csa", "full"), "\n", hh("nested", "fixed_cs"), "\n")
}
cat("\n===== domains chosen by inner LOCO (frequency across held-out x outcome x target) =====\n")
tab <- sort(table(unlist(strsplit(CH$chosen, " \\| "))), decreasing = TRUE)
print(round(100 * tab / nrow(CH), 1))
cat(sprintf("\nmean number of domains chosen: %.1f\n", mean(CH$n_chosen)))
cat("\nDONE\n")
