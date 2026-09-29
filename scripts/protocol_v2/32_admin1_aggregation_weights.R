# =============================================================================
# scripts/protocol_v2/32_admin1_aggregation_weights.R   [AG-01]
#
# IS THE 12/12 REGIONAL TRANSPORT RESULT SENSITIVE TO HOW REGIONS ARE BUILT?
#
# The headline regional result (script 16) aggregates district outcomes to
# regions weighting by effective survey n, and aggregates predictors by a
# simple mean of districts. Both are choices. A programme would think of a
# region's prevalence as population-weighted; a remote-sensing analyst might
# average predictors over area or population. If the 12/12 depends on the
# choice it is fragile; if it does not, it can be quoted without a footnote.
#
# Three aggregation schemes on identical LOCO folds, domain_index arm:
#   neff      outcome weighted by effective n, predictors simple mean  (record)
#   pop       outcome weighted by target-group population, predictors simple mean
#   pop_both  outcome AND predictors weighted by population
#
#   Rscript scripts/protocol_v2/32_admin1_aggregation_weights.R
# -> results/tables/protocol_v2/admin1_aggregation_weights.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; set.seed(20260903L)
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
POP <- readRDS("dashboard/data/admin2_population.rds"); BND <- readRDS("dashboard/data/admin2_boundaries.rds")
POP$country <- gsub(" ", "", POP$country)  # (label also normalised inside admin2_population_v2(), JK-01)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w); if (!any(ok)) NA_real_ else stats::weighted.mean(x[ok], w[ok]) }

build_a1 <- function(cn, on, target, scheme) {
  ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]; if (!nrow(t)) return(NULL)
  t <- join_admin2_v2(t, admin2_population_v2(POP, cn, pop_for(on)), what = paste("pop", cn, on), quiet = TRUE)   # JK-01 pair key, no fan
  t <- t[is.finite(t$pop) & t$pop > 0, ]; if (!nrow(t)) return(NULL)
  wy <- if (scheme == "neff") t[[wcol]] else t$pop
  a1 <- t |> mutate(wy = wy) |> group_by(Admin1) |> summarise(y = wmean(.data[[ycol]], wy), w = sum(.data[[wcol]]), .groups = "drop") |> filter(is.finite(y))
  # Predictor side. `pp` was undefined here until 2026-09-29, so build_a1() threw
  # "object 'pp' not found" for every country, the tryCatch below swallowed it and
  # the script wrote an EMPTY table over its own result. It is the same
  # outcome-specific population the outcome side uses above, and it joins on the
  # pair key: S is the 243-row Malawi predictor set, where TA Lundu, TA Ngabu,
  # TA Pemba and TA Malemia each occur in two Admin1 regions, so a name-only join
  # fans 243 rows to 251.
  pp <- admin2_population_v2(POP, cn, pop_for(on))
  sc <- join_admin2_v2(S[S$country == cn, c("Admin1", "Admin2", PREDS)], pp,
                       what = paste("preds", cn, on), quiet = TRUE)
  x1 <- if (scheme == "pop_both") sc |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ wmean(.x, pop)), .groups = "drop")
        else sc |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
  m <- a1 |> inner_join(x1, by = "Admin1") |> inner_join(c1, by = "Admin1"); m <- m[is.finite(m$lon), ]; if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), y_nat = m$y, y_mod = if (target == "prev") .v2_logit(m$y) else m$y, X = Xr, w = m$w, Admin1 = m$Admin1, lon = m$lon, lat = m$lat)
}

rows <- list()
for (scheme in c("neff", "pop", "pop_both")) for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  # A bare error -> NULL hid the undefined `pp` for three weeks: every country
  # failed, every combination was skipped, and the run still "succeeded". Report.
  cl <- list(); for (cn in COUNTRIES) {
    z <- tryCatch(build_a1(cn, on, target, scheme),
                  error = function(e) { message(sprintf("  [skip] %s / %s / %s / %s: %s",
                                                        scheme, target, on, cn, conditionMessage(e))); NULL })
    if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X))); if (length(common) < 20) next
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod)))); Xm <- do.call(rbind, lapply(cl, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); ynat <- unlist(lapply(cl, function(z) z$y_nat)); wv <- unlist(lapply(cl, function(z) z$w))
  aux <- list(lon = unlist(lapply(cl, function(z) z$lon)), lat = unlist(lapply(cl, function(z) z$lat)),
              Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$Admin1))), y_nat = Y)
  folds <- as.integer(factor(ctry)); pred <- rep(NA_real_, length(Y))
  for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, aux), error = function(e) rep(NA_real_, length(te)))
    if (length(p) == length(te)) pred[te] <- p }
  for (cn in unique(ctry)) { k <- which(ctry == cn); s <- score_v2(ynat[k], pred[k], wv[k], scale = if (target == "prev") "prev" else "level")
    rows[[length(rows) + 1L]] <- data.frame(scheme = scheme, target = target, outcome = on, country = cn, n_units = length(k), spearman = s$spearman, thin = length(k) < 8, stringsAsFactors = FALSE) }
  cat("done", scheme, target, on, "\n")
}
R <- bind_rows(rows)
# Never overwrite a real result with an empty frame. bind_rows(list()) is 0 rows,
# so with the old blanket tryCatch a total failure wrote an empty CSV over the
# previous run and reported nothing.
if (!nrow(R)) stop("[AG-01] no cells were built - refusing to overwrite ",
                   file.path(OUTDIR, "admin1_aggregation_weights.csv"),
                   ". See the [skip] lines above for why each country failed.")
write.csv(R, file.path(OUTDIR, "admin1_aggregation_weights.csv"), row.names = FALSE)
cat("\n===== AG-01: regional transport (LOCO Spearman, domain_index) by aggregation scheme =====\n")
print(as.data.frame(R |> group_by(target, scheme) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
  median_rho = round(median(spearman, na.rm = TRUE), 3), cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |> arrange(target, scheme)), row.names = FALSE)
cat("\n(excluding thin countries, n_units < 8)\n")
print(as.data.frame(R[!R$thin, ] |> group_by(target, scheme) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3),
  cells_positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop") |> arrange(target, scheme)), row.names = FALSE)
for (tg in unique(R$target)) { W <- pivot_wider(R[R$target == tg, c("outcome", "country", "scheme", "spearman")], names_from = scheme, values_from = spearman)
  for (a in c("pop", "pop_both")) { d <- W[[a]] - W$neff
    cat(sprintf("  %-5s %-8s vs neff: better in %2d of %2d | median %+.3f | max |diff| %.3f\n", tg, a, sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE), max(abs(d), na.rm = TRUE))) } }
cat("\nDONE\n")
