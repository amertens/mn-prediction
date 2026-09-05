# =============================================================================
# scripts/protocol_v2/35_anchor_and_rank.R   [AR-01]
#
# A NATIONAL ANCHOR PLUS A TRANSPORTED RANKING, COSTED AGAINST A DISTRICT SURVEY
#
# Two earlier results point at a survey design nobody has simulated. District
# RANKINGS transport across countries (LOCO Spearman 0.28-0.31 at the district
# rung, 0.50-0.56 regional) but LEVELS do not; and anchoring the transported
# surface to a country's own national estimate halves the level error while
# an external anchor triples it (LV-01). So the cheapest survey that makes a
# transported ranking usable is one that measures only the NATIONAL (or
# regional) prevalence and lets the covariates order the districts. This
# costs that design against the survey a country would otherwise field, at
# matched effective sample size.
#
# Designs, each given a fraction f of the country's full-survey effective
# sample (f x N_eff respondents, spread proportionally where relevant):
#   A1  national anchor + transported ranking:
#         district_d = expit(logit(anchor) + rho_train * sd_train * z_d)
#       z_d = leave-one-country-out domain_index prediction, standardised
#       within the held-out country; sd_train = between-district SD of logit
#       prevalence in the training countries; rho_train = nested LOCO Spearman
#       among the training countries (best-linear-predictor shrinkage). The
#       held-out country contributes nothing but the anchor.
#   A2  regional anchors + transported ranking WITHIN region.
#   B   district survey: the full survey's district estimate plus the
#       incremental noise of a fraction-f survey, var = (1/f - 1) p(1-p)/n_eff
#       (script 26's convention; zero at f = 1).
#   C   regional survey: regional estimates from the same sample, flat across
#       districts -- what the survey is designed to deliver.
# Two rankings (from the prevalence and from the biomarker-level target) and
# two domain sets (full; climate + soil) are run.
#
# Scored against the full survey's district prevalences on the consistent
# district rung (146 units; Malawi at its 27 districts in 3 real regions).
# B converges to the truth at f = 1 by construction and A never does, so the
# crossover -- the smallest f at which a district survey beats the anchored
# ranking -- is CONSERVATIVE for the anchored design.
#
#   Rscript scripts/protocol_v2/35_anchor_and_rank.R
# -> results/tables/protocol_v2/anchor_and_rank.csv, _summary.csv, _crossover.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"
FRACTIONS <- c(0.05, 0.10, 0.15, 0.25, 0.40, 0.60, 0.80, 1.00)
REPS <- as.integer(Sys.getenv("AR_REPS", "20")); TOPFRAC <- 0.20; MIN_TRAIN <- 20L; set.seed(20260904L)
MALAWI_REGION <- c(
  Chitipa = "Northern", Karonga = "Northern", Likoma = "Northern", Mzimba = "Northern",
  `Nkhata Bay` = "Northern", Rumphi = "Northern",
  Dedza = "Central", Dowa = "Central", Kasungu = "Central", Lilongwe = "Central", Mchinji = "Central",
  Nkhotakota = "Central", Ntcheu = "Central", Ntchisi = "Central", Salima = "Central",
  Balaka = "Southern", Blantyre = "Southern", Chikwawa = "Southern", Chiradzulu = "Southern",
  Machinga = "Southern", Mangochi = "Southern", Mulanje = "Southern", Mwanza = "Southern",
  Neno = "Southern", Nsanje = "Southern", Phalombe = "Southern", Thyolo = "Southern", Zomba = "Southern")
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds"); POP <- readRDS("dashboard/data/admin2_population.rds")
# the population file spells Sierra Leone with a space; targets_v2 does not.
# Without this line the join silently drops the country (first run: 12 cells, no Sierra Leone)
POP$country <- gsub(" ", "", POP$country)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
domains <- sort(unique(stats::na.omit(domain_of[PREDS])))
prefix_of <- stats::setNames(domains, make.names(substr(domains, 1, 12)))
col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
SETS <- list(full = domains, climate_soil = c("Climate and weather", "Soil characteristics"))
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"
wm <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
clamp <- function(p, eps = 0.005) pmin(pmax(p, eps), 1 - eps)

build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff) & t$n_eff > 0, ]; if (!nrow(t)) return(NULL)
  pc <- pop_for(on); pp <- POP[POP$country == cn, c("Admin2", pc)]; names(pp)[2] <- "pop"
  t <- left_join(t, pp, by = "Admin2"); t <- t[is.finite(t$pop) & t$pop > 0, ]; if (!nrow(t)) return(NULL)
  if (cn == "Malawi") {
    a <- t |> group_by(Admin1) |> summarise(y = wm(y_prev, n_eff), yl = wm(y_level, n_eff_cont), w = sum(n_eff), n_raw = sum(n_raw), pop = sum(pop), .groups = "drop")
    x <- S[S$country == cn, c("Admin1", PREDS)] |> group_by(Admin1) |> summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1") |> inner_join(c1, by = "Admin1"); unit <- m$Admin1; region <- unname(MALAWI_REGION[m$Admin1])
  } else {
    m <- t[, c("Admin1", "Admin2", "y_prev", "y_level", "n_eff", "n_raw", "pop")]; names(m)[3:5] <- c("y", "yl", "w")
    m <- m |> inner_join(S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2")) |>
      inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2")); unit <- m$Admin2; region <- m$Admin1
  }
  keep <- is.finite(m$lon) & is.finite(m$y) & !is.na(region); m <- m[keep, ]; unit <- unit[keep]; region <- region[keep]
  if (nrow(m) < 8 || dplyr::n_distinct(region) < 2) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), unit = unit, region = region, y = m$y, yl = m$yl, w = m$w, n_raw = m$n_raw, pop = m$pop, X = Xr)
}
yv <- function(z, target) if (target == "prev") .v2_logit(clamp(z$y)) else z$yl
# pooled fit on tr_names, domain_index predictions for every row of te_name; NULL if not possible
fit_pred <- function(cl, tr_names, te_name, target, set) {
  use <- c(tr_names, te_name)
  if (any(vapply(cl[use], function(z) sum(is.finite(yv(z, target))) < 5, TRUE))) return(NULL)
  common <- Reduce(intersect, lapply(cl[use], function(z) colnames(z$X))); if (length(common) < 20) return(NULL)
  Y <- unlist(lapply(cl[use], function(z) { v <- yv(z, target); v[!is.finite(v)] <- mean(v, na.rm = TRUE); as.numeric(scale(v)) }))
  Xm <- do.call(rbind, lapply(cl[use], function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(use, vapply(cl[use], function(z) z$n, 0L))
  tr <- which(ctry %in% tr_names); te <- which(ctry == te_name); if (length(tr) < MIN_TRAIN) return(NULL)
  Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  keep <- which(col_domain(colnames(Dm)) %in% SETS[[set]]); if (!length(keep)) return(NULL)
  aux <- list(Admin1 = paste(ctry, unlist(lapply(cl[use], function(z) z$region))), y_nat = Y)
  p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, Dm[, keep, drop = FALSE], aux), error = function(e) NULL)
  if (is.null(p) || length(p) != length(te) || !all(is.finite(p)) || stats::sd(p) == 0) return(NULL)
  p
}

rows <- list()
for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on), error = function(e) { cat("  build error", cn, on, conditionMessage(e), "\n"); NULL }); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  for (h in names(cl)) { pool <- setdiff(names(cl), h); Z <- cl[[h]]
    p <- Z$y; w <- Z$w; n <- Z$n; reg <- as.character(Z$region); pv <- clamp(p)
    N <- sum(w); p_nat <- sum(p * w) / N
    p_r <- tapply(p * w, reg, sum) / tapply(w, reg, sum); N_r <- tapply(w, reg, sum); ri <- match(reg, names(p_r))
    burden <- p * Z$pop; k_top <- max(1L, round(TOPFRAC * n))
    capture <- function(s) { sel <- order(s, decreasing = TRUE)[seq_len(k_top)]; sum(burden[sel]) / sum(burden) }
    for (target in c("prev", "level")) for (set in names(SETS)) {
      pred <- fit_pred(cl, pool, h, target, set); if (is.null(pred)) next
      z <- as.numeric(scale(pred)); zc <- z - stats::ave(z, reg)
      rt <- vapply(pool, function(t) { p2 <- fit_pred(cl, setdiff(pool, t), t, target, set)
        if (is.null(p2)) NA_real_ else suppressWarnings(stats::cor(cl[[t]]$y, p2, method = "spearman")) }, numeric(1))
      rho_tr <- mean(rt, na.rm = TRUE); if (!is.finite(rho_tr) || rho_tr < 0) rho_tr <- 0
      sd_tr <- mean(vapply(pool, function(t) stats::sd(.v2_logit(clamp(cl[[t]]$y))), numeric(1)), na.rm = TRUE)
      rho_obs <- suppressWarnings(stats::cor(p, z, method = "spearman"))
      for (f in FRACTIONS) for (r in seq_len(REPS)) {
        ex <- 1 / f - 1
        a_nat <- clamp(p_nat + stats::rnorm(1, 0, sqrt(ex * p_nat * (1 - p_nat) / N)))
        a_r <- clamp(p_r + stats::rnorm(length(p_r), 0, sqrt(ex * p_r * (1 - p_r) / N_r)))
        est <- list(
          A1_anchor_rank        = .v2_expit(.v2_logit(a_nat) + rho_tr * sd_tr * z),
          A2_region_anchor_rank = .v2_expit(.v2_logit(a_r[ri]) + rho_tr * sd_tr * zc),
          B_district_survey     = pmin(pmax(p + stats::rnorm(n, 0, sqrt(ex * pv * (1 - pv) / w)), 0), 1),
          C_regional_survey     = unname(a_r[ri]))
        for (nm in names(est)) { e <- est[[nm]]
          rows[[length(rows) + 1L]] <- data.frame(outcome = on, country = h, rank_from = target, set = set, fraction = f, rep = r, design = nm,
            n_units = n, n_raw_total = sum(Z$n_raw), n_eff_total = round(N, 1), rho_train = round(rho_tr, 3), rho_obs = round(rho_obs, 3),
            mae = 100 * mean(abs(p - e)), spearman = suppressWarnings(stats::cor(p, e, method = "spearman")), capture = capture(e), stringsAsFactors = FALSE) }
      }
    }
    cat("done", on, h, "\n")
  }
}
R <- bind_rows(rows); if (!nrow(R)) stop("no cells")
write.csv(R, file.path(OUTDIR, "anchor_and_rank.csv"), row.names = FALSE)
CELL <- R |> group_by(outcome, country, rank_from, set, fraction, design, n_units, n_raw_total, n_eff_total, rho_train, rho_obs) |>
  summarise(mae = mean(mae), spearman = mean(spearman, na.rm = TRUE), capture = mean(capture, na.rm = TRUE), .groups = "drop")
SUMM <- CELL |> group_by(rank_from, set, fraction, design) |> summarise(cells = dplyr::n(), mae = round(median(mae), 2), spearman = round(median(spearman, na.rm = TRUE), 3), capture = round(median(capture, na.rm = TRUE), 3), .groups = "drop")
write.csv(SUMM, file.path(OUTDIR, "anchor_and_rank_summary.csv"), row.names = FALSE)
# crossover: smallest f at which the district survey (B) beats each anchored design on MAE
W <- CELL |> select(outcome, country, rank_from, set, fraction, design, mae, n_raw_total) |> pivot_wider(names_from = design, values_from = mae)
XO <- W |> group_by(outcome, country, rank_from, set, n_raw_total) |> arrange(fraction, .by_group = TRUE) |>
  summarise(f_cross_A1 = { i <- which(B_district_survey < A1_anchor_rank); if (length(i)) fraction[min(i)] else NA_real_ },
            f_cross_A2 = { i <- which(B_district_survey < A2_region_anchor_rank); if (length(i)) fraction[min(i)] else NA_real_ },
            f_cross_C_A1 = { i <- which(C_regional_survey < A1_anchor_rank); if (length(i)) fraction[min(i)] else NA_real_ },
            mae_A1_f05 = A1_anchor_rank[fraction == 0.05], mae_B_f05 = B_district_survey[fraction == 0.05], mae_C_f05 = C_regional_survey[fraction == 0.05],
            mae_A1_f25 = A1_anchor_rank[fraction == 0.25], mae_B_f25 = B_district_survey[fraction == 0.25], .groups = "drop") |>
  mutate(resp_cross_A1 = round(f_cross_A1 * n_raw_total), resp_cross_A2 = round(f_cross_A2 * n_raw_total))
write.csv(XO, file.path(OUTDIR, "anchor_and_rank_crossover.csv"), row.names = FALSE)

cat("\n===== AR-01: national/regional anchor + transported ranking vs district survey (consistent rung) =====\n")
for (tg in unique(SUMM$rank_from)) for (st in names(SETS)) { s <- SUMM[SUMM$rank_from == tg & SUMM$set == st, ]; if (!nrow(s)) next
  cat(sprintf("\n-- ranking from %s target | domain set %s | %d cells --\n", tg, st, max(s$cells)))
  cat("median district MAE (pp):\n"); print(as.data.frame(pivot_wider(s[, c("fraction", "design", "mae")], names_from = design, values_from = mae)), row.names = FALSE)
  cat("median Spearman:\n"); print(as.data.frame(pivot_wider(s[, c("fraction", "design", "spearman")], names_from = design, values_from = spearman)), row.names = FALSE)
  cat("median burden captured in worst-ranked fifth:\n"); print(as.data.frame(pivot_wider(s[, c("fraction", "design", "capture")], names_from = design, values_from = capture)), row.names = FALSE)
  x <- XO[XO$rank_from == tg & XO$set == st, ]
  cat(sprintf("crossover (district survey beats A1 on MAE): median f %.2f | never within grid in %d of %d cells | median respondents at crossover %s\n",
              median(x$f_cross_A1, na.rm = TRUE), sum(is.na(x$f_cross_A1)), nrow(x), format(round(median(x$resp_cross_A1, na.rm = TRUE)))))
  cat(sprintf("crossover (district survey beats A2 on MAE): median f %.2f | never within grid in %d of %d cells\n", median(x$f_cross_A2, na.rm = TRUE), sum(is.na(x$f_cross_A2)), nrow(x)))
  cat(sprintf("at f = 0.05: A1 beats B in %d of %d cells; A1 beats C in %d of %d; at f = 0.25: A1 beats B in %d of %d\n",
              sum(x$mae_A1_f05 < x$mae_B_f05), nrow(x), sum(x$mae_A1_f05 < x$mae_C_f05), nrow(x), sum(x$mae_A1_f25 < x$mae_B_f25), nrow(x)))
}
cat("\n-- per cell (prev ranking, climate_soil): transported rho in the held-out country vs the shrinkage used --\n")
print(as.data.frame(CELL |> filter(rank_from == "prev", set == "climate_soil", fraction == 0.05, design == "A1_anchor_rank") |>
  select(outcome, country, n_units, n_raw_total, rho_train, rho_obs, mae_A1_f05 = mae) |> mutate(mae_A1_f05 = round(mae_A1_f05, 2))), row.names = FALSE)
cat("\nDONE\n")
