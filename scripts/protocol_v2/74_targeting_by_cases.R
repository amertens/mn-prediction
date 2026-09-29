# =============================================================================
# scripts/protocol_v2/74_targeting_by_cases.R   [TC-01]
#
# TARGETING BY PREDICTED CASES, NOT RATE, UNDER TWO BUDGET RULES
#
# Why. Script 12 ranks districts by predicted RATE and counts the deficient
# people in the worst-ranked fifth. With the country's own survey removed the
# flagged fifth is more deficient than the country (22.0 against 19.1 per cent,
# higher in 13 of 22) but holds only 16.8 per cent of its deficient people,
# less than a random fifth, because high-rate districts are small (28 Sep).
# A programme with a budget in DISTRICTS should rank by expected cases (rate x
# population); one with a budget in PEOPLE reached should rank by rate. This
# scores both, against population alone, which needs no model at all.
#
# PRE-REGISTERED DESIGN (written 28 Sep 2026 before any result was seen)
# Units, cells, folds and the rate ranking are script 12's exactly (prevalence
# target, headline tiers, Malawi at Traditional Authorities, in-fill 5-fold x 10
# draws, country transport pooled as 02b); capture is among surveyed districts.
#   Framing 1, a budget of districts (k = round(0.2 n)): share of the deficient
#     people inside the k selected districts. Arms:
#       pop_only        rank by target-group population (no model)
#       model_rate      rank by the index (script 12's arm; reproduction check)
#       model_cases     rank by predicted prevalence x population. In-fill: the
#                       calibrated index (domain_index_cal), expit of its logit.
#                       Transport: the held-out country anchored at its own
#                       national prevalence (n_eff-weighted mean of its surveyed
#                       districts) and spread by the training countries, as
#                       AR-01: expit(logit(p_nat) + rho_train * sd_train * z),
#                       z the index standardised within the held-out country,
#                       rho_train the mean nested leave-one-country-out Spearman
#                       among the training countries (floored at 0), sd_train
#                       the mean SD of logit prevalence in the training countries.
#       region_rate, region_cases   the survey's regional figure (region_mean_jk,
#                       fold-based), in-fill only
#       oracle_rate, oracle_cases   the true district prevalence / cases (ceilings)
#       random          k / n (expected)
#   Framing 2, a budget of people (20% of the target population; districts taken
#     in order until the budget is spent, the last one in part): share of the
#     deficient people reached. Arms: model_rate, region_rate, oracle_rate,
#     random = 0.20.
#   Primary reading, per estimand (in-fill, transport), over the measurable
#   combinations (results/figures/mnf15/cell_master.csv keep):
#     F1: model_cases beats pop_only on the mean AND in a majority of cells.
#     F2: model_rate beats 0.20 on the mean AND in a majority of cells.
#   Both framings and both estimands are reported whatever they show; all cells
#   (not only measurable) in the cell table.
#
#   Rscript -e "source('scripts/protocol_v2/74_targeting_by_cases.R')"
# -> results/tables/protocol_v2/tc01_targeting_by_cases.csv        (cell x estimand x arm x framing x rep)
# -> results/tables/protocol_v2/tc01_targeting_by_cases_cells.csv  (mean over draws)
# -> results/tables/protocol_v2/tc01_targeting_by_cases_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R")

OUTDIR <- "results/tables/protocol_v2"
TOPFRAC <- 0.20
REPS <- as.integer(Sys.getenv("NCE_REPS", "10"))
set.seed(20260991L)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
POP <- readRDS("dashboard/data/admin2_population.rds")
POP$country <- gsub(" ", "", POP$country)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
NCE <- read.csv(file.path(OUTDIR, "nce_targeting_metrics.csv"))
KEEP <- read.csv("results/figures/mnf15/cell_master.csv") |> filter(keep) |> transmute(country, outcome, measurable = TRUE)

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"
clamp <- function(p, eps = 0.005) pmin(pmax(p, eps), 1 - eps)

# as script 12
build <- function(cn, on) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_prev) & is.finite(t$n_eff), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2")) |>
    inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2"))
  m <- join_admin2_v2(m, admin2_population_v2(POP, cn, pop_for(on)), what = paste("pop", cn, on), quiet = TRUE)
  m <- m[is.finite(m$pop) & m$pop > 0, ]
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  list(country = cn, outcome = on, n = nrow(m), y = m$y_prev, pop = m$pop,
       X = Xr, D = domain_representation_v2(Xr, domain_of),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = m$y_prev),
       Admin1 = m$Admin1, w = m$n_eff)
}

# Framing 1: share of burden in the k highest-scored districts (script 12's rule)
cap_districts <- function(y, pop, score) {
  ok <- is.finite(y) & is.finite(pop) & is.finite(score)
  if (sum(ok) < 5) return(c(NA, NA, NA))
  y <- y[ok]; pop <- pop[ok]; score <- score[ok]; burden <- y * pop
  k <- max(1, round(TOPFRAC * length(y)))
  sel <- order(score, decreasing = TRUE)[seq_len(k)]
  c(sum(burden[sel]) / sum(burden), sum(burden[sel]) / sum(pop[sel]), sum(burden) / sum(pop))
}
# Framing 2: share of burden reached with TOPFRAC of the population, the last district in part
cap_people <- function(y, pop, score) {
  ok <- is.finite(y) & is.finite(pop) & is.finite(score)
  if (sum(ok) < 5) return(NA_real_)
  y <- y[ok]; pop <- pop[ok]; score <- score[ok]; burden <- y * pop
  o <- order(score, decreasing = TRUE); budget <- TOPFRAC * sum(pop)
  cp <- cumsum(pop[o]); full <- which(cp <= budget)
  got <- sum(burden[o][full]); used <- if (length(full)) cp[max(full)] else 0
  nxt <- if (length(full)) max(full) + 1 else 1
  if (nxt <= length(o) && used < budget) got <- got + burden[o][nxt] * (budget - used) / pop[o][nxt]
  got / sum(burden)
}

rows <- list()
add <- function(cl, estimand, arm, framing, rep, v, prev_sel = NA, prev_nat = NA)
  rows[[length(rows) + 1L]] <<- data.frame(country = cl$country, outcome = cl$outcome, estimand = estimand, arm = arm,
                                           framing = framing, rep = rep, n_areas = cl$n, capture = v,
                                           prev_sel = prev_sel, prev_nat = prev_nat)
both <- function(cl, estimand, arm, rep, score, people = TRUE) {
  f1 <- cap_districts(cl$y, cl$pop, score); add(cl, estimand, arm, "districts", rep, f1[1], f1[2], f1[3])
  if (people) add(cl, estimand, arm, "people", rep, cap_people(cl$y, cl$pop, score))
}

cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
built <- list()
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL)
  if (is.null(cl)) next
  built[[paste(cn, on)]] <- cl
  ymod <- .v2_logit(cl$y)
  auxc <- c(cl$aux, list(target = "prev"))
  # ---- in-fill, replicated ----
  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
    P <- list(domain_index = rep(NA_real_, cl$n), domain_index_cal = rep(NA_real_, cl$n), region_mean_jk = rep(NA_real_, cl$n))
    for (a in names(P)) for (f in unique(folds)) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 12) next
      p <- tryCatch(ARMS_V2[[a]](tr, te, ymod, cl$X, cl$D, if (a == "domain_index_cal") auxc else cl$aux),
                    error = function(e) rep(NA_real_, length(te)))
      if (length(p) == length(te)) P[[a]][te] <- p
    }
    if (all(is.na(P$domain_index))) next
    both(cl, "infill", "model_rate", r, P$domain_index)
    both(cl, "infill", "model_cases", r, .v2_expit(P$domain_index_cal) * cl$pop, people = FALSE)
    both(cl, "infill", "region_rate", r, P$region_mean_jk)
    both(cl, "infill", "region_cases", r, .v2_expit(P$region_mean_jk) * cl$pop, people = FALSE)
  }
  if (any(vapply(rows, function(z) z$country == cn && z$outcome == on, NA))) {
    both(cl, "infill", "pop_only", 1L, cl$pop, people = FALSE)
    both(cl, "infill", "oracle_rate", 1L, cl$y)
    both(cl, "infill", "oracle_cases", 1L, cl$y * cl$pop, people = FALSE)
  }
  cat("in-fill done", cn, on, "\n")
}

# ---- country transport (script 12, section B), plus the anchored rate ----
loco_pred <- function(cl, common, tr_names, te_name) {
  use <- c(tr_names, te_name)
  Y <- unlist(lapply(cl[use], function(z) as.numeric(scale(.v2_logit(z$y)))))
  Dm <- do.call(rbind, lapply(cl[use], function(z) z$D[, common, drop = FALSE]))
  ctry <- rep(use, vapply(cl[use], function(z) z$n, 0L))
  aux <- list(Admin1 = paste(ctry, unlist(lapply(cl[use], function(z) z$Admin1))), y_nat = Y)
  tr <- which(ctry %in% tr_names); te <- which(ctry == te_name)
  if (length(tr) < 20) return(NULL)
  p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, Dm, Dm, aux), error = function(e) NULL)
  if (is.null(p) || length(p) != length(te)) NULL else p
}
RHO <- list()
for (on in unique(cells$outcome)) {
  cl <- built[paste(COUNTRIES, on)]; cl <- cl[!vapply(cl, is.null, NA)]
  if (length(cl) < 3) next
  names(cl) <- vapply(cl, function(z) z$country, "")
  common <- Reduce(intersect, lapply(cl, function(z) colnames(z$D)))
  if (length(common) < 5) next
  for (h in names(cl)) {
    pool <- setdiff(names(cl), h); Z <- cl[[h]]
    pred <- loco_pred(cl, common, pool, h); if (is.null(pred)) next
    rt <- vapply(pool, function(t) { p2 <- loco_pred(cl, common, setdiff(pool, t), t)
      if (is.null(p2)) NA_real_ else suppressWarnings(stats::cor(cl[[t]]$y, p2, method = "spearman")) }, numeric(1))
    rho_tr <- mean(rt, na.rm = TRUE); if (!is.finite(rho_tr) || rho_tr < 0) rho_tr <- 0
    sd_tr <- mean(vapply(pool, function(t) stats::sd(.v2_logit(clamp(cl[[t]]$y))), numeric(1)), na.rm = TRUE)
    p_nat <- sum(Z$y * Z$w) / sum(Z$w)
    z <- as.numeric(scale(pred))
    rate_hat <- .v2_expit(.v2_logit(clamp(p_nat)) + rho_tr * sd_tr * z)
    both(Z, "country", "model_rate", 1L, pred)
    both(Z, "country", "model_cases", 1L, rate_hat * Z$pop, people = FALSE)
    both(Z, "country", "pop_only", 1L, Z$pop, people = FALSE)
    both(Z, "country", "oracle_rate", 1L, Z$y)
    both(Z, "country", "oracle_cases", 1L, Z$y * Z$pop, people = FALSE)
    RHO[[length(RHO) + 1L]] <- data.frame(country = h, outcome = on, rho_train = rho_tr, sd_train = sd_tr, p_nat = p_nat)
    cat(sprintf("loco done %s %s  rho_train %.2f sd_train %.2f\n", on, h, rho_tr, sd_tr))
  }
}

R <- bind_rows(rows)
# random: k / n districts (framing 1), 0.20 of people (framing 2)
RND <- R |> distinct(country, outcome, estimand, n_areas) |>
  mutate(k = pmax(1, round(TOPFRAC * n_areas)))
R <- bind_rows(R,
  RND |> transmute(country, outcome, estimand, arm = "random", framing = "districts", rep = 1L, n_areas, capture = k / n_areas),
  RND |> transmute(country, outcome, estimand, arm = "random", framing = "people", rep = 1L, n_areas, capture = TOPFRAC))

# ---- reproduction: the rate ranking must equal script 12's published capture ----
rep_in <- R |> filter(estimand == "infill", arm == "model_rate", framing == "districts") |>
  inner_join(NCE |> filter(estimand == "infill", arm == "domain_index") |> select(country, outcome, rep, cap12 = capture_top20),
             by = c("country", "outcome", "rep"))
rep_lc <- R |> filter(estimand == "country", arm == "model_rate", framing == "districts") |>
  inner_join(NCE |> filter(estimand == "country", arm == "domain_index") |> select(country, outcome, cap12 = capture_top20),
             by = c("country", "outcome"))
chk(identical(is.na(rep_in$capture), is.na(rep_in$cap12)), "the same in-fill cells are missing in both (Sierra Leone: too few districts)")
rep_in <- rep_in[is.finite(rep_in$capture), ]
cat(sprintf("reproduction: in-fill %d rows, max |diff| %.2e; transport %d rows, max |diff| %.2e\n",
            nrow(rep_in), max(abs(rep_in$capture - rep_in$cap12)), nrow(rep_lc), max(abs(rep_lc$capture - rep_lc$cap12))))
chk(nrow(rep_in) > 100 && max(abs(rep_in$capture - rep_in$cap12)) < 1e-6, "in-fill rate capture reproduces script 12")
chk(nrow(rep_lc) >= 20 && max(abs(rep_lc$capture - rep_lc$cap12)) < 1e-6, "transport rate capture reproduces script 12")

write.csv(R, file.path(OUTDIR, "tc01_targeting_by_cases.csv"), row.names = FALSE)
CELLS <- R |> group_by(country, outcome, estimand, arm, framing, n_areas) |>
  summarise(capture = mean(capture), prev_sel = mean(prev_sel), prev_nat = mean(prev_nat), .groups = "drop") |>
  left_join(KEEP, by = c("country", "outcome")) |> mutate(measurable = !is.na(measurable))
CELLS <- CELLS |> left_join(bind_rows(RHO) |> mutate(estimand = "country"), by = c("country", "outcome", "estimand"))
write.csv(CELLS, file.path(OUTDIR, "tc01_targeting_by_cases_cells.csv"), row.names = FALSE)

# Sierra Leone has no in-country test (14 districts leave fewer than 12 to train on), so its two
# measurable combinations drop out of the in-fill rows, as in script 34 (summary only, design unchanged)
ok_in <- CELLS |> filter(estimand == "infill", arm == "model_rate", is.finite(capture)) |> transmute(key = paste(country, outcome))
CELLS <- CELLS |> filter(is.finite(capture), estimand == "country" | paste(country, outcome) %in% ok_in$key)
W <- CELLS |> filter(measurable) |> select(country, outcome, estimand, framing, arm, capture) |>
  tidyr::pivot_wider(names_from = arm, values_from = capture)
SUMM <- CELLS |> filter(measurable) |> group_by(estimand, framing, arm) |>
  summarise(cells = dplyr::n(), mean_capture = round(mean(capture), 3), median_capture = round(median(capture), 3), .groups = "drop")
write.csv(SUMM, file.path(OUTDIR, "tc01_targeting_by_cases_summary.csv"), row.names = FALSE)
cat("\n== summary over measurable combinations ==\n"); print(as.data.frame(SUMM), row.names = FALSE)
for (e in c("infill", "country")) {
  a <- W |> filter(estimand == e, framing == "districts"); b <- W |> filter(estimand == e, framing == "people")
  cat(sprintf("\n[%s] F1 districts: model_cases %.3f vs pop_only %.3f; better in %d of %d  -> %s\n", e,
              mean(a$model_cases), mean(a$pop_only), sum(a$model_cases > a$pop_only), nrow(a),
              ifelse(mean(a$model_cases) > mean(a$pop_only) && sum(a$model_cases > a$pop_only) > nrow(a) / 2, "PASS", "FAIL")))
  cat(sprintf("[%s] F2 people: model_rate %.3f vs 0.20; above in %d of %d  -> %s\n", e,
              mean(b$model_rate), sum(b$model_rate > TOPFRAC), nrow(b),
              ifelse(mean(b$model_rate) > TOPFRAC && sum(b$model_rate > TOPFRAC) > nrow(b) / 2, "PASS", "FAIL")))
}
