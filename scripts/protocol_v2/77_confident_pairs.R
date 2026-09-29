# =============================================================================
# scripts/protocol_v2/77_confident_pairs.R   [SEP-01]
#
# ACCURACY ON THE PAIRS THE SURVEY CAN ACTUALLY TELL APART
#
# The talk's "pairs" metric (share of district pairs the prediction puts in the
# survey's order; 50% = coin toss) counts every pair, including pairs whose
# survey difference is well inside the survey's own sampling noise. Those pairs
# pull every method toward 50%. SEP-01 scores the same held-out predictions on
# the pairs the survey itself separates with confidence.
#
# DESIGN (fixed 2026-09-28, before any SEP-01 result was computed)
# ----------------------------------------------------------------
# Per-district standard error
#   level target      SE_d = sd_level_d / sqrt(n_eff_cont_d). sd_level is the
#                     respondent-level survey-weighted SD of the district's
#                     negated log biomarker (01_build_targets_v2.R:
#                     sqrt(.v2_wvar(ycont, w))), NOT a standard error, so it is
#                     divided by the root of the design-effect-adjusted n.
#                     Districts with one respondent have sd_level = NA (8 of
#                     1,342 rows); their SE is set to Inf, so none of their
#                     pairs is confident and each gets probability weight 0.
#   prevalence target SE_d = sqrt(p_d (1 - p_d) / n_eff_d).
# Confident pair      |y_i - y_j| > 1.96 * sqrt(SE_i^2 + SE_j^2). Single,
#                     pre-specified threshold; not tuned.
# Arms (each district's held-out prediction averaged over the ten in-fill draws,
#   as scripts/policy_deck/28_mnf15_v3_percell_pairs.R does; 5-fold district
#   folds from make_folds_v2("kfold_district", n, k = 5, rep_id = 1..10); on the
#   prevalence target each draw is back-transformed with expit before
#   averaging, as 73_predicted_vs_observed.R does)
#   domain_index   the model
#   spatial        the neighbour-map GAM on centroids (ARMS_V2$spatial)
#   regional_fair  the survey's regional figure scored fairly, exactly as
#                  scripts/policy_deck/39_fair_regional_pairs.py: each
#                  district's figure is the mean of the OTHER surveyed
#                  districts in its region (all other districts if it is alone
#                  in its region); two districts in the same region, or with
#                  equal figures, count as a tie = half a pair
#   coin           50%
# Pair scoring        survey-tied pairs (y_i == y_j) are skipped everywhere.
#                     domain_index and spatial: ties in the prediction are
#                     skipped (as the deck's concord()). regional_fair: ties
#                     count as half (as script 39).
# Metrics per cell
#   1 share_confident  confident pairs / all n(n-1)/2 pairs
#   2 pairs_confident  pairs in the survey's order among confident pairs
#   3 pairs_all        the same over all pairs (reproduces the deck)
#   4 pairs_probw      every pair weighted by 2*Phi(|d| / SE_diff) - 1, the
#                      probability the survey's order is right minus the
#                      probability it is wrong (SE_diff = 0 -> weight 1,
#                      SE_diff = Inf -> weight 0)
# Targets             level primary, prevalence secondary.
# Cells               the 14 measurable in-country combinations of the deck's
#                     per-nutrient slide (cell_master.csv keep == TRUE, Sierra
#                     Leone excluded: no in-country test) primary; all 18
#                     in-country combinations of the benchmark also reported.
# Summary             unweighted mean over cells (the deck's convention),
#                     plus pairs pooled over cells as a secondary view, and
#                     the number of cells in which the model beats the fair
#                     regional figure on confident pairs. A cell with no
#                     confident pair (prevalence, Ghana and Malawi women's
#                     vitamin A; neither is in the 14) has no confident-pairs
#                     value and drops out of that mean and count
#                     (cells_with_confident says how many remain). This NA
#                     handling was added after the first run, whose summary
#                     showed NA for those counts; nothing else changed.
# Reproduction checks (the script stops if any fails)
#   - domain_index and spatial: mean per-draw Spearman within 0.005 of
#     benchmarks_v2_cells.csv (infill), both targets;
#   - level, domain_index: pairs_all within 0.005 of
#     results/tables/policy_deck/v3_percell_pairs.csv;
#   - level, regional_fair: pairs_all within 0.005 of
#     results/tables/policy_deck/v6_fair_regional_pairs.csv (14 cells).
#
#   Rscript -e "source('scripts/protocol_v2/77_confident_pairs.R')"
# -> results/tables/protocol_v2/sep01_cells.csv          per cell x target x arm
# -> results/tables/protocol_v2/sep01_summary.csv        mean over cells / pooled
# -> results/tables/protocol_v2/sep01_districts.csv      per-district SE and averaged predictions
# -> results/tables/protocol_v2/sep01_reproduction.csv   the checks above
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
stopifnot(Sys.getenv("V2_PREDICTOR_TIERS") == "open,survey_public")
source("R/protocol_v2.R")

P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"; REPS <- 10L; TOL <- 0.005; Z <- 1.96
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND  <- readRDS("dashboard/data/admin2_boundaries.rds")
CELL <- read.csv(file.path(P2, "benchmarks_v2_cells.csv")) |> filter(estimand == "infill")
V3   <- read.csv(file.path(PD, "v3_percell_pairs.csv")) |> filter(estimand == "infill", arm == "domain_index")
V6   <- read.csv(file.path(PD, "v6_fair_regional_pairs.csv"))
CM   <- read.csv("results/figures/mnf15/cell_master.csv")
KEEP14 <- paste(CM$country, CM$outcome)[CM$keep %in% TRUE & CM$country != "SierraLeone"]
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

# as 73_predicted_vs_observed.R build_cell (itself verbatim from 02_run_benchmarks_v2.R), plus the SE
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
  se <- if (target == "prev") sqrt(y_nat * (1 - y_nat) / m$n_eff) else m$sd_level / sqrt(m$n_eff_cont)
  se[!is.finite(se)] <- Inf                                   # one-respondent districts: never confident
  list(country = cn, outcome = on, target = target, y_nat = y_nat, y_mod = y_mod, X = Xr, D = D,
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = y_nat, target = target),
       se = se, Admin1 = m$Admin1, Admin2 = m$Admin2, n = nrow(m))
}

# all pairs i < j with the survey difference, its SE, the confident flag and the probability weight
pair_frame <- function(y, se, a1) {
  n <- length(y); ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]; d <- y[i] - y[j]; sed <- sqrt(se[i]^2 + se[j]^2)
  conf <- is.finite(sed) & abs(d) > Z * sed
  w <- ifelse(!is.finite(sed), 0, ifelse(sed > 0, 2 * pnorm(abs(d) / sed) - 1, as.numeric(d != 0)))
  data.frame(i = i, j = j, d = d, sed = sed, conf = conf, w = w, same = a1[i] == a1[j])
}
# per-pair score in [0, 1] and whether the pair counts; survey ties never count
score_model <- function(pf, p) {           # prediction ties skipped (the deck's concord())
  s <- sign(pf$d) * sign(p[pf$i] - p[pf$j])
  list(score = as.numeric(s > 0), use = s != 0)
}
score_regional_fair <- function(pf, fig) { # same region or equal figures = half a pair (script 39)
  dp <- ifelse(pf$same, 0, fig[pf$i] - fig[pf$j])
  sc <- ifelse(dp == 0, 0.5, as.numeric(sign(pf$d) == sign(dp)))
  list(score = sc, use = pf$d != 0)
}
fair_region_figure <- function(y, a1) {    # verbatim logic of script 39's loo
  n <- length(y)
  vapply(seq_len(n), function(i) {
    same <- which(a1 == a1[i] & seq_len(n) != i)
    if (length(same)) mean(y[same]) else mean(y[-i])
  }, 0)
}
metrics <- function(pf, sc) {
  u <- sc$use; cu <- u & pf$conf
  data.frame(n_used_all = sum(u), n_used_conf = sum(cu),
             agree_all = sum(sc$score[u]), agree_conf = sum(sc$score[cu]),
             pairs_all = if (any(u)) mean(sc$score[u]) else NA_real_,
             pairs_confident = if (any(cu)) mean(sc$score[cu]) else NA_real_,
             w_sum = sum(pf$w[u]), w_agree = sum(pf$w[u] * sc$score[u]),
             pairs_probw = if (sum(pf$w[u]) > 0) sum(pf$w[u] * sc$score[u]) / sum(pf$w[u]) else NA_real_)
}

ARMS <- c("domain_index", "spatial", "region_mean_jk")   # region_mean_jk only for the reproduction check
cell_rows <- list(); dist_rows <- list(); chk_rows <- list()
chk <- function(what, cn, on, target, got, ref) {
  ok <- length(ref) == 1 && is.finite(ref) && is.finite(got) && abs(got - ref) <= TOL
  chk_rows[[length(chk_rows) + 1]] <<- data.frame(check = what, country = cn, outcome = on, target = target,
                                                  value = got, reference = if (length(ref) == 1) ref else NA_real_, pass = ok)
  if (!ok) stop(sprintf("REPRODUCTION FAILED: %s %s %s %s: %.4f vs %s", what, cn, on, target, got, paste(ref, collapse = ",")))
}

for (target in c("level", "prev")) {
  for (cn in COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) {
    cc <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (is.null(cc)) next
    bref <- CELL[CELL$country == cn & CELL$outcome == on & CELL$target == target, ]
    if (!any(bref$arm == "domain_index" & is.finite(bref$spearman))) { cat("skipped (not scored in the benchmark):", cn, on, target, "\n"); next }
    n <- cc$n; scale <- if (target == "prev") "prev" else "level"
    P <- setNames(lapply(ARMS, function(a) matrix(NA_real_, n, REPS)), ARMS)
    sp <- setNames(lapply(ARMS, function(a) numeric(REPS)), ARMS)
    for (r in seq_len(REPS)) {
      folds <- make_folds_v2("kfold_district", n, k = 5, rep_id = r)
      for (a in ARMS) {
        pm <- rep(NA_real_, n)
        for (f in unique(folds)) {
          te <- which(folds == f); tr <- which(folds != f)
          if (length(tr) < 12) next
          pm[te] <- ARMS_V2[[a]](tr, te, cc$y_mod, cc$X, cc$D, cc$aux)
        }
        P[[a]][, r] <- if (target == "prev") .v2_expit(pm) else pm
        sp[[a]][r] <- score_v2(cc$y_nat, P[[a]][, r], NULL, scale = scale)$spearman
      }
    }
    for (a in c("domain_index", "spatial"))
      chk(paste0("spearman_", a), cn, on, target, mean(sp[[a]]), bref$spearman[bref$arm == a])

    pbar <- lapply(P, rowMeans)
    fig  <- fair_region_figure(cc$y_nat, cc$Admin1)
    pf   <- pair_frame(cc$y_nat, cc$se, cc$Admin1)
    M <- list(domain_index  = metrics(pf, score_model(pf, pbar$domain_index)),
              spatial       = metrics(pf, score_model(pf, pbar$spatial)),
              regional_fair = metrics(pf, score_regional_fair(pf, fig)),
              regional_jk_as_plotted = metrics(pf, score_model(pf, pbar$region_mean_jk)))
    if (target == "level") {
      chk("pairs_all_domain_index_vs_v3", cn, on, target, M$domain_index$pairs_all, V3$pairs[V3$country == cn & V3$outcome == on])
      if (paste(cn, on) %in% KEEP14)
        chk("pairs_all_regional_fair_vs_v6", cn, on, target, M$regional_fair$pairs_all, V6$regional_fair[V6$country == cn & V6$outcome == on])
    }
    npairs <- nrow(pf)
    for (a in names(M)) cell_rows[[length(cell_rows) + 1]] <- data.frame(
      country = cn, outcome = on, target = target, keep14 = paste(cn, on) %in% KEEP14, n = n, n_pairs = npairs,
      n_confident = sum(pf$conf), share_confident = sum(pf$conf) / npairs, mean_weight = mean(pf$w),
      arm = a, M[[a]])
    cell_rows[[length(cell_rows) + 1]] <- data.frame(
      country = cn, outcome = on, target = target, keep14 = paste(cn, on) %in% KEEP14, n = n, n_pairs = npairs,
      n_confident = sum(pf$conf), share_confident = sum(pf$conf) / npairs, mean_weight = mean(pf$w), arm = "coin",
      n_used_all = sum(pf$d != 0), n_used_conf = sum(pf$conf), agree_all = 0.5 * sum(pf$d != 0), agree_conf = 0.5 * sum(pf$conf),
      pairs_all = 0.5, pairs_confident = 0.5, w_sum = sum(pf$w), w_agree = 0.5 * sum(pf$w), pairs_probw = 0.5)
    dist_rows[[length(dist_rows) + 1]] <- data.frame(country = cn, outcome = on, target = target, Admin1 = cc$Admin1,
      Admin2 = cc$Admin2, observed = cc$y_nat, se = cc$se, pred_domain_index = pbar$domain_index,
      pred_spatial = pbar$spatial, regional_fair_figure = fig)
    cat(sprintf("%-5s %-11s %-12s n=%2d  confident %4.0f%% of %4d pairs | confident: index %.3f spatial %.3f regional %.3f | all: index %.3f\n",
                target, cn, on, n, 100 * sum(pf$conf) / npairs, npairs, M$domain_index$pairs_confident,
                M$spatial$pairs_confident, M$regional_fair$pairs_confident, M$domain_index$pairs_all))
  }
}

CE  <- bind_rows(cell_rows)
CHK <- bind_rows(chk_rows)
stopifnot(all(CHK$pass))
stopifnot(sum(CHK$check == "pairs_all_domain_index_vs_v3") == 18, sum(CHK$check == "pairs_all_regional_fair_vs_v6") == 14)
stopifnot(n_distinct(paste(CE$country, CE$outcome)[CE$target == "level" & CE$keep14]) == 14)

summ <- function(d, set) d |> group_by(target, arm) |> summarise(
  cells = n(), cells_with_confident = sum(is.finite(pairs_confident)),
  share_confident = mean(share_confident), n_confident_median = median(n_confident),
  pairs_confident = mean(pairs_confident, na.rm = TRUE), pairs_all = mean(pairs_all), pairs_probw = mean(pairs_probw),
  pooled_pairs_confident = sum(agree_conf) / sum(n_used_conf), pooled_pairs_all = sum(agree_all) / sum(n_used_all),
  pooled_pairs_probw = sum(w_agree) / sum(w_sum), .groups = "drop") |> mutate(cellset = set, .before = 1)
SU <- bind_rows(summ(CE[CE$keep14, ], "14 measurable"), summ(CE, "all 18"))
wins <- CE |> filter(arm %in% c("domain_index", "regional_fair", "spatial")) |>
  select(target, keep14, country, outcome, arm, pairs_confident) |>
  tidyr::pivot_wider(names_from = arm, values_from = pairs_confident)
WN <- bind_rows(
  wins |> filter(keep14) |> group_by(target) |> summarise(cellset = "14 measurable", index_beats_regional_conf = sum(domain_index > regional_fair, na.rm = TRUE),
    index_beats_spatial_conf = sum(domain_index > spatial, na.rm = TRUE), cells = n(), .groups = "drop"),
  wins |> group_by(target) |> summarise(cellset = "all 18", index_beats_regional_conf = sum(domain_index > regional_fair, na.rm = TRUE),
    index_beats_spatial_conf = sum(domain_index > spatial, na.rm = TRUE), cells = n(), .groups = "drop"))
SU <- left_join(SU, WN |> select(-cells), by = c("cellset", "target"))
SU$index_beats_regional_conf[SU$arm != "domain_index"] <- NA; SU$index_beats_spatial_conf[SU$arm != "domain_index"] <- NA

write.csv(CE, file.path(P2, "sep01_cells.csv"), row.names = FALSE)
write.csv(SU, file.path(P2, "sep01_summary.csv"), row.names = FALSE)
write.csv(bind_rows(dist_rows), file.path(P2, "sep01_districts.csv"), row.names = FALSE)
write.csv(CHK, file.path(P2, "sep01_reproduction.csv"), row.names = FALSE)

cat("\nreproduction: all", nrow(CHK), "checks pass (max abs diff", signif(max(abs(CHK$value - CHK$reference)), 3), ")\n")
options(width = 220)
print(as.data.frame(SU |> filter(arm != "regional_jk_as_plotted") |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
