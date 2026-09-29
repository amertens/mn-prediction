# =============================================================================
# scripts/protocol_v2/76_nutrient_specificity_swap.R   [NX-01, 2026-09-28]
#
# IS EACH TRANSPORTED MAP SPECIFIC TO ITS NUTRIENT? A WEIGHT-SWAP NEGATIVE CONTROL
#
# The transported maps for different nutrients agree strongly with each other
# (Cote d'Ivoire: about 0.93 between outcomes). If so, the talk's "which
# deficiencies can it map" table may describe one general deprivation map
# rather than nutrient maps. This script tests that directly.
#
# DESIGN (leave-one-country-out, estimand C, level target primary)
# ---------------------------------------------------------------
# For each held-out country H and each outcome B measured in H (the 22 LOCO
# cells of benchmarks_v2_cells.csv, estimand "country", arm "domain_index"),
# score B's districts in H with index weights from different sources:
#   own                      the index trained on B in the countries other than
#                            H. The standard LOCO domain_index; GUARD 1 below
#                            requires it to reproduce benchmarks_v2_cells.csv
#                            (country, domain_index, level) to 0.005 or the
#                            script stops before writing anything.
#   swap (one per A)         for every outcome A of a DIFFERENT nutrient that is
#                            measured in >= 2 countries other than H, the index
#                            trained on A in the countries other than H only
#                            (never H's own survey), applied to H's B districts
#                            and scored against B.
#   same_nutrient_other_pop  the same nutrient in the other population (e.g.
#                            child_iron for women_iron), same eligibility rule.
#                            NOT a swap; reported separately.
#   generic                  one "deficiency in general" index: the sum of the
#                            z-weight vectors of every eligible outcome (own,
#                            same-nutrient-other-population and all swaps; i.e.
#                            every outcome measured in >= 2 countries other
#                            than H), each z computed as in `own` (Fisher-z of
#                            Spearman x sqrt(n-3) on that outcome's training
#                            rows, outcome standardised within country). The z
#                            vectors can only be summed on one shared basis, so
#                            the generic's domain PCs are learned once on the
#                            stacked training rows of all eligible outcomes, on
#                            the columns common to all of them and to H.
# Nutrients: vitA, iron, folate, b12, zinc, selenium, iodine. In the current
# targets_v2.csv zinc is Malawi-only and selenium / iodine are absent, so none
# of them is ever eligible (>= 2 countries other than H); the eligible sources
# are child/women vitA, child/women iron, women folate and women B12.
# Every fit mirrors 02b_merge_and_loco.R estimand C: the same cell builder
# (within-country rank-normalisation of the cell's own districts), pooling on
# the columns common to the training countries and H, domain PCs learned on
# the training rows only (sign_rows), outcome standardised within country
# (logit for prevalence), scored within the held-out country. Tiers
# open,survey_public; V2_INDEX_SHRINK / V2_DOMAIN_REP / V2_PREP_SCALE /
# V2_KEEP_NATIONAL / V2_DROP_MODELLED unset so the protocol defaults apply.
#
# PRE-REGISTERED CLAIM (fixed before any result was seen)
# -------------------------------------------------------
#   Test       a cell is "nutrient-specific" iff  own - max(swap_A) >= 0.03
#              in level Spearman.
#   Also       own - mean(swap_A) and own - generic, per cell.
#   Summary    count of specific cells out of 22, and by nutrient
#              (iron 8 cells, vitA 8, folate 3, B12 3). B12 and iron are the
#              cases the talk leans on.
#   Verdict    the talk's nutrient table is DEFENDED for a nutrient if a
#              strict majority of its LOCO cells are nutrient-specific (by the
#              max test above). Otherwise the honest wording is "one general
#              deprivation map, which ranks this nutrient about as well as any
#              other".
#   Bias       a max over several noisy swaps is biased upward, so
#              own - mean(swap) is the fairer second reading. Both readings are
#              reported (the mean reading with the same >= 0.03 bar and the
#              same majority rule); neither is chosen after seeing them. The
#              verdict above is the max reading, as registered.
#   Secondary  the same on the prevalence target, reported, not used for the
#              verdict.
# Known asymmetry, stated in advance: a swap is trained on A's countries other
# than H, which can differ in number from B's (folate / B12 exist in 3
# countries, iron / vitA in 4). So for iron / vitA cells in Ghana, Malawi and
# Sierra Leone the folate / B12 swaps train on 2 countries against own's 3,
# and for folate / B12 cells the iron / vitA swaps train on 3 against own's 2.
# n_train_countries is recorded on every row.
#
#   Rscript -e "source('scripts/protocol_v2/76_nutrient_specificity_swap.R')"
# -> results/tables/protocol_v2/nx01_per_source.csv   cell x weight source
#    results/tables/protocol_v2/nx01_per_cell.csv     cell: own / max / mean swap / generic / same-nutrient, deltas, flags
#    results/tables/protocol_v2/nx01_summary.csv      counts and verdict by nutrient, per target
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
Sys.unsetenv(c("V2_INDEX_SHRINK", "V2_DOMAIN_REP", "V2_PREP_SCALE", "V2_KEEP_NATIONAL", "V2_DROP_MODELLED"))
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R")

OUTDIR    <- "results/tables/protocol_v2"
PASS_BAR  <- 0.03    # pre-registered: own - max(swap) >= 0.03
REPRO_TOL <- 0.005   # GUARD 1 tolerance against benchmarks_v2_cells.csv

TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
OUTCOMES <- unique(TG$outcome)

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

# identical to build_cell() in 02b_merge_and_loco.R
build_cell <- function(cn, on, target) {
  t <- TG[TG$country == cn & TG$outcome == on, ]
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
  m <- inner_join(t, sc, by = c("Admin1", "Admin2"))
  m <- inner_join(m, CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")],
                  by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 20) return(NULL)
  y_nat <- m[[ycol]]
  list(country = cn, n = nrow(m), y_nat = y_nat,
       y_mod = if (target == "prev") .v2_logit(y_nat) else y_nat,
       X = Xr, w = m[[wcol]])
}

nutrient_of   <- function(o) sub("^(child|women)_", "", o)
population_of <- function(o) sub("_.*$", "", o)

#' pooled index on the training cells, applied to the held-out cell (02b estimand C)
fit_apply <- function(train_cells, test_cell) {
  common <- Reduce(intersect, c(lapply(train_cells, function(z) colnames(z$X)), list(colnames(test_cell$X))))
  if (length(common) < 20) return(NULL)
  Xtr <- do.call(rbind, lapply(train_cells, function(z) z$X[, common, drop = FALSE]))
  ntr <- nrow(Xtr)
  if (ntr < 20) return(NULL)
  Xm <- rbind(Xtr, test_cell$X[, common, drop = FALSE])
  tr <- seq_len(ntr); te <- ntr + seq_len(test_cell$n)
  Y  <- c(unlist(lapply(train_cells, function(z) as.numeric(scale(z$y_mod)))), rep(NA_real_, test_cell$n))
  Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  p  <- ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, NULL)
  list(pred = p, n_common = length(common), n_axes = ncol(Dm), n_train_areas = ntr)
}

#' generic: sum of every eligible outcome's z vector on one shared basis
fit_generic <- function(src_cells, test_cell) {
  all_tr <- unlist(src_cells, recursive = FALSE)
  common <- Reduce(intersect, c(lapply(all_tr, function(z) colnames(z$X)), list(colnames(test_cell$X))))
  if (length(common) < 20) return(NULL)
  Xtr <- do.call(rbind, lapply(all_tr, function(z) z$X[, common, drop = FALSE]))
  ntr <- nrow(Xtr)
  Xm <- rbind(Xtr, test_cell$X[, common, drop = FALSE])
  tr <- seq_len(ntr); te <- ntr + seq_len(test_cell$n)
  Y  <- unlist(lapply(all_tr, function(z) as.numeric(scale(z$y_mod))))
  grp <- rep(rep(names(src_cells), vapply(src_cells, length, 0L)), vapply(all_tr, function(z) z$n, 0L))
  stopifnot(length(grp) == ntr, length(Y) == ntr)
  Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
  zsum <- rep(0, ncol(Dm))
  for (o in names(src_cells)) {
    k <- which(grp == o)
    zsum <- zsum + .index_weights_v2(Dm[k, , drop = FALSE], Y[k])
  }
  list(pred = as.numeric(Dm[te, , drop = FALSE] %*% zsum), n_common = length(common),
       n_axes = ncol(Dm), n_train_areas = ntr)
}

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

rows <- list(); n_fail <- 0L
for (target in c("level", "prev")) {
  CELLS <- list()
  for (on in OUTCOMES) for (cn in COUNTRIES) {
    cc <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
    if (!is.null(cc)) CELLS[[paste(cn, on, sep = "|")]] <- cc
  }
  have <- function(on) unname(COUNTRIES[paste(COUNTRIES, on, sep = "|") %in% names(CELLS)])
  loco_outcomes <- OUTCOMES[vapply(OUTCOMES, function(o) length(have(o)) >= 3, NA)]
  for (B in loco_outcomes) for (H in have(B)) {
    test <- CELLS[[paste(H, B, sep = "|")]]
    srcs <- OUTCOMES[vapply(OUTCOMES, function(o) length(setdiff(have(o), H)) >= 2, NA)]
    stopifnot(B %in% srcs)
    preds <- list(); src_cells <- list()
    for (A in srcs) {
      role <- if (A == B) "own" else if (nutrient_of(A) == nutrient_of(B)) "same_nutrient_other_pop" else "swap"
      trc_names <- setdiff(have(A), H)
      trc <- CELLS[paste(trc_names, A, sep = "|")]
      src_cells[[A]] <- trc
      f <- tryCatch(fit_apply(trc, test), error = function(e) { message("FIT ERROR ", target, " ", H, " ", B, " <- ", A, ": ", conditionMessage(e)); NULL })
      if (is.null(f) || length(f$pred) != test$n) { n_fail <- n_fail + 1L; next }
      preds[[A]] <- f$pred
      s <- score_v2(test$y_nat, f$pred, test$w, scale = if (target == "prev") "prev" else "level")
      rows[[length(rows) + 1L]] <- data.frame(
        target = target, heldout = H, outcome = B, nutrient = nutrient_of(B),
        source = A, role = role, n_train_countries = length(trc_names),
        train_countries = paste(trc_names, collapse = ";"), n_train_areas = f$n_train_areas,
        n_areas = test$n, n_common_cols = f$n_common, n_axes = f$n_axes,
        spearman = s$spearman, pearson = s$pearson, topk = s$topk, stringsAsFactors = FALSE)
    }
    g <- tryCatch(fit_generic(src_cells, test), error = function(e) { message("GENERIC ERROR ", target, " ", H, " ", B, ": ", conditionMessage(e)); NULL })
    if (is.null(g)) { n_fail <- n_fail + 1L } else {
      preds[["generic"]] <- g$pred
      s <- score_v2(test$y_nat, g$pred, test$w, scale = if (target == "prev") "prev" else "level")
      rows[[length(rows) + 1L]] <- data.frame(
        target = target, heldout = H, outcome = B, nutrient = nutrient_of(B),
        source = "generic", role = "generic", n_train_countries = NA_integer_,
        train_countries = paste(srcs, collapse = ";"), n_train_areas = g$n_train_areas,
        n_areas = test$n, n_common_cols = g$n_common, n_axes = g$n_axes,
        spearman = s$spearman, pearson = s$pearson, topk = s$topk, stringsAsFactors = FALSE)
    }
    # descriptive: how much each source's map agrees with own's map in H
    k0 <- length(rows) - length(preds) + 1L
    for (i in k0:length(rows)) {
      src <- rows[[i]]$source
      rows[[i]]$map_agreement_with_own <- if (!is.null(preds[[B]]) && !is.null(preds[[src]])) sp(preds[[src]], preds[[B]]) else NA_real_
    }
    own_i <- which(vapply(rows[k0:length(rows)], function(r) r$role == "own", NA))
    cat(sprintf("  %-5s %-12s %-12s own %.3f  sources %d\n", target, H, B,
                if (length(own_i)) rows[[k0 + own_i[1] - 1L]]$spearman else NA_real_, length(preds)))
  }
}
PS <- bind_rows(rows)
if (n_fail > 0) stop(sprintf("%d fits failed; nothing written", n_fail))

# ── GUARD 1: own reproduces benchmarks_v2_cells.csv (country, domain_index) ──
BM <- read.csv(file.path(OUTDIR, "benchmarks_v2_cells.csv"), stringsAsFactors = FALSE) |>
  filter(estimand == "country", arm == "domain_index") |>
  select(target, heldout = country, outcome, bench_spearman = spearman)
REPRO <- PS |> filter(role == "own") |> select(target, heldout, outcome, own = spearman) |>
  full_join(BM, by = c("target", "heldout", "outcome")) |> mutate(repro_diff = own - bench_spearman)
rl <- REPRO |> filter(target == "level"); rp <- REPRO |> filter(target == "prev")
cat(sprintf("\nGUARD 1 level: %d cells (benchmark %d), max |own - benchmark| = %.2e\n",
            sum(is.finite(rl$own)), sum(is.finite(rl$bench_spearman)), max(abs(rl$repro_diff))))
cat(sprintf("GUARD 1 prev : %d cells (benchmark %d), max |own - benchmark| = %.2e (reported, not a stop rule)\n",
            sum(is.finite(rp$own)), sum(is.finite(rp$bench_spearman)), max(abs(rp$repro_diff))))
if (nrow(rl) != 22 || anyNA(rl$repro_diff) || max(abs(rl$repro_diff)) > REPRO_TOL) {
  print(as.data.frame(rl), digits = 3)
  stop("GUARD 1 failed: own does not reproduce benchmarks_v2_cells.csv (country, domain_index, level) to ", REPRO_TOL)
}

# ── per cell ────────────────────────────────────────────────────────────────
PC <- PS |> group_by(target, heldout, outcome, nutrient) |>
  summarise(n_areas = first(n_areas),
            own = spearman[role == "own"][1],
            max_swap = if (any(role == "swap")) max(spearman[role == "swap"], na.rm = TRUE) else NA_real_,
            max_swap_source = if (any(role == "swap")) source[role == "swap"][which.max(spearman[role == "swap"])] else NA_character_,
            mean_swap = mean(spearman[role == "swap"], na.rm = TRUE),
            n_swaps = sum(role == "swap"),
            generic = spearman[role == "generic"][1],
            same_nutrient_other_pop = if (any(role == "same_nutrient_other_pop")) spearman[role == "same_nutrient_other_pop"][1] else NA_real_,
            mean_map_agreement_swap = mean(map_agreement_with_own[role == "swap"], na.rm = TRUE),
            .groups = "drop") |>
  mutate(own_minus_max_swap = own - max_swap,
         own_minus_mean_swap = own - mean_swap,
         own_minus_generic = own - generic,
         specific = own_minus_max_swap >= PASS_BAR,
         specific_mean_reading = own_minus_mean_swap >= PASS_BAR) |>
  left_join(REPRO |> select(target, heldout, outcome, bench_spearman, repro_diff),
            by = c("target", "heldout", "outcome")) |>
  arrange(target, nutrient, outcome, heldout)

# ── summary by nutrient ─────────────────────────────────────────────────────
summ_one <- function(d, label) {
  n <- nrow(d); k <- sum(d$specific); km <- sum(d$specific_mean_reading)
  data.frame(nutrient = label, cells = n, specific = k, specific_mean_reading = km,
             defended = k > n / 2, defended_mean_reading = km > n / 2,
             mean_own = mean(d$own), mean_max_swap = mean(d$max_swap), mean_mean_swap = mean(d$mean_swap),
             mean_generic = mean(d$generic), mean_same_nutrient_other_pop = mean(d$same_nutrient_other_pop, na.rm = TRUE),
             mean_own_minus_max = mean(d$own_minus_max_swap), mean_own_minus_mean = mean(d$own_minus_mean_swap),
             mean_own_minus_generic = mean(d$own_minus_generic),
             mean_map_agreement_swap = mean(d$mean_map_agreement_swap, na.rm = TRUE),
             stringsAsFactors = FALSE)
}
SUMM <- bind_rows(lapply(c("level", "prev"), function(tg) {
  d <- PC[PC$target == tg, ]
  out <- bind_rows(lapply(c("iron", "vitA", "folate", "b12"), function(nu) summ_one(d[d$nutrient == nu, ], nu)),
                   summ_one(d, "all"))
  out$target <- tg
  out$role <- if (tg == "level") "primary (verdict)" else "secondary (not used for verdict)"
  out$verdict <- ifelse(out$nutrient == "all", NA_character_,
    ifelse(out$defended, "nutrient table defended",
           "one general deprivation map, which ranks this nutrient about as well as any other"))
  out
})) |> select(target, role, nutrient, everything())

write.csv(PS,   file.path(OUTDIR, "nx01_per_source.csv"), row.names = FALSE)
write.csv(PC,   file.path(OUTDIR, "nx01_per_cell.csv"),   row.names = FALSE)
write.csv(SUMM, file.path(OUTDIR, "nx01_summary.csv"),    row.names = FALSE)

cat("\n=== NX-01 per cell (level; specific = own - max swap >= 0.03) ===\n")
print(as.data.frame(PC |> filter(target == "level") |>
  select(nutrient, outcome, heldout, n_areas, own, max_swap, max_swap_source, mean_swap, generic,
         same_nutrient_other_pop, own_minus_max_swap, own_minus_mean_swap, own_minus_generic, specific,
         mean_map_agreement_swap) |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== NX-01 per cell (prev, secondary) ===\n")
print(as.data.frame(PC |> filter(target == "prev") |>
  select(nutrient, outcome, heldout, own, max_swap, mean_swap, generic, same_nutrient_other_pop,
         own_minus_max_swap, own_minus_mean_swap, specific) |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\n=== NX-01 summary ===\n")
print(as.data.frame(SUMM |> select(target, nutrient, cells, specific, specific_mean_reading, defended,
  defended_mean_reading, mean_own, mean_max_swap, mean_mean_swap, mean_generic, mean_same_nutrient_other_pop,
  mean_map_agreement_swap) |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nDONE\n")
