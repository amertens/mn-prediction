# =============================================================================
# scripts/accuracy_impact/wsc3_iron_bridge.R
#
# WS-C3. The iron-anaemia bridge, and the transport test that matters.
#
# WHY THIS MODEL AND NOT ANOTHER
# ------------------------------
# Every covariate result in this project so far has been a socioeconomic
# gradient wearing a nutrition label: WS4b measured education at 0.679 to 0.808
# and stunting at -0.122 from the same predictors. The bridge is the one
# candidate with a physiological mechanism rather than a correlational one.
# Haemoglobin responds to iron status; malaria and the inherited haemoglobin
# disorders are the two large confounders of that relationship at the population
# level. A model built from those three, and nothing else, is a mechanism test.
#
# It is also the transport test that matters. A socioeconomic proxy cannot be
# expected to carry across countries with different economies. A physiological
# relationship should, and leave-one-country-out is where that is decided.
#
# WHERE THE COVARIATES COME FROM, AND WHY NOT FROM THE 294
# --------------------------------------------------------
# Stage 0 established that the harmonised 294 contains no haemoglobin and no
# sickle-cell layers. They exist in two other places and both are used here:
#
#   haemoglobin, anaemia, iron in pregnancy   data/covariates/nutrition_proximal.csv
#                                             (WS-F2, from the prior DHS round)
#   MAP malaria and inherited disorders       the individual-level merge, where
#                                             they are district-level values
#                                             replicated across respondents, so
#                                             one value per district recovers the
#                                             Admin-2 layer
#
# The prior-round DHS is a DIFFERENT SURVEY from the biomarker survey, so a
# haemoglobin predictor here is not the outcome measured twice. It is also the
# deployable case: a country contemplating a micronutrient survey usually has a
# recent DHS with haemoglobin in it.
#
#   PROFILE=smoke   Ghana only
#   Rscript scripts/accuracy_impact/wsc3_iron_bridge.R
# -> results/tables/iron_bridge.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
SEED <- 20260921L; K_SCREEN <- 9L
OUTCOMES <- c("child_iron", "women_iron")
set.seed(SEED)
nz <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS294 <- setdiff(names(H), c("country", "Admin1", "Admin2"))
NP <- suppressMessages(readr::read_csv(
  here("data", "covariates", "nutrition_proximal.csv"), show_col_types = FALSE)) |>
  as.data.frame()
NPC <- grep("^np_", names(NP), value = TRUE); NPC <- NPC[!grepl("_n$", NPC)]

# THE BRIDGE IS PRE-SPECIFIED AND SMALL, WHICH IS THE POINT.
#
# One column per mechanistic concept: haemoglobin and anaemia in each
# population, iron supplementation in pregnancy, malaria parasite rate, the two
# inherited haemoglobin disorders that confound anaemia in West Africa, and
# bednet coverage. The first version of this script took every MAP layer that
# passed a filter and arrived at 31 covariates including Duffy negativity,
# G6PD deficiency, temperature suitability and non-malarial fever. A 31-column
# model with a top-8 screen is an exploratory model, not the mechanism test the
# workstream specifies, so the set is now named rather than filtered.
MAP_CORE <- c(
  "map2_explorer_2020_global_pfpr",
  "map2_blood_disorders_201201_global_sickle_haemoglobin_hbs_al",
  "map2_blood_disorders_201201_africa_hbc_allele_frequency")
NP_CORE <- c("np_child_hb_gdl", "np_women_hb_gdl", "np_child_anaemia_any",
             "np_women_anaemia_any", "np_iron_pregnancy")
ITN_CORE <- c("dhs_hh_itn_any")

# ── recover the MAP Admin-2 layers from the individual-level merge ──────────
# The layers are extracted per respondent, at the district in some countries and
# at the cluster in others, so the district MEAN is taken. How many columns
# actually varied within a district is counted and reported rather than assumed.
map_admin2 <- function(cn, cc) {
  ocn <- names(cc$outcomes)[1]
  od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                       store = STORE), error = function(e) NULL)
  if (is.null(od)) return(NULL)
  d <- od$data
  full <- od$Xvars_full %||% od$Xvars
  cand <- intersect(intersect(MAP_CORE, full), names(d))
  if (!length(cand)) return(NULL)
  a1 <- trimws(as.character(d[[cc$admin1_col]])); a2 <- trimws(as.character(d[[cc$admin2_col]]))
  ok <- !is.na(a1) & !is.na(a2)
  key <- paste(a1, a2, sep = "|::|")[ok]
  out <- unique(data.frame(country = cc$country, Admin1 = a1[ok], Admin2 = a2[ok],
                           stringsAsFactors = FALSE))
  kept <- character(0); nvary <- 0L
  for (v in cand) {
    x <- nz(d[[v]])[ok]
    # DISTRICT MEAN, not a constancy requirement. The first version of this
    # function accepted a column only where it was constant within every
    # district, which silently kept 22 columns for three countries and ONE for
    # Malawi, because Malawi's layers are extracted at the cluster and vary
    # inside a district. Averaging is correct either way: for a district-level
    # layer the mean is the value, and for a cluster-level layer it is the
    # district aggregate.
    nd <- tapply(x, key, function(z) length(unique(z[is.finite(z)])))
    if (max(nd, na.rm = TRUE) > 1) nvary <- nvary + 1L
    val <- tapply(x, key, function(z) { z <- z[is.finite(z)]
      if (length(z)) mean(z) else NA_real_ })
    out[[v]] <- as.numeric(val[paste(out$Admin1, out$Admin2, sep = "|::|")])
    if (mean(is.finite(out[[v]])) > 0.5) kept <- c(kept, v) else out[[v]] <- NULL
  }
  if (!length(kept)) return(NULL)
  attr(out, "cols") <- kept; attr(out, "n_vary_within_district") <- nvary
  out
}

cfgs <- get_country_configs()
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]
frames <- list(); map_cols <- character(0)
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]
  mp <- map_admin2(cn, cc)
  if (!is.null(mp)) map_cols <- union(map_cols, attr(mp, "cols"))
  for (ocn in intersect(OUTCOMES, names(cc$outcomes))) {
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(sv) || nrow(sv) < 10) next
    hc <- H[H$country == cn, , drop = FALSE]
    m <- dplyr::inner_join(
      sv[, intersect(c("Admin1","Admin2","svy_prev","n_svy"), names(sv)), drop = FALSE],
      hc, by = admin2_join_by(sv, hc))
    m <- m[is.finite(m$svy_prev), , drop = FALSE]
    if (!nrow(m)) next
    npc <- NP[NP$country == cc$country, , drop = FALSE]
    n_before <- nrow(m)
    m <- dplyr::left_join(m, npc[, c("Admin1","Admin2", NPC)],
                          by = admin2_join_by(m, npc))
    if (!is.null(mp)) m <- dplyr::left_join(m, mp[, c("Admin1","Admin2", attr(mp,"cols"))],
                                            by = c("Admin1","Admin2"))
    m$country <- cc$country; m$outcome <- ocn
    frames[[paste(cn, ocn)]] <- m
    cat(sprintf("  %-13s %-11s districts=%3d  np matched=%3d (%.0f%%)  map cols=%d\n",
                cn, ocn, nrow(m), sum(is.finite(m[[NPC[1]]])),
                100 * mean(is.finite(m[[NPC[1]]])),
                if (is.null(mp)) 0L else length(attr(mp, "cols"))))
  }
}
if (!length(frames)) stop("No cells assembled.")
BRIDGE <- unique(c(NP_CORE, map_cols, intersect(ITN_CORE, COVS294)))
cat(sprintf("\nbridge covariates (%d): %s\n", length(BRIDGE), paste(BRIDGE, collapse = ", ")))

fit_score <- function(m, cols, arm, kscr) {
  cols <- intersect(cols, names(m))
  X <- as.matrix(m[, cols, drop = FALSE]); X[!is.finite(X)] <- NA
  ok <- colMeans(is.finite(X)) > 0.6
  X <- X[, ok, drop = FALSE]
  if (ncol(X) < 2) return(NULL)
  for (j in seq_len(ncol(X))) { md <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(md)) md else 0 }
  sdok <- apply(X, 2, stats::sd) > 1e-12; X <- X[, sdok, drop = FALSE]
  if (ncol(X) < 2 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  oof <- rep(NA_real_, nrow(m))
  for (rg in unique(m$Admin1)) {
    i <- which(m$Admin1 == rg); tr <- setdiff(seq_len(nrow(m)), i)
    if (length(tr) < 6) next
    f <- .ds_fit(X[tr, , drop = FALSE], m$svy_prev[tr], k_screen = min(kscr, ncol(X)))
    p <- .ds_predict(f, X[i, , drop = FALSE]); if (!is.null(p)) oof[i] <- p
  }
  fin <- is.finite(oof)
  if (sum(fin) < 8 || stats::sd(m$svy_prev[fin]) == 0) return(NULL)
  data.frame(arm = arm, eval = "LORO", n_units = sum(fin), n_pred = ncol(X),
             r = round(suppressWarnings(stats::cor(m$svy_prev[fin], oof[fin])), 4),
             mae_pp = round(100 * mean(abs(oof[fin] - m$svy_prev[fin])), 2),
             bias_pp = round(100 * mean(oof[fin] - m$svy_prev[fin]), 2),
             stringsAsFactors = FALSE)
}
# The jackknifed regional mean, the WS-B1 baseline, on the identical rows.
flat_jk <- function(m) {
  p <- rep(NA_real_, nrow(m))
  for (i in seq_len(nrow(m))) {
    j <- which(m$Admin1 == m$Admin1[i] & seq_len(nrow(m)) != i)
    if (length(j) >= 2) p[i] <- stats::weighted.mean(m$svy_prev[j], pmax(m$n_svy[j], 1))
  }
  fin <- is.finite(p)
  if (sum(fin) < 8) return(NULL)
  data.frame(arm = "flat regional mean (jackknife)", eval = "LORO", n_units = sum(fin),
             n_pred = 0L,
             r = round(suppressWarnings(stats::cor(m$svy_prev[fin], p[fin])), 4),
             mae_pp = round(100 * mean(abs(p[fin] - m$svy_prev[fin])), 2),
             bias_pp = round(100 * mean(p[fin] - m$svy_prev[fin]), 2),
             stringsAsFactors = FALSE)
}

rows <- list()
for (nm in names(frames)) {
  m <- frames[[nm]]
  for (a in list(list("bridge", BRIDGE, K_SCREEN),
                 list("full 294", COVS294, 20L))) {
    r <- fit_score(m, a[[2]], a[[1]], a[[3]])
    if (!is.null(r)) { r$country <- m$country[1]; r$outcome <- m$outcome[1]
      rows[[length(rows)+1L]] <- r }
  }
  r <- flat_jk(m)
  if (!is.null(r)) { r$country <- m$country[1]; r$outcome <- m$outcome[1]
    rows[[length(rows)+1L]] <- r }
}

# ── LOCO: fit on three countries, predict the fourth ────────────────────────
for (ocn in OUTCOMES) {
  fs <- frames[grepl(paste0(" ", ocn, "$"), names(frames))]
  if (length(fs) < 3) next
  for (a in list(list("bridge", BRIDGE, K_SCREEN), list("full 294", COVS294, 20L))) {
    cols <- Reduce(intersect, c(list(a[[2]]), lapply(fs, names)))
    for (held in names(fs)) {
      tr <- dplyr::bind_rows(fs[setdiff(names(fs), held)])
      te <- fs[[held]]
      X <- as.matrix(tr[, cols, drop = FALSE]); Z <- as.matrix(te[, cols, drop = FALSE])
      ok <- colMeans(is.finite(X)) > 0.6 & colMeans(is.finite(Z)) > 0.6
      X <- X[, ok, drop = FALSE]; Z <- Z[, ok, drop = FALSE]
      if (ncol(X) < 2) next
      for (j in seq_len(ncol(X))) { md <- stats::median(X[, j], na.rm = TRUE)
        X[!is.finite(X[, j]), j] <- if (is.finite(md)) md else 0
        Z[!is.finite(Z[, j]), j] <- if (is.finite(md)) md else 0 }
      sdok <- apply(X, 2, stats::sd) > 1e-12
      X <- X[, sdok, drop = FALSE]; Z <- Z[, sdok, drop = FALSE]
      if (ncol(X) < 2) next
      f <- .ds_fit(X, tr$svy_prev, k_screen = min(a[[3]], ncol(X)))
      p <- .ds_predict(f, Z); if (is.null(p)) next
      fin <- is.finite(p) & is.finite(te$svy_prev)
      if (sum(fin) < 8 || stats::sd(te$svy_prev[fin]) == 0) next
      rows[[length(rows)+1L]] <- data.frame(
        arm = a[[1]], eval = "LOCO", n_units = sum(fin), n_pred = ncol(X),
        r = round(suppressWarnings(stats::cor(te$svy_prev[fin], p[fin])), 4),
        mae_pp = round(100 * mean(abs(p[fin] - te$svy_prev[fin])), 2),
        bias_pp = round(100 * mean(p[fin] - te$svy_prev[fin]), 2),
        country = te$country[1], outcome = ocn, stringsAsFactors = FALSE)
    }
  }
}
res <- dplyr::bind_rows(rows)
front <- c("country","outcome","eval","arm")
res <- res[, c(front, setdiff(names(res), front))]
readr::write_csv(res, here("results","tables", sprintf("iron_bridge%s.csv", SUF)))

cat("\n=== WS-C3: iron-anaemia bridge ===\n")
print(as.data.frame(res |> group_by(eval, arm) |>
  summarise(cells = dplyr::n(), mean_r = round(mean(r, na.rm=TRUE),3),
            med_r = round(stats::median(r, na.rm=TRUE),3),
            mae = round(mean(mae_pp, na.rm=TRUE),2),
            med_np = stats::median(n_pred), .groups="drop")), row.names = FALSE)
cat("\n=== per cell ===\n")
print(as.data.frame(res |> arrange(eval, country, outcome, arm) |>
  select(eval, country, outcome, arm, n_units, n_pred, r, mae_pp)), row.names = FALSE)
