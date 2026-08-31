# =============================================================================
# scripts/accuracy_impact/ws2_anchor_controls.R
#
# WS2. Controls and circularity tests for the Admin-1 anchor.
#
# WHAT SECTION 4 CLAIMS, AND WHAT IT DOES NOT TEST
# ------------------------------------------------
# Anchoring an Admin-2 model to its region's design-based survey estimate takes
# mean r from 0.164 to 0.413 (source:
# results/tables/frozen_2026-09/admin1_arms.csv). Section 4.4 argues this is not
# circular, on the grounds that the anchor supplies only the region's mean and
# says nothing about which districts inside it are worse.
#
# The argument has a gap. The regional anchor is computed from ALL respondents
# in the region, INCLUDING the district being scored. With three to four
# districts per region in two of these countries, a district contributes roughly
# a quarter to a third of its own anchor, so each district's own sampling error
# is imported into the number used to correct it. That is a real channel and
# Section 4.4 does not close it.
#
# Section 12 proposes checking whether the gain survives a within-region
# correlation analysis. That check is vacuous: the shift is constant within a
# region and monotone, so within-region rank correlation is unchanged BY
# CONSTRUCTION. The tests below replace it.
#
# THE ARMS
#   flat regional mean     every district predicted by its region's survey
#                          estimate, no covariates at all. If the anchored
#                          covariate model does not beat this, the survey's
#                          regional estimates carry the information and the
#                          covariate pattern adds little.
#   flat national mean     the same at national resolution.
#   hard / shrunk anchor   reproductions of the published arms.
#   jackknife anchor       each district's regional anchor is computed from the
#                          region's OTHER districts only, so no district
#                          contributes to its own correction.
#   split-sample anchor    within each region clusters are split in two; the
#                          anchor comes from one half and the scoring from the
#                          other. Repeated over splits.
#
# All arms are scored on the identical cells and the identical leave-one-region-
# out folds as the published table, so every comparison is like for like.
#
# MAE and signed bias are reported as the primary metrics, per WS2d. Section 4.4
# itself concedes that r is the wrong summary for an arm that moves only the
# level.
#
#   PROFILE=smoke   Ghana only
#   Rscript scripts/accuracy_impact/ws2_anchor_controls.R
# -> results/tables/anchor_controls.csv
# -> results/tables/anchor_implied_shifts.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(readr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
SEED <- 20260905L; N_SPLIT <- 25L; K_SCREEN <- 20L
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]

# Scoring is identical to scripts/covariates/08_admin1_arms.R except that
# r_share now uses the EMPIRICAL ceiling and applies the r_max > 0.05 guard that
# add_reliability_columns() uses. The published script divided by an analytic
# r_max as small as 0.030, which is how a cell reached r_share 13.61.
REL <- tryCatch(read.csv(here("results", "tables", "reliability_empirical.csv"),
                         stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(REL)) REL <- REL[REL$scheme == "within", ]
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))

score <- function(svy, pred_df, col, country, outcome, arm) {
  jb <- admin2_join_by(svy, pred_df)
  keep_s <- intersect(c("Admin1", "Admin2", "svy_prev", "n_svy"), names(svy))
  keep_p <- intersect(c("Admin1", "Admin2", col), names(pred_df))
  m <- dplyr::inner_join(svy[, keep_s, drop = FALSE],
                         pred_df[, keep_p, drop = FALSE], by = jb)
  names(m)[names(m) == col] <- "p"
  m <- m[is.finite(m$svy_prev) & is.finite(m$p), , drop = FALSE]
  if (nrow(m) < 5) return(NULL)
  r <- suppressWarnings(stats::cor(m$svy_prev, m$p))
  rme <- if (!is.null(REL)) {
    i <- which(kk(REL$country) == kk(country) & REL$outcome == outcome)
    if (length(i)) REL$r_max_emp[i[1]] else NA_real_
  } else NA_real_
  data.frame(country = country, outcome = outcome, arm = arm, n_admin2 = nrow(m),
             pearson_r = round(r, 4),
             mae_pp = round(100 * mean(abs(m$svy_prev - m$p)), 2),
             bias_pp = round(100 * mean(m$p - m$svy_prev), 2),
             abs_bias_pp = round(abs(100 * mean(m$p - m$svy_prev)), 2),
             r_max_emp = round(rme, 4),
             r_share_emp = if (is.finite(rme) && rme > 0.05) round(r / rme, 3) else NA_real_,
             stringsAsFactors = FALSE)
}

# Design-based regional prevalence, optionally excluding a set of districts.
# This is admin1_design_based() with a hold-out, which is what the jackknife
# needs and what the published arms never computed.
a1_prev_excluding <- function(d, cc, oc, drop_admin2 = character(0)) {
  y <- suppressWarnings(as.numeric(d[[oc$binary]]))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  g  <- trimws(as.character(d[[cc$admin1_col]]))
  a2 <- trimws(as.character(d[[cc$admin2_col]]))
  ok <- is.finite(y) & !is.na(g) & !(a2 %in% drop_admin2)
  if (sum(ok) < 10) return(NULL)
  y <- y[ok]; w <- w[ok]; g <- g[ok]
  rg <- sort(unique(g))
  data.frame(Admin1 = rg,
             prev = vapply(rg, function(r) stats::weighted.mean(y[g == r], w[g == r]), numeric(1)),
             n    = vapply(rg, function(r) sum(g == r), numeric(1)),
             stringsAsFactors = FALSE)
}

rows <- list(); shifts <- list()
set.seed(SEED)
for (ctry in names(cfgs)) {
  cc <- cfgs[[ctry]]
  key <- H[H$country == ctry, c("Admin1", "Admin2")]
  key <- key[!duplicated(key[, c("Admin1", "Admin2")]), ]
  if (!nrow(key)) next
  hc <- H[H$country == ctry, ]
  pop <- data.frame(Admin1 = hc$Admin1, Admin2 = hc$Admin2,
                    pop = if ("ghs_pop" %in% names(hc)) hc$ghs_pop else NA_real_)

  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    sfx <- paste0(tolower(ctry), "_", ocn)
    svy <- rd(paste0("svy_admin2_", sfx)); od <- rd(paste0("outcome_data_", sfx))
    if (is.null(svy) || is.null(od) || nrow(svy) < 5) next
    d <- od$data
    if (!all(c(cc$admin1_col, cc$admin2_col, oc$binary, cc$weight_col) %in% names(d))) next

    # ── the published base model: admin-2 fit, leave-one-region-out ──────────
    tr <- dplyr::inner_join(
      svy[, intersect(c("Admin1", "Admin2", "svy_prev", "n_svy"), names(svy)), drop = FALSE],
      hc, by = admin2_join_by(svy, hc))
    tr <- tr[is.finite(tr$svy_prev), , drop = FALSE]
    if (nrow(tr) < 10 || dplyr::n_distinct(tr$Admin1) < 3) next
    Xtr <- as.matrix(tr[, COVS, drop = FALSE])
    oof <- rep(NA_real_, nrow(tr))
    for (r in unique(tr$Admin1)) {
      i <- which(tr$Admin1 == r)
      m <- .ds_fit(Xtr[-i, , drop = FALSE], tr$svy_prev[-i], k_screen = K_SCREEN)
      p <- .ds_predict(m, Xtr[i, , drop = FALSE])
      if (!is.null(p)) oof[i] <- p
    }
    p2 <- data.frame(Admin2 = tr$Admin2, Admin1 = tr$Admin1, pred = oof)
    p2 <- p2[is.finite(p2$pred), , drop = FALSE]
    if (nrow(p2) < 5) next
    message("=== ", ctry, " / ", ocn, " (", nrow(p2), " districts) ===")

    nat <- national_design_based(od, cc, oc)
    a1t <- admin1_design_based(od, cc, oc)
    if (is.null(nat) || is.null(a1t)) next

    add <- function(df, col, arm) {
      r <- score(svy, df, col, ctry, ocn, arm)
      if (!is.null(r)) rows[[length(rows) + 1L]] <<- r
    }

    add(p2, "pred", "1 no anchor (LORO)")

    # ── 2a. anchor-only arms: no covariates at all ──────────────────────────
    fr <- data.frame(Admin2 = p2$Admin2, Admin1 = p2$Admin1,
                     pred = a1t$prev[match(p2$Admin1, a1t$Admin1)])
    add(fr, "pred", "2a flat REGIONAL mean (no covariates)")
    fn <- data.frame(Admin2 = p2$Admin2, Admin1 = p2$Admin1,
                     pred = rep(nat$prev, nrow(p2)))
    add(fn, "pred", "2a flat NATIONAL mean (no covariates)")

    # ── published anchors, reproduced ───────────────────────────────────────
    a1h <- a1t; a1h$target <- a1h$prev
    ph <- benchmark_admin2_to_admin1(p2, "pred", a1h, national = nat$prev,
                                     admin1_map = key, pop = pop)
    add(ph, "pred", "3 ADMIN-1 anchor (hard)")
    lg <- attr(ph, "benchmark_a1")
    if (!is.null(lg) && nrow(lg)) {
      lg$country <- ctry; lg$outcome <- ocn; lg$arm <- "hard"
      shifts[[length(shifts) + 1L]] <- lg
    }
    a1s <- shrink_admin1_targets(a1t, nat$prev)
    ps <- benchmark_admin2_to_admin1(p2, "pred", a1s, national = nat$prev,
                                     admin1_map = key, pop = pop)
    add(ps, "pred", "3 ADMIN-1 anchor (shrunk)")

    pn <- benchmark_admin2_table(p2, "pred", nat$prev, pop)
    add(pn, "pred", "4 national anchor")
    # 2c: the shift the national anchor applies, and the un-anchored national
    # aggregate it is correcting. Section 3 reports the pipeline recovers the
    # national prevalence to 0.96 pp, so a large shift here needs explaining.
    wpop <- pop$pop[match(paste0(p2$Admin1, "\r", p2$Admin2),
                          paste0(pop$Admin1, "\r", pop$Admin2))]
    agg_un <- area_aggregate(p2$pred, wpop)
    shifts[[length(shifts) + 1L]] <- data.frame(
      Admin1 = "(national)", n_areas = nrow(p2), target = nat$prev,
      delta = NA_real_, before = agg_un, after = nat$prev, used = "national",
      country = ctry, outcome = ocn, arm = "national", stringsAsFactors = FALSE)

    # ── 2b. jackknife anchor: a district never contributes to its own ───────
    pj <- p2; pj$pred <- NA_real_
    for (i in seq_len(nrow(p2))) {
      a1x <- a1_prev_excluding(d, cc, oc, drop_admin2 = p2$Admin2[i])
      if (is.null(a1x)) next
      reg <- p2$Admin1[i]
      j <- which(p2$Admin1 == reg & p2$Admin2 != p2$Admin2[i])
      t_r <- a1x$prev[match(reg, a1x$Admin1)]
      if (!length(j) || !is.finite(t_r)) next
      # Solve the shift on the region's OTHER districts, then apply it to this
      # one. No quantity here has seen district i.
      b <- benchmark_area_predictions(p2$pred[j], t_r, wpop[j], method = "logit_shift")
      if (!is.finite(b$delta)) next
      pj$pred[i] <- stats::plogis(stats::qlogis(
        pmin(pmax(p2$pred[i], 1e-4), 1 - 1e-4)) + b$delta)
    }
    add(pj, "pred", "5 ADMIN-1 anchor (hard, JACKKNIFE)")

    # ── 2b. split-sample anchor: anchor on half the clusters, score on all ──
    clv <- trimws(as.character(d[[cc$cluster_id]]))
    acc <- matrix(NA_real_, nrow(p2), N_SPLIT)
    for (s in seq_len(N_SPLIT)) {
      cu <- unique(clv[!is.na(clv)])
      half <- sample(rep_len(c(TRUE, FALSE), length(cu))); names(half) <- cu
      dh <- d[!is.na(clv) & half[clv], , drop = FALSE]
      a1x <- a1_prev_excluding(dh, cc, oc)
      if (is.null(a1x)) next
      pp <- p2
      a1x$target <- a1x$prev
      pb <- benchmark_admin2_to_admin1(pp, "pred", a1x, national = nat$prev,
                                       admin1_map = key, pop = pop)
      acc[, s] <- pb$pred
    }
    pss <- p2; pss$pred <- rowMeans(acc, na.rm = TRUE)
    add(pss, "pred", "6 ADMIN-1 anchor (hard, SPLIT-SAMPLE)")
  }
}

res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
readr::write_csv(res, here("results", "tables", sprintf("anchor_controls%s.csv", SUF)))
sh <- dplyr::bind_rows(shifts)
readr::write_csv(sh, here("results", "tables", sprintf("anchor_implied_shifts%s.csv", SUF)))

cat("\n=== WS2: arms scored on identical cells and folds ===\n")
s <- res |> group_by(arm) |>
  summarise(cells = dplyr::n(),
            mean_r = round(mean(pearson_r, na.rm = TRUE), 3),
            med_r  = round(stats::median(pearson_r, na.rm = TRUE), 3),
            MAE_pp = round(mean(mae_pp, na.rm = TRUE), 2),
            mean_abs_bias = round(mean(abs_bias_pp, na.rm = TRUE), 2),
            mean_signed_bias = round(mean(bias_pp, na.rm = TRUE), 2),
            med_r_share_emp = round(stats::median(r_share_emp, na.rm = TRUE), 2),
            .groups = "drop") |> arrange(arm)
print(as.data.frame(s), row.names = FALSE)

cat("\n=== paired against 'no anchor', same cells ===\n")
w <- res |> select(country, outcome, arm, pearson_r, mae_pp, abs_bias_pp) |>
  tidyr::pivot_longer(c(pearson_r, mae_pp, abs_bias_pp))
base <- res[res$arm == "1 no anchor (LORO)", c("country","outcome","pearson_r","mae_pp","abs_bias_pp")]
names(base)[3:5] <- c("r0","mae0","bias0")
cmp <- dplyr::inner_join(res, base, by = c("country","outcome"))
p <- cmp |> filter(arm != "1 no anchor (LORO)") |> group_by(arm) |>
  summarise(cells = dplyr::n(),
            d_r = round(mean(pearson_r - r0, na.rm = TRUE), 3),
            better_r = sprintf("%d/%d", sum(pearson_r > r0, na.rm = TRUE), dplyr::n()),
            d_MAE = round(mean(mae_pp - mae0, na.rm = TRUE), 2),
            better_MAE = sprintf("%d/%d", sum(mae_pp < mae0, na.rm = TRUE), dplyr::n()),
            d_absbias = round(mean(abs_bias_pp - bias0, na.rm = TRUE), 2),
            .groups = "drop")
print(as.data.frame(p), row.names = FALSE)

cat("\n=== 2c: implied national shift, un-anchored aggregate against the survey ===\n")
nsh <- sh[sh$arm == "national", ]
nsh$gap_pp <- round(100 * (nsh$before - nsh$target), 2)
print(as.data.frame(nsh[, c("country","outcome","before","target","gap_pp")]), row.names = FALSE)
cat(sprintf("\nmean |gap| between the un-anchored population-weighted national aggregate\nand the design-based survey national estimate: %.2f pp over %d cells\n",
            mean(abs(nsh$gap_pp), na.rm = TRUE), nrow(nsh)))
