# =============================================================================
# scripts/protocol_v2/78_flat_national_comparator.R   [LV-02, 2026-09-28]
#
# DOES THE TRANSPORTED RANKING ADD TO LEVEL ACCURACY, OR IS THE NATIONAL
# FIGURE DOING ALL THE WORK?
#
# AR-01 (script 35) scores "national anchor + transported ranking" (A1) but
# has no design that hands every district the national anchor itself. So the
# ranking's own contribution to level error was never measured. This adds the
# flat comparator and nothing else.
#
# (a) INTERNAL, the 22 AR-01 cells. Script 35's machinery is COPIED below
#     (lines marked "verbatim from 35"); script 35 itself is not edited or run.
#     A0  flat national anchor: every district gets a_nat, the SAME draw A1
#         uses, at the same fraction f, same units, same scoring (unweighted
#         district MAE in points against the full survey's district prevalence).
#     The random stream is consumed exactly as in script 35 (all four of its
#     designs are still drawn, in order), so A1 here reproduces
#     anchor_and_rank.csv draw for draw; A0 adds no draws. The script stops if
#     the A1 summary is more than 0.05 MAE from anchor_and_rank_summary.csv.
#
# (b) EXTERNAL, the WHO VMNIS sub-national cells of XV-01/02 (prevalence only:
#     21 African cells on iSDAsoil, 8 South Asian cells on SoilGrids).
#     z per admin-1 unit: xv_transport.csv keeps only per-cell Spearman, not
#     per-unit predictions, so the domain_index predictions are REBUILT here
#     with script 04's own code and inputs (copied, marked "verbatim from 04")
#     and checked against xv_transport.csv: every cell's Spearman must match to
#     1e-6 or the script stops. z = prediction standardised within the country.
#     A1_u = expit(logit(anchor) + rho_train * sd_train * z_u);  A0_u = anchor.
#     rho_train, pre-specified: as in script 35, the mean nested LOCO Spearman
#       among the four training countries (each held out in turn, the other
#       three trained, script 04's admin-1 pipeline, same soil block and
#       outcome), floored at 0.
#     sd_train, pre-specified: as in script 35, the mean over the training
#       countries of the between-unit SD of logit prevalence, at the rung of
#       the units being predicted (admin-1, script 04's panel units).
#       SECONDARY (declared before any result): sd_train at the district rung,
#       exactly as the Cote d'Ivoire figure computes it (script 29, Malawi at
#       district), same rho.
#     anchor, pre-specified order, recorded per cell:
#       1. the deposit's own national estimate for the same survey (VMNIS
#          export: same country, begin year, indicator, population group,
#          "national", both urban and rural, no education/wealth breakdown,
#          Gender All/Female, same inflammation adjustment as the admin-1 rows;
#          the row with the largest sample size, then the widest age range);
#       2. else the population-weighted mean of the admin-1 values: NOT
#          AVAILABLE, no population data on disk for these six countries;
#       3. else the unweighted mean of the admin-1 values scored in the cell.
#       Options 2 and 3 are mildly circular: the anchor is built from the very
#       values being scored, which favours the flat A0 most.
#     Scored as MAE (points) over the admin-1 units of each cell.
#
# Pre-registered reading: the ranking adds to level accuracy only where A1's
# MAE is below A0's, on average and in a majority of cells.
#
#   Rscript -e "source('scripts/protocol_v2/78_flat_national_comparator.R')"
# -> results/tables/protocol_v2/lv02_internal_cells.csv      per cell x rank_from x set x f
#    results/tables/protocol_v2/lv02_internal_summary.csv    A1 vs A0 by f, set, rank_from
#    results/tables/protocol_v2/lv02_reproduction.csv        A1 here vs anchor_and_rank_summary.csv
#    results/tables/protocol_v2/lv02_external_cells.csv      29 VMNIS cells
#    results/tables/protocol_v2/lv02_external_summary.csv    A1 vs A0 by arm group
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")   # the tiers script 35 ran at (rerun_downstream.sh)
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"

# ============================================================================ (a) INTERNAL
# ---- verbatim from 35 (lines 47-115) ----------------------------------------
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
  t <- join_admin2_v2(t, admin2_population_v2(POP, cn, pop_for(on)), what = paste("pop", cn, on), quiet = TRUE)
  t <- t[is.finite(t$pop) & t$pop > 0, ]; if (!nrow(t)) return(NULL)
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
# ---- end verbatim ------------------------------------------------------------

rows <- list()
for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on), error = function(e) { cat("  build error", cn, on, conditionMessage(e), "\n"); NULL }); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  for (h in names(cl)) { pool <- setdiff(names(cl), h); Z <- cl[[h]]
    p <- Z$y; w <- Z$w; n <- Z$n; reg <- as.character(Z$region); pv <- clamp(p)
    N <- sum(w); p_nat <- sum(p * w) / N
    p_r <- tapply(p * w, reg, sum) / tapply(w, reg, sum); N_r <- tapply(w, reg, sum); ri <- match(reg, names(p_r))
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
        # the random stream exactly as script 35 consumes it: a_nat, a_r, then B's district noise
        a_nat <- clamp(p_nat + stats::rnorm(1, 0, sqrt(ex * p_nat * (1 - p_nat) / N)))
        a_r <- clamp(p_r + stats::rnorm(length(p_r), 0, sqrt(ex * p_r * (1 - p_r) / N_r)))
        invisible(stats::rnorm(n, 0, sqrt(ex * pv * (1 - pv) / w)))   # B_district_survey's draw, discarded
        est <- list(
          A1_anchor_rank = .v2_expit(.v2_logit(a_nat) + rho_tr * sd_tr * z),
          A0_flat_anchor = rep(a_nat, n))                               # LV-02: no ranking, same anchor draw
        rows[[length(rows) + 1L]] <- data.frame(outcome = on, country = h, rank_from = target, set = set, fraction = f, rep = r,
          n_units = n, rho_train = round(rho_tr, 3), rho_obs = round(rho_obs, 3), rs_used = rho_tr * sd_tr,
          mae_A1 = 100 * mean(abs(p - est$A1_anchor_rank)), mae_A0 = 100 * mean(abs(p - est$A0_flat_anchor)), stringsAsFactors = FALSE)
      }
    }
    cat("done", on, h, "\n")
  }
}
R <- bind_rows(rows)

# ---- reproduction check against script 35's published outputs ----------------
AR <- read.csv(file.path(OUTDIR, "anchor_and_rank.csv"), stringsAsFactors = FALSE) |> filter(design == "A1_anchor_rank")
chk_rows <- R |> inner_join(AR |> select(outcome, country, rank_from, set, fraction, rep, mae_pub = mae),
                            by = c("outcome", "country", "rank_from", "set", "fraction", "rep"))
cat(sprintf("\nreproduction, draw level: %d of %d published A1 draws matched; max |diff| %.2e MAE\n",
            nrow(chk_rows), nrow(AR), max(abs(chk_rows$mae_A1 - chk_rows$mae_pub))))
CELL <- R |> group_by(outcome, country, rank_from, set, fraction, n_units, rho_train, rho_obs, rs_used) |>
  summarise(mae_A1 = mean(mae_A1), mae_A0 = mean(mae_A0), .groups = "drop") |>
  mutate(diff_A1_minus_A0 = mae_A1 - mae_A0, A1_better = mae_A1 < mae_A0)
SUMA1 <- CELL |> group_by(rank_from, set, fraction) |> summarise(cells = dplyr::n(), mae_A1 = round(median(mae_A1), 2), .groups = "drop")
PUB <- read.csv(file.path(OUTDIR, "anchor_and_rank_summary.csv"), stringsAsFactors = FALSE) |> filter(design == "A1_anchor_rank")
REP <- PUB |> select(rank_from, set, fraction, cells_pub = cells, mae_pub = mae) |>
  full_join(SUMA1, by = c("rank_from", "set", "fraction")) |> mutate(abs_diff = abs(mae_A1 - mae_pub))
write.csv(REP, file.path(OUTDIR, "lv02_reproduction.csv"), row.names = FALSE)
cat(sprintf("reproduction, summary level: %d rows, max |A1 here - published| = %.3f MAE, cells here %s vs published %s\n",
            nrow(REP), max(REP$abs_diff, na.rm = TRUE), paste(unique(REP$cells), collapse = "/"), paste(unique(REP$cells_pub), collapse = "/")))
if (any(!is.finite(REP$abs_diff)) || max(REP$abs_diff) > 0.05) stop("LV-02: A1 does not reproduce anchor_and_rank_summary.csv to 0.05 MAE; stopping")

write.csv(CELL, file.path(OUTDIR, "lv02_internal_cells.csv"), row.names = FALSE)
SUMM <- CELL |> group_by(rank_from, set, fraction) |>
  summarise(cells = dplyr::n(), median_mae_A1 = median(mae_A1), median_mae_A0 = median(mae_A0),
            median_diff = median(diff_A1_minus_A0), mean_diff = mean(diff_A1_minus_A0),
            cells_A1_better = sum(A1_better), .groups = "drop") |>
  mutate(ranking_adds = mean_diff < 0 & cells_A1_better > cells / 2)
write.csv(SUMM, file.path(OUTDIR, "lv02_internal_summary.csv"), row.names = FALSE)
cat("\n===== LV-02 (a): A1 (anchor + ranking) minus A0 (flat anchor), district MAE in points, 22 cells =====\n")
print(as.data.frame(SUMM |> mutate(across(c(median_mae_A1, median_mae_A0, median_diff, mean_diff), ~ round(.x, 2)))), row.names = FALSE)

# ============================================================================ (b) EXTERNAL
# ---- verbatim from 04 (setup, build_panel, build_new) ------------------------
XV <- "data/external_validation"
PANEL <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
KEY <- c("country", "Admin1", "Admin2")
VM <- read.csv(file.path(XV, "vmnis_admin1_targets.csv"), stringsAsFactors = FALSE)
XW <- read.csv(file.path(XV, "gadm_to_vmnis_crosswalk.csv"), stringsAsFactors = FALSE)
CL <- read.csv(file.path(XV, "gee_clim_admin2.csv"), check.names = FALSE)
clim_cols <- setdiff(names(CL), KEY)
SOILS <- list(isda = "gee_isda_admin2.csv", sgrid = "gee_sgrid_admin2.csv")
# the pre-specified XV prevalence cells: Africa on iSDA (pre-registered block), South Asia on SoilGrids
ARM_OF_SOIL <- c(isda = "africa", sgrid = "offcontinent")
XVT <- read.csv("results/tables/external_validation/xv_transport.csv", stringsAsFactors = FALSE) |>
  filter(arm == "domain_index", target == "prev", arm_group == ARM_OF_SOIL[soil], is.finite(spearman))
stopifnot(sum(XVT$arm_group == "africa") == 21, sum(XVT$arm_group == "offcontinent") == 8)

# ---- VMNIS national rows (anchor option 1) -----------------------------------
POPGROUP <- c(`Preschool-age children` = "child", `Non-pregnant women (NPW)` = "women", `Women of reproductive age` = "women")
vm_zip <- "data/RA_2026-09/VMNIS.zip"; tmpd <- file.path(tempdir(), "lv02_vmnis"); dir.create(tmpd, showWarnings = FALSE)
read_ind <- function(ind) {
  f <- grep(paste0("VMNISIndicator_", ind, "_"), unzip(vm_zip, list = TRUE)$Name, fixed = TRUE, value = TRUE)
  f <- f[grepl("\\.xlsx$", f)][1]; unzip(vm_zip, files = f, exdir = tmpd, overwrite = TRUE)
  d <- as.data.frame(readxl::read_excel(file.path(tmpd, f), sheet = "Export", guess_max = 100000))
  nm <- names(d); lo <- tolower(nm)
  pc <- nm[grepl("prevalence", lo) & !grepl("cut-off", lo) & !grepl("marginal|severe|elevated|excess", lo)]
  # script 01's rule: first non-missing prevalence column, left to right
  d$prev_any <- if (length(pc)) apply(d[, pc, drop = FALSE], 1, function(r) { r <- suppressWarnings(as.numeric(r)); r <- r[is.finite(r)]; if (length(r)) r[1] else NA_real_ }) else NA_real_
  for (cc in c("Mothers education", "Wealth quantile", "Gender", "Data adjusted for")) if (!cc %in% nm) d[[cc]] <- NA
  d$pg <- unname(POPGROUP[d$Population])
  d$age_months <- (suppressWarnings(as.numeric(d$`Age to`)) - suppressWarnings(as.numeric(d$`Age from`))) * ifelse(d$`Age unit` %in% "Year", 12, 1)
  d
}
IND <- list(); for (ind in unique(VM$indicator)) IND[[ind]] <- read_ind(ind)
national_anchor <- function(cn, yr, on, ind) {
  d <- IND[[ind]]; pg <- sub("_.*$", "", on)
  x <- d[d$Country == cn & d$`Begin year` %in% yr & d$pg %in% pg & is.finite(d$prev_any), ]
  sub_adj <- unique(x$`Data adjusted for`[x$Representativeness. != "national"])
  nat <- x[x$Representativeness. == "national" & x$`Area covered` %in% "both urban and rural" &
             is.na(x$`Mothers education`) & is.na(x$`Wealth quantile`) & (is.na(x$Gender) | x$Gender %in% c("All", "Female")), ]
  if (nrow(nat) && length(sub_adj) == 1 && any(nat$`Data adjusted for` %in% sub_adj)) nat <- nat[nat$`Data adjusted for` %in% sub_adj, ]
  if (!nrow(nat)) return(NULL)
  ss <- suppressWarnings(as.numeric(nat$`Sample size`)); ss[!is.finite(ss)] <- -Inf
  am <- nat$age_months; am[!is.finite(am)] <- -Inf
  i <- order(-ss, -am)[1]
  list(p = nat$prev_any[i] / 100, detail = sprintf("VMNIS national row: %s, %s-%s %s, n %s, adjusted for %s",
       nat$Population[i], nat$`Age from`[i], nat$`Age to`[i], nat$`Age unit`[i], nat$`Sample size`[i], nat$`Data adjusted for`[i]))
}

# ---- district-rung sd_train (secondary; script 29's computation) -------------
lg <- function(p) { p <- pmin(pmax(p, 0.005), 0.995); log(p / (1 - p)) }
sd_district <- function(on) {
  TG |> filter(outcome == on, is.finite(y_prev), is.finite(n_eff), n_eff > 0) |>
    mutate(unit = ifelse(country == "Malawi", Admin1, paste(Admin1, Admin2))) |> group_by(country, unit) |>
    summarise(p = weighted.mean(y_prev, n_eff), .groups = "drop") |> group_by(country) |>
    summarise(sd = sd(lg(p)), .groups = "drop") |> pull(sd) |> mean(na.rm = TRUE)
}

# one fit, script 04's steps: prep per country, common columns, Y z-scored within country, domain PCs signed on training rows
xv_fit <- function(cl, te_name, domain_of_x) {
  Xs <- lapply(cl, function(q) prep_predictors_v2(q$X))
  common <- Reduce(intersect, lapply(Xs, colnames)); if (length(common) < 20) return(NULL)
  Xm <- do.call(rbind, lapply(Xs, function(q) q[, common, drop = FALSE]))
  Y <- unlist(lapply(cl, function(q) as.numeric(scale(q$y))))
  ctry <- rep(names(cl), vapply(cl, function(q) length(q$y), 0L))
  te <- which(ctry == te_name); tr <- which(ctry != te_name)
  Dm <- domain_representation_v2(Xm, domain_of_x, sign_rows = tr)
  aux <- list(lon = rep(0, length(Y)), lat = rep(0, length(Y)), y_nat = Y)
  p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, Xm, Dm, aux), error = function(e) rep(NA_real_, length(te)))
  if (length(p) != length(te)) p <- rep(NA_real_, length(te))
  list(p = p, y = Y[te])
}

ext <- list()
for (soil in names(SOILS)) {
  IS <- read.csv(file.path(XV, SOILS[[soil]]), check.names = FALSE)
  soil_cols <- setdiff(names(IS), KEY); PREDS_X <- c(clim_cols, soil_cols)
  COV <- dplyr::inner_join(CL, IS, by = KEY) |> group_by(across(all_of(KEY))) |>
    summarise(across(all_of(PREDS_X), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |> as.data.frame()
  XWu <- unique(XW[, c("country", "Admin1", "Admin2", "vmnis_unit")])
  domain_of_x <- stats::setNames(c(rep("Climate and weather", length(clim_cols)), rep("Soil characteristics", length(soil_cols))), PREDS_X)
  build_panel <- function(cn, on, target) {           # verbatim from 04
    t <- TG[TG$country == cn & TG$outcome == on, ]
    ycol <- if (target == "prev") "y_prev" else "y_level"; wcol <- if (target == "prev") "n_eff" else "n_eff_cont"
    t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]) & t[[wcol]] > 0, ]; if (!nrow(t)) return(NULL)
    a1 <- t |> group_by(unit = Admin1) |> summarise(y = stats::weighted.mean(.data[[ycol]], .data[[wcol]]), w = sum(.data[[wcol]]), .groups = "drop")
    x1 <- COV[COV$country == cn, ] |> group_by(unit = Admin1) |> summarise(across(all_of(PREDS_X), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    m <- dplyr::inner_join(a1, x1, by = "unit"); if (nrow(m) < 3) return(NULL)
    list(country = cn, unit = m$unit, y = m$y, w = m$w, X = as.matrix(m[, PREDS_X, drop = FALSE]))
  }
  build_new <- function(cn, on, target) {             # verbatim from 04
    v <- VM[VM$country == cn & VM$outcome == on, ]; if (!nrow(v)) return(NULL)
    v$y <- if (target == "prev") v$prev / 100 else -v$level
    v <- v[is.finite(v$y), c("unit", "y")]; if (nrow(v) < 5) return(NULL)
    cv <- COV[COV$country == cn, ] |> dplyr::inner_join(XWu, by = KEY, relationship = "many-to-one")
    x1 <- cv |> group_by(unit = vmnis_unit) |> summarise(across(all_of(PREDS_X), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    m <- dplyr::inner_join(v, x1, by = "unit"); if (nrow(m) < 5) return(NULL)
    list(country = cn, unit = m$unit, y = m$y, w = rep(1, nrow(m)), X = as.matrix(m[, PREDS_X, drop = FALSE]))
  }
  cells_here <- XVT[XVT$soil == soil, ]
  for (on in sort(unique(cells_here$outcome))) {
    tr_list <- list(); for (cn in PANEL) { z <- tryCatch(build_panel(cn, on, "prev"), error = function(e) NULL); if (!is.null(z)) tr_list[[cn]] <- z }
    # rho_train: nested LOCO among the training countries, script 04's pipeline; sd_train: admin-1 rung
    rt <- vapply(names(tr_list), function(t) { f <- xv_fit(c(tr_list[setdiff(names(tr_list), t)], tr_list[t]), t, domain_of_x)
      if (is.null(f) || !all(is.finite(f$p)) || stats::sd(f$p) == 0) NA_real_ else suppressWarnings(stats::cor(tr_list[[t]]$y, f$p, method = "spearman")) }, numeric(1))
    rho_tr <- mean(rt, na.rm = TRUE); if (!is.finite(rho_tr) || rho_tr < 0) rho_tr <- 0
    sd_a1 <- mean(vapply(tr_list, function(q) stats::sd(.v2_logit(clamp(q$y))), numeric(1)), na.rm = TRUE)
    sd_d <- sd_district(on)
    for (cn in cells_here$country[cells_here$outcome == on]) {
      z <- build_new(cn, on, "prev"); f <- xv_fit(c(tr_list, stats::setNames(list(z), cn)), cn, domain_of_x)
      pub <- cells_here$spearman[cells_here$country == cn & cells_here$outcome == on]
      rho_here <- suppressWarnings(stats::cor(f$y, f$p, method = "spearman"))
      if (!is.finite(rho_here) || abs(rho_here - pub) > 1e-6) stop(sprintf("LV-02 (b): %s %s %s Spearman %.6f does not reproduce xv_transport.csv %.6f; stopping", soil, cn, on, rho_here, pub))
      zz <- as.numeric(scale(f$p)); y <- z$y
      yr <- unique(VM$survey_year[VM$country == cn]); ind <- unique(VM$indicator[VM$country == cn & VM$outcome == on])
      na <- national_anchor(cn, yr, on, ind)
      if (!is.null(na)) { anc <- na$p; src <- "1_deposit_national"; det <- na$detail
      } else { anc <- mean(y); src <- "3_unweighted_admin1_mean"; det <- "no national row in the deposit; no population for option 2" }
      a <- clamp(anc)
      A1 <- .v2_expit(.v2_logit(a) + rho_tr * sd_a1 * zz); A1d <- .v2_expit(.v2_logit(a) + rho_tr * sd_d * zz); A0 <- rep(a, length(y))
      ext[[length(ext) + 1L]] <- data.frame(soil = soil, arm_group = ARM_OF_SOIL[[soil]], country = cn, outcome = on, n_units = length(y),
        spearman = rho_here, spearman_published = pub, rho_train = rho_tr, rho_train_by_country = paste(sprintf("%s %.2f", names(rt), rt), collapse = "; "),
        sd_train_admin1 = sd_a1, sd_train_district = sd_d, anchor = anc, anchor_source = src, anchor_detail = det,
        admin1_unweighted_mean = mean(y), mae_A0 = 100 * mean(abs(y - A0)), mae_A1 = 100 * mean(abs(y - A1)),
        mae_A1_sd_district = 100 * mean(abs(y - A1d)), stringsAsFactors = FALSE)
    }
  }
  cat("  external", soil, "done\n")
}
E <- bind_rows(ext) |> mutate(diff_A1_minus_A0 = mae_A1 - mae_A0, A1_better = mae_A1 < mae_A0,
                              diff_sd_district = mae_A1_sd_district - mae_A0, A1_sd_district_better = mae_A1_sd_district < mae_A0)
stopifnot(nrow(E) == 29)
write.csv(E, file.path(OUTDIR, "lv02_external_cells.csv"), row.names = FALSE)
esum <- function(d, label) data.frame(set = label, cells = nrow(d), mean_mae_A0 = mean(d$mae_A0), mean_mae_A1 = mean(d$mae_A1),
  mean_diff = mean(d$diff_A1_minus_A0), median_diff = median(d$diff_A1_minus_A0), cells_A1_better = sum(d$A1_better),
  ranking_adds = mean(d$diff_A1_minus_A0) < 0 && sum(d$A1_better) > nrow(d) / 2,
  mean_diff_sd_district = mean(d$diff_sd_district), cells_A1_sd_district_better = sum(d$A1_sd_district_better), stringsAsFactors = FALSE)
ES <- bind_rows(esum(E[E$arm_group == "africa", ], "Africa (iSDA), 21"), esum(E[E$arm_group == "offcontinent", ], "South Asia (SoilGrids), 8"),
                esum(E, "all 29"), esum(E[E$anchor_source == "1_deposit_national", ], "anchor = deposit national"),
                esum(E[E$anchor_source != "1_deposit_national", ], "anchor = unweighted admin-1 mean"))
write.csv(ES, file.path(OUTDIR, "lv02_external_summary.csv"), row.names = FALSE)
cat("\n===== LV-02 (b): WHO VMNIS admin-1 prevalence, A1 minus A0 (MAE, points) =====\n")
print(as.data.frame(E |> transmute(arm_group, country, outcome, n_units, spearman = round(spearman, 2), rho_train = round(rho_train, 2),
  sd_a1 = round(sd_train_admin1, 2), anchor = round(100 * anchor, 1), anchor_source, mae_A0 = round(mae_A0, 2), mae_A1 = round(mae_A1, 2),
  diff = round(diff_A1_minus_A0, 2))), row.names = FALSE)
print(ES |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
cat("\nDONE\n")
