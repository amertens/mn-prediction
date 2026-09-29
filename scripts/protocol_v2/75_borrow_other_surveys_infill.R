# =============================================================================
# scripts/protocol_v2/75_borrow_other_surveys_infill.R   [BO-01, 2026-09-28]
#
# BORROW THE OTHER SURVEYS WHEN FILLING IN A SURVEYED COUNTRY
#
# The in-fill index learns its weights from the 24-70 training districts of one
# survey. The same outcome was measured in up to three other surveys (about
# 130-180 districts). Does adding their evidence to the weights help inside a
# surveyed country?
#
# ── PRE-REGISTRATION (written before any result was seen; not changed after) ──
#
# ARMS (nothing tuned)
#   domain_index         the deployed arm, exactly as 02_run_benchmarks_v2.R
#                        runs it: per-country domain components over all of the
#                        country's columns, ARMS_V2$domain_index, rho = 1.
#   index_borrow         z_total = z_own + sum_C z_C  (precision-weighted
#                        fixed-effect sum; the index's own construction)
#       z_own = .index_weights_v2(D_T[train], y_T[train])
#       z_C   = .index_weights_v2(D_C, y_C) on ALL of donor C's districts
#       donors = every other country with outcome O whose cell builds under
#                02b's build_cell (>= 12 districts, >= 3 regions, >= 20 columns)
#       D = the components 02b uses to make domain PCs comparable across
#           countries: each country's predictors rank-normalised WITHIN country
#           (prep_predictors_v2), restricted to the columns common to T and all
#           donors, stacked, and ONE principal-component basis per domain learned
#           with sign_rows = the fold's training rows (T's training districts +
#           every donor district), then applied to every row. This mirrors
#           02b line 107 exactly (sign_rows = the training rows). A per-country
#           basis is NOT used for borrowing: its components are different linear
#           combinations in each country (the majority-positive flip fixes a
#           sign, not an axis), so summing weights over countries would add
#           incomparable axes. bo01_borrow_alignment.csv quantifies this.
#       prediction = D_T %*% z_total, rescaled exactly as arm_domain_index_v2
#           (training-index mean and sd, then the training outcome's sd and
#           mean, rho = 1).
#       No donor (Malawi child/women zinc), or fewer than 20 common columns:
#           nothing to borrow, index_borrow = domain_index; a tie, not a win.
#   domain_index_common  z_own alone on the SAME fold-specific pooled basis.
#                        Diagnostic: domain_index_common - domain_index is the
#                        change of representation (common columns, pooled
#                        basis); index_borrow - domain_index_common is the
#                        borrowing itself. Not part of the verdict.
#
# TARGET: level = y_level (native scale); prevalence = .v2_logit(y_prev),
#   back-transformed with expit before scoring, exactly as 02.
#
# LEAKAGE: T's held-out districts never enter a weight. z_own uses T's training
#   rows; the basis sees T's training rows only (and no outcome); the donors are
#   separate surveys, used in full. Checked at run time: on the first draw of
#   every cell (and every region fold) T's held-out outcomes are replaced by
#   arbitrary values and the index_borrow predictions must not move (1e-10).
#
# PRIMARY TEST
#   index_borrow against domain_index, paired on identical folds; level target;
#   in-fill, make_folds_v2("kfold_district", n, k = 5, rep_id = 1..10); per-cell
#   Spearman = mean over the 10 draws; the 18 cells with a finite in-fill level
#   domain_index Spearman in benchmarks_v2_cells.csv.
#   PASS iff the mean gain >= +0.03 AND index_borrow is strictly better in
#   >= 12 of the 18. Otherwise FAIL.
#   REPRODUCTION FIRST: domain_index must match benchmarks_v2_cells.csv to 0.005
#   in every primary cell, or the script stops before writing a verdict.
#
# SECONDARY (reported, not used for the verdict)
#   prevalence target (logit) in-fill; region hold-out (exhaustive leave-one-
#   Admin1-out), level and prevalence; pairs (share of district pairs in the
#   survey's order; in-fill uses each district's held-out prediction averaged
#   over the ten draws, as scripts/policy_deck/28); gain by country.
#   Prior expectation: the largest gain in The Gambia and under region hold-out.
#
# Settings forced to the protocol defaults: V2_PREDICTOR_TIERS=open,survey_public
# (the headline tiers), V2_PREP_SCALE=rank, V2_DOMAIN_REP=pcvar,
# V2_INDEX_SHRINK=none. One R process, no parallel workers.
#
#   Rscript -e "source('scripts/protocol_v2/75_borrow_other_surveys_infill.R')"
#   PROFILE=smoke  Gambia only, 2 draws, outputs to BO01_OUTDIR (default tempdir)
# -> results/tables/protocol_v2/bo01_borrow_raw.csv        one row per draw x arm
#    results/tables/protocol_v2/bo01_borrow_cells.csv      per cell x target x estimand
#    results/tables/protocol_v2/bo01_borrow_summary.csv    verdict + secondary summaries
#    results/tables/protocol_v2/bo01_borrow_alignment.csv  per-country PC alignment check
#    results/tables/protocol_v2/bo01_borrow_reproduction.csv  domain_index vs the benchmark
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public", V2_PREP_SCALE = "rank",
           V2_DOMAIN_REP = "pcvar", V2_INDEX_SHRINK = "none")
source("R/protocol_v2.R")
t_start <- Sys.time()
stamp <- function(...) cat(sprintf("[%s | %5.1f min] ", format(Sys.time(), "%H:%M:%S"),
                                   as.numeric(difftime(Sys.time(), t_start, units = "mins"))), ..., "\n", sep = "")

PROFILE <- Sys.getenv("PROFILE", "full")
P2      <- "results/tables/protocol_v2"
OUTDIR  <- if (PROFILE == "smoke") Sys.getenv("BO01_OUTDIR", tempdir()) else P2
REPS    <- if (PROFILE == "smoke") 2L else 10L
TOL_REP <- 0.005

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND  <- readRDS("dashboard/data/admin2_boundaries.rds")
BENCH <- read.csv(file.path(P2, "benchmarks_v2_cells.csv"), stringsAsFactors = FALSE)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
TARGET_COUNTRIES <- if (PROFILE == "smoke") "Gambia" else unname(COUNTRIES)

CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]],
             Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2),
             lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE)
}))

# ── one cell: 02's build_cell (identical rows and X to 02b's build_cell) ─────
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
  D  <- domain_representation_v2(Xr, domain_of)
  y_nat <- m[[ycol]]
  y_mod <- if (target == "prev") .v2_logit(y_nat) else y_nat
  list(country = cn, outcome = on, target = target, y_nat = y_nat, y_mod = y_mod,
       X = Xr, D = D, w = m[[wcol]], Admin1 = m$Admin1, Admin2 = m$Admin2, n = nrow(m),
       aux = list(lon = m$lon, lat = m$lat, Admin1 = m$Admin1, y_nat = y_nat, target = target))
}

# ── the pooled design for (T, O, target): 02b's common columns, stacked ──────
make_pool <- function(cell, donors) {
  if (!length(donors)) return(NULL)
  common <- Reduce(intersect, c(list(colnames(cell$X)), lapply(donors, function(d) colnames(d$X))))
  if (length(common) < 20) return(NULL)                      # 02b's threshold
  Xm <- do.call(rbind, c(list(cell$X[, common, drop = FALSE]),
                         lapply(donors, function(d) d$X[, common, drop = FALSE])))
  off <- cell$n; drows <- list()
  for (d in donors) { drows[[d$country]] <- off + seq_len(d$n); off <- off + d$n }
  list(Xm = Xm, common = common, rows_T = seq_len(cell$n), donor_rows = drows,
       donor_y = lapply(donors, function(d) d$y_mod), donors = names(donors),
       donor_n = vapply(donors, function(d) d$n, 0L))
}

# arm_domain_index_v2's rescale with rho = 1, for a given weight vector
idx_predict <- function(D, tr, te, z, y) {
  ytr <- y[tr]
  itr <- as.numeric(D[tr, , drop = FALSE] %*% z)
  ite <- as.numeric(D[te, , drop = FALSE] %*% z)
  if (stats::sd(itr) == 0) return(rep(mean(ytr), length(te)))
  ((ite - mean(itr)) / stats::sd(itr)) * 1 * stats::sd(ytr) + mean(ytr)
}

# one fold of index_borrow and domain_index_common on the pooled basis
borrow_fold <- function(pool, tr, te, y) {
  sr <- c(pool$rows_T[tr], unlist(pool$donor_rows, use.names = FALSE))   # training rows only (02b)
  Dm <- domain_representation_v2(pool$Xm, domain_of, sign_rows = sr)
  DT <- Dm[pool$rows_T, , drop = FALSE]
  z_own <- .index_weights_v2(DT[tr, , drop = FALSE], y[tr])
  z_C   <- Map(function(r, yy) .index_weights_v2(Dm[r, , drop = FALSE], yy), pool$donor_rows, pool$donor_y)
  z_don <- Reduce(`+`, z_C)
  list(borrow = idx_predict(DT, tr, te, z_own + z_don, y),
       common = idx_predict(DT, tr, te, z_own, y),
       ncomp = ncol(Dm),
       r_own_don = suppressWarnings(stats::cor(z_own, z_don)),
       norm_ratio = sqrt(sum(z_don^2)) / max(sqrt(sum(z_own^2)), 1e-12))
}

concord <- function(o, p) {          # scripts/policy_deck/28, ties in either skipped
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (length(o) < 5) return(NA_real_)
  so <- sign(outer(o, o, "-")); sp <- sign(outer(p, p, "-"))
  m <- so * sp; ut <- upper.tri(m)
  agree <- sum(m[ut] > 0); tot <- sum(m[ut] != 0)
  if (tot == 0) NA_real_ else agree / tot
}

ARMS <- c("domain_index", "domain_index_common", "index_borrow")

# run all three arms over one fold assignment; returns natural-scale predictions
run_folds <- function(cell, pool, folds, leak_check) {
  n <- cell$n; y <- cell$y_mod
  P <- matrix(NA_real_, n, 3, dimnames = list(NULL, ARMS))
  diag <- list()
  for (f in unique(folds)) {
    te <- which(folds == f); tr <- which(folds != f)
    if (length(tr) < 12 || !length(te)) next                              # 02's rule
    p_di <- ARMS_V2$domain_index(tr, te, y, cell$X, cell$D, cell$aux)
    # the rescale here must be the arm's own, bit for bit
    p_chk <- idx_predict(cell$D, tr, te, .index_weights_v2(cell$D[tr, , drop = FALSE], y[tr]), y)
    if (max(abs(p_chk - p_di)) > 1e-10) stop("idx_predict does not mirror arm_domain_index_v2")
    P[te, "domain_index"] <- p_di
    if (is.null(pool)) {
      P[te, "domain_index_common"] <- p_di; P[te, "index_borrow"] <- p_di
    } else {
      b <- borrow_fold(pool, tr, te, y)
      P[te, "domain_index_common"] <- b$common; P[te, "index_borrow"] <- b$borrow
      diag[[length(diag) + 1]] <- c(ncomp = b$ncomp, r_own_don = b$r_own_don, norm_ratio = b$norm_ratio)
      if (leak_check) {
        y_bad <- y; y_bad[te] <- y[te] + 1e3 * seq_along(te) - 7
        b2 <- borrow_fold(pool, tr, te, y_bad)
        if (max(abs(b2$borrow - b$borrow)) > 1e-10 || max(abs(b2$common - b$common)) > 1e-10)
          stop(sprintf("LEAK: %s %s %s held-out outcomes move the prediction", cell$country, cell$outcome, cell$target))
        LEAK_CHECKS <<- LEAK_CHECKS + 1L
      }
    }
  }
  if (cell$target == "prev") P <- .v2_expit(P)
  list(P = P, diag = if (length(diag)) do.call(rbind, diag) else NULL)
}

score_row <- function(cell, pred, estimand, arm, rep) {
  s <- score_v2(cell$y_nat, pred, cell$w, scale = if (cell$target == "prev") "prev" else "level")
  data.frame(country = cell$country, outcome = cell$outcome, target = cell$target,
             estimand = estimand, arm = arm, rep = rep, n_areas = cell$n,
             spearman = s$spearman, pearson = s$pearson, topk = s$topk)
}

bench_ref <- function(cn, on, tg, es) {
  r <- BENCH$spearman[BENCH$country == cn & BENCH$outcome == on & BENCH$target == tg &
                        BENCH$estimand == es & BENCH$arm == "domain_index"]
  if (length(r) != 1) NA_real_ else r
}

# ── build every cell once ────────────────────────────────────────────────────
OUTCOMES <- unique(TG$outcome)
B <- list()
for (tg in c("level", "prev")) for (cn in COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) {
  cc <- tryCatch(build_cell(cn, on, tg), error = function(e) { message("build failed ", cn, " ", on, ": ", conditionMessage(e)); NULL })
  if (!is.null(cc)) B[[paste(cn, on, tg)]] <- cc
}
stamp("built ", length(B), " cells")

# ── per-country bases are not aligned: the check that justifies the pooled basis ──
al <- list()
for (on in OUTCOMES) {
  cl <- Filter(Negate(is.null), lapply(COUNTRIES, function(cn) B[[paste(cn, on, "level")]]))
  if (length(cl) < 2) next
  names(cl) <- vapply(cl, function(z) z$country, "")
  prs <- utils::combn(names(cl), 2)
  for (j in seq_len(ncol(prs))) {
    Ba <- attr(cl[[prs[1, j]]]$D, "basis"); Bb <- attr(cl[[prs[2, j]]]$D, "basis")
    for (dm in intersect(names(Ba), names(Bb))) {
      sc <- intersect(Ba[[dm]]$cols, Bb[[dm]]$cols)
      if (length(sc) < 2) next
      for (k in intersect(colnames(Ba[[dm]]$rot), colnames(Bb[[dm]]$rot))) {
        va <- Ba[[dm]]$rot[sc, k]; vb <- Bb[[dm]]$rot[sc, k]
        den <- sqrt(sum(va^2) * sum(vb^2)); if (den == 0) next
        al[[length(al) + 1]] <- data.frame(outcome = on, country_a = prs[1, j], country_b = prs[2, j],
          domain = dm, component = k, pc_rank = as.integer(sub(".*_PC", "", k)), n_shared_cols = length(sc),
          cosine = sum(va * vb) / den)
      }
    }
  }
}
AL <- bind_rows(al)
AL_SUM <- AL |> mutate(which_pc = ifelse(pc_rank == 1, "PC1", "PC2+")) |>
  group_by(outcome, which_pc) |>
  summarise(n_compared = n(), median_cosine = median(cosine), share_negative = mean(cosine < 0),
            share_abs_ge_0.7 = mean(abs(cosine) >= 0.7), .groups = "drop")
write.csv(AL_SUM, file.path(OUTDIR, "bo01_borrow_alignment.csv"), row.names = FALSE)
cat("\nPer-country bases, same-named components compared across countries (loadings on shared columns):\n")
print(as.data.frame(AL_SUM |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# ── main loop ────────────────────────────────────────────────────────────────
LEAK_CHECKS <- 0L
raw <- list(); cellrows <- list(); repro <- list()
same_or_close <- function(ours, ref) (!is.finite(ref) && !is.finite(ours)) ||
  (is.finite(ref) && is.finite(ours) && abs(ours - ref) <= TOL_REP)

for (tg in c("level", "prev")) for (cn in TARGET_COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) {
  cell <- B[[paste(cn, on, tg)]]; if (is.null(cell)) next
  donors <- Filter(Negate(is.null), lapply(setdiff(COUNTRIES, cn), function(c2) B[[paste(c2, on, tg)]]))
  names(donors) <- vapply(donors, function(d) d$country, "")
  stopifnot(!cn %in% names(donors))
  pool <- make_pool(cell, donors)
  n <- cell$n

  # A. in-fill, 5-fold x REPS draws (identical folds for every arm)
  PA <- array(NA_real_, c(n, REPS, 3), dimnames = list(NULL, NULL, ARMS))
  SPA <- matrix(NA_real_, REPS, 3, dimnames = list(NULL, ARMS)); dg <- list()
  for (r in seq_len(REPS)) {
    folds <- make_folds_v2("kfold_district", n, k = 5, rep_id = r)
    out <- run_folds(cell, pool, folds, leak_check = (r == 1))
    PA[, r, ] <- out$P; dg[[r]] <- out$diag
    for (a in ARMS) {
      s <- score_row(cell, out$P[, a], "infill", a, r)
      SPA[r, a] <- s$spearman; raw[[length(raw) + 1]] <- s
    }
  }
  spA <- colMeans(SPA, na.rm = TRUE)                      # the benchmark's own summary
  spA[!is.finite(spA)] <- NA_real_
  ref <- bench_ref(cn, on, tg, "infill"); ok <- same_or_close(spA[["domain_index"]], ref)
  repro[[length(repro) + 1]] <- data.frame(country = cn, outcome = on, target = tg, estimand = "infill",
                                           ours = spA[["domain_index"]], bench = ref, reproduced = ok)
  if (PROFILE != "smoke" && tg == "level" && is.finite(ref) && !ok) {
    write.csv(bind_rows(repro), file.path(OUTDIR, "bo01_borrow_reproduction.csv"), row.names = FALSE)
    stop(sprintf("REPRODUCTION FAILED %s %s level in-fill: %.4f vs benchmark %.4f", cn, on, spA[["domain_index"]], ref))
  }
  # pairs: each district's held-out prediction averaged over the draws (script 28)
  pairs_in <- vapply(ARMS, function(a) concord(cell$y_nat, rowMeans(matrix(PA[, , a], nrow = n))), 0)
  D_ALL <- do.call(rbind, dg)

  # B. region hold-out, exhaustive leave-one-Admin1-out (deterministic)
  folds <- make_folds_v2("loro", n, blocks = cell$Admin1)
  outB <- run_folds(cell, pool, folds, leak_check = TRUE)
  spB <- setNames(rep(NA_real_, 3), ARMS)
  for (a in ARMS) {
    s <- score_row(cell, outB$P[, a], "region", a, 1L)
    spB[[a]] <- s$spearman; raw[[length(raw) + 1]] <- s
  }
  refB <- bench_ref(cn, on, tg, "region"); okB <- same_or_close(spB[["domain_index"]], refB)
  repro[[length(repro) + 1]] <- data.frame(country = cn, outcome = on, target = tg, estimand = "region",
                                           ours = spB[["domain_index"]], bench = refB, reproduced = okB)
  pairs_rg <- vapply(ARMS, function(a) concord(cell$y_nat, outB$P[, a]), 0)

  info <- data.frame(
    donors = if (is.null(pool)) "" else paste(pool$donors, collapse = "+"),
    donor_districts = if (is.null(pool)) 0L else sum(pool$donor_n),
    n_x_own = ncol(cell$X), n_x_common = if (is.null(pool)) ncol(cell$X) else length(pool$common),
    n_comp_own = ncol(cell$D),
    n_comp_pooled = if (is.null(D_ALL)) NA_real_ else mean(D_ALL[, "ncomp"]),
    r_own_donor_weights = if (is.null(D_ALL)) NA_real_ else mean(D_ALL[, "r_own_don"], na.rm = TRUE),
    donor_to_own_weight_norm = if (is.null(D_ALL)) NA_real_ else mean(D_ALL[, "norm_ratio"], na.rm = TRUE),
    stringsAsFactors = FALSE)
  for (es in c("infill", "region")) {
    sp <- if (es == "infill") spA else spB
    pr <- if (es == "infill") pairs_in else pairs_rg
    cellrows[[length(cellrows) + 1]] <- data.frame(
      country = cn, outcome = on, target = tg, estimand = es, n_areas = n,
      borrowed = !is.null(pool), info,
      bench_domain_index = if (es == "infill") ref else refB,
      reproduced = if (es == "infill") ok else okB,
      sp_domain_index = sp[["domain_index"]], sp_domain_index_common = sp[["domain_index_common"]],
      sp_index_borrow = sp[["index_borrow"]],
      gain = sp[["index_borrow"]] - sp[["domain_index"]],
      gain_representation = sp[["domain_index_common"]] - sp[["domain_index"]],
      gain_borrowing_same_basis = sp[["index_borrow"]] - sp[["domain_index_common"]],
      pairs_domain_index = pr[["domain_index"]], pairs_domain_index_common = pr[["domain_index_common"]],
      pairs_index_borrow = pr[["index_borrow"]], stringsAsFactors = FALSE)
  }
  stamp(sprintf("%-11s %-12s %-5s n=%2d donors=%-26s infill %.3f / %.3f / %.3f | region %.3f / %.3f / %.3f  (index / common / borrow)",
                cn, on, tg, n, info$donors, spA[1], spA[2], spA[3], spB[1], spB[2], spB[3]))
}

RAW   <- bind_rows(raw)
CELLS <- bind_rows(cellrows)
REPRO <- bind_rows(repro)
write.csv(RAW,   file.path(OUTDIR, "bo01_borrow_raw.csv"), row.names = FALSE)
write.csv(CELLS, file.path(OUTDIR, "bo01_borrow_cells.csv"), row.names = FALSE)
write.csv(REPRO, file.path(OUTDIR, "bo01_borrow_reproduction.csv"), row.names = FALSE)

# ── verdict and summaries ────────────────────────────────────────────────────
# a cell is scored where the benchmark scored domain_index (the 18 primary cells
# on the level in-fill; Sierra Leone's 14 districts leave < 12 training rows)
summ <- function(d, setting, group = "all") {
  d <- d[is.finite(d$bench_domain_index) & is.finite(d$sp_index_borrow) & is.finite(d$sp_domain_index), ]
  if (!nrow(d)) return(NULL)
  data.frame(setting = setting, group = group, cells = nrow(d),
             mean_domain_index = mean(d$sp_domain_index), mean_domain_index_common = mean(d$sp_domain_index_common),
             mean_index_borrow = mean(d$sp_index_borrow),
             mean_gain = mean(d$gain), median_gain = median(d$gain),
             wins = sum(d$gain > 0), losses = sum(d$gain < 0), ties = sum(d$gain == 0),
             mean_gain_representation = mean(d$gain_representation),
             mean_gain_borrowing_same_basis = mean(d$gain_borrowing_same_basis),
             wins_borrow_vs_common = sum(d$gain_borrowing_same_basis > 0),
             mean_pairs_domain_index = mean(d$pairs_domain_index),
             mean_pairs_domain_index_common = mean(d$pairs_domain_index_common),
             mean_pairs_index_borrow = mean(d$pairs_index_borrow),
             verdict = "", stringsAsFactors = FALSE)
}
S_ALL <- list()
for (tg in c("level", "prev")) for (es in c("infill", "region")) {
  d <- CELLS[CELLS$target == tg & CELLS$estimand == es, ]
  st <- paste(es, tg, sep = "_")
  s0 <- summ(d, st); if (is.null(s0)) next
  if (st == "infill_level") {
    if (PROFILE != "smoke") stopifnot(s0$cells == 18, all(d$reproduced[is.finite(d$bench_domain_index)]))
    s0$verdict <- sprintf("%s (pre-registered: mean gain >= +0.03 AND wins >= 12 of 18; observed %+.3f and %d of %d)",
                          if (s0$mean_gain >= 0.03 && s0$wins >= 12) "PASS" else "FAIL",
                          s0$mean_gain, s0$wins, s0$cells)
  } else s0$verdict <- "secondary"
  S_ALL[[length(S_ALL) + 1]] <- s0
  for (cn in unique(d$country)) S_ALL[[length(S_ALL) + 1]] <- summ(d[d$country == cn, ], st, cn)
}
SUMM <- bind_rows(S_ALL)
write.csv(SUMM, file.path(OUTDIR, "bo01_borrow_summary.csv"), row.names = FALSE)

cat("\n=============== BO-01 PRIMARY: level, in-fill, index_borrow vs domain_index ===============\n")
pc <- CELLS[CELLS$target == "level" & CELLS$estimand == "infill" & is.finite(CELLS$bench_domain_index), ]
print(pc |> transmute(country, outcome, donors, x_own = n_x_own, x_common = n_x_common, comp_own = n_comp_own,
                      comp_pool = round(n_comp_pooled, 1), bench = round(bench_domain_index, 3),
                      index = round(sp_domain_index, 3), common = round(sp_domain_index_common, 3),
                      borrow = round(sp_index_borrow, 3), gain = round(gain, 3),
                      r_w = round(r_own_donor_weights, 2), w_ratio = round(donor_to_own_weight_norm, 2)), row.names = FALSE)
cat("\n"); print(SUMM |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
cat(sprintf("\nReproduction: domain_index matches benchmarks_v2_cells.csv in %d of %d cell-settings (max |diff| %.4f); %d leakage checks passed\n",
            sum(REPRO$reproduced), nrow(REPRO), max(abs(REPRO$ours - REPRO$bench), na.rm = TRUE), LEAK_CHECKS))
if (!all(REPRO$reproduced)) print(REPRO[!REPRO$reproduced, ], row.names = FALSE)
stamp("DONE")
