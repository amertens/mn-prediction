# =============================================================================
# scripts/signal_probes/p5_omnibus_dense_signal.R
#
# P5. Is there ANY covariate signal in a cell? An omnibus test with power
#     against MANY SMALL EFFECTS, plus a spatial-confound control.
#
# WHY THIS TEST DOES NOT EXIST YET
# --------------------------------
# The project's headline evidence for "no covariate signal" is
#   "0 of 294 predictors survive FDR control, in all 24 cells"
# (scripts/covariates/16_bivariate_fdr.R). That is a PER-PREDICTOR test with a
# multiplicity correction over ~300 tests at n = 14-87 areas. Its detectable
# effect size is |r| ~ 0.45-0.6. If the truth is "150 predictors each with
# |r| ~ 0.20, all pointing the same way", that design has essentially zero
# power and returns 0/294 with near-certainty. "0 of 294 survive" is then a
# statement about POWER, not about the absence of signal.
#
# The right instrument for a dense weak alternative is an OMNIBUS statistic
# that pools evidence across predictors instead of penalising for looking at
# them. This script runs two, side by side, on identical data and folds:
#
#   Q    = mean_j r_j^2      dense alternative  (Goeman globaltest / SKAT with a
#                            linear kernel is monotone in this; it is the
#                            locally most powerful statistic when many
#                            predictors carry small independent effects)
#   Tmax = max_j |r_j|       sparse alternative (what the FDR scan effectively
#                            tests, restated as a single statistic)
#
# Both are calibrated by the SAME permutation null, so the CONTRAST between
# them is interpretable:
#   Q significant, Tmax not  -> signal is real, dense and weak; per-predictor
#                               screening cannot see it, and no amount of FDR
#                               correction will make it appear.
#   Tmax significant, Q not  -> a few strong predictors; screening is right.
#   neither                  -> no detectable association in this cell.
#
# THE SPATIAL CONFOUND, AND WHY IT IS CONTROLLED HERE
# ---------------------------------------------------
# A free permutation of the outcome across areas destroys the outcome's spatial
# autocorrelation, so it counts "deficiency is smooth in space AND covariates
# are smooth in space" as signal. That is the standing caveat on P1/P4 and the
# check (p3b) that was written but never ran. Three nested variants are run:
#
#   S0_none    no spatial control (comparable to P1/P4)
#   S1_latlon  outcome and every predictor residualised on centroid lat+lon
#              (2 df) - the p3b design, generalised to Admin-2
#   S2_mem     residualised on the leading Moran eigenvector maps (MEM) of a
#              k-nearest-neighbour graph - a flexible smooth spatial field, so
#              this asks the strict question: does the covariate block carry
#              information BEYOND ANY SMOOTH SPATIAL SURFACE?
#
# S2_mem is the honest inferential analogue of the project's covariate-free
# spatial smoother comparator. If Q survives S2_mem, the covariates carry
# information that geography alone does not.
# Permutation under S1/S2 is Freedman-Lane: the reduced-model residual vector
# is permuted, which is the correct null for a partial association.
#
# INTERPRETING A NEGATIVE RESULT UNDER S2_mem
# -------------------------------------------
# Losing significance under S2_mem does NOT mean the covariates are useless.
# For a country with no survey, no spatial smoother can be fitted (it needs
# observed outcomes), so a covariate signal that is collinear with geography is
# still the only deployable signal there. S2_mem answers a scientific question
# ("is this more than geography?"), not the deployment question (P3a does that).
#
# DESIGN
# ------
# unit         admin2 (Gambia 30, Ghana 75, Malawi 87, SierraLeone 14) and
#              admin1 (6, 16, 27, 4). Admin-2 is joined on the PAIR KEY
#              (Admin1, Admin2) - Malawi has Admin-2 names that repeat across
#              regions, and a bare-name join fans rows (documented defect).
# outcome      svy_prev; at admin1, the n-weighted mean of member districts
# predictors   the 373 shared harmonised predictors, rank-normalised within
#              country over the analysis rows (outcome-independent)
# statistic    r_j = Pearson on rank-transformed, spatially-residualised values
#              (= Spearman when the spatial basis is empty)
# null         permute the (residualised) outcome across areas within country;
#              ONE draw per country per replicate, reused across that country's
#              outcome cells, so the pooled test respects cross-outcome
#              dependence
# pooling      per-cell Q standardised by its own permutation mean/sd -> Qz;
#              pooled statistic = mean Qz over cells, with the global null built
#              from the same country-level draws
# families     the full 373-predictor block, and each of the 18 domains
# min_n        primary keeps all areas; sensitivity drops areas with n_svy < 15
#
#   Rscript scripts/signal_probes/p5_omnibus_dense_signal.R
# -> results/tables/signal_probes/p5_omnibus_cells.csv
# -> results/tables/signal_probes/p5_omnibus_domains.csv
# -> results/tables/signal_probes/p5_omnibus_pooled.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")

STORE  <- "_targets_full"
B      <- 2000L
SEED   <- 20260941L
# P5_IMPUTE=1 keeps every predictor with >= 80% coverage on the analysis rows and
# median-imputes the rest, instead of dropping any column with a single NA. This
# matters because the 153-column DHS-derived block has 0 columns fully present on
# Ghana's 75 survey districts (15% of districts missing, scattered), so a
# complete-case column rule costs Ghana 174 of 373 predictors - including the
# demographic/nutritional domains. Imputation is outcome-independent, so the
# permutation null remains exact.
IMPUTE <- identical(Sys.getenv("P5_IMPUTE"), "1")
SUF    <- if (IMPUTE) "_imputed" else ""
MIN_COV <- 0.80
OUTDIR <- "results/tables/signal_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS     <- intersect(MD$column, names(S))
domain_of <- setNames(MD$domain, MD$column)
domains   <- sort(unique(domain_of[PREDS]))
cat("shared predictors:", length(PREDS), "| domains:", length(domains), "\n")

BND2 <- readRDS("dashboard/data/admin2_boundaries.rds")
BND1 <- readRDS("dashboard/data/admin1_boundaries.rds")

# ── centroids, keyed the same way the analysis is keyed ─────────────────────
centroids <- function(bnd, keycols) {
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(bnd))))
  out <- as.data.frame(sf::st_drop_geometry(bnd)[, keycols, drop = FALSE])
  out$lon <- xy[, 1]; out$lat <- xy[, 2]
  out
}
CENT2 <- lapply(names(COUNTRIES), function(lc) centroids(BND2[[lc]], c("Admin1", "Admin2")))
CENT1 <- lapply(names(COUNTRIES), function(lc) centroids(BND1[[lc]], c("Admin1")))
names(CENT2) <- names(CENT1) <- COUNTRIES

# ── outcome cells ───────────────────────────────────────────────────────────
rknorm <- function(x) {
  ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && stats::sd(x[ok]) > 0)
    out[ok] <- stats::qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}

meta_t <- targets::tar_meta(store = STORE)
svy_names <- grep("^svy_admin2_(gambia|ghana|malawi|sierraleone)_",
                  meta_t$name, value = TRUE)

build_cells <- function(unit, min_n) {
  out <- list()
  for (tn in svy_names) {
    parts <- sub("^svy_admin2_", "", tn)
    lc <- sub("_(child|women)_.*$", "", parts)
    oc <- sub("^[a-z]+_", "", parts)
    cn <- COUNTRIES[[lc]]
    sv <- tryCatch(targets::tar_read_raw(tn, store = STORE), error = function(e) NULL)
    if (is.null(sv)) next
    sv <- sv[is.finite(sv$svy_prev) & is.finite(sv$n_svy) & sv$n_svy > 0, ]
    if (min_n > 0) sv <- sv[sv$n_svy >= min_n, ]
    if (sum(sv$svy_prev * sv$n_svy) < 10) next

    Sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)]
    if (unit == "admin2") {
      # PAIR KEY join - never join Admin-2 on the bare name
      m <- dplyr::inner_join(sv[, c("Admin1", "Admin2", "svy_prev", "n_svy")],
                             Sc, by = c("Admin1", "Admin2"))
      m <- dplyr::inner_join(m, CENT2[[cn]], by = c("Admin1", "Admin2"))
      y <- m$svy_prev; nn <- m$n_svy
      X <- as.matrix(m[, PREDS, drop = FALSE])
    } else {
      a1 <- sv |> group_by(Admin1) |>
        summarise(y = sum(svy_prev * n_svy) / sum(n_svy), n = sum(n_svy),
                  .groups = "drop")
      # predictors aggregated over the covariate table's OWN districts
      xa <- Sc |> group_by(Admin1) |>
        summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
      m <- dplyr::inner_join(a1, xa, by = "Admin1")
      m <- dplyr::inner_join(m, CENT1[[cn]], by = "Admin1")
      y <- m$y; nn <- m$n
      X <- as.matrix(m[, PREDS, drop = FALSE])
    }
    if (length(y) < 8 || stats::sd(y) == 0) next

    # rank-normalise predictors within country over the analysis rows, drop
    # columns that are constant or unusable HERE (they carry no information)
    Xr <- apply(X, 2, rknorm)
    if (IMPUTE) {
      cov_j <- colMeans(is.finite(Xr))
      keep  <- cov_j >= MIN_COV &
        apply(Xr, 2, function(z) sum(is.finite(z)) > 2 &&
                stats::sd(z[is.finite(z)]) > 0)
      Xr[!is.finite(Xr)] <- 0          # 0 = median on the rank-normal scale
    } else {
      keep <- apply(Xr, 2, function(z) all(is.finite(z)) && stats::sd(z) > 0)
    }
    if (sum(keep) < 20) next
    out[[paste(lc, oc, sep = "|")]] <- list(
      country = cn, outcome = oc, unit = unit, min_n = min_n,
      y = y, n = nn, yr = rknorm(y),
      X = Xr[, keep, drop = FALSE], cols = PREDS[keep],
      lat = m$lat, lon = m$lon, n_areas = length(y),
      blk = if (unit == "admin2") as.character(m$Admin1) else as.character(m$Admin1))
  }
  out
}

# ── spatial bases ───────────────────────────────────────────────────────────
mem_basis <- function(lat, lon, k = 4L, m_max = 8L) {
  n <- length(lat)
  if (n < 12) return(NULL)
  D <- as.matrix(stats::dist(cbind(lon, lat)))
  diag(D) <- Inf
  kk <- min(k, n - 1L)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) A[i, order(D[i, ])[seq_len(kk)]] <- 1
  A <- (A + t(A)) / 2                      # symmetrise
  P <- diag(n) - matrix(1 / n, n, n)
  M <- P %*% A %*% P
  ev <- eigen(M, symmetric = TRUE)
  pos <- which(ev$values > 1e-8)
  if (!length(pos)) return(NULL)
  m <- min(m_max, length(pos), max(2L, floor(n / 8)))
  ev$vectors[, pos[seq_len(m)], drop = FALSE]
}

# within-Admin-1 restricted permutation: the SAME null 16_bivariate_fdr.R uses
# for its headline "0 of 294 survive FDR" result. Reused here so that the only
# difference from that scan is the STATISTIC (omnibus vs per-predictor), which
# isolates how much of "0 of 294" is a power limitation.
block_perm_matrix <- function(blk, B) {
  n <- length(blk); idx <- split(seq_len(n), blk)
  usable <- vapply(idx, length, 1L) > 1L
  if (!any(usable)) return(NULL)
  replicate(B, {
    p <- seq_len(n)
    for (g in idx[usable]) p[g] <- g[sample.int(length(g))]
    p
  })
}

basis_for <- function(cell, variant) {
  if (variant == "S0_none")      return(NULL)
  if (variant == "S3_blockperm") return(NULL)
  if (variant == "S1_latlon") return(cbind(lat = cell$lat, lon = cell$lon))
  mem_basis(cell$lat, cell$lon)
}

# residualise a matrix (or vector) on a basis; returns centred columns
resid_on <- function(Z, Bm) {
  Z <- as.matrix(Z)
  Z <- sweep(Z, 2, colMeans(Z), "-")
  if (is.null(Bm)) return(Z)
  Bm <- sweep(as.matrix(Bm), 2, colMeans(Bm), "-")
  qq <- qr(Bm)
  Z - qr.fitted(qq, Z)
}
unitcol <- function(Z) {                    # scale columns to unit norm
  s <- sqrt(colSums(Z^2)); s[s <= 0] <- NA_real_
  sweep(Z, 2, s, "/")
}

# ── the test, vectorised over B permutations ─────────────────────────────────
# Xs: n x p, unit-norm columns.  Ys: n x (1+B), unit-norm columns
#   column 1 = observed, 2..(B+1) = permuted.  R = t(Xs) %*% Ys is p x (1+B).
cell_stats <- function(cell, variant, perm_index) {
  Bm <- basis_for(cell, variant)
  if (!variant %in% c("S0_none", "S3_blockperm")) {
    if (is.null(Bm)) return(NULL)
    if (cell$n_areas - ncol(Bm) - 1 < 8) return(NULL)
  }
  Xs <- unitcol(resid_on(cell$X, Bm))
  ok <- which(is.finite(colSums(Xs)))
  if (length(ok) < 20) return(NULL)
  Xs <- Xs[, ok, drop = FALSE]; cols <- cell$cols[ok]

  yres <- as.numeric(resid_on(cell$yr, Bm))
  Ymat <- cbind(yres, matrix(yres[perm_index], nrow = length(yres)))
  Ys <- unitcol(Ymat)
  R  <- crossprod(Xs, Ys)                   # p x (1+B)
  R[!is.finite(R)] <- 0

  R2  <- R^2
  Qall  <- colMeans(R2)
  Tmax  <- apply(abs(R), 2, max)
  dcols <- domain_of[cols]
  Qdom <- t(vapply(domains, function(dm) {
    idx <- which(dcols == dm)
    if (length(idx) < 3) return(rep(NA_real_, ncol(R2)))
    colMeans(R2[idx, , drop = FALSE])
  }, numeric(ncol(R2))))
  rownames(Qdom) <- domains
  list(p_used = length(cols), n_basis = if (is.null(Bm)) 0L else ncol(Bm),
       Qall = Qall, Tmax = Tmax, Qdom = Qdom,
       ndom = vapply(domains, function(dm) sum(dcols == dm), 1L))
}

pval <- function(v) (1 + sum(v[-1] >= v[1], na.rm = TRUE)) /
                    (1 + sum(is.finite(v[-1])))
zscore <- function(v) {
  mu <- mean(v[-1], na.rm = TRUE); sd0 <- stats::sd(v[-1], na.rm = TRUE)
  if (!is.finite(sd0) || sd0 <= 0) return(NA_real_)
  (v[1] - mu) / sd0
}

# ── run ─────────────────────────────────────────────────────────────────────
VARIANTS <- c("S0_none", "S1_latlon", "S2_mem", "S3_blockperm")
rows_cell <- list(); rows_dom <- list(); pooled_rows <- list()

for (unit in c("admin2", "admin1")) {
  for (min_n in c(0L, 15L)) {
    cells <- build_cells(unit, min_n)
    if (!length(cells)) next
    cat("\n== unit =", unit, "| min_n =", min_n, "| cells =", length(cells), "\n")

    # one permutation matrix per (country, n_areas), reused across that
    # country's outcome cells so the pooled null respects dependence
    perm_cache <- list(); bperm_cache <- list()
    for (kx in names(cells)) {
      cl <- cells[[kx]]
      pk <- paste(cl$country, cl$n_areas)
      if (is.null(perm_cache[[pk]]))
        perm_cache[[pk]] <- replicate(B, sample.int(cl$n_areas))
      bk <- paste(cl$country, paste(cl$blk, collapse = "|"))
      if (!bk %in% names(bperm_cache))
        bperm_cache[[bk]] <- block_perm_matrix(cl$blk, B)
    }
    perm_for <- function(cl, variant) {
      if (variant == "S3_blockperm")
        return(bperm_cache[[paste(cl$country, paste(cl$blk, collapse = "|"))]])
      perm_cache[[paste(cl$country, cl$n_areas)]]
    }

    Qz_by_variant <- list()
    for (variant in VARIANTS) {
      if (variant == "S3_blockperm" && unit == "admin1") next
      qz <- list()
      for (kx in names(cells)) {
        cl <- cells[[kx]]
        pmx <- perm_for(cl, variant)
        if (is.null(pmx)) next
        st <- cell_stats(cl, variant, pmx)
        if (is.null(st)) next
        rows_cell[[paste(unit, min_n, variant, kx)]] <- data.frame(
          unit = unit, min_n = min_n, spatial = variant,
          country = cl$country, outcome = cl$outcome,
          n_areas = cl$n_areas, p_used = st$p_used, n_basis = st$n_basis,
          Q_obs = st$Qall[1], Q_null_mean = mean(st$Qall[-1]),
          Q_ratio = st$Qall[1] / mean(st$Qall[-1]),
          Q_z = zscore(st$Qall), p_perm_Q = pval(st$Qall),
          Tmax_obs = st$Tmax[1], p_perm_Tmax = pval(st$Tmax),
          row.names = NULL)
        qz[[kx]] <- st$Qall                  # keep the full vector for pooling
        for (dm in domains) {
          v <- st$Qdom[dm, ]
          if (!is.finite(v[1])) next
          rows_dom[[paste(unit, min_n, variant, kx, dm)]] <- data.frame(
            unit = unit, min_n = min_n, spatial = variant,
            country = cl$country, outcome = cl$outcome, domain = dm,
            n_pred = st$ndom[[dm]], Q_obs = v[1],
            Q_ratio = v[1] / mean(v[-1]), Q_z = zscore(v),
            p_perm = pval(v), row.names = NULL)
        }
      }
      if (length(qz) < 3) next
      # pooled: standardise each cell's Q trajectory by its own null, then average
      Z <- do.call(rbind, lapply(qz, function(v) {
        mu <- mean(v[-1], na.rm = TRUE); s <- stats::sd(v[-1], na.rm = TRUE)
        if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(v)))
        (v - mu) / s
      }))
      Z <- Z[stats::complete.cases(Z), , drop = FALSE]
      if (nrow(Z) < 3) next
      pooled <- colMeans(Z)
      pooled_rows[[paste(unit, min_n, variant)]] <- data.frame(
        unit = unit, min_n = min_n, spatial = variant, n_cells = nrow(Z),
        pooled_Qz = pooled[1],
        p_global = (1 + sum(pooled[-1] >= pooled[1])) / length(pooled),
        n_cells_p05 = sum(vapply(qz, function(v) pval(v) < 0.05, TRUE)),
        row.names = NULL)
      Qz_by_variant[[variant]] <- pooled[1]
    }
    cat("   pooled Qz:",
        paste(sprintf("%s=%.2f", names(Qz_by_variant),
                      unlist(Qz_by_variant)), collapse = "  "), "\n")
  }
}

cellsdf <- bind_rows(rows_cell)
domdf   <- bind_rows(rows_dom)
pooldf  <- bind_rows(pooled_rows)

# BH across cells within (unit, min_n, spatial)
cellsdf <- cellsdf |> group_by(unit, min_n, spatial) |>
  mutate(q_bh_Q = p.adjust(p_perm_Q, "BH"),
         q_bh_Tmax = p.adjust(p_perm_Tmax, "BH")) |> ungroup()
domdf <- domdf |> group_by(unit, min_n, spatial, outcome) |>
  mutate(q_bh = p.adjust(p_perm, "BH")) |> ungroup()

write.csv(cellsdf, file.path(OUTDIR, paste0("p5_omnibus_cells", SUF, ".csv")), row.names = FALSE)
write.csv(domdf,   file.path(OUTDIR, paste0("p5_omnibus_domains", SUF, ".csv")), row.names = FALSE)
write.csv(pooldf,  file.path(OUTDIR, paste0("p5_omnibus_pooled", SUF, ".csv")), row.names = FALSE)

cat("\n================ POOLED (all cells) ================\n")
print(as.data.frame(pooldf), row.names = FALSE)

cat("\n======== Q vs Tmax: how many cells reach p < 0.05 ========\n")
print(as.data.frame(cellsdf |> group_by(unit, min_n, spatial) |>
  summarise(cells = n(),
            Q_p05 = sum(p_perm_Q < 0.05), Q_q05 = sum(q_bh_Q < 0.05),
            Tmax_p05 = sum(p_perm_Tmax < 0.05), Tmax_q05 = sum(q_bh_Tmax < 0.05),
            median_Q_ratio = round(median(Q_ratio), 3), .groups = "drop")),
  row.names = FALSE)

cat("\n======== primary cells (admin2, min_n=0) ========\n")
print(as.data.frame(cellsdf |>
  filter(unit == "admin2", min_n == 0) |>
  select(spatial, country, outcome, n_areas, p_used, Q_ratio, Q_z,
         p_perm_Q, p_perm_Tmax) |>
  arrange(spatial, p_perm_Q)), row.names = FALSE, digits = 3)

cat("\nDONE\n")
