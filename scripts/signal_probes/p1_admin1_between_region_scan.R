# =============================================================================
# scripts/signal_probes/p1_admin1_between_region_scan.R
#
# P1. Is there BETWEEN-REGION covariate signal, and which predictors carry it?
#
# WHY THIS TEST DOES NOT EXIST YET
# --------------------------------
# The pooled consistency scan (WS-A) permutes outcomes WITHIN Admin-1 region,
# so its null retains every predictor's between-region association in every
# draw: it can only ever detect within-region signal. The project's own WS4a
# finding is that most of what signal exists sits BETWEEN regions — so the
# existing scan is structurally blind to the signal class most likely to be
# real. This scan tests that class directly.
#
# WHY IT IS LINKAGE-ROBUST
# ------------------------
# The covariate tables sit on a different Admin-2 vintage than the survey
# outcomes (Ghana 260 vs 75 rows, Malawi 243 vs 87, Gambia 37 vs 30), so any
# Admin-2 name join risks silent polygon mismatch. Here each side is
# aggregated to Admin-1 on its OWN geography and joined on Admin-1 names,
# which match exactly in all four countries (verified upstream). No Admin-2
# join is performed at all.
#
# DESIGN
# ------
# unit         Admin-1 region (Gambia 6, Ghana 16, Malawi 27, SierraLeone 4)
# outcome      n-weighted mean of district svy_prev within region
# predictor    unweighted mean over the covariate table's own districts
# statistic    Spearman r per country-outcome cell -> Fisher z * sqrt(n-3)
# pooling      per outcome: DerSimonian-Laird meta across countries
#              (one cell per country per outcome, country = cluster)
# null         permute the REGION outcome vector across regions within
#              country; ONE permutation per country per draw, applied to all
#              that country's cells, preserving cross-outcome dependence
# multiplicity family-wise max-|meta z| calibration over predictors (the same
#              calibration WS-A used), plus per-predictor permutation p
# families     373 shared predictors, and 18 domain scores (mean of
#              sign-aligned rank-normalized members; signs from pooled PC1)
#
# CAVEAT stated up front: a between-region association cannot distinguish
# "covariate carries information" from "covariate tracks a smooth spatial /
# socioeconomic gradient that also tracks deficiency". For the question being
# asked — do these covariates carry PREDICTIVE information, transportable
# across countries — that distinction does not matter; cross-country sign
# replication is the evidence standard, and the meta-analysis enforces it.
#
#   Rscript scripts/signal_probes/p1_admin1_between_region_scan.R
# -> results/tables/signal_probes/p1_admin1_scan_predictors.csv
# -> results/tables/signal_probes/p1_admin1_scan_domains.csv
# -> results/tables/signal_probes/p1_admin1_scan_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd('C:/Users/andre/OneDrive/Documents/mn-prediction')

STORE  <- "_targets_full"
B      <- 2000L
SEED   <- 20260931L
OUTDIR <- "results/tables/signal_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

# ── outcome cells: aggregate survey districts to Admin-1 ────────────────────
meta_t <- targets::tar_meta(store = STORE)
svy_names <- grep("^svy_admin2_(gambia|ghana|malawi|sierraleone)_",
                  meta_t$name, value = TRUE)
cells <- list()
for (tn in svy_names) {
  parts <- sub("^svy_admin2_", "", tn)
  cc <- sub("_(child|women)_.*$", "", parts)
  oc <- sub("^[a-z]+_", "", parts)
  sv <- tryCatch(targets::tar_read_raw(tn, store = STORE), error = function(e) NULL)
  if (is.null(sv)) next
  sv <- sv[is.finite(sv$svy_prev) & is.finite(sv$n_svy) & sv$n_svy > 0, ]
  # guard against degenerate cells (e.g. SierraLeone women_b12: 4 cases total)
  n_cases <- sum(sv$svy_prev * sv$n_svy)
  a1 <- sv |>
    group_by(Admin1) |>
    summarise(y = sum(svy_prev * n_svy) / sum(n_svy), n = sum(n_svy),
              .groups = "drop")
  if (n_cases < 10 || nrow(a1) < 4 || sd(a1$y) == 0) next
  cells[[paste(cc, oc, sep = "|")]] <- list(country = COUNTRIES[[cc]],
                                            outcome = oc, a1 = a1)
}
cat("cells retained:", length(cells), "\n")

# ── predictors: aggregate covariate table to Admin-1 on its own geography ───
X_a1 <- list()
for (cn in COUNTRIES) {
  sc <- S[S$country == cn, c("Admin1", PREDS)]
  X_a1[[cn]] <- sc |> group_by(Admin1) |>
    summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

# ── domain scores: sign-align via pooled PC1 of within-country rank-normal ──
rknorm <- function(x) {
  ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && sd(x[ok]) > 0)
    out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}
# stack rank-normalized within-country matrices over the covariate geography
RN <- do.call(rbind, lapply(COUNTRIES, function(cn) {
  sc <- S[S$country == cn, PREDS]
  as.data.frame(lapply(sc, rknorm))
}))
domain_of <- setNames(MD$domain, MD$column)
domains <- sort(unique(MD$domain))
sign_of <- setNames(rep(1, length(PREDS)), PREDS)
for (dm in domains) {
  cols <- PREDS[domain_of[PREDS] == dm]
  Mx <- RN[, cols, drop = FALSE]
  keep <- colSums(is.finite(as.matrix(Mx))) > 0.5 * nrow(Mx)
  cols <- cols[keep]
  if (length(cols) < 2) next
  Mx <- as.matrix(RN[, cols, drop = FALSE])
  Mx[!is.finite(Mx)] <- 0
  pc <- tryCatch(prcomp(Mx, center = TRUE, scale. = FALSE), error = function(e) NULL)
  if (is.null(pc)) next
  ld <- pc$rotation[, 1]
  if (mean(ld > 0) < 0.5) ld <- -ld        # majority-positive orientation
  sign_of[cols] <- ifelse(ld >= 0, 1, -1)
}
# domain score at Admin-1: mean of sign-aligned rank-normalized members,
# aggregated over the covariate table's own districts
D_a1 <- list()
for (cn in COUNTRIES) {
  sc <- S[S$country == cn, c("Admin1", PREDS)]
  rn <- as.data.frame(lapply(sc[PREDS], rknorm))
  rn <- sweep(as.matrix(rn), 2, sign_of[PREDS], "*")
  ds <- sapply(domains, function(dm) {
    cols <- which(domain_of[PREDS] == dm)
    if (!length(cols)) return(rep(NA_real_, nrow(rn)))
    rowMeans(rn[, cols, drop = FALSE], na.rm = TRUE)
  })
  ds <- as.data.frame(ds); names(ds) <- domains
  ds$Admin1 <- sc$Admin1
  D_a1[[cn]] <- ds |> group_by(Admin1) |>
    summarise(across(all_of(domains), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

# ── assemble per-cell aligned matrices (region x predictor), pre-ranked ─────
prep_cell <- function(cell, Xtab, cols) {
  m <- inner_join(cell$a1, Xtab, by = "Admin1")
  y <- m$y; X <- as.matrix(m[, cols, drop = FALSE])
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z) > 0)
  list(yr = rank(y), Xr = apply(X[, keep, drop = FALSE], 2, rank),
       cols = cols[keep], n = nrow(m))
}
cellsP <- lapply(cells, function(cl) prep_cell(cl, X_a1[[cl$country]], PREDS))
cellsD <- lapply(cells, function(cl) prep_cell(cl, D_a1[[cl$country]], domains))

spearman_block <- function(Xr, yr) {
  Xc <- sweep(Xr, 2, colMeans(Xr), "-"); yc <- yr - mean(yr)
  r <- as.numeric(crossprod(Xc, yc)) / (sqrt(colSums(Xc^2)) * sqrt(sum(yc^2)))
  r[!is.finite(r)] <- NA_real_; r
}
fisher_z <- function(r, n) {
  r <- pmin(pmax(r, -0.999), 0.999)
  0.5 * log((1 + r) / (1 - r)) * sqrt(pmax(n - 3, 1))
}

# per-cell z for all predictors under a given region permutation (NULL = obs)
cell_z <- function(cp, perm = NULL) {
  yr <- if (is.null(perm)) cp$yr else cp$yr[perm]
  setNames(fisher_z(spearman_block(cp$Xr, yr), cp$n), cp$cols)
}

# DL meta across countries for a matrix of country z-scores (rows=country)
dl_meta <- function(zmat) {           # zmat: country x predictor, NA allowed
  k  <- colSums(is.finite(zmat))
  zm <- ifelse(is.finite(zmat), zmat, 0)
  fin <- is.finite(zmat) * 1
  mu_fe <- colSums(zm) / pmax(k, 1)
  Q  <- colSums(fin * (sweep(zmat, 2, mu_fe, "-"))^2, na.rm = TRUE)
  tau2 <- pmax(0, (Q - (k - 1)) / pmax(k - 1e-9, 1))   # vi = 1 for a z-score
  w  <- fin / (1 + matrix(tau2, nrow(zmat), ncol(zmat), byrow = TRUE))
  mu <- colSums(w * zm, na.rm = TRUE) / colSums(w)
  se <- 1 / sqrt(colSums(w))
  list(z = mu / se, mu = mu, tau2 = tau2, k = k)
}

run_family <- function(cellsX, cols, pool) {
  # pool: named list outcome-group -> character vector of cell keys
  cellkeys <- names(cellsX)
  ctry <- vapply(cells[cellkeys], function(x) x$country, "")
  obs_z <- lapply(cellsX, cell_z)
  # observed meta per group
  make_zmat <- function(zlist, keys) {
    zm <- matrix(NA_real_, length(COUNTRIES), length(cols),
                 dimnames = list(COUNTRIES, cols))
    for (cn in COUNTRIES) {
      ks <- keys[ctry[keys] == cn]
      if (!length(ks)) next
      acc <- matrix(NA_real_, length(ks), length(cols),
                    dimnames = list(ks, cols))
      for (kx in ks) acc[kx, names(zlist[[kx]])] <- zlist[[kx]]
      zm[cn, ] <- colMeans(acc, na.rm = TRUE)
    }
    zm
  }
  obs_meta <- lapply(pool, function(keys) dl_meta(make_zmat(obs_z, keys)))
  # permutation null
  perm_absz <- lapply(pool, function(...) matrix(NA_real_, B, length(cols)))
  for (b in seq_len(B)) {
    perms <- lapply(COUNTRIES, function(cn) NULL)
    z_b <- list()
    for (kx in cellkeys) {
      cn <- ctry[kx]; np <- cellsX[[kx]]$n
      if (is.null(perms[[cn]]) || length(perms[[cn]]) != np)
        perms[[cn]] <- sample(np)               # one draw per country per b
      z_b[[kx]] <- cell_z(cellsX[[kx]], perms[[cn]])
    }
    for (g in names(pool))
      perm_absz[[g]][b, ] <- abs(dl_meta(make_zmat(z_b, pool[[g]]))$z)
  }
  out <- list()
  for (g in names(pool)) {
    om <- obs_meta[[g]]
    pm <- perm_absz[[g]]
    p_perm <- vapply(seq_along(cols), function(j)
      (1 + sum(pm[, j] >= abs(om$z[j]), na.rm = TRUE)) / (B + 1), 0)
    maxdist <- apply(pm, 1, max, na.rm = TRUE)
    p_fwer <- vapply(abs(om$z), function(zz)
      (1 + sum(maxdist >= zz)) / (B + 1), 0)
    zmat <- make_zmat(obs_z, pool[[g]])
    sgn  <- colSums(sign(zmat) == rep(sign(om$mu), each = nrow(zmat)),
                    na.rm = TRUE)
    out[[g]] <- data.frame(group = g, predictor = cols, meta_z = om$z,
                           mu = om$mu, tau2 = om$tau2, k_countries = om$k,
                           sign_agree = sgn, p_perm = p_perm,
                           p_fwer = p_fwer, row.names = NULL)
  }
  bind_rows(out)
}

outc <- vapply(cells, function(x) x$outcome, "")
pool <- c(
  split(names(cells), outc),                          # per-outcome groups
  list(iron = names(cells)[grepl("iron", outc)],
       vitA = names(cells)[grepl("vitA", outc)],
       all  = names(cells))
)
pool <- pool[vapply(pool, length, 1L) >= 3]

cat("running predictor family (", length(PREDS), "predictors, B =", B, ")\n")
resP <- run_family(cellsP, PREDS, pool)
resP$domain <- domain_of[resP$predictor]
cat("running domain family (", length(domains), "domains )\n")
resD <- run_family(cellsD, domains, pool)

resP <- resP |> group_by(group) |>
  mutate(q_bh_perm = p.adjust(p_perm, "BH")) |> ungroup() |> arrange(p_fwer)
resD <- resD |> group_by(group) |>
  mutate(q_bh_perm = p.adjust(p_perm, "BH")) |> ungroup() |> arrange(p_fwer)

write.csv(resP, file.path(OUTDIR, "p1_admin1_scan_predictors.csv"), row.names = FALSE)
write.csv(resD, file.path(OUTDIR, "p1_admin1_scan_domains.csv"), row.names = FALSE)

# summary: overdispersion of meta z per group (distributed-weak-signal check)
summ <- resP |> group_by(group) |>
  summarise(n_pred = n(),
            n_absz_gt2 = sum(abs(meta_z) > 2, na.rm = TRUE),
            n_absz_gt3 = sum(abs(meta_z) > 3, na.rm = TRUE),
            exp_gt2 = round(0.0455 * n(), 1), exp_gt3 = round(0.0027 * n(), 1),
            sd_metaz = sd(meta_z, na.rm = TRUE),
            n_fwer_05 = sum(p_fwer < 0.05, na.rm = TRUE),
            n_perm_01 = sum(p_perm < 0.01, na.rm = TRUE), .groups = "drop")
write.csv(summ, file.path(OUTDIR, "p1_admin1_scan_summary.csv"), row.names = FALSE)
print(as.data.frame(summ))
cat("\nTop 15 predictor rows by FWER p (group=all):\n")
print(head(resP[resP$group == "all",
                c("predictor", "domain", "meta_z", "k_countries", "sign_agree",
                  "p_perm", "p_fwer")], 15), row.names = FALSE)
cat("\nDomain results (group=all):\n")
print(resD[resD$group == "all",
           c("predictor", "meta_z", "k_countries", "sign_agree", "p_perm",
             "p_fwer", "q_bh_perm")], row.names = FALSE)
cat("\nDONE\n")
