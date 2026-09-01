# =============================================================================
# scripts/signal_probes/p4_admin1_continuous_scan.R
#
# P4. Between-region scan with the CONTINUOUS biomarker outcome.
#
# Same design as p1_admin1_between_region_scan.R (Admin-1 aggregation on each
# side's own geography; region-permutation null; DL meta across countries;
# max-|z| FWER calibration), but the outcome is the survey-weighted Admin-1
# mean of the log biomarker concentration rather than the deficiency
# prevalence. The project's own measurements (continuous_vs_binary.csv) show
# the continuous target carries about 3.5x the predictive skill of the binary
# one (median r 0.206 vs 0.058), so this is the higher-power version of P1.
#
# Biomarker direction: LOWER concentration = worse status for all outcomes
# here (RBP, ferritin, folate, B12). To keep signs comparable with P1 (where
# positive meta_z = covariate associated with MORE deficiency), the outcome is
# NEGATED before correlation.
#
#   Rscript scripts/signal_probes/p4_admin1_continuous_scan.R
# -> results/tables/signal_probes/p4_admin1_continuous_predictors.csv
# -> results/tables/signal_probes/p4_admin1_continuous_domains.csv
# -> results/tables/signal_probes/p4_admin1_continuous_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd('C:/Users/andre/OneDrive/Documents/mn-prediction')

STORE  <- "_targets_full"
B      <- 2000L
SEED   <- 20260932L
OUTDIR <- "results/tables/signal_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

source("R/config.R")
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }

S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")

# ── outcome cells: Admin-1 weighted mean of -log(biomarker) ─────────────────
cfgs <- get_country_configs()
cells <- list()
for (cn in names(cfgs)) {
  lc <- tolower(cn)
  if (!lc %in% names(COUNTRIES)) next
  cc <- cfgs[[cn]]
  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    if (!all(c(oc$continuous, "Admin1") %in% names(d))) next
    y <- num(d[[oc$continuous]])
    w <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d))
           num(d[[cc$weight_col]]) else rep(1, nrow(d))
    w[!is.finite(w) | w <= 0] <- 1
    t <- if (identical(oc$cutoff_scale, "log")) y else
         { if (any(y <= 0, na.rm = TRUE)) y[y <= 0] <- NA; log(y) }
    keep <- is.finite(t) & is.finite(w)
    if (sum(keep) < 100) next
    a1 <- data.frame(Admin1 = as.character(d$Admin1)[keep],
                     t = t[keep], w = w[keep]) |>
      group_by(Admin1) |>
      summarise(y = -wmean(t, w), n = dplyr::n(), .groups = "drop") |>
      filter(is.finite(y), n >= 20)
    if (nrow(a1) < 4 || sd(a1$y) == 0) next
    cells[[paste(lc, on, sep = "|")]] <- list(country = COUNTRIES[[lc]],
                                              outcome = on, a1 = a1)
  }
}
cat("cells retained:", length(cells), "\n")
for (kx in names(cells)) cat(" ", kx, "regions:", nrow(cells[[kx]]$a1), "\n")

# ── predictors and domain scores at Admin-1 (same as P1) ────────────────────
X_a1 <- list()
for (cn in COUNTRIES) {
  sc <- S[S$country == cn, c("Admin1", PREDS)]
  X_a1[[cn]] <- sc |> group_by(Admin1) |>
    summarise(across(all_of(PREDS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}
rknorm <- function(x) {
  ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && sd(x[ok]) > 0)
    out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}
RN <- do.call(rbind, lapply(COUNTRIES, function(cn) {
  as.data.frame(lapply(S[S$country == cn, PREDS], rknorm))
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
  Mx <- as.matrix(RN[, cols, drop = FALSE]); Mx[!is.finite(Mx)] <- 0
  pc <- tryCatch(prcomp(Mx, center = TRUE), error = function(e) NULL)
  if (is.null(pc)) next
  ld <- pc$rotation[, 1]
  if (mean(ld > 0) < 0.5) ld <- -ld
  sign_of[cols] <- ifelse(ld >= 0, 1, -1)
}
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

# ── scan machinery (identical to P1) ────────────────────────────────────────
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
cell_z <- function(cp, perm = NULL) {
  yr <- if (is.null(perm)) cp$yr else cp$yr[perm]
  setNames(fisher_z(spearman_block(cp$Xr, yr), cp$n), cp$cols)
}
dl_meta <- function(zmat) {
  k  <- colSums(is.finite(zmat))
  zm <- ifelse(is.finite(zmat), zmat, 0)
  fin <- is.finite(zmat) * 1
  mu_fe <- colSums(zm) / pmax(k, 1)
  Q  <- colSums(fin * (sweep(zmat, 2, mu_fe, "-"))^2, na.rm = TRUE)
  tau2 <- pmax(0, (Q - (k - 1)) / pmax(k - 1e-9, 1))
  w  <- fin / (1 + matrix(tau2, nrow(zmat), ncol(zmat), byrow = TRUE))
  mu <- colSums(w * zm, na.rm = TRUE) / colSums(w)
  se <- 1 / sqrt(colSums(w))
  list(z = mu / se, mu = mu, tau2 = tau2, k = k)
}
run_family <- function(cellsX, cols, pool) {
  cellkeys <- names(cellsX)
  ctry <- vapply(cells[cellkeys], function(x) x$country, "")
  obs_z <- lapply(cellsX, cell_z)
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
  perm_absz <- lapply(pool, function(...) matrix(NA_real_, B, length(cols)))
  for (b in seq_len(B)) {
    perms <- list(); z_b <- list()
    for (kx in cellkeys) {
      cn <- ctry[kx]; np <- cellsX[[kx]]$n
      pk <- paste(cn, np)
      if (is.null(perms[[pk]])) perms[[pk]] <- sample(np)
      z_b[[kx]] <- cell_z(cellsX[[kx]], perms[[pk]])
    }
    for (g in names(pool))
      perm_absz[[g]][b, ] <- abs(dl_meta(make_zmat(z_b, pool[[g]]))$z)
  }
  out <- list()
  for (g in names(pool)) {
    om <- obs_meta[[g]]; pm <- perm_absz[[g]]
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
pool <- c(split(names(cells), outc),
          list(iron = names(cells)[grepl("iron", outc)],
               vitA = names(cells)[grepl("vitA", outc)],
               all  = names(cells)))
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
write.csv(resP, file.path(OUTDIR, "p4_admin1_continuous_predictors.csv"), row.names = FALSE)
write.csv(resD, file.path(OUTDIR, "p4_admin1_continuous_domains.csv"), row.names = FALSE)

summ <- resP |> group_by(group) |>
  summarise(n_pred = n(),
            n_absz_gt2 = sum(abs(meta_z) > 2, na.rm = TRUE),
            n_absz_gt3 = sum(abs(meta_z) > 3, na.rm = TRUE),
            exp_gt2 = round(0.0455 * n(), 1), exp_gt3 = round(0.0027 * n(), 1),
            sd_metaz = sd(meta_z, na.rm = TRUE),
            n_fwer_05 = sum(p_fwer < 0.05, na.rm = TRUE),
            n_perm_01 = sum(p_perm < 0.01, na.rm = TRUE), .groups = "drop")
write.csv(summ, file.path(OUTDIR, "p4_admin1_continuous_summary.csv"), row.names = FALSE)
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
