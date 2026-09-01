# =============================================================================
# scripts/signal_probes/p3_loco_domain_index.R
#
# P3a. Does the between-region signal TRANSPORT? A zero-tuning check.
#
# For each held-out country, a linear index over the 18 domain scores is built
# using ONLY the other three countries: per-predictor sign alignment (pooled
# PC1 of the training countries' rank-normalized covariate rows), and a weight
# per domain equal to the training-countries meta-analytic mean z for the
# held-out outcome. No hyperparameters, no fitting on the held-out country;
# the held-out country contributes only its own covariates (rank-normalized
# within itself — no outcome information). The index is scored on Spearman
# rank correlation against the held-out country's observed Admin-1 prevalence.
# For a country with no survey, the honest comparator is a flat map (r = 0),
# so any consistent positive rank correlation is pure transported information.
#
# P3b. Robustness of P1 to smooth spatial gradients: the same per-cell
# association scan after residualising outcome and predictors on region
# centroid latitude+longitude, for the two countries with enough regions to
# support it (Ghana 16, Malawi 27). Freedman-Lane style permutation on
# reduced-model residuals.
#
#   Rscript scripts/signal_probes/p3_loco_domain_index.R
# -> results/tables/signal_probes/p3_loco_domain_index.csv
# -> results/tables/signal_probes/p3_loco_summary.csv
# -> results/tables/signal_probes/p3b_latlon_partialed_domains.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd('C:/Users/andre/OneDrive/Documents/mn-prediction')

STORE  <- "_targets_full"
B      <- 2000L
SEED   <- 20260933L
OUTDIR <- "results/tables/signal_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- intersect(MD$column, names(S))
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")
domain_of <- setNames(MD$domain, MD$column)
domains <- sort(unique(MD$domain))

# ── binary outcome cells at Admin-1 (as in P1) ──────────────────────────────
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
  n_cases <- sum(sv$svy_prev * sv$n_svy)
  a1 <- sv |> group_by(Admin1) |>
    summarise(y = sum(svy_prev * n_svy) / sum(n_svy), n = sum(n_svy),
              .groups = "drop")
  if (n_cases < 10 || nrow(a1) < 4 || sd(a1$y) == 0) next
  cells[[paste(cc, oc, sep = "|")]] <- list(country = COUNTRIES[[cc]],
                                            outcome = oc, a1 = a1)
}
ctry <- vapply(cells, function(x) x$country, "")
outc <- vapply(cells, function(x) x$outcome, "")
cat("cells:", length(cells), "\n")

rknorm <- function(x) {
  ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && sd(x[ok]) > 0)
    out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}
# rank-normalized covariates per country over the covariate geography
RN_by_country <- lapply(COUNTRIES, function(cn)
  as.data.frame(lapply(S[S$country == cn, PREDS], rknorm)))
A1_by_country <- lapply(COUNTRIES, function(cn) S$Admin1[S$country == cn])

sign_from <- function(train_countries) {
  RN <- do.call(rbind, RN_by_country[names(COUNTRIES)[match(train_countries, COUNTRIES)]])
  sg <- setNames(rep(1, length(PREDS)), PREDS)
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
    sg[cols] <- ifelse(ld >= 0, 1, -1)
  }
  sg
}
domain_scores_a1 <- function(cn, sg) {          # cn = pretty country name
  lc <- names(COUNTRIES)[match(cn, COUNTRIES)]
  rn <- sweep(as.matrix(RN_by_country[[lc]]), 2, sg[PREDS], "*")
  ds <- sapply(domains, function(dm) {
    cols <- which(domain_of[PREDS] == dm)
    if (!length(cols)) return(rep(NA_real_, nrow(rn)))
    rowMeans(rn[, cols, drop = FALSE], na.rm = TRUE)
  })
  ds <- as.data.frame(ds); names(ds) <- domains
  ds$Admin1 <- A1_by_country[[lc]]
  ds |> group_by(Admin1) |>
    summarise(across(all_of(domains), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}
fisher_z <- function(r, n) {
  r <- pmin(pmax(r, -0.999), 0.999)
  0.5 * log((1 + r) / (1 - r)) * sqrt(pmax(n - 3, 1))
}
dl_mu <- function(z) {                          # z: vector over countries
  z <- z[is.finite(z)]; k <- length(z)
  if (!k) return(NA_real_)
  mu_fe <- mean(z)
  Q <- sum((z - mu_fe)^2)
  tau2 <- max(0, (Q - (k - 1)) / max(k - 1e-9, 1))
  w <- rep(1 / (1 + tau2), k)
  sum(w * z) / sum(w)
}

# ── P3a: LOCO index ─────────────────────────────────────────────────────────
rows <- list()
held_scores <- list()   # per cell: list(y=..., idx=...) for global permutation
for (hold in COUNTRIES) {
  train <- setdiff(COUNTRIES, hold)
  sg <- sign_from(train)
  DS <- lapply(COUNTRIES, function(cn) domain_scores_a1(cn, sg))
  names(DS) <- COUNTRIES
  hold_cells <- names(cells)[ctry == hold]
  for (kx in hold_cells) {
    on <- outc[kx]
    # training weights: per training country, z of each domain for this outcome
    zmat <- sapply(domains, function(dm) {
      vapply(train, function(tc) {
        tk <- names(cells)[ctry == tc & outc == on]
        if (!length(tk)) return(NA_real_)
        cl <- cells[[tk]]
        m <- inner_join(cl$a1, DS[[tc]], by = "Admin1")
        x <- m[[dm]]
        if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) return(NA_real_)
        ok <- is.finite(x)
        fisher_z(cor(m$y[ok], x[ok], method = "spearman"), sum(ok))
      }, 0)
    })
    w_d <- apply(zmat, 2, dl_mu)                 # meta mean z per domain
    w_d[!is.finite(w_d)] <- 0
    # index on held-out regions
    cl <- cells[[kx]]
    m <- inner_join(cl$a1, DS[[hold]], by = "Admin1")
    Xd <- as.matrix(m[, domains, drop = FALSE]); Xd[!is.finite(Xd)] <- 0
    idx <- as.numeric(Xd %*% w_d)
    if (sd(idx) == 0) next
    r_obs <- cor(m$y, idx, method = "spearman")
    p_perm <- (1 + sum(replicate(B, abs(cor(sample(m$y), idx,
                 method = "spearman"))) >= abs(r_obs))) / (B + 1)
    held_scores[[kx]] <- list(y = m$y, idx = idx, country = hold)
    rows[[kx]] <- data.frame(held_out = hold, outcome = on,
                             n_regions = nrow(m), n_train_countries =
                               sum(apply(zmat, 1, function(z) any(is.finite(z)))),
                             spearman = r_obs, p_perm = p_perm)
  }
}
res <- bind_rows(rows)
write.csv(res, file.path(OUTDIR, "p3_loco_domain_index.csv"), row.names = FALSE)
print(res, row.names = FALSE)

# global permutation: mean Spearman across all held-out cells, one region
# permutation per country per draw (preserves cross-outcome dependence)
obs_mean <- mean(res$spearman)
perm_means <- numeric(B)
for (b in seq_len(B)) {
  perms <- list()
  vals <- vapply(names(held_scores), function(kx) {
    hs <- held_scores[[kx]]; np <- length(hs$y)
    pk <- paste(hs$country, np)
    if (is.null(perms[[pk]])) perms[[pk]] <<- sample(np)
    cor(hs$y[perms[[pk]]], hs$idx, method = "spearman")
  }, 0)
  perm_means[b] <- mean(vals)
}
p_global <- (1 + sum(perm_means >= obs_mean)) / (B + 1)
summ <- data.frame(
  n_cells = nrow(res), mean_spearman = obs_mean,
  median_spearman = median(res$spearman),
  n_positive = sum(res$spearman > 0),
  n_perm_p_lt_05 = sum(res$p_perm < 0.05),
  p_global_onesided = p_global)
write.csv(summ, file.path(OUTDIR, "p3_loco_summary.csv"), row.names = FALSE)
cat("\nGLOBAL: mean spearman =", round(obs_mean, 3), " one-sided p =",
    signif(p_global, 3), "\n")

# ── P3b: lat/lon-partialed domain scan, Ghana + Malawi ──────────────────────
bl <- readRDS("dashboard/data/admin1_boundaries.rds")
cent <- list()
for (lc in c("ghana", "malawi")) {
  bb <- bl[[lc]]
  ct <- suppressWarnings(sf::st_coordinates(sf::st_centroid(bb)))
  nmcol <- intersect(c("Admin1", "NAME_1", "shapeName", "name"), names(bb))[1]
  cent[[lc]] <- data.frame(Admin1 = as.character(bb[[nmcol]]),
                           lon = ct[, 1], lat = ct[, 2])
}
sg_all <- sign_from(COUNTRIES)     # full-sample signs fine for a robustness check
DS_all <- lapply(COUNTRIES, function(cn) domain_scores_a1(cn, sg_all))
names(DS_all) <- COUNTRIES
p3b <- list()
for (lc in c("ghana", "malawi")) {
  cn <- COUNTRIES[[lc]]
  for (kx in names(cells)[ctry == cn]) {
    cl <- cells[[kx]]
    m <- inner_join(cl$a1, DS_all[[cn]], by = "Admin1") |>
      inner_join(cent[[lc]], by = "Admin1")
    if (nrow(m) < 8) next
    ry <- resid(lm(y ~ lat + lon, data = m))
    for (dm in domains) {
      x <- m[[dm]]
      if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) next
      x[!is.finite(x)] <- mean(x, na.rm = TRUE)
      rx <- resid(lm(x ~ lat + lon, data = m))
      r_obs <- cor(ry, rx)
      # Freedman-Lane: permute reduced-model residuals of y
      pd <- replicate(B, abs(cor(sample(ry), rx)))
      p3b[[paste(kx, dm)]] <- data.frame(
        country = cn, outcome = cl$outcome, domain = dm, n = nrow(m),
        r_partial = r_obs,
        p_perm = (1 + sum(pd >= abs(r_obs))) / (B + 1))
    }
  }
}
p3b <- bind_rows(p3b)
# pool per domain across cells within country (mean z), then across the two
# countries (Stouffer)
pool3b <- p3b |>
  mutate(z = fisher_z(r_partial, n)) |>
  group_by(domain, country) |> summarise(zc = mean(z), .groups = "drop") |>
  group_by(domain) |>
  summarise(k = n(), stouffer_z = sum(zc) / sqrt(n()), .groups = "drop") |>
  arrange(desc(abs(stouffer_z)))
write.csv(p3b, file.path(OUTDIR, "p3b_latlon_partialed_cells.csv"), row.names = FALSE)
write.csv(pool3b, file.path(OUTDIR, "p3b_latlon_partialed_domains.csv"), row.names = FALSE)
cat("\nP3b pooled (lat/lon partialed, Ghana+Malawi):\n")
print(as.data.frame(pool3b), row.names = FALSE)
cat("\nDONE\n")
