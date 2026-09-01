# =============================================================================
# scripts/accuracy_impact/wsf7_confirmatory.R
#
# WS-F7. The confirmatory round: the SAME meta-analysis restricted to the
# pre-specified nutrition-proximal family.
#
# WHY THIS IS REPORTED SEPARATELY FROM THE 294-WIDE SCAN
# ------------------------------------------------------
# WS-A searched 294 predictors and found two surviving a family-wise permutation
# correction, neither of them nutrition-proximal. That is an exploratory scan and
# its multiplicity burden is large.
#
# The mechanism hypothesis deserves a fair test, which means a SMALL family fixed
# in advance: the nutrition-proximal columns built in WS-F1 and WS-F2 plus the
# curated sets over the harmonised covariates. A family of tens rather than
# hundreds carries a far smaller correction, so an association of the same size
# has a real chance of surviving here that it does not have in the wide scan.
#
# If nothing survives in a pre-specified family of this size, that is a much
# stronger negative than the wide scan can deliver.
#
#   Rscript scripts/accuracy_impact/wsf7_confirmatory.R
# -> results/tables/nutrition_proximal_confirmatory.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
NPERM <- as.integer(Sys.getenv("WSA_PERM", if (PROFILE == "smoke") "50" else "500"))
SEED <- 20260926L
SHARED <- c("child_iron", "child_vitA", "women_iron", "women_vitA")
TDIR <- here("results", "tables"); FDIR <- here("results", "figures")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

HH <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
NP <- suppressMessages(readr::read_csv(
  here("data", "covariates", "nutrition_proximal.csv"),
  show_col_types = FALSE)) |> as.data.frame()
SETS <- suppressMessages(readr::read_csv(
  here("data", "covariates", "nutrition_proximal_sets.csv"),
  show_col_types = FALSE)) |> as.data.frame()

# The pre-specified family: every column named by a curated set, from either
# source, deduplicated. Joining NP onto H keeps one row per district pair.
FAMILY <- unique(SETS$column)
npc <- intersect(grep("^np_", FAMILY, value = TRUE), names(NP))
H <- dplyr::left_join(HH, NP[, c("Admin1", "Admin2", npc)],
                      by = admin2_join_by(HH, NP))
if (nrow(H) != nrow(HH)) stop("nutrition-proximal join fanned rows")
COVS <- intersect(FAMILY, names(H))
cat(sprintf("[F7] pre-specified family: %d columns (%d nutrition-proximal, %d harmonised)
",
            length(COVS), sum(grepl("^np_", COVS)), sum(!grepl("^np_", COVS))))
REL <- read.csv(file.path(TDIR, "reliability_empirical.csv"), stringsAsFactors = FALSE)
REL <- REL[REL$scheme == "within", ]
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))

.logit <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
# Residualise a matrix on region indicators. Returns the residuals; a column
# that is constant within every region becomes all zeros and is dropped later.
.resid_on <- function(M, g) {
  g <- as.factor(g)
  for (lv in levels(g)) { i <- which(g == lv)
    if (length(i) > 1) M[i, ] <- sweep(M[i, , drop = FALSE], 2,
                                       colMeans(M[i, , drop = FALSE], na.rm = TRUE), "-")
    else M[i, ] <- 0 }
  M
}
# Spearman correlation of every column of M with y, vectorised through ranks.
.cor_cols <- function(M, y, method = c("spearman", "pearson")) {
  method <- match.arg(method)
  if (method == "spearman") { M <- apply(M, 2, rank); y <- rank(y) }
  Mc <- sweep(M, 2, colMeans(M), "-"); yc <- y - mean(y)
  sd_m <- sqrt(colSums(Mc^2)); sd_y <- sqrt(sum(yc^2))
  out <- as.numeric(crossprod(Mc, yc)) / (sd_m * sd_y)
  out[!is.finite(out)] <- NA_real_
  out
}
.fisher <- function(r) { r <- pmin(pmax(r, -0.999), 0.999); 0.5 * log((1 + r) / (1 - r)) }

# ── assemble one district-level frame per cell ──────────────────────────────
cfgs <- get_country_configs()
cells <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- H[H$country == cn, , drop = FALSE]
  if (!nrow(hc)) next
  ocs <- names(cc$outcomes)
  if (PROFILE == "smoke") ocs <- intersect(ocs, "child_iron")
  for (ocn in ocs) {
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(sv) || nrow(sv) < 10) next
    m <- dplyr::inner_join(
      sv[, intersect(c("Admin1", "Admin2", "svy_prev", "n_svy"), names(sv)), drop = FALSE],
      hc, by = admin2_join_by(sv, hc))
    m <- m[is.finite(m$svy_prev), , drop = FALSE]
    if (nrow(m) < 10 || dplyr::n_distinct(m$Admin1) < 2) next
    X <- as.matrix(m[, COVS, drop = FALSE])
    # The wide scan required a column to be complete in a cell. That rule is too
    # strict for this family: the nutrition-proximal columns come from the prior
    # DHS round through a name crosswalk that leaves 6 of Ghana's 75 districts
    # unmatched, which dropped the family to FOUR usable columns for Ghana on the
    # first run. A column is admitted here at 80 percent completeness and the
    # remainder is median-imputed within the cell, with the rate recorded.
    frac_ok <- colMeans(is.finite(X))
    keep <- frac_ok >= 0.8 & apply(X, 2, function(z) {
      sdv <- stats::sd(z, na.rm = TRUE); is.finite(sdv) && sdv > 0 })
    keep[is.na(keep)] <- FALSE
    for (j in which(keep)) {
      md <- stats::median(X[, j], na.rm = TRUE)
      X[!is.finite(X[, j]), j] <- if (is.finite(md)) md else 0
    }
    imp_rate <- if (any(keep)) round(1 - mean(frac_ok[keep]), 3) else NA_real_
    rmx <- REL$r_max_emp[match(paste(kk(cc$country), ocn),
                               paste(kk(REL$country), REL$outcome))]
    cells[[paste(cc$country, ocn)]] <- list(
      country = cc$country, outcome = ocn, n = nrow(m),
      y = .logit(m$svy_prev), region = m$Admin1, X = X, keep = keep,
      r_max = if (length(rmx) && is.finite(rmx)) rmx else NA_real_,
      shared = ocn %in% SHARED)
    cat(sprintf("  cell %-13s %-13s districts=%3d usable=%2d imputed=%-6s r_max=%s\n",
                cc$country, ocn, nrow(m), sum(keep),
                ifelse(is.finite(imp_rate), sprintf("%.1f%%", 100 * imp_rate), "NA"),
                ifelse(is.finite(rmx), sprintf("%.3f", rmx), "NA")))
  }
}
if (!length(cells)) stop("No cells assembled.")

# ── per-cell associations, both families ────────────────────────────────────
assoc <- function(cell, family, y = NULL) {
  X <- cell$X; yy <- if (is.null(y)) cell$y else y
  if (family == "region_partialed") {
    X <- .resid_on(X, cell$region); yy <- as.numeric(.resid_on(matrix(yy, ncol = 1), cell$region))
  }
  # A column carrying any NA returns NA from sd(), and a column that is constant
  # WITHIN every region becomes exactly zero after residualisation. Both must be
  # excluded before subsetting, and `keep` (complete cases in this cell) is
  # applied first so the two masks cannot disagree.
  sdv  <- suppressWarnings(apply(X, 2, stats::sd))
  sdok <- cell$keep & is.finite(sdv) & sdv > 1e-12
  sdok[is.na(sdok)] <- FALSE
  r <- rep(NA_real_, ncol(X))
  if (any(sdok) && is.finite(stats::sd(yy)) && stats::sd(yy) > 1e-12)
    r[sdok] <- .cor_cols(X[, sdok, drop = FALSE], yy, "spearman")
  r
}
FAM <- c("marginal", "region_partialed")
obs <- list()
for (f in FAM) for (nm in names(cells))
  obs[[f]][[nm]] <- assoc(cells[[nm]], f)

# ── the meta statistic: within-country mean z, then random effects over 4 ────
meta_stat <- function(rlist, cellset, disatt = TRUE) {
  ctys <- unique(vapply(cellset, function(z) z$country, character(1)))
  P <- length(COVS)
  zc <- matrix(NA_real_, P, length(ctys), dimnames = list(COVS, ctys))
  wc <- matrix(0,        P, length(ctys), dimnames = list(COVS, ctys))
  for (nm in names(cellset)) {
    cl <- cellset[[nm]]; r <- rlist[[nm]]
    if (is.null(r)) next
    rr <- r
    if (disatt && is.finite(cl$r_max) && cl$r_max > 0.05)
      rr <- pmin(pmax(r / cl$r_max, -0.99), 0.99)
    z <- .fisher(rr); w <- max(cl$n - 3, 1)
    j <- match(cl$country, ctys)
    ok <- is.finite(z)
    zc[ok, j] <- ifelse(is.na(zc[ok, j]), 0, zc[ok, j]) + z[ok] * w
    wc[ok, j] <- wc[ok, j] + w
  }
  zbar <- zc / ifelse(wc > 0, wc, NA)          # weighted mean z within country
  vbar <- 1 / ifelse(wc > 0, wc, NA)           # variance of that mean
  # DerSimonian-Laird across countries
  out <- matrix(NA_real_, P, 3, dimnames = list(COVS, c("z", "tau", "k")))
  for (i in seq_len(P)) {
    zi <- zbar[i, ]; vi <- vbar[i, ]
    ok <- is.finite(zi) & is.finite(vi)
    k <- sum(ok); if (k < 2) next
    zi <- zi[ok]; vi <- vi[ok]
    wf <- 1 / vi; zf <- sum(wf * zi) / sum(wf)
    Q <- sum(wf * (zi - zf)^2)
    C <- sum(wf) - sum(wf^2) / sum(wf)
    tau2 <- max(0, (Q - (k - 1)) / C)
    wr <- 1 / (vi + tau2); zr <- sum(wr * zi) / sum(wr)
    out[i, ] <- c(zr / sqrt(1 / sum(wr)), sqrt(tau2), k)   # z as a test statistic
  }
  list(stat = out, zbar = zbar)
}

cellsets <- list(shared = cells[vapply(cells, function(z) z$shared, logical(1))],
                 all    = cells)
res <- list()
for (f in FAM) for (cs in names(cellsets)) {
  st <- meta_stat(obs[[f]], cellsets[[cs]], disatt = TRUE)
  stn <- meta_stat(obs[[f]], cellsets[[cs]], disatt = FALSE)
  sign_ct <- rowSums(vapply(names(cellsets[[cs]]),
    function(nm) sign(obs[[f]][[nm]]), numeric(length(COVS))) > 0, na.rm = TRUE)
  n_cells <- rowSums(vapply(names(cellsets[[cs]]),
    function(nm) is.finite(obs[[f]][[nm]]), numeric(length(COVS))))
  res[[paste(f, cs)]] <- data.frame(
    predictor = COVS, family = f, cellset = cs,
    z = round(st$stat[, "z"], 4), tau = round(st$stat[, "tau"], 4),
    k_countries = st$stat[, "k"],
    z_undisattenuated = round(stn$stat[, "z"], 4),
    n_cells = n_cells, n_positive = sign_ct,
    p_analytic = 2 * stats::pnorm(-abs(st$stat[, "z"])),
    stringsAsFactors = FALSE)
  for (j in seq_len(ncol(st$zbar)))
    res[[paste(f, cs)]][[paste0("z_", kk(colnames(st$zbar)[j]))]] <- round(st$zbar[, j], 4)
}

# ── permutation calibration ─────────────────────────────────────────────────
cat(sprintf("\n[perm] %d permutations, outcomes shuffled within region within country\n", NPERM))
perm_max <- list()
for (f in FAM) for (cs in names(cellsets)) {
  key <- paste(f, cs); cset <- cellsets[[cs]]
  nullz <- matrix(NA_real_, NPERM, length(COVS))
  for (b in seq_len(NPERM)) {
    rl <- list()
    for (nm in names(cset)) {
      cl <- cset[[nm]]
      yp <- cl$y
      for (lv in unique(cl$region)) { i <- which(cl$region == lv)
        if (length(i) > 1) yp[i] <- sample(yp[i]) }
      rl[[nm]] <- assoc(cl, f, y = yp)
    }
    nullz[b, ] <- meta_stat(rl, cset, disatt = TRUE)$stat[, "z"]
    if (b %% 100 == 0) cat(sprintf("   %s %d/%d\n", key, b, NPERM))
  }
  # Per-predictor two-sided permutation p, and the max-|z| null for a
  # family-wise statement that needs no independence assumption.
  d <- res[[key]]
  d$p_perm <- vapply(seq_along(COVS), function(i) {
    nz <- nullz[, i]; nz <- nz[is.finite(nz)]
    if (!length(nz) || !is.finite(d$z[i])) return(NA_real_)
    (1 + sum(abs(nz) >= abs(d$z[i]))) / (1 + length(nz)) }, numeric(1))
  mx <- apply(abs(nullz), 1, function(v) suppressWarnings(max(v, na.rm = TRUE)))
  mx <- mx[is.finite(mx)]
  d$p_perm_fwer <- vapply(seq_along(COVS), function(i)
    if (!is.finite(d$z[i])) NA_real_ else (1 + sum(mx >= abs(d$z[i]))) / (1 + length(mx)),
    numeric(1))
  d$q_analytic <- stats::p.adjust(d$p_analytic, "BH")
  d$q_perm     <- stats::p.adjust(d$p_perm, "BH")
  d$n_permutations <- NPERM
  res[[key]] <- d
  perm_max[[key]] <- mx
}

out <- dplyr::bind_rows(res)
readr::write_csv(out, file.path(TDIR, sprintf("nutrition_proximal_confirmatory%s.csv", SUF)))

cat("\n=== WS-A: survivors by family and cell set ===\n")
s <- out |> group_by(family, cellset) |>
  summarise(predictors = dplyr::n(),
            testable = sum(is.finite(z)),
            q_analytic_lt_05 = sum(q_analytic < 0.05, na.rm = TRUE),
            q_perm_lt_05 = sum(q_perm < 0.05, na.rm = TRUE),
            fwer_lt_05 = sum(p_perm_fwer < 0.05, na.rm = TRUE),
            max_abs_z = round(max(abs(z), na.rm = TRUE), 3),
            median_tau = round(stats::median(tau, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(s), row.names = FALSE)

top <- out |> filter(cellset == "shared") |> group_by(family) |>
  arrange(desc(abs(z))) |> slice_head(n = 20) |> ungroup()
readr::write_csv(top, file.path(TDIR, sprintf("nutrition_proximal_confirmatory_top%s.csv", SUF)))
cat("\n=== top 10 per family, shared-outcome cells ===\n")
print(as.data.frame(top |> group_by(family) |> slice_head(n = 10) |>
  select(family, predictor, z, tau, k_countries, n_positive, n_cells,
         q_analytic, q_perm, p_perm_fwer)), row.names = FALSE)

png(file.path(FDIR, sprintf("nutrition_proximal_confirmatory%s.png", SUF)),
    width = 1200, height = 800, res = 130)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (f in FAM) {
  d <- out[out$family == f & out$cellset == "shared", ]
  z <- d$z[is.finite(d$z)]
  hist(z, breaks = 30, col = "grey88", border = "grey40",
       main = sub("_", " ", f), xlab = "meta z")
  abline(v = c(-1.96, 1.96), lty = 2, col = "grey50")
  mxk <- perm_max[[paste(f, "shared")]]
  if (length(mxk)) abline(v = c(-1, 1) * stats::quantile(mxk, 0.95), col = "firebrick", lwd = 2)
}
par(op); dev.off()
cat(sprintf("\n-> %s\n", file.path("results", "figures",
                                   sprintf("nutrition_proximal_confirmatory%s.png", SUF))))
