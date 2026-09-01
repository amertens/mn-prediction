# =============================================================================
# scripts/accuracy_impact/wso1_split_half_arena.R
#
# WS-O. Does shrinkage cost ranking ON THE REAL SURVEYS, and how do the
# empirical-Bayes estimators compare with covariate models and an ensemble?
#
# THE PROBLEM THIS SOLVES
# -----------------------
# On real data there is no known truth. Scoring an estimator against the
# survey's own district estimates gives the DIRECT estimate a correlation of
# exactly 1 by construction, so it can never be compared with anything. That is
# why the estimator tournament simulates: it needs a truth to score against.
# Everything the tournament concluded about shrinkage costing ranking is
# therefore a simulation result, and it has never been checked on the surveys.
#
# THE DESIGN
# ----------
# Split each district's respondents into halves A and B. Build every estimator
# from HALF A only; score against HALF B, which no estimator has seen.
#
# Half B is a noisy target, so every correlation is attenuated by roughly
# sqrt(reliability of a half sample). But that factor is the SAME for every
# estimator, because none of them saw B, so the COMPARISON between estimators is
# unbiased even though the absolute numbers are depressed. Absolute values here
# are lower bounds and should not be quoted against the tournament's.
#
# Splitting is done WITHIN CLUSTER so both halves inherit the same cluster
# structure; splitting by cluster would give the halves different designs and
# confound the comparison with a design difference.
#
# THE ARMS
#   direct        the district's own half-A estimate. No borrowing at all.
#   flat_region   its region's half-A mean, jackknifed.
#   eb_blend      empirical Bayes between the two. Shrinkage, NO covariates.
#   ridge         covariates, no shrinkage, fitted out of region.
#   rf            random forest on covariates, out of region.
#   spatial       a covariate-free smooth of district centroids, out of region.
#   ensemble      non-negative least squares stack of ridge, rf and spatial,
#                 weights fitted out of region. The SuperLearner-style arm.
#   eb_ensemble   the ensemble used as the shrinkage target.
#
#   Rscript scripts/accuracy_impact/wso1_split_half_arena.R
# -> results/tables/split_half_arena.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here); library(sf)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")
SEED <- 20260902L
SPLITS <- as.integer(Sys.getenv("WSO1_SPLITS", "20"))
MIN_AREAS <- 12L; HOLDOUT_FRAC <- 0.30; TOPFRAC <- 0.20

P <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_shared.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(P), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# --- out-of-region prediction helpers --------------------------------------
oof_generic <- function(X, y, reg, fitfun) {
  ok <- is.finite(y); regs <- unique(reg[ok]); if (length(regs) < 3) return(NULL)
  k <- max(2L, round(1 / HOLDOUT_FRAC))
  fold <- split(sample(regs), rep(seq_len(k), length.out = length(regs)))
  out <- rep(NA_real_, length(y))
  for (f in seq_len(k)) {
    te <- which(reg %in% fold[[f]] & ok); tr <- which(!(reg %in% fold[[f]]) & ok)
    if (!length(te) || length(tr) < 8) next
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE); if (!is.finite(mu)) mu <- 0
      Xtr[!is.finite(Xtr[, j]), j] <- mu; Xte[!is.finite(Xte[, j]), j] <- mu
    }
    # The floor is 3 for a covariate matrix and 2 for the spatial arm, whose
    # design matrix is just (lon, lat). A flat floor of 3 skipped every fold of
    # the spatial arm, returning all-NA, which then silently emptied the
    # ensemble that tried to stack it.
    keep <- apply(Xtr, 2, function(z) stats::sd(z) > 0)
    if (sum(keep) < min(3L, ncol(Xtr))) next
    pp <- tryCatch(fitfun(Xtr[, keep, drop = FALSE], y[tr], Xte[, keep, drop = FALSE]),
                   error = function(e) NULL)
    if (!is.null(pp) && length(pp) == length(te)) out[te] <- pp
  }
  out
}
fit_ridge <- function(Xtr, ytr, Xte) {
  f <- .ds_fit(Xtr, ytr, k_screen = min(20L, ncol(Xtr)), link = "logit")
  if (is.null(f)) return(NULL); .ds_predict(f, Xte)
}
fit_rf <- function(Xtr, ytr, Xte) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  # screen first: 370 predictors on ~40 training areas is not a forest, it is
  # a lottery. Same in-fold screen the ridge uses.
  r <- abs(suppressWarnings(stats::cor(Xtr, ytr))); r[!is.finite(r)] <- 0
  sel <- colnames(Xtr)[order(r, decreasing = TRUE)][seq_len(min(20L, ncol(Xtr)))]
  m <- ranger::ranger(y = ytr, x = as.data.frame(Xtr[, sel, drop = FALSE]),
                      num.trees = 500, min.node.size = 3, seed = SEED)
  as.numeric(stats::predict(m, as.data.frame(Xte[, sel, drop = FALSE]))$predictions)
}

run_cell <- function(d, cc, oc, hc) {
  need <- c(cc$admin1_col, cc$admin2_col, cc$cluster_id, cc$weight_col, oc$binary)
  if (!all(need %in% names(d))) return(NULL)
  y <- num(d[[oc$binary]]); w <- num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- 1
  a1 <- trimws(as.character(d[[cc$admin1_col]])); a2 <- trimws(as.character(d[[cc$admin2_col]]))
  cl <- as.character(d[[cc$cluster_id]])
  ok <- is.finite(y) & !is.na(a1) & !is.na(a2) & !is.na(cl)
  y <- y[ok]; w <- w[ok]; a1 <- a1[ok]; a2 <- a2[ok]; cl <- cl[ok]
  key <- paste(a1, a2, sep = "||"); ukey <- sort(unique(key))
  if (length(ukey) < MIN_AREAS) return(NULL)
  spine <- data.frame(Admin1 = sub("\\|\\|.*$", "", ukey),
                      Admin2 = sub("^.*\\|\\|", "", ukey), stringsAsFactors = FALSE)
  m <- dplyr::inner_join(spine, hc, by = admin2_join_by(spine, hc))
  if (nrow(m) != nrow(spine)) return(NULL)
  covs <- COVS[vapply(COVS, function(v) mean(is.finite(m[[v]])) > 0.5 &&
                        stats::sd(m[[v]], na.rm = TRUE) > 0, logical(1))]
  if (length(covs) < 20) return(NULL)
  X <- as.matrix(m[, covs, drop = FALSE])
  ctr <- tryCatch({
    cen <- sf::st_drop_geometry(load_admin2_centroids(cc$gadm_code))
    j <- dplyr::left_join(spine, cen[, intersect(c("Admin1","Admin2","lon","lat"), names(cen))],
                          by = admin2_join_by(spine, cen))
    j <- j[!duplicated(paste(j$Admin1, j$Admin2)), ]
    if (nrow(j) == nrow(spine)) as.matrix(j[, c("lon","lat")]) else NULL
  }, error = function(e) NULL)

  wprev <- function(idx) vapply(ukey, function(k) {
    i <- idx[key[idx] == k]
    if (!length(i)) NA_real_ else stats::weighted.mean(y[i], w[i])
  }, numeric(1))
  nof <- function(idx) vapply(ukey, function(k) sum(key[idx] == k), numeric(1))

  acc <- list()
  for (b in seq_len(SPLITS)) {
    set.seed(SEED + b)
    # split WITHIN cluster so both halves keep the same design
    half <- unlist(lapply(split(seq_along(y), cl), function(i)
      sample(rep(c(TRUE, FALSE), length.out = length(i)))), use.names = FALSE)
    ord <- unlist(split(seq_along(y), cl), use.names = FALSE)
    isA <- logical(length(y)); isA[ord] <- half
    A <- which(isA); B <- which(!isA)
    pA <- wprev(A); pB <- wprev(B); nA <- nof(A)
    fin <- is.finite(pA) & is.finite(pB) & nA >= 2
    if (sum(fin) < MIN_AREAS) next

    est <- list(direct = pA)
    fr <- rep(NA_real_, length(ukey))
    for (i in seq_along(ukey)) {
      j <- which(spine$Admin1 == spine$Admin1[i] & seq_along(ukey) != i & fin)
      if (length(j) >= 2) fr[i] <- stats::weighted.mean(pA[j], pmax(nA[j], 1))
    }
    est$flat_region <- fr
    v_d <- area_sampling_var_shrunk(pA, pmax(nA, 1), deff = 1.5)
    okb <- fin & is.finite(fr)
    tau2 <- if (sum(okb) >= 10) max(stats::var(pA[okb] - fr[okb]) - mean(v_d[okb]), 1e-8) else 1e-8
    lam <- tau2 / (tau2 + v_d); lam[!is.finite(lam)] <- 0
    est$eb_blend <- ifelse(is.finite(fr), lam * pA + (1 - lam) * fr, pA)

    est$ridge <- oof_generic(X, ifelse(fin, pA, NA_real_), spine$Admin1, fit_ridge)
    est$rf    <- oof_generic(X, ifelse(fin, pA, NA_real_), spine$Admin1, fit_rf)
    if (!is.null(ctr) && requireNamespace("mgcv", quietly = TRUE)) {
      est$spatial <- oof_generic(ctr, ifelse(fin, pA, NA_real_), spine$Admin1,
        function(Xtr, ytr, Xte) {
          dd <- as.data.frame(Xtr); names(dd) <- c("lon","lat")
          kk <- min(20L, max(5L, floor(nrow(dd)/4)))
          g <- mgcv::gam(stats::qlogis(pmin(pmax(ytr,.005),.995)) ~ s(lon,lat,k=kk),
                         data = dd, method = "REML")
          de <- as.data.frame(Xte); names(de) <- c("lon","lat")
          stats::plogis(as.numeric(stats::predict(g, newdata = de)))
        })
    }
    # NNLS ensemble of the covariate arms, weights fitted on half A only
    # An all-NA arm is not a candidate. It passes a NULL check and then poisons
    # every row of the stack.
    cand <- Filter(function(z) !is.null(z) && any(is.finite(z)),
                   est[c("ridge","rf","spatial")])
    if (b == 1L) cat(sprintf("      [arms] covariate candidates: %s
",
                             paste(names(cand), collapse = ", ")))
    if (length(cand) >= 2 && requireNamespace("nnls", quietly = TRUE)) {
      M <- do.call(cbind, cand)
      o2 <- fin & apply(M, 1, function(z) all(is.finite(z)))
      if (sum(o2) >= 10) {
        fitw <- tryCatch(nnls::nnls(M[o2, , drop = FALSE], pA[o2]), error = function(e) NULL)
        if (!is.null(fitw)) {
          wt <- fitw$x; if (!any(wt > 0)) wt <- rep(1/ncol(M), ncol(M)); wt <- wt/sum(wt)
          est$ensemble <- as.numeric(M %*% wt)
          tg <- est$ensemble; okt <- fin & is.finite(tg)
          t2 <- if (sum(okt) >= 10) max(stats::var(pA[okt]-tg[okt]) - mean(v_d[okt]), 1e-8) else 1e-8
          l2 <- t2/(t2+v_d); l2[!is.finite(l2)] <- 1
          est$eb_ensemble <- ifelse(is.finite(tg), l2*pA + (1-l2)*tg, est$eb_blend)
        }
      }
    }
    for (nm in names(est)) {
      p <- est[[nm]]; if (is.null(p)) next
      o <- fin & is.finite(p); if (sum(o) < MIN_AREAS) next
      kk <- max(1L, round(TOPFRAC * sum(o)))
      sel <- order(p[o], decreasing = TRUE)[seq_len(kk)]
      bv <- pB[o]; ov <- mean(bv)
      acc[[nm]] <- rbind(acc[[nm]], data.frame(
        r = if (stats::sd(p[o]) > 0 && stats::sd(bv) > 0)
              suppressWarnings(stats::cor(p[o], bv)) else NA_real_,
        mae = 100 * mean(abs(p[o] - bv)),
        lift = if (ov > 0) mean(bv[sel]) / ov else NA_real_))
    }
  }
  if (!length(acc)) return(NULL)
  dplyr::bind_rows(lapply(names(acc), function(nm) {
    A <- acc[[nm]]
    data.frame(estimator = nm, splits = nrow(A),
               r = round(mean(A$r, na.rm = TRUE), 4),
               mae_pp = round(mean(A$mae, na.rm = TRUE), 3),
               lift = round(mean(A$lift, na.rm = TRUE), 4),
               stringsAsFactors = FALSE)
  }))
}

rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- P[P$country == cn, , drop = FALSE]; if (!nrow(hc)) next
  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]; if (is.null(oc$binary)) next
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od)) next
    r <- tryCatch(run_cell(od$data, cc, oc, hc), error = function(e) NULL)
    if (is.null(r)) next
    r$country <- cc$country; r$outcome <- ocn
    rows[[length(rows) + 1L]] <- r
    cat(sprintf("  [ok] %-13s %-13s %d arms\n", cn, ocn, nrow(r)))
  }
}
if (!length(rows)) stop("[wso1] nothing scored")
out <- dplyr::bind_rows(rows)
front <- c("country", "outcome", "estimator")
out <- out[, c(front, setdiff(names(out), front))]
readr::write_csv(out, file.path(TDIR, "split_half_arena.csv"))

cat("\n=== scored on a HELD-OUT HALF of the same survey, pooled over cells ===\n")
s <- out |> group_by(estimator) |>
  summarise(cells = dplyr::n(),
            r = round(mean(r, na.rm = TRUE), 3),
            mae_pp = round(mean(mae_pp, na.rm = TRUE), 2),
            lift = round(mean(lift, na.rm = TRUE), 3), .groups = "drop") |>
  arrange(desc(lift))
print(as.data.frame(s), row.names = FALSE)
cat("\nAbsolute values are attenuated by the noise in the held-out half and are\n")
cat("lower bounds; the COMPARISON between estimators is unbiased.\n")
cat(sprintf("-> %s\n", file.path("results","tables","split_half_arena.csv")))
