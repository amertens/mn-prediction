# =============================================================================
# scripts/signal_probes/p6_within_country_omics_scan.R
#
# P6. WITHIN-COUNTRY signal, on the PRE-HARMONIZATION (native) vocabulary,
#     with -omics-style machinery.
#
# MOTIVATION (PI request)
# -----------------------
# The harmonization step cuts the native GEE vocabulary (~296 columns per
# country) down to the 67 GEE columns shared across countries; every pooled
# scan so far tested only the shared 373. This scan asks, per country and per
# cell, the low-N-high-P question the way genomics asks it:
#   (a) mass bivariate tests with BH correction and a Storey pi0 estimate of
#       the fraction of truly-null predictors,
#   (b) an omnibus dense-signal statistic (mean r^2, permutation-calibrated),
#   (c) cross-validated penalized regression (elastic net) under BOTH fold
#       schemes — random folds (the interpolation estimand) and
#       leave-one-Admin1-out folds (the extrapolation estimand the pipeline
#       uses) — so the fold-scheme cost is measured on identical data,
#   (d) the highly adaptive lasso (hal9001) on a fold-internal prescreen,
#       random folds, core cells only (it is slow).
#
# Decision rule stated by the PI: if NO predictor shows within-country signal
# anywhere, suspect data linkage/cleaning. (Spoiler from P5: the omnibus
# already rejects that for Gambia, Ghana, Malawi.)
#
# PREDICTOR SET per country = native gee_admin2 columns (pre-harmonization)
#   UNION the country's rows of the shared 373 set (which carries the
#   non-GEE domains: DHS prior-round, SoilGrids/iSDA, MAP, MapSPAM, ...).
#   Joined to outcomes on the PAIR KEY (Admin1, Admin2). Rank-normalized
#   within country; columns with < 70% coverage on analysis rows dropped,
#   remainder median-imputed (outcome-independent).
#
# OUTCOMES per cell: binary svy_prev (logit for fitting, Spearman for
#   scoring) and the continuous district mean of -log(biomarker).
#
#   Rscript scripts/signal_probes/p6_within_country_omics_scan.R
# -> results/tables/signal_probes/p6_cells.csv        one row per cell x target
# -> results/tables/signal_probes/p6_bivariate_top.csv top hits with q-values
# -> results/tables/signal_probes/p6_model_skill.csv   glmnet/HAL by fold scheme
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(glmnet)})
setwd('C:/Users/andre/OneDrive/Documents/mn-prediction')

STORE  <- "_targets_full"
B      <- 2000L
SEED   <- 20260951L
OUTDIR <- "results/tables/signal_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
RUN_HAL <- !identical(Sys.getenv("P6_SKIP_HAL"), "1")

COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi",
               sierraleone = "SierraLeone")
CORE <- c("child_iron", "child_vitA", "women_iron", "women_vitA")

S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv",
               check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
SHARED_NONGEE <- MD$column[MD$source != "GEE"]
source("R/config.R")
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
rknorm <- function(x) {
  ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (sum(ok) > 2 && sd(x[ok]) > 0) out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
  out
}
.logit <- function(p, eps = 0.005) qlogis(pmin(pmax(p, eps), 1 - eps))

# ── native + shared predictor matrix per country, on covariate geography ────
native <- list()
for (lc in names(COUNTRIES)) {
  cv <- targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE)
  g <- cv$gee_admin2
  keys <- intersect(c("Admin1", "Admin2"), names(g))
  gee_cols <- setdiff(names(g), c(keys, "country",
                                  grep("^(id|ID|geometry)", names(g), value = TRUE)))
  gee_cols <- gee_cols[vapply(g[gee_cols], is.numeric, TRUE)]
  gtab <- g[, c(keys, gee_cols)]
  sh <- S[S$country == COUNTRIES[[lc]],
          c("Admin1", "Admin2", intersect(SHARED_NONGEE, names(S)))]
  by <- intersect(keys, c("Admin1", "Admin2"))
  m <- if (length(by) == 2) full_join(gtab, sh, by = by) else
       full_join(gtab, sh, by = "Admin2")
  native[[lc]] <- m
  cat(lc, ": native GEE", length(gee_cols), "cols + shared non-GEE",
      ncol(sh) - 2, "-> matrix", nrow(m), "x", ncol(m) - 2, "\n")
}

# ── outcome frames per cell (binary + continuous at Admin-2) ────────────────
cfgs <- get_country_configs()
cells <- list()
for (cn in names(cfgs)) {
  lc <- tolower(cn)
  if (!lc %in% names(COUNTRIES)) next
  cc <- cfgs[[cn]]
  for (on in names(cc$outcomes)) {
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on),
                                         store = STORE), error = function(e) NULL)
    if (is.null(sv)) next
    sv <- sv[is.finite(sv$svy_prev) & is.finite(sv$n_svy) & sv$n_svy > 0, ]
    if (sum(sv$svy_prev * sv$n_svy) < 10) next
    # continuous district means from individual data
    oc <- cc$outcomes[[on]]
    cont <- NULL
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on),
                                         store = STORE), error = function(e) NULL)
    if (!is.null(od) && all(c(oc$continuous, "Admin1", "Admin2") %in% names(od$data))) {
      d <- od$data
      y <- num(d[[oc$continuous]])
      w <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d))
             num(d[[cc$weight_col]]) else rep(1, nrow(d))
      w[!is.finite(w) | w <= 0] <- 1
      t <- if (identical(oc$cutoff_scale, "log")) y else
           { y[!is.finite(y) | y <= 0] <- NA; log(y) }
      keep <- is.finite(t)
      cont <- data.frame(Admin1 = as.character(d$Admin1)[keep],
                         Admin2 = as.character(d$Admin2)[keep],
                         t = t[keep], w = w[keep]) |>
        group_by(Admin1, Admin2) |>
        summarise(y_cont = -wmean(t, w), n_cont = dplyr::n(), .groups = "drop") |>
        filter(is.finite(y_cont), n_cont >= 5)
    }
    cells[[paste(lc, on, sep = "|")]] <- list(lc = lc, outcome = on,
                                              sv = sv, cont = cont)
  }
}
cat("cells:", length(cells), "\n")

# ── assemble per-cell X and y, run tests ────────────────────────────────────
spearman_block <- function(Xr, yr) {
  Xc <- sweep(Xr, 2, colMeans(Xr), "-"); yc <- yr - mean(yr)
  r <- as.numeric(crossprod(Xc, yc)) / (sqrt(colSums(Xc^2)) * sqrt(sum(yc^2)))
  r[!is.finite(r)] <- 0; r
}
storey_pi0 <- function(p) min(1, mean(p > 0.5, na.rm = TRUE) / 0.5)

cell_rows <- list(); biv_rows <- list(); skill_rows <- list()
for (kx in names(cells)) {
  cl <- cells[[kx]]; lc <- cl$lc; cn <- COUNTRIES[[lc]]
  Xtab <- native[[lc]]
  by <- intersect(c("Admin1", "Admin2"), intersect(names(cl$sv), names(Xtab)))
  base <- inner_join(cl$sv[, c("Admin1", "Admin2", "svy_prev", "n_svy")],
                     Xtab, by = by)
  if (!is.null(cl$cont))
    base <- left_join(base, cl$cont, by = c("Admin1", "Admin2"))
  pcols <- setdiff(names(base), c("Admin1", "Admin2", "svy_prev", "n_svy",
                                  "y_cont", "n_cont"))
  for (target in c("binary", "continuous")) {
    y_raw <- if (target == "binary") base$svy_prev else base$y_cont
    if (is.null(y_raw) || all(is.na(y_raw))) next
    rows_ok <- is.finite(y_raw)
    y <- y_raw[rows_ok]
    if (length(y) < 10 || sd(y) == 0) next
    X0 <- as.matrix(base[rows_ok, pcols, drop = FALSE])
    Xr <- apply(X0, 2, rknorm)
    cov_j <- colMeans(is.finite(Xr))
    keep <- cov_j >= 0.70 &
      apply(Xr, 2, function(z) sum(is.finite(z)) > 2 && sd(z[is.finite(z)]) > 0)
    Xr <- Xr[, keep, drop = FALSE]; Xr[!is.finite(Xr)] <- 0
    if (ncol(Xr) < 20) next
    cols <- pcols[keep]
    yr <- rank(y); n <- length(y)

    # (a)+(b): mass bivariate + omnibus under free permutation within country
    P <- replicate(B, sample(n))
    Ymat <- cbind(yr, matrix(yr[P], nrow = n))
    Ys <- sweep(Ymat, 2, colMeans(Ymat), "-")
    Ys <- sweep(Ys, 2, sqrt(colSums(Ys^2)), "/")
    Xs <- sweep(Xr, 2, colMeans(Xr), "-")
    sdx <- sqrt(colSums(Xs^2)); sdx[sdx <= 0] <- NA
    Xs <- sweep(Xs, 2, sdx, "/")
    R <- crossprod(Xs, Ys); R[!is.finite(R)] <- 0     # p x (1+B)
    r_obs <- R[, 1]
    p_biv <- (1 + rowSums(abs(R[, -1, drop = FALSE]) >= abs(r_obs))) / (B + 1)
    q_biv <- p.adjust(p_biv, "BH")
    pi0 <- storey_pi0(p_biv)
    Q <- colMeans(R^2)
    p_Q <- (1 + sum(Q[-1] >= Q[1])) / (B + 1)
    Tm <- apply(abs(R), 2, max)
    p_T <- (1 + sum(Tm[-1] >= Tm[1])) / (B + 1)

    cell_rows[[paste(kx, target)]] <- data.frame(
      country = cn, outcome = cl$outcome, target = target, n_areas = n,
      p_native = ncol(Xr), n_q10 = sum(q_biv < 0.10), n_q05 = sum(q_biv < 0.05),
      min_q = min(q_biv), pi0 = pi0, est_n_nonnull = round((1 - pi0) * ncol(Xr)),
      Q_obs = Q[1], Q_null = mean(Q[-1]),
      Q_z = (Q[1] - mean(Q[-1])) / sd(Q[-1]), p_omnibus_Q = p_Q, p_Tmax = p_T)

    topk <- order(p_biv)[seq_len(min(10, sum(q_biv < 0.10) + 3))]
    biv_rows[[paste(kx, target)]] <- data.frame(
      country = cn, outcome = cl$outcome, target = target,
      predictor = cols[topk], r = round(r_obs[topk], 3),
      p_perm = p_biv[topk], q_bh = round(q_biv[topk], 4))

    # (c) elastic net under two fold schemes
    yfit <- if (target == "binary") .logit(y) else y
    blk <- base$Admin1[rows_ok]
    for (scheme in c("random10", "loro")) {
      folds <- if (scheme == "random10") sample(rep(seq_len(min(10, n)), length.out = n))
               else as.integer(factor(blk))
      nf <- length(unique(folds))
      if (nf < 3) next
      oof <- rep(NA_real_, n)
      for (f in unique(folds)) {
        tr <- folds != f
        if (sum(tr) < 15) next
        cvf <- tryCatch(cv.glmnet(Xr[tr, , drop = FALSE], yfit[tr], alpha = 0.5,
                                  nfolds = 5, standardize = FALSE),
                        error = function(e) NULL)
        if (is.null(cvf)) next
        oof[!tr] <- as.numeric(predict(cvf, Xr[!tr, , drop = FALSE],
                                       s = "lambda.min"))
      }
      ok <- is.finite(oof)
      sk <- if (sum(ok) > 5 && sd(oof[ok]) > 0)
              cor(oof[ok], y[ok], method = "spearman") else NA_real_
      skill_rows[[paste(kx, target, "enet", scheme)]] <- data.frame(
        country = cn, outcome = cl$outcome, target = target, model = "enet",
        scheme = scheme, n_areas = n, oof_spearman = sk)
    }

    # (d) HAL, random folds, core cells, binary target only (runtime)
    if (RUN_HAL && target == "binary" && cl$outcome %in% CORE) {
      folds <- sample(rep(seq_len(5), length.out = n))
      oof <- rep(NA_real_, n)
      for (f in unique(folds)) {
        tr <- folds != f
        rtr <- spearman_block(apply(Xr[tr, , drop = FALSE], 2, rank), rank(yfit[tr]))
        top <- order(-abs(rtr))[seq_len(min(15, ncol(Xr)))]
        fit <- tryCatch(suppressWarnings(hal9001::fit_hal(
                 X = Xr[tr, top, drop = FALSE], Y = yfit[tr],
                 max_degree = 2, smoothness_orders = 1,
                 num_knots = c(10, 3), family = "gaussian",
                 fit_control = list(cv_select = TRUE, nfolds = 5))),
               error = function(e) NULL)
        if (is.null(fit)) next
        oof[!tr] <- tryCatch(predict(fit, new_data = Xr[!tr, top, drop = FALSE]),
                             error = function(e) NA_real_)
      }
      ok <- is.finite(oof)
      sk <- if (sum(ok) > 5 && sd(oof[ok]) > 0)
              cor(oof[ok], y[ok], method = "spearman") else NA_real_
      skill_rows[[paste(kx, "hal")]] <- data.frame(
        country = cn, outcome = cl$outcome, target = "binary", model = "hal",
        scheme = "random5", n_areas = n, oof_spearman = sk)
      cat("  HAL", kx, ":", round(sk, 3), "\n")
    }
  }
  cat("done", kx, "\n")
}

cellT <- bind_rows(cell_rows); bivT <- bind_rows(biv_rows)
skillT <- bind_rows(skill_rows)
write.csv(cellT, file.path(OUTDIR, "p6_cells.csv"), row.names = FALSE)
write.csv(bivT, file.path(OUTDIR, "p6_bivariate_top.csv"), row.names = FALSE)
write.csv(skillT, file.path(OUTDIR, "p6_model_skill.csv"), row.names = FALSE)

cat("\n=== per-cell summary (binary) ===\n")
print(cellT |> filter(target == "binary") |>
  select(country, outcome, n_areas, p_native, n_q10, min_q, pi0,
         est_n_nonnull, Q_z, p_omnibus_Q) |> as.data.frame(), row.names = FALSE)
cat("\n=== per-cell summary (continuous) ===\n")
print(cellT |> filter(target == "continuous") |>
  select(country, outcome, n_areas, p_native, n_q10, min_q, pi0,
         est_n_nonnull, Q_z, p_omnibus_Q) |> as.data.frame(), row.names = FALSE)
cat("\n=== model skill by fold scheme ===\n")
print(skillT |> group_by(model, scheme, target) |>
  summarise(cells = n(), mean_sp = round(mean(oof_spearman, na.rm = TRUE), 3),
            n_pos = sum(oof_spearman > 0, na.rm = TRUE), .groups = "drop") |>
  as.data.frame(), row.names = FALSE)
cat("\nDONE\n")
