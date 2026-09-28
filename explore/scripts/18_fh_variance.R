# =============================================================================
# explore/scripts/18_fh_variance.R   [probe VF-01]
#
# HOW MUCH DOES THE FAY-HERRIOT SAMPLING-VARIANCE ASSUMPTION COST?
#
# The Bayesian SAE preprint (arXiv 2604.14971) does two things this project's
# Fay-Herriot does not: it MODELS the sampling variance instead of plugging one
# in, and it augments single-cluster areas with synthetic observations so a
# variance can be estimated at all.
#
# This project's FH (R/benchmark_models.R:229-270) instead:
#   v_i = p_i (1 - p_i) / (n_i / 1.5)        design effect FIXED at 1.5
#   sv  <- pmax(sv, 1e-8)                    single-cluster areas FLOORED
#
# Both assumptions are measurably wrong on this data:
#   - the district design effect (DE-01, already in targets_v2 as deff_binary)
#     has median 2.57, IQR 1.56-3.46, max 6.1 - not 1.5. Understating v_i makes
#     FH shrink too little, i.e. trust noisy direct estimates too much.
#   - 85% of Malawi, 83% of Ghana and 57% of Gambia district-rows are
#     SINGLE-CLUSTER (n_psu == 1). Flooring their variance at 1e-8 gives them
#     near-maximal weight when beta is estimated.
#
# WHAT IS ACTUALLY BEING TESTED. Under in-fill CV the held-out district's own
# direct estimate is hidden, so its prediction is the synthetic part x'beta.
# The sampling variances therefore act through the WEIGHTS on the training
# districts when beta is fitted. Understate them unevenly and beta is pulled
# toward the noisiest areas.
#
# ARMS (prevalence target, which is the scale FH is used on here)
#   fh_deff15     production: deff fixed at 1.5, variance floored at 1e-8
#   fh_deff_meas  each district's OWN measured design effect (n_eff column)
#   fh_varsmooth  the preprint's joint variance-mean smoothing:
#                 log v_i = g0 + g1 log(p(1-p)) + g2 log(n) fitted on the
#                 TRAINING districts, fitted values used as the variances
#   fh_phantom    the preprint's phantom-cluster idea: a single-cluster
#                 district does not get a floor, it borrows the pooled
#                 within-Admin1 variance of the multi-cluster districts
#   + domain_index and null_train_mean on identical folds
#
# FH is implemented here directly so that ONLY the variance differs between
# arms; fh_deff15 is checked against sae::eblupFH on one cell.
#
#   Rscript explore/scripts/18_fh_variance.R
# -> explore/out/18_fh_variance.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
TG <- E$TG
ix <- exp_cell_index(E)

# ── Fay-Herriot, variance supplied by the caller ────────────────────────────
#' @param y   direct estimates (training areas)
#' @param X   covariate matrix including intercept
#' @param D   sampling variances for the training areas
#' @return beta and the estimated between-area variance A
fh_fit <- function(y, X, D, maxit = 100L, tol = 1e-8) {
  n <- length(y); p <- ncol(X)
  A <- max(stats::var(y) - mean(D), 1e-6)     # moment start
  for (it in seq_len(maxit)) {
    w <- 1 / (A + D)
    XtW <- t(X * w)
    b <- tryCatch(solve(XtW %*% X, XtW %*% y), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    r <- as.numeric(y - X %*% b)
    # Fay-Herriot moment update for A
    A_new <- max((sum(w * (r^2 - D)) / sum(w)), 1e-8)
    if (abs(A_new - A) < tol * max(A, 1e-8)) { A <- A_new; break }
    A <- A_new
  }
  list(beta = b, A = A)
}

#' An arm factory: the variance rule is the only thing that changes.
#' `vrule(p, n_raw, n_eff, n_psu, admin1, tr)` returns variances for ALL rows.
make_fh_arm <- function(vrule, n_var_cap = 5L) {
  function(tr, te, y, X, D, aux) {
    p  <- pmin(pmax(aux$y_nat, 1e-4), 1 - 1e-4)
    v  <- vrule(p, aux$n_raw, aux$n_eff, aux$n_psu, aux$Admin1, tr)
    v[!is.finite(v) | v <= 0] <- stats::median(v[is.finite(v) & v > 0], na.rm = TRUE)
    # production restricts FH to the top-correlated covariates; do the same
    cc <- apply(D[tr, , drop = FALSE], 2, function(z)
      suppressWarnings(abs(stats::cor(z, y[tr]))))
    cc[!is.finite(cc)] <- 0
    keep <- order(cc, decreasing = TRUE)[seq_len(min(n_var_cap, ncol(D)))]
    Xd <- cbind(1, D[, keep, drop = FALSE])
    f <- fh_fit(y[tr], Xd[tr, , drop = FALSE], v[tr])
    if (is.null(f)) return(rep(mean(y[tr]), length(te)))
    as.numeric(Xd[te, , drop = FALSE] %*% f$beta)   # synthetic part: area unseen
  }
}

v_deff15 <- function(p, n_raw, n_eff, n_psu, a1, tr)
  pmax(p * (1 - p) / pmax(n_raw / 1.5, 1), 1e-8)
v_measured <- function(p, n_raw, n_eff, n_psu, a1, tr)
  p * (1 - p) / pmax(n_eff, 1)
v_varsmooth <- function(p, n_raw, n_eff, n_psu, a1, tr) {
  v_obs <- p * (1 - p) / pmax(n_eff, 1)
  d <- data.frame(lv = log(pmax(v_obs, 1e-10)),
                  lp = log(pmax(p * (1 - p), 1e-10)), ln = log(pmax(n_raw, 1)))
  fit <- tryCatch(stats::lm(lv ~ lp + ln, d[tr, , drop = FALSE]),
                  error = function(e) NULL)
  if (is.null(fit)) return(v_obs)
  exp(as.numeric(stats::predict(fit, newdata = d)))
}
v_phantom <- function(p, n_raw, n_eff, n_psu, a1, tr) {
  v <- p * (1 - p) / pmax(n_eff, 1)
  single <- !is.finite(n_psu) | n_psu < 2
  multi_tr <- tr[!single[tr]]
  if (!length(multi_tr)) return(v)
  # a single-cluster district borrows the variance of the multi-cluster
  # districts in its own Admin1, scaled to its own n; national pool as fallback
  for (i in which(single)) {
    same <- multi_tr[a1[multi_tr] == a1[i]]
    pool <- if (length(same) >= 2) same else multi_tr
    ratio <- stats::median(v[pool] * pmax(n_eff[pool], 1) /
                             pmax(p[pool] * (1 - p[pool]), 1e-10), na.rm = TRUE)
    v[i] <- p[i] * (1 - p[i]) * ratio / pmax(n_eff[i], 1)
  }
  v
}

ARMS <- c(exp_baseline_arms()[c("null_train_mean", "domain_index")],
          list(fh_deff15    = make_fh_arm(v_deff15),
               fh_deff_meas = make_fh_arm(v_measured),
               fh_varsmooth = make_fh_arm(v_varsmooth),
               fh_phantom   = make_fh_arm(v_phantom)))

# ── attach the design columns the variance rules need ───────────────────────
with_design <- function(cell) {
  t <- TG[TG$country == cell$country & TG$outcome == cell$outcome, ]
  i <- match(paste(cell$Admin1, cell$Admin2), paste(t$Admin1, t$Admin2))
  cell$aux$n_raw <- as.numeric(t$n_raw[i])
  cell$aux$n_eff <- as.numeric(t$n_eff[i])
  cell$aux$n_psu <- as.numeric(t$n_psu[i])
  cell$aux$Admin1 <- cell$Admin1
  cell
}

rows <- list()
for (i in seq_len(nrow(ix))) {
  cell <- tryCatch(with_design(exp_cell(E, ix$country[i], ix$outcome[i], "prev")),
                   error = function(e) NULL)
  if (is.null(cell) || !any(is.finite(cell$aux$n_eff))) next
  rows[[i]] <- exp_infill(cell, ARMS, reps = REPS)
  message("  ", ix$country[i], " ", ix$outcome[i],
          sprintf("  single-cluster %.0f%%", 100 * mean(cell$aux$n_psu < 2, na.rm = TRUE)))
}
SM <- exp_summarise(dplyr::bind_rows(rows)); exp_write(SM, "18_fh_variance")

cat("\n== in-fill, prevalence target: median Spearman over cells ==\n")
print(aggregate(spearman ~ arm, data = SM,
                FUN = function(z) round(median(z, na.rm = TRUE), 3)), row.names = FALSE)
cat("\n== and by MAE (points) ==\n")
print(aggregate(mae ~ arm, data = SM,
                FUN = function(z) round(median(z, na.rm = TRUE), 2)), row.names = FALSE)

cat("\n== paired against the production variance rule (fh_deff15) ==\n")
w <- reshape(SM[, c("country", "outcome", "arm", "spearman")],
             idvar = c("country", "outcome"), timevar = "arm", direction = "wide")
names(w) <- sub("^spearman[.]", "", names(w))
for (a in c("fh_deff_meas", "fh_varsmooth", "fh_phantom")) {
  g <- w[[a]] - w$fh_deff15; ok <- is.finite(g)
  b <- stats::aggregate(list(gain = g[ok]), by = list(b = paste(w$country[ok])), FUN = mean)
  cat(sprintf("  %-13s mean %+.4f  cells %2d/%2d  countries %d/%d  sign p=%.3f\n",
      a, mean(g[ok]), sum(g[ok] > 0), sum(ok), sum(b$gain > 0), nrow(b),
      stats::binom.test(sum(g[ok] > 0), sum(ok), 0.5)$p.value))
}
cat("\n== is FH competitive with the index at all? ==\n")
for (a in c("fh_deff15", "fh_deff_meas", "fh_varsmooth", "fh_phantom")) {
  g <- w[[a]] - w$domain_index; ok <- is.finite(g)
  cat(sprintf("  %-13s vs domain_index: mean %+.4f  cells %2d/%2d\n",
      a, mean(g[ok]), sum(g[ok] > 0), sum(ok)))
}
