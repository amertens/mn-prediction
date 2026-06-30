# =============================================================================
# R_corrected/p2_p6_methods.R — corrected analogues for P2,P3,P4,P5,P6
# All consume the P1 fit object (honest OOF) and production survey aggregates.
# Joins are keyed on (country, Admin2) throughout (P8).
# =============================================================================

# ── P2: out-of-fold calibration (Platt), honest vs in-sample ────────────────
# Mirrors R/benchmark_models.R::calibrate_binary_oof (fold-aware) and contrasts
# it with the in-sample calibration the production diagnostics actually use.
.platt_fit <- function(p, y) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  fit <- suppressWarnings(stats::glm(y ~ stats::qlogis(p),
                                     family = stats::binomial()))
  function(pnew) {
    pnew <- pmin(pmax(pnew, 1e-6), 1 - 1e-6)
    as.numeric(stats::predict(fit, newdata = data.frame(p = pnew),
                              type = "response"))
  }
}
.calib_slope <- function(p, y) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  fit <- tryCatch(suppressWarnings(stats::glm(y ~ stats::qlogis(p),
                                              family = stats::binomial())),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(intercept = NA, slope = NA))
  co <- stats::coef(fit); c(intercept = unname(co[1]), slope = unname(co[2]))
}

diagnostics_calibrated_oof <- function(slfit) {
  if (!is.null(slfit$error)) return(NULL)
  d <- slfit$honest_cluster$oof
  y <- d$Y; p <- d$yhat_full; fold <- d$fold
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]; fold <- fold[ok]

  # in-sample Platt (fit & eval on all — the production anti-pattern)
  fn_in <- .platt_fit(p, y); p_in <- fn_in(p)
  # OOF Platt: fit on out-of-fold, apply in-fold (honest)
  p_oof <- rep(NA_real_, length(p))
  for (f in unique(fold)) {
    inf <- which(fold == f); ouf <- which(fold != f)
    if (length(unique(y[ouf])) < 2) { p_oof[inf] <- mean(y[ouf]); next }
    fn <- .platt_fit(p[ouf], y[ouf]); p_oof[inf] <- fn(p[inf])
  }
  mk <- function(tag, pp) {
    cs <- .calib_slope(pp, y)
    data.frame(country = slfit$country, outcome = slfit$outcome,
               calibration = tag,
               brier = round(mean((pp - y)^2), 4),
               calib_intercept = round(cs["intercept"], 3),
               calib_slope = round(cs["slope"], 3), row.names = NULL)
  }
  rbind(mk("raw (uncalibrated)", p),
        mk("in-sample Platt (optimistic)", p_in),
        mk("out-of-fold Platt (honest)", p_oof))
}

# ── helper: corrected Admin-2 predictions from honest OOF ────────────────────
# Aggregate individual honest OOF predicted probabilities to district means and
# join survey-observed prevalence on (country, Admin2).
corrected_admin2 <- function(slfit, svy_admin2) {
  d <- slfit$honest_cluster$oof
  d <- d[!is.na(d$yhat_full) & !is.na(d$Admin2), , drop = FALSE]
  # (Issue 2) survey-WEIGHTED district mean of the honest OOF predictions, so the
  # predicted prevalence shares a weighting basis with the design-weighted
  # `svy_prev` truth. `w` is threaded in by fit_corrected_sl()/attach_meta();
  # if absent (older slfit objects) fall back to an unweighted mean.
  if (is.null(d$w)) d$w <- 1
  d$w[!is.finite(d$w) | d$w <= 0] <- 1
  agg <- do.call(rbind, lapply(split(seq_len(nrow(d)), d$Admin2), function(ix)
    data.frame(Admin2    = d$Admin2[ix][1],
               pred_prev = stats::weighted.mean(d$yhat_full[ix], d$w[ix], na.rm = TRUE),
               stringsAsFactors = FALSE)))
  rownames(agg) <- NULL
  agg$country <- slfit$country
  sv <- svy_admin2
  if (is.null(sv) || !"Admin2" %in% colnames(sv)) return(agg)
  sv$country <- slfit$country               # single-country slice
  keep <- intersect(c("country", "Admin2", "svy_prev", "svy_prev_se",
                      "svy_prev_low", "svy_prev_upp", "n_svy"), colnames(sv))
  merge(agg, sv[, keep, drop = FALSE], by = c("country", "Admin2"),
        all.x = TRUE, sort = FALSE)          # (P8) keyed on (country, Admin2)
}

# ── P5: sampling-error-aware Admin-2 error metrics ──────────────────────────
admin2_error_corrected <- function(slfit, svy_admin2) {
  if (!is.null(slfit$error)) return(NULL)
  m <- corrected_admin2(slfit, svy_admin2)
  m <- m[is.finite(m$pred_prev) & is.finite(m$svy_prev), , drop = FALSE]
  if (nrow(m) < 3) return(NULL)
  resid2 <- (m$pred_prev - m$svy_prev)^2
  # sampling variance of the "truth": prefer reported SE, else p(1-p)/n
  sv <- m$svy_prev_se^2
  if (all(is.na(sv)) && "n_svy" %in% colnames(m))
    sv <- m$svy_prev * (1 - m$svy_prev) / pmax(m$n_svy, 1)
  naive_mse <- mean(resid2)
  samp_mse  <- mean(sv, na.rm = TRUE)
  # (Issue 5) Truncating at 0 makes adjusted_rmse_pp a CONSERVATIVE (upward-biased)
  # point estimate of the noise-removed error: when the true adjusted MSE is near
  # zero the max(0, .) floor cannot return the (occasionally negative) unbiased
  # difference, so the reported value is biased toward positive. Footnote in any
  # table that surfaces adjusted_rmse_pp.
  adj_mse   <- max(0, naive_mse - samp_mse)   # remove truth's sampling noise
  # (Issue 6) bootstrap CI + permutation null on the district correlation; the
  # resampling unit is the Admin-2 area (rows of m).
  cin <- if (exists("metric_ci_null"))
    metric_ci_null(m$pred_prev, m$svy_prev, "pearson", seed = 202L)
  else data.frame(pearson_ci_lo = NA_real_, pearson_ci_hi = NA_real_,
                  pearson_perm_p = NA_real_, n_boot = 0L)
  data.frame(
    country = slfit$country, outcome = slfit$outcome, n_admin2 = nrow(m),
    naive_rmse_pp = round(sqrt(naive_mse) * 100, 2),
    sampling_rmse_pp = round(sqrt(samp_mse) * 100, 2),
    adjusted_rmse_pp = round(sqrt(adj_mse) * 100, 2),
    naive_mae_pp = round(mean(abs(m$pred_prev - m$svy_prev)) * 100, 2),
    pearson_r = round(suppressWarnings(stats::cor(m$pred_prev, m$svy_prev)), 3),
    pearson_ci_lo = cin$pearson_ci_lo, pearson_ci_hi = cin$pearson_ci_hi,
    pearson_perm_p = cin$pearson_perm_p, n_boot = cin$n_boot
  )
}

# ── P3: decision-value targets (targeting accuracy / lift) ──────────────────
decision_value_corrected <- function(slfit, svy_admin2) {
  if (!is.null(slfit$error)) return(NULL)
  m <- corrected_admin2(slfit, svy_admin2)
  m <- m[is.finite(m$pred_prev) & is.finite(m$svy_prev), , drop = FALSE]
  if (nrow(m) < 3) return(NULL)
  w <- if ("n_svy" %in% colnames(m)) m$n_svy else rep(1, nrow(m))
  w[!is.finite(w) | w <= 0] <- stats::median(w[w > 0], na.rm = TRUE) %||% 1
  burden <- m$svy_prev * w; totB <- sum(burden); totW <- sum(w)
  om <- order(m$pred_prev, decreasing = TRUE)
  cx <- c(0, cumsum(w[om]) / totW); cy <- c(0, cumsum(burden[om]) / totB)
  reach20 <- stats::approx(cx, cy, xout = 0.20, ties = "ordered", rule = 2)$y
  k <- max(1L, round(0.20 * nrow(m)))
  hit20 <- length(intersect(order(m$pred_prev, decreasing = TRUE)[seq_len(k)],
                            order(m$svy_prev,  decreasing = TRUE)[seq_len(k)])) / k
  overall <- stats::weighted.mean(m$svy_prev, w)
  sel_prev <- stats::weighted.mean(m$svy_prev[om[seq_len(k)]], w[om[seq_len(k)]])
  data.frame(
    country = slfit$country, outcome = slfit$outcome, n_admin2 = nrow(m),
    reach_at_20pct = round(reach20, 3),
    lift_vs_no_targeting = round(if (overall > 0) sel_prev / overall else NA, 2),
    hit_rate_top20 = round(hit20, 3),
    note = "within-country corrected SL admin-2 predictions vs survey truth; weight = n_svy"
  )
}

# ── P4: split-conformal (held-out calibration split) + design-based CIs ─────
# District-level intervals; NO Admin-1 -> Admin-2 broadcast. Split-conformal
# reserves a calibration subset of districts (never used to set the model
# prediction) to estimate the interval half-width, then validates coverage on
# the remaining districts. Design-based CIs come straight from the survey
# (srvyr) estimates per district.
intervals_corrected <- function(slfit, svy_admin2, alpha = 0.10, seed = 7L) {
  if (!is.null(slfit$error)) return(NULL)
  m <- corrected_admin2(slfit, svy_admin2)
  m <- m[is.finite(m$pred_prev) & is.finite(m$svy_prev), , drop = FALSE]
  if (nrow(m) < 6) return(list(districts = m, note = "too few districts"))
  set.seed(seed)
  cal <- sample(seq_len(nrow(m)), size = floor(nrow(m) / 2))
  res <- abs(m$pred_prev[cal] - m$svy_prev[cal])           # calibration residuals
  # Split-conformal level with the finite-sample (n+1)/n correction, clamped to
  # 1: for a small calibration set (n+1)(1-alpha) can exceed n, in which case the
  # conformal quantile is the maximum residual (prob = 1) — the standard behaviour.
  qprob <- min(1, (1 - alpha) * (length(cal) + 1) / length(cal))
  qh  <- as.numeric(stats::quantile(res, probs = qprob, names = FALSE, type = 1))
  qh <- min(qh, 1)
  m$conf_pi_lo <- pmax(0, m$pred_prev - qh)
  m$conf_pi_hi <- pmin(1, m$pred_prev + qh)
  # design-based CI from survey (district-level; no broadcast)
  m$design_lo <- if ("svy_prev_low" %in% colnames(m)) m$svy_prev_low else NA_real_
  m$design_hi <- if ("svy_prev_upp" %in% colnames(m)) m$svy_prev_upp else NA_real_
  m$design_available <- !is.na(m$design_lo) & !is.na(m$design_hi)
  # empirical coverage of split-conformal on the held-out districts
  test <- setdiff(seq_len(nrow(m)), cal)
  cov <- mean(m$svy_prev[test] >= m$conf_pi_lo[test] &
              m$svy_prev[test] <= m$conf_pi_hi[test])
  list(
    districts = m[, c("country", "Admin2", "pred_prev", "svy_prev",
                      "conf_pi_lo", "conf_pi_hi", "design_lo", "design_hi",
                      "design_available")],
    summary = data.frame(
      country = slfit$country, outcome = slfit$outcome,
      target_coverage = 1 - alpha,
      split_conformal_halfwidth_pp = round(qh * 100, 2),
      empirical_coverage_heldout = round(cov, 3),
      n_districts = nrow(m),
      n_design_ci = sum(m$design_available))
  )
}

# ── P6: out-of-support / trust flags ────────────────────────────────────────
# Per-district reliability from (a) covariate out-of-support distance to the
# surveyed (training) districts, (b) survey sample size, (c) prediction-interval
# width. gee_admin2 supplies covariates for ALL polygons (surveyed + not).
trust_flags_corrected <- function(slfit, svy_admin2, gee_admin2,
                                  interval_obj = NULL) {
  if (!is.null(slfit$error)) return(NULL)
  if (is.null(gee_admin2) || !"Admin2" %in% colnames(gee_admin2)) return(NULL)
  gee_cols <- grep("^gee_", colnames(gee_admin2), value = TRUE)
  if (length(gee_cols) < 3) return(NULL)
  # use a compact, well-behaved covariate subset
  G <- gee_admin2[, c("Admin2", gee_cols), drop = FALSE]
  X <- as.matrix(G[, gee_cols, drop = FALSE])
  X[!is.finite(X)] <- NA
  keepc <- which(apply(X, 2, function(z) mean(is.na(z)) < 0.3 &&
                         stats::sd(z, na.rm = TRUE) > 0))
  X <- X[, keepc, drop = FALSE]
  if (ncol(X) > 30) {                                    # cap dims, take top-var
    v <- apply(X, 2, stats::var, na.rm = TRUE)
    X <- X[, order(v, decreasing = TRUE)[1:30], drop = FALSE]
  }
  ctr <- colMeans(X, na.rm = TRUE); sds <- apply(X, 2, stats::sd, na.rm = TRUE)
  sds[sds < 1e-8] <- 1
  Z <- sweep(sweep(X, 2, ctr, "-"), 2, sds, "/")
  for (j in seq_len(ncol(Z))) Z[is.na(Z[, j]), j] <- 0
  # training (surveyed) districts define in-support region
  surveyed <- if (!is.null(svy_admin2)) svy_admin2$Admin2 else character(0)
  is_train <- G$Admin2 %in% surveyed
  if (sum(is_train) < 3) is_train <- rep(TRUE, nrow(G))
  tr_ctr <- colMeans(Z[is_train, , drop = FALSE])
  dist <- sqrt(rowMeans(sweep(Z, 2, tr_ctr, "-")^2))     # mean standardized dist
  thr <- stats::quantile(dist[is_train], 0.95, names = FALSE)
  out <- data.frame(country = slfit$country, outcome = slfit$outcome,
                    Admin2 = G$Admin2,
                    support_distance = round(dist, 3),
                    out_of_support = dist > thr,
                    surveyed = G$Admin2 %in% surveyed,
                    stringsAsFactors = FALSE)
  if (!is.null(svy_admin2) && "n_svy" %in% colnames(svy_admin2)) {
    out <- merge(out, data.frame(Admin2 = svy_admin2$Admin2,
                                 n_svy = svy_admin2$n_svy),
                 by = "Admin2", all.x = TRUE, sort = FALSE)
    out$small_sample <- !is.na(out$n_svy) & out$n_svy < 15
  } else out$small_sample <- NA
  if (!is.null(interval_obj) && !is.null(interval_obj$districts)) {
    w <- interval_obj$districts
    w$pi_width <- w$conf_pi_hi - w$conf_pi_lo
    out <- merge(out, w[, c("Admin2", "pi_width")], by = "Admin2",
                 all.x = TRUE, sort = FALSE)
    out$wide_interval <- !is.na(out$pi_width) &
      out$pi_width > stats::median(out$pi_width, na.rm = TRUE) * 1.5
  } else out$wide_interval <- NA
  # traffic-light verdict
  red <- out$out_of_support | (isTRUE_vec(out$small_sample) & isTRUE_vec(out$wide_interval))
  amber <- !red & (isTRUE_vec(out$small_sample) | isTRUE_vec(out$wide_interval) |
                     out$support_distance > thr * 0.8)
  out$trust <- ifelse(red, "Do not rely",
                      ifelse(amber, "Use with caution", "OK to use"))
  out
}
isTRUE_vec <- function(x) !is.na(x) & x
