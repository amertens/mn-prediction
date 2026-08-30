# =============================================================================
# sandbox_parsimony/R/03_core.R  --  shared fitting / CV / metric machinery
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(glmnet)})

.logit  <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.ilogit <- stats::plogis

#' Median-impute + drop-constant, using TRAIN statistics only
prep_X <- function(train, test, vars) {
  Xtr <- as.matrix(train[, vars, drop = FALSE])
  Xte <- as.matrix(test[,  vars, drop = FALSE])
  storage.mode(Xtr) <- "double"; storage.mode(Xte) <- "double"
  med <- apply(Xtr, 2, function(z) stats::median(z[is.finite(z)]))
  med[!is.finite(med)] <- 0
  for (j in seq_along(vars)) {
    Xtr[!is.finite(Xtr[, j]), j] <- med[j]
    Xte[!is.finite(Xte[, j]), j] <- med[j]
  }
  sdv <- apply(Xtr, 2, stats::sd); keep <- is.finite(sdv) & sdv > 1e-8
  list(Xtr = Xtr[, keep, drop = FALSE], Xte = Xte[, keep, drop = FALSE],
       vars = vars[keep])
}

# -- Model fitters. Each returns predicted PREVALENCE on [0,1] for `test`. ---
# `w` = observation weights (survey n); p_tr / n_tr = training prevalence and n.

fit_glmnet_logit <- function(Xtr, Xte, p_tr, n_tr, alpha = 0.5, lambda = "lambda.min",
                             w = NULL, seed = 1L) {
  if (ncol(Xtr) < 2) return(rep(mean(p_tr), nrow(Xte)))
  y <- .logit(p_tr)
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, y, alpha = alpha, weights = w,
                                   nfolds = min(5, length(y))), error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(p_tr), nrow(Xte)))
  pmin(pmax(.ilogit(as.numeric(stats::predict(cv, newx = Xte, s = lambda))), 0), 1)
}

# Binomial glmnet on the actual COUNTS.
# This is the correct likelihood for a prevalence built from n_svy Bernoulli
# draws: it weights each area by its own information instead of treating
# logit(p) as homoskedastic, and it needs no 0.5%/99.5% clamping of the areas
# that came back all-deficient or none-deficient.
fit_glmnet_binom <- function(Xtr, Xte, p_tr, n_tr, alpha = 0.5, lambda = "lambda.min",
                             deff = 1.5, w = NULL, seed = 1L) {
  if (ncol(Xtr) < 2) return(rep(mean(p_tr), nrow(Xte)))
  n_eff <- pmax(round(n_tr / deff), 1)            # design-effect-deflated trials
  succ  <- pmin(pmax(round(p_tr * n_eff), 0), n_eff)
  ymat  <- cbind(failure = n_eff - succ, success = succ)
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, ymat, family = "binomial", alpha = alpha,
                                   nfolds = min(5, nrow(Xtr))), error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(p_tr), nrow(Xte)))
  as.numeric(stats::predict(cv, newx = Xte, s = lambda, type = "response"))
}

fit_ranger <- function(Xtr, Xte, p_tr, n_tr, w = NULL, seed = 1L) {
  if (!requireNamespace("ranger", quietly = TRUE) || ncol(Xtr) < 2)
    return(rep(mean(p_tr), nrow(Xte)))
  df <- data.frame(y = .logit(p_tr), Xtr, check.names = TRUE)
  rf <- tryCatch(ranger::ranger(y ~ ., data = df, num.trees = 800, min.node.size = 5,
                                case.weights = w, seed = seed), error = function(e) NULL)
  if (is.null(rf)) return(rep(mean(p_tr), nrow(Xte)))
  te <- data.frame(Xte, check.names = TRUE)
  .ilogit(stats::predict(rf, data = te)$predictions)
}

# Intercept-only benchmark (the country mean). Any model that cannot beat this
# on RMSE is adding nothing.
fit_null <- function(Xtr, Xte, p_tr, n_tr, w = NULL, ...) rep(mean(p_tr), nrow(Xte))

# Pure spatial smoother on lon/lat -- NO covariates. Tells you how much of a
# model apparent skill is just "neighbouring districts look alike", which the
# pipeline currently never benchmarks against.
fit_spatial <- function(coords_tr, coords_te, p_tr, w = NULL, k = 20) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(rep(mean(p_tr), nrow(coords_te)))
  ok <- is.finite(coords_tr[, 1]) & is.finite(coords_tr[, 2])
  if (sum(ok) < 15) return(rep(mean(p_tr), nrow(coords_te)))
  d <- data.frame(y = .logit(p_tr), lon = coords_tr[, 1], lat = coords_tr[, 2])[ok, ]
  kk <- min(k, max(5, floor(nrow(d) / 3)))
  g <- tryCatch(mgcv::gam(y ~ s(lon, lat, k = kk), data = d,
                          weights = if (is.null(w)) NULL else w[ok]),
                error = function(e) NULL)
  if (is.null(g)) return(rep(mean(p_tr), nrow(coords_te)))
  nd <- data.frame(lon = coords_te[, 1], lat = coords_te[, 2])
  out <- .ilogit(as.numeric(stats::predict(g, newdata = nd)))
  out[!is.finite(out)] <- mean(p_tr); out
}

# -- Metrics ----------------------------------------------------------------
# `prec` = precision weights (1/sampling variance) for the precision-weighted
# correlation: a district estimated from 6 children should not carry the same
# weight in the evaluation as one estimated from 200.
metrics <- function(pred, obs, prec = NULL, r_max = NA_real_) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 4) return(NULL)
  p <- pred[ok]; o <- obs[ok]
  wr <- NA_real_
  if (!is.null(prec)) {
    w <- prec[ok]; w[!is.finite(w) | w <= 0] <- stats::median(w[is.finite(w) & w > 0])
    if (all(is.finite(w)) && stats::sd(p) > 0) {
      mp <- stats::weighted.mean(p, w); mo <- stats::weighted.mean(o, w)
      wr <- sum(w * (p - mp) * (o - mo)) /
        sqrt(sum(w * (p - mp)^2) * sum(w * (o - mo)^2))
    }
  }
  r <- suppressWarnings(stats::cor(p, o))
  data.frame(
    n = sum(ok),
    pearson  = r,
    spearman = suppressWarnings(stats::cor(p, o, method = "spearman")),
    pearson_w = wr,
    r_share  = if (is.finite(r_max) && r_max > .05) r / r_max else NA_real_,
    rmse_pp  = sqrt(mean((p - o)^2)) * 100,
    mae_pp   = mean(abs(p - o)) * 100,
    bias_pp  = mean(p - o) * 100,
    calib    = tryCatch(unname(stats::coef(stats::lm(o ~ p))[2]), error = function(e) NA_real_),
    stringsAsFactors = FALSE)
}

# Sampling variance of an area prevalence.
# The raw plug-in p(1-p)/n collapses to EXACTLY ZERO for the many districts that
# came back 0% (or 100%) deficient -- which at a median n_svy of 6-14 is 11-89%
# of districts. Inverting that gives those districts infinite precision. Use the
# Agresti-Coull style shrunk proportion for the variance instead, so a district
# observed at 0/6 is treated as roughly as informative as 6 observations allow,
# not as if it were measured perfectly.
sampling_var <- function(p, n, deff = 1.5) {
  n <- pmax(n, 1)
  p_s <- (p * n + 0.5) / (n + 1)
  deff * p_s * (1 - p_s) / n
}

# Sampling variance for AVERAGING over districts.
# E[p(1-p)] = pi(1-pi)(n-1)/n for a binomial proportion, so p(1-p)/(n-1) is
# unbiased for pi(1-pi)/n. Individual districts observed at 0% still return 0,
# but the MEAN over districts -- which is all the variance decomposition needs
# -- is unbiased, whereas the Agresti-Coull version above would over-correct it.
sampling_var_unbiased <- function(p, n, deff = 1.5) {
  n <- pmax(n, 2)
  deff * p * (1 - p) / (n - 1)
}

# Reliability (share of between-area variance that is signal) and the implied
# ceiling on Pearson r. See sandbox_parsimony/R/01_noise_audit.R.
reliability <- function(p, n, deff = 1.5) {
  v_obs <- stats::var(p, na.rm = TRUE)
  v_samp <- mean(sampling_var_unbiased(p, n, deff), na.rm = TRUE)
  lam <- max(0, v_obs - v_samp) / v_obs
  list(lambda = lam, r_max = sqrt(lam), v_samp = v_samp)
}

# Repeated K-fold CV over one country Admin-2 rows.
# spec = list(name, features = function(Xtr, ytr, vars) -> character,
#             transform = function(Xtr, Xte) -> list(Xtr, Xte),
#             fit = function(Xtr, Xte, p_tr, n_tr, w, seed) -> prevalence,
#             kind = optional "spatial")
run_cv <- function(dat, vars, spec, n_rep = 10L, k = 5L, deff = 1.5, base_seed = 7L) {
  n <- nrow(dat)
  if (n < 10) return(NULL)
  rel <- reliability(dat$svy_prev, dat$n_svy, deff)
  prec <- 1 / sampling_var(dat$svy_prev, dat$n_svy, deff)
  coords <- as.matrix(dat[, c("lon", "lat")])

  out <- list()
  for (rep_i in seq_len(n_rep)) {
    set.seed(base_seed * 1000L + rep_i)
    fold <- sample(rep_len(seq_len(min(k, n)), n))
    oof <- rep(NA_real_, n)
    for (f in unique(fold)) {
      te <- which(fold == f); tr <- setdiff(seq_len(n), te)
      if (length(tr) < 8) next
      pp <- prep_X(dat[tr, , drop = FALSE], dat[te, , drop = FALSE], vars)
      Xtr <- pp$Xtr; Xte <- pp$Xte; vv <- pp$vars
      if (!is.null(spec$features) && length(vv) > 1) {
        sel <- tryCatch(spec$features(Xtr, .logit(dat$svy_prev[tr]), vv),
                        error = function(e) vv)
        if (length(sel) >= 1) { Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE] }
      }
      if (!is.null(spec$transform)) {
        tf <- tryCatch(spec$transform(Xtr, Xte), error = function(e) NULL)
        if (!is.null(tf)) { Xtr <- tf$Xtr; Xte <- tf$Xte }
      }
      w <- pmax(dat$n_svy[tr], 1)
      pr <- tryCatch(
        if (identical(spec$kind, "spatial"))
          fit_spatial(coords[tr, , drop = FALSE], coords[te, , drop = FALSE],
                      dat$svy_prev[tr], w)
        else spec$fit(Xtr, Xte, dat$svy_prev[tr], dat$n_svy[tr], w = w,
                      seed = base_seed * 1000L + rep_i),
        error = function(e) rep(mean(dat$svy_prev[tr]), length(te)))
      if (length(pr) == length(te)) oof[te] <- pr
    }
    m <- metrics(oof, dat$svy_prev, prec, rel$r_max)
    if (!is.null(m)) { m$rep <- rep_i; out[[rep_i]] <- m }
  }
  if (!length(out)) return(NULL)
  res <- bind_rows(out)
  res$n_pred <- length(vars); res$reliability <- rel$lambda; res$r_max <- rel$r_max
  res
}

summarise_cv <- function(res) {
  if (is.null(res)) return(NULL)
  res |> summarise(
    n = dplyr::first(n), n_pred = dplyr::first(n_pred),
    reliability = round(dplyr::first(reliability), 3),
    r_max = round(dplyr::first(r_max), 3),
    pearson_sd = round(stats::sd(pearson, na.rm = TRUE), 3),
    pearson = round(mean(pearson, na.rm = TRUE), 3),
    spearman = round(mean(spearman, na.rm = TRUE), 3),
    pearson_w = round(mean(pearson_w, na.rm = TRUE), 3),
    r_share = round(mean(r_share, na.rm = TRUE), 3),
    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
    mae_pp = round(mean(mae_pp, na.rm = TRUE), 2),
    calib = round(mean(calib, na.rm = TRUE), 2))
}
