# =============================================================================
# sandbox_parsimony/R/13_combine_space.R
#
# The binomial GAM (08_spatial_models.R) did NOT beat the simpler models -- at
# 30-170 districts the extra machinery costs more than the correct likelihood
# buys. But the bake-off also showed a covariate-free lon/lat smoother matching
# the 237-predictor production recipe, which says geography carries signal the
# covariates are not capturing.
#
# So: take the models that actually won (random forest on decorrelated
# representatives; ridge on the curated 16) and give them geography, the two
# cheap ways:
#   *_xy       add lon/lat as ordinary predictors
#   *_sstack   fit the covariate model, then smooth its residuals over lon/lat
#              with a thin-plate spline and add the two together
#
# Same repeated 5-fold CV as 04_within_country.R, so rows are comparable.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages(library(mgcv))

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES  <- c("child_vitA", "child_iron", "women_iron", "women_folate", "women_vitA")
MIN_AREAS <- 25L

# Smooth the training residuals over space and add the correction back.
# The spline is fitted to residuals only, so it cannot simply relearn the mean,
# and its own smoothing parameter decides how much spatial structure is left
# once the covariates have had their turn.
spatial_residual_stack <- function(pred_tr, pred_te, y_tr, coords_tr, coords_te) {
  ok <- is.finite(coords_tr[, 1]) & is.finite(coords_tr[, 2])
  if (sum(ok) < 25) return(pred_te)
  r <- .logit(y_tr) - .logit(pmin(pmax(pred_tr, .005), .995))
  d <- data.frame(r = r, lon = coords_tr[, 1], lat = coords_tr[, 2])[ok, ]
  kk <- min(20L, max(5L, floor(nrow(d) / 4)))
  g <- tryCatch(mgcv::gam(r ~ s(lon, lat, k = kk), data = d, method = "REML"),
                error = function(e) NULL)
  if (is.null(g)) return(pred_te)
  adj <- tryCatch(as.numeric(stats::predict(
    g, newdata = data.frame(lon = coords_te[, 1], lat = coords_te[, 2]))),
    error = function(e) rep(0, nrow(coords_te)))
  adj[!is.finite(adj)] <- 0
  pmin(pmax(.ilogit(.logit(pmin(pmax(pred_te, .005), .995)) + adj), 0), 1)
}

# run_cv has no slot for coordinates, so this driver mirrors it and threads them
# through explicitly.
run_cv_xy <- function(dat, vars, spec, n_rep = 5L, k = 5L, deff = 1.5, base_seed = 7L) {
  n <- nrow(dat); if (n < 10) return(NULL)
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
      sel <- tryCatch(spec$features(Xtr, .logit(dat$svy_prev[tr]), vv),
                      error = function(e) vv)
      Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
      if (isTRUE(spec$add_xy)) {
        Xtr <- cbind(Xtr, lon = coords[tr, 1], lat = coords[tr, 2])
        Xte <- cbind(Xte, lon = coords[te, 1], lat = coords[te, 2])
      }
      w <- pmax(dat$n_svy[tr], 1)
      sd_ <- base_seed * 1000L + rep_i
      pr <- tryCatch(spec$fit(Xtr, Xte, dat$svy_prev[tr], dat$n_svy[tr], w = w, seed = sd_),
                     error = function(e) rep(mean(dat$svy_prev[tr]), length(te)))
      if (isTRUE(spec$sstack)) {
        pr_tr <- tryCatch(spec$fit(Xtr, Xtr, dat$svy_prev[tr], dat$n_svy[tr],
                                   w = w, seed = sd_),
                          error = function(e) rep(mean(dat$svy_prev[tr]), length(tr)))
        # NB pr_tr is an in-sample fit, so the residual field is estimated on
        # optimistic residuals and the correction is conservative. Honest
        # nested-CV residuals would be better; this is the cheap version.
        pr <- spatial_residual_stack(pr_tr, pr, dat$svy_prev[tr],
                                     coords[tr, , drop = FALSE],
                                     coords[te, , drop = FALSE])
      }
      if (length(pr) == length(te)) oof[te] <- pr
    }
    m <- metrics(oof, dat$svy_prev, prec, rel$r_max)
    if (!is.null(m)) { m$rep <- rep_i; out[[rep_i]] <- m }
  }
  if (!length(out)) return(NULL)
  res <- dplyr::bind_rows(out)
  res$n_pred <- length(vars); res$reliability <- rel$lambda; res$r_max <- rel$r_max
  res
}

f_decorr20 <- function(X, y, v) decorr_reps(X, v, k = 20L)
f_curated  <- function(X, y, v) intersect(curated_vars(v), v)
rf  <- fit_ranger
rdg <- function(...) fit_glmnet_logit(..., alpha = 0)

SPECS <- list(
  list(name = "decorr20_rf_xy",       features = f_decorr20, fit = rf,  add_xy = TRUE),
  list(name = "decorr20_rf_sstack",   features = f_decorr20, fit = rf,  sstack = TRUE),
  list(name = "curated16_ridge_xy",   features = f_curated,  fit = rdg, add_xy = TRUE),
  list(name = "curated16_ridge_sstack", features = f_curated, fit = rdg, sstack = TRUE),
  list(name = "curated16_rf_xy",      features = f_curated,  fit = rf,  add_xy = TRUE)
)

res_all <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  for (ctry in P$countries) {
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < MIN_AREAS) next
    message(sprintf("\n== %s / %s (n=%d) ==", oc, ctry, nrow(dat)))
    for (sp in SPECS) {
      r <- run_cv_xy(dat, P$predictors, sp)
      if (is.null(r)) next
      s <- summarise_cv(r)
      s$outcome <- oc; s$country <- ctry; s$model <- sp$name
      res_all[[paste(oc, ctry, sp$name)]] <- s
      message(sprintf("   %-24s r=%+.3f (sd %.3f) rho=%+.3f rmse=%5.2f",
                      sp$name, s$pearson, s$pearson_sd, s$spearman, s$rmse_pp))
    }
  }
}

res <- dplyr::bind_rows(res_all)
write.csv(res, "sandbox_parsimony/out/combine_space.csv", row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/combine_space.csv")
