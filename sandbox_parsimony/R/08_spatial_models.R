# =============================================================================
# sandbox_parsimony/R/08_spatial_models.R
#
# The headline candidate: parsimonious covariates + a spatial field, fitted
# with the BINOMIAL likelihood on the actual survey counts.
#
#   y_i ~ Binomial(n_i / deff, p_i)
#   logit(p_i) = beta0 + X_i beta + f(lon_i, lat_i)
#
# Why this and not the current area-level elastic net on logit(prevalence):
#   * binomial likelihood weights each district by its own information, so the
#     n=6 districts stop shouting as loudly as the n=200 ones, and the areas
#     that came back 0% or 100% need no arbitrary clamping;
#   * the spatial field captures the district-to-district structure that the
#     covariates demonstrably miss -- in the bake-off a covariate-free lon/lat
#     smoother already beats the 237-predictor production recipe in several
#     country x outcome cells;
#   * shrinking the covariate set first is what makes a GAM feasible at all at
#     n = 30-170 districts.
#
# Scored on the SAME repeated 5-fold CV as 04_within_country.R, so the rows can
# be pasted straight into that table.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages(library(mgcv))

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
OUTCOMES  <- c("child_vitA", "child_iron", "women_iron", "women_folate", "women_vitA")
MIN_AREAS <- 25L
DEFF <- 1.5

#' Binomial GAM: covariates (linear, ridge-penalised) + optional 2-D spatial field
#' @param sp_k basis size for s(lon,lat); NULL disables the spatial term
fit_binom_gam <- function(Xtr, Xte, p_tr, n_tr, coords_tr, coords_te,
                          sp_k = 20L, deff = DEFF) {
  n_eff <- pmax(round(n_tr / deff), 1)
  succ  <- pmin(pmax(round(p_tr * n_eff), 0), n_eff)
  Ztr <- scale(Xtr)
  ctr <- attr(Ztr, "scaled:center"); scl <- attr(Ztr, "scaled:scale")
  scl[!is.finite(scl) | scl < 1e-8] <- 1
  Ztr <- as.data.frame(sweep(sweep(Xtr, 2, ctr, "-"), 2, scl, "/"))
  Zte <- as.data.frame(sweep(sweep(Xte, 2, ctr, "-"), 2, scl, "/"))
  nmz <- paste0("z", seq_len(ncol(Ztr))); names(Ztr) <- names(Zte) <- nmz

  Ztr$lon <- coords_tr[, 1]; Ztr$lat <- coords_tr[, 2]
  Zte$lon <- coords_te[, 1]; Zte$lat <- coords_te[, 2]
  Ztr$Y <- cbind(succ, n_eff - succ)

  sp_ok <- !is.null(sp_k) && all(is.finite(coords_tr)) && nrow(Ztr) >= 25
  kk <- if (sp_ok) min(sp_k, max(5, floor(nrow(Ztr) / 4))) else NA
  # ridge-penalise the linear block (paraPen) so p ~ 12-25 covariates stay
  # estimable at n ~ 30 districts; the spatial field gets its own GCV-tuned
  # smoothing parameter, so covariates and space compete on equal footing.
  form <- if (sp_ok)
    stats::as.formula(paste("Y ~", paste(nmz, collapse = " + "),
                            sprintf("+ s(lon, lat, k = %d)", kk)))
  else stats::as.formula(paste("Y ~", paste(nmz, collapse = " + ")))

  S <- list(diag(length(nmz)))
  # women_vitA in Ghana/Malawi is ~85-90% zero-prevalence districts, which makes
  # the binomial fit near-separable; without an iteration cap mgcv can spin for
  # tens of minutes on a cell that has no signal to find anyway.
  ctrl <- mgcv::gam.control(maxit = 50L, mgcv.tol = 1e-5, nthreads = 1L)
  g <- tryCatch(
    mgcv::gam(form, family = stats::binomial(), data = Ztr, method = "REML",
              paraPen = list(X = S), select = TRUE, control = ctrl),
    error = function(e) NULL)
  if (is.null(g))
    g <- tryCatch(mgcv::gam(form, family = stats::binomial(), data = Ztr,
                            method = "REML", control = ctrl),
                  error = function(e) NULL)
  if (is.null(g)) return(rep(mean(p_tr), nrow(Xte)))
  out <- tryCatch(as.numeric(stats::predict(g, newdata = Zte, type = "response")),
                  error = function(e) rep(mean(p_tr), nrow(Xte)))
  out[!is.finite(out)] <- mean(p_tr)
  pmin(pmax(out, 0), 1)
}

# run_cv passes `spec$fit(Xtr, Xte, p_tr, n_tr, w, seed)` and has no slot for
# coordinates, so carry them through a spec-level closure instead.
make_gam_spec <- function(name, feat, sp_k) {
  list(name = name, features = feat, kind = "gam", sp_k = sp_k)
}

run_cv_gam <- function(dat, vars, spec, n_rep = 5L, k = 5L, deff = DEFF, base_seed = 7L) {
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
      if (length(sel) >= 1) { Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE] }
      pr <- tryCatch(fit_binom_gam(Xtr, Xte, dat$svy_prev[tr], dat$n_svy[tr],
                                   coords[tr, , drop = FALSE],
                                   coords[te, , drop = FALSE], sp_k = spec$sp_k),
                     error = function(e) rep(mean(dat$svy_prev[tr]), length(te)))
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

f_decorr20  <- function(X, y, v) decorr_reps(X, v, k = 20L)
f_decorr12  <- function(X, y, v) decorr_reps(X, v, k = 12L)
f_curated   <- function(X, y, v) intersect(curated_vars(v), v)
f_d40s12    <- function(X, y, v) decorr_then_screen(X, y, v, 40L, 12L)

SPECS <- list(
  make_gam_spec("gam_spatial_only",     function(X, y, v) v[1], 25L),
  make_gam_spec("gam_decorr12",         f_decorr12, NULL),
  make_gam_spec("gam_decorr12_spatial", f_decorr12, 25L),
  make_gam_spec("gam_decorr20_spatial", f_decorr20, 25L),
  make_gam_spec("gam_curated",          f_curated,  NULL),
  make_gam_spec("gam_curated_spatial",  f_curated,  25L),
  make_gam_spec("gam_d40s12_spatial",   f_d40s12,   25L)
)

res_all <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  for (ctry in P$countries) {
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < MIN_AREAS) next
    # Skip cells the noise audit says have nothing to find. Besides saving time,
    # it avoids the near-separable binomial fits (women_vitA is 85-90% zero-
    # prevalence districts in Ghana and Malawi) that make mgcv spin.
    rel0 <- reliability(dat$svy_prev, dat$n_svy, DEFF)
    if (rel0$r_max < 0.15) {
      message(sprintf("\n== %s / %s (n=%d) -- SKIPPED, r_max = %.2f, no signal to fit",
                      oc, ctry, nrow(dat), rel0$r_max))
      next
    }
    message(sprintf("\n== %s / %s (n=%d, r_max=%.2f) ==",
                    oc, ctry, nrow(dat), rel0$r_max))
    for (sp in SPECS) {
      t0 <- Sys.time()
      r <- run_cv_gam(dat, P$predictors, sp)
      if (is.null(r)) next
      s <- summarise_cv(r)
      s$outcome <- oc; s$country <- ctry; s$model <- sp$name
      s$secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
      res_all[[paste(oc, ctry, sp$name)]] <- s
      message(sprintf("   %-22s r=%+.3f (sd %.3f) rho=%+.3f r_share=%+.2f rmse=%5.2f [%.0fs]",
                      sp$name, s$pearson, s$pearson_sd, s$spearman,
                      ifelse(is.na(s$r_share), NA, s$r_share), s$rmse_pp, s$secs))
      # checkpoint after every fit so a hang later does not lose the run
      utils::write.csv(dplyr::bind_rows(res_all),
                       "sandbox_parsimony/out/spatial_bakeoff_partial.csv",
                       row.names = FALSE)
    }
  }
}

res <- dplyr::bind_rows(res_all) |>
  dplyr::select(outcome, country, model, n, n_pred, reliability, r_max,
                pearson, pearson_sd, spearman, pearson_w, r_share,
                rmse_pp, mae_pp, calib, secs)
write.csv(res, "sandbox_parsimony/out/spatial_bakeoff.csv", row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/spatial_bakeoff.csv")
