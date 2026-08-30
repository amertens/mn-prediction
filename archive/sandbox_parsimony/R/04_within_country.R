# =============================================================================
# sandbox_parsimony/R/04_within_country.R
#
# Within-country bake-off for the GEE/proxy-only Admin-2 model.
#
# Every model is scored by the SAME repeated 5-fold CV (10 repeats, fixed
# seeds), so differences are not fold-luck. Reported alongside:
#   r_max    -- ceiling imposed by survey sampling noise in the outcome
#   r_share  -- pearson / r_max, i.e. share of the ACHIEVABLE signal captured
#   pearson_w-- precision-weighted r (tiny-n districts down-weighted)
#   null / spatial-only rows -- the benchmarks the pipeline currently lacks
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")

OUTCOMES  <- c("child_vitA", "child_iron", "women_iron", "women_folate", "women_vitA")
MIN_AREAS <- 25L

SPECS <- list(
  list(name = "null_mean",         fit = fit_null,          features = NULL),
  list(name = "spatial_only",      fit = fit_null,          features = NULL, kind = "spatial"),

  # --- current production recipe and its stated max-performance alternative --
  list(name = "PROD enet_screen30",
       features = function(X, y, v) screen_topK(X, y, v, K = 30L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "PROD ridge_all",
       features = NULL,
       fit = function(...) fit_glmnet_logit(..., alpha = 0)),

  # --- likelihood fix: binomial on counts instead of logit-of-proportion -----
  list(name = "binom_screen30",
       features = function(X, y, v) screen_topK(X, y, v, K = 30L),
       fit = function(...) fit_glmnet_binom(..., alpha = 0.5)),
  list(name = "binom_ridge_all",
       features = NULL,
       fit = function(...) fit_glmnet_binom(..., alpha = 0)),

  # --- parsimony: decorrelated representatives instead of top-K duplicates ---
  list(name = "decorr20_enet",
       features = function(X, y, v) decorr_reps(X, v, k = 20L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "decorr40screen12_enet",
       features = function(X, y, v) decorr_then_screen(X, y, v, k_clust = 40L, K = 12L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "decorr40screen12_binom",
       features = function(X, y, v) decorr_then_screen(X, y, v, k_clust = 40L, K = 12L),
       fit = function(...) fit_glmnet_binom(..., alpha = 0.5)),

  # --- parsimony: a priori mechanistic set (~16 constructs) -----------------
  list(name = "curated16_enet",
       features = function(X, y, v) intersect(curated_vars(v), v),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "curated16_binom",
       features = function(X, y, v) intersect(curated_vars(v), v),
       fit = function(...) fit_glmnet_binom(..., alpha = 0.5)),
  list(name = "curated16_ridge",
       features = function(X, y, v) intersect(curated_vars(v), v),
       fit = function(...) fit_glmnet_logit(..., alpha = 0)),

  # --- parsimony: one PC per conceptual domain (~11 features) ---------------
  list(name = "domainPC1_enet",
       features = NULL,
       transform = function(Xtr, Xte) domain_pcs(Xtr, Xte, colnames(Xtr), npc = 1L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0.5)),
  list(name = "domainPC2_ridge",
       features = NULL,
       transform = function(Xtr, Xte) domain_pcs(Xtr, Xte, colnames(Xtr), npc = 2L),
       fit = function(...) fit_glmnet_logit(..., alpha = 0)),

  # --- nonlinear on a parsimonious set -------------------------------------
  list(name = "decorr20_rf",
       features = function(X, y, v) decorr_reps(X, v, k = 20L),
       fit = fit_ranger)
)

res_all <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  vars <- P$predictors
  for (ctry in P$countries) {
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < MIN_AREAS) next
    message(sprintf("\n== %s / %s  (n=%d areas, p=%d) ==", oc, ctry, nrow(dat), length(vars)))
    for (sp in SPECS) {
      t0 <- Sys.time()
      r <- run_cv(dat, vars, sp, n_rep = 5L, k = 5L)
      if (is.null(r)) next
      s <- summarise_cv(r)
      s$outcome <- oc; s$country <- ctry; s$model <- sp$name
      s$secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
      res_all[[paste(oc, ctry, sp$name)]] <- s
      message(sprintf("   %-24s r=%+.3f (sd %.3f)  rho=%+.3f  r_share=%+.2f  rmse=%5.2f  [%.0fs]",
                      sp$name, s$pearson, s$pearson_sd, s$spearman,
                      ifelse(is.na(s$r_share), NA, s$r_share), s$rmse_pp, s$secs))
    }
  }
}

res <- dplyr::bind_rows(res_all) |>
  dplyr::select(outcome, country, model, n, n_pred, reliability, r_max,
                pearson, pearson_sd, spearman, pearson_w, r_share,
                rmse_pp, mae_pp, calib, secs)
write.csv(res, "sandbox_parsimony/out/within_country_bakeoff.csv", row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/within_country_bakeoff.csv")
