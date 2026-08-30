# =============================================================================
# scripts/covariates/13_resolution_comparison.R
#
# DOES THE COVARIATE MODEL BEAT A NULL AND A SPATIAL SMOOTHER AT ADMIN-1
# WHEN IT CANNOT AT ADMIN-2?
#
# That is a claim about WHERE THE SIGNAL LIVES, and it is only worth making if
# the two resolutions are tested identically. So this script runs the SAME
# three arms, through the SAME leave-one-unit-out loop, on the SAME countries
# and outcomes, at both resolutions:
#
#   null      predict the (precision-weighted) mean of the training units
#   spatial   thin-plate spline on unit centroids, NO covariates
#   covariate the pipeline's NNLS ensemble on the harmonised covariates
#
# Everything that could differ between the two resolutions -- learner, folds,
# screening, weighting, scoring -- is held fixed. The only thing that changes
# is the geographic unit.
#
# WHAT THIS CANNOT DO, AND WHY IT IS STATED UP FRONT
# --------------------------------------------------
# Admin-1 buys sample size per unit but destroys the number of units:
#   Gambia 30 -> 6, Ghana 75 -> 16, Sierra Leone 14 -> 4, Malawi 87 -> 27.
# Leave-one-out on 4 or 6 units estimates a correlation from 4 or 6 points; the
# result is dominated by which unit is held out. So:
#   - Sierra Leone (4) is NOT ESTIMABLE at admin-1. Reported as such, not
#     silently dropped.
#   - Gambia (6) is reported but flagged unreliable, and gets no spatial arm.
#   - Ghana (16) and Malawi (27) carry the claim.
# A "we beat the null at admin-1" headline resting on two countries must say
# so. MIN_UNITS_* below encode these thresholds explicitly.
#
#   Rscript scripts/covariates/13_resolution_comparison.R
# -> results/tables/resolution_comparison.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")

MIN_UNITS_CV      <- 8L    # below this, leave-one-out r is not interpretable
MIN_UNITS_SPATIAL <- 12L   # below this, a 2-D spline has nothing to fit
K_SCREEN          <- 20L

# ── One leave-one-unit-out evaluation of all three arms on one frame ─────────
# `d` must carry: unit, prev, n, lon, lat, + numeric covariate columns.
eval_arms <- function(d, covs, label) {
  n <- nrow(d)
  y <- stats::qlogis(pmin(pmax(d$prev, 0.005), 0.995))
  w <- area_precision_weights(d$prev, d$n)
  X <- as.matrix(d[, covs, drop = FALSE])
  X[!is.finite(X)] <- NA
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  keep <- apply(X, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]

  oof_null <- oof_spat <- oof_cov <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    # null: precision-weighted mean of the TRAINING units only
    oof_null[i] <- stats::weighted.mean(d$prev[tr], w[tr])
    # spatial: thin-plate spline on centroids, no covariates
    if (n >= MIN_UNITS_SPATIAL && requireNamespace("mgcv", quietly = TRUE)) {
      kk <- min(20L, max(4L, floor(n / 3)))
      g <- tryCatch(mgcv::gam(y[tr] ~ s(lon, lat, k = kk), data = d[tr, ],
                              weights = pmax(d$n[tr], 1), method = "REML"),
                    error = function(e) NULL)
      if (!is.null(g))
        oof_spat[i] <- stats::plogis(as.numeric(
          stats::predict(g, newdata = d[i, , drop = FALSE])))
    }
    # covariate ensemble: screen INSIDE the fold, then the pipeline's stack
    if (ncol(X) >= 2) {
      sel <- .awsl_screen(X[tr, , drop = FALSE], y[tr], K_SCREEN)
      s <- tryCatch(.awsl_stack(X[tr, sel, drop = FALSE], y[tr], w[tr]),
                    error = function(e) NULL)
      if (!is.null(s))
        oof_cov[i] <- stats::plogis(.awsl_predict(s, X[i, sel, drop = FALSE]))
    }
  }

  met <- function(p, arm) {
    ok <- is.finite(p) & is.finite(d$prev)
    data.frame(arm = arm, n_units = n, n_scored = sum(ok),
               # a constant prediction has no correlation; NA, not 0
               r = if (sum(ok) > 3 && stats::sd(p[ok]) > 0)
                     round(suppressWarnings(stats::cor(d$prev[ok], p[ok])), 3)
                   else NA_real_,
               mae_pp = if (any(ok)) round(100 * mean(abs(d$prev[ok] - p[ok])), 2)
                        else NA_real_,
               stringsAsFactors = FALSE)
  }
  out <- rbind(met(oof_null, "null (train mean)"),
               met(oof_spat, "spatial only"),
               met(oof_cov,  "covariates"))
  out$level <- label
  # Keep the per-unit out-of-fold predictions: the WHO risk-category accuracy
  # below (the quantity the Ghana deck reports at admin-1 vs admin-2) needs
  # predictions, not summary metrics.
  attr(out, "preds") <- data.frame(
    unit = d$unit, obs = d$prev, null = oof_null,
    spatial = oof_spat, covariates = oof_cov, stringsAsFactors = FALSE)
  out
}

# ── WHO vitamin-A / severity banding, for the deck's "exact match" table ─────
# Thresholds are the WHO public-health-significance bands for VAD.
who_band <- function(p) cut(100 * p, breaks = c(-Inf, 2, 10, 20, Inf),
                            labels = c("None", "Mild", "Moderate", "Severe"),
                            right = FALSE)
band_acc <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 3) return(c(exact = NA_real_, within1 = NA_real_, n = sum(ok)))
  a <- as.integer(who_band(obs[ok])); b <- as.integer(who_band(pred[ok]))
  c(exact = mean(a == b), within1 = mean(abs(a - b) <= 1), n = sum(ok))
}

configs <- get_country_configs()
rows <- list()
preds <- list()

for (cn in names(configs)) {
  cc <- configs[[cn]]; lc <- tolower(cn)
  acov <- tryCatch(targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE),
                   error = function(e) NULL)
  if (is.null(acov)) acov <- tryCatch(
    targets::tar_read_raw(paste0("gee_admin2_", lc), store = STORE),
    error = function(e) NULL)
  cen <- tryCatch(sf::st_drop_geometry(load_admin2_centroids(cc$gadm_code)),
                  error = function(e) NULL)
  # area_covariates_* is a LIST(gee_admin2, polygons); unwrap before use.
  if (is.list(acov) && !is.data.frame(acov) && "gee_admin2" %in% names(acov))
    acov <- acov$gee_admin2
  if (is.null(acov) || is.null(cen)) { message("skip ", cn, ": no covariates/centroids"); next }

  for (on in names(cc$outcomes)) {
    oc  <- cc$outcomes[[on]]
    svy <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", on), store = STORE),
                    error = function(e) NULL)
    od  <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE),
                    error = function(e) NULL)
    if (is.null(svy) || !nrow(svy)) next

    # ---- ADMIN-2 ----------------------------------------------------------
    by2 <- admin2_join_by(svy, acov)
    d2  <- dplyr::inner_join(
      svy[, c(by2, "svy_prev", "n_svy")], acov, by = by2)
    cby <- admin2_join_by(d2, cen)
    d2  <- dplyr::left_join(d2, cen[, c(cby, "lon", "lat")], by = cby)
    d2  <- d2[is.finite(d2$svy_prev) & is.finite(d2$lon), , drop = FALSE]
    covs2 <- setdiff(names(d2), c("Admin1", "Admin2", "svy_prev", "n_svy", "lon", "lat"))
    covs2 <- covs2[vapply(covs2, function(v) is.numeric(d2[[v]]), logical(1))]
    if (nrow(d2) >= MIN_UNITS_CV && length(covs2) >= 2) {
      r2 <- eval_arms(data.frame(unit = d2$Admin2, prev = d2$svy_prev, n = d2$n_svy,
                                 lon = d2$lon, lat = d2$lat,
                                 d2[, covs2, drop = FALSE], check.names = FALSE),
                      covs2, "admin-2")
      rows[[length(rows) + 1L]] <- cbind(country = cn, outcome = on, r2)
      pr <- attr(r2, "preds")
      preds[[length(preds) + 1L]] <- cbind(country = cn, outcome = on,
                                           level = "admin-2", pr)
    }

    # ---- ADMIN-1 ----------------------------------------------------------
    a1 <- tryCatch(admin1_design_based(od, cc, oc), error = function(e) NULL)
    if (is.null(a1) || !nrow(a1)) next
    names(a1)[names(a1) == "prev"] <- "svy_prev"
    x1 <- tryCatch(aggregate_covariates_to_admin1(acov), error = function(e) NULL)
    c1 <- cen %>% group_by(Admin1) %>%
      summarise(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE),
                .groups = "drop")
    if (is.null(x1)) next
    d1 <- a1 %>% inner_join(x1, by = "Admin1") %>% left_join(c1, by = "Admin1")
    nm <- intersect(c("svy_prev", "prev"), names(d1))[1]
    nn <- intersect(c("n", "n_svy"), names(d1))[1]
    d1 <- d1[is.finite(d1[[nm]]) & is.finite(d1$lon), , drop = FALSE]
    covs1 <- setdiff(names(d1), c("Admin1", nm, nn, "lon", "lat", "se", "n"))
    covs1 <- covs1[vapply(covs1, function(v) is.numeric(d1[[v]]), logical(1))]

    if (nrow(d1) < MIN_UNITS_CV) {
      rows[[length(rows) + 1L]] <- data.frame(
        country = cn, outcome = on, arm = "NOT ESTIMABLE",
        n_units = nrow(d1), n_scored = 0L, r = NA_real_, mae_pp = NA_real_,
        level = "admin-1", stringsAsFactors = FALSE)
      next
    }
    r1 <- eval_arms(data.frame(unit = d1$Admin1, prev = d1[[nm]], n = d1[[nn]],
                               lon = d1$lon, lat = d1$lat,
                               d1[, covs1, drop = FALSE], check.names = FALSE),
                    covs1, "admin-1")
    rows[[length(rows) + 1L]] <- cbind(country = cn, outcome = on, r1)
    pr <- attr(r1, "preds")
    preds[[length(preds) + 1L]] <- cbind(country = cn, outcome = on,
                                         level = "admin-1", pr)
  }
  message("done: ", cn)
}

res <- dplyr::bind_rows(rows)
dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(res, here("results", "tables", "resolution_comparison.csv"))

cat("\n================ COVARIATES vs NULL vs SPATIAL, BY RESOLUTION ==========\n")
for (lv in c("admin-2", "admin-1")) {
  s <- res[res$level == lv & res$arm != "NOT ESTIMABLE", ]
  if (!nrow(s)) next
  cat(sprintf("\n--- %s ---\n", lv))
  print(s %>% group_by(arm) %>%
          summarise(cells = n(), mean_r = round(mean(r, na.rm = TRUE), 3),
                    med_r = round(stats::median(r, na.rm = TRUE), 3),
                    mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                    .groups = "drop") %>% as.data.frame(), row.names = FALSE)
  w <- s %>% select(country, outcome, arm, mae_pp) %>%
    tidyr::pivot_wider(names_from = arm, values_from = mae_pp)
  if (all(c("covariates", "null (train mean)") %in% names(w)))
    cat(sprintf("  covariates beat null:    %d/%d cells\n",
                sum(w$covariates < w$`null (train mean)`, na.rm = TRUE),
                sum(is.finite(w$covariates) & is.finite(w$`null (train mean)`))))
  if (all(c("covariates", "spatial only") %in% names(w)))
    cat(sprintf("  covariates beat spatial: %d/%d cells\n",
                sum(w$covariates < w$`spatial only`, na.rm = TRUE),
                sum(is.finite(w$covariates) & is.finite(w$`spatial only`))))
}
ne <- res[res$arm == "NOT ESTIMABLE", ]
if (nrow(ne)) cat(sprintf("\nNOT ESTIMABLE at admin-1 (< %d units): %s\n", MIN_UNITS_CV,
                          paste(unique(ne$country), collapse = ", ")))
cat("\n-> results/tables/resolution_comparison.csv\n")

# ── HYBRID: pick the level PER COUNTRY, a priori ─────────────────────────────
# The rule is fixed BEFORE looking at performance: use admin-1 wherever it has
# at least MIN_UNITS_CV units, otherwise fall back to admin-2. It therefore
# depends only on how a country is subdivided -- Gambia (6 regions) and Sierra
# Leone (4) get admin-2, Ghana (16) and Malawi (27) get admin-1.
#
# THIS IS THE POINT: choosing the level by which one SCORED BETTER would be
# selection on the test data -- precisely the optimism this project measured at
# +0.182 r across 20/20 cells. A hybrid picked that way would inherit exactly
# that bias. Unit count is observable without the outcome, so it is a
# legitimate basis for the choice.
pr <- dplyr::bind_rows(preds)
readr::write_csv(pr, here("results", "tables", "resolution_predictions.csv"))

units <- res %>% filter(arm != "NOT ESTIMABLE") %>%
  select(country, level, n_units) %>% distinct()
a1n <- units %>% filter(level == "admin-1") %>% select(country, n1 = n_units)
choice <- units %>% distinct(country) %>% left_join(a1n, by = "country") %>%
  mutate(chosen = ifelse(!is.na(n1) & n1 >= MIN_UNITS_CV, "admin-1", "admin-2"))
cat("\n================ HYBRID (level chosen a priori by unit count) =========\n")
print(as.data.frame(choice), row.names = FALSE)

hyb <- res %>% filter(arm != "NOT ESTIMABLE") %>%
  inner_join(choice %>% select(country, chosen), by = "country") %>%
  filter(level == chosen)
cat("\n--- hybrid, pooled over countries at their chosen level ---\n")
print(hyb %>% group_by(arm) %>%
        summarise(cells = n(), mean_r = round(mean(r, na.rm = TRUE), 3),
                  med_r = round(stats::median(r, na.rm = TRUE), 3),
                  mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)
wh <- hyb %>% select(country, outcome, arm, mae_pp) %>%
  tidyr::pivot_wider(names_from = arm, values_from = mae_pp)
if (all(c("covariates", "spatial only") %in% names(wh)))
  cat(sprintf("  covariates beat spatial: %d/%d cells\n",
              sum(wh$covariates < wh$`spatial only`, na.rm = TRUE),
              sum(is.finite(wh$covariates) & is.finite(wh$`spatial only`))))

# ── WHO risk-category accuracy: the Ghana deck reports this admin-1 vs -2 ────
cat("\n================ WHO RISK-CATEGORY ACCURACY (child vitamin A) =========\n")
cat("deck claimed: admin-2 exact 44%/31%, within-1 75%/49%;",
    "admin-1 exact 56%/56%, within-1 100%/69% (full model / proxies only)\n\n")
cv <- pr %>% filter(outcome == "child_vitA")
if (nrow(cv)) {
  bands <- lapply(split(cv, cv$level, drop = TRUE), function(z) {
    a <- band_acc(z$obs, z$covariates); b <- band_acc(z$obs, z$null)
    data.frame(level = z$level[1], n = as.integer(a[["n"]]),
               cov_exact = round(100 * a[["exact"]]),
               cov_within1 = round(100 * a[["within1"]]),
               null_exact = round(100 * b[["exact"]]),
               null_within1 = round(100 * b[["within1"]]))
  })
  print(dplyr::bind_rows(bands), row.names = FALSE)
} else cat("  no child_vitA predictions available\n")
cat("\n-> results/tables/resolution_predictions.csv\n")
