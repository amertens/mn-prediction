# =============================================================================
# sandbox_parsimony/R/17_admin1_layer.R
#
# NEXT STEP 1: an Admin-1 layer alongside Admin-2.
#
# FINDINGS.md section 12 showed the ceiling roughly doubles going from Admin-2
# to Admin-1 (r_max 0.31 -> 0.59) because each unit carries ~3x the sample. But
# a higher ceiling is not automatically a better map: Admin-1 also leaves only
# 8-32 units per country to fit on. This script measures whether the reliability
# gain survives the loss of design points, rather than assuming it does.
#
# It also answers the second question: does letting the data CHOOSE the
# predictors beat the 16 pre-specified constructs? At n = 8-32 that is exactly
# where selection overfits, so selection is done strictly INSIDE each training
# fold -- and a deliberately leaky variant is included to show what the same
# procedure looks like when selection happens once on all the data, which is
# the mistake it is easy to make.
#
# Everything runs on the EXTENDED covariate set from 21_extend_covariates.R.
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr); library(ranger); library(mgcv)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi", "tanzania")
CFG_KEY <- c(gambia = "Gambia", ghana = "Ghana", sierraleone = "SierraLeone",
             malawi = "Malawi", tanzania = "Tanzania")
`%|%` <- function(a, b) if (is.null(a)) b else a

sink(stdout(), type = "message")
source("R/config.R")
cfg <- get_country_configs()
cov_ext <- readRDS("sandbox_parsimony/out/cov_ext.rds")

WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"

# ---------------------------------------------------------------------------
# Aggregate the Admin-2 covariate table up to Admin-1.
# Population-weighted where a population proxy exists: a region's covariate
# should represent where its people are, not the unweighted average of its
# districts (which lets a huge empty district dominate).
# ---------------------------------------------------------------------------
source("sandbox_parsimony/R/17_admin1_layer_fns.R")

# ---------------------------------------------------------------------------
# Feature strategies. Every one of these is called with TRAINING data only.
# ---------------------------------------------------------------------------
f_curated  <- function(X, y, v) intersect(curated_vars(v), v)
f_decorr   <- function(k) function(X, y, v) decorr_reps(X, v, k = k)
f_topK     <- function(k) function(X, y, v) screen_topK(X, y, v, K = k)
f_lasso    <- function(kmax) function(X, y, v) {
  Xs <- X[, v, drop = FALSE]
  if (ncol(Xs) < 2 || nrow(Xs) < 6) return(v)
  cv <- tryCatch(glmnet::cv.glmnet(Xs, y, alpha = 1, nfolds = min(5, nrow(Xs))),
                 error = function(e) NULL)
  if (is.null(cv)) return(utils::head(v, kmax))
  b <- as.matrix(stats::coef(cv, s = "lambda.min"))
  sel <- rownames(b)[-1][b[-1, 1] != 0]
  if (!length(sel)) return(utils::head(v, kmax))
  utils::head(sel, kmax)
}

fit_rf  <- function(Xtr, Xte, p_tr, n_tr, seed = 1L) {
  df <- data.frame(y = .logit(p_tr), Xtr, check.names = TRUE)
  rf <- tryCatch(ranger::ranger(y ~ ., data = df, num.trees = 1000, min.node.size = 3,
                                case.weights = pmax(n_tr, 1), seed = seed),
                 error = function(e) NULL)
  if (is.null(rf)) return(rep(mean(p_tr), nrow(Xte)))
  .ilogit(stats::predict(rf, data = data.frame(Xte, check.names = TRUE))$predictions)
}
fit_rdg <- function(Xtr, Xte, p_tr, n_tr, seed = 1L) {
  if (ncol(Xtr) < 2) return(rep(mean(p_tr), nrow(Xte)))
  set.seed(seed)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, .logit(p_tr), alpha = 0,
                                   weights = pmax(n_tr, 1),
                                   nfolds = min(5, nrow(Xtr))), error = function(e) NULL)
  if (is.null(cv)) return(rep(mean(p_tr), nrow(Xte)))
  pmin(pmax(.ilogit(as.numeric(stats::predict(cv, newx = Xte, s = "lambda.min"))), 0), 1)
}

SPECS <- list(
  list(name = "null_mean",          feat = NULL,          fit = NULL),
  list(name = "curated_xy",         feat = f_curated,     fit = fit_rf,  xy = TRUE),
  list(name = "curated_ridge",      feat = f_curated,     fit = fit_rdg, xy = FALSE),
  list(name = "decorr20_xy",        feat = f_decorr(20),  fit = fit_rf,  xy = TRUE),
  # --- data-adaptive, selection strictly inside the training fold ------------
  list(name = "adaptive_top10_xy",  feat = f_topK(10),    fit = fit_rf,  xy = TRUE),
  list(name = "adaptive_top20_xy",  feat = f_topK(20),    fit = fit_rf,  xy = TRUE),
  list(name = "adaptive_lasso_xy",  feat = f_lasso(20),   fit = fit_rf,  xy = TRUE),
  list(name = "adaptive_top20_ridge", feat = f_topK(20),  fit = fit_rdg, xy = FALSE),
  # --- the same top-20 screen, but chosen ONCE on all the data --------------
  list(name = "LEAKY_top20_xy",     feat = f_topK(20),    fit = fit_rf,  xy = TRUE,
       leaky = TRUE)
)

#' Leave-one-out CV (Admin-1 has too few units for 5-fold) or repeated K-fold
run_level_cv <- function(dat, vars, spec, level, n_rep = 10L) {
  n <- nrow(dat); if (n < 8) return(NULL)
  loo <- level == "Admin1" || n < 30
  rel <- reliability(dat$svy_prev, dat$n_svy, 1.5)
  coords <- as.matrix(dat[, c("lon", "lat")])

  # For the leaky variant, choose the features ONCE using every row, then hand
  # the same fixed set to every fold. This is the optimism being quantified.
  fixed_sel <- NULL
  if (isTRUE(spec$leaky) && !is.null(spec$feat)) {
    pp <- prep_X(dat, dat, vars)
    fixed_sel <- tryCatch(spec$feat(pp$Xtr, .logit(dat$svy_prev), pp$vars),
                          error = function(e) pp$vars)
  }

  one_pass <- function(fold) {
    oof <- rep(NA_real_, n)
    for (f in unique(fold)) {
      te <- which(fold == f); tr <- setdiff(seq_len(n), te)
      if (length(tr) < 5) next
      if (is.null(spec$fit)) { oof[te] <- stats::weighted.mean(dat$svy_prev[tr],
                                                               pmax(dat$n_svy[tr], 1)); next }
      pp <- prep_X(dat[tr, , drop = FALSE], dat[te, , drop = FALSE], vars)
      sel <- if (!is.null(fixed_sel)) intersect(fixed_sel, pp$vars) else
        tryCatch(spec$feat(pp$Xtr, .logit(dat$svy_prev[tr]), pp$vars),
                 error = function(e) pp$vars)
      if (!length(sel)) sel <- pp$vars
      Xtr <- pp$Xtr[, sel, drop = FALSE]; Xte <- pp$Xte[, sel, drop = FALSE]
      if (isTRUE(spec$xy)) {
        Xtr <- cbind(Xtr, lon = coords[tr, 1], lat = coords[tr, 2])
        Xte <- cbind(Xte, lon = coords[te, 1], lat = coords[te, 2])
      }
      pr <- tryCatch(spec$fit(Xtr, Xte, dat$svy_prev[tr], dat$n_svy[tr], seed = f),
                     error = function(e) rep(mean(dat$svy_prev[tr]), length(te)))
      if (length(pr) == length(te)) oof[te] <- pr
    }
    oof
  }

  res <- list()
  reps <- if (loo) 1L else n_rep
  for (i in seq_len(reps)) {
    fold <- if (loo) seq_len(n) else { set.seed(900L + i); sample(rep_len(1:5, n)) }
    oof <- one_pass(fold)
    m <- metrics(oof, dat$svy_prev, NULL, rel$r_max)
    if (!is.null(m)) res[[i]] <- m
  }
  if (!length(res)) return(NULL)
  out <- bind_rows(res)
  out$r_max <- rel$r_max; out$n_units <- n; out$cv <- if (loo) "LOO" else "5-fold x10"
  out
}

# ---------------------------------------------------------------------------
rows <- list(); preds <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  cov2 <- cov_ext[[ctry]]; if (is.null(cov2)) next
  cov1 <- covariates_admin1(cov2)
  for (level in c("Admin1", "Admin2")) {
    cv <- if (level == "Admin1") cov1 else cov2
    if (is.null(cv)) next
    for (oc in cc$outcomes) {
      sv <- outcome_at(ctry, oc, cc, level); if (is.null(sv)) next
      dat <- dplyr::inner_join(sv, cv, by = "Admin2")
      dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy) &
                   is.finite(dat$lon) & is.finite(dat$lat), , drop = FALSE]
      if (nrow(dat) < 8) next
      vars <- setdiff(names(dat), c("Admin2", "Admin1", "svy_prev", "n_svy",
                                    "lon", "lat", "n_admin2", ".w"))
      vars <- vars[vapply(vars, function(v) is.numeric(dat[[v]]) &&
                            sum(is.finite(dat[[v]])) > 2 &&
                            stats::sd(dat[[v]], na.rm = TRUE) > 1e-8, logical(1))]
      if (length(vars) < 5) next
      for (sp in SPECS) {
        r <- run_level_cv(dat, vars, sp, level); if (is.null(r)) next
        rows[[paste(ctry, oc$tag, level, sp$name)]] <- data.frame(
          country = ctry, outcome = oc$tag, level = level, model = sp$name,
          n_units = r$n_units[1], n_pred = length(vars), cv = r$cv[1],
          r_max = round(r$r_max[1], 3),
          pearson = round(mean(r$pearson, na.rm = TRUE), 3),
          spearman = round(mean(r$spearman, na.rm = TRUE), 3),
          rmse_pp = round(mean(r$rmse_pp, na.rm = TRUE), 2),
          mae_pp = round(mean(r$mae_pp, na.rm = TRUE), 2),
          stringsAsFactors = FALSE)
      }
      message(sprintf("  %-11s %-13s %-7s n=%3d p=%4d done", ctry, oc$tag, level,
                      nrow(dat), length(vars)))
    }
  }
}

res <- bind_rows(rows)
write.csv(res, "sandbox_parsimony/out/admin1_vs_admin2.csv", row.names = FALSE)

cat("\n=== Admin-1 vs Admin-2, same models, same countries ===\n")
print(as.data.frame(res |> group_by(level, model) |>
  summarise(cells = n(), median_units = round(median(n_units)),
            r_max = round(mean(r_max), 3),
            pearson = round(mean(pearson, na.rm = TRUE), 3),
            spearman = round(mean(spearman, na.rm = TRUE), 3),
            mae_pp = round(mean(mae_pp, na.rm = TRUE), 2),
            .groups = "drop") |> arrange(level, desc(spearman))), row.names = FALSE)

cat("\n=== Paired by cell: does Admin-1 beat Admin-2? (curated_xy) ===\n")
p <- res |> filter(model == "curated_xy") |>
  select(country, outcome, level, n_units, r_max, spearman, mae_pp) |>
  tidyr::pivot_wider(names_from = level, values_from = c(n_units, r_max, spearman, mae_pp))
print(as.data.frame(p), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/admin1_vs_admin2.csv")
