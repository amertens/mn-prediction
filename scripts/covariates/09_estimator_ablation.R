# =============================================================================
# scripts/covariates/09_estimator_ablation.R
#
# Settle the estimator choices BEFORE the pipeline rebuild, so the rebuild
# produces the final method rather than an intermediate one.
#
# FACTORS TESTED (one at a time from a baseline that mirrors fit_area_sl)
#   learners     the current 5-learner area stack vs trimmed alternatives
#   weighting    equal area weights vs precision weights (inverse sampling var)
#   AlphaEarth   64 free dimensions vs an 8-PC block vs excluded
#   screening    all covariates vs a top-k prescreen, screened on the ADMIN-1
#                target rather than the noisy admin-2 one
#
# HONEST EVALUATION
# -----------------
# Leave-one-REGION-out at admin-2: every fold holds out a whole admin-1 region,
# so no district is predicted by a model that saw its neighbours. Screening and
# ensemble weighting happen strictly inside the training folds. Metrics are
# reported against the admin-2 direct estimates with r_max alongside, because
# r_share -- not r -- is the quantity that says whether a change helped.
#
# FIDELITY CAVEAT
# ---------------
# This harness re-implements the area-level ensemble (same five learner
# families, NNLS-stacked on inner out-of-fold predictions) rather than calling
# mlr3superlearner, because the production path costs hours across this many
# configurations. The winning configuration is re-checked against the real
# fit_area_sl in 10_confirm_ablation.R before anything is adopted.
#
# Output: results/tables/estimator_ablation.csv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
targets::tar_source("R/")

STORE <- here("_targets_full")
rd <- function(nm) tryCatch(targets::tar_read_raw(nm, store = STORE), error = function(e) NULL)
H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))
AEF  <- grep("^aef_A", COVS, value = TRUE)
cfgs <- get_country_configs()
set.seed(20260829L)

# ── learners ────────────────────────────────────────────────────────────────
.fit_glmnet <- function(X, y, w, alpha) {
  f <- tryCatch(glmnet::cv.glmnet(X, y, alpha = alpha, weights = w,
                                  nfolds = max(3L, min(8L, floor(length(y) / 4)))),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, newx = Z, s = "lambda.min")))
}
.fit_ranger <- function(X, y, w) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  d <- data.frame(.y = y, X)
  f <- tryCatch(ranger::ranger(.y ~ ., data = d, num.trees = 500, min.node.size = 5,
                               case.weights = w, num.threads = 2),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, data = as.data.frame(Z))$predictions))
}
.fit_xgb <- function(X, y, w) {
  if (!requireNamespace("xgboost", quietly = TRUE)) return(NULL)
  f <- tryCatch(xgboost::xgboost(data = X, label = y, weight = w, max_depth = 2,
                                 eta = 0.05, nrounds = 200, min_child_weight = 3,
                                 subsample = 0.8, verbose = 0, nthread = 2),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(predict = function(Z) as.numeric(stats::predict(f, Z)))
}
LEARNERS <- list(
  lasso = function(X, y, w) .fit_glmnet(X, y, w, 1),
  ridge = function(X, y, w) .fit_glmnet(X, y, w, 0),
  enet  = function(X, y, w) .fit_glmnet(X, y, w, 0.5),
  ranger = .fit_ranger,
  xgb    = .fit_xgb)

#' NNLS-stacked ensemble over a named learner set, weights from inner OOF.
.fit_stack <- function(X, y, w, which) {
  ls <- LEARNERS[which]
  n <- length(y)
  if (length(ls) == 1L) {
    m <- ls[[1]](X, y, w)
    if (is.null(m)) return(NULL)
    return(list(models = list(m), wts = 1))
  }
  k <- max(3L, min(5L, floor(n / 5)))
  fold <- sample(rep_len(seq_len(k), n))
  P <- matrix(NA_real_, n, length(ls), dimnames = list(NULL, names(ls)))
  for (f in seq_len(k)) {
    i <- which(fold == f)
    for (j in seq_along(ls)) {
      m <- ls[[j]](X[-i, , drop = FALSE], y[-i], w[-i])
      if (!is.null(m)) P[i, j] <- m$predict(X[i, , drop = FALSE])
    }
  }
  ok <- stats::complete.cases(P)
  wts <- rep(1 / length(ls), length(ls))
  if (sum(ok) > length(ls) + 2) {
    nn <- if (requireNamespace("nnls", quietly = TRUE))
      tryCatch(nnls::nnls(P[ok, , drop = FALSE], y[ok]), error = function(e) NULL) else NULL
    if (!is.null(nn) && sum(nn$x) > 0) {
      wts <- nn$x / sum(nn$x)
    } else {
      # Non-negative least squares is the SuperLearner convention; clipping an
      # OLS stack at zero is the standard fallback when nnls is unavailable.
      cf <- tryCatch(stats::coef(stats::lm.fit(P[ok, , drop = FALSE], y[ok])),
                     error = function(e) NULL)
      if (!is.null(cf)) { cf[!is.finite(cf) | cf < 0] <- 0
                          if (sum(cf) > 0) wts <- cf / sum(cf) }
    }
  }
  models <- lapply(ls, function(L) L(X, y, w))
  keep <- !vapply(models, is.null, logical(1))
  if (!any(keep)) return(NULL)
  list(models = models[keep], wts = wts[keep] / max(sum(wts[keep]), 1e-12))
}
.pred_stack <- function(s, Z) {
  if (is.null(s)) return(rep(NA_real_, nrow(Z)))
  P <- vapply(s$models, function(m) m$predict(Z), numeric(nrow(Z)))
  if (is.null(dim(P))) P <- matrix(P, nrow = nrow(Z))
  as.numeric(P %*% s$wts)
}

# ── design matrix construction per configuration ────────────────────────────
.build_X <- function(df, covs, aef_mode, pcs = NULL) {
  use <- covs
  if (aef_mode == "drop") use <- setdiff(use, AEF)
  X <- as.matrix(df[, use, drop = FALSE])
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  if (aef_mode == "pca") {
    A <- as.matrix(df[, AEF, drop = FALSE]); A[!is.finite(A)] <- 0
    if (is.null(pcs)) pcs <- stats::prcomp(A, center = TRUE, scale. = FALSE)
    S <- stats::predict(pcs, A)[, seq_len(min(8L, ncol(pcs$rotation))), drop = FALSE]
    colnames(S) <- paste0("aefPC", seq_len(ncol(S)))
    X <- cbind(X[, setdiff(colnames(X), AEF), drop = FALSE], S)
  }
  list(X = X, pcs = pcs)
}

CONFIGS <- list(
  list(id = "baseline (current stack)",      learners = names(LEARNERS), wt = FALSE, aef = "all",  screen = 0),
  list(id = "enet only",                     learners = "enet",           wt = FALSE, aef = "all",  screen = 0),
  list(id = "enet + ranger",                 learners = c("enet","ranger"), wt = FALSE, aef = "all", screen = 0),
  list(id = "linear only (lasso/ridge/enet)", learners = c("lasso","ridge","enet"), wt = FALSE, aef = "all", screen = 0),
  list(id = "baseline + precision weights",  learners = names(LEARNERS), wt = TRUE,  aef = "all",  screen = 0),
  list(id = "baseline + AEF as 8 PCs",       learners = names(LEARNERS), wt = FALSE, aef = "pca",  screen = 0),
  list(id = "baseline, AEF dropped",         learners = names(LEARNERS), wt = FALSE, aef = "drop", screen = 0),
  list(id = "baseline + screen top-30",      learners = names(LEARNERS), wt = FALSE, aef = "all",  screen = 30)
)

rows <- list()
for (ctry in names(cfgs)) {
  cc <- cfgs[[ctry]]; hc <- H[H$country == ctry, ]
  for (ocn in names(cc$outcomes)) {
    svy <- rd(paste0("svy_admin2_", tolower(ctry), "_", ocn))
    if (is.null(svy) || nrow(svy) < 10) next
    tr <- dplyr::inner_join(svy[, c("Admin2", "svy_prev", "n_svy")], hc, by = "Admin2")
    tr <- tr[is.finite(tr$svy_prev) & !is.na(tr$Admin1), , drop = FALSE]
    if (nrow(tr) < 12 || dplyr::n_distinct(tr$Admin1) < 3) next
    rel <- admin2_reliability(svy, deff = 1.5, boot = 0); rmax <- rel$r_max %||% NA_real_
    y  <- stats::qlogis(pmin(pmax(tr$svy_prev, 0.005), 0.995))
    # precision weight on the LOGIT scale (delta method), normalised to mean 1
    vp <- area_sampling_var_shrunk(tr$svy_prev, tr$n_svy)
    p  <- pmin(pmax(tr$svy_prev, 0.005), 0.995)
    vl <- vp / (p * (1 - p))^2
    wprec <- 1 / pmax(vl, .Machine$double.eps); wprec <- wprec / mean(wprec)
    regions <- unique(tr$Admin1)
    message(sprintf("%s / %s  (n=%d, %d regions, r_max=%.3f)", ctry, ocn,
                    nrow(tr), length(regions), rmax))

    for (cfg in CONFIGS) {
      oof <- rep(NA_real_, nrow(tr))
      for (r in regions) {
        i <- which(tr$Admin1 == r)
        if (length(i) == nrow(tr)) next
        covs <- COVS
        if (cfg$screen > 0) {
          # Screen on the ADMIN-1 target built from TRAINING regions only.
          a1 <- tr[-i, ] |> group_by(Admin1) |>
            summarise(p1 = weighted.mean(svy_prev, pmax(n_svy, 1)), .groups = "drop")
          A <- aggregate_covariates_to_admin1(tr[-i, c("Admin1", "Admin2", COVS)])
          A <- dplyr::inner_join(a1, A, by = "Admin1")
          rr <- vapply(COVS, function(v) {
            x <- suppressWarnings(as.numeric(A[[v]]))
            if (sum(is.finite(x)) < 4 || stats::sd(x, na.rm = TRUE) == 0) return(0)
            abs(suppressWarnings(stats::cor(x, A$p1, use = "complete.obs")))
          }, numeric(1))
          rr[!is.finite(rr)] <- 0
          covs <- names(sort(rr, decreasing = TRUE))[seq_len(min(cfg$screen, sum(rr > 0)))]
          if (length(covs) < 2) covs <- COVS
        }
        btr <- .build_X(tr[-i, ], covs, cfg$aef)
        bte <- .build_X(tr[i, , drop = FALSE], covs, cfg$aef, pcs = btr$pcs)
        common <- intersect(colnames(btr$X), colnames(bte$X))
        s <- .fit_stack(btr$X[, common, drop = FALSE], y[-i],
                        if (cfg$wt) wprec[-i] else rep(1, nrow(tr) - length(i)),
                        cfg$learners)
        oof[i] <- .pred_stack(s, bte$X[, common, drop = FALSE])
      }
      pr <- stats::plogis(oof)
      ok <- is.finite(pr) & is.finite(tr$svy_prev)
      if (sum(ok) < 5) next
      r <- suppressWarnings(stats::cor(tr$svy_prev[ok], pr[ok]))
      rows[[length(rows) + 1L]] <- data.frame(
        country = ctry, outcome = ocn, config = cfg$id, n = sum(ok),
        pearson_r = round(r, 3), r_max = round(rmax, 3),
        r_share = if (is.finite(rmax) && rmax > 0) round(r / rmax, 2) else NA_real_,
        mae_pp = round(100 * mean(abs(tr$svy_prev[ok] - pr[ok])), 2),
        bias_pp = round(100 * mean(pr[ok] - tr$svy_prev[ok]), 2),
        stringsAsFactors = FALSE)
    }
  }
}

out <- dplyr::bind_rows(rows)
readr::write_csv(out, here("results", "tables", "estimator_ablation.csv"))
message("\n-> results/tables/estimator_ablation.csv (", nrow(out), " rows)")

base <- out %>% filter(config == "baseline (current stack)") %>%
  select(country, outcome, r0 = pearson_r, mae0 = mae_pp)
summ <- out %>% left_join(base, by = c("country", "outcome")) %>%
  group_by(config) %>%
  summarise(cells = dplyr::n(),
            median_r = round(stats::median(pearson_r, na.rm = TRUE), 3),
            median_r_share = round(stats::median(r_share, na.rm = TRUE), 2),
            median_mae_pp = round(stats::median(mae_pp, na.rm = TRUE), 2),
            beats_baseline_mae = sum(mae_pp < mae0, na.rm = TRUE),
            beats_baseline_r = sum(pearson_r > r0, na.rm = TRUE),
            .groups = "drop") %>% arrange(median_mae_pp)
print(as.data.frame(summ), row.names = FALSE)
