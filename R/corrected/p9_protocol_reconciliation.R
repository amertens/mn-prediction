# =============================================================================
# R/corrected/p9_protocol_reconciliation.R
#
# Puts the individual-level SuperLearner (aggregated to Admin-2) and the
# area-level model on ONE protocol so their district-ranking skill (Admin-2
# Pearson r vs survey-weighted prevalence) is directly comparable.
#
# Motivation: the previously reported individual-model r (~0.53) and area-model
# r (~0.12) are NOT comparable — they differ in fold scheme (cluster- vs
# region-blocked), aggregation weighting (unweighted vs survey-weighted), and
# model class. The dominant confound is the fold scheme: cluster-blocked folds
# do NOT hold out whole districts, so district-level signal leaks across folds
# and the aggregated individual predictions look artificially good. Under one
# honest protocol (region/spatial-block CV, survey-weighted target) the two are
# essentially tied; a principled unit-level small-area model (nested-error
# GLMM) is the preferred district-ranking estimator.
#
# This target REUSES the already-fitted corrected SL object (P1), so it does NOT
# refit the SuperLearner — it only re-aggregates its existing out-of-fold (OOF)
# predictions under each protocol and fits the cheap area-enet / SAE comparators.
#
# Columns returned (Admin-2 Pearson r unless noted):
#   indiv_optim_unwt    individual SL, random CV + full-data preproc, unweighted
#                       aggregation  (≈ the optimistic / production path)
#   indiv_cluster_unwt  individual SL, cluster-block CV, unweighted  (the old "0.53")
#   indiv_cluster_wt    individual SL, cluster-block CV, survey-weighted
#   indiv_region_unwt   individual SL, region/spatial-block CV, unweighted
#   indiv_region_wt     individual SL, region/spatial-block CV, survey-weighted (MATCHED)
#   area_enet_region    area elastic-net, region-block leave-one-region-out (≈ the old "0.12")
#   sae_glmm_region     unit-level SAE (nested-error GLMM), region-block (honest, RE=0 OOS)
#   sae_glmm_insample   unit-level SAE in-sample EB (ceiling where the district HAS data)
# =============================================================================

# weighted mean ignoring non-finite
.rc_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
# guarded Pearson correlation (NA if degenerate)
.rc_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok]))
}

# Aggregate an OOF frame (row, yhat_full, Admin2) to Admin-2 and return the merged
# district (pred, svy_prev) pairs. `w_ind` is the per-row survey weight.
.rc_agg_pairs <- function(oof, w_ind, sv, weighted) {
  d <- data.frame(Admin2 = oof$Admin2, yhat = oof$yhat_full, w = w_ind,
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$yhat) & !is.na(d$Admin2), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  if (weighted) {
    idx <- split(seq_len(nrow(d)), d$Admin2)
    ag  <- vapply(idx, function(i) .rc_wmean(d$yhat[i], d$w[i]), numeric(1))
  } else {
    ag  <- tapply(d$yhat, d$Admin2, mean)
  }
  merge(data.frame(Admin2 = names(ag), pred = as.numeric(ag),
                   stringsAsFactors = FALSE), sv, by = "Admin2")
}
.rc_agg_r <- function(oof, w_ind, sv, weighted) {
  m <- .rc_agg_pairs(oof, w_ind, sv, weighted)
  if (is.null(m) || nrow(m) == 0) return(NA_real_)
  .rc_cor(m$pred, m$svy_prev)
}

# (Issue 6) bootstrap CI + permutation null for a matched district correlation;
# returns a 1-row frame of <prefix>_ci_lo/_ci_hi/_perm_p (NA if helper absent).
.rc_ci <- function(pred, obs, prefix, seed) {
  cols <- paste0(prefix, c("_ci_lo", "_ci_hi", "_perm_p"))
  if (!exists("metric_ci_null")) {
    z <- data.frame(NA_real_, NA_real_, NA_real_); names(z) <- cols; return(z)
  }
  metric_ci_null(pred, obs, "pearson", seed = seed, prefix = prefix)[, cols, drop = FALSE]
}

#' Reconcile individual-vs-area district-ranking skill on one protocol.
#'
#' @param slfit  a fit_corrected_sl() object (P1) — reused, never refit.
#' @param outcome_data build_outcome_dataset() output (for weights + admin map).
#' @param svy_admin2 compute_svy_admin2() output (Admin2, svy_prev, n_svy).
#' @param gee_admin2 area covariate frame (Admin2 + gee_ columns).
#' @param cc,oc country / outcome config.
#' @return one-row data.frame of matched-protocol Pearson r's (see file header).
reconcile_protocols <- function(slfit, outcome_data, svy_admin2, gee_admin2,
                                cc, oc, sae_k = 8L) {
  if (is.null(slfit) || !is.null(slfit$error)) return(NULL)
  if (is.null(svy_admin2) || !"svy_prev" %in% colnames(svy_admin2)) return(NULL)
  if (is.null(outcome_data$data)) return(NULL)
  sv <- svy_admin2[is.finite(svy_admin2$svy_prev),
                   intersect(c("Admin2", "svy_prev", "n_svy"), colnames(svy_admin2)),
                   drop = FALSE]
  if (!"n_svy" %in% colnames(sv)) sv$n_svy <- 1

  a2c <- cc$admin2_col %||% "Admin2"
  a1c <- cc$admin1_col %||% "Admin1"
  d   <- outcome_data$data
  ybin <- oc$binary
  if (is.null(ybin) || !ybin %in% colnames(d)) return(NULL)
  y    <- as.numeric(d[[ybin]]); keep <- !is.na(y)        # SAME filter as P1
  wcol <- cc$weight_col
  w_ind <- if (!is.null(wcol) && wcol %in% colnames(d))
    suppressWarnings(as.numeric(d[[wcol]]))[keep] else rep(1, sum(keep))
  w_ind[!is.finite(w_ind) | w_ind <= 0] <- 1

  res <- data.frame(
    country = slfit$country, outcome = slfit$outcome,
    n_ind = slfit$n, n_area = nrow(sv),
    indiv_optim_unwt   = .rc_agg_r(slfit$optimistic$oof,     w_ind, sv, FALSE),
    indiv_cluster_unwt = .rc_agg_r(slfit$honest_cluster$oof, w_ind, sv, FALSE),
    indiv_cluster_wt   = .rc_agg_r(slfit$honest_cluster$oof, w_ind, sv, TRUE),
    indiv_region_unwt  = .rc_agg_r(slfit$honest_region$oof,  w_ind, sv, FALSE),
    indiv_region_wt    = .rc_agg_r(slfit$honest_region$oof,  w_ind, sv, TRUE),
    area_enet_region   = NA_real_, sae_glmm_region = NA_real_,
    sae_glmm_insample  = NA_real_, stringsAsFactors = FALSE)

  # (Issue 6) bootstrap CI + permutation null on the three matched headline
  # district correlations. indiv_region_wt is available now; area_enet_region and
  # sae_glmm_region are filled below once those comparators are fit.
  pr_rw <- .rc_agg_pairs(slfit$honest_region$oof, w_ind, sv, TRUE)
  res <- cbind(res,
    .rc_ci(if (is.null(pr_rw)) numeric(0) else pr_rw$pred,
           if (is.null(pr_rw)) numeric(0) else pr_rw$svy_prev,
           "indiv_region_wt", 301L),
    data.frame(area_enet_region_ci_lo = NA_real_, area_enet_region_ci_hi = NA_real_,
               area_enet_region_perm_p = NA_real_,
               sae_glmm_region_ci_lo = NA_real_, sae_glmm_region_ci_hi = NA_real_,
               sae_glmm_region_perm_p = NA_real_))

  # Admin2 -> Admin1 crosswalk (region blocking key for the area/SAE models)
  cw <- unique(data.frame(Admin2 = as.character(d[[a2c]]),
                          Admin1 = as.character(d[[a1c]]), stringsAsFactors = FALSE))
  cw <- cw[!is.na(cw$Admin2) & !is.na(cw$Admin1), , drop = FALSE]
  gcols <- if (!is.null(gee_admin2)) grep("^gee_", colnames(gee_admin2), value = TRUE) else character(0)
  if (length(gcols) < 2 || nrow(cw) == 0) return(res)

  area <- merge(merge(sv, cw, by = "Admin2"),
                gee_admin2[, c("Admin2", gcols)], by = "Admin2")
  area <- area[is.finite(area$svy_prev), , drop = FALSE]
  if (nrow(area) < 8 || length(unique(area$Admin1)) < 2) return(res)

  # ── Area elastic-net, leave-one-region-out ─────────────────────────────────
  if (requireNamespace("glmnet", quietly = TRUE)) {
    pa <- rep(NA_real_, nrow(area))
    X0 <- as.matrix(area[, gcols]); X0[!is.finite(X0)] <- NA
    for (rg in unique(area$Admin1)) {
      tri <- area$Admin1 != rg; tei <- area$Admin1 == rg
      if (sum(tri) < 6) next
      X <- X0
      mu <- apply(X[tri, , drop = FALSE], 2, function(z) stats::median(z, na.rm = TRUE))
      for (j in seq_along(gcols)) {
        z <- X[, j]; z[!is.finite(z)] <- if (is.finite(mu[j])) mu[j] else 0; X[, j] <- z
      }
      kv <- which(apply(X[tri, , drop = FALSE], 2, stats::sd) > 0)
      if (length(kv) < 1) next
      cvf <- tryCatch(glmnet::cv.glmnet(X[tri, kv, drop = FALSE], area$svy_prev[tri],
                        alpha = 0.5, weights = pmax(area$n_svy[tri], 1)),
                      error = function(e) NULL)
      if (!is.null(cvf))
        pa[tei] <- as.numeric(stats::predict(cvf, X[tei, kv, drop = FALSE], s = "lambda.min"))
    }
    res$area_enet_region <- .rc_cor(pa, area$svy_prev)
    ci <- .rc_ci(pa, area$svy_prev, "area_enet_region", 302L)
    res$area_enet_region_ci_lo  <- ci$area_enet_region_ci_lo
    res$area_enet_region_ci_hi  <- ci$area_enet_region_ci_hi
    res$area_enet_region_perm_p <- ci$area_enet_region_perm_p
  }

  # ── Unit-level SAE: nested-error GLMM y_i ~ area covars + (1|Admin2) ───────
  if (requireNamespace("lme4", quietly = TRUE)) {
    di <- data.frame(y = y[keep],
                     Admin2 = as.character(d[[a2c]])[keep],
                     Admin1 = as.character(d[[a1c]])[keep], stringsAsFactors = FALSE)
    di <- di[is.finite(di$y) & !is.na(di$Admin2) & !is.na(di$Admin1), , drop = FALSE]
    gee_std <- gee_admin2[, c("Admin2", gcols)]
    ctrl <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
    pick_top <- function(areas_sub) {            # screen top-K covars on a subset
      cc2 <- vapply(gcols, function(v) {
        r <- .rc_cor(areas_sub[[v]], areas_sub$svy_prev); if (is.na(r)) 0 else abs(r)
      }, numeric(1))
      names(sort(cc2, decreasing = TRUE))[seq_len(min(sae_k, sum(cc2 > 0)))]
    }
    fit_sae <- function(train_regions) {
      top <- pick_top(area[area$Admin1 %in% train_regions, , drop = FALSE])
      if (length(top) < 1) return(NULL)
      av <- gee_std[, c("Admin2", top)]
      for (v in top) av[[v]] <- as.numeric(scale(suppressWarnings(as.numeric(av[[v]]))))
      dd <- merge(di, av, by = "Admin2")
      dd <- dd[stats::complete.cases(dd[, top, drop = FALSE]), , drop = FALSE]
      list(formula = stats::as.formula(paste("y ~", paste(top, collapse = " + "), "+ (1|Admin2)")),
           data = dd, top = top)
    }
    # region-blocked (honest): predict held-out region via fixed effects (RE=0)
    pg <- numeric(0); pa2 <- character(0)
    for (rg in unique(di$Admin1)) {
      sp <- fit_sae(setdiff(unique(di$Admin1), rg)); if (is.null(sp)) next
      tr <- sp$data[sp$data$Admin1 != rg, , drop = FALSE]
      te <- sp$data[sp$data$Admin1 == rg, , drop = FALSE]
      if (nrow(tr) < 50 || nrow(te) < 1 || length(unique(tr$y)) < 2) next
      g <- tryCatch(lme4::glmer(sp$formula, data = tr, family = stats::binomial, control = ctrl),
                    error = function(e) NULL)
      if (is.null(g)) next
      pg  <- c(pg, stats::predict(g, newdata = te, type = "response",
                                  re.form = NA, allow.new.levels = TRUE))
      pa2 <- c(pa2, te$Admin2)
    }
    if (length(pg) > 0) {
      ag <- tapply(pg, pa2, mean, na.rm = TRUE)
      m  <- merge(data.frame(Admin2 = names(ag), pred = as.numeric(ag)), sv, by = "Admin2")
      res$sae_glmm_region <- .rc_cor(m$pred, m$svy_prev)
      ci <- .rc_ci(m$pred, m$svy_prev, "sae_glmm_region", 303L)
      res$sae_glmm_region_ci_lo  <- ci$sae_glmm_region_ci_lo
      res$sae_glmm_region_ci_hi  <- ci$sae_glmm_region_ci_hi
      res$sae_glmm_region_perm_p <- ci$sae_glmm_region_perm_p
    }
    # in-sample EB (ceiling): full fit, predictions include the area random effect
    sp0 <- fit_sae(unique(di$Admin1))
    if (!is.null(sp0) && length(unique(sp0$data$y)) >= 2) {
      g0 <- tryCatch(lme4::glmer(sp0$formula, data = sp0$data, family = stats::binomial, control = ctrl),
                     error = function(e) NULL)
      if (!is.null(g0)) {
        pin <- stats::predict(g0, type = "response")
        ag  <- tapply(pin, sp0$data$Admin2, mean, na.rm = TRUE)
        m   <- merge(data.frame(Admin2 = names(ag), pred = as.numeric(ag)), sv, by = "Admin2")
        res$sae_glmm_insample <- .rc_cor(m$pred, m$svy_prev)
      }
    }
  }
  res
}

#' Bundle per-slice reconciliations: bind, write CSV, and report medians.
build_protocol_reconciliation <- function(rows_list) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, rows_list)
  if (length(rows) == 0) return(NULL)
  per_slice <- do.call(rbind, rows)
  num <- setdiff(names(per_slice), c("country", "outcome", "n_ind", "n_area"))
  medians <- sapply(per_slice[num], stats::median, na.rm = TRUE)
  out_dir <- file.path("results", "tables", "corrected")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(per_slice, file.path(out_dir, "protocol_reconciliation.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(metric = names(medians), median_r = round(medians, 3)),
                   file.path(out_dir, "protocol_reconciliation_medians.csv"),
                   row.names = FALSE)
  list(per_slice = per_slice, medians = round(medians, 3))
}
