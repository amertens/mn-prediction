# =============================================================================
# R/corrected/area_recipe.R  — area-level prevalence estimator ("area recipe")
#
# The deployable district-prevalence recipe arrived at by the methods-development
# sandbox (iterations 1-3 + 2b; see docs/AREA_LEVEL_RECIPE_SPEC.md). Models the
# AREA prevalence directly (not pseudo-replicated individual binaries), which is
# the dominant lever at effective n = number of districts.
#
# Recipe:
#   - target  : area survey-weighted prevalence; BINOMIAL (count) loss, weight n_svy
#   - features: regime-specific. mode="enriched" = gee_ (soil/climate/veg/...) +
#               MAP_ malaria + ihme_ anaemia/stunting  (within-country gap-filling);
#               mode="universal" = gee_ only (cross-country transport — MAP/ihme are
#               country-calibrated and hurt transport).
#   - select  : mild GLOBAL near-null pre-filter (drop bottom `prefilter_drop` by
#               marginal |cor|, validated leakage-free by the honest re-screen null)
#               -> in-fold supervised top-K screen -> penalised regression (enet).
#   - eval    : honest district-block AND region-block CV; metrics carry bootstrap
#               CI + permutation null (metric_ci_null).  Transport = leave-one-
#               country-out, scored by within-held-out Spearman (level-invariant).
#
# SAE note: do NOT Fay-Herriot-shrink the level toward the covariate mean — the
# covariates carry ranking, not level (sandbox 2b). Anchor the level to known
# national totals; use SAE/conformal for intervals only.
#
# All functions take the SAME objects the existing corrected targets produce
# (outcome_data, svy_admin2, gee_admin2, cc, oc) so wiring is a drop-in.
# Depends on 00_corrected_utils.R (%||%, make_group_folds, metric_ci_null).
# =============================================================================

# ---- regime covariate-domain patterns --------------------------------------
AR_UNIVERSAL_PREFIX <- "^gee_"
AR_ENRICH_PREFIX    <- "^MAP_|^map2|^map_|^ihme_"

#' Build the one-row-per-Admin2 area frame (svy_prev, n_svy, Admin1, covariates).
#' @param mode "enriched" (gee + MAP + IHME) or "universal" (gee only).
ar_build_frame <- function(outcome_data, svy_admin2, gee_admin2, cc, oc,
                           mode = c("enriched", "universal")) {
  mode <- match.arg(mode)
  if (is.null(svy_admin2) || !"svy_prev" %in% colnames(svy_admin2)) return(NULL)
  if (!"n_svy" %in% colnames(svy_admin2)) return(NULL)
  sv <- as.data.frame(svy_admin2); sv$Admin2 <- as.character(sv$Admin2)
  gee <- as.data.frame(gee_admin2); gee$Admin2 <- as.character(gee$Admin2)
  gcols <- grep(AR_UNIVERSAL_PREFIX, colnames(gee), value = TRUE)
  area <- merge(sv, gee[, c("Admin2", gcols)], by = "Admin2", sort = FALSE)

  # Admin1 crosswalk + (enriched) MAP/IHME aggregated to Admin2, from the
  # individual data the slice already loaded.
  a2c <- cc$admin2_col %||% "Admin2"; a1c <- cc$admin1_col %||% "Admin1"
  d <- outcome_data$data; d$.a2 <- as.character(d[[a2c]])
  cw <- unique(data.frame(Admin2 = d$.a2, Admin1 = as.character(d[[a1c]]),
                          stringsAsFactors = FALSE))
  cw <- cw[!is.na(cw$Admin2) & !is.na(cw$Admin1), ]; cw <- cw[!duplicated(cw$Admin2), ]
  area <- merge(area, cw, by = "Admin2", all.x = TRUE, sort = FALSE)

  covcols <- gcols
  if (mode == "enriched") {
    ex <- grep(AR_ENRICH_PREFIX, colnames(d), value = TRUE, ignore.case = TRUE)
    ex <- ex[vapply(ex, function(c) is.numeric(d[[c]]) || is.integer(d[[c]]), logical(1))]
    if (length(ex)) {
      agg <- stats::aggregate(d[, ex, drop = FALSE], list(Admin2 = d$.a2),
                              function(z) mean(z, na.rm = TRUE))
      area <- merge(area, agg, by = "Admin2", all.x = TRUE, sort = FALSE)
      covcols <- c(covcols, colnames(agg)[-1])
    }
  }
  area <- area[is.finite(area$svy_prev) & is.finite(area$n_svy) & area$n_svy >= 5, ]
  # keep covariates that are >50% finite across the surveyed areas
  covcols <- covcols[colSums(is.finite(as.matrix(area[, covcols, drop = FALSE]))) >
                       nrow(area) * 0.5]
  attr(area, "covcols") <- covcols
  attr(area, "mode") <- mode
  area
}

.ar_il  <- function(p) stats::qlogis(pmin(pmax(p, .005), .995))
.ar_imp <- function(M, med) { for (j in seq_len(ncol(M))) { z <- M[, j]; z[!is.finite(z)] <- med[j]; M[, j] <- z }; M }

# global near-null pre-filter: keep covariates whose |marginal cor with the target|
# is in the top (1 - drop) fraction. Honest (leakage-free) per the re-screen null.
ar_prefilter_target <- function(X, y, drop = 0.70) {
  r <- abs(suppressWarnings(as.numeric(stats::cor(X, y)))); r[!is.finite(r)] <- 0
  which(r >= stats::quantile(r, drop, na.rm = TRUE))
}
ar_prefilter <- function(X, prev, drop = 0.70) ar_prefilter_target(X, .ar_il(prev), drop)

# ---- (c) opt-in variants: domain map, late fusion, distributional ----------
.ar_domain <- function(cn) {
  d <- rep("other", length(cn))
  d[grepl("gee_soil", cn)] <- "soil"
  d[grepl("precip|chirps|tavg|lst|temp|tmax|tmin|prec|fldas", cn, ignore.case = TRUE)] <- "climate"
  d[grepl("ndvi|evi|fpar|lai|grassland|cropland|landcover|vegetation", cn, ignore.case = TRUE)] <- "veg"
  d[grepl("elevation|srtm|slope|terrain", cn, ignore.case = TRUE)] <- "terrain"
  d[grepl("pop|built|ghs|smod|accessib|nightlight|density", cn, ignore.case = TRUE)] <- "human"
  d[grepl("MAP_|map2|map_|pf_|pv_|parasite|malaria", cn, ignore.case = TRUE)] <- "malaria"
  d[grepl("ihme_", cn, ignore.case = TRUE)] <- "ihme"
  d[grepl("livestock|glw|cattle|gleam", cn, ignore.case = TRUE)] <- "livestock"
  d
}
.ar_enet_g <- function(Xtr, y, w, Xte) {
  f <- tryCatch(glmnet::cv.glmnet(Xtr, y, family = "gaussian", weights = w, nfolds = 5),
                error = function(e) NULL)
  if (is.null(f)) return(rep(mean(y), nrow(Xte)))
  as.numeric(stats::predict(f, Xte, s = "lambda.min"))
}
# late fusion: one regularised learner per domain (>=3 covariates), NNLS-stacked
# on the logit-prevalence scale. Wins on distributed signal; can lose sparse
# signal (see spec) -> opt-in per outcome.
.ar_late_fusion <- function(Xtr, prev_tr, w_tr, Xte, dom) {
  y <- .ar_il(prev_tr); doms <- unique(dom); doms <- doms[vapply(doms, function(z) sum(dom == z) >= 3, logical(1))]
  if (length(doms) < 2) return(.ar_enet_g(scale_safe(Xtr), y, w_tr, scale_safe(Xte, Xtr)) |> stats::plogis())
  inner <- make_random_folds(length(y), min(4L, length(y)), 7L)
  Z <- matrix(NA_real_, length(y), length(doms), dimnames = list(NULL, doms))
  for (ii in sort(unique(inner))) {
    itr <- inner != ii; ite <- inner == ii; if (sum(itr) < 8) next
    for (dd in doms) { dc <- which(dom == dd)
      Z[ite, dd] <- .ar_enet_g(scale_safe(Xtr[itr, dc, drop = FALSE]),
                               y[itr], w_tr[itr], scale_safe(Xtr[ite, dc, drop = FALSE], Xtr[itr, dc, drop = FALSE])) }
  }
  ok <- which(colSums(is.finite(Z)) == length(y)); if (!length(ok)) return(rep(mean(prev_tr), nrow(Xte)))
  Z <- Z[, ok, drop = FALSE]; doms <- doms[ok]
  cf <- tryCatch(nnls::nnls(Z * sqrt(w_tr), y * sqrt(w_tr))$x, error = function(e) rep(1, ncol(Z)))
  if (sum(cf) <= 0) cf <- rep(1, ncol(Z)); cf <- cf / sum(cf)
  Pte <- sapply(doms, function(dd) { dc <- which(dom == dd)
    .ar_enet_g(scale_safe(Xtr[, dc, drop = FALSE]), y, w_tr, scale_safe(Xte[, dc, drop = FALSE], Xtr[, dc, drop = FALSE])) })
  stats::plogis(as.numeric(Pte %*% cf))
}
# standardise Xte by Xtr (or self) column stats; robust to constants
scale_safe <- function(X, ref = X) {
  mu <- colMeans(ref); sdv <- apply(ref, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  m <- scale(X, mu, sdv); m[!is.finite(m)] <- 0; m
}
# per-area continuous biomarker mean + pooled within-area SD, on the cut-point scale
.ar_area_cont <- function(outcome_data, cc, oc) {
  d <- outcome_data$data; a2 <- as.character(d[[cc$admin2_col %||% "Admin2"]])
  cont <- suppressWarnings(as.numeric(d[[oc$continuous]]))
  logscale <- identical(oc$cutoff_scale, "log")
  ok <- is.finite(cont) & (if (logscale) cont > 0 else TRUE)
  val <- if (logscale) log(cont) else cont
  df <- data.frame(a2 = a2[ok], v = val[ok])
  mu <- tapply(df$v, df$a2, mean); resid <- df$v - mu[df$a2]
  list(mean = mu, sdpool = sqrt(mean(resid^2, na.rm = TRUE)),
       thr = if (logscale) log(oc$cutoff) else oc$cutoff)
}

# binomial enet on (prescreened, standardised) area X; weight = n_svy counts.
.ar_fit_predict <- function(Xtr, prev_tr, n_tr, Xte) {
  n <- pmax(round(n_tr), 1)
  ymat <- cbind(no = pmax(round((1 - prev_tr) * n), 0L), yes = pmax(round(prev_tr * n), 0L))
  f <- tryCatch(glmnet::cv.glmnet(Xtr, ymat, family = "binomial", nfolds = 5),
                error = function(e) NULL)
  if (is.null(f)) return(rep(mean(prev_tr), nrow(Xte)))
  as.numeric(stats::predict(f, Xte, s = "lambda.min", type = "response"))
}

# one honest CV pass: predict each held-out fold's area prevalence.
ar_cv_predict <- function(area, scheme = c("district", "region"),
                          drop = 0.70, screen_k = NULL, V = 10L, seed = 12345L,
                          estimator = c("binomial", "late_fusion", "distributional")) {
  scheme <- match.arg(scheme); estimator <- match.arg(estimator)
  covcols <- attr(area, "covcols"); prev <- area$svy_prev; w <- pmax(area$n_svy, 1)
  X0 <- as.matrix(area[, covcols, drop = FALSE]); dom <- .ar_domain(covcols)
  acont <- attr(area, "acont")              # distributional: list(mean, sdpool, thr)
  grp <- if (scheme == "district") area$Admin2 else area$Admin1
  Veff <- if (scheme == "district") V else min(V, length(unique(grp[!is.na(grp)])))
  folds <- make_group_folds(grp, Veff, seed)
  pred <- rep(NA_real_, nrow(area))
  for (f in sort(unique(folds))) {
    te <- which(folds == f); tr <- which(folds != f)
    if (length(tr) < 8) next
    med <- apply(X0[tr, , drop = FALSE], 2, stats::median, na.rm = TRUE); med[!is.finite(med)] <- 0
    Xtr <- .ar_imp(X0[tr, , drop = FALSE], med); Xte <- .ar_imp(X0[te, , drop = FALSE], med)
    k <- screen_k %||% max(4L, min(15L, floor(length(tr) / 4)))
    if (estimator == "late_fusion") {
      pred[te] <- .ar_late_fusion(Xtr, prev[tr], w[tr], Xte, dom)
    } else if (estimator == "distributional" && !is.null(acont)) {
      ymv <- as.numeric(acont$mean[match(area$Admin2, names(acont$mean))])
      ok <- tr[is.finite(ymv[tr])]; if (length(ok) < 8) next
      keep <- ar_prefilter_target(Xtr[match(ok, tr), , drop = FALSE], ymv[ok], drop)
      rr <- abs(suppressWarnings(as.numeric(stats::cor(Xtr[match(ok, tr), keep, drop = FALSE], ymv[ok]))))
      rr[!is.finite(rr)] <- 0; sel <- keep[order(rr, decreasing = TRUE)[seq_len(min(k, length(keep)))]]
      Ztr <- scale_safe(Xtr[match(ok, tr), sel, drop = FALSE])
      Zte <- scale_safe(Xte[, sel, drop = FALSE], Xtr[match(ok, tr), sel, drop = FALSE])
      pmu <- .ar_enet_g(Ztr, ymv[ok], w[ok], Zte)
      pred[te] <- stats::pnorm((acont$thr - pmu) / acont$sdpool)
    } else {                                # binomial (default)
      keep <- ar_prefilter(Xtr, prev[tr], drop)
      r <- abs(suppressWarnings(as.numeric(stats::cor(Xtr[, keep, drop = FALSE], .ar_il(prev[tr])))))
      r[!is.finite(r)] <- 0
      sel <- keep[order(r, decreasing = TRUE)[seq_len(min(k, length(keep)))]]
      Ztr <- scale_safe(Xtr[, sel, drop = FALSE]); Zte <- scale_safe(Xte[, sel, drop = FALSE], Xtr[, sel, drop = FALSE])
      pred[te] <- .ar_fit_predict(Ztr, prev[tr], w[tr], Zte)
    }
  }
  pred
}

#' Within-country area-recipe evaluation for one slice: district + region block,
#' each with rank metrics + bootstrap CI + permutation null. One row per scheme.
ar_within_country <- function(outcome_data, svy_admin2, gee_admin2, cc, oc,
                              mode = "enriched", drop = 0.70, seed = 12345L,
                              estimator = c("binomial", "late_fusion", "distributional")) {
  estimator <- match.arg(estimator)
  area <- ar_build_frame(outcome_data, svy_admin2, gee_admin2, cc, oc, mode)
  if (is.null(area) || nrow(area) < 8) return(NULL)
  # distributional (recommended for vitamin A): attach per-area continuous mean.
  if (estimator == "distributional") {
    if (is.null(oc$continuous) || !oc$continuous %in% colnames(outcome_data$data)) {
      estimator <- "binomial"           # no continuous biomarker -> fall back
    } else attr(area, "acont") <- .ar_area_cont(outcome_data, cc, oc)
  }
  out <- lapply(c("district", "region"), function(sch) {
    p <- tryCatch(ar_cv_predict(area, sch, drop, seed = seed, estimator = estimator),
                  error = function(e) NULL)
    if (is.null(p)) return(NULL)
    m <- data.frame(pred = p, obs = area$svy_prev)
    m <- m[is.finite(m$pred) & is.finite(m$obs), ]
    if (nrow(m) < 4 || stats::sd(m$pred) == 0) return(NULL)
    cin <- if (exists("metric_ci_null")) metric_ci_null(m$pred, m$obs, "pearson", seed = 501L)
           else data.frame(pearson_ci_lo = NA, pearson_ci_hi = NA, pearson_perm_p = NA, n_boot = 0L)
    slope <- tryCatch(stats::coef(stats::lm(obs ~ pred, m))[2], error = function(e) NA)
    data.frame(country = cc$country, outcome = oc$tag %||% oc$label, mode = mode,
               estimator = estimator, scheme = sch, n_area = nrow(m),
               pearson_r = round(stats::cor(m$pred, m$obs), 3),
               spearman_r = round(stats::cor(m$pred, m$obs, method = "spearman"), 3),
               calib_slope = round(as.numeric(slope), 2),
               mae_pp = round(mean(abs(m$pred - m$obs)) * 100, 2),
               pearson_ci_lo = cin$pearson_ci_lo, pearson_ci_hi = cin$pearson_ci_hi,
               pearson_perm_p = cin$pearson_perm_p, n_boot = cin$n_boot,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

#' Full-coverage district predictions for deployment (e.g. the dashboard): train
#' the recipe on the surveyed areas and predict EVERY admin-2 polygon, surveyed
#' or not. Uses the UNIVERSAL (gee, polygon-zonal) feature set, which is the only
#' one available at unsurveyed polygons (MAP/IHME are aggregated from survey
#' clusters). Returns one row per polygon with the predicted prevalence and, for
#' surveyed polygons, the direct survey estimate.
ar_full_coverage <- function(outcome_data, svy_admin2, gee_admin2, cc, oc, drop = 0.70) {
  train <- ar_build_frame(outcome_data, svy_admin2, gee_admin2, cc, oc, mode = "universal")
  if (is.null(train) || nrow(train) < 8) return(NULL)
  covcols <- attr(train, "covcols")
  gee <- as.data.frame(gee_admin2); gee$Admin2 <- as.character(gee$Admin2)
  poly <- gee[, c("Admin2", intersect(covcols, colnames(gee))), drop = FALSE]
  covcols <- intersect(covcols, colnames(poly))
  Xtr0 <- as.matrix(train[, covcols, drop = FALSE]); Xpe0 <- as.matrix(poly[, covcols, drop = FALSE])
  med <- apply(Xtr0, 2, stats::median, na.rm = TRUE); med[!is.finite(med)] <- 0
  Xtr <- .ar_imp(Xtr0, med); Xpe <- .ar_imp(Xpe0, med)
  keep <- ar_prefilter(Xtr, train$svy_prev, drop)
  Ztr <- scale_safe(Xtr[, keep, drop = FALSE]); Zpe <- scale_safe(Xpe[, keep, drop = FALSE], Xtr[, keep, drop = FALSE])
  pr <- .ar_fit_predict(Ztr, train$svy_prev, pmax(train$n_svy, 1), Zpe)
  out <- data.frame(country = cc$country, outcome = oc$tag %||% oc$label,
                    Admin2 = poly$Admin2, pred_prev = round(pr, 4),
                    surveyed = poly$Admin2 %in% train$Admin2, stringsAsFactors = FALSE)
  sv <- as.data.frame(svy_admin2); sv$Admin2 <- as.character(sv$Admin2)
  out <- merge(out, sv[, c("Admin2", "svy_prev", "n_svy")], by = "Admin2", all.x = TRUE, sort = FALSE)
  names(out)[names(out) == "svy_prev"] <- "direct_prev"
  out
}

#' Cross-country LOCO transport with the area recipe (universal covariates by
#' default). `frames` = named list of ar_build_frame() outputs (one per country),
#' all built on a COMMON covariate set. Scores within-held-out Spearman + null.
ar_transport_loco <- function(frames, drop = 0.70, seed = 7L) {
  frames <- Filter(function(a) !is.null(a) && nrow(a) >= 8, frames)
  if (length(frames) < 3) return(NULL)
  common <- Reduce(intersect, lapply(frames, function(a) attr(a, "covcols")))
  common <- common[vapply(common, function(c)
    all(vapply(frames, function(a) mean(is.finite(a[[c]])) > 0.5, logical(1))), logical(1))]
  if (length(common) < 5) return(NULL)
  rows <- lapply(names(frames), function(ho) {
    tr <- do.call(rbind, lapply(frames[setdiff(names(frames), ho)],
                                function(a) a[, c("svy_prev", "n_svy", common)]))
    te <- frames[[ho]]
    Xtr0 <- as.matrix(tr[, common, drop = FALSE]); Xte0 <- as.matrix(te[, common, drop = FALSE])
    med <- apply(Xtr0, 2, stats::median, na.rm = TRUE); med[!is.finite(med)] <- 0
    Xtr <- .ar_imp(Xtr0, med); Xte <- .ar_imp(Xte0, med)
    keep <- ar_prefilter(Xtr, tr$svy_prev, drop)
    mu <- colMeans(Xtr[, keep, drop = FALSE]); sdv <- apply(Xtr[, keep, drop = FALSE], 2, stats::sd)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    Ztr <- scale(Xtr[, keep, drop = FALSE], mu, sdv); Zte <- scale(Xte[, keep, drop = FALSE], mu, sdv)
    pr <- .ar_fit_predict(Ztr, tr$svy_prev, pmax(tr$n_svy, 1), Zte)
    m <- data.frame(pred = pr, obs = te$svy_prev); m <- m[is.finite(m$pred) & is.finite(m$obs), ]
    if (nrow(m) < 4 || stats::sd(m$pred) == 0) return(NULL)
    cin <- if (exists("metric_ci_null")) metric_ci_null(m$pred, m$obs, "spearman", seed = 601L)
           else data.frame(spearman_ci_lo = NA, spearman_ci_hi = NA, spearman_perm_p = NA, n_boot = 0L)
    data.frame(held_out = ho, n_area = nrow(m),
               spearman_r = round(stats::cor(m$pred, m$obs, method = "spearman"), 3),
               pearson_r = round(stats::cor(m$pred, m$obs), 3),
               spearman_ci_lo = cin$spearman_ci_lo, spearman_ci_hi = cin$spearman_ci_hi,
               spearman_perm_p = cin$spearman_perm_p, n_boot = cin$n_boot,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

#' Roll up per-slice area-recipe results, write CSVs, return the bundle.
#' @param within_list   list of ar_within_country() frames (one per slice)
#' @param transport_list named list (by outcome) of ar_transport_loco() frames
build_area_recipe_tables <- function(within_list, transport_list = list()) {
  od <- file.path("results", "tables", "corrected")
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  w <- do.call(rbind, Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0,
                             within_list))
  if (!is.null(w)) utils::write.csv(w, file.path(od, "area_recipe_within.csv"), row.names = FALSE)
  tr <- do.call(rbind, Filter(Negate(is.null), lapply(names(transport_list), function(nm) {
    x <- transport_list[[nm]]
    if (is.null(x) || !nrow(x)) NULL else cbind(outcome = nm, x, stringsAsFactors = FALSE)
  })))
  if (!is.null(tr)) utils::write.csv(tr, file.path(od, "area_recipe_transport.csv"), row.names = FALSE)
  message(sprintf("[area_recipe] wrote within (%s rows) + transport (%s rows)",
                  if (is.null(w)) 0 else nrow(w), if (is.null(tr)) 0 else nrow(tr)))
  list(within = w, transport = tr)
}
