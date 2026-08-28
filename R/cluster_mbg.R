# =============================================================================
# R/cluster_mbg.R  —  WS5: point-referenced (cluster-level) spatial model
#
# WHY THIS EXISTS
# ---------------
# The strongest structure this project has found is spatial, not covariate.
# sandbox_parsimony/out/spatial_plus_soil_rescored.csv has the covariate-free
# smoother at mean LOCO Pearson 0.2861 against 0.1426 for honestly-selected soil
# covariates, and WS1 showed the covariate-based transport signal is close to
# chance once selection is nested. So it is worth asking what a spatial model
# does when it is allowed to work at the resolution the survey actually has,
# namely the GPS cluster, instead of at the Admin-2 polygon.
#
# WHAT THIS IS NOT
# ----------------
# The original plan gated this on a displacement-integrated Earth Engine
# extraction of covariates at cluster points. That was dropped, on evidence:
#
#   * sandbox_parsimony/FINDINGS.md section 12 measured the noise ceiling
#     by spatial level and found it FALLS going finer (Admin-1 0.59, Admin-2
#     0.31, cluster 0.22). Finer units are noisier faster than they multiply.
#   * results/cluster_vs_admin2_comparison.csv has cluster-level correlation
#     worse than Admin-2 in 6 of 10 comparable cells and MAE worse in 8 of 10.
#   * The spatial smoother needs COORDINATES, which already exist for all four
#     countries (70, 90, 60 and 105 distinct cluster locations), not cluster
#     COVARIATES.
#
# So the covariate arm here uses ADMIN-2 covariates joined to clusters. Those
# are constant within district (Ghana: 75 distinct values across 2,450 rows, a
# between-district variance share of 1.00, FINDINGS.md section 12). That is a
# real limitation of the covariate arm and it is recorded per row in the
# `covariates_district_constant` column rather than buried. It does not affect
# the spatial arms, which is the point of the comparison.
#
# THE EVALUATION
# --------------
# Leave-one-Admin-2-out. Every cluster in a district is held out together, the
# model is fitted on the rest, and the district is predicted. That is the
# unsurveyed-district scenario: the smoother has to extrapolate to a location
# where it has no data.
#
# The design makes this a demanding test. The median district has ONE cluster in
# Gambia, Ghana and Malawi (62 of 75 Ghana districts and 76 of 89 Malawi
# districts have exactly one), so holding out a district usually means holding
# out a single point. Only Sierra Leone has 2 to 8 clusters per district. This
# is the same fact that puts Admin-2 below the survey's design resolution.
#
# AGGREGATION
# -----------
# Predictions are made at held-out cluster coordinates and aggregated to the
# district weighted by cluster biomarker count. This is the plan's documented
# approximation: population-weighted aggregation over a fine grid needs the
# WorldPop rasters, which are a separate workstream. The approximation is
# labelled in `aggregation` on every row so it cannot be mistaken for a
# population-weighted figure.
# =============================================================================

CLUSTER_MBG_ARMS <- c("national_mean", "covariates_only",
                      "spatial_only", "spatial_plus_covariates", "matern_spamm")

#' Covariates for the covariate arms. Same fixed a-priori set as WS4, so the
#' spatial-versus-covariate contrast is not confounded with covariate count.
#'
#' A function rather than a constant: tar_source() loads R/ alphabetically, so
#' R/cluster_mbg.R is sourced before R/corrected/p12_distributional.R and a
#' top-level reference to DIST_COVARIATE_PATTERNS would not resolve. Deferring
#' to call time keeps one definition instead of two that can drift.
cluster_mbg_covariates <- function() DIST_COVARIATE_PATTERNS

#' Collapse individual records to one row per GPS cluster.
#'
#' @param merged per-country merged data.frame
#' @param cc,oc country and outcome configs
#' @param value continuous biomarker vector aligned to `merged`
#' @param threshold deficiency cut-point
#' @return data.frame with cluster, lon, lat, Admin2, n, k, p_w, or NULL
build_cluster_dataset <- function(merged, cc, oc, value, threshold) {
  d <- merged
  pop_col <- cc$child_flag
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    keep  <- !is.na(d[[pop_col]]) & d[[pop_col]] == oc$child_flag_val
    d     <- d[keep, , drop = FALSE]
    value <- value[keep]
  }
  lon <- suppressWarnings(as.numeric(d$lon))
  lat <- suppressWarnings(as.numeric(d$lat))
  w   <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  v   <- suppressWarnings(as.numeric(value))

  ok <- is.finite(v) & v > 0 & is.finite(lon) & is.finite(lat) &
        !is.na(d[[cc$cluster_id]]) & !is.na(d[[cc$admin2_col]])
  if (sum(ok) < 100) return(NULL)

  # Admin1 is carried so the covariate join can key on the PAIR: the area
  # covariate table is not deduplicated and its Admin-2 names are not unique.
  a1 <- if (!is.null(cc$admin1_col) && cc$admin1_col %in% colnames(d))
    as.character(d[[cc$admin1_col]])[ok] else NA_character_
  ind <- data.frame(
    cluster = as.character(d[[cc$cluster_id]])[ok],
    Admin1  = a1,
    Admin2  = as.character(d[[cc$admin2_col]])[ok],
    lon = lon[ok], lat = lat[ok], w = w[ok],
    def = as.integer(v[ok] < threshold),
    stringsAsFactors = FALSE
  )
  if (length(unique(ind$def)) < 2) return(NULL)

  cl <- ind |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      Admin1 = .data$Admin1[1],
      Admin2 = .data$Admin2[1],
      lon = stats::median(.data$lon), lat = stats::median(.data$lat),
      n = dplyr::n(),
      k = sum(.data$def),
      # Survey-weighted cluster prevalence. mgcv is given this as the response
      # with `weights = n`, rather than cbind(k, n - k), because that is how the
      # survey weight enters: cbind() counts cannot carry a non-integer weight.
      # With equal weights the two are identical.
      p_w = stats::weighted.mean(.data$def, .data$w),
      .groups = "drop"
    ) |>
    as.data.frame()
  cl <- cl[is.finite(cl$p_w) & cl$n > 0, , drop = FALSE]
  if (nrow(cl) < 20) return(NULL)
  cl
}

#' Basis dimension for the spatial smoother.
#'
#' With 60 to 105 clusters, mgcv's default of 30 basis functions for s(lon, lat)
#' is too rich. Capped at one third of the training clusters so the smoother
#' cannot interpolate the data it is fitted on.
.cm_k <- function(n_train) max(5L, min(30L, floor(n_train / 3)))

#' Fit one arm on training clusters and predict at test cluster coordinates.
#'
#' @return numeric vector of predicted prevalences, or NULL when the arm cannot
#'   be fitted (never a silent fallback to another arm)
fit_cluster_arm <- function(arm, train, test, covs, seed = 12345L) {
  set.seed(seed)
  wt <- train$n
  ok_cov <- length(covs) > 0 && all(covs %in% names(train))

  # Standardize covariates on the TRAINING clusters only.
  if (ok_cov) {
    mu <- vapply(covs, function(c) mean(train[[c]], na.rm = TRUE), numeric(1))
    sd <- vapply(covs, function(c) stats::sd(train[[c]], na.rm = TRUE), numeric(1))
    sd[!is.finite(sd) | sd == 0] <- 1
    for (c in covs) {
      train[[c]] <- (train[[c]] - mu[[c]]) / sd[[c]]
      test[[c]]  <- (test[[c]]  - mu[[c]]) / sd[[c]]
      train[[c]][!is.finite(train[[c]])] <- 0
      test[[c]][!is.finite(test[[c]])]   <- 0
    }
  }

  if (arm == "national_mean")
    return(rep(stats::weighted.mean(train$p_w, wt), nrow(test)))

  cov_term <- if (ok_cov) paste(covs, collapse = " + ") else NULL

  if (arm == "covariates_only") {
    if (!ok_cov) return(NULL)
    f <- stats::as.formula(paste("p_w ~", cov_term))
    m <- tryCatch(stats::glm(f, data = train, family = stats::binomial(), weights = wt),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    return(as.numeric(stats::predict(m, newdata = test, type = "response")))
  }

  if (arm %in% c("spatial_only", "spatial_plus_covariates")) {
    kk <- .cm_k(nrow(train))
    rhs <- sprintf("s(lon, lat, bs = 'tp', k = %d)", kk)
    if (arm == "spatial_plus_covariates") {
      if (!ok_cov) return(NULL)
      rhs <- paste(rhs, "+", cov_term)
    }
    f <- stats::as.formula(paste("p_w ~", rhs))
    m <- tryCatch(mgcv::gam(f, data = train, family = stats::binomial(),
                            weights = wt, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    p <- tryCatch(as.numeric(mgcv::predict.gam(m, newdata = test, type = "response")),
                  error = function(e) NULL)
    return(p)
  }

  if (arm == "matern_spamm") {
    # Optional engine. Never a hard dependency; the arm is simply absent when
    # spaMM is not installed, and that absence is recorded by the caller.
    if (!requireNamespace("spaMM", quietly = TRUE)) return(NULL)
    # Matern() must resolve UNQUALIFIED inside the formula. Writing
    # spaMM::Matern(1 | lon + lat) fails with "the condition has length > 1",
    # because spaMM inspects the formula terms rather than evaluating them.
    # Setting the formula's environment to spaMM's namespace makes Matern
    # visible without attaching the package.
    f <- stats::as.formula("cbind(k, n - k) ~ 1 + Matern(1 | lon + lat)")
    environment(f) <- asNamespace("spaMM")
    m <- tryCatch(spaMM::fitme(f, data = train, family = stats::binomial()),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    p <- tryCatch(as.numeric(stats::predict(m, newdata = test, type = "response")),
                  error = function(e) NULL)
    return(p)
  }

  NULL
}

#' Leave-one-Admin-2-out evaluation of every arm, for one country-outcome.
#'
#' @param cl build_cluster_dataset() output, with covariates already joined
#' @param svy compute_svy_admin2() output for this cell
#' @param covs covariate column names present on `cl`
#' @param arms which arms to run
#' @return data.frame, one row per arm, or NULL
run_cluster_mbg <- function(cl, svy, covs, cc, oc,
                            arms = CLUSTER_MBG_ARMS, seed = 12345L) {
  districts <- sort(unique(cl$Admin2))
  if (length(districts) < 8) return(NULL)

  preds <- matrix(NA_real_, nrow = nrow(cl), ncol = length(arms),
                  dimnames = list(NULL, arms))

  for (dd in districts) {
    te <- which(cl$Admin2 == dd)
    tr <- which(cl$Admin2 != dd)
    if (!length(te) || length(tr) < 20) next
    for (a in arms) {
      p <- fit_cluster_arm(a, cl[tr, , drop = FALSE], cl[te, , drop = FALSE],
                           covs, seed = seed)
      if (!is.null(p) && length(p) == length(te)) preds[te, a] <- p
    }
  }

  sv <- as.data.frame(svy); sv$Admin2 <- as.character(sv$Admin2)
  truth <- sv$svy_prev[match(districts, sv$Admin2)]
  n_svy <- sv$n_svy[match(districts, sv$Admin2)]
  evalable <- is.finite(truth) & is.finite(n_svy) & n_svy >= DIST_MIN_NSVY

  out <- list()
  for (a in arms) {
    # District prediction = cluster predictions weighted by cluster biomarker
    # count. Documented approximation to population weighting; see file header.
    ap <- vapply(districts, function(dd) {
      i <- cl$Admin2 == dd & is.finite(preds[, a])
      if (!any(i)) NA_real_ else stats::weighted.mean(preds[i, a], cl$n[i])
    }, numeric(1))

    ok <- evalable & is.finite(ap)
    if (sum(ok) < 5) next
    cal <- tryCatch(unname(stats::coef(stats::lm(truth[ok] ~ ap[ok]))[2]),
                    error = function(e) NA_real_)

    # Two reasons a correlation here can be an artifact rather than performance.
    #
    # 1. MECHANICAL. The national-mean arm predicts the leave-one-out training
    #    mean, which is anti-correlated with the held-out value BY CONSTRUCTION:
    #    a high-prevalence district drags the mean down when it is removed. Its
    #    correlation is a property of the CV scheme, not of the model, so it is
    #    never used as the correlation baseline. MAE and RMSE are the metrics the
    #    null comparison is read on.
    # 2. DEGENERATE. Any arm whose predictions are near-constant relative to the
    #    truth yields a correlation dominated by rounding. Flagged and set to NA.
    pred_sd  <- stats::sd(ap[ok])
    truth_sd <- stats::sd(truth[ok])
    degenerate <- !is.finite(pred_sd) || !is.finite(truth_sd) ||
                  truth_sd <= 0 || pred_sd < 0.02 * truth_sd
    r_p <- if (degenerate) NA_real_ else suppressWarnings(stats::cor(ap[ok], truth[ok]))
    r_s <- if (degenerate) NA_real_ else
      suppressWarnings(stats::cor(ap[ok], truth[ok], method = "spearman"))
    m <- data.frame(
      country = cc$country, outcome = oc$tag, arm = a,
      n_clusters = nrow(cl), n_districts = length(districts),
      n_districts_eval = sum(ok),
      median_clusters_per_district = stats::median(as.numeric(table(cl$Admin2))),
      districts_single_cluster = sum(as.numeric(table(cl$Admin2)) == 1),
      pearson_r  = r_p,
      spearman_r = r_s,
      pred_sd_over_truth_sd = if (is.finite(truth_sd) && truth_sd > 0) round(pred_sd / truth_sd, 4) else NA_real_,
      correlation_degenerate = degenerate,
      correlation_is_mechanical = identical(a, "national_mean"),
      mae_pp     = 100 * mean(abs(ap[ok] - truth[ok])),
      rmse_pp    = 100 * sqrt(mean((ap[ok] - truth[ok])^2)),
      bias_pp    = 100 * (mean(ap[ok]) - mean(truth[ok])),
      calib_slope = cal,
      aggregation = "cluster-count weighted (not population weighted)",
      covariates_district_constant = a %in% c("covariates_only", "spatial_plus_covariates"),
      seed = seed,
      stringsAsFactors = FALSE
    )
    m <- add_reliability_columns(m, sv[match(districts[ok], sv$Admin2), , drop = FALSE])
    out[[a]] <- m
  }
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out)
}
