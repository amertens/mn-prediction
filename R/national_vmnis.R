# =============================================================================
# R/national_vmnis.R
#
# THE NATIONAL LEVEL, AND WHAT IT IS WORTH TO A DISTRICT MAP.
#
# Promoted 2026-08-31 from sandbox_parsimony/R/30-33 and /36. Four pieces:
#
#   fit_national_loco()      leave-one-COUNTRY-out prevalence prediction from
#                            the WHO VMNIS panel. The unit of observation is the
#                            country-survey, so a ~70-country panel gives an
#                            order of magnitude more design points for the
#                            BETWEEN-country relationship than aggregating the
#                            four-country Admin-2 model ever could.
#   predict_national_level() the same model pointed at a country the panel does
#                            not cover -- the deployment case.
#   compose_level_pattern()  the composition this exists for.
#   national_noise_ceiling() variance-component ceiling for the panel.
#
# WHY COMPOSE. The subnational work established that the proxy covariates carry
# a district PATTERN but not a LEVEL, which is why every production map is
# anchored to a design-based survey total (R/area_level_comparison.R; anchoring
# lifts r from 0.164 to 0.413 and cuts absolute bias from 3.24 to 1.59 pp). A
# country with no survey has no such anchor, and its transported map is
# therefore a ranking, not a set of prevalences. The VMNIS model supplies a
# level for exactly that case: predict the national prevalence from national
# covariates, then shift the transported pattern onto it. Level from VMNIS,
# pattern from the subnational model.
#
# WHAT IS TESTABLE. VMNIS carries Folate, Vitamin A, Vitamin B12, Vitamin D and
# Zinc -- it has NO iron or ferritin panel -- and the LOCO transport predictions
# cover iron and vitamin A. The intersection is vitamin A, so the composition is
# scored on 2 outcomes x 4 countries. That is a best case rather than a
# representative one: Vitamin A / preschool children is also the strongest VMNIS
# panel (69 countries). Read it as an upper bound on what the composition buys,
# not as an average over nutrients.
#
# THE CEILING, AND A CORRECTION. Two earlier ceilings for this panel (0.66 and
# 0.505) differed only in whether repeat surveys were grouped by country-year or
# by country-year-METHOD -- an arbitrary choice deciding what counts as error.
# national_noise_ceiling() estimates the variance components instead, using
# VMNIS's recorded Samplesize for the sampling term. For Vitamin A / preschool
# children (528 surveys, 70 countries) it gives r_max = 0.818 for what a survey
# will REPORT and 0.866 for a standardised prevalence; the best LOCO model
# reaches r = 0.655, so 80% of the ceiling. There is real headroom here, unlike
# at Admin-2 where the honest arms are near-saturated -- but see
# scripts/covariates/19_national_composition.R for why closing it would not
# help a district map.
# =============================================================================

.nat_logit  <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.nat_ilogit <- stats::plogis

# Deficiency outcome -> VMNIS (micronutrient, population) panel. Iron is absent
# from VMNIS by design, so it has no entry rather than a wrong one.
NAT_PANEL_MAP <- list(
  child_vitA   = c("Vitamin A",   "Preschool-age children"),
  child_zinc   = c("Zinc",        "Preschool-age children"),
  women_vitA   = c("Vitamin A",   "Non-pregnant women (NPW)"),
  women_folate = c("Folate",      "Non-pregnant women (NPW)"),
  women_b12    = c("Vitamin B12", "Non-pregnant women (NPW)")
)

national_panel_for <- function(outcome) NAT_PANEL_MAP[[outcome]]

vmnis_national <- function(path = here::here("data", "national",
                                             "vmnis_national_rep.rds")) {
  if (!file.exists(path)) stop("VMNIS national panel not found at: ", path)
  nat <- readRDS(path)
  nat$iso3c <- suppressWarnings(countrycode::countrycode(
    nat$country, "country.name", "iso3c", warn = FALSE))
  nat$pop <- as.character(nat$Population)
  nat[!is.na(nat$iso3c), , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Leave-one-country-out national prevalence.
#
# Every survey from the held-out country is removed, so the model has never
# seen that country's prevalence at any date -- the situation the product is
# for. Returns the scored metrics AND the per-country-year predictions; the
# sandbox version returned only metrics, and the composition needs the
# predictions.
#
# Surveys are never dropped for missing covariates. With up to ~1,600 candidate
# columns a complete-case rule would discard nearly all of them, and "this
# country never reported this indicator" is signal in a DHS-derived panel, not
# noise to be papered over with a median. Missingness indicators, kNN
# imputation and screening all happen INSIDE the fold, on training rows only.
# ---------------------------------------------------------------------------
fit_national_loco <- function(d, vars, log_vars, label, source = "wdi",
                              screen_k = 150L, knn_k = 5L, min_countries = 12L) {
  d <- d[is.finite(d$prev), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  .m <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  d <- d |>
    dplyr::group_by(iso3c, year) |>
    dplyr::summarise(prev = mean(prev), n_svy = .m(Samplesize),
                     dplyr::across(dplyr::all_of(vars), .m), .groups = "drop") |>
    as.data.frame()
  ctys <- unique(d$iso3c)
  if (length(ctys) < min_countries) return(NULL)

  X <- as.matrix(d[, vars, drop = FALSE]); storage.mode(X) <- "double"
  X[!is.finite(X)] <- NA_real_
  keep <- colSums(!is.na(X)) > 0L &
    apply(X, 2, function(x) length(unique(x[!is.na(x)])) > 1L)
  X <- X[, keep, drop = FALSE]
  for (v in intersect(log_vars, colnames(X))) X[, v] <- log1p(pmax(X[, v], 0))
  y <- .nat_logit(d$prev)

  preds <- matrix(NA_real_, nrow(d), 3,
                  dimnames = list(NULL, c("null", "ridge", "rf")))
  n_used <- n_mi <- integer(0)
  for (ho in ctys) {
    te <- which(d$iso3c == ho); tr <- setdiff(seq_len(nrow(d)), te)
    if (length(tr) < 10) next
    preds[te, "null"] <- .nat_ilogit(mean(y[tr]))

    MI <- missing_indicators(X, tr)
    Xi <- knn_impute_fold(X, tr, k = knn_k)
    Xf <- if (is.null(MI)) Xi else cbind(Xi, MI)
    n_mi <- c(n_mi, if (is.null(MI)) 0L else ncol(MI))

    sel <- prescreen_cols(Xf[tr, , drop = FALSE], y[tr], screen_k)
    Xtr <- Xf[tr, sel, drop = FALSE]; Xte <- Xf[te, sel, drop = FALSE]
    n_used <- c(n_used, ncol(Xtr))

    cvm <- tryCatch(glmnet::cv.glmnet(Xtr, y[tr], alpha = 0, nfolds = 5),
                    error = function(e) NULL)
    if (!is.null(cvm))
      preds[te, "ridge"] <- .nat_ilogit(as.numeric(
        stats::predict(cvm, newx = Xte, s = "lambda.min")))
    rf <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = y[tr], Xtr),
                                  num.trees = 800, min.node.size = 5, seed = 1),
                   error = function(e) NULL)
    if (!is.null(rf))
      preds[te, "rf"] <- .nat_ilogit(stats::predict(
        rf, data = data.frame(Xte))$predictions)
  }

  metrics <- dplyr::bind_rows(lapply(colnames(preds), function(k) {
    p <- preds[, k]; ok <- is.finite(p)
    if (sum(ok) < 10) return(NULL)
    data.frame(panel = label, source = source, model = k,
               n_cov_pool = ncol(X),
               n_mi_added = if (length(n_mi)) round(mean(n_mi)) else NA_integer_,
               n_cov_used = if (length(n_used)) round(mean(n_used)) else NA_integer_,
               n_surveys = sum(ok),
               n_countries = dplyr::n_distinct(d$iso3c[ok]),
               mean_prev_pp = round(100 * mean(d$prev[ok]), 1),
               mae_pp  = round(100 * mean(abs(p[ok] - d$prev[ok])), 2),
               rmse_pp = round(100 * sqrt(mean((p[ok] - d$prev[ok])^2)), 2),
               pearson = round(suppressWarnings(stats::cor(p[ok], d$prev[ok])), 3),
               spearman = round(suppressWarnings(stats::cor(
                 p[ok], d$prev[ok], method = "spearman")), 3),
               stringsAsFactors = FALSE)
  }))

  list(metrics = metrics,
       predictions = data.frame(panel = label, iso3c = d$iso3c, year = d$year,
                                observed = d$prev, n_svy = d$n_svy,
                                null = preds[, "null"], ridge = preds[, "ridge"],
                                rf = preds[, "rf"], stringsAsFactors = FALSE))
}


# ---------------------------------------------------------------------------
# Predict a national level for a country the panel need not cover.
#
# fit_national_loco() can only score countries VMNIS holds. The deployment case
# is the opposite one: a country with no national survey of its own. This trains
# on the whole panel minus that country -- so the exclusion stays honest even
# where VMNIS does happen to hold it -- and predicts at the requested year.
# ---------------------------------------------------------------------------
predict_national_level <- function(d, vars, log_vars, iso3c, year,
                                   cov_df, screen_k = 150L, knn_k = 5L) {
  .m <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  tr <- d[is.finite(d$prev) & d$iso3c != iso3c, , drop = FALSE] |>
    dplyr::group_by(iso3c, year) |>
    dplyr::summarise(prev = mean(prev),
                     dplyr::across(dplyr::all_of(vars), .m), .groups = "drop") |>
    as.data.frame()
  if (nrow(tr) < 20) return(NA_real_)

  # The target row is the country's covariates at the survey year; if that exact
  # year is absent, take the nearest year the panel does hold.
  cd <- cov_df[cov_df$iso3c == iso3c, , drop = FALSE]
  if (!nrow(cd)) return(NA_real_)
  te <- cd[which.min(abs(cd$year - year)), , drop = FALSE]

  X <- rbind(as.matrix(tr[, vars, drop = FALSE]),
             as.matrix(te[, vars, drop = FALSE]))
  storage.mode(X) <- "double"; X[!is.finite(X)] <- NA_real_
  itr <- seq_len(nrow(tr)); ite <- nrow(X)
  keep <- colSums(!is.na(X[itr, , drop = FALSE])) > 0L &
    apply(X[itr, , drop = FALSE], 2, function(x) length(unique(x[!is.na(x)])) > 1L)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(NA_real_)
  for (v in intersect(log_vars, colnames(X))) X[, v] <- log1p(pmax(X[, v], 0))
  y <- .nat_logit(tr$prev)

  MI <- missing_indicators(X, itr)
  Xi <- knn_impute_fold(X, itr, k = knn_k)
  Xf <- if (is.null(MI)) Xi else cbind(Xi, MI)
  sel <- prescreen_cols(Xf[itr, , drop = FALSE], y, screen_k)
  cvm <- tryCatch(glmnet::cv.glmnet(Xf[itr, sel, drop = FALSE], y,
                                    alpha = 0, nfolds = 5),
                  error = function(e) NULL)
  if (is.null(cvm)) return(NA_real_)
  .nat_ilogit(as.numeric(stats::predict(
    cvm, newx = Xf[ite, sel, drop = FALSE], s = "lambda.min")))
}


# ---------------------------------------------------------------------------
# THE COMPOSITION: level from VMNIS, pattern from the subnational model.
#
# The shift is a single constant on the LOGIT scale, solved so the
# sample-size-weighted mean of the shifted district predictions equals the
# supplied level. A monotone shift leaves the district RANKING untouched by
# construction, so Spearman is identical across arms -- and that is the point:
# the arms differ only in level, which is the quantity under test.
# ---------------------------------------------------------------------------
shift_to_level <- function(p, target, w = NULL) {
  ok <- is.finite(p)
  if (!any(ok) || !is.finite(target)) return(p)
  if (is.null(w)) w <- rep(1, length(p))
  w[!is.finite(w)] <- 0
  if (sum(w[ok]) <= 0) w <- rep(1, length(p))
  z <- .nat_logit(p)
  f <- function(dz) sum(w[ok] * .nat_ilogit(z[ok] + dz)) / sum(w[ok]) - target
  if (f(-25) > 0 || f(25) < 0) return(p)
  dz <- stats::uniroot(f, c(-25, 25), tol = 1e-8)$root
  .nat_ilogit(z + dz)
}

compose_level_pattern <- function(pattern_df, level, weight_col = "n_svy") {
  w <- if (weight_col %in% names(pattern_df)) pattern_df[[weight_col]] else NULL
  pattern_df$composed <- shift_to_level(pattern_df$modeled_prev, level, w)
  pattern_df
}


# ---------------------------------------------------------------------------
# Variance-component noise ceiling for a VMNIS panel.
#
#   logit(prev)_s = country_c + method_m + year trend + residual_s + sampling_s
#
# r_max_report is the ceiling for predicting what a survey will REPORT, with the
# method unknown in advance; r_max_standardised is the ceiling for a country's
# prevalence at one fixed method. The second is higher precisely because method
# variance has moved out of the error term -- which is what outcome
# standardisation actually buys, stated as a number rather than as a side
# effect of how repeat surveys were grouped.
# ---------------------------------------------------------------------------
# WS6a, 2026-08-31. THE SAMPLING TERM WAS WRONG BY A FACTOR OF 4.7, AND THE
# "DEGENERATE lmer FIT" WAS NOT A DEGENERATE lmer FIT.
#
# The published version reported sd_sampling 0.816 for the Vitamin A preschool
# panel. On the logit scale at a prevalence near 25 percent that implies an
# effective sample size of about THIRTEEN, and national nutrition surveys do not
# have thirteen respondents. Three causes, each measured
# (results/tables/vmnis_sampling_audit.csv):
#
#   1. v_bar was the ARITHMETIC MEAN of a reciprocal. The mean of 1/n is set by
#      the smallest surveys, not the typical one. That panel holds 34 surveys
#      under n = 50 and a minimum of n = 8, against a median of 373.5. Using the
#      median of the same per-survey variances gives 0.172 instead of 0.816.
#   2. The prevalence clamp sat at 0.005, so p(1-p) could fall to 0.005 and a
#      single near-zero-prevalence survey could contribute up to 301/(n-1). That
#      panel holds 37 surveys with prevalence under 2 percent.
#   3. Rows with no Samplesize dropped out of v_bar but stayed in the lmer fit,
#      so the variance components and the sampling term described different sets
#      of surveys.
#
# CONSEQUENCE FOR THE RECORD. Section 6.3 of the review document flags sd_resid
# at exactly 0.000 in two panels and attributes it to a degenerate lmer fit,
# concluding r_max is untrustworthy there. That diagnosis was wrong. The raw
# residual variance is non-zero in all four panels; the reported zero was the
# over-large sampling term being subtracted from a healthy residual and floored
# at zero. With the corrected term no panel sits at a boundary and all four are
# usable, so the refit with weakly informative priors that Section 12 recommends
# is not needed.
#
# WHAT CHANGES IN THE OUTPUT. Vitamin A preschool moves from r_max_report 0.818
# to 0.869 and sd_resid from 0.000 to 0.703, so the model's saturation against
# the ceiling falls from about 80 percent to about 75 percent. There is MORE
# headroom than the published table reports, not less. Three diagnostic columns
# are added (resid_at_boundary, sampling_exceeds_resid, usable); every previously
# reported column keeps its name.
#
# RE-RUNNING scripts/covariates/19_national_composition.R WILL THEREFORE CHANGE
# results/tables/national_vmnis_ceiling.csv. The pre-change values are preserved
# in results/tables/frozen_2026-09/national_vmnis_ceiling.csv and the
# side-by-side comparison is in results/tables/national_vmnis_ceiling_revised.csv.
#
# @param centre "median" (default) or "mean" for the per-survey sampling
#   variances. "mean" reproduces the published behaviour and is kept so the
#   difference stays inspectable rather than becoming folklore.
# @param clamp prevalence clamp; 0.02 by default, 0.005 in the published version
# @param require_n drop rows with no recorded Samplesize before fitting, so the
#   variance components and the sampling term describe the same surveys
national_noise_ceiling <- function(d, deff = 1.5, clamp = 0.02,
                                   centre = c("median", "mean"),
                                   require_n = TRUE) {
  centre <- match.arg(centre)
  need <- c("prev", "iso3c", "year", "Samplesize")
  if (!all(need %in% names(d))) return(NULL)
  d <- d[is.finite(d$prev) & !is.na(d$iso3c), , drop = FALSE]
  d$.n <- suppressWarnings(as.numeric(d$Samplesize))
  if (require_n) d <- d[is.finite(d$.n), , drop = FALSE]
  if (nrow(d) < 30 || dplyr::n_distinct(d$iso3c) < 10) return(NULL)
  d$y <- .nat_logit(d$prev)
  if (!"method" %in% names(d)) d$method <- "unspecified"
  d$yr_c <- d$year - mean(d$year, na.rm = TRUE)

  fm <- tryCatch(lme4::lmer(y ~ yr_c + (1 | iso3c) + (1 | method), data = d,
                            REML = TRUE), error = function(e) NULL)
  if (is.null(fm)) return(NULL)
  vc <- as.data.frame(lme4::VarCorr(fm))
  gv <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v[1] else 0 }
  s2_country <- gv("iso3c"); s2_method <- gv("method")
  s2_resid   <- vc$vcov[vc$grp == "Residual"][1]

  # Binomial sampling variance on the logit scale, design-effect inflated. The
  # residual already contains it, so it is subtracted out rather than added.
  p <- pmin(pmax(d$prev, clamp), 1 - clamp)
  v_s <- deff / (pmax(d$.n, 2) - 1) / (p * (1 - p))
  v_s <- v_s[is.finite(v_s)]
  if (!length(v_s)) return(NULL)
  v_bar <- if (centre == "median") stats::median(v_s) else mean(v_s)

  # Flags rather than a silent floor: a reader should be able to see when the
  # subtraction, not the fit, produced a zero.
  resid_boundary <- s2_resid <= 1e-8
  sampling_exceeds_resid <- v_bar > s2_resid
  s2_resid_adj <- max(s2_resid - v_bar, 0)

  lam_raw <- s2_country / (s2_country + s2_method + s2_resid_adj + v_bar)
  lam_std <- s2_country / (s2_country + s2_resid_adj + v_bar)
  data.frame(surveys = nrow(d), countries = dplyr::n_distinct(d$iso3c),
             sd_country = round(sqrt(s2_country), 3),
             sd_method  = round(sqrt(s2_method), 3),
             sd_resid   = round(sqrt(s2_resid_adj), 3),
             sd_sampling = round(sqrt(v_bar), 3),
             r_max_report = round(sqrt(lam_raw), 3),
             r_max_standardised = round(sqrt(lam_std), 3),
             resid_at_boundary = resid_boundary,
             sampling_exceeds_resid = sampling_exceeds_resid,
             usable = !resid_boundary & !sampling_exceeds_resid,
             stringsAsFactors = FALSE)
}
