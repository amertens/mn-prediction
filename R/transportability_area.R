# =============================================================================
# R/transportability_area.R
#
# Universal, parsimonious AREA-LEVEL (Admin-2) transportability models for
# mapping micronutrient deficiency from NON-SURVEY proxy predictors only
# (GEE remote sensing + IHME + Malaria Atlas + food security).
#
# Selected recipe (see scripts/transportability_experiments.R for the model
# bake-off): elastic net on the top-30 correlation-screened predictors, fit on
# the logit-prevalence scale and weighted by survey sample size (~8 predictors
# retained). On the full ~1.2k-covariate set with the harmonized outcome,
# heavy correlation screening (top-15) overfits and within-country centering
# does not help; un-screened ridge is marginally better but uses all
# predictors, so the screened elastic net is the parsimonious default
# (AREA_TRANSPORT_RECIPE_RIDGE holds the max-performance alternative).
#
# Outcome comparability: the Admin-2 survey prevalence is derived from a uniform
# WHO cutoff on the adjusted biomarker (compute_svy_admin2 / UNIFORM_TRANSPORT_
# TAGS), so every country uses the SAME definition (ferritin<12/15 = ID,
# RBP<0.70 = VAD) rather than the heterogeneous survey ID/IDA binaries.
#
# Key empirical finding: under the comparable outcome, child iron and child
# vitamin A spatial patterns transfer modestly across borders (within-country
# Spearman ~0.30); women's iron does not (≈0, a data limitation).
# =============================================================================

# Default recommended recipe ------------------------------------------------
# Re-selected (scripts/transportability_experiments.R) AFTER outcome
# harmonization on the full covariate set (~1.2k common predictors). On that
# set, screening to 15 overfits and within-country centering no longer helps;
# the parsimonious sweet spot is elastic net on the top-30 correlation-screened
# predictors (~8 retained), at ~83% of the best (un-screened ridge) recipe.
# For maximum performance over parsimony, use ridge with screen_K = NULL.
AREA_TRANSPORT_RECIPE <- list(
  model    = "enet",         # elastic net (glmnet alpha = 0.5)
  alpha    = 0.5,
  screen_K = 30,             # correlation prescreen: keep top-K candidate predictors
  center   = FALSE,          # centering did not help on the harmonized/full set
  weight   = TRUE,           # weight Admin-2 units by survey n
  scale    = "logit",
  lambda   = "lambda.min",   # lambda.1se over-shrinks small folds to the null
  seed     = 12345           # fixes cv.glmnet folds for reproducible maps
)

# Maximum-performance (less parsimonious) alternative: ridge on all predictors.
AREA_TRANSPORT_RECIPE_RIDGE <- list(
  model = "ridge", alpha = 0, screen_K = NULL, center = FALSE,
  weight = TRUE, scale = "logit", lambda = "lambda.min", seed = 12345
)

# Domains treated as universal, non-survey proxies (present in all countries).
# Pricing (wfp_) and chirps/soil/worldpop/ntl are NOT universal (e.g. Sierra
# Leone lacks WFP price data) and are therefore excluded from the universal
# model; they can be explored as country-subset sensitivities.
AREA_TRANSPORT_DOMAINS <- "^(ihme_|MAP_|fsec_)"

.tr_logit  <- function(p) stats::qlogis(pmin(pmax(p, 0.005), 0.995))
.tr_ilogit <- stats::plogis

.tr_as_num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  if (is.factor(v)) return(suppressWarnings(as.numeric(as.character(v))))
  suppressWarnings(as.numeric(v))
}

#' Harmonize year-suffixed GEE variable names into common stems
#' Strips YYYY and YYYYMMDD tokens so multi-year layers collapse to one stem
#' (this lifts the cross-country common-GEE set from ~13 to ~158 variables).
harmonize_gee_names <- function(x) {
  x <- gsub("_(19|20)[0-9]{2}([0-9]{4})?", "", x)
  x <- gsub("__+", "_", x)
  sub("_$", "", x)
}

#' Build enriched Admin-2 covariate table for one country
#'
#' GEE comes from the clean Admin-2 polygon extraction (year-harmonized and
#' averaged across years); IHME / MAP / food-security come from the merged
#' cluster-level table aggregated to Admin-2 means.
#'
#' @param gee_admin2 data.frame with Admin2 + gee_* columns
#' @param merged     per-country merged cluster-level data.frame (has Admin2)
#' @return data.frame keyed on Admin2 with harmonized proxy covariates
build_admin2_covariates <- function(gee_admin2, merged) {
  gee_cols <- grep("^gee_", names(gee_admin2), value = TRUE)
  stems <- harmonize_gee_names(gee_cols)
  gmat <- sapply(gee_cols, function(c) .tr_as_num(gee_admin2[[c]]))
  gee_df <- as.data.frame(sapply(unique(stems), function(s) {
    cols <- gee_cols[stems == s]
    if (length(cols) == 1) gmat[, cols] else rowMeans(gmat[, cols, drop = FALSE], na.rm = TRUE)
  }))
  names(gee_df) <- unique(stems)
  gee_df$Admin2 <- as.character(gee_admin2$Admin2)

  dom_cols <- grep(AREA_TRANSPORT_DOMAINS, names(merged), value = TRUE)
  agg <- data.frame(Admin2 = as.character(merged$Admin2))
  for (c in dom_cols) agg[[c]] <- .tr_as_num(merged[[c]])
  agg <- agg |>
    dplyr::group_by(Admin2) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    as.data.frame()

  dplyr::full_join(gee_df, agg, by = "Admin2")
}

#' Assemble the pooled multi-country Admin-2 dataset for one outcome
#'
#' @param svy_admin2_list named list of svy_admin2 data.frames (Admin2, svy_prev, n_svy)
#' @param cov_list        named list of build_admin2_covariates() outputs
#' @param outcome         outcome tag (for labelling)
#' @return list(data, predictors, countries, outcome) or NULL if <2 countries
assemble_area_transport <- function(svy_admin2_list, cov_list, outcome = NA) {
  frames <- list()
  for (cn in names(svy_admin2_list)) {
    s <- svy_admin2_list[[cn]]
    if (is.null(s) || !any(is.finite(s$svy_prev))) next
    sv <- data.frame(
      Admin2   = as.character(s$Admin2),
      svy_prev = s$svy_prev,
      n_svy    = if ("n_svy" %in% names(s)) s$n_svy else NA_real_,
      stringsAsFactors = FALSE
    )
    sv <- sv[is.finite(sv$svy_prev), , drop = FALSE]
    d <- dplyr::inner_join(sv, cov_list[[cn]], by = "Admin2")
    if (nrow(d) == 0) next
    d$country <- cn
    frames[[cn]] <- d
  }
  if (length(frames) < 2) return(NULL)
  covnames <- function(df) {
    cv <- setdiff(names(df), c("Admin2", "svy_prev", "n_svy", "country"))
    cv[sapply(cv, function(c) any(is.finite(df[[c]])))]
  }
  common <- Reduce(intersect, lapply(frames, covnames))
  # Optional covariate hygiene, applied after the intersection so every country
  # is pruned identically. Active by default since WS2b; set
  # GEE_COVARIATE_HYGIENE=false to keep the unpruned set.
  common <- setdiff(common, prune_gee_covariates(common))
  pooled <- dplyr::bind_rows(lapply(frames, function(df)
    df[, c("country", "Admin2", "svy_prev", "n_svy", common), drop = FALSE]))
  list(data = pooled, predictors = common,
       countries = names(frames), outcome = outcome)
}

# --- internal: matrix prep, screening, centering, fitting ------------------
.tr_prep_X <- function(train, test, preds) {
  Xtr <- as.matrix(train[, preds, drop = FALSE])
  Xte <- as.matrix(test[,  preds, drop = FALSE])
  med <- apply(Xtr, 2, function(z) stats::median(z[is.finite(z)]))
  med[!is.finite(med)] <- 0
  for (j in seq_along(preds)) {
    Xtr[!is.finite(Xtr[, j]), j] <- med[j]
    Xte[!is.finite(Xte[, j]), j] <- med[j]
  }
  sdv <- apply(Xtr, 2, stats::sd); keep <- is.finite(sdv) & sdv > 1e-8
  list(Xtr = Xtr[, keep, drop = FALSE], Xte = Xte[, keep, drop = FALSE],
       preds = preds[keep])
}
.tr_screen <- function(Xtr, ytr, K) {
  r <- abs(suppressWarnings(apply(Xtr, 2, function(z) stats::cor(z, ytr))))
  r[!is.finite(r)] <- 0
  sort(order(r, decreasing = TRUE)[seq_len(min(K, ncol(Xtr)))])
}
.tr_center_by <- function(X, grp) {
  for (g in unique(grp)) {
    idx <- grp == g
    X[idx, ] <- sweep(X[idx, , drop = FALSE], 2, colMeans(X[idx, , drop = FALSE]), "-")
  }
  X
}
.tr_fit_predict <- function(Xtr, ytr, Xte, recipe, w = NULL) {
  if (recipe$model == "enet" || recipe$model == "lasso" || recipe$model == "ridge") {
    a <- switch(recipe$model, enet = recipe$alpha %||% 0.5, lasso = 1, ridge = 0)
    s <- recipe$lambda %||% "lambda.1se"
    if (ncol(Xtr) < 2) return(list(pred = rep(mean(ytr), nrow(Xte)), vars = colnames(Xtr)))
    if (!is.null(recipe$seed)) set.seed(recipe$seed)   # reproducible CV folds
    cv <- tryCatch(glmnet::cv.glmnet(Xtr, ytr, alpha = a, weights = w, nfolds = 5),
                   error = function(e) NULL)
    if (is.null(cv)) return(list(pred = rep(mean(ytr), nrow(Xte)), vars = character()))
    b <- as.matrix(stats::coef(cv, s = s))
    sel <- rownames(b)[-1][b[-1, 1] != 0]
    list(pred = as.numeric(stats::predict(cv, newx = Xte, s = s)), vars = sel)
  } else if (recipe$model == "rf") {
    df <- data.frame(y = ytr, Xtr)
    rf <- ranger::ranger(y ~ ., data = df, num.trees = 800, min.node.size = 5,
                         case.weights = w, seed = recipe$seed %||% 1L)
    list(pred = stats::predict(rf, data.frame(Xte))$predictions, vars = colnames(Xtr))
  } else stop("unknown model")
}
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
# ^ NA-aware; kept byte-identical to the canonical definition in
#   R/corrected/00_corrected_utils.R so behaviour is independent of source order.

#' Run leave-one-country-out area-level transportability CV
#'
#' For each held-out country, trains the recipe on the remaining countries and
#' predicts the held-out country's Admin-2 prevalences (out-of-sample). Returns
#' both per-country metrics and per-Admin-2 predictions (for difference maps).
#'
#' @param pooled output of assemble_area_transport()
#' @param recipe list like AREA_TRANSPORT_RECIPE
#' @return list(metrics, predictions, selected_vars)
run_area_transport_loco <- function(pooled, recipe = AREA_TRANSPORT_RECIPE) {
  d <- pooled$data; preds0 <- pooled$predictors
  metrics <- list(); preds_out <- list(); selvars <- list()
  for (ho in pooled$countries) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 10 || nrow(te) < 4) next
    ytr <- .tr_logit(tr$svy_prev); yte <- te$svy_prev
    w <- if (isTRUE(recipe$weight)) pmax(tr$n_svy, 1, na.rm = TRUE) else NULL

    pp <- .tr_prep_X(tr, te, preds0); Xtr <- pp$Xtr; Xte <- pp$Xte
    ytr_fit <- ytr; level <- mean(ytr)
    if (isTRUE(recipe$center)) {
      # (Issue 3) Center the held-out country on the TRAINING column means, never
      # on its own (test) means: `colMeans(Xte)` would leak held-out information
      # into its own features. NB group-centering train + pooled-mean centering
      # test is an approximation; the principled transportable analogue is
      # fit_predict_two_stage() in R/benchmark_models.R.
      tr_means <- colMeans(Xtr)                        # pooled train means (pre-centering)
      Xtr <- .tr_center_by(Xtr, tr$country)
      Xte <- sweep(Xte, 2, tr_means, "-")
      for (g in unique(tr$country)) {
        idx <- tr$country == g; ytr_fit[idx] <- ytr[idx] - mean(ytr[idx])
      }
    }
    if (!is.null(recipe$screen_K) && is.finite(recipe$screen_K)) {
      sel <- .tr_screen(Xtr, ytr_fit, recipe$screen_K)
      Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
    }
    fp <- .tr_fit_predict(Xtr, ytr_fit, Xte, recipe, w)
    plogit <- fp$pred
    if (isTRUE(recipe$center)) plogit <- plogit + level
    pred <- pmin(pmax(.tr_ilogit(plogit), 0), 1)

    ok <- is.finite(pred) & is.finite(yte)
    if (sum(ok) < 4) next
    selvars[[ho]] <- fp$vars
    preds_out[[ho]] <- data.frame(
      outcome = pooled$outcome, country = ho, Admin2 = te$Admin2,
      n_svy = te$n_svy, survey_prev = yte, modeled_prev = pred,
      diff_pp = (pred - yte) * 100, stringsAsFactors = FALSE)
    pr <- suppressWarnings(stats::cor(yte[ok], pred[ok]))
    sr <- suppressWarnings(stats::cor(yte[ok], pred[ok], method = "spearman"))
    cal <- tryCatch(stats::coef(stats::lm(yte[ok] ~ pred[ok]))[2], error = function(e) NA)
    # (Issue 6) area-bootstrap CI + within-held-out-country permutation null on the
    # transport correlation, so a reader can tell r ~ 0.30 from sampling noise.
    cin <- if (exists("metric_ci_null"))
      metric_ci_null(pred[ok], yte[ok], "pearson", seed = 101L)
    else data.frame(pearson_ci_lo = NA_real_, pearson_ci_hi = NA_real_,
                    pearson_perm_p = NA_real_, n_boot = 0L)
    metrics[[ho]] <- data.frame(
      outcome = pooled$outcome, held_out = ho, n_train = nrow(tr), n_test = sum(ok),
      n_pred = ncol(Xtr), n_selected = length(fp$vars),
      pearson_r = round(pr, 3),
      pearson_ci_lo = cin$pearson_ci_lo, pearson_ci_hi = cin$pearson_ci_hi,
      pearson_perm_p = cin$pearson_perm_p,
      spearman_r = round(sr, 3),
      rmse_pp = round(sqrt(mean((yte[ok] - pred[ok])^2)) * 100, 2),
      mae_pp  = round(mean(abs(yte[ok] - pred[ok])) * 100, 2),
      nat_bias_pp = round((mean(pred[ok]) - mean(yte[ok])) * 100, 2),
      calib_slope = round(as.numeric(cal), 2), n_boot = cin$n_boot,
      stringsAsFactors = FALSE)
  }
  list(
    metrics      = if (length(metrics)) dplyr::bind_rows(metrics) else NULL,
    predictions  = if (length(preds_out)) dplyr::bind_rows(preds_out) else NULL,
    selected_vars = selvars
  )
}
