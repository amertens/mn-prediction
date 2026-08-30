# =============================================================================
# sandbox_parsimony/R/20_recommended_recipe.R
#
# The two recipes the bake-offs point to, written so they can be lifted into
# R/ with minimal edits. Nothing here is wired into _targets.R.
#
#   area_model_v2()      within-country Admin-2 map from proxy covariates
#   area_transport_v2()  cross-country (LOCO) transport of the spatial pattern
#
# Both take the survey COUNTS, not a pre-computed proportion, because the whole
# point is to stop treating a district measured on 6 children the same as one
# measured on 200.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(glmnet)})

# ---------------------------------------------------------------------------
# 1. Within-country Admin-2 model
#
# Random forest on ~16 a-priori constructs PLUS longitude and latitude, fitted
# on the logit scale with survey-n case weights. Beat the production screened
# elastic net in 14 of 16 country x outcome cells (paired mean +0.063 Spearman,
# t = 4.9) using 18 predictors instead of 237-487.
#
# Changes relative to fit_area_sl() / the area-level elastic net:
#   * a parsimonious predictor set instead of a top-K correlation screen over
#     ~1,200 mostly-duplicated columns, which spends its budget on near-copies
#     of one or two underlying constructs;
#   * geography IN the model. A covariate-free lon/lat smoother already matches
#     the production recipe, so leaving coordinates out of both the model and
#     the benchmark table hides where the skill is coming from;
#   * a forest rather than a linear fit, so the covariate-by-geography
#     interactions are available without hand-specifying them.
#
# Things that were tried and did NOT help, so they are deliberately absent:
#   * a binomial likelihood on the survey counts (statistically the right model
#     for a proportion built from n Bernoulli draws, and it needs no clamping of
#     the 0%-prevalence districts) -- at 30-170 districts the extra parameters
#     cost more than the likelihood buys. See 08_spatial_models.R.
#   * a thin-plate spline field instead of raw coordinates -- same story.
#
# CAVEAT ON lon/lat. With random CV folds a district's neighbours are in the
# training set, so coordinates let the forest interpolate. That matches the
# real use case (predicting an unsurveyed district surrounded by surveyed ones)
# but flatters the model when a whole REGION is unsurveyed. Stress-test with
# leave-one-Admin1-out CV before promising accuracy for an unsurveyed region.
# ---------------------------------------------------------------------------

#' @param train data.frame with Admin2, svy_prev, n_svy, lon, lat + covariates
#' @param newdata districts to predict (all of them, surveyed or not)
#' @param vars    candidate covariate names
#' @param curated a-priori construct names; falls back to correlation-cluster
#'   representatives for whatever is missing
#' @param k_rep  target predictor count before coordinates are added
#' @param use_xy include lon/lat as predictors
#' @return list(pred, se, vars_used, fit)
area_model_v2 <- function(train, newdata, vars, curated = NULL, k_rep = 20L,
                          use_xy = TRUE, num_trees = 1000L, seed = 1L) {
  stopifnot(all(c("svy_prev", "n_svy") %in% names(train)))
  pp <- prep_X(train, newdata, vars)
  Xtr <- pp$Xtr; Xne <- pp$Xte; vv <- pp$vars

  sel <- if (!is.null(curated)) intersect(curated, vv) else character(0)
  if (length(sel) < k_rep)
    sel <- unique(c(sel, decorr_reps(Xtr, setdiff(vv, sel), k = k_rep - length(sel))))
  Xtr <- Xtr[, sel, drop = FALSE]; Xne <- Xne[, sel, drop = FALSE]

  if (isTRUE(use_xy) && all(is.finite(train$lon)) && all(is.finite(newdata$lon))) {
    Xtr <- cbind(Xtr, lon = train$lon, lat = train$lat)
    Xne <- cbind(Xne, lon = newdata$lon, lat = newdata$lat)
  }

  df <- data.frame(y = .logit(train$svy_prev), Xtr, check.names = TRUE)
  fit <- tryCatch(
    ranger::ranger(y ~ ., data = df, num.trees = num_trees, min.node.size = 5,
                   case.weights = pmax(train$n_svy, 1), keep.inbag = TRUE,
                   seed = seed),
    error = function(e) NULL)
  if (is.null(fit))
    return(list(pred = rep(mean(train$svy_prev), nrow(newdata)),
                se = rep(NA_real_, nrow(newdata)), vars_used = colnames(Xtr),
                fit = NULL))

  nd <- data.frame(Xne, check.names = TRUE)
  p  <- stats::predict(fit, data = nd)$predictions
  # Infinitesimal-jackknife SE on the logit scale, mapped back through the
  # logistic derivative. This is FOREST uncertainty only -- it excludes the
  # district sampling error entirely, so do not present it as a credible
  # interval for the underlying prevalence.
  # The infinitesimal jackknife cannot be calibrated on fewer than ~20 rows, so
  # skip it rather than return an uncalibrated (sometimes NaN) SE. In production
  # newdata is every district in the country, so this branch is for CV folds.
  se_link <- if (nrow(nd) > 20)
    tryCatch(stats::predict(fit, data = nd, type = "se", se.method = "infjack")$se,
             error = function(e) rep(NA_real_, nrow(nd)))
  else rep(NA_real_, nrow(nd))
  pred <- .ilogit(p)
  list(pred = pmin(pmax(pred, 0), 1),
       se = se_link * pred * (1 - pred),
       vars_used = colnames(Xtr), fit = fit, xy_used = isTRUE(use_xy))
}

# ---------------------------------------------------------------------------
# 2. Cross-country transport
#
# The premise: a cross-survey biomarker LEVEL offset does not transport (raw
# ferritin differs several-fold between these surveys) but the within-country
# spatial PATTERN partly does. So the model is trained on the within-country
# z-score of the transformed prevalence, which makes level and spread
# structurally unlearnable, and the held-out country level is supplied by an
# anchor -- normally its published national prevalence.
#
# Report the two error components separately: rank agreement for the pattern,
# and bias/RMSE for the level. Collapsing them into one Pearson r is what makes
# a 44 pp level miss look like r = 0.52.
# ---------------------------------------------------------------------------

#' @param train pooled Admin-2 rows from the training countries; needs `country`
#' @param test  the held-out country rows (covariates only are used)
#' @param national_anchor prevalence to centre the held-out country on. NULL
#'   falls back to the training-country mean, which is what production does and
#'   is exactly where the large level biases come from.
#' @return list(pred, zhat, vars_used)
area_transport_v2 <- function(train, test, vars, national_anchor = NULL,
                              alpha = 0, scale = c("raw", "logit")) {
  scale <- match.arg(scale)
  tf  <- if (scale == "logit") function(p) .logit(p) else identity
  itf <- if (scale == "logit") stats::plogis else identity

  pp <- prep_X(train, test, vars); Xtr <- pp$Xtr; Xte <- pp$Xte

  # centre AND scale covariates within country, including the held-out country
  # on its own means. Only its COVARIATES are used, and those are observable
  # everywhere without a survey -- that is the premise of the method, so it is
  # not leakage of the held-out outcome.
  for (g in unique(train$country)) {
    i <- train$country == g
    m <- colMeans(Xtr[i, , drop = FALSE])
    s <- apply(Xtr[i, , drop = FALSE], 2, stats::sd)
    s[!is.finite(s) | s < 1e-8] <- 1
    Xtr[i, ] <- sweep(sweep(Xtr[i, , drop = FALSE], 2, m, "-"), 2, s, "/")
  }
  m_te <- colMeans(Xte); s_te <- apply(Xte, 2, stats::sd)
  s_te[!is.finite(s_te) | s_te < 1e-8] <- 1
  Xte <- sweep(sweep(Xte, 2, m_te, "-"), 2, s_te, "/")

  y <- tf(train$svy_prev); yz <- y
  for (g in unique(train$country)) {
    i <- train$country == g
    s <- stats::sd(y[i]); if (!is.finite(s) || s < 1e-8) s <- 1
    yz[i] <- (y[i] - mean(y[i])) / s
  }

  set.seed(12345L)
  cv <- tryCatch(glmnet::cv.glmnet(Xtr, yz, alpha = alpha,
                                   weights = pmax(train$n_svy, 1), nfolds = 5),
                 error = function(e) NULL)
  if (is.null(cv))
    return(list(pred = rep(mean(train$svy_prev), nrow(test)), zhat = NULL,
                vars_used = character(0)))
  zhat <- as.numeric(stats::predict(cv, newx = Xte, s = "lambda.min"))

  lev <- if (!is.null(national_anchor)) tf(national_anchor) else mean(y)
  spread <- mean(vapply(unique(train$country),
                        function(g) stats::sd(y[train$country == g]), numeric(1)),
                 na.rm = TRUE)
  if (!is.finite(spread) || spread < 1e-8) spread <- 1
  sz <- stats::sd(zhat); if (!is.finite(sz) || sz < 1e-8) sz <- 1

  b <- as.matrix(stats::coef(cv, s = "lambda.min"))
  list(pred = pmin(pmax(itf(lev + spread * (zhat - mean(zhat)) / sz), 0), 1),
       zhat = zhat,
       vars_used = rownames(b)[-1][b[-1, 1] != 0])
}

# ---------------------------------------------------------------------------
# Demo: run both on child_iron and print what they give.
# ---------------------------------------------------------------------------
if (!isTRUE(get0(".RECIPE_FNS_ONLY", ifnotfound = FALSE))) {
  source("sandbox_parsimony/R/02_features.R")
  source("sandbox_parsimony/R/03_core.R")
  suppressPackageStartupMessages(library(ranger))
  pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")

  P <- pooled_all$child_iron
  d <- P$data[is.finite(P$data$svy_prev) & is.finite(P$data$n_svy), ]
  cur <- curated_vars(P$predictors)

  message("\n--- area_model_v2, within Ghana, 5-fold CV ---")
  gh <- d[d$country == "ghana", ]
  set.seed(1); fold <- sample(rep_len(1:5, nrow(gh)))
  oof <- rep(NA_real_, nrow(gh))
  for (f in 1:5) {
    te <- which(fold == f)
    m <- area_model_v2(gh[-te, ], gh[te, ], P$predictors, curated = cur)
    oof[te] <- m$pred
  }
  rel <- reliability(gh$svy_prev, gh$n_svy)
  message(sprintf("  r = %.3f, rho = %.3f, RMSE = %.1f pp | ceiling r_max = %.2f (%.0f%% of it)",
                  cor(oof, gh$svy_prev), cor(oof, gh$svy_prev, method = "spearman"),
                  sqrt(mean((oof - gh$svy_prev)^2)) * 100, rel$r_max,
                  100 * cor(oof, gh$svy_prev) / rel$r_max))

  message("\n--- area_transport_v2, LOCO, with and without a national anchor ---")
  for (ho in unique(d$country)) {
    tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
    if (nrow(tr) < 30 || nrow(te) < 8) next
    anchor <- stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1))  # published national figure
    a <- area_transport_v2(tr, te, P$predictors, national_anchor = anchor)
    b <- area_transport_v2(tr, te, P$predictors, national_anchor = NULL)
    message(sprintf("  %-12s rho = %+.3f | RMSE anchored %5.1f pp vs unanchored %5.1f pp | %d vars",
                    ho, cor(a$pred, te$svy_prev, method = "spearman"),
                    sqrt(mean((a$pred - te$svy_prev)^2)) * 100,
                    sqrt(mean((b$pred - te$svy_prev)^2)) * 100,
                    length(a$vars_used)))
  }
}
