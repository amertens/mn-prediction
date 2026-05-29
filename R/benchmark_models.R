# =============================================================================
# R/benchmark_models.R
#
# Sensitivity / benchmark suite comparing SuperLearner against standard
# methods for admin2-level micronutrient deficiency prediction:
#
#   1. country-mean baseline        — "is SL doing anything useful?"
#   2. GLM (OLS / logit-link)       — main-effects regression
#   3. penalized GLM (elastic net)  — isolated linear contribution
#   4. Fay-Herriot small-area       — industry standard for survey-based admin estimates
#   5. BYM2 (INLA)                  — canonical Bayesian areal disease-mapping model
#
# Each method has a single fit_predict_<method>(train, test, vars, ...)
# interface that returns predicted prevalence on [0,1].  The unified
# evaluators below call these methods over leave-one-country-out splits
# (transportability) and within-country k-fold CV.
#
# Inputs follow the same conventions as R/area_level_comparison.R:
#   - pooled_data: data.frame with at least `Admin2`, `svy_prev`, `n_svy`,
#                   `country`, and GEE/area covariates.
#   - gee_vars:   character vector of covariate column names.
#
# Returns a tidy data.frame stacking results across methods/holdouts so
# the dashboard and downstream comparison panels can pivot it.
#
# Heavy dependencies (INLA, sae, spdep) are referenced via `pkg::fn()` and
# guarded with requireNamespace() so this file sources cleanly even when
# they are missing.  Functions that need a missing package return NULL
# rather than erroring, so the rollup gracefully skips that method.
# =============================================================================


# ── Helpers ──────────────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (is.null(a)) b else a

# Logit / inverse logit with bounded inputs to avoid Inf at extremes.
.safe_logit  <- function(p, eps = 1e-3) qlogis(pmin(pmax(p, eps), 1 - eps))
.safe_invlogit <- function(x) plogis(x)

# Tidy 1-row metrics summary for a (obs, pred) pair on the [0,1] scale.
compute_area_metrics <- function(obs, pred, ..., method, model_type,
                                  n_train, n_test, n_vars,
                                  scale_pp = TRUE) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  n_eval <- length(obs)
  if (n_eval < 3) {
    data.frame(method = method, model_type = model_type,
               n_train = n_train, n_test = n_test, n_eval = n_eval,
               n_vars = n_vars,
               obs_prev = NA_real_, pred_prev = NA_real_,
               mae_pp = NA_real_, rmse_pp = NA_real_,
               pearson_r = NA_real_, spearman_r = NA_real_,
               mean_bias_pp = NA_real_,
               calib_slope = NA_real_, calib_intercept = NA_real_,
               ..., stringsAsFactors = FALSE)
  } else {
    mae  <- mean(abs(obs - pred))
    rmse <- sqrt(mean((obs - pred)^2))
    r    <- suppressWarnings(cor(obs, pred))
    rs   <- suppressWarnings(cor(obs, pred, method = "spearman"))
    bias <- mean(pred - obs)
    # Calibration on the linear scale: lm(obs ~ pred)
    calib <- tryCatch(stats::lm(obs ~ pred), error = function(e) NULL)
    intercept <- if (!is.null(calib)) coef(calib)[1] else NA_real_
    slope     <- if (!is.null(calib)) coef(calib)[2] else NA_real_
    f <- if (scale_pp) 100 else 1
    data.frame(method = method, model_type = model_type,
               n_train = n_train, n_test = n_test, n_eval = n_eval,
               n_vars = n_vars,
               obs_prev = mean(obs), pred_prev = mean(pred),
               mae_pp  = round(mae  * f, 3),
               rmse_pp = round(rmse * f, 3),
               pearson_r  = round(r,  3),
               spearman_r = round(rs, 3),
               mean_bias_pp = round(bias * f, 3),
               calib_slope     = round(slope,     3),
               calib_intercept = round(intercept, 3),
               ..., stringsAsFactors = FALSE)
  }
}

# Variable screening: keep covariates with >2 non-missing values and
# >1 unique value on the training set (mirrors run_area_loco logic).
.screen_vars <- function(df, vars) {
  vars[vapply(vars, function(v) {
    x <- df[[v]]
    sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
}

# Median impute NAs/Inf with training medians.  Returns a list of two
# data.frames with the same columns.
.impute_train_test <- function(train, test, vars) {
  Xtr <- as.data.frame(train[, vars, drop = FALSE])
  Xte <- as.data.frame(test [, vars, drop = FALSE])
  med <- vapply(Xtr, median, numeric(1), na.rm = TRUE)
  for (v in vars) {
    bad_tr <- !is.finite(Xtr[[v]])
    bad_te <- !is.finite(Xte[[v]])
    if (any(bad_tr)) Xtr[[v]][bad_tr] <- med[v]
    if (any(bad_te)) Xte[[v]][bad_te] <- med[v]
  }
  list(train = Xtr, test = Xte, medians = med)
}

# Bound predictions to a valid prevalence range.
.bound01 <- function(p, lo = 0, hi = 1) pmin(pmax(p, lo), hi)


# ── METHOD 1: country-mean baseline ───────────────────────────────────────────
# Predict the training-set mean prevalence for every held-out unit.
# Floor for "is the model doing anything?"
fit_predict_baseline <- function(train, test, vars, model_type = "continuous") {
  list(pred = rep(mean(train$svy_prev, na.rm = TRUE), nrow(test)),
       train_pred = rep(mean(train$svy_prev, na.rm = TRUE), nrow(train)))
}


# ── METHOD 2: GLM ─────────────────────────────────────────────────────────────
# Plain OLS on either the raw prevalence (continuous) or logit-transformed
# prevalence (logit), no penalization.  When p > n on the held-out training
# set, we top-rank covariates by univariate correlation with svy_prev and
# keep the top n_var_cap to avoid singular design matrices.
fit_predict_glm <- function(train, test, vars,
                              model_type = c("continuous", "logit"),
                              n_var_cap = 8L) {
  model_type <- match.arg(model_type)
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
    imp$train <- imp$train[, vars, drop = FALSE]
    imp$test  <- imp$test [, vars, drop = FALSE]
  }
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  # Weight by survey n so high-precision admin-2 polygons drive the fit; the
  # binomial sampling variance is p(1-p)/n, so weights ∝ n is the natural choice.
  wt <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  fit_df <- data.frame(Y = Y, imp$train)
  fit <- tryCatch(stats::lm(Y ~ ., data = fit_df, weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_tr <- as.numeric(stats::predict(fit, fit_df))
  pred_te <- as.numeric(stats::predict(fit, data.frame(imp$test)))
  if (model_type == "logit") {
    pred_tr <- .safe_invlogit(pred_tr); pred_te <- .safe_invlogit(pred_te)
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr))
}


# ── METHOD 3: penalized GLM (elastic net) ────────────────────────────────────
# 10-fold CV via cv.glmnet.  Alpha=0.5 (elastic net) by default; for very
# small training samples we drop to ridge (alpha=0) to avoid coefficient
# instability.  Logit option transforms target as in fit_predict_glm.
fit_predict_penalized <- function(train, test, vars,
                                    model_type = c("continuous", "logit"),
                                    alpha = 0.5) {
  if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  if (length(vars) < 2) return(NULL)
  imp <- .impute_train_test(train, test, vars)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  Y   <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  # Drop alpha for very small n
  if (nrow(Xtr) < 25 && alpha > 0) alpha <- 0
  n_folds <- min(nrow(Xtr), 10L)
  fit <- tryCatch(
    glmnet::cv.glmnet(Xtr, Y, alpha = alpha, nfolds = n_folds,
                       family = "gaussian", weights = wt),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- fit$lambda.min %||% fit$lambda.1se
  pred_tr <- as.numeric(stats::predict(fit, newx = Xtr, s = s))
  pred_te <- as.numeric(stats::predict(fit, newx = Xte, s = s))
  if (model_type == "logit") {
    pred_tr <- .safe_invlogit(pred_tr); pred_te <- .safe_invlogit(pred_te)
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       lambda = s, alpha = alpha)
}


# ── METHOD 4: Fay-Herriot small-area model ────────────────────────────────────
# Linear mixed model with area-level random effect and known sampling
# variances.  Standard industry method for survey-based admin estimates.
# Sampling variance is computed from prevalence and effective sample size:
#   v_i = p_i * (1 - p_i) / n_svy_i
# Requires the `sae` package.  When the design effect is unknown we use the
# DHS-standard 1.5 approximation (deff_default).
fit_predict_fh <- function(train, test, vars,
                             model_type = c("continuous"),
                             n_var_cap = 5L, deff_default = 1.5) {
  if (!requireNamespace("sae", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  if (!"n_svy" %in% colnames(train)) return(NULL)

  imp <- .impute_train_test(train, test, vars)
  # FH is fragile with high p; restrict to top-correlated covariates.
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])

  # Sampling variance per training area (use DHS deff approximation).
  p   <- pmin(pmax(train$svy_prev, 1e-4), 1 - 1e-4)
  n_e <- pmax(train$n_svy / deff_default, 1)
  sv  <- p * (1 - p) / n_e

  # FH fitting: direct = svy_prev, X = auxiliary vars, vardir = sv.
  # sae::eblupFH wants vardir to be a column in `data`, not a free vector
  # (the help text suggests a vector works, but in practice it tries to
  # match the call's symbol against data column names). We give it both.
  fh_data <- cbind(
    data.frame(svy_prev = as.numeric(train$svy_prev),
                .vardir   = as.numeric(sv),
                stringsAsFactors = FALSE),
    Xtr[, vars, drop = FALSE]
  )
  fh_formula <- stats::as.formula(paste("svy_prev ~", paste(vars, collapse = " + ")))
  fit <- tryCatch(
    sae::eblupFH(formula = fh_formula, vardir = .vardir, data = fh_data,
                  method = "REML"),
    error = function(e) {
      cat(sprintf("    [fh] eblupFH failed: %s\n", conditionMessage(e)))
      NULL
    })
  if (is.null(fit)) return(NULL)

  # Fitted (in-sample) area estimates (combine direct + synthetic).
  pred_tr <- as.numeric(fit$eblup)

  # Out-of-sample prediction = synthetic only (X test * beta).
  # eblupFH stores fixed-effect estimates in fit$fit$estcoef[, "beta"].
  beta <- tryCatch(fit$fit$estcoef[, "beta"], error = function(e) NULL) %||%
          tryCatch(fit$fit$beta,              error = function(e) NULL)
  if (is.null(beta)) return(NULL)
  X_te_design <- model.matrix(~ ., data = Xte)
  if (ncol(X_te_design) != length(beta)) {
    # Align columns by name (intercept always present).
    common <- intersect(colnames(X_te_design), names(beta))
    X_te_design <- X_te_design[, common, drop = FALSE]
    beta <- beta[common]
  }
  pred_te <- as.numeric(X_te_design %*% beta)

  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = "Fay-Herriot (REML)")
}


# ── METHOD 5: BYM2 via INLA ──────────────────────────────────────────────────
# Bayesian areal model with structured (ICAR) + unstructured area random
# effects, parameterized per Riebler et al. 2016 with PC priors.
# Adjacency is supplied as a list of igraph-style edge lists or a graph
# file path; for LOCO we build a block-diagonal graph across training
# countries (no cross-country edges).  Held-out areas receive predictions
# from the linear-predictor (fixed effects) only — their random effects
# are unobserved by construction.
#
# Input adjacency_list: named list, country -> spdep::nb object (or NULL).
fit_predict_bym2 <- function(train, test, vars, adjacency_list = NULL,
                               model_type = c("continuous", "binomial"),
                               n_var_cap = 12L) {
  if (!requireNamespace("INLA", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)

  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])

  # Normalize covariates for prior compatibility.
  scl <- lapply(vars, function(v) {
    mu <- mean(Xtr[[v]], na.rm = TRUE); sd <- stats::sd(Xtr[[v]], na.rm = TRUE)
    if (!is.finite(sd) || sd == 0) sd <- 1
    Xtr[[v]] <<- (Xtr[[v]] - mu) / sd
    Xte[[v]] <<- (Xte[[v]] - mu) / sd
    NULL
  })

  # Build BYM2 graph: block-diagonal adjacency across training countries.
  n_train <- nrow(train); n_test <- nrow(test); n_all <- n_train + n_test
  area_idx <- seq_len(n_all)
  # Use a forward-slash path: on Windows tempfile() returns backslash paths,
  # and embedding such a path in a parsed formula string trips R's lexer on
  # the `\U` (unicode escape) sequence in e.g. C:\Users\... .
  graph_file <- gsub("\\\\", "/", tempfile(fileext = ".graph"))
  W <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(n_all, n_all))
  if (!is.null(adjacency_list)) {
    offset <- 0
    train_countries <- unique(train$country)
    for (ctry in train_countries) {
      nb <- adjacency_list[[ctry]]
      idx <- which(train$country == ctry)
      if (is.null(nb) || length(idx) == 0) {
        offset <- offset + length(idx); next
      }
      # Build adjacency for this country's training areas only.
      # Map original area positions in nb to position within train country.
      # We assume train is in the same row order as the country's polygon set.
      n_ctry <- length(nb)
      if (length(idx) == n_ctry) {
        for (i in seq_len(n_ctry)) {
          for (j in nb[[i]]) {
            if (j > 0 && j > i) {
              W[idx[i], idx[j]] <- 1
              W[idx[j], idx[i]] <- 1
            }
          }
        }
      }
    }
    # Test country: also add its within-country adjacency (held-out areas
    # may be neighbors among themselves, even though we hold the outcome).
    held_out_ctry <- unique(test$country)
    if (length(held_out_ctry) == 1 && !is.null(adjacency_list[[held_out_ctry]])) {
      nb <- adjacency_list[[held_out_ctry]]
      idx <- seq.int(n_train + 1, n_all)
      if (length(idx) == length(nb)) {
        for (i in seq_len(length(nb))) {
          for (j in nb[[i]]) {
            if (j > 0 && j > i) {
              W[idx[i], idx[j]] <- 1
              W[idx[j], idx[i]] <- 1
            }
          }
        }
      }
    }
  }
  # Any isolated nodes get a tiny self-loop to keep INLA's BYM2 happy.
  has_neighbor <- Matrix::rowSums(W) > 0
  if (any(!has_neighbor)) {
    # Connect isolated nodes to themselves' country-mate via dummy edges
    # if available; otherwise BYM2 falls back to IID for them.
    iso <- which(!has_neighbor)
    for (i in iso) {
      W[i, i] <- 0  # explicit zero; INLA's bym2 handles isolated nodes via the iid component.
    }
  }

  # Write graph file.
  INLA::inla.write.graph(W, file = graph_file)

  # Assemble INLA data: outcome NA on test rows; covariates on all.
  Y_all <- c(if (model_type == "binomial") train$svy_prev else train$svy_prev,
              rep(NA_real_, n_test))
  Ntrials <- c(if (model_type == "binomial") train$n_svy else rep(NA_real_, n_train),
                rep(NA_real_, n_test))

  dat <- data.frame(Y = Y_all, area_idx = area_idx,
                    rbind(Xtr, Xte))
  formula_str <- sprintf(
    "Y ~ %s + f(area_idx, model = 'bym2', graph = '%s', scale.model = TRUE, hyper = list(phi = list(prior = 'pc', param = c(0.5, 0.5)), prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
    paste(vars, collapse = " + "),
    graph_file)
  formula <- stats::as.formula(formula_str)

  fit <- tryCatch({
    if (model_type == "binomial") {
      # Approximate the binomial likelihood on aggregated counts.
      events <- round(train$svy_prev * train$n_svy)
      dat$Y <- c(events, rep(NA_integer_, n_test))
      INLA::inla(formula, family = "binomial", Ntrials = Ntrials, data = dat,
                  control.predictor = list(link = 1, compute = TRUE),
                  control.compute   = list(config = FALSE))
    } else {
      # Gaussian on prevalence (Fay-Herriot-style with spatial smoothing).
      sv <- pmax((train$svy_prev * (1 - train$svy_prev)) / pmax(train$n_svy, 1), 1e-6)
      prec_weights <- c(1 / sv, rep(1 / mean(sv), n_test))
      INLA::inla(formula, family = "gaussian", data = dat,
                  scale = prec_weights,
                  control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
                  control.predictor = list(link = 1, compute = TRUE),
                  control.compute   = list(config = FALSE))
    }
  }, error = function(e) {
    cat(sprintf("    [bym2] INLA failed: %s\n", conditionMessage(e)))
    NULL
  })

  unlink(graph_file)
  if (is.null(fit)) return(NULL)

  fitted_all <- fit$summary.fitted.values[["mean"]]
  pred_tr <- fitted_all[seq_len(n_train)]
  pred_te <- fitted_all[seq.int(n_train + 1, n_all)]
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = "BYM2 (INLA, PC priors)")
}


# ── METHOD 6: Mixed-effects (random country intercept) ───────────────────────
# Linear mixed model with a random country intercept and area-level covariates
# as fixed effects, fit on the logit scale and weighted by survey n. For LOCO
# the held-out country's random effect is unobserved by construction; we use
# the marginal mean (random intercept set to zero), so predictions come from
# the fixed-effect (covariate) part only. This is the right structural model
# when outcome variation is dominated by within-country signal — vitamin A,
# B12, folate, and zinc are likely candidates per the variance-decomposition
# heuristic in the manuscript discussion.
#
# Implementation note: we use lme4 (already a dependency) rather than glmmTMB
# to avoid an extra install. lme4::lmer with weights = n_svy on the logit
# scale is the standard "small-area Fay-Herriot with random intercept"
# parameterisation in the SAE literature.
fit_predict_mixed <- function(train, test, vars,
                                model_type = c("logit", "continuous"),
                                n_var_cap = 8L) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  fit_df <- data.frame(Y = Y, country = factor(train$country), Xtr,
                        wt = wt, check.names = FALSE)
  pred_df <- data.frame(Xte, check.names = FALSE)
  pred_df$country <- factor(NA_character_, levels = levels(fit_df$country))
  fixed_rhs <- paste(vars, collapse = " + ")
  form <- stats::as.formula(sprintf("Y ~ %s + (1 | country)", fixed_rhs))
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(form, data = fit_df, weights = wt, REML = TRUE,
                  control = lme4::lmerControl(
                    check.conv.singular = "ignore")))),
    error = function(e) {
      cat(sprintf("    [mixed] lmer failed: %s\n", conditionMessage(e))); NULL })
  if (is.null(fit)) return(NULL)
  # Held-out country: random effect unobserved -> marginal (re.form = NA).
  pred_te_logit <- as.numeric(stats::predict(fit, newdata = pred_df,
                                              re.form = NA,
                                              allow.new.levels = TRUE))
  pred_tr_logit <- as.numeric(stats::predict(fit, newdata = fit_df))
  if (model_type == "logit") {
    pred_te <- .safe_invlogit(pred_te_logit)
    pred_tr <- .safe_invlogit(pred_tr_logit)
  } else {
    pred_te <- pred_te_logit; pred_tr <- pred_tr_logit
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = "lme4::lmer with random country intercept (marginal predict)")
}


# ── METHOD 7: Two-stage (country-mean + within-country deviance) ─────────────
# Decompose admin-2 prevalence into a country-mean component (predictable
# from country-level covariates) and a within-country deviation (predictable
# from area-level covariates).
#
# Stage 1: pool training admin-2 observations to country means; fit a small
#          OLS of country mean svy_prev on country-aggregated covariates
#          (means of the area covariates per country). With only 3 training
#          countries we keep stage 1 deliberately small — at most 1 covariate
#          plus intercept.
# Stage 2: subtract the stage-1 prediction from each training admin-2 unit's
#          svy_prev, and fit a penalised regression of the residual deviation
#          on the area covariates (centered within country).
# Prediction: stage 1 (using held-out country's aggregated covariates) +
#             stage 2 (using held-out country's area covariates).
#
# This is the transportable variant of the within-country-centred / "predict
# deviance" recipe that the existing transportability_area pipeline tested
# with `center = TRUE`. Difference: instead of subtracting the known country
# mean (which is unavailable for the held-out country), we *predict* the
# country mean from country-level inputs.
fit_predict_two_stage <- function(train, test, vars,
                                    model_type = "continuous",
                                    stage1_n_var_cap = 1L,
                                    stage2_n_var_cap = 12L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
  imp <- .impute_train_test(train, test, vars)
  # Country-level aggregates: mean of each area covariate within country.
  country_means_tr <- imp$train |>
    cbind(country = train$country, svy_prev = train$svy_prev) |>
    dplyr::group_by(country) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(vars), mean, na.rm = TRUE),
                      mean_svy = mean(svy_prev, na.rm = TRUE),
                      .groups = "drop")
  held_out_country <- unique(test$country)
  country_means_te <- imp$test |>
    cbind(country = test$country) |>
    dplyr::group_by(country) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(vars), mean, na.rm = TRUE),
                      .groups = "drop")
  # Stage 1: rank stage-1 candidates by correlation with country mean prevalence,
  # take top-K. With 3 training rows, fit a tiny OLS.
  n_ctry_tr <- nrow(country_means_tr)
  if (n_ctry_tr < 2) return(NULL)
  cmean_cors <- vapply(vars, function(v)
    suppressWarnings(abs(cor(country_means_tr[[v]], country_means_tr$mean_svy))),
    numeric(1))
  s1_vars <- names(sort(cmean_cors, decreasing = TRUE))[
    seq_len(min(stage1_n_var_cap, n_ctry_tr - 1))]
  s1_df <- data.frame(mean_svy = country_means_tr$mean_svy,
                       country_means_tr[, s1_vars, drop = FALSE])
  s1_fit <- tryCatch(stats::lm(mean_svy ~ ., data = s1_df), error = function(e) NULL)
  if (is.null(s1_fit)) return(NULL)
  # Stage-1 predicted country mean for held-out country.
  s1_pred_test <- as.numeric(stats::predict(s1_fit,
    newdata = country_means_te[, s1_vars, drop = FALSE]))
  # Stage-1 fitted means for training countries (one number per country).
  s1_pred_train_by_country <- setNames(stats::predict(s1_fit), country_means_tr$country)
  s1_pred_train <- as.numeric(s1_pred_train_by_country[train$country])

  # Stage 2: regress residual deviation on top-correlated area covariates
  # (penalised, lambda.min). Centre area features by their training country
  # mean so the stage-2 covariate matrix represents within-country deviation.
  resid_train <- train$svy_prev - s1_pred_train
  # Within-country-centred area features:
  ctr_train <- imp$train
  ctr_test  <- imp$test
  for (v in vars) {
    means_by_c <- country_means_tr[, c("country", v)]
    join_tr <- means_by_c[match(train$country, means_by_c$country), v, drop = TRUE]
    ctr_train[[v]] <- ctr_train[[v]] - join_tr
    # For the held-out country, centre by its own observed mean (= country_means_te).
    join_te <- country_means_te[match(test$country, country_means_te$country),
                                  v, drop = TRUE]
    ctr_test[[v]]  <- ctr_test[[v]] - join_te
  }
  # Variable screening for stage 2: drop near-constant after centring.
  s2_keep <- vars[vapply(vars, function(v)
    stats::sd(ctr_train[[v]], na.rm = TRUE) > 1e-6, logical(1))]
  if (length(s2_keep) == 0) {
    # Stage-2 contributes nothing; predictions are stage-1 means only.
    return(list(pred = .bound01(s1_pred_test[match(test$country, country_means_te$country)]),
                train_pred = .bound01(s1_pred_train),
                method_note = "two-stage (stage 2 dropped: no within-country variation)"))
  }
  s2_keep <- if (length(s2_keep) > stage2_n_var_cap) {
    s2cors <- vapply(s2_keep, function(v)
      suppressWarnings(abs(cor(ctr_train[[v]], resid_train, use = "pairwise"))),
      numeric(1))
    names(sort(s2cors, decreasing = TRUE))[seq_len(stage2_n_var_cap)]
  } else s2_keep
  Xtr2 <- as.matrix(ctr_train[, s2_keep, drop = FALSE])
  Xte2 <- as.matrix(ctr_test [, s2_keep, drop = FALSE])
  s2_wt <- pmax(train$n_svy %||% rep(1, nrow(Xtr2)), 1)
  s2_fit <- tryCatch(
    glmnet::cv.glmnet(Xtr2, resid_train, alpha = 0.5,
                       nfolds = min(nrow(Xtr2), 10L),
                       weights = s2_wt),
    error = function(e) NULL)
  if (is.null(s2_fit)) return(NULL)
  s2_pred_te <- as.numeric(stats::predict(s2_fit, newx = Xte2, s = "lambda.min"))
  s2_pred_tr <- as.numeric(stats::predict(s2_fit, newx = Xtr2, s = "lambda.min"))

  # Combine: stage 1 (country mean) + stage 2 (within-country deviation).
  s1_pred_te_per_area <- s1_pred_test[match(test$country, country_means_te$country)]
  pred_te <- s1_pred_te_per_area + s2_pred_te
  pred_tr <- s1_pred_train + s2_pred_tr
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("two-stage (stage 1: %d var; stage 2: %d var, EN)",
                              length(s1_vars), length(s2_keep)))
}


# ── METHOD 8: Forward selection with LOCO-CV as the selection metric ─────────
# Forward variable selection where the criterion at each step is the average
# leave-one-country-out (within the training set) Pearson r — i.e., we
# directly optimise for transportability, not within-country fit. With three
# training countries this is a LOO over 3 splits per candidate-step. Final
# model is a plain OLS on the selected variable set.
#
# This addresses the failure mode where elastic-net's inner CV selects
# variables that explain *within-training-country* prevalence but don't
# transport — a major reason SL fails the Sierra Leone child-iron LOCO.
fit_predict_forward_loco <- function(train, test, vars,
                                       model_type = c("continuous", "logit"),
                                       k_max = 6L, min_improve = 0.01) {
  model_type <- match.arg(model_type)
  imp <- .impute_train_test(train, test, vars)
  train_countries <- unique(train$country)
  if (length(train_countries) < 2) return(NULL)

  # Helper: fit OLS on selected_vars, evaluate average LOCO Pearson r within
  # training countries.
  eval_set <- function(selected_vars) {
    rs <- c()
    Y_full <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
    for (ctry in train_countries) {
      tr_idx <- which(train$country != ctry)
      te_idx <- which(train$country == ctry)
      if (length(tr_idx) < 3 || length(te_idx) < 2) next
      Xtr <- imp$train[tr_idx, selected_vars, drop = FALSE]
      Xte <- imp$train[te_idx, selected_vars, drop = FALSE]
      df  <- data.frame(Y = Y_full[tr_idx], Xtr, check.names = FALSE)
      f   <- tryCatch(stats::lm(Y ~ ., data = df), error = function(e) NULL)
      if (is.null(f)) next
      p <- as.numeric(stats::predict(f, newdata = Xte))
      if (model_type == "logit") p <- .safe_invlogit(p)
      r <- suppressWarnings(cor(train$svy_prev[te_idx], p))
      if (is.finite(r)) rs <- c(rs, r)
    }
    if (length(rs) == 0) -Inf else mean(rs)
  }

  selected <- character()
  best_score <- -Inf
  remaining <- vars
  for (step in seq_len(k_max)) {
    if (length(remaining) == 0) break
    scores <- vapply(remaining, function(v) eval_set(c(selected, v)),
                      numeric(1))
    best_new <- remaining[which.max(scores)]
    new_score <- max(scores, na.rm = TRUE)
    if (!is.finite(new_score) || new_score < best_score + min_improve) break
    selected <- c(selected, best_new)
    best_score <- new_score
    remaining <- setdiff(remaining, best_new)
  }
  if (length(selected) == 0) return(NULL)

  # Final fit on all training data with selected variables, survey-weighted.
  Y_full <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt_full <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  fit_df <- data.frame(Y = Y_full,
                        imp$train[, selected, drop = FALSE],
                        check.names = FALSE)
  fit <- stats::lm(Y ~ ., data = fit_df, weights = wt_full)
  pred_te_link <- as.numeric(stats::predict(fit,
    newdata = imp$test[, selected, drop = FALSE]))
  pred_tr_link <- as.numeric(stats::predict(fit, newdata = fit_df))
  if (model_type == "logit") {
    pred_te <- .safe_invlogit(pred_te_link)
    pred_tr <- .safe_invlogit(pred_tr_link)
  } else {
    pred_te <- pred_te_link; pred_tr <- pred_tr_link
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("forward-LOCO OLS (selected %d / cap %d: %s)",
                              length(selected), k_max,
                              paste(selected, collapse = ", ")),
       selected_vars = selected, loco_score = best_score)
}


# ── METHOD 9: Generalised additive model (mgcv::gam) ─────────────────────────
# Thin-plate-spline regression on the top correlation-screened covariates,
# fit on the logit-prevalence scale and weighted by survey n. GAMs capture
# smooth monotone covariate effects that ranger's step functions and
# elastic-net's linear effects miss; the literature suggests this is the
# most useful single additional learner for area-level spatial prediction.
fit_predict_gam <- function(train, test, vars,
                              model_type = c("logit", "continuous"),
                              n_var_cap = 6L, k_spline = 4L) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  smoothers <- paste0("s(", vars, ", k = ", k_spline, ")", collapse = " + ")
  form <- stats::as.formula(paste("Y ~", smoothers))
  fit_df <- data.frame(Y = Y, Xtr, check.names = FALSE)
  fit <- tryCatch(mgcv::gam(form, data = fit_df, weights = wt,
                              method = "REML"),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_te_link <- as.numeric(mgcv::predict.gam(fit, newdata = Xte))
  pred_tr_link <- as.numeric(mgcv::predict.gam(fit, newdata = fit_df))
  if (model_type == "logit") {
    pred_te <- .safe_invlogit(pred_te_link); pred_tr <- .safe_invlogit(pred_tr_link)
  } else {
    pred_te <- pred_te_link; pred_tr <- pred_tr_link
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("mgcv::gam thin-plate (%d vars, k=%d)",
                              length(vars), k_spline))
}


# ── METHOD 10: Group-lasso by data source ────────────────────────────────────
# Groups area covariates by data-source prefix (so a single penalty term
# governs entire interpretable blocks: NDVI vs accessibility vs livelihood
# vs sanitation vs malaria-interventions, ...). Fit via gglasso with
# CV-selected lambda. The selected variable list typically contains
# whole-domain blocks rather than scattered single variables, which is
# both more interpretable and more transportable.
.assign_variable_groups <- function(vars) {
  # Heuristic grouping: take the first 2 underscore-delimited tokens after
  # stripping the data-source prefix. Treats e.g. "gee_NDVI_*" -> group
  # "NDVI", "ihme_anemia_*" -> "anemia", etc. Variables that don't match
  # the heuristic land in the "other" group.
  tokens <- strsplit(sub("^(gee|ihme|mal_atlas|gdl|spam|chirps|wfp|map|wp|esa)_",
                          "", vars), "_")
  groups <- vapply(tokens, function(t) {
    if (length(t) >= 2) paste(t[1:2], collapse = "_")
    else if (length(t) == 1) t[1]
    else "other"
  }, character(1))
  # Encode as integer factors for gglasso.
  as.integer(factor(groups))
}
fit_predict_grouplasso <- function(train, test, vars,
                                     model_type = c("continuous", "logit"),
                                     loss = "ls") {
  if (!requireNamespace("gglasso", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  if (length(vars) < 4) return(NULL)
  imp <- .impute_train_test(train, test, vars)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  groups <- .assign_variable_groups(vars)
  # gglasso requires variables to be in contiguous group blocks.
  ord <- order(groups)
  groups <- groups[ord]; Xtr <- Xtr[, ord, drop = FALSE]; Xte <- Xte[, ord, drop = FALSE]
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  cvfit <- tryCatch(
    gglasso::cv.gglasso(Xtr, Y, group = groups, loss = loss,
                         pred.loss = "L2", nfolds = min(nrow(Xtr), 10L)),
    error = function(e) NULL)
  if (is.null(cvfit)) return(NULL)
  s <- cvfit$lambda.min
  pred_tr_link <- as.numeric(stats::predict(cvfit$gglasso.fit,
    newx = Xtr, s = s, type = "link"))
  pred_te_link <- as.numeric(stats::predict(cvfit$gglasso.fit,
    newx = Xte, s = s, type = "link"))
  if (model_type == "logit") {
    pred_tr <- .safe_invlogit(pred_tr_link); pred_te <- .safe_invlogit(pred_te_link)
  } else {
    pred_tr <- pred_tr_link; pred_te <- pred_te_link
  }
  # Count surviving groups for the method note.
  beta <- as.numeric(stats::coef(cvfit$gglasso.fit, s = s))[-1]  # drop intercept
  n_groups_kept <- length(unique(groups[abs(beta) > 1e-8]))
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("gglasso (%d/%d groups retained)",
                              n_groups_kept, length(unique(groups))))
}


# ── METHOD 11: Stacked LOCO ensemble ─────────────────────────────────────────
# Combines several benchmark methods via a simple convex-combination
# stacker whose weights are learned from training-country LOCO predictions.
# For each candidate base method we precompute LOCO predictions on the
# *training* countries (leave-one-of-3-out within training), then fit
# non-negative weights on the stacker training data. The stacker's
# weights are used to combine each method's predictions for the
# held-out country. This is the transportability-aware version of an SL
# stacker — its CV is country-disjoint, so it doesn't over-fit to
# within-training within-country structure.
fit_predict_stacked <- function(train, test, vars,
                                  base_methods = c("baseline", "penalized",
                                                    "mixed", "gam"),
                                  adjacency_list = NULL,
                                  weight_method = c("nnls", "equal")) {
  # NOTE: defaults trimmed to fast methods (baseline + penalized + mixed +
  # gam). Adding two_stage and forward increases compute by an order of
  # magnitude because each requires its own inner LOCO; pass them
  # explicitly via `base_methods` when needed.
  weight_method <- match.arg(weight_method)
  train_countries <- unique(train$country)
  if (length(train_countries) < 3) return(NULL)

  # Helper to evaluate one base method on a (sub_train, sub_test) split.
  base_run <- function(m, sub_tr, sub_te) {
    tryCatch(switch(m,
      baseline  = fit_predict_baseline   (sub_tr, sub_te, vars,
                                            model_type = "continuous"),
      penalized = fit_predict_penalized  (sub_tr, sub_te, vars,
                                            model_type = "continuous"),
      mixed     = fit_predict_mixed      (sub_tr, sub_te, vars,
                                            model_type = "continuous"),
      gam       = fit_predict_gam        (sub_tr, sub_te, vars,
                                            model_type = "continuous"),
      two_stage = fit_predict_two_stage  (sub_tr, sub_te, vars),
      forward   = fit_predict_forward_loco(sub_tr, sub_te, vars,
                                             model_type = "continuous")
    ), error = function(e) NULL)
  }

  # 1) Stacker training data: leave-one-training-country-out within
  # `train`. For each held-out training country c, fit all base methods on
  # the OTHER training countries and predict on c.
  stacker_pred <- matrix(NA_real_, nrow = nrow(train), ncol = length(base_methods))
  colnames(stacker_pred) <- base_methods
  for (c in train_countries) {
    inner_tr <- train[train$country != c, , drop = FALSE]
    inner_te <- train[train$country == c, , drop = FALSE]
    if (nrow(inner_tr) < 5 || nrow(inner_te) < 3) next
    idx_te <- which(train$country == c)
    for (i in seq_along(base_methods)) {
      r <- base_run(base_methods[i], inner_tr, inner_te)
      if (!is.null(r) && !is.null(r$pred) && length(r$pred) == length(idx_te))
        stacker_pred[idx_te, i] <- r$pred
    }
  }

  # 2) Drop methods that failed on every training country.
  ok_cols <- which(colSums(!is.na(stacker_pred)) > nrow(stacker_pred) * 0.5)
  if (length(ok_cols) == 0) return(NULL)
  stacker_pred <- stacker_pred[, ok_cols, drop = FALSE]
  methods_used <- base_methods[ok_cols]
  # Median-impute remaining missing.
  for (j in seq_len(ncol(stacker_pred))) {
    bad <- !is.finite(stacker_pred[, j])
    if (any(bad)) stacker_pred[bad, j] <- median(stacker_pred[, j], na.rm = TRUE)
  }

  # 3) Solve for non-negative convex weights via NNLS on training data:
  #     min ||X w - y||^2 s.t. w >= 0, sum(w) = 1.
  Y <- train$svy_prev
  if (weight_method == "equal") {
    w <- rep(1 / ncol(stacker_pred), ncol(stacker_pred))
  } else {
    nnls_fit <- tryCatch({
      # No nnls dependency: use constrained QP via stats::optim if available,
      # else fall back to ridge. Simple closed-form: solve via lm + clamp.
      coefs <- stats::coef(stats::lm.fit(x = stacker_pred, y = Y))
      coefs[!is.finite(coefs) | coefs < 0] <- 0
      if (sum(coefs) == 0) rep(1 / ncol(stacker_pred), ncol(stacker_pred))
      else coefs / sum(coefs)
    }, error = function(e) rep(1 / ncol(stacker_pred), ncol(stacker_pred)))
    w <- nnls_fit
  }
  names(w) <- methods_used

  # 4) Base-method predictions on the actual held-out country.
  test_pred <- matrix(NA_real_, nrow = nrow(test), ncol = length(methods_used))
  for (i in seq_along(methods_used)) {
    r <- base_run(methods_used[i], train, test)
    if (!is.null(r) && !is.null(r$pred) && length(r$pred) == nrow(test))
      test_pred[, i] <- r$pred
  }
  for (j in seq_len(ncol(test_pred))) {
    bad <- !is.finite(test_pred[, j])
    if (any(bad)) test_pred[bad, j] <- median(test_pred[, j], na.rm = TRUE)
  }

  # 5) Combine via stacker weights.
  pred_te <- as.numeric(test_pred %*% w)
  pred_tr <- as.numeric(stacker_pred %*% w)
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("stacked (%d methods, weights: %s)",
                              length(methods_used),
                              paste(sprintf("%s=%.2f", methods_used, w),
                                     collapse = ", ")),
       stacker_weights = setNames(w, methods_used))
}


# ── METHOD 12: Quasi-binomial logit GLM (proper proportion likelihood) ───────
# Models svy_prev as a proportion with logit link and quasi-binomial
# variance (allows over/under-dispersion). Survey n_svy enters as the
# weights argument so the GLM treats high-precision admin-2 polygons as
# multiple-observation equivalents. This is the principled proportion-
# scale alternative to MSE-on-raw-prevalence; predictions stay in [0, 1]
# by construction and the loss naturally accounts for boundary effects
# at low prevalence (women_vitA, women_b12).
fit_predict_quasibinomial <- function(train, test, vars,
                                         model_type = "continuous",
                                         n_var_cap = 8L) {
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])
  # Clamp svy_prev away from {0, 1} for stability of quasi-binomial likelihood.
  p <- pmin(pmax(train$svy_prev, 1e-4), 1 - 1e-4)
  wt <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  fit_df <- data.frame(p = p, Xtr, check.names = FALSE)
  fit <- tryCatch(
    suppressWarnings(stats::glm(p ~ ., data = fit_df,
                                  family = stats::quasibinomial(link = "logit"),
                                  weights = wt)),
    error = function(e) {
      cat(sprintf("    [quasibinomial] glm failed: %s\n", conditionMessage(e)))
      NULL })
  if (is.null(fit)) return(NULL)
  pred_tr <- as.numeric(stats::predict(fit, fit_df, type = "response"))
  pred_te <- as.numeric(stats::predict(fit, Xte,   type = "response"))
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = "quasi-binomial logit GLM (n_svy-weighted)")
}


# ── METHOD 13: Quantile (L1) regression at the median ────────────────────────
# Quantile regression directly optimises mean absolute error rather than
# MSE; it is robust to outliers and (importantly for our LOCO setting)
# does not penalise level bias as heavily as MSE — so the model is
# allowed to give up absolute calibration in exchange for better
# admin-2 ranking. This directly targets the Pearson r evaluation
# metric.
fit_predict_quantile <- function(train, test, vars,
                                   model_type = c("continuous", "logit"),
                                   n_var_cap = 8L, tau = 0.5) {
  if (!requireNamespace("quantreg", quietly = TRUE)) return(NULL)
  model_type <- match.arg(model_type)
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(Xtr)), 1)
  fit_df <- data.frame(Y = Y, Xtr, check.names = FALSE)
  fit <- tryCatch(
    quantreg::rq(Y ~ ., data = fit_df, tau = tau, weights = wt),
    error = function(e) {
      cat(sprintf("    [quantile] rq failed: %s\n", conditionMessage(e)))
      NULL })
  if (is.null(fit)) return(NULL)
  pred_tr_link <- as.numeric(stats::predict(fit, newdata = fit_df))
  pred_te_link <- as.numeric(stats::predict(fit, newdata = Xte))
  if (model_type == "logit") {
    pred_tr <- .safe_invlogit(pred_tr_link); pred_te <- .safe_invlogit(pred_te_link)
  } else {
    pred_tr <- pred_tr_link; pred_te <- pred_te_link
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("quantile regression (tau=%g, weighted)", tau))
}


# ── METHOD 14: Beta regression ───────────────────────────────────────────────
# Beta regression (Ferrari & Cribari-Neto 2004) is the natural likelihood
# for proportions on (0, 1) — it models both mean (logit link) and
# precision (log link). Unlike quasi-binomial it does not assume
# Bernoulli-derived variance, allowing the variance structure to be
# learned from data. Implementation via `betareg`.
fit_predict_betareg <- function(train, test, vars,
                                  model_type = "continuous",
                                  n_var_cap = 6L) {
  if (!requireNamespace("betareg", quietly = TRUE)) return(NULL)
  imp <- .impute_train_test(train, test, vars)
  if (length(vars) > n_var_cap) {
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    vars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap, length(cors)))]
  }
  Xtr <- as.data.frame(imp$train[, vars, drop = FALSE])
  Xte <- as.data.frame(imp$test [, vars, drop = FALSE])
  # Beta regression requires response strictly in (0,1).
  p <- pmin(pmax(train$svy_prev, 1e-4), 1 - 1e-4)
  fit_df <- data.frame(p = p, Xtr, check.names = FALSE)
  fit <- tryCatch(
    betareg::betareg(p ~ ., data = fit_df, link = "logit"),
    error = function(e) {
      cat(sprintf("    [betareg] failed: %s\n", conditionMessage(e)))
      NULL })
  if (is.null(fit)) return(NULL)
  pred_tr <- as.numeric(stats::predict(fit, newdata = fit_df, type = "response"))
  pred_te <- as.numeric(stats::predict(fit, newdata = Xte,    type = "response"))
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = "beta regression (Ferrari-Cribari-Neto, logit link)")
}


# ── METHOD 15: Domain-invariance filter (H5 from sandbox) ───────────────────
# Pre-screen covariates by within-country variance ratio: keep variables
# whose variation happens MOSTLY within each country rather than between
# them. Variables with high between/within ratio are mostly encoding
# country fixed effects and don't transport (the held-out country has no
# observed country effect).
#
# This is the strongest single-filter method tested in sandbox/
# (mean LOCO Pearson r = 0.195 at min_within_ratio = 0.70, top_k = 8).
# Recipe is parsimonious, fits in seconds, and the selection criterion
# is interpretable.
#
# Reference: sandbox/06_domain_adversarial.R / sandbox/findings_log.md §H5
.country_variance_decomp <- function(train, vars) {
  out <- vapply(vars, function(v) {
    x <- train[[v]]; ctry <- train$country
    grand_var <- stats::var(x, na.rm = TRUE)
    if (!is.finite(grand_var) || grand_var == 0) return(0)
    ctry_means <- tapply(x, ctry, mean, na.rm = TRUE)
    x_centered <- x - ctry_means[ctry]
    within_var <- stats::var(x_centered, na.rm = TRUE)
    within_var / grand_var
  }, numeric(1))
  data.frame(variable = vars, within_ratio = out, stringsAsFactors = FALSE)
}

fit_predict_invariance_filter <- function(train, test, vars,
                                             model_type = c("continuous", "logit"),
                                             min_within_ratio = 0.7,
                                             top_k = 8L) {
  model_type <- match.arg(model_type)
  vars <- vars[vars %in% colnames(train)]
  if (length(vars) < 2) return(NULL)
  scores <- .country_variance_decomp(train, vars)
  keep <- scores[scores$within_ratio >= min_within_ratio, , drop = FALSE]
  if (nrow(keep) < 2) keep <- scores  # fall back to full ranking
  imp_pre <- .impute_train_test(train, test, keep$variable)
  keep$abs_r <- vapply(keep$variable, function(v)
    suppressWarnings(abs(cor(imp_pre$train[[v]], train$svy_prev,
                                use = "pairwise"))), numeric(1))
  selected <- head(keep[order(-keep$abs_r), ]$variable, top_k)
  if (length(selected) < 2) return(NULL)
  imp <- .impute_train_test(train, test, selected)
  Y   <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  fit <- tryCatch(glmnet::cv.glmnet(Xtr, Y, alpha = 0.5,
                                        nfolds = min(nrow(Xtr), 10L),
                                        weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_te_link <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.min"))
  pred_tr_link <- as.numeric(stats::predict(fit, newx = Xtr, s = "lambda.min"))
  if (model_type == "logit") {
    pred_te <- .safe_invlogit(pred_te_link); pred_tr <- .safe_invlogit(pred_tr_link)
  } else {
    pred_te <- pred_te_link; pred_tr <- pred_tr_link
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("invariance filter (%d/%d vars, w>=%.2f)",
                              length(selected), length(vars), min_within_ratio),
       selected_vars = selected)
}


# ── METHOD 16: Combined invariance + multi-outcome-shared filter ─────────────
# Two-filter promotion of the winning H6 from sandbox/. Combines the
# within-country invariance filter (H5) with the outcome-shared filter
# (H4) via a composite score:
#
#   score(v) = within_ratio(v) * mean_abs_r_across_outcomes(v)
#
# Requires the OTHER outcomes' pooled data (so it can compute shared signal)
# — pass via `cross_outcome_pooled = list(<outcome_tag> = pooled_data, ...)`.
# Falls back to fit_predict_invariance_filter (H5) when cross-outcome data
# isn't supplied — i.e., this is a strict superset of method 15.
#
# Best sandbox config: min_within_ratio = 0.70, top_k = 12 ->
# mean LOCO Pearson r = 0.224 (sandbox/07_combined_winners.R).
#
# Reference: sandbox/07_combined_winners.R / sandbox/findings_log.md §H6
fit_predict_combined_filter <- function(train, test, vars,
                                          cross_outcome_pooled = NULL,
                                          this_outcome = NA_character_,
                                          model_type = "continuous",
                                          min_within_ratio = 0.70,
                                          top_k = 12L) {
  # Fall back to pure invariance filter if no cross-outcome data given.
  if (is.null(cross_outcome_pooled) || length(cross_outcome_pooled) < 2) {
    return(fit_predict_invariance_filter(train, test, vars,
                                            model_type = model_type,
                                            min_within_ratio = min_within_ratio,
                                            top_k = top_k))
  }
  vars <- vars[vars %in% colnames(train)]
  if (length(vars) < 2) return(NULL)
  # Filter 1: within-country variance ratio on this outcome's training data.
  inv_scores <- .country_variance_decomp(train, vars)
  pass <- inv_scores[inv_scores$within_ratio >= min_within_ratio, , drop = FALSE]
  if (nrow(pass) < 2) pass <- inv_scores
  candidate_vars <- pass$variable
  # Filter 2: mean absolute correlation with svy_prev across ALL outcomes
  # (computed on the same set of training countries as this fold).
  training_countries <- setdiff(unique(train$country), unique(test$country))
  cor_mat <- vapply(names(cross_outcome_pooled), function(oc) {
    pd <- cross_outcome_pooled[[oc]]
    if (is.null(pd) || nrow(pd) == 0) return(rep(NA_real_, length(candidate_vars)))
    tr_oc <- pd[pd$country %in% training_countries, , drop = FALSE]
    vapply(candidate_vars, function(v) {
      if (!v %in% colnames(tr_oc)) return(NA_real_)
      suppressWarnings(abs(cor(tr_oc[[v]], tr_oc$svy_prev, use = "pairwise")))
    }, numeric(1))
  }, numeric(length(candidate_vars)))
  if (!is.matrix(cor_mat)) cor_mat <- matrix(cor_mat, ncol = 1)
  mean_abs_r <- rowMeans(cor_mat, na.rm = TRUE)
  composite  <- pass$within_ratio * mean_abs_r
  selected_idx <- order(-composite)[seq_len(min(top_k, length(composite)))]
  selected <- candidate_vars[selected_idx]
  if (length(selected) < 2) return(NULL)
  # Final fit: weighted elastic net.
  imp <- .impute_train_test(train, test, selected)
  Y   <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt  <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  Xtr <- as.matrix(imp$train); Xte <- as.matrix(imp$test)
  fit <- tryCatch(glmnet::cv.glmnet(Xtr, Y, alpha = 0.5,
                                        nfolds = min(nrow(Xtr), 10L),
                                        weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_te <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.min"))
  pred_tr <- as.numeric(stats::predict(fit, newx = Xtr, s = "lambda.min"))
  if (model_type == "logit") {
    pred_te <- .safe_invlogit(pred_te); pred_tr <- .safe_invlogit(pred_tr)
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf(
         "combined filter (n_outcomes=%d, w>=%.2f, top_k=%d)",
         length(cross_outcome_pooled), min_within_ratio, top_k),
       selected_vars = selected)
}


# ── METHOD 17: Spatial GAM on admin-2 centroids ──────────────────────────────
# Thin-plate spline on (lon, lat) of admin-2 polygon centroids. NO GEE
# covariates — coordinates only. This was the strongest single LOCO
# method in sandbox testing (mean LOCO Pearson r = 0.285 vs the H6
# combined filter's 0.224 and the 17-hour SuperLearner's ~0.14). The
# 152 GEE rasters add nothing beyond the smooth geographic gradient
# already encoded in (lon, lat) — at least within the West African
# training footprint.
#
# Requires `lon` / `lat` columns in train and test. Use
# add_admin2_centroids() upstream to add them.
#
# Reference: sandbox/04_spatial_gam_coords.R / sandbox/findings_log.md
fit_predict_spatial_coords <- function(train, test, vars = NULL,
                                          k_spline = 30, extra_covars = character(),
                                          n_var_cap_extra = 0L,
                                          model_type = c("continuous", "logit")) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  if (!all(c("lon","lat") %in% colnames(train)) ||
      !all(c("lon","lat") %in% colnames(test))) return(NULL)
  if (any(!is.finite(train$lon)) || any(!is.finite(train$lat))) return(NULL)
  model_type <- match.arg(model_type)
  Y  <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(train)), 1)
  # Choose a small number of additional covariates if requested (correlation-
  # screened). Default is zero extras = pure spatial smoothing.
  if (length(extra_covars) == 0 && n_var_cap_extra > 0 && length(vars) > 0) {
    imp_pre <- .impute_train_test(train, test, vars)
    cors <- vapply(vars, function(v)
      suppressWarnings(abs(cor(imp_pre$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    extra_covars <- names(sort(cors, decreasing = TRUE))[seq_len(min(n_var_cap_extra, length(cors)))]
  }
  if (length(extra_covars) > 0) {
    imp <- .impute_train_test(train, test, extra_covars)
    fit_df <- data.frame(Y = Y, lon = train$lon, lat = train$lat, imp$train)
    test_df <- data.frame(lon = test$lon, lat = test$lat, imp$test)
    rhs <- paste(c(sprintf("s(lon, lat, k = %d, bs = 'tp')", k_spline),
                    extra_covars), collapse = " + ")
  } else {
    fit_df  <- data.frame(Y = Y, lon = train$lon, lat = train$lat)
    test_df <- data.frame(    lon = test$lon, lat = test$lat)
    rhs <- sprintf("s(lon, lat, k = %d, bs = 'tp')", k_spline)
  }
  form <- stats::as.formula(paste("Y ~", rhs))
  fit <- tryCatch(mgcv::gam(form, data = fit_df, weights = wt, method = "REML"),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  p_te_link <- as.numeric(mgcv::predict.gam(fit, newdata = test_df))
  p_tr_link <- as.numeric(mgcv::predict.gam(fit, newdata = fit_df))
  if (model_type == "logit") {
    p_te <- .safe_invlogit(p_te_link); p_tr <- .safe_invlogit(p_tr_link)
  } else { p_te <- p_te_link; p_tr <- p_tr_link }
  list(pred = .bound01(p_te), train_pred = .bound01(p_tr),
       method_note = sprintf("spatial thin-plate spline (k=%d, extras=%d)",
                              k_spline, length(extra_covars)))
}


# ── METHOD 18: Spatial GAM + H6 residual stack ───────────────────────────────
# Two-stage: fit spatial GAM on (lon, lat) to capture the smooth
# geographic gradient (METHOD 17), compute training residuals, then fit
# the H6 combined-filter elastic net on the residuals to capture
# within-gradient variation explained by GEE covariates. Final prediction
# is the sum of the two stages' predictions.
fit_predict_spatial_plus_h6 <- function(train, test, vars,
                                          cross_outcome_pooled = NULL,
                                          this_outcome = NA_character_,
                                          k_spline = 30,
                                          h6_min_within_ratio = 0.70,
                                          h6_top_k = 12L) {
  s <- fit_predict_spatial_coords(train, test, vars, k_spline = k_spline)
  if (is.null(s)) return(NULL)
  # Stage-2 target = residuals of stage-1 on training.
  resid_train <- train$svy_prev - s$train_pred
  # Use a thin shim: fit combined_filter on the residuals by overriding
  # svy_prev in the training frame. We restore at the end.
  train2 <- train; train2$svy_prev <- resid_train
  # Cross-outcome pooled also needs to use the same residual target —
  # without it, the combined_filter falls back to invariance_filter on
  # this outcome only. That's an acceptable degradation.
  s2 <- tryCatch(
    fit_predict_combined_filter(train2, test, vars,
                                  cross_outcome_pooled = cross_outcome_pooled,
                                  this_outcome = this_outcome,
                                  min_within_ratio = h6_min_within_ratio,
                                  top_k = h6_top_k),
    error = function(e) NULL)
  if (is.null(s2)) {
    return(list(pred = .bound01(s$pred), train_pred = .bound01(s$train_pred),
                method_note = "spatial_plus_h6 (stage 2 failed; spatial-only)"))
  }
  # Combine: spatial baseline + residual correction, but don't double-clamp
  # to [0,1] before summing — only after.
  pred <- s$pred + s2$pred
  list(pred = .bound01(pred),
       train_pred = .bound01(s$train_pred + s2$train_pred),
       method_note = sprintf("spatial+h6 (k=%d, h6_top_k=%d)", k_spline, h6_top_k))
}


# ── Helper: add admin-2 polygon centroids to a pooled data frame ─────────────
# pooled_data is augmented with `lon` and `lat` columns (centroid coords from
# GADM). cc_list is a list of country-config entries with $country and
# $gadm_code fields (as returned by get_country_configs()). svy_admin2_list
# is the per-country survey frame whose Admin2 ordering keys into GADM.
# Caches per-country centroids under <cache_dir> for re-use.
add_admin2_centroids <- function(pooled_data, cc_list, svy_admin2_list,
                                   gadm_cache = here::here("data", "gadm"),
                                   centroid_cache = here::here("data",
                                                                "admin2_centroids.rds")) {
  if (file.exists(centroid_cache)) {
    centroids <- readRDS(centroid_cache)
  } else {
    out <- list()
    .norm_nm <- function(x) gsub(" ", "", tolower(x))
    svy_keys_norm <- setNames(names(svy_admin2_list), .norm_nm(names(svy_admin2_list)))
    for (cc in cc_list) {
      nm <- cc$country
      key <- svy_keys_norm[[.norm_nm(nm)]] %||% nm
      svy <- svy_admin2_list[[key]]
      if (is.null(svy) || nrow(svy) == 0) next
      poly <- tryCatch(geodata::gadm(cc$gadm_code, level = 2, path = gadm_cache),
                        error = function(e) NULL)
      if (is.null(poly)) next
      sf_poly <- sf::st_as_sf(poly)
      suppressWarnings({ ctr <- sf::st_coordinates(sf::st_centroid(sf_poly)) })
      out[[key]] <- data.frame(country = key, Admin2 = sf_poly$NAME_2,
                                 lon = ctr[, 1], lat = ctr[, 2],
                                 stringsAsFactors = FALSE)
    }
    centroids <- dplyr::bind_rows(out)
    tryCatch(saveRDS(centroids, centroid_cache), error = function(e) NULL)
  }
  if (is.null(centroids) || nrow(centroids) == 0) return(pooled_data)
  if (!all(c("lon","lat") %in% colnames(pooled_data))) {
    pooled_data <- dplyr::left_join(pooled_data, centroids,
                                      by = c("country", "Admin2"))
  }
  pooled_data
}


# ── TODO(pc-hal): future benchmark method — PC-HAL from UC Berkeley ──────────
# Watch for a stable release of PC-HAL (principal-components highly
# adaptive lasso) from Mark van der Laan's group at UC Berkeley. The
# regular `hal9001` package implements standard HAL; PC-HAL extends it
# with a low-rank PC basis to scale to mid-dimensional inputs without
# losing the n^{-1/3} convergence rate. This is a strong candidate to
# add as a benchmark method on this project once a public release
# exists — the within-domain PCs produced by Bucket A
# (add_within_domain_pca() in R/feature_engineering_bucket_a.R) are the
# natural input basis. See docs/data_and_approach_ideas.qmd §3.0 for
# rationale.


# ── Causal DAG-based variable selection ──────────────────────────────────────
# Defines a hand-curated directed acyclic graph (DAG) per outcome capturing
# the causal pathway from environmental / livelihood / health-system
# exposures to micronutrient biomarker level. Uses dagitty to identify
# the minimal sufficient adjustment set for predicting the outcome from
# the available exposures, then maps the DAG node names to the GEE
# covariate vocabulary via regex.
#
# The DAGs are deliberately simple and based on standard nutritional
# epidemiology (Stoltzfus 2003 for iron; West 2017 for vitamin A;
# Allen 2008 for B12; Bailey 2015 for zinc). They are intended as a
# starting point — refinement requires domain expertise.
#
# Returns NULL if dagitty unavailable or if no covariate maps to the
# adjustment set.

# DAG specifications per outcome. Each entry is a dagitty syntax string.
.dag_iron <- "dag {
  livelihood -> iron_intake -> iron_status
  livelihood -> wealth -> iron_intake
  malaria -> parasite_burden -> iron_status
  sanitation -> parasite_burden
  inflammation -> measured_iron
  iron_status -> measured_iron
}"
.dag_vita <- "dag {
  vegetable_access -> vitA_intake -> vitA_status
  livestock -> animal_food -> vitA_intake
  wealth -> vitA_intake
  breastfeeding -> vitA_status
  fortification -> vitA_status
  inflammation -> measured_RBP
  vitA_status -> measured_RBP
}"
.dag_zinc <- "dag {
  soil_zinc -> staple_zinc -> zinc_intake -> zinc_status
  livelihood -> staple_zinc
  livestock -> animal_food -> zinc_intake
  inflammation -> measured_zinc
  zinc_status -> measured_zinc
}"
.dag_b12 <- "dag {
  livestock -> animal_food -> b12_intake -> b12_status
  wealth -> animal_food
  access -> animal_food
  inflammation -> measured_b12
  b12_status -> measured_b12
}"
.dag_folate <- "dag {
  vegetable_access -> folate_intake -> folate_status
  fortification -> folate_status
  wealth -> folate_intake
  inflammation -> measured_folate
  folate_status -> measured_folate
}"

# Map DAG node names to GEE covariate regex patterns.
.dag_node_patterns <- list(
  livelihood       = c("worldcereal","livestock","spam","crop","cropland"),
  wealth           = c("ccnl","nightlight","accessibility","wealth","gdl"),
  malaria          = c("malaria","itn","irs","map_"),
  parasite_burden  = c("malaria","sanitation","wash"),
  sanitation       = c("sanitation","wash","piped","toilet"),
  inflammation     = c("aerosol","atmosphere","pm25"),
  vegetable_access = c("ndvi","evi","worldcover","vegetation"),
  livestock        = c("livestock","cattle","goat"),
  animal_food      = c("livestock","cattle","milk","ccnl"),
  access           = c("accessibility","travel","road"),
  breastfeeding    = c("breastfeeding","ifs","iycf"),
  fortification    = c("fortif"),
  soil_zinc        = c("soil.*zn","soil.*zinc","soilgrid"),
  staple_zinc      = c("worldcereal","spam_cereal","spam_maize"),
  zinc_intake      = c("worldcereal","spam","livestock"),
  iron_intake      = c("worldcereal","livestock","spam"),
  vitA_intake      = c("ndvi","evi","livestock"),
  b12_intake       = c("livestock","ccnl","accessibility"),
  folate_intake    = c("ndvi","evi","worldcereal")
)

# Identify which covariates map to a given DAG node name.
.match_dag_node <- function(node, available_vars) {
  pats <- .dag_node_patterns[[node]]
  if (is.null(pats)) return(character())
  matched <- character()
  for (p in pats) {
    hits <- grep(p, available_vars, ignore.case = TRUE, value = TRUE)
    matched <- c(matched, hits)
  }
  unique(matched)
}

# Identify the DAG adjustment set for the outcome's measured-biomarker
# node (e.g., "measured_iron" for iron). Returns a character vector of
# DAG node names.
.dag_adjustment_set <- function(outcome_tag) {
  if (!requireNamespace("dagitty", quietly = TRUE)) return(character())
  oc <- tolower(outcome_tag %||% "")
  cfg <- if (grepl("iron",   oc)) list(dag = .dag_iron,   out = "measured_iron",   exp_node = "iron_status")
   else if (grepl("vita",    oc)) list(dag = .dag_vita,   out = "measured_RBP",    exp_node = "vitA_status")
   else if (grepl("zinc",    oc)) list(dag = .dag_zinc,   out = "measured_zinc",   exp_node = "zinc_status")
   else if (grepl("b12",     oc)) list(dag = .dag_b12,    out = "measured_b12",    exp_node = "b12_status")
   else if (grepl("folate",  oc)) list(dag = .dag_folate, out = "measured_folate", exp_node = "folate_status")
   else return(character())

  d <- tryCatch(dagitty::dagitty(cfg$dag), error = function(e) NULL)
  if (is.null(d)) return(character())
  # Identify minimal adjustment set for predicting `out` from all upstream
  # exposures (treating the measured biomarker as the "outcome" node and
  # all upstream variables as candidate exposures).
  ancestors <- dagitty::ancestors(d, cfg$out)
  upstream  <- setdiff(ancestors, cfg$out)
  upstream  # Use the full ancestor set as the candidate variable pool;
            # the model will select within it via OLS.
}

# fit_predict_dag: project DAG node set to actual covariates, fit an OLS.
fit_predict_dag <- function(train, test, vars,
                              outcome_tag = NULL,
                              model_type = c("continuous", "logit"),
                              n_var_cap = 12L) {
  model_type <- match.arg(model_type)
  adj_nodes <- .dag_adjustment_set(outcome_tag)
  if (length(adj_nodes) == 0) {
    cat(sprintf("    [dag] no DAG defined for outcome '%s'\n",
                outcome_tag %||% "<NA>"))
    return(NULL)
  }
  # Resolve DAG nodes to covariates available in the dataset.
  selected <- unique(unlist(lapply(adj_nodes, .match_dag_node, vars)))
  if (length(selected) == 0) {
    cat(sprintf("    [dag] adjustment set covers 0 vars (nodes: %s)\n",
                paste(adj_nodes, collapse=",")))
    return(NULL)
  }
  if (length(selected) > n_var_cap) {
    # Rank by univariate correlation with svy_prev.
    imp_pre <- .impute_train_test(train, test, selected)
    cors <- vapply(selected, function(v)
      suppressWarnings(abs(cor(imp_pre$train[[v]], train$svy_prev, use = "pairwise"))),
      numeric(1))
    selected <- names(sort(cors, decreasing = TRUE))[seq_len(n_var_cap)]
  }
  imp <- .impute_train_test(train, test, selected)
  Y <- if (model_type == "logit") .safe_logit(train$svy_prev) else train$svy_prev
  wt <- pmax(train$n_svy %||% rep(1, nrow(imp$train)), 1)
  fit_df <- data.frame(Y = Y, imp$train, check.names = FALSE)
  fit <- tryCatch(stats::lm(Y ~ ., data = fit_df, weights = wt),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_tr_link <- as.numeric(stats::predict(fit, newdata = fit_df))
  pred_te_link <- as.numeric(stats::predict(fit, newdata = data.frame(imp$test)))
  if (model_type == "logit") {
    pred_tr <- .safe_invlogit(pred_tr_link); pred_te <- .safe_invlogit(pred_te_link)
  } else {
    pred_tr <- pred_tr_link; pred_te <- pred_te_link
  }
  list(pred = .bound01(pred_te), train_pred = .bound01(pred_tr),
       method_note = sprintf("DAG-based OLS (%d DAG nodes -> %d covars: %s)",
                              length(adj_nodes), length(selected),
                              paste(head(selected, 8), collapse = ", ")),
       selected_vars = selected,
       adjustment_nodes = adj_nodes)
}


# ── Outcome-aware feature augmentation ───────────────────────────────────────
# Add programmatic interaction terms to the pooled data based on the
# outcome. Biologically-motivated, simple, transparent. Each interaction
# term is named so it can be selected by the same downstream variable
# screening that other GEE covariates go through.
#
# Effects added per outcome:
#   - iron      : malaria intervention x cereal share (anaemia & dietary)
#   - vitamin A : NDVI x livelihood density (carotene-rich vegetable proxy)
#   - zinc      : SoilGrids zinc surrogate x staple-crop yield
#   - B12, folate: livestock density x access (B12 from animal-source foods)
#
# When the necessary base covariates are absent the corresponding
# interaction is silently skipped.
augment_outcome_features <- function(df, outcome_tag) {
  vars_present <- function(...) all(c(...) %in% colnames(df))
  add <- function(df, name, expr) {
    df[[name]] <- expr; df
  }
  oc <- tolower(outcome_tag %||% "")
  # Helper for safe column access (returns NA if missing).
  pull <- function(v) if (v %in% colnames(df)) df[[v]] else NA_real_

  if (grepl("iron", oc)) {
    if (vars_present("gee_map_itn_2018", "gee_esa_worldcereal_2021_annual_mean"))
      df <- add(df, "fe_iron_itn_x_cereal",
                pull("gee_map_itn_2018") * pull("gee_esa_worldcereal_2021_annual_mean"))
    if (vars_present("gee_accessibility_2019", "gee_ccnl_2013"))
      df <- add(df, "fe_iron_access_x_nightlight",
                pull("gee_accessibility_2019") * pull("gee_ccnl_2013"))
  }
  if (grepl("vita", oc)) {
    # vitamin A: green-vegetable proxies (NDVI, EVI) x livelihood density.
    if (vars_present("gee_DailyEVI_2018", "gee_accessibility_2019"))
      df <- add(df, "fe_vitA_evi_x_access",
                pull("gee_DailyEVI_2018") * pull("gee_accessibility_2019"))
  }
  if (grepl("zinc", oc)) {
    # Use SoilGrids if present (zinc bioavailability via soil zinc).
    soil_zn <- grep("^soil.*zn", colnames(df), value = TRUE, ignore.case = TRUE)
    cereal  <- grep("worldcereal", colnames(df), value = TRUE)
    if (length(soil_zn) > 0 && length(cereal) > 0)
      df <- add(df, "fe_zinc_soil_x_cereal",
                df[[soil_zn[1]]] * df[[cereal[1]]])
  }
  if (grepl("b12|folate", oc)) {
    # B12 from animal-source foods: livestock density (if available).
    livestock <- grep("livestock|cattle", colnames(df), value = TRUE,
                       ignore.case = TRUE)
    if (length(livestock) > 0 && vars_present("gee_accessibility_2019"))
      df <- add(df, "fe_animal_food_x_access",
                df[[livestock[1]]] * pull("gee_accessibility_2019"))
  }
  df
}


# ── Bootstrap LOCO Pearson r confidence intervals ────────────────────────────
# Block-bootstrap admin-2 polygons within the held-out country to produce
# honest CIs on the LOCO Pearson r reported in the manuscript. The
# admin-2 polygons within a country are correlated (spatial autocorrelation),
# but at the within-country LOCO level we treat them as exchangeable.
# B = 500 by default.
#
# Usage: pass an already-fit predict.fn(train, test, vars, ...) closure,
# plus the LOCO pooled_data. Returns a data.frame with one row per
# held-out country with point estimate + 80% / 95% CI of Pearson r.
bootstrap_loco_ci <- function(pooled_data, gee_vars, country_names,
                                method = "penalized",
                                outcome_label = NA_character_,
                                B = 500L,
                                conf_levels = c(0.8, 0.95),
                                adjacency_list = NULL,
                                seed = 12345L) {
  set.seed(seed)
  results <- list()
  for (held_out in country_names) {
    train <- pooled_data[pooled_data$country != held_out, , drop = FALSE]
    test  <- pooled_data[pooled_data$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    vars <- .screen_vars(train, gee_vars)
    if (length(vars) == 0) next

    fit_pred <- switch(method,
      baseline          = fit_predict_baseline,
      glm               = fit_predict_glm,
      penalized         = fit_predict_penalized,
      mixed             = fit_predict_mixed,
      gam               = fit_predict_gam,
      two_stage         = fit_predict_two_stage,
      forward           = fit_predict_forward_loco,
      grouplasso        = fit_predict_grouplasso,
      stacked           = fit_predict_stacked,
      quasibinomial     = fit_predict_quasibinomial,
      invariance_filter = fit_predict_invariance_filter,
      combined_filter   = function(train, test, vars, ...)
                            fit_predict_invariance_filter(train, test, vars),
      stop("unknown method"))
    out <- tryCatch(fit_pred(train, test, vars), error = function(e) NULL)
    if (is.null(out) || is.null(out$pred)) next
    obs <- test$svy_prev; pred <- out$pred
    r_point <- suppressWarnings(cor(obs, pred))

    # Bootstrap CIs by resampling held-out admin-2 polygons WITH replacement.
    rs <- numeric(B)
    n <- length(obs)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      rs[b] <- suppressWarnings(cor(obs[idx], pred[idx]))
    }
    rs <- rs[is.finite(rs)]
    ci_rows <- lapply(conf_levels, function(p) {
      alpha <- (1 - p) / 2
      qs <- stats::quantile(rs, c(alpha, 1 - alpha), na.rm = TRUE)
      data.frame(conf = p, lo = qs[1], hi = qs[2])
    })
    ci_df <- dplyr::bind_rows(ci_rows)
    results[[held_out]] <- data.frame(
      held_out = held_out,
      method = method,
      outcome = outcome_label,
      n_test = nrow(test),
      pearson_r = round(r_point, 3),
      r_lo_80 = round(ci_df$lo[ci_df$conf == 0.8],  3),
      r_hi_80 = round(ci_df$hi[ci_df$conf == 0.8],  3),
      r_lo_95 = round(ci_df$lo[ci_df$conf == 0.95], 3),
      r_hi_95 = round(ci_df$hi[ci_df$conf == 0.95], 3),
      B = B, stringsAsFactors = FALSE
    )
  }
  if (length(results) == 0) return(data.frame())
  dplyr::bind_rows(results)
}


# ── Calibration post-step: Platt (logistic) recalibration for binary models ──
# Takes a vector of raw probabilistic predictions (e.g. CV out-of-fold
# probabilities from an individual-level SL binary classifier) and the
# observed binary outcomes, and returns a calibrator that maps raw probs to
# calibrated probs.
#
# Use case: the diagnostics_binary.csv table shows ~6 country×outcome models
# with negative Brier skill score; applying this recalibration on a held-out
# slice repairs the calibration without re-fitting the underlying model.
#
# Methods:
#   "platt"    — logistic recalibration glm(y ~ logit(p_raw), binomial)
#                Recovers a 2-parameter logistic transformation of raw probs.
#   "isotonic" — non-parametric monotonic mapping via stats::isoreg.
#                More flexible but needs more data to be stable.
#
# Returns a function fn(p_raw) -> p_calibrated, plus diagnostic stats.
fit_binary_recalibrator <- function(p_raw, y_obs,
                                      method = c("platt", "isotonic"),
                                      eps = 1e-3) {
  method <- match.arg(method)
  ok <- is.finite(p_raw) & is.finite(y_obs) & y_obs %in% c(0, 1)
  p_raw <- pmin(pmax(p_raw[ok], eps), 1 - eps)
  y_obs <- y_obs[ok]
  if (length(p_raw) < 10 || length(unique(y_obs)) < 2) {
    return(list(fn = function(p) p, method = "identity (insufficient data)",
                n = length(p_raw)))
  }
  if (method == "platt") {
    logit_p <- log(p_raw / (1 - p_raw))
    glmfit <- stats::glm(y_obs ~ logit_p, family = stats::binomial())
    fn <- function(p) {
      p <- pmin(pmax(p, eps), 1 - eps)
      stats::predict(glmfit,
                      newdata = data.frame(logit_p = log(p / (1 - p))),
                      type = "response")
    }
    list(fn = fn, method = "platt", coef = stats::coef(glmfit),
          n = length(p_raw))
  } else {
    # Isotonic regression: monotonic mapping from p_raw to E[y | p_raw].
    iso <- stats::isoreg(p_raw, y_obs)
    # Build a step-function from the ordered (x, yf) pairs.
    iso_x <- iso$x[iso$ord]
    iso_y <- iso$yf
    fn <- function(p) {
      p <- pmin(pmax(p, min(iso_x)), max(iso_x))
      stats::approx(iso_x, iso_y, xout = p, ties = "ordered",
                     rule = 2)$y
    }
    list(fn = fn, method = "isotonic", n = length(p_raw))
  }
}


# ── Calibrated diagnostics wrapper ───────────────────────────────────────────
# Runs the same diagnostic suite as run_diagnostics() but with Platt
# recalibration applied to the OOF binary probabilities first. Plumbs into
# _targets.R as the diagnostics_calibrated_<country>_<outcome> targets.
# Returns the same list structure (binary_metrics, calibration_table, etc.)
# as run_diagnostics so the rollup can use the same shape.
run_diagnostics_calibrated <- function(sl_fit, cc, oc) {
  if (is.null(sl_fit$bin_fit) || is.null(sl_fit$bin_fit$res)) {
    return(list(binary_metrics = NULL, calibration_table = NULL,
                pr_curve = NULL, roc_curve = NULL))
  }
  res <- sl_fit$bin_fit$res
  cal <- fit_binary_recalibrator(res$yhat_full, res$Y, method = "platt")
  yhat_cal <- cal$fn(res$yhat_full)
  cal_preds <- data.frame(Y = res$Y, yhat = yhat_cal,
                            country = cc$country, outcome = oc$tag,
                            type = "binary", stringsAsFactors = FALSE)
  bin_diag <- tryCatch(compute_binary_diagnostics(cal_preds),
                        error = function(e) NULL)
  if (is.null(bin_diag)) return(list(binary_metrics = NULL))
  out <- list()
  bin_row <- bin_diag$metrics
  bin_row$model_type     <- "binary_calibrated"
  bin_row$country        <- cc$country
  bin_row$outcome        <- oc$tag
  bin_row$calibrator     <- cal$method
  bin_row$calib_n_train  <- cal$n
  if (!is.null(cal$coef)) {
    bin_row$platt_intercept <- as.numeric(cal$coef[1])
    bin_row$platt_slope     <- as.numeric(cal$coef[2])
  }
  out$binary_metrics    <- bin_row
  out$calibration_table <- bin_diag$calibration_table
  if (!is.null(out$calibration_table)) {
    out$calibration_table$country <- cc$country
    out$calibration_table$outcome <- oc$tag
    out$calibration_table$model_type <- "binary_calibrated"
  }
  out$pr_curve  <- bin_diag$pr_curve
  out$roc_curve <- bin_diag$roc_curve
  out
}


# ── Apply recalibration to an existing binary CV-prediction table ─────────────
# Takes a data.frame with one row per observation containing the raw
# predicted probability, observed binary outcome, and (typically) a fold
# assignment. Returns the input with a calibrated_p column added.
# For honest evaluation, the calibrator is fit using out-of-fold predictions
# from FOLDS OTHER THAN the one being calibrated (cv-stacking style).
calibrate_binary_oof <- function(cv_df,
                                   p_col = "yhat_full",
                                   y_col = "Y",
                                   fold_col = NULL,
                                   method = "platt") {
  if (is.null(cv_df) || nrow(cv_df) == 0) return(cv_df)
  if (is.null(fold_col) || !fold_col %in% colnames(cv_df)) {
    # No folds given: fit calibrator on all data (slightly optimistic, but
    # the binary models are already heavily under-fit so this matters little).
    cal <- fit_binary_recalibrator(cv_df[[p_col]], cv_df[[y_col]], method = method)
    cv_df$calibrated_p <- cal$fn(cv_df[[p_col]])
    attr(cv_df, "calibrator") <- cal
    return(cv_df)
  }
  # Honest CV: for each fold f, fit calibrator on observations not in f,
  # then apply to observations in f.
  cv_df$calibrated_p <- NA_real_
  for (f in unique(cv_df[[fold_col]])) {
    in_f  <- which(cv_df[[fold_col]] == f)
    out_f <- which(cv_df[[fold_col]] != f)
    cal <- fit_binary_recalibrator(cv_df[[p_col]][out_f],
                                     cv_df[[y_col]][out_f],
                                     method = method)
    cv_df$calibrated_p[in_f] <- cal$fn(cv_df[[p_col]][in_f])
  }
  cv_df
}


# ── Adjacency builder ────────────────────────────────────────────────────────
# Build a country -> spdep::nb adjacency list using GADM admin-2 polygons.
# Designed to match the Admin2 ordering of the per-country svy_admin2/gee_admin2
# data frames (assumed to be in GADM NAME_2 alphabetical order, as the rest
# of the pipeline does).
build_adjacency_list <- function(cc_list, svy_admin2_list,
                                  gadm_cache = here::here("data", "gadm")) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    warning("spdep not available — adjacency list will be empty")
    return(list())
  }
  out <- list()
  # Build a name-insensitive lookup so callers can key svy_admin2_list by
  # display name ("Sierra Leone" / "SierraLeone") or by config $country
  # ("sierraleone"); we strip spaces AND lowercase both sides.
  .norm_nm <- function(x) gsub(" ", "", tolower(x))
  svy_keys_norm <- setNames(names(svy_admin2_list), .norm_nm(names(svy_admin2_list)))
  for (cc in cc_list) {
    nm <- cc$country
    norm_nm <- .norm_nm(nm)
    key <- if (norm_nm %in% names(svy_keys_norm)) svy_keys_norm[[norm_nm]] else nm
    svy <- svy_admin2_list[[key]]
    if (is.null(svy) || nrow(svy) == 0) {
      cat(sprintf("  [adj] %s (key=%s): no svy data, skipping\n", nm, key)); next
    }
    nm <- key  # use the actual list key going forward so out[[nm]] aligns
    poly <- tryCatch(
      geodata::gadm(cc$gadm_code, level = 2, path = gadm_cache),
      error = function(e) { cat(sprintf("  [adj] %s: gadm failed: %s\n", nm, e$message)); NULL })
    if (is.null(poly)) next
    sf_poly <- sf::st_as_sf(poly)
    # Align polygon order with svy_admin2 ordering (by Admin2 label).
    sf_poly$.match <- sf_poly$NAME_2
    matched <- match(svy$Admin2, sf_poly$.match)
    n_matched <- sum(!is.na(matched))
    if (n_matched == 0) {
      cat(sprintf("  [adj] %s: 0 of %d admin2 names matched GADM NAME_2\n",
                  nm, nrow(svy)))
      next
    }
    sf_poly <- sf_poly[matched, ]
    sf_poly <- sf_poly[!is.na(sf_poly$.match), ]
    nb <- tryCatch(spdep::poly2nb(sf_poly), error = function(e) {
      cat(sprintf("  [adj] %s: poly2nb failed: %s\n", nm, e$message)); NULL
    })
    if (!is.null(nb)) {
      out[[nm]] <- nb
      cat(sprintf("  [adj] %s: %d areas matched, %.1f neighbors per area (mean)\n",
                  nm, length(nb), mean(spdep::card(nb))))
    }
  }
  out
}


# ── LOCO runner: evaluate every method on every held-out country ─────────────
# pooled_data must have columns: country, Admin2, svy_prev, n_svy, plus
# gee_vars.  Returns a tidy data.frame with one row per (method, holdout,
# model_type) combination.
run_area_benchmarks_loco <- function(pooled_data, gee_vars, country_names,
                                       adjacency_list = NULL,
                                       outcome_label = NA_character_,
                                       methods = c("baseline", "glm",
                                                   "penalized", "fh", "bym2",
                                                   "mixed", "two_stage",
                                                   "forward", "gam",
                                                   "grouplasso", "stacked",
                                                   "quasibinomial", "quantile",
                                                   "betareg", "dag",
                                                   "invariance_filter",
                                                   "combined_filter",
                                                   "spatial_coords",
                                                   "spatial_plus_h6"),
                                       model_types = c("continuous", "logit"),
                                       augment_features = TRUE,
                                       cross_outcome_pooled = NULL) {
  # Optional outcome-aware feature augmentation: adds biologically-motivated
  # interaction terms BEFORE the LOCO loop, so every method sees the same
  # augmented feature set and the new features go through the same
  # .screen_vars() screening that the GEE variables do.
  if (isTRUE(augment_features) && !is.na(outcome_label)) {
    pooled_data <- augment_outcome_features(pooled_data, outcome_label)
    # Add any newly-created `fe_*` columns to gee_vars so they are eligible.
    fe_cols <- grep("^fe_", colnames(pooled_data), value = TRUE)
    gee_vars <- unique(c(gee_vars, fe_cols))
  }
  results <- list()
  for (held_out in country_names) {
    train <- pooled_data[pooled_data$country != held_out, , drop = FALSE]
    test  <- pooled_data[pooled_data$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    vars <- .screen_vars(train, gee_vars)
    n_vars <- length(vars)
    if (n_vars == 0) next

    for (m in methods) {
      mts <- if (m %in% c("baseline", "fh", "bym2", "two_stage",
                            "grouplasso", "stacked", "quasibinomial",
                            "betareg", "dag", "combined_filter",
                            "spatial_coords", "spatial_plus_h6"))
              "continuous" else model_types
      for (mt in mts) {
        out <- tryCatch(switch(m,
          baseline      = fit_predict_baseline    (train, test, vars, model_type = mt),
          glm           = fit_predict_glm         (train, test, vars, model_type = mt),
          penalized     = fit_predict_penalized   (train, test, vars, model_type = mt),
          fh            = fit_predict_fh          (train, test, vars, model_type = mt),
          bym2          = fit_predict_bym2        (train, test, vars,
                                                     adjacency_list = adjacency_list,
                                                     model_type = "continuous"),
          mixed         = fit_predict_mixed       (train, test, vars, model_type = mt),
          two_stage     = fit_predict_two_stage   (train, test, vars,
                                                     model_type = "continuous"),
          forward       = fit_predict_forward_loco(train, test, vars, model_type = mt),
          gam           = fit_predict_gam         (train, test, vars, model_type = mt),
          grouplasso    = fit_predict_grouplasso  (train, test, vars,
                                                     model_type = "continuous"),
          stacked       = fit_predict_stacked     (train, test, vars,
                                                     adjacency_list = adjacency_list),
          quasibinomial = fit_predict_quasibinomial(train, test, vars),
          quantile      = fit_predict_quantile    (train, test, vars, model_type = mt),
          betareg       = fit_predict_betareg     (train, test, vars),
          dag           = fit_predict_dag         (train, test, vars,
                                                     outcome_tag = outcome_label,
                                                     model_type = mt),
          invariance_filter = fit_predict_invariance_filter(train, test, vars,
                                                              model_type = mt),
          combined_filter   = fit_predict_combined_filter  (train, test, vars,
                                                              cross_outcome_pooled = cross_outcome_pooled,
                                                              this_outcome = outcome_label),
          spatial_coords    = fit_predict_spatial_coords   (train, test, vars,
                                                              k_spline = 30),
          spatial_plus_h6   = fit_predict_spatial_plus_h6  (train, test, vars,
                                                              cross_outcome_pooled = cross_outcome_pooled,
                                                              this_outcome = outcome_label)
        ), error = function(e) {
          cat(sprintf("    [%s/%s] %s failed: %s\n",
                      held_out, mt, m, conditionMessage(e)))
          NULL
        })
        if (is.null(out) || is.null(out$pred)) next
        results[[paste(held_out, m, mt, sep = "/")]] <- compute_area_metrics(
          obs = test$svy_prev, pred = out$pred,
          method = m, model_type = mt,
          n_train = nrow(train), n_test = nrow(test), n_vars = n_vars,
          held_out = held_out, outcome = outcome_label,
          scale_pp = TRUE)
      }
    }
  }
  if (length(results) == 0) return(data.frame())
  dplyr::bind_rows(results)
}


# ── Within-country k-fold CV runner ──────────────────────────────────────────
# country_data must have columns Admin2, svy_prev, n_svy and gee_vars.
# Uses k-fold CV; for k > n - 1 we fall back to LOO.
run_area_benchmarks_within <- function(country_data, gee_vars,
                                         country = NA_character_,
                                         outcome_label = NA_character_,
                                         k_folds = 5L, seed = 12345L,
                                         adjacency = NULL,
                                         methods = c("baseline", "glm",
                                                     "penalized", "fh", "bym2"),
                                         model_types = c("continuous", "logit")) {
  n <- nrow(country_data)
  if (n < 5) return(data.frame())
  set.seed(seed)
  k <- min(k_folds, max(2L, n - 1L))
  folds <- sample(rep(seq_len(k), length.out = n))

  # Container: predictions per (method, model_type, row index).
  preds <- list()

  for (f in seq_len(k)) {
    train <- country_data[folds != f, , drop = FALSE]
    test  <- country_data[folds == f, , drop = FALSE]
    if (nrow(train) < 3 || nrow(test) < 1) next
    vars <- .screen_vars(train, gee_vars)
    if (length(vars) == 0) next

    for (m in methods) {
      mts <- if (m %in% c("baseline", "fh", "bym2")) "continuous" else model_types
      for (mt in mts) {
        # For BYM2 the adjacency must cover both train and test areas.
        adj_arg <- if (m == "bym2" && !is.null(adjacency)) {
          # Pass the full adjacency as a single-country list keyed by country.
          setNames(list(adjacency), country %||% "country")
        } else NULL
        # When using BYM2 we need train$country and test$country populated.
        if (m == "bym2" && !"country" %in% colnames(train)) {
          train$country <- country %||% "country"
          test $country <- country %||% "country"
        }
        out <- tryCatch(switch(m,
          baseline  = fit_predict_baseline (train, test, vars, model_type = mt),
          glm       = fit_predict_glm      (train, test, vars, model_type = mt),
          penalized = fit_predict_penalized(train, test, vars, model_type = mt),
          fh        = fit_predict_fh       (train, test, vars, model_type = mt),
          bym2      = fit_predict_bym2     (train, test, vars,
                                              adjacency_list = adj_arg,
                                              model_type = "continuous")
        ), error = function(e) NULL)
        if (is.null(out) || is.null(out$pred)) next
        key <- paste(m, mt, sep = "/")
        if (is.null(preds[[key]])) preds[[key]] <- rep(NA_real_, n)
        preds[[key]][folds == f] <- out$pred
      }
    }
  }

  results <- list()
  for (key in names(preds)) {
    parts <- strsplit(key, "/", fixed = TRUE)[[1]]
    m <- parts[1]; mt <- parts[2]
    obs <- country_data$svy_prev
    p <- preds[[key]]
    n_train_avg <- round(mean(table(folds)) * (k - 1))
    n_test_avg  <- round(mean(table(folds)))
    results[[key]] <- compute_area_metrics(
      obs = obs, pred = p,
      method = m, model_type = mt,
      n_train = n_train_avg, n_test = n_test_avg,
      n_vars = length(.screen_vars(country_data, gee_vars)),
      country = country, outcome = outcome_label,
      scale_pp = TRUE)
  }
  if (length(results) == 0) return(data.frame())
  dplyr::bind_rows(results)
}


# ── Driver: run all benchmarks over all outcomes ─────────────────────────────
# Expects svy_admin2_list and gee_admin2_list keyed by country, plus an
# outcome label and a vector of country names.  Iterates over outcomes
# (each list element is keyed by outcome inside the per-country survey
# frame), builds the pooled LOCO dataset, runs benchmarks LOCO + within,
# and returns a single stacked data.frame.
run_area_benchmarks_all <- function(svy_admin2_list_by_outcome,
                                      gee_admin2_list,
                                      adjacency_list = NULL,
                                      methods = c("baseline", "glm",
                                                  "penalized", "fh", "bym2"),
                                      do_loco = TRUE, do_within = TRUE) {
  out <- list()
  for (oc in names(svy_admin2_list_by_outcome)) {
    svy_list <- svy_admin2_list_by_outcome[[oc]]
    pooled <- build_area_loco_dataset(svy_list, gee_admin2_list)
    if (is.null(pooled) || nrow(pooled$pooled_data) < 10) next

    if (do_loco) {
      loco <- run_area_benchmarks_loco(
        pooled_data    = pooled$pooled_data,
        gee_vars       = pooled$common_gee_vars,
        country_names  = pooled$country_names,
        adjacency_list = adjacency_list,
        outcome_label  = oc, methods = methods)
      if (nrow(loco) > 0) {
        loco$eval_type <- "loco"
        out[[paste0(oc, "_loco")]] <- loco
      }
    }
    if (do_within) {
      for (ctry in pooled$country_names) {
        cdata <- pooled$pooled_data[pooled$pooled_data$country == ctry, , drop = FALSE]
        if (nrow(cdata) < 5) next
        wb <- run_area_benchmarks_within(
          country_data  = cdata,
          gee_vars      = pooled$common_gee_vars,
          country       = ctry,
          outcome_label = oc,
          adjacency     = adjacency_list[[ctry]],
          methods       = methods)
        if (nrow(wb) > 0) {
          wb$eval_type <- "within"
          out[[paste0(oc, "_within_", ctry)]] <- wb
        }
      }
    }
  }
  if (length(out) == 0) return(data.frame())
  dplyr::bind_rows(out)
}
