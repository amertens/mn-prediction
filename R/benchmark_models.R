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
  fit_df <- data.frame(Y = Y, imp$train)
  fit <- tryCatch(stats::lm(Y ~ ., data = fit_df), error = function(e) NULL)
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
  # Drop alpha for very small n
  if (nrow(Xtr) < 25 && alpha > 0) alpha <- 0
  n_folds <- min(nrow(Xtr), 10L)
  fit <- tryCatch(
    glmnet::cv.glmnet(Xtr, Y, alpha = alpha, nfolds = n_folds, family = "gaussian"),
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
  # Build a case-insensitive lookup so callers can key svy_admin2_list by
  # display name ("SierraLeone") or by config $country ("sierraleone").
  svy_keys_lower <- setNames(names(svy_admin2_list), tolower(names(svy_admin2_list)))
  for (cc in cc_list) {
    nm <- cc$country
    key <- svy_keys_lower[[tolower(nm)]] %||% nm
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
                                                   "penalized", "fh", "bym2"),
                                       model_types = c("continuous", "logit")) {
  results <- list()
  for (held_out in country_names) {
    train <- pooled_data[pooled_data$country != held_out, , drop = FALSE]
    test  <- pooled_data[pooled_data$country == held_out, , drop = FALSE]
    if (nrow(train) < 5 || nrow(test) < 3) next
    vars <- .screen_vars(train, gee_vars)
    n_vars <- length(vars)
    if (n_vars == 0) next

    for (m in methods) {
      mts <- if (m %in% c("baseline", "fh", "bym2")) "continuous" else model_types
      for (mt in mts) {
        out <- tryCatch(switch(m,
          baseline  = fit_predict_baseline (train, test, vars, model_type = mt),
          glm       = fit_predict_glm      (train, test, vars, model_type = mt),
          penalized = fit_predict_penalized(train, test, vars, model_type = mt),
          fh        = fit_predict_fh       (train, test, vars, model_type = mt),
          bym2      = fit_predict_bym2     (train, test, vars,
                                              adjacency_list = adjacency_list,
                                              model_type = "continuous")
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
