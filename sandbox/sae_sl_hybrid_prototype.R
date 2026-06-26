# =============================================================================
# sandbox/sae_sl_hybrid_prototype.R
#
# PROTOTYPE — Option 1 hybrid: feed area-level SuperLearner out-of-fold
# predictions as the single covariate into a small-area-estimation model
# (Fay-Herriot and BYM2), so the SAE layer adds variance-weighted shrinkage
# and credible intervals on top of the ML mean.
#
# Question: do shrinkage + intervals improve WITHIN-COUNTRY calibration without
# hurting the area-level SL's strong rank correlation?
#
# Scope: child_iron, Ghana (75 areas) + Malawi (87 areas). MSE-loss SL variant.
#
# Honest design notes:
#  - The SL covariate is genuine k-fold OUT-OF-FOLD (fit_area_sl's own loo_preds
#    are in-sample refits, so we regenerate them here).
#  - "FH synthetic (OOF)" is a held-out, variance-weighted recalibration of the
#    SL prediction: in each outer fold, FH is fit on training areas and the
#    synthetic mean (a + b * sl_oof) is predicted for held-out areas. This is
#    NOT circular — it never sees the held-out area's direct estimate.
#  - "FH EBLUP (in-sample)" and "BYM2 (in-sample)" DO use each area's direct
#    estimate (that is what shrinkage means); we report them mainly for
#    interval width + coverage, not for "beating" r.
#
# Output: results/prototype_sae_sl/{metrics.csv, calibration_<country>.png}
# =============================================================================

suppressWarnings(suppressMessages({
  library(targets); library(dplyr); library(ggplot2); library(here)
  library(mlr3superlearner); library(mlr3extralearners)
}))

STORE   <- "_targets_full"
OUTCOME <- "child_iron"
DEFF    <- 1.5                 # DHS design-effect approximation
K       <- 10L                 # outer CV folds for OOF predictions
SEED    <- 20260615            # fixed (no Date/random in workflow context)
OUT     <- here::here("results", "prototype_sae_sl")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Compact area-level SL library (matches fit_area_sl defaults; no BART/GP)
SL_LIB <- list(
  list("glmnet", alpha = 1,   id = "lasso"),
  list("glmnet", alpha = 0,   id = "ridge"),
  list("glmnet", alpha = 0.5, id = "elastic_net"),
  list("ranger", num.trees = 500, min.node.size = 5, id = "ranger"),
  list("xgboost", max_depth = 2, eta = 0.05, nrounds = 200,
       min_child_weight = 3, subsample = 0.8, id = "xgb")
)

# ── helpers ──────────────────────────────────────────────────────────────────
load_country <- function(cn) {
  svy <- tar_read_raw(paste0("svy_admin2_", cn, "_", OUTCOME), store = STORE)
  gee <- tar_read_raw(paste0("gee_admin2_", cn), store = STORE)
  df <- dplyr::inner_join(
    svy |> dplyr::select(Admin2, svy_prev, n_svy),
    gee, by = "Admin2") |>
    dplyr::filter(!is.na(svy_prev), is.finite(svy_prev), !is.na(n_svy))
  gee_vars <- setdiff(names(gee), c("Admin1", "Admin2"))
  valid <- gee_vars[vapply(gee_vars, function(v) {
    x <- df[[v]]; sum(!is.na(x)) > 2 && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
  list(df = df, vars = valid)
}

# median impute using TRAIN-fold medians, applied to a target frame
impute_with <- function(target, train, vars) {
  for (v in vars) {
    x <- train[[v]]; if (inherits(x, "haven_labelled")) x <- as.double(unclass(x))
    med <- median(x[is.finite(x)], na.rm = TRUE)
    tx <- target[[v]]; if (inherits(tx, "haven_labelled")) tx <- as.double(unclass(tx))
    tx[!is.finite(tx)] <- med
    target[[v]] <- tx
  }
  target
}

# genuine k-fold OOF area-level SL (continuous / MSE loss)
oof_area_sl <- function(df, vars) {
  set.seed(SEED)
  n <- nrow(df); fold <- sample(rep_len(seq_len(K), n))
  oof <- rep(NA_real_, n)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    tr_i <- impute_with(df[tr, , drop = FALSE], df[tr, , drop = FALSE], vars)
    te_i <- impute_with(df[te, , drop = FALSE], df[tr, , drop = FALSE], vars)
    fit <- tryCatch(mlr3superlearner::mlr3superlearner(
      data = data.frame(Y = tr_i$svy_prev, tr_i[, vars, drop = FALSE]),
      target = "Y", library = SL_LIB, outcome_type = "continuous",
      folds = min(5L, length(tr))), error = function(e) NULL)
    if (is.null(fit)) next
    p <- tryCatch(as.numeric(predict(
      fit, data.frame(te_i[, vars, drop = FALSE]))), error = function(e) NULL)
    if (!is.null(p)) oof[te] <- pmin(pmax(p, 0), 1)
  }
  oof
}

# metric helpers
cal_slope <- function(y, x) tryCatch(unname(coef(lm(y ~ x))[2]), error=function(e) NA)
cal_int   <- function(y, x) tryCatch(unname(coef(lm(y ~ x))[1]), error=function(e) NA)
metric_row <- function(method, y, pred, lwr = NULL, upr = NULL) {
  ok <- is.finite(y) & is.finite(pred)
  cov <- NA_real_; width <- NA_real_
  if (!is.null(lwr)) {
    cov   <- round(mean(y[ok] >= lwr[ok] & y[ok] <= upr[ok]) * 100, 1)
    width <- round(mean((upr[ok] - lwr[ok])) * 100, 1)
  }
  data.frame(method = method, n = sum(ok),
             pearson_r   = round(cor(y[ok], pred[ok]), 3),
             spearman    = round(cor(y[ok], pred[ok], method = "spearman"), 3),
             calib_slope = round(cal_slope(y[ok], pred[ok]), 3),
             calib_int   = round(cal_int(y[ok], pred[ok]) * 100, 2),
             mae_pp      = round(mean(abs(y[ok] - pred[ok])) * 100, 2),
             pi95_width_pp = width, pi95_cover_pct = cov,
             stringsAsFactors = FALSE)
}

run_country <- function(cn, label) {
  cat(sprintf("\n========== %s (%s) ==========\n", label, OUTCOME))
  d <- load_country(cn); df <- d$df; vars <- d$vars
  cat(sprintf("  %d areas, %d GEE covariates\n", nrow(df), length(vars)))

  # 1) genuine OOF area-level SL
  df$sl_oof <- oof_area_sl(df, vars)
  df <- df |> dplyr::filter(is.finite(sl_oof))

  # sampling variance (DHS deff)
  p   <- pmin(pmax(df$svy_prev, 1e-4), 1 - 1e-4)
  n_e <- pmax(df$n_svy / DEFF, 1)
  df$sv <- p * (1 - p) / n_e

  res <- list()
  # ── raw SL (no uncertainty) ────────────────────────────────────────────────
  res[["raw_SL_oof"]] <- metric_row("raw_SL_oof", df$svy_prev, df$sl_oof)

  # ── FH synthetic, held-out (variance-weighted recalibration of SL) ──────────
  set.seed(SEED); fold <- sample(rep_len(seq_len(K), nrow(df)))
  fh_syn <- rep(NA_real_, nrow(df))
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    fit <- tryCatch(sae::eblupFH(
      svy_prev ~ sl_oof,
      vardir = sv,
      data = data.frame(svy_prev = df$svy_prev[tr], sl_oof = df$sl_oof[tr],
                        sv = df$sv[tr]), method = "REML"),
      error = function(e) NULL)
    if (is.null(fit)) next
    b <- tryCatch(fit$fit$estcoef[, "beta"], error = function(e) NULL)
    if (is.null(b) || length(b) != 2) next
    fh_syn[te] <- b[1] + b[2] * df$sl_oof[te]
  }
  fh_syn <- pmin(pmax(fh_syn, 0), 1)
  res[["FH_synthetic_oof"]] <- metric_row("FH_synthetic_oof", df$svy_prev, fh_syn)

  # ── FH EBLUP in-sample + intervals (shrinkage; uses direct estimate) ────────
  fh <- tryCatch(sae::mseFH(svy_prev ~ sl_oof, vardir = sv,
    data = data.frame(svy_prev = df$svy_prev, sl_oof = df$sl_oof, sv = df$sv),
    method = "REML"), error = function(e) NULL)
  if (!is.null(fh)) {
    eblup <- pmin(pmax(as.numeric(fh$est$eblup), 0), 1)
    se    <- sqrt(pmax(as.numeric(fh$mse), 0))
    # predictive interval for a NEW direct estimate: add sampling variance back
    pse <- sqrt(pmax(as.numeric(fh$mse), 0) + df$sv)
    res[["FH_eblup_insample"]] <- metric_row("FH_eblup_insample", df$svy_prev,
      eblup, lwr = eblup - 1.96 * pse, upr = eblup + 1.96 * pse)
    df$fh_eblup <- eblup; df$fh_lwr <- eblup - 1.96*pse; df$fh_upr <- eblup + 1.96*pse
  }

  # ── BYM2 via INLA: spatial smoothing on top of SL mean + credible intervals ─
  bym2 <- NULL
  if (requireNamespace("INLA", quietly = TRUE)) {
    adj <- tryCatch(tar_read_raw("area_adjacency_list", store = STORE),
                    error = function(e) NULL)
    nb <- adj[[label]]
    if (!is.null(nb)) {
      # align nb (built on svy ordering) to df rows by Admin2 if possible
      g <- tryCatch(INLA::inla.read.graph(spdep::nb2mat(nb, style = "B",
            zero.policy = TRUE)), error = function(e) NULL)
      if (!is.null(g) && g$n == nrow(df)) {
        dat <- data.frame(y = df$svy_prev, sl_oof = df$sl_oof,
                          area = seq_len(nrow(df)), scale = 1 / df$sv)
        bym2 <- tryCatch(INLA::inla(
          y ~ sl_oof + f(area, model = "bym2", graph = g,
                         scale.model = TRUE,
                         hyper = list(
                           prec = list(prior = "pc.prec", param = c(1, 0.01)),
                           phi  = list(prior = "pc",      param = c(0.5, 0.5)))),
          family = "gaussian", data = dat, scale = dat$scale,
          control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
          control.predictor = list(compute = TRUE),
          control.compute = list(return.marginals.predictor = TRUE)),
          error = function(e) { cat("    [bym2] INLA failed:", conditionMessage(e), "\n"); NULL })
      } else {
        cat(sprintf("    [bym2] adjacency size %s != %d areas; skipping spatial\n",
            if (is.null(g)) "NA" else g$n, nrow(df)))
      }
    }
  }
  if (!is.null(bym2)) {
    fitted <- pmin(pmax(bym2$summary.fitted.values$mean, 0), 1)
    lwr <- pmin(pmax(bym2$summary.fitted.values$`0.025quant`, 0), 1)
    upr <- pmin(pmax(bym2$summary.fitted.values$`0.975quant`, 0), 1)
    res[["BYM2_insample"]] <- metric_row("BYM2_insample", df$svy_prev,
      fitted, lwr = lwr, upr = upr)
  }

  out <- dplyr::bind_rows(res); out$country <- label; out$outcome <- OUTCOME
  print(out[, c("method","pearson_r","spearman","calib_slope","calib_int",
                "mae_pp","pi95_width_pp","pi95_cover_pct")], row.names = FALSE)

  # calibration plot: raw SL vs FH EBLUP
  if (!is.null(res[["FH_eblup_insample"]]) && "fh_eblup" %in% names(df)) {
    pdat <- rbind(
      data.frame(method = "raw SL (OOF)", pred = df$sl_oof, svy = df$svy_prev,
                 lwr = NA, upr = NA),
      data.frame(method = "SL -> FH EBLUP", pred = df$fh_eblup, svy = df$svy_prev,
                 lwr = df$fh_lwr, upr = df$fh_upr))
    g <- ggplot(pdat, aes(pred, svy)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
      { if (any(!is.na(pdat$lwr))) geom_errorbarh(
          aes(xmin = lwr, xmax = upr), height = 0, alpha = .3, na.rm = TRUE) } +
      geom_point(colour = "#2c7bb6", alpha = .8) +
      geom_smooth(method = "lm", se = FALSE, colour = "#d7191c", linewidth = .6) +
      facet_wrap(~method) +
      labs(title = sprintf("%s — child iron: SL vs SL->FH hybrid", label),
           x = "predicted prevalence", y = "survey (direct) prevalence") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUT, sprintf("calibration_%s.png", label)), g,
           width = 8, height = 4, dpi = 150)
  }
  out
}

all_out <- dplyr::bind_rows(
  run_country("ghana",  "Ghana"),
  run_country("malawi", "Malawi"))
write.csv(all_out, file.path(OUT, "metrics.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUT, "metrics.csv"), "\n")
