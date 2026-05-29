# =============================================================================
# scripts/run_calibration_recalibration.R
#
# Applies Platt (logistic) recalibration to the out-of-fold (OOF) binary
# predictions stored in each `sl_fit_<country>_<outcome>` target. Repairs
# the negative-Brier-skill / inverted-calibration-slope failures
# identified in diagnostics_binary.csv (Malawi child_vitA, Ghana women_vitA,
# Sierra Leone women_b12, etc.).
#
# Outputs:
#   results/tables/diagnostics_binary_calibrated.csv
#       Per-(country, outcome) row with both RAW and CALIBRATED diagnostics
#       side-by-side. Platt calibration is fit two ways:
#         "in_sample"  — fit on all OOF preds, apply to same. Simple but
#                        mildly optimistic (2 params estimated from same
#                        data they're calibrating — usually ≤0.001 AUC
#                        bias at N≥200 per outcome).
#         "cv5"        — 5-fold CV: fit calibrator on 4/5 of OOF preds,
#                        apply to the held-out 5th. Honest, slightly
#                        noisier.
#   results/tables/calibrators.rds
#       Named list of Platt coefficient pairs (intercept, slope on
#       logit-of-raw-prob) per (country, outcome). Use these to
#       recalibrate predictions on NEW data:
#           p_cal <- plogis(intercept + slope * qlogis(p_raw))
#
# Reads sl_fit_<country>_<outcome> objects from the _targets_full store
# WITHOUT tar_make() to avoid the package-version cascade.
# =============================================================================

suppressMessages({
  library(targets); library(dplyr); library(here)
})

source(here::here("R", "benchmark_models.R"))
source(here::here("R", "diagnostics.R"))

STORE <- here::here("_targets_full")
OUT_DIR <- here::here("results", "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Country / outcome inventory ──────────────────────────────────────────
COUNTRIES <- c(Gambia = "gambia", Ghana = "ghana",
                SierraLeone = "sierraleone", Malawi = "malawi")
OUTCOMES <- list(
  Gambia      = c("child_vitA","women_vitA","child_iron","women_iron"),
  Ghana       = c("child_vitA","women_vitA","child_iron","women_iron",
                   "women_b12","women_folate"),
  SierraLeone = c("child_vitA","women_vitA","child_iron","women_iron",
                   "women_b12","women_folate"),
  Malawi      = c("child_vitA","women_vitA","child_iron","women_iron",
                   "women_b12","women_folate","child_zinc","women_zinc")
)

# ── Helper: fit Platt calibrator on logit(p_raw) → y ─────────────────────
.platt_fit <- function(p_raw, y, eps = 1e-3) {
  ok <- is.finite(p_raw) & is.finite(y) & y %in% c(0, 1)
  p_raw <- pmin(pmax(p_raw[ok], eps), 1 - eps); y <- y[ok]
  if (length(p_raw) < 10 || length(unique(y)) < 2) return(NULL)
  fit <- stats::glm(y ~ I(log(p_raw / (1 - p_raw))),
                     family = stats::binomial())
  list(intercept = as.numeric(stats::coef(fit)[1]),
        slope     = as.numeric(stats::coef(fit)[2]))
}

.platt_apply <- function(p_raw, fit, eps = 1e-3) {
  if (is.null(fit)) return(p_raw)
  p_raw <- pmin(pmax(p_raw, eps), 1 - eps)
  stats::plogis(fit$intercept + fit$slope * log(p_raw / (1 - p_raw)))
}

# ── Helper: 5-fold CV calibration (honest) ───────────────────────────────
.platt_cv5 <- function(p_raw, y, K = 5L, seed = 12345L) {
  set.seed(seed)
  n <- length(y); folds <- sample(rep(seq_len(K), length.out = n))
  out <- rep(NA_real_, n)
  for (k in seq_len(K)) {
    tr <- folds != k; te <- folds == k
    fit <- .platt_fit(p_raw[tr], y[tr])
    if (is.null(fit)) out[te] <- p_raw[te]
    else              out[te] <- .platt_apply(p_raw[te], fit)
  }
  out
}

# ── Helper: compute binary diagnostics (mirrors R/diagnostics.R) ─────────
.bin_metrics <- function(Y, yhat) {
  if (length(Y) < 5) return(NULL)
  yhat <- pmin(pmax(yhat, 1e-6), 1 - 1e-6)
  prev <- mean(Y == 1, na.rm = TRUE)
  brier <- mean((yhat - Y)^2, na.rm = TRUE)
  # No-skill Brier = prevalence * (1 - prevalence) (Bernoulli variance).
  brier_skill <- 1 - brier / (prev * (1 - prev))
  # AUC via Mann-Whitney U
  pos <- yhat[Y == 1]; neg <- yhat[Y == 0]
  auc <- tryCatch({
    if (length(pos) < 1 || length(neg) < 1) NA_real_
    else {
      # Mean over all positive-negative pairs.
      mean(outer(pos, neg, FUN = ">") +
            0.5 * outer(pos, neg, FUN = "==") )
    }
  }, error = function(e) NA_real_)
  # Calibration slope from logistic recalibration check
  cal <- tryCatch({
    logit_yhat <- log(yhat / (1 - yhat))
    f <- stats::glm(Y ~ logit_yhat, family = stats::binomial())
    list(intercept = as.numeric(stats::coef(f)[1]),
          slope     = as.numeric(stats::coef(f)[2]))
  }, error = function(e) list(intercept = NA_real_, slope = NA_real_))
  data.frame(
    n        = length(Y),
    n_events = sum(Y == 1, na.rm = TRUE),
    prevalence  = round(prev, 4),
    auc         = round(auc, 4),
    brier       = round(brier, 4),
    brier_skill = round(brier_skill, 4),
    calib_intercept = round(cal$intercept, 4),
    calib_slope     = round(cal$slope, 4)
  )
}

# ── Main loop ────────────────────────────────────────────────────────────
rows <- list()
calibrators <- list()
t0 <- Sys.time()

for (ctry_disp in names(COUNTRIES)) {
  ctry <- COUNTRIES[[ctry_disp]]
  for (oc in OUTCOMES[[ctry_disp]]) {
    tgt <- paste0("sl_fit_", ctry, "_", oc)
    sl <- tryCatch(targets::tar_read_raw(tgt, store = STORE),
                    error = function(e) NULL)
    if (is.null(sl) || is.null(sl$bin_fit) || is.null(sl$bin_fit$res)) next
    res <- sl$bin_fit$res
    Y  <- as.numeric(res$Y)
    p_raw <- as.numeric(res$yhat_full)
    if (length(Y) < 5) next
    cat(sprintf("[cal] %s %s (n=%d) ...", ctry_disp, oc, length(Y)))

    # Raw diagnostics
    raw <- .bin_metrics(Y, p_raw)
    # In-sample Platt
    fit_in <- .platt_fit(p_raw, Y)
    p_cal_in <- .platt_apply(p_raw, fit_in)
    cal_in   <- .bin_metrics(Y, p_cal_in)
    # 5-fold CV Platt (honest)
    p_cal_cv <- .platt_cv5(p_raw, Y)
    cal_cv   <- .bin_metrics(Y, p_cal_cv)

    row <- data.frame(
      country = ctry_disp, outcome = oc,
      n = raw$n, n_events = raw$n_events, prevalence = raw$prevalence,
      # raw
      auc_raw = raw$auc, brier_raw = raw$brier,
      brier_skill_raw = raw$brier_skill,
      calib_slope_raw = raw$calib_slope,
      calib_intercept_raw = raw$calib_intercept,
      # in-sample Platt
      auc_platt_insample = cal_in$auc, brier_platt_insample = cal_in$brier,
      brier_skill_platt_insample = cal_in$brier_skill,
      calib_slope_platt_insample = cal_in$calib_slope,
      # 5-fold CV Platt
      auc_platt_cv5 = cal_cv$auc, brier_platt_cv5 = cal_cv$brier,
      brier_skill_platt_cv5 = cal_cv$brier_skill,
      calib_slope_platt_cv5 = cal_cv$calib_slope,
      # Platt coefficients (in-sample fit)
      platt_intercept = if (!is.null(fit_in)) fit_in$intercept else NA_real_,
      platt_slope     = if (!is.null(fit_in)) fit_in$slope     else NA_real_,
      stringsAsFactors = FALSE
    )
    rows[[paste(ctry_disp, oc)]] <- row
    calibrators[[paste(ctry_disp, oc)]] <- fit_in

    cat(sprintf(" brier_skill %+.3f -> %+.3f  slope %+.2f -> %+.2f\n",
                raw$brier_skill, cal_cv$brier_skill,
                raw$calib_slope, cal_cv$calib_slope))
  }
}

result <- dplyr::bind_rows(rows)
csv_path <- file.path(OUT_DIR, "diagnostics_binary_calibrated.csv")
rds_path <- file.path(OUT_DIR, "calibrators.rds")
readr::write_csv(result, csv_path)
saveRDS(calibrators, rds_path)

# ── Side-by-side summary ─────────────────────────────────────────────────
cat("\n===== Brier skill score: raw vs calibrated (5-fold CV) =====\n")
summary_tbl <- result |>
  dplyr::transmute(
    country, outcome,
    raw       = round(brier_skill_raw,        3),
    insample  = round(brier_skill_platt_insample, 3),
    cv5       = round(brier_skill_platt_cv5,  3),
    `delta_cv5_minus_raw` = round(brier_skill_platt_cv5 - brier_skill_raw, 3)
  ) |>
  dplyr::arrange(dplyr::desc(`delta_cv5_minus_raw`))
print(as.data.frame(summary_tbl), row.names = FALSE)

cat("\n===== calibration slope: raw vs calibrated (5-fold CV) =====\n")
slope_tbl <- result |>
  dplyr::transmute(country, outcome,
                    raw = round(calib_slope_raw, 2),
                    cv5 = round(calib_slope_platt_cv5, 2)) |>
  dplyr::arrange(abs(raw - 1))   # furthest from ideal (slope = 1) first
print(as.data.frame(slope_tbl), row.names = FALSE)

cat(sprintf("\nwritten:\n  %s\n  %s\n", csv_path, rds_path))
cat(sprintf("\nTotal elapsed: %.1f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))
