# =============================================================================
# scripts/accuracy_impact/wsl1_within_country.R
#
# WS-L. The consolidated within-country run: PREVALENCE and BIOMARKER LEVEL,
# one protocol, one predictor set, one table.
#
# WHY THIS EXISTS
# ---------------
# Within-country prediction was scattered across three tables using three
# different nulls and three different fold schemes, which is how the same
# project came to report both "covariates beat the null in 21 of 24 cells" and
# "covariates lose to the null in every arm". Both were true of their own
# table. Neither was labelled well enough to reconcile.
#
# This script fixes all three choices and reports both outcome types together:
#
#   FOLD    region-blocked. Holding out a district while keeping its neighbours
#           in training reads 0.416 where region-blocking reads 0.154, on
#           identical data, so the permissive version is not reported at all.
#   NULL    the mean of the TRAINING areas only. The full-data national mean
#           contains the held-out district's own respondents and is therefore
#           contaminated with the answer; scoring against it swings the
#           conclusion by 2.80 pp.
#   SET     data/covariates/harmonized/predictors_admin2_shared.csv, with the
#           coverage and domain metadata beside it.
#
# CONTINUOUS OUTCOMES ARE NOT INTERCHANGEABLE ACROSS COUNTRIES
# ------------------------------------------------------------
# All 24 cells define one, but they are not on a common footing, and three
# differences were measured rather than assumed:
#
#   SCALE      Gambia's iron variable is a LOG ferritin (median 2.14) where
#              Ghana, Sierra Leone and Malawi are linear ug/L (median 24-48).
#              Neither exp() nor 10^ reconciles it against the others, so no
#              transformation is applied and cross-country pooling of the
#              continuous iron outcome is refused rather than guessed.
#   CENSORING  Ghana and Sierra Leone share assay ceilings exactly -- 1476 for
#              B12 and 45.40 for folate, with 7 women at the folate cap in each.
#              Malawi is uncapped. So the same nutrient is right-censored in two
#              countries and not in a third.
#   SENTINEL   Malawi child zinc contains one value of -100, which is not a
#              concentration. Dropped, and the count reported.
#
# Within a country each definition is self-consistent, so within-country
# modelling is sound. It is the pooling that is not, and this script scores
# within country only.
#
#   Rscript scripts/accuracy_impact/wsl1_within_country.R
# -> results/tables/within_country_consolidated.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")
SEED <- 20260929L
MIN_AREAS <- 10L
HOLDOUT_FRAC <- 0.30   # share of regions held out per fold

SHARED <- here("data", "covariates", "harmonized", "predictors_admin2_shared.csv")
P <- suppressMessages(readr::read_csv(
  if (file.exists(SHARED)) SHARED else
    here("data","covariates","harmonized","predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(P), c("country", "Admin1", "Admin2"))
cat(sprintf("[wsl1] predictor set: %s, %d predictors\n", basename(SHARED), length(COVS)))
cfgs <- get_country_configs()
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# ---------------------------------------------------------------------------
# Area-level outcome: prevalence (binary) and mean level (continuous), both
# survey-weighted, from the same individual rows so the two are paired.
# ---------------------------------------------------------------------------
area_outcomes <- function(d, cc, oc) {
  a1 <- trimws(as.character(d[[cc$admin1_col]]))
  a2 <- trimws(as.character(d[[cc$admin2_col]]))
  w  <- num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- 1
  key <- paste(a1, a2, sep = "||")
  out <- data.frame(Admin1 = tapply(a1, key, function(z) z[1]),
                    Admin2 = tapply(a2, key, function(z) z[1]),
                    stringsAsFactors = FALSE)
  wm <- function(v) {
    x <- num(v)
    tapply(seq_along(x), key, function(i) {
      ok <- is.finite(x[i]) & is.finite(w[i])
      if (!any(ok)) NA_real_ else stats::weighted.mean(x[i][ok], w[i][ok])
    })
  }
  nn <- function(v) { x <- num(v); tapply(seq_along(x), key, function(i) sum(is.finite(x[i]))) }
  if (!is.null(oc$binary) && oc$binary %in% names(d)) {
    out$prev <- as.numeric(wm(d[[oc$binary]])); out$n_prev <- as.integer(nn(d[[oc$binary]]))
  }
  if (!is.null(oc$continuous) && oc$continuous %in% names(d)) {
    x <- num(d[[oc$continuous]])
    # A biomarker concentration cannot be negative or zero. Malawi zinc carries
    # one -100; dropping it here rather than letting it enter a mean.
    bad <- is.finite(x) & x <= 0
    if (any(bad)) {
      cat(sprintf("      [clean] %d non-positive %s values dropped\n",
                  sum(bad), oc$continuous))
      x[bad] <- NA_real_
      d[[oc$continuous]] <- x
    }
    out$level <- as.numeric(wm(d[[oc$continuous]])); out$n_level <- as.integer(nn(d[[oc$continuous]]))
  }
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Region-blocked folds, and the training-mean null computed inside each fold.
# ---------------------------------------------------------------------------
score_cell <- function(m, y, n_y, covs, label, link = "logit") {
  ok <- is.finite(y) & is.finite(m$.reg_i)
  if (sum(ok) < MIN_AREAS) return(NULL)
  regs <- unique(m$Admin1[ok])
  if (length(regs) < 3) return(NULL)
  set.seed(SEED)
  k <- max(2L, round(1 / HOLDOUT_FRAC))
  assign_r <- split(sample(regs), rep(seq_len(k), length.out = length(regs)))
  pred <- rep(NA_real_, nrow(m)); nullp <- rep(NA_real_, nrow(m))
  X <- as.matrix(m[, covs, drop = FALSE])
  for (f in seq_len(k)) {
    te <- which(m$Admin1 %in% assign_r[[f]] & ok)
    tr <- which(!(m$Admin1 %in% assign_r[[f]]) & ok)
    if (length(te) < 1 || length(tr) < 8) next
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      Xtr[!is.finite(Xtr[, j]), j] <- mu; Xte[!is.finite(Xte[, j]), j] <- mu
    }
    keep <- apply(Xtr, 2, function(z) stats::sd(z) > 0)
    if (sum(keep) < 3) next
    # THE NULL: mean of the TRAINING areas only, weighted by their sample size.
    nullp[te] <- stats::weighted.mean(y[tr], pmax(n_y[tr], 1), na.rm = TRUE)
    f1 <- tryCatch(.ds_fit(Xtr[, keep, drop = FALSE], y[tr],
                           k_screen = min(20L, sum(keep)), link = link),
                   error = function(e) NULL)
    if (is.null(f1)) next
    pp <- tryCatch(.ds_predict(f1, Xte[, keep, drop = FALSE]), error = function(e) NULL)
    if (!is.null(pp)) pred[te] <- pp
  }
  o <- is.finite(y) & is.finite(pred) & is.finite(nullp)
  if (sum(o) < MIN_AREAS) return(NULL)
  sc <- function(p) {
    r <- if (stats::sd(p[o]) > 0 && stats::sd(y[o]) > 0)
      suppressWarnings(stats::cor(p[o], y[o])) else NA_real_
    c(r = r, mae = mean(abs(p[o] - y[o])), rmse = sqrt(mean((p[o] - y[o])^2)))
  }
  a <- sc(pred); b <- sc(nullp)
  # Targeting lift and its ceiling, on the same held-out predictions.
  kk <- max(1L, round(0.20 * sum(o)))
  yv <- y[o]; ov <- mean(yv)
  lift <- if (is.finite(ov) && ov > 0)
    mean(yv[order(pred[o], decreasing = TRUE)[seq_len(kk)]]) / ov else NA_real_
  orc <- if (is.finite(ov) && ov > 0)
    mean(yv[order(yv, decreasing = TRUE)[seq_len(kk)]]) / ov else NA_real_
  data.frame(target = label, n_areas = sum(o),
             r = round(a[["r"]], 4), mae = round(a[["mae"]], 5),
             null_mae = round(b[["mae"]], 5),
             mae_vs_null = round(a[["mae"]] - b[["mae"]], 5),
             beats_null = a[["mae"]] < b[["mae"]],
             lift = round(lift, 3), lift_oracle = round(orc, 3),
             lift_share = round(if (is.finite(orc) && orc > 1) (lift - 1)/(orc - 1) else NA_real_, 3),
             stringsAsFactors = FALSE)
}

rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- P[P$country == cn, , drop = FALSE]
  if (!nrow(hc)) next
  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    need <- c(cc$admin1_col, cc$admin2_col, cc$weight_col)
    if (!all(need %in% names(d))) next
    cat(sprintf("  %-13s %-13s\n", cn, ocn))
    ao <- area_outcomes(d, cc, oc)
    m <- dplyr::inner_join(ao, hc, by = admin2_join_by(ao, hc))
    if (nrow(m) < MIN_AREAS) next
    m$.reg_i <- as.integer(factor(m$Admin1))
    covs <- COVS[vapply(COVS, function(v)
      mean(is.finite(m[[v]])) > 0.5 && stats::sd(m[[v]], na.rm = TRUE) > 0, logical(1))]
    if (length(covs) < 5) next
    for (tgt in c("prev", "level")) {
      if (!tgt %in% names(m)) next
      nc <- if (tgt == "prev") "n_prev" else "n_level"
      # Prevalences are proportions; biomarker levels are positive and
      # right-skewed (measured skew 0.29 to 9.28), so they take a log link.
      r <- score_cell(m, m[[tgt]], m[[nc]], covs, tgt,
                      link = if (tgt == "prev") "logit" else "log")
      if (is.null(r)) next
      r$country <- cc$country; r$outcome <- ocn; r$n_cov <- length(covs)
      r$variable <- if (tgt == "prev") oc$binary %||% NA_character_ else
        oc$continuous %||% NA_character_
      rows[[length(rows) + 1L]] <- r
    }
  }
}

if (!length(rows)) stop("[wsl1] nothing scored")
out <- dplyr::bind_rows(rows)
front <- c("country", "outcome", "target", "variable")
out <- out[, c(front, setdiff(names(out), front))]
readr::write_csv(out, file.path(TDIR, "within_country_consolidated.csv"))

cat("\n=== within country, region-blocked folds, training-mean null ===\n")
s <- out |> group_by(target) |>
  summarise(cells = dplyr::n(),
            median_r = round(stats::median(r, na.rm = TRUE), 3),
            beats_null = sprintf("%d/%d", sum(beats_null, na.rm = TRUE), dplyr::n()),
            median_lift = round(stats::median(lift, na.rm = TRUE), 2),
            median_lift_share = round(stats::median(lift_share, na.rm = TRUE), 2),
            .groups = "drop")
print(as.data.frame(s), row.names = FALSE)

cat("\n=== per cell ===\n")
print(as.data.frame(out[, c("country","outcome","target","n_areas","r",
                            "mae","null_mae","beats_null","lift","lift_oracle")]),
      row.names = FALSE)
cat(sprintf("\n-> %s\n", file.path("results","tables","within_country_consolidated.csv")))
