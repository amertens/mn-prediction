# =============================================================================
# scripts/covariates/15_continuous_vs_binary.R
#
# MODEL THE BIOMARKER, OR MODEL THE BINARY?
#
# Every district outcome in this project is a PREVALENCE: the share of people
# below a cutoff. Dichotomising throws information away, and it throws away
# most of it when the outcome is rare -- Malawi women's vitamin A is 11 events
# in 836. The alternative is to model the district MEAN of the continuous
# biomarker and convert to prevalence afterwards.
#
#   binary arm      y_j = logit(p_j)                     -> p_hat directly
#   continuous arm  y_j = mean_j( log biomarker )        -> p_hat = Phi((c - m_hat)/sigma)
#
# WHY A CONVERSION CEILING ARM IS INCLUDED
# ----------------------------------------
# The continuous arm can fail for two quite different reasons: the model may
# predict the district mean badly, or the mean-to-prevalence conversion may be
# a poor approximation even with a perfect mean. Arm 4 feeds the TRUE district
# mean through the same conversion, so it measures the second in isolation. If
# arm 4 is already worse than the binary arm, no improvement in prediction can
# rescue the continuous approach, and that is worth knowing before anyone
# rebuilds the pipeline around it.
#
# THE ARMS
#   1 binary (direct)          the status quo
#   2 continuous -> prevalence conversion with the POOLED training-district SD
#   3 continuous, oracle SD    conversion with the held-out district's OWN SD;
#                              uses held-out information, so an UPPER BOUND
#   4 conversion ceiling       TRUE district mean + pooled SD, no prediction
#   5 null                     training-mean prevalence
#
# Scored on: RMSE of the district mean on the modelling scale (arms 2-3 share
# one prediction, so one number), and prevalence MAE / RMSE / rank r for every
# arm, which is the comparison that actually decides the question.
#
#   Rscript scripts/covariates/15_continuous_vs_binary.R
# -> results/tables/continuous_vs_binary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")
K_SCREEN <- 20L
MIN_AREAS <- 12L

num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  m <- wmean(x, w); v <- sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
  sqrt(v * sum(ok) / max(1, sum(ok) - 1)) }

rows <- list(); cont_rows <- list()
cfgs <- get_country_configs()

for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  cov <- tryCatch(targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE),
                  error = function(e) NULL)
  if (is.list(cov) && !is.data.frame(cov) && "gee_admin2" %in% names(cov))
    cov <- cov$gee_admin2
  if (is.null(cov)) next

  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE),
                   error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    need <- c(oc$continuous, oc$binary, "Admin2")
    if (!all(need %in% names(d))) next

    y   <- num(d[[oc$continuous]])
    def <- num(d[[oc$binary]])
    w   <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d))
             num(d[[cc$weight_col]]) else rep(1, nrow(d))
    w[!is.finite(w) | w <= 0] <- 1

    # Modelling scale. Biomarker concentrations are right-skewed, so a normal
    # approximation only makes sense in logs. Gambia's iron marker is ALREADY
    # log (cutoff_scale == "log"), so it must not be logged twice.
    if (identical(oc$cutoff_scale, "log")) {
      t <- y; cut_t <- oc$cutoff
    } else {
      if (any(y <= 0, na.rm = TRUE)) next
      t <- log(y); cut_t <- log(oc$cutoff)
    }
    keep <- is.finite(t) & is.finite(def)
    if (sum(keep) < 100) next
    dd <- data.frame(Admin1 = if ("Admin1" %in% names(d)) as.character(d$Admin1) else NA_character_,
                     Admin2 = as.character(d$Admin2), t = t, def = def, w = w,
                     stringsAsFactors = FALSE)[keep, ]

    agg <- dd %>% group_by(Admin1, Admin2) %>%
      summarise(n = dplyr::n(), m = wmean(t, w), s = wsd(t, w),
                p = wmean(def, w), .groups = "drop") %>%
      filter(is.finite(m), is.finite(p), n >= 5)
    if (nrow(agg) < MIN_AREAS) next

    by <- admin2_join_by(agg, cov)
    fr <- dplyr::inner_join(agg, cov, by = by)
    fr <- fr[is.finite(fr$m), , drop = FALSE]
    if (nrow(fr) < MIN_AREAS) next

    covs <- setdiff(names(fr), c("Admin1", "Admin2", "n", "m", "s", "p"))
    covs <- covs[vapply(covs, function(v) is.numeric(fr[[v]]), logical(1))]
    X <- as.matrix(fr[, covs, drop = FALSE]); X[!is.finite(X)] <- NA
    for (j in seq_len(ncol(X))) { mu <- stats::median(X[, j], na.rm = TRUE)
      X[!is.finite(X[, j]), j] <- if (is.finite(mu)) mu else 0 }
    X <- X[, apply(X, 2, function(z) stats::sd(z) > 0), drop = FALSE]
    if (ncol(X) < 2) next

    folds <- if (sum(!is.na(fr$Admin1)) == nrow(fr) &&
                 dplyr::n_distinct(fr$Admin1) >= 3) as.character(fr$Admin1) else
             as.character(sample(rep_len(seq_len(min(10, nrow(fr))), nrow(fr))))
    wt <- area_precision_weights(fr$p, fr$n)

    oof_p_bin <- oof_m <- rep(NA_real_, nrow(fr))
    sigma_hat <- rep(NA_real_, nrow(fr))
    for (f in unique(folds)) {
      i <- which(folds == f); tr <- setdiff(seq_len(nrow(fr)), i)
      if (!length(tr) || length(i) >= nrow(fr)) next
      fit <- function(yv) {
        sel <- .awsl_screen(X[tr, , drop = FALSE], yv[tr], K_SCREEN)
        s <- tryCatch(.awsl_stack(X[tr, sel, drop = FALSE], yv[tr], wt[tr]),
                      error = function(e) NULL)
        if (is.null(s)) return(rep(NA_real_, length(i)))
        .awsl_predict(s, X[i, sel, drop = FALSE])
      }
      oof_p_bin[i] <- stats::plogis(fit(stats::qlogis(pmin(pmax(fr$p, .002), .998))))
      oof_m[i]     <- fit(fr$m)
      # pooled within-district SD from TRAINING districts only
      sigma_hat[i] <- stats::weighted.mean(fr$s[tr], fr$n[tr], na.rm = TRUE)
    }

    conv <- function(mhat, sig) stats::pnorm((cut_t - mhat) / sig)
    p_bin    <- oof_p_bin
    p_cont   <- conv(oof_m, sigma_hat)
    p_oracle <- conv(oof_m, fr$s)
    p_ceil   <- conv(fr$m, sigma_hat)          # true mean, no prediction
    p_null   <- rep(stats::weighted.mean(fr$p, fr$n), nrow(fr))

    met <- function(p, arm) { ok <- is.finite(p) & is.finite(fr$p)
      data.frame(country = cn, outcome = on, arm = arm, n_areas = sum(ok),
        prev_pct = round(100 * stats::weighted.mean(fr$p, fr$n), 2),
        mae_pp  = round(100 * mean(abs(fr$p[ok] - p[ok])), 2),
        rmse_pp = round(100 * sqrt(mean((fr$p[ok] - p[ok])^2)), 2),
        r = if (sum(ok) > 3 && stats::sd(p[ok]) > 0)
              round(suppressWarnings(stats::cor(fr$p[ok], p[ok])), 3) else NA_real_,
        stringsAsFactors = FALSE) }

    rows[[length(rows) + 1L]] <- rbind(
      met(p_bin,    "1 binary (direct)"),
      met(p_cont,   "2 continuous -> prevalence"),
      met(p_oracle, "3 continuous, oracle SD"),
      met(p_ceil,   "4 conversion ceiling (true mean)"),
      met(p_null,   "5 null (train mean)"))

    ok <- is.finite(oof_m)
    cont_rows[[length(cont_rows) + 1L]] <- data.frame(
      country = cn, outcome = on, n_areas = sum(ok),
      sd_of_mean = round(stats::sd(fr$m), 3),
      rmse_mean = round(sqrt(mean((fr$m[ok] - oof_m[ok])^2)), 3),
      r_mean = if (sum(ok) > 3) round(suppressWarnings(stats::cor(fr$m[ok], oof_m[ok])), 3) else NA_real_,
      pooled_sd = round(mean(sigma_hat, na.rm = TRUE), 3), stringsAsFactors = FALSE)
    cat(sprintf("  %-12s %-13s n=%3d prev=%5.1f%% | bin MAE=%5.2f r=%+.2f | cont MAE=%5.2f r=%+.2f | ceil MAE=%5.2f\n",
        cn, on, nrow(fr), 100*stats::weighted.mean(fr$p, fr$n),
        met(p_bin,"")$mae_pp, met(p_bin,"")$r, met(p_cont,"")$mae_pp,
        met(p_cont,"")$r, met(p_ceil,"")$mae_pp))
  }
}

res <- dplyr::bind_rows(rows); cont <- dplyr::bind_rows(cont_rows)
readr::write_csv(res, here("results", "tables", "continuous_vs_binary.csv"))
readr::write_csv(cont, here("results", "tables", "continuous_mean_rmse.csv"))

cat("\n================ PREVALENCE ACCURACY BY ARM ==========================\n")
print(res %>% group_by(arm) %>%
        summarise(cells = n(), mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                  med_mae = round(stats::median(mae_pp, na.rm = TRUE), 2),
                  mean_rmse = round(mean(rmse_pp, na.rm = TRUE), 2),
                  mean_r = round(mean(r, na.rm = TRUE), 3), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

w <- res %>% select(country, outcome, arm, mae_pp) %>%
  tidyr::pivot_wider(names_from = arm, values_from = mae_pp)
cmp <- function(a, b, lab) cat(sprintf("%-46s %d/%d cells\n", lab,
  sum(w[[a]] < w[[b]], na.rm = TRUE), sum(is.finite(w[[a]]) & is.finite(w[[b]]))))
cat("\n")
cmp("2 continuous -> prevalence", "1 binary (direct)", "continuous beats binary (MAE)")
cmp("2 continuous -> prevalence", "5 null (train mean)", "continuous beats null")
cmp("1 binary (direct)", "5 null (train mean)", "binary beats null")
cmp("4 conversion ceiling (true mean)", "1 binary (direct)",
    "conversion ceiling beats binary  <- is the ceiling even high enough?")

cat("\n================ DOES RARITY CHANGE THE ANSWER? ======================\n")
rr <- res %>% filter(arm %in% c("1 binary (direct)", "2 continuous -> prevalence")) %>%
  mutate(band = cut(prev_pct, c(-Inf, 5, 20, Inf),
                    labels = c("rare (<5%)", "mid (5-20%)", "common (>20%)")))
print(rr %>% group_by(band, arm) %>%
        summarise(cells = n(), mae = round(mean(mae_pp), 2),
                  r = round(mean(r, na.rm = TRUE), 3), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\n================ DISTRICT MEAN, CONTINUOUS SCALE =====================\n")
print(as.data.frame(cont), row.names = FALSE)
cat("\n-> results/tables/continuous_vs_binary.csv\n")
cat("-> results/tables/continuous_mean_rmse.csv\n")
