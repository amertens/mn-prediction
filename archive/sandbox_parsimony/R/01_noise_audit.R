# =============================================================================
# sandbox_parsimony/R/01_noise_audit.R
#
# How much of the between-Admin-2 variance in survey prevalence is REAL signal
# and how much is survey sampling noise? This sets the ceiling that any
# GEE-only or LOCO model can possibly hit.
#
# Method-of-moments decomposition per country x outcome:
#   V_obs  = var(p_i)                       observed spread of area prevalences
#   V_samp = mean(sampling_var_unbiased(p_i, n_i, deff))  sampling variance
#   lambda = (V_obs - V_samp) / V_obs       reliability
#   r_max  = sqrt(lambda)                   ceiling on Pearson r vs ANY predictor
#
# NOTE ON deff. The stored `svy_prev_se` cannot be used to estimate the design
# effect: most Admin-2 areas contain a single PSU, so the between-cluster
# variance estimate collapses and survey::svymean returns SE ~ 1e-17. We
# therefore sweep deff over a plausible range for DHS/MICS biomarker
# subsamples and report the reliability as a band.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
source("sandbox_parsimony/R/03_core.R")   # sampling_var()

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
DEFFS <- c(1.0, 1.5, 2.0)

rel_at <- function(p, n, deff) {
  v_obs  <- stats::var(p, na.rm = TRUE)
  v_samp <- mean(sampling_var_unbiased(p, n, deff), na.rm = TRUE)
  max(0, v_obs - v_samp) / v_obs
}

audit <- list()
for (oc in names(pooled_all)) {
  d <- pooled_all[[oc]]$data
  for (ctry in unique(d$country)) {
    s <- d[d$country == ctry, ]
    p <- s$svy_prev; n <- s$n_svy; se <- s$svy_prev_se
    ok <- is.finite(p) & is.finite(n)
    if (sum(ok) < 5) next
    p <- p[ok]; n <- n[ok]; se <- se[ok]

    lam <- vapply(DEFFS, function(dd) rel_at(p, n, dd), numeric(1))
    # lambda is a difference of two variances estimated from 14-168 districts,
    # so it carries real uncertainty. Bootstrap districts to show how much.
    set.seed(4242L)
    bs <- replicate(500L, {
      i <- sample.int(length(p), replace = TRUE)
      rel_at(p[i], n[i], 1.5)
    })
    bs <- bs[is.finite(bs)]                    # a resample can be degenerate
    if (!length(bs)) bs <- NA_real_
    audit[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, n_areas = length(p),
      median_n = stats::median(n), min_n = min(n),
      pct_n_lt_25 = round(100 * mean(n < 25), 0),
      mean_prev_pp = round(100 * mean(p), 1),
      sd_obs_pp = round(100 * stats::sd(p), 2),
      sd_samp_pp_d15 = round(100 * sqrt(mean(sampling_var_unbiased(p, n, 1.5))), 2),
      pct_se_degenerate = round(100 * mean(!is.finite(se) | se < 1e-6), 0),
      lambda_d10 = round(lam[1], 3), lambda_d15 = round(lam[2], 3),
      lambda_d20 = round(lam[3], 3),
      lambda_d15_lo = round(stats::quantile(bs, .025, names = FALSE, na.rm = TRUE), 3),
      lambda_d15_hi = round(stats::quantile(bs, .975, names = FALSE, na.rm = TRUE), 3),
      r_max_d15 = round(sqrt(lam[2]), 3),
      r_max_d15_hi = round(sqrt(stats::quantile(bs, .975, names = FALSE, na.rm = TRUE)), 3),
      stringsAsFactors = FALSE)
  }
}
audit <- bind_rows(audit) |> arrange(outcome, country)
write.csv(audit, "sandbox_parsimony/out/noise_audit.csv", row.names = FALSE)

cat("\n=== Admin-2 outcome reliability (deff sweep) ===\n")
cat("lambda = share of between-district variance that is real signal.\n")
cat("r_max  = sqrt(lambda) = best Pearson r ANY model could achieve.\n")
cat("pct_se_degenerate = %% of districts where the stored design SE collapsed\n")
cat("                    to ~0 because the district holds a single PSU.\n\n")
print(audit |> mutate(cell = paste(outcome, country)) |>
        select(cell, median_n, pct_n_lt_25, mean_prev_pp, sd_obs_pp,
               sd_samp_pp_d15, pct_se_degenerate, lambda_d10, lambda_d15,
               lambda_d20, lambda_d15_lo, lambda_d15_hi, r_max_d15),
      row.names = FALSE)

# ---- published within-country CV r vs. the ceiling -------------------------
ceil_f <- "results/tables/transportability_within_country_ceiling.csv"
if (file.exists(ceil_f)) {
  ceil <- read.csv(ceil_f, stringsAsFactors = FALSE)
  cmp <- ceil |>
    mutate(country = tolower(gsub(" ", "", country))) |>
    inner_join(audit |> select(outcome, country, lambda_d15, r_max_d15,
                              r_max_d15_hi),
               by = c("outcome", "country")) |>
    mutate(pct_of_ceiling = ifelse(r_max_d15 > .1,
                                   round(100 * cv_pearson / r_max_d15, 0), NA),
           pct_of_ceiling_optimistic = ifelse(r_max_d15_hi > .1,
                                   round(100 * cv_pearson / r_max_d15_hi, 0), NA))
  cat("\n=== Published within-country CV r vs. the noise ceiling ===\n")
  print(cmp |> select(outcome, country, n, cv_pearson, r_max_d15, r_max_d15_hi,
                      pct_of_ceiling, pct_of_ceiling_optimistic), row.names = FALSE)
  write.csv(cmp, "sandbox_parsimony/out/ceiling_vs_published.csv", row.names = FALSE)
}

# ---- what the tiny district samples do to a Fay-Herriot shrinkage weight ---
cat("\n=== Consequence for Fay-Herriot (dashboard/data-raw/_build_fh_layer.R) ===\n")
cat("FH shrinks the direct estimate by gamma = A/(A+D), D = sampling variance.\n")
cat("That file takes D from svy_prev_se^2, treats only D<=0 as bad, then floors\n")
cat("at 1e-5. A degenerate SE of 1e-17 therefore survives as D = 1e-5\n")
cat("(SE = 0.32 pp) instead of the true p(1-p)/n.\n\n")
fh <- list()
for (oc in names(pooled_all)) {
  d <- pooled_all[[oc]]$data
  for (ctry in unique(d$country)) {
    s <- d[d$country == ctry, ]
    ok <- is.finite(s$svy_prev) & is.finite(s$n_svy); s <- s[ok, ]
    if (nrow(s) < 5) next
    D_used  <- pmax(ifelse(is.finite(s$svy_prev_se) & s$svy_prev_se^2 > 0,
                           s$svy_prev_se^2, s$svy_prev * (1 - s$svy_prev) /
                             pmax(s$n_svy / 1.5, 1)), 1e-5)
    D_true  <- sampling_var(s$svy_prev, s$n_svy, 1.5)  # per-district: AC-shrunk
    A <- max(stats::var(s$svy_prev) -
               mean(sampling_var_unbiased(s$svy_prev, s$n_svy, 1.5)), 1e-4)
    fh[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, n_areas = nrow(s),
      mean_gamma_as_coded = round(mean(A / (A + D_used)), 3),
      mean_gamma_correct  = round(mean(A / (A + D_true)), 3),
      pct_gamma_gt_95_as_coded = round(100 * mean(A / (A + D_used) > .95), 0),
      pct_gamma_gt_95_correct  = round(100 * mean(A / (A + D_true) > .95), 0),
      stringsAsFactors = FALSE)
  }
}
fh <- bind_rows(fh) |> arrange(outcome, country)
print(fh, row.names = FALSE)
write.csv(fh, "sandbox_parsimony/out/fh_shrinkage_check.csv", row.names = FALSE)
cat("\ngamma = 1 means 'no shrinkage, keep the raw direct estimate'.\n")
