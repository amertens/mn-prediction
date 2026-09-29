# =============================================================================
# scripts/protocol_v2/79_prevalence_noise_split.R   [NZ-01, 2026-09-28]
#
# HOW MUCH OF THE DISTRICT PREVALENCE ERROR IS THE SURVEY'S OWN NOISE?
#
# The district prevalence error is always scored against the survey's district
# figure, which is itself a small-sample estimate. Method of moments, per cell:
#
#   E[(p_svy - p_hat)^2] = E[(p_true - p_hat)^2] + Var(p_svy | p_true)
#
# because p_hat is out-of-fold (it never saw the district's own respondents),
# so the survey's sampling error is independent of the model's error.
#   observed MSE   mean over districts of (p_svy - p_hat)^2
#   survey noise   mean over districts of p_d (1 - p_d) / n_eff_d, plug-in p_d
#   true-error MSE max(MSE - noise, 0); where noise exceeds MSE it is floored
#                  at 0 and flagged, and the cell is kept
# n_eff (targets_v2.csv, script 01, V2_DEFF_METHOD = "district", the default):
#   n_eff_d = n_kish_d / (1 + (b_d - 1) rho), b_d = n_raw_d / n_psu_d (the
#   district's respondents per cluster), rho = the intra-cluster correlation
#   implied by the NATIONAL design effect of that country x outcome after
#   removing the weighting part (deff / (n_raw / n_kish)), floored at 1; the
#   national deff is the fallback where rho is not estimable. So it carries
#   both the weighting and the clustering design effect.
#
# Residuals: the post-fix rerun of script 66 (CP-01): the calibrated index
# (domain_index_cal), in-fill (5-fold by district, 10 draws), anchored on
# population to the national figure; p_hat is the mean of the out-of-fold
# predictions. The Malawi selenium / iodine cells take n_eff from their own
# table (script 61), as dashboard/data-raw/00_read_targets.R does (logic
# copied here, the dashboard file is not sourced).
#
# Headline, pre-specified: pooled over every district of every cell (each
# district counts once per outcome); the median over cells is reported too.
# Also written: the CP-01 rerun against its pre-fix snapshot.
#
#   Rscript -e "source('scripts/protocol_v2/79_prevalence_noise_split.R')"
# -> results/tables/protocol_v2/nz01_cells.csv              per-cell split
#    results/tables/protocol_v2/nz01_overall.csv            pooled and median-of-cells
#    results/tables/protocol_v2/nz01_cp01_rerun_compare.csv  CP-01 bands, pre-fix vs post-fix
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
P2 <- "results/tables/protocol_v2"
SNAP <- "C:/Users/andre/mn-prediction-snapshots/pre_rebuild_2026-09-27/results/tables/protocol_v2"
nk <- function(x) gsub(" ", "", x)

# ---- targets with n_eff (00_read_targets.R's logic, copied) ------------------
TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
se <- read.csv(file.path(P2, "malawi_selenium_iodine_targets.csv"), stringsAsFactors = FALSE)
se <- se[is.finite(se$p_low) & is.finite(se$n_eff) & se$n_eff > 0, ]
TGX <- bind_rows(TG[, c("country", "outcome", "Admin1", "Admin2", "y_prev", "n_eff", "n_raw", "deff_binary")],
                 data.frame(country = "Malawi", outcome = se$outcome, Admin1 = se$Admin1, Admin2 = se$Admin2,
                            y_prev = se$p_low, n_eff = se$n_eff, n_raw = se$n, deff_binary = NA_real_, stringsAsFactors = FALSE)) |>
  mutate(country = nk(country), Admin1 = trimws(Admin1), Admin2 = trimws(Admin2))
stopifnot(!anyDuplicated(TGX[, c("country", "outcome", "Admin1", "Admin2")]))

# ---- post-fix residuals from script 66 ---------------------------------------
D <- read.csv(file.path(P2, "conformal_prev_districts.csv"), stringsAsFactors = FALSE) |>
  mutate(ckey = nk(country), Admin1 = trimws(Admin1), Admin2 = trimws(Admin2))
M <- D |> left_join(TGX, by = c(ckey = "country", "outcome", "Admin1", "Admin2"), suffix = c("", ".tg"))
cat(sprintf("districts %d | joined to an n_eff %d | max |p_survey - y_prev| %.2e\n", nrow(M), sum(is.finite(M$n_eff)),
            max(abs(M$p_survey - M$y_prev), na.rm = TRUE)))
if (any(!is.finite(M$n_eff))) stop("NZ-01: districts without an n_eff; stopping")
if (max(abs(M$p_survey - M$y_prev)) > 1e-8) stop("NZ-01: conformal_prev_districts.csv p_survey does not match the targets; stopping")
M <- M |> mutate(sq = resid^2, noise = p_survey * (1 - p_survey) / n_eff)

# national figure script 66 anchored on vs the post-fix district data (diagnostic only)
NE <- read.csv("results/tables/national_estimates_all.csv", stringsAsFactors = FALSE) |> transmute(ckey = nk(country), outcome, p_nat_anchor_used = obs_prev)
PN <- TGX |> group_by(ckey = country, outcome) |> summarise(p_nat_postfix_nraw = stats::weighted.mean(y_prev, n_raw), .groups = "drop")

CELLS <- M |> group_by(country, outcome) |>
  summarise(n = dplyr::n(), mean_p = mean(p_survey), median_n_eff = median(n_eff), deff_national = dplyr::first(deff_binary),
            mean_resid_pp = 100 * mean(resid), mse = mean(sq), noise = mean(noise), .groups = "drop") |>
  mutate(ckey = nk(country), noise_share_raw = noise / mse, noise_exceeds_mse = noise > mse, noise_share = pmin(noise_share_raw, 1),
         rmse_obs_pp = 100 * sqrt(mse), rmse_noise_pp = 100 * sqrt(noise), rmse_true_pp = 100 * sqrt(pmax(mse - noise, 0))) |>
  left_join(NE, by = c("ckey", "outcome")) |> left_join(PN, by = c("ckey", "outcome")) |> select(-ckey)
write.csv(CELLS, file.path(P2, "nz01_cells.csv"), row.names = FALSE)

pool <- function(d, label) { mse <- mean(d$sq); nz <- mean(d$noise)
  data.frame(scope = label, cells = dplyr::n_distinct(paste(d$country, d$outcome)), districts = nrow(d), mse = mse, noise = nz,
             noise_share = min(nz / mse, 1), rmse_obs_pp = 100 * sqrt(mse), rmse_true_pp = 100 * sqrt(max(mse - nz, 0)), stringsAsFactors = FALSE) }
PANEL6 <- c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12")
OVR <- bind_rows(pool(M, "pooled, all cells"), pool(M[M$outcome %in% PANEL6, ], "pooled, the 22 panel cells (no zinc, selenium, iodine)"),
  data.frame(scope = "median over cells", cells = nrow(CELLS), districts = nrow(M), mse = median(CELLS$mse), noise = median(CELLS$noise),
             noise_share = median(CELLS$noise_share), rmse_obs_pp = median(CELLS$rmse_obs_pp), rmse_true_pp = median(CELLS$rmse_true_pp)))
write.csv(OVR, file.path(P2, "nz01_overall.csv"), row.names = FALSE)

# ---- CP-01 rerun vs the pre-fix snapshot ---------------------------------------
old <- read.csv(file.path(SNAP, "conformal_prev_cells.csv"), stringsAsFactors = FALSE)
new <- read.csv(file.path(P2, "conformal_prev_cells.csv"), stringsAsFactors = FALSE)
CMP <- full_join(old |> select(country, outcome, n_old = n, half_width_old = half_width_pp, mae_med_old = median_abs_err_pp, loo_old = loo_coverage_90),
                 new |> select(country, outcome, n_new = n, half_width_new = half_width_pp, mae_med_new = median_abs_err_pp, loo_new = loo_coverage_90),
                 by = c("country", "outcome")) |> mutate(half_width_change = half_width_new - half_width_old)
write.csv(CMP, file.path(P2, "nz01_cp01_rerun_compare.csv"), row.names = FALSE)

cat("\n===== CP-01 rerun: 90% band half-width (pp) and leave-one-out coverage, pre-fix -> post-fix =====\n")
print(as.data.frame(CMP |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
cat(sprintf("cells %d -> %d | median half-width %.1f -> %.1f pp | LOO coverage %.2f-%.2f -> %.2f-%.2f (mean %.3f -> %.3f)\n",
            nrow(old), nrow(new), median(old$half_width_pp), median(new$half_width_pp), min(old$loo_coverage_90), max(old$loo_coverage_90),
            min(new$loo_coverage_90), max(new$loo_coverage_90), mean(old$loo_coverage_90), mean(new$loo_coverage_90)))
cat("\n===== NZ-01: share of the observed squared error that is survey sampling noise =====\n")
print(as.data.frame(CELLS |> transmute(country, outcome, n, median_n_eff = round(median_n_eff, 1), mean_p = round(100 * mean_p, 1),
  bias_pp = round(mean_resid_pp, 1), rmse_obs_pp = round(rmse_obs_pp, 1), rmse_noise_pp = round(rmse_noise_pp, 1),
  noise_share = round(noise_share, 2), flag = ifelse(noise_exceeds_mse, "noise > MSE, floored", ""), rmse_true_pp = round(rmse_true_pp, 1),
  nat_used = round(100 * p_nat_anchor_used, 1), nat_postfix = round(100 * p_nat_postfix_nraw, 1))), row.names = FALSE)
print(OVR |> mutate(across(where(is.numeric), ~ round(.x, 4))), row.names = FALSE)
o <- OVR[1, ]
cat(sprintf("\nHEADLINE: About %.0f%% of the gap between the model and the survey's district figure is the survey's own sampling noise; against the true prevalence the typical error is about %.1f points rather than %.1f.\n",
            100 * o$noise_share, o$rmse_true_pp, o$rmse_obs_pp))
cat(sprintf("cells where the noise estimate exceeds the observed MSE: %d of %d\n", sum(CELLS$noise_exceeds_mse), nrow(CELLS)))
cat("\nDONE\n")
