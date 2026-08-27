# =============================================================================
# sandbox_parsimony/R/99_summary.R -- print every headline table in one place
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
O <- function(f) file.path("sandbox_parsimony/out", f)
rd <- function(f) if (file.exists(O(f))) read.csv(O(f), stringsAsFactors = FALSE) else NULL
hdr <- function(x) cat("\n\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- 1. how much signal is there ------------------------------------------
a <- rd("noise_audit.csv")
if (!is.null(a)) {
  hdr("1. Admin-2 outcome reliability (deff = 1.5, 95% bootstrap band)")
  print(a |> mutate(cell = paste(outcome, country),
                    lambda = sprintf("%.2f [%.2f, %.2f]", lambda_d15,
                                     lambda_d15_lo, lambda_d15_hi)) |>
          select(cell, n_areas, median_n, pct_n_lt_25, mean_prev_pp, sd_obs_pp,
                 sd_samp_pp_d15, lambda, r_max_d15, r_max_d15_hi) |>
          arrange(desc(r_max_d15)), row.names = FALSE)
}
c2 <- rd("ceiling_vs_published.csv")
if (!is.null(c2)) {
  hdr("2. Published within-country CV r as a share of the achievable ceiling")
  print(c2 |> select(outcome, country, n, cv_pearson, r_max_d15, r_max_d15_hi,
                     pct_of_ceiling, pct_of_ceiling_optimistic), row.names = FALSE)
}
fh <- rd("fh_shrinkage_check.csv")
if (!is.null(fh)) {
  hdr("3. Fay-Herriot shrinkage weight: as coded vs. with a real sampling variance")
  print(fh |> summarise(
    cells = n(),
    mean_gamma_as_coded = round(mean(mean_gamma_as_coded), 3),
    mean_gamma_correct  = round(mean(mean_gamma_correct), 3),
    pct_no_shrinkage_as_coded = round(mean(pct_gamma_gt_95_as_coded), 0),
    pct_no_shrinkage_correct  = round(mean(pct_gamma_gt_95_correct), 0)),
    row.names = FALSE)
  cat("\ngamma near 1 = the EBLUP just returns the raw direct estimate.\n")
}

# --- 4. within-country bake-off -------------------------------------------
w <- bind_rows(rd("within_country_bakeoff.csv"), rd("combine_space.csv"))
if (!nrow(w)) w <- NULL
# 08_spatial_models.R deliberately skips cells with no detectable signal, so the
# binomial-GAM rows cover a SUBSET of the cells. Keeping them in the main mean
# would compare models on different data, so they get their own paired table.
sp <- rd("spatial_bakeoff.csv")
if (!is.null(w)) {
  # The bake-off CSVs carry whatever r_max was current when they were written.
  # Always re-join the ceiling from the (single, authoritative) noise audit so
  # r_share and the "signal-bearing" filter cannot drift between runs.
  if (!is.null(a)) {
    w <- w |> select(-any_of(c("r_max", "reliability", "r_share"))) |>
      left_join(a |> select(outcome, country, r_max = r_max_d15,
                            reliability = lambda_d15),
                by = c("outcome", "country")) |>
      mutate(r_share = ifelse(r_max > .05, round(pearson / r_max, 2), NA_real_))
  }
  hdr("4. Within-country bake-off (repeated 5-fold CV, mean over 16 cells)")
  print(w |> group_by(model) |>
          summarise(cells = n(),
                    rho = round(mean(spearman, na.rm = TRUE), 3),
                    r = round(mean(pearson, na.rm = TRUE), 3),
                    r_share = round(mean(r_share, na.rm = TRUE), 2),
                    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
                    med_fold_sd = round(median(pearson_sd, na.rm = TRUE), 3),
                    .groups = "drop") |> arrange(desc(rho)) |> as.data.frame(),
        row.names = FALSE)
  hdr("4b. Same, restricted to cells with detectable signal (r_max >= 0.35)")
  print(w |> filter(r_max >= 0.35) |> group_by(model) |>
          summarise(cells = n(),
                    rho = round(mean(spearman, na.rm = TRUE), 3),
                    r = round(mean(pearson, na.rm = TRUE), 3),
                    r_share = round(mean(r_share, na.rm = TRUE), 2),
                    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
                    .groups = "drop") |> arrange(desc(rho)) |> as.data.frame(),
        row.names = FALSE)

  if (!is.null(sp)) {
    cells <- unique(paste(sp$outcome, sp$country))
    both <- bind_rows(w, sp) |> filter(paste(outcome, country) %in% cells)
    hdr(sprintf("4c. Binomial GAM (counts + spatial field) vs the rest, on its %d cells",
                length(cells)))
    print(both |> group_by(model) |>
            summarise(cells = n(),
                      rho = round(mean(spearman, na.rm = TRUE), 3),
                      r = round(mean(pearson, na.rm = TRUE), 3),
                      rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2),
                      .groups = "drop") |> arrange(desc(rho)) |> as.data.frame(),
          row.names = FALSE)
  }
}

# --- 5. LOCO --------------------------------------------------------------
l <- rd("loco_headline.csv")
if (!is.null(l)) {
  hdr("5. LOCO bake-off (4 outcomes x 4 held-out countries, as in benchmarks_all.csv)")
  print(l |> group_by(variant, anchor) |>
          summarise(cells = n(),
                    rho = round(mean(spearman, na.rm = TRUE), 3),
                    rho_sd = round(sd(spearman, na.rm = TRUE), 3),
                    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 1),
                    pattern_rmse_pp = round(mean(pattern_rmse_pp, na.rm = TRUE), 1),
                    abs_level_bias_pp = round(mean(abs(level_bias_pp), na.rm = TRUE), 1),
                    .groups = "drop") |>
          arrange(anchor, desc(rho)) |> as.data.frame(), row.names = FALSE)
  hdr("5b. What a known national anchor buys, per cell (best model each way)")
  print(l |> group_by(outcome, held_out, anchor) |>
          summarise(best_rmse = round(min(rmse_pp, na.rm = TRUE), 1), .groups = "drop") |>
          tidyr::pivot_wider(names_from = anchor, values_from = best_rmse) |>
          mutate(rmse_pp_saved = round(train_mean - oracle_national, 1)) |>
          arrange(desc(rmse_pp_saved)) |> as.data.frame(), row.names = FALSE)
}

sc <- rd("scale_sweep.csv")
if (!is.null(sc)) {
  hdr("6. Response scale: logit (production) vs raw prevalence, LOCO")
  print(sc |> group_by(model, scale) |>
          summarise(rho = round(mean(spearman, na.rm = TRUE), 3), .groups = "drop") |>
          tidyr::pivot_wider(names_from = scale, values_from = rho) |>
          mutate(gain_from_raw_scale = round(continuous - logit, 3)) |>
          arrange(desc(gain_from_raw_scale)) |> as.data.frame(), row.names = FALSE)
}

pc <- rd("pool_composition.csv")
if (!is.null(pc)) {
  hdr("7. Cost of pooling with a country whose GEE extraction is incomplete")
  print(pc |> group_by(outcome, pool, n_pred_pool) |>
          summarise(cells = n(), rho = round(mean(spearman, na.rm = TRUE), 3),
                    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 1), .groups = "drop") |>
          as.data.frame(), row.names = FALSE)
}

hy <- rd("hygiene_effect.csv")
if (!is.null(hy)) {
  hdr("8. GEE_COVARIATE_HYGIENE on vs off (paired on the same cells)")
  print(hy |> group_by(eval, model, hygiene) |>
          summarise(cells = n(), rho = round(mean(spearman, na.rm = TRUE), 3),
                    rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 2), .groups = "drop") |>
          tidyr::pivot_wider(names_from = hygiene, values_from = c(rho, rmse_pp)) |>
          mutate(rho_gain = round(rho_on - rho_off, 3)) |> as.data.frame(),
        row.names = FALSE)
}

jd <- rd("join_diagnostics.csv")
if (!is.null(jd)) {
  hdr("9. Join sanity (Moran I on covariates should be well above 0)")
  print(jd |> group_by(country) |>
          summarise(cells = n(), moran_cov = round(mean(moran_cov), 3),
                    moran_outcome = round(mean(moran_out), 3), .groups = "drop") |>
          as.data.frame(), row.names = FALSE)
}
cat("\n")
