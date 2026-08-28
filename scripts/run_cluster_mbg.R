# =============================================================================
# scripts/run_cluster_mbg.R
#
# WS5. Cluster-level spatial model, leave-one-Admin-2-out, over every
# country-outcome. See R/cluster_mbg.R for the design and for why the Earth
# Engine cluster-covariate extraction was dropped.
#
# Profiles (ANALYSIS_PROFILE):
#   smoke  child_iron only, first two countries
#   full   every country, every outcome
#
# Run:
#   Rscript scripts/run_cluster_mbg.R
#
# Writes results/tables/corrected/cluster_mbg_within_country.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(tidyr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

params    <- get_pipeline_params()
PROFILE   <- params$analysis_profile
SEED      <- params$seed
configs   <- get_country_configs()
countries <- names(configs)

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

if (PROFILE == "smoke") {
  countries <- countries[seq_len(min(2L, length(countries)))]
  only_outcomes <- "child_iron"
} else {
  only_outcomes <- NULL
}

have_spamm <- requireNamespace("spaMM", quietly = TRUE)
arms <- if (have_spamm) CLUSTER_MBG_ARMS else setdiff(CLUSTER_MBG_ARMS, "matern_spamm")

cat(sprintf("\n[ws5] profile=%s seed=%d spaMM=%s\n[ws5] arms: %s\n",
            PROFILE, SEED, have_spamm, paste(arms, collapse = ", ")))

rows <- list(); skipped <- list()

for (cn in countries) {
  cc <- configs[[cn]]
  merged <- read_target(paste0("merged_", tolower(cn)))
  acov   <- read_target(paste0("area_covariates_", tolower(cn)))
  if (is.null(merged) || is.null(acov) || is.null(acov$gee_admin2)) {
    cat(sprintf("[skip] %s: merged or area_covariates absent\n", cn)); next
  }
  area_cov <- acov$gee_admin2
  covs <- intersect(cluster_mbg_covariates(), names(area_cov))
  # area_covariates_*$gee_admin2 is kept in polygon order and NOT deduplicated,
  # so its Admin-2 names are not unique. Joining it on the name alone multiplied
  # Malawi's 103 clusters into 107. area_covariate_lookup() drops water polygons
  # and collapses on the pair key. See R/admin2_key_hygiene.R.
  ac <- area_covariate_lookup(area_cov, covs, sprintf("%s area covariates", cn))

  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    if (!is.null(only_outcomes) && !(oc$tag %in% only_outcomes)) next
    svy <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), oc$tag))
    if (is.null(svy)) next

    bio <- dist_continuous_biomarker(merged, cc, oc)
    if (is.null(bio)) {
      skipped[[length(skipped) + 1L]] <- data.frame(
        country = cc$country, outcome = oc$tag,
        reason = "no usable continuous biomarker", stringsAsFactors = FALSE)
      next
    }
    cl <- build_cluster_dataset(merged, cc, oc, bio$value, bio$threshold)
    if (is.null(cl)) {
      skipped[[length(skipped) + 1L]] <- data.frame(
        country = cc$country, outcome = oc$tag,
        reason = "too few clusters or no outcome variation", stringsAsFactors = FALSE)
      next
    }
    # Admin-2 covariates joined onto clusters. Constant within district by
    # construction; the arm rows record that.
    by_cov <- admin2_join_by(cl, ac)
    n_cl <- nrow(cl)
    cl <- dplyr::left_join(cl, ac, by = by_cov)
    cl <- warn_if_join_multiplied(n_cl, cl, sprintf("%s cluster-covariate join", cn))

    t0 <- Sys.time()
    r <- tryCatch(run_cluster_mbg(cl, svy, covs, cc, oc, arms = arms, seed = SEED),
                  error = function(e) { message("  ERROR ", cn, " ", oc$tag, ": ",
                                                conditionMessage(e)); NULL })
    if (is.null(r)) {
      skipped[[length(skipped) + 1L]] <- data.frame(
        country = cc$country, outcome = oc$tag,
        reason = "too few districts or no arm fitted", stringsAsFactors = FALSE)
      cat(sprintf("[skip] %s %s\n", cn, oc$tag)); next
    }
    r$elapsed_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    rows[[paste(cn, oc$tag)]] <- r
    cat(sprintf("[ok]   %-13s %-13s clusters=%3d districts=%3d (1-cluster: %3d)  %.1f min\n",
                cn, oc$tag, r$n_clusters[1], r$n_districts[1],
                r$districts_single_cluster[1], r$elapsed_min[1]))
  }
}

if (!length(rows)) stop("No cluster-MBG results produced.", call. = FALSE)
out <- dplyr::bind_rows(rows)
out$analysis_profile <- PROFILE

dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
suffix <- if (PROFILE == "smoke") "_SMOKE" else ""
readr::write_csv(out, here("results", "tables", "corrected",
                           sprintf("cluster_mbg_within_country%s.csv", suffix)))
if (length(skipped))
  readr::write_csv(dplyr::bind_rows(skipped),
                   here("results", "tables", "corrected",
                        sprintf("cluster_mbg_gaps%s.csv", suffix)))

cat("\n=== by arm, over all cells ===\n")
cat("national_mean correlation is mechanical: the leave-one-out training mean is\n",
    "anti-correlated with the held-out value by construction, so it is blanked in\n",
    "the correlation columns. MAE and RMSE are the metrics the null is read on.\n",
    sep = "")
print(as.data.frame(out |>
  dplyr::mutate(pearson_r = ifelse(correlation_is_mechanical, NA_real_, pearson_r),
                r_share   = ifelse(correlation_is_mechanical, NA_real_, r_share)) |>
  dplyr::group_by(arm) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   mean_pearson = round(mean(pearson_r, na.rm = TRUE), 4),
                   median_pearson = round(stats::median(pearson_r, na.rm = TRUE), 4),
                   mean_mae_pp = round(mean(mae_pp, na.rm = TRUE), 3),
                   mean_rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 3),
                   mean_abs_bias_pp = round(mean(abs(bias_pp), na.rm = TRUE), 3),
                   mean_r_share = round(mean(r_share, na.rm = TRUE), 3),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_pearson))), row.names = FALSE)

# Head-to-head against the national-mean null, which is what a district with no
# survey clusters would otherwise be given.
cat("\n=== each arm minus the national-mean null, paired within cell ===\n")
paired <- out |>
  dplyr::select(country, outcome, arm, pearson_r, mae_pp, rmse_pp) |>
  tidyr::pivot_longer(c(pearson_r, mae_pp, rmse_pp), names_to = "metric",
                      values_to = "value") |>
  tidyr::pivot_wider(names_from = arm, values_from = value)
if ("national_mean" %in% names(paired)) {
  cmp <- paired |>
    tidyr::pivot_longer(dplyr::any_of(setdiff(arms, "national_mean")),
                        names_to = "arm", values_to = "value") |>
    # Correlation against the null is dropped for the reason above. Only the
    # error metrics answer "is this better than giving an unsurveyed district
    # the national average".
    dplyr::filter(metric != "pearson_r") |>
    dplyr::filter(!is.na(value), !is.na(national_mean)) |>
    dplyr::mutate(delta = value - national_mean,
                  lower_better = metric %in% c("mae_pp", "rmse_pp"),
                  better = ifelse(lower_better, delta < 0, delta > 0))
  readr::write_csv(cmp, here("results", "tables", "corrected",
                             sprintf("cluster_mbg_vs_null%s.csv", suffix)))
  print(as.data.frame(cmp |> dplyr::group_by(metric, arm) |>
    dplyr::summarise(n = dplyr::n(), cells_better = sum(better, na.rm = TRUE),
                     mean_delta = round(mean(delta), 4),
                     .groups = "drop") |>
    dplyr::arrange(metric, dplyr::desc(cells_better))), row.names = FALSE)
}

if (length(skipped)) {
  cat("\n=== skipped cells ===\n")
  print(as.data.frame(dplyr::bind_rows(skipped)), row.names = FALSE)
}
cat("\nwrote results/tables/corrected/cluster_mbg_within_country", suffix, ".csv\n", sep = "")
