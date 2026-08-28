# =============================================================================
# scripts/run_subsample.R
#
# WS6. Rebuild the survey-subsample cost-of-accuracy simulation from its prose
# specification, add the spatial variant, and trace the calibration learning
# curve. See R/corrected/p14_subsample.R.
#
# Profiles (ANALYSIS_PROFILE):
#   smoke  child_iron only, first two countries, 5 replicates
#   full   every country, the four primary outcomes, 25 replicates
#
# Run:
#   Rscript scripts/run_subsample.R
#
# Writes
#   results/tables/corrected/subsample_cost_of_accuracy.csv
#   results/tables/corrected/subsample_summary.csv
#   results/tables/corrected/calibration_learning_curve.csv
#   results/figures/calibration_learning_curve.png
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
OUTCOMES  <- c("child_vitA", "women_vitA", "child_iron", "women_iron")
REPS      <- SUBSAMPLE_REPLICATES

if (PROFILE == "smoke") {
  countries <- countries[seq_len(min(2L, length(countries)))]
  OUTCOMES  <- "child_iron"
  REPS      <- 5L
}

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

cat(sprintf("\n[ws6] profile=%s seed=%d reps=%d\n[ws6] retention: %s\n",
            PROFILE, SEED, REPS, paste(SUBSAMPLE_RETENTION, collapse = ", ")))

# ── 6a and 6b: retention simulation ─────────────────────────────────────────
rows <- list()
cluster_cache <- list()

for (cn in countries) {
  cc <- configs[[cn]]
  merged <- read_target(paste0("merged_", tolower(cn)))
  acov   <- read_target(paste0("area_covariates_", tolower(cn)))
  if (is.null(merged) || is.null(acov) || is.null(acov$gee_admin2)) next
  area_cov <- acov$gee_admin2
  covs <- intersect(cluster_mbg_covariates(), names(area_cov))
  ac <- area_cov[, c("Admin2", covs), drop = FALSE]
  ac$Admin2 <- as.character(ac$Admin2)
  for (k in covs) ac[[k]] <- suppressWarnings(as.numeric(ac[[k]]))

  for (otag in OUTCOMES) {
    oc <- cc$outcomes[[otag]]
    if (is.null(oc)) next
    svy <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), otag))
    if (is.null(svy)) next
    bio <- dist_continuous_biomarker(merged, cc, oc)
    if (is.null(bio)) next
    cl <- build_cluster_dataset(merged, cc, oc, bio$value, bio$threshold)
    if (is.null(cl)) next
    cl <- dplyr::left_join(cl, ac, by = "Admin2")
    cluster_cache[[paste(cn, otag)]] <- list(cl = cl, svy = svy, cc = cc, oc = oc)

    t0 <- Sys.time()
    r <- tryCatch(run_subsample_simulation(cl, svy, covs, reps = REPS, seed = SEED),
                  error = function(e) { message("  ERROR ", cn, " ", otag, ": ",
                                                conditionMessage(e)); NULL })
    if (is.null(r)) { cat(sprintf("[skip] %s %s\n", cn, otag)); next }
    r$country <- cc$country; r$outcome <- otag
    rows[[paste(cn, otag)]] <- r
    cat(sprintf("[ok]   %-13s %-12s clusters=%3d districts=%3d  %.1f min\n",
                cn, otag, r$n_clusters_total[1], r$n_districts_total[1],
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

if (!length(rows)) stop("No subsample results.", call. = FALSE)
sub <- dplyr::bind_rows(rows)
sub$analysis_profile <- PROFILE
dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
suffix <- if (PROFILE == "smoke") "_SMOKE" else ""
readr::write_csv(sub, here("results", "tables", "corrected",
                           sprintf("subsample_cost_of_accuracy%s.csv", suffix)))

summ <- sub |>
  dplyr::group_by(stratum, strategy) |>
  dplyr::summarise(n_rows = dplyr::n(),
                   mean_rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 3),
                   mean_mae_pp = round(mean(mae_pp, na.rm = TRUE), 3),
                   mean_districts = round(mean(n_districts), 1),
                   .groups = "drop") |>
  dplyr::arrange(stratum, mean_rmse_pp)
readr::write_csv(summ, here("results", "tables", "corrected",
                            sprintf("subsample_summary%s.csv", suffix)))

cat("\n=== RMSE against the full survey, by coverage stratum ===\n")
print(as.data.frame(summ), row.names = FALSE)

# The two headline contrasts the August note states, recomputed here.
.delta <- function(stratum, a, b) {
  d <- sub |> dplyr::filter(stratum == !!stratum, strategy %in% c(a, b)) |>
    dplyr::select(country, outcome, retention, replicate, strategy, rmse_pp) |>
    tidyr::pivot_wider(names_from = strategy, values_from = rmse_pp) |>
    dplyr::filter(is.finite(.data[[a]]), is.finite(.data[[b]]))
  if (!nrow(d)) return(NULL)
  data.frame(stratum = stratum, comparison = paste(a, "minus", b),
             n_replicates = nrow(d),
             mean_delta_pp = round(mean(d[[a]] - d[[b]]), 4),
             median_delta_pp = round(stats::median(d[[a]] - d[[b]]), 4),
             pct_replicates_a_better = round(100 * mean(d[[a]] < d[[b]]), 1),
             stringsAsFactors = FALSE)
}
contrasts <- dplyr::bind_rows(
  .delta("has_clusters",  "model_cv",   "direct_or_national"),
  .delta("has_clusters",  "spatial_cv", "direct_or_national"),
  .delta("zero_clusters", "model_cv",   "direct_or_national"),
  .delta("zero_clusters", "spatial_cv", "direct_or_national"),
  .delta("zero_clusters", "spatial_cv", "model_cv")
)
readr::write_csv(contrasts, here("results", "tables", "corrected",
                                 sprintf("subsample_contrasts%s.csv", suffix)))
cat("\n=== headline contrasts (negative favours the first strategy) ===\n")
print(as.data.frame(contrasts), row.names = FALSE)

# ── 6c: calibration learning curve ──────────────────────────────────────────
# Transported deployment: for each target country, the area model is trained on
# the OTHER countries, then calibrated on k of the target's own districts.
lc_rows <- list()
for (otag in OUTCOMES) {
  svy <- list()
  for (cn in countries) {
    s <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), otag))
    if (!is.null(s) && nrow(s) > 0) svy[[cn]] <- s
  }
  if (length(svy) < 2) next
  gee <- list()
  for (cn in names(svy)) {
    g <- read_target(paste0("gee_admin2_", tolower(cn)))
    if (!is.null(g) && ncol(g) > 2) gee[[cn]] <- g
  }
  p <- build_area_loco_dataset(svy, gee[names(svy)])
  d <- p$pooled_data
  if (is.null(d) || !nrow(d)) next
  if (!"n_svy" %in% names(d)) d$n_svy <- 1
  pool <- list(data = d, predictors = p$common_gee_vars,
               countries = p$country_names, outcome = otag)
  res <- tryCatch(run_area_transport_loco(pool, AREA_TRANSPORT_RECIPE),
                  error = function(e) NULL)
  if (is.null(res) || is.null(res$predictions)) next
  for (ho in unique(res$predictions$country)) {
    dd <- res$predictions[res$predictions$country == ho, , drop = FALSE]
    lc <- calibration_learning_curve(dd$modeled_prev, dd$survey_prev,
                                     reps = REPS, seed = SEED)
    if (is.null(lc)) next
    lc$outcome <- otag; lc$country <- ho
    lc_rows[[paste(otag, ho)]] <- lc
  }
}

if (length(lc_rows)) {
  lc <- dplyr::bind_rows(lc_rows)
  lc$analysis_profile <- PROFILE
  readr::write_csv(lc, here("results", "tables", "corrected",
                            sprintf("calibration_learning_curve%s.csv", suffix)))
  cat("\n=== calibration learning curve: held-out RMSE by sentinel count ===\n")
  lcs <- lc |> dplyr::group_by(k_sentinel) |>
    dplyr::summarise(n = dplyr::n(),
                     mean_rmse_pp = round(mean(rmse_pp), 3),
                     lo = round(stats::quantile(rmse_pp, .025), 3),
                     hi = round(stats::quantile(rmse_pp, .975), 3),
                     .groups = "drop")
  print(as.data.frame(lcs), row.names = FALSE)

  dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
  png(here("results", "figures", sprintf("calibration_learning_curve%s.png", suffix)),
      width = 1000, height = 700, res = 130)
  plot(lcs$k_sentinel, lcs$mean_rmse_pp, type = "b", pch = 19,
       ylim = range(c(lcs$lo, lcs$hi), na.rm = TRUE),
       xlab = "sentinel districts used for calibration (k)",
       ylab = "held-out district RMSE (pp)",
       main = "Calibration learning curve, transported area model")
  graphics::arrows(lcs$k_sentinel, lcs$lo, lcs$k_sentinel, lcs$hi,
                   angle = 90, code = 3, length = 0.04, col = "grey40")
  graphics::grid()
  invisible(dev.off())
  cat("wrote results/figures/calibration_learning_curve", suffix, ".png\n", sep = "")
}

cat("\nwrote results/tables/corrected/subsample_cost_of_accuracy", suffix, ".csv\n", sep = "")
