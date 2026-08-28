# =============================================================================
# scripts/run_distributional.R
#
# WS4. Run the distributional estimator against the binary-classifier path for
# every country-outcome that has a continuous biomarker.
# See R/corrected/p12_distributional.R.
#
# Profiles (ANALYSIS_PROFILE):
#   smoke  child_iron only, first two countries
#   full   every country, every outcome with a continuous biomarker
#
# Run:
#   Rscript scripts/run_distributional.R
#
# Writes results/tables/corrected/distributional_comparison.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
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

cat(sprintf("\n[ws4] profile=%s seed=%d countries=%s\n",
            PROFILE, SEED, paste(countries, collapse = ", ")))

rows <- list()
skipped <- list()

for (cn in countries) {
  cc <- configs[[cn]]
  merged <- read_target(paste0("merged_", tolower(cn)))
  acov   <- read_target(paste0("area_covariates_", tolower(cn)))
  if (is.null(merged) || is.null(acov)) {
    cat(sprintf("[skip] %s: merged or area_covariates absent from the store\n", cn)); next
  }
  area_cov <- acov$gee_admin2
  if (is.null(area_cov)) { cat(sprintf("[skip] %s: no gee_admin2 in area_covariates\n", cn)); next }

  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    if (!is.null(only_outcomes) && !(oc$tag %in% only_outcomes)) next
    svy <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), oc$tag))
    if (is.null(svy)) {
      skipped[[length(skipped) + 1L]] <- data.frame(
        country = cc$country, outcome = oc$tag,
        reason = "no svy_admin2_* target in the store", stringsAsFactors = FALSE)
      next
    }
    t0 <- Sys.time()
    r <- tryCatch(run_distributional_cell(merged, area_cov, svy, cc, oc, seed = SEED),
                  error = function(e) { message("  ERROR ", cn, " ", oc$tag, ": ",
                                                conditionMessage(e)); NULL })
    if (is.null(r)) {
      # Recorded, not imputed. The brief is explicit that a missing continuous
      # biomarker is documented rather than filled in.
      skipped[[length(skipped) + 1L]] <- data.frame(
        country = cc$country, outcome = oc$tag,
        reason = "no usable continuous biomarker, too few areas, or no outcome variation",
        stringsAsFactors = FALSE)
      cat(sprintf("[skip] %s %s\n", cn, oc$tag))
      next
    }
    r$elapsed_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    rows[[paste(cn, oc$tag)]] <- r
    cat(sprintf("[ok]   %-14s %-13s prev=%.1f%%  areas=%d  (%.1f min)\n",
                cn, oc$tag, 100 * r$national_prev[1], r$n_areas_total[1],
                r$elapsed_min[1]))
  }
}

if (!length(rows)) stop("No distributional results produced.", call. = FALSE)
out <- dplyr::bind_rows(rows)
out$analysis_profile <- PROFILE

dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
suffix <- if (PROFILE == "smoke") "_SMOKE" else ""
readr::write_csv(out, here("results", "tables", "corrected",
                           sprintf("distributional_comparison%s.csv", suffix)))
if (length(skipped))
  readr::write_csv(dplyr::bind_rows(skipped),
                   here("results", "tables", "corrected",
                        sprintf("distributional_gaps%s.csv", suffix)))

cat("\n=== per-cell comparison ===\n")
print(as.data.frame(out |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))) |>
  dplyr::select(country, outcome, scheme, national_prev, auc, brier_skill,
                pearson_r, r_max, r_share, calib_slope)), row.names = FALSE)

cat("\n=== paired, distributional minus binary ===\n")
paired <- out |>
  dplyr::select(country, outcome, scheme, auc, brier_skill, pearson_r, mae_pp,
                calib_slope, r_share) |>
  tidyr::pivot_longer(c(auc, brier_skill, pearson_r, mae_pp, calib_slope, r_share),
                      names_to = "metric", values_to = "value") |>
  tidyr::pivot_wider(names_from = scheme, values_from = value) |>
  dplyr::filter(!is.na(binary_classifier), !is.na(distributional_default)) |>
  dplyr::mutate(delta = distributional_default - binary_classifier,
                lower_better = metric %in% c("mae_pp"),
                dist_better = ifelse(lower_better, delta < 0, delta > 0))
readr::write_csv(paired, here("results", "tables", "corrected",
                              sprintf("distributional_paired%s.csv", suffix)))

# A cell whose reliability ceiling is indistinguishable from zero has no
# predictable between-district variation at all (R/area_reliability.R:98), so a
# correlation computed on it is noise and averaging it into a recommendation
# would be misleading. The `signal` flag uses the OPTIMISTIC upper bound, so a
# cell is excluded only when even that says there is nothing there.
has_signal <- out |>
  dplyr::group_by(country, outcome) |>
  dplyr::summarise(has_signal = any(signal %in% TRUE), .groups = "drop")
paired_sig <- paired |>
  dplyr::inner_join(has_signal, by = c("country", "outcome")) |>
  dplyr::filter(has_signal)
readr::write_csv(paired_sig, here("results", "tables", "corrected",
                                  sprintf("distributional_paired_signal%s.csv", suffix)))

cat("\n=== restricted to cells with signal ===\n")
print(as.data.frame(paired_sig |> dplyr::group_by(metric) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   cells_dist_better = sum(dist_better, na.rm = TRUE),
                   mean_binary = round(mean(binary_classifier), 4),
                   mean_dist   = round(mean(distributional_default), 4),
                   mean_delta  = round(mean(delta), 4),
                   min_delta   = round(min(delta), 4),
                   max_delta   = round(max(delta), 4),
                   lower_better = dplyr::first(lower_better),
                   .groups = "drop")), row.names = FALSE)

print(as.data.frame(paired |> dplyr::group_by(metric) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   cells_dist_better = sum(dist_better, na.rm = TRUE),
                   mean_binary = round(mean(binary_classifier), 4),
                   mean_dist   = round(mean(distributional_default), 4),
                   mean_delta  = round(mean(delta), 4),
                   min_delta   = round(min(delta), 4),
                   max_delta   = round(max(delta), 4),
                   lower_better = dplyr::first(lower_better),
                   .groups = "drop")), row.names = FALSE)

if (length(skipped)) {
  cat("\n=== cells with no usable continuous biomarker (documented, not imputed) ===\n")
  print(as.data.frame(dplyr::bind_rows(skipped)), row.names = FALSE)
}
cat("\nwrote results/tables/corrected/distributional_comparison", suffix, ".csv\n", sep = "")
