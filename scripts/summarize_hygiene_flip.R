# =============================================================================
# scripts/summarize_hygiene_flip.R
#
# WS2b. Turn the v1/v2 covariate-hygiene comparison into the paired, fold-level
# evidence the default-flip decision is made on.
#
# scripts/compare_covariate_hygiene.R reports means across held-out countries.
# A mean over 4 folds hides whether the two variants move together, so this
# script pairs every (outcome, held_out) fold and reports the within-fold delta.
#
# With 4 held-out countries per outcome, inference is limited to the paired
# deltas themselves and their range. No p-value is reported: 4 folds do not
# support one.
#
# Run:
#   Rscript scripts/summarize_hygiene_flip.R
#
# Writes results/sensitivity/covariate_hygiene_paired_<scheme>.csv (one row per
# outcome x held-out country x metric) and covariate_hygiene_paired_summary_<scheme>.csv
# (one row per metric), so every figure quoted in the findings file is read from
# a produced CSV rather than from console output.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(tidyr)
})

SCHEME <- Sys.getenv("HYGIENE_SCHEME", "ws2b_2026-08")
IN_PATH <- here("results", "sensitivity",
                sprintf("covariate_hygiene_comparison_%s.csv", SCHEME))
OUT_PATH <- here("results", "sensitivity",
                 sprintf("covariate_hygiene_paired_%s.csv", SCHEME))

if (!file.exists(IN_PATH))
  stop("Run scripts/compare_covariate_hygiene.R first; missing ", IN_PATH)

cmp <- readr::read_csv(IN_PATH, show_col_types = FALSE)

# Metrics where a LOWER value is better, so the sign of "improvement" is
# recorded per metric rather than assumed.
LOWER_BETTER <- c("rmse_pp", "mae_pp", "abs_nat_bias_pp")
METRICS <- c("pearson_r", "spearman_r", "rmse_pp", "mae_pp", "abs_nat_bias_pp",
             "n_selected")

cmp$abs_nat_bias_pp <- abs(cmp$nat_bias_pp)
METRICS <- intersect(METRICS, names(cmp))

paired <- cmp |>
  dplyr::select(outcome, held_out, variant, dplyr::all_of(METRICS)) |>
  tidyr::pivot_longer(dplyr::all_of(METRICS), names_to = "metric",
                      values_to = "value") |>
  tidyr::pivot_wider(names_from = variant, values_from = value) |>
  dplyr::filter(!is.na(v1_current), !is.na(v2_hygiene)) |>
  dplyr::mutate(
    delta       = v2_hygiene - v1_current,
    lower_better = metric %in% LOWER_BETTER,
    v2_better   = ifelse(lower_better, delta < 0, delta > 0),
    scheme      = SCHEME
  )

# Fold-level summary: how many of the 16 folds each metric improves in, and the
# spread of the paired delta. `folds_better` is the headline; the mean delta is
# reported alongside it because a mean can be carried by one fold.
summ <- paired |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    n_folds      = dplyr::n(),
    folds_v2_better = sum(v2_better, na.rm = TRUE),
    mean_delta   = round(mean(delta, na.rm = TRUE), 4),
    median_delta = round(stats::median(delta, na.rm = TRUE), 4),
    min_delta    = round(min(delta, na.rm = TRUE), 4),
    max_delta    = round(max(delta, na.rm = TRUE), 4),
    lower_better = dplyr::first(lower_better),
    .groups = "drop"
  ) |>
  dplyr::arrange(metric)

readr::write_csv(paired, OUT_PATH)
readr::write_csv(summ, sub("_paired_", "_paired_summary_", OUT_PATH))

cat("\n=== paired fold-level deltas (v2_hygiene - v1_current) ===\n")
cat(sprintf("scheme: %s   folds: %d outcomes x %d held-out countries\n",
            SCHEME, dplyr::n_distinct(paired$outcome),
            dplyr::n_distinct(paired$held_out)))
print(as.data.frame(summ), row.names = FALSE)
cat("\nwrote ", OUT_PATH, "\n", sep = "")
