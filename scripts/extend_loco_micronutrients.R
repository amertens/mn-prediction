# =============================================================================
# scripts/extend_loco_micronutrients.R
#
# Extend the individual-level LOCO transportability table to folate and B12.
# Uses the SAME inputs as the pipeline's loco_* targets (merged_<country>,
# sl_learners, pipeline_params) so the new rows are comparable to the existing
# four outcomes. build_pooled_dataset() auto-skips countries lacking the
# outcome, so folate/B12 pool Ghana + Sierra Leone + Malawi (Gambia has neither).
#
# NOTE: zinc is Malawi-only — LOCO needs >=2 countries to hold one out, so zinc
# cannot have transport folds and is intentionally excluded.
#
# Appends to results/tables/transportability_loco_results.csv (idempotent:
# drops any pre-existing rows for these outcomes before appending).
# =============================================================================

suppressWarnings(suppressMessages({ library(targets); targets::tar_source("R/") }))
suppressPackageStartupMessages({ library(dplyr) })

STORE <- "_targets_full"
# Outcomes to (re)compute; override on the command line, e.g.
#   Rscript scripts/extend_loco_micronutrients.R women_folate
NEW_OUTCOMES <- c("women_folate", "women_b12")
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) > 0) NEW_OUTCOMES <- .args
CSV <- "results/tables/transportability_loco_results.csv"

cfgs <- get_country_configs()
all_merged <- list(
  Gambia      = tar_read_raw("merged_gambia",      store = STORE),
  Ghana       = tar_read_raw("merged_ghana",       store = STORE),
  SierraLeone = tar_read_raw("merged_sierraleone", store = STORE),
  Malawi      = tar_read_raw("merged_malawi",      store = STORE)
)
sl_learners <- tar_read_raw("sl_learners",     store = STORE)
params      <- tar_read_raw("pipeline_params", store = STORE)

# Harmonize Malawi folate to the cross-country definition (serum folate
# < 10 nmol/L), matching Ghana/Sierra Leone (gw_wFolate < 10). Malawi's native
# folate_def uses a stricter < 3 ng/mL (~6.8 nmol/L) rule that yields ~0.1%
# deficiency — not comparable. fol_nmol is serum folate in nmol/L
# (= fol[ng/mL] x 2.266), so the same < 10 cutoff applies directly.
if ("fol_nmol" %in% names(all_merged$Malawi)) {
  fn <- suppressWarnings(as.numeric(all_merged$Malawi$fol_nmol))
  all_merged$Malawi$folate_def <- ifelse(is.na(fn), NA_integer_, as.integer(fn < 10))
  cat(sprintf("[harmonize] Malawi folate redefined as fol_nmol < 10 nmol/L: %.1f%% deficient (n=%d)\n",
              100 * mean(all_merged$Malawi$folate_def, na.rm = TRUE),
              sum(!is.na(all_merged$Malawi$folate_def))))
}

new_rows <- list()
for (otag in NEW_OUTCOMES) {
  cat(sprintf("\n================ LOCO: %s ================\n", otag))
  pooled <- tryCatch(build_pooled_dataset(all_merged, cfgs, otag),
                     error = function(e) { cat("  pool err:", conditionMessage(e), "\n"); NULL })
  if (is.null(pooled)) next
  res <- tryCatch(run_loco_cv(pooled, sl_learners, params),
                  error = function(e) { cat("  loco err:", conditionMessage(e), "\n"); NULL })
  if (is.null(res) || nrow(res) == 0) next
  res$outcome <- otag
  new_rows[[otag]] <- res
  # Write incrementally so partial progress survives a crash.
  saveRDS(bind_rows(new_rows), "results/tables/_loco_new_partial.rds")
  cat(sprintf("  -> %d folds for %s\n", nrow(res), otag))
}

if (length(new_rows) == 0) { cat("\nNo new LOCO rows produced.\n"); quit(status = 1) }
new_df <- bind_rows(new_rows)

old <- read.csv(CSV, stringsAsFactors = FALSE)
keep_cols <- names(old)
# Align new rows to the existing schema (fill any missing cols with NA).
for (cc in setdiff(keep_cols, names(new_df))) new_df[[cc]] <- NA
new_df <- new_df[, keep_cols, drop = FALSE]
# Idempotent: drop prior rows for these outcomes, then append.
old <- old[!old$outcome %in% NEW_OUTCOMES, , drop = FALSE]
combined <- rbind(old, new_df)
write.csv(combined, CSV, row.names = FALSE)

cat(sprintf("\nAppended %d rows (%s). Table now %d rows.\n",
            nrow(new_df), paste(NEW_OUTCOMES, collapse = ", "), nrow(combined)))
print(new_df[, c("held_out", "n_test", "prev_test", "pred_prev", "auc", "calib_slope", "outcome")])

# Hard exit so any package/teardown hook from tar_source can't mask success.
quit(save = "no", status = 0)
