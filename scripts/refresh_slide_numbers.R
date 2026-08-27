# =============================================================================
# scripts/refresh_slide_numbers.R
#
# Print the figures the slide deck's Summary and Results slides quote, read
# from the CURRENT targets store. Run this AFTER a completed pipeline rebuild,
# then paste the numbers into docs/mn_prediction_slides.qmd.
#
# WHY THIS EXISTS. The Summary slide carries hard-coded claims -- "Ghana child
# iron AUC = 0.828", "LOCO transportability is moderate (AUC 0.47-0.73)" -- that
# were written before the 2026-08 methods review and before the fixes that
# change model inputs (water-body snapping, Tanzania's restored covariates,
# duplicate-key collapsing, benchmarking). They cannot be updated by hand
# without re-deriving them, and re-deriving them by hand is how stale numbers
# get perpetuated.
#
# It prints nothing it cannot read. A missing target is reported as missing
# rather than filled with the old value.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

STORE <- "_targets_full/objects"
rd <- function(nm) {
  p <- file.path(STORE, nm)
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
}
hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

cat(sprintf("store: %s   (%d objects, newest %s)\n", STORE,
            length(list.files(STORE)),
            format(max(file.info(list.files(STORE, full.names = TRUE))$mtime),
                   "%Y-%m-%d %H:%M")))

# ---- best within-country binary performance --------------------------------
hdr("1. Best within-country AUC  (Summary slide quotes Ghana child iron 0.828)")
cvp <- rd("cv_perf")
if (is.null(cvp)) {
  cat("cv_perf not in the store - rebuild has not reached it yet\n")
} else {
  auc_col <- intersect(c("auc", "AUC", "roc_auc"), names(cvp))
  if (!length(auc_col)) {
    cat("no AUC column in cv_perf; columns are:\n"); print(names(cvp))
  } else {
    top <- cvp |> filter(is.finite(.data[[auc_col[1]]])) |>
      arrange(desc(.data[[auc_col[1]]])) |>
      select(any_of(c("country", "outcome", "model_type")), all_of(auc_col[1])) |>
      head(8)
    print(as.data.frame(top), row.names = FALSE)
  }
}

# ---- LOCO transportability --------------------------------------------------
hdr("2. LOCO transportability  (Summary slide quotes AUC 0.47-0.73)")
lo <- rd("transportability_area_loco_metrics")
if (is.null(lo)) {
  f <- "results/tables/transportability_area_loco_metrics.csv"
  lo <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}
if (is.null(lo)) cat("LOCO metrics unavailable\n") else {
  print(as.data.frame(lo |> group_by(outcome) |>
    summarise(cells = n(),
              spearman = round(mean(spearman_r, na.rm = TRUE), 3),
              pearson = round(mean(pearson_r, na.rm = TRUE), 3),
              rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 1),
              nat_bias_pp = round(mean(nat_bias_pp, na.rm = TRUE), 1),
              .groups = "drop")), row.names = FALSE)
  cat("\nNOTE: quote Spearman (spatial pattern) and national bias (level)\n")
  cat("SEPARATELY. A single AUC or r hides which of the two failed -- see\n")
  cat("FINDINGS.md sections 6 and 19.\n")
}

# ---- area-level comparison, with the new benchmark columns ------------------
hdr("3. Area-level comparison  (now carries r_max / r_share / signal / national_err)")
ac <- rd("area_comparison_all")
if (is.null(ac)) {
  f <- "results/tables/area_comparison_all.csv"
  ac <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}
if (is.null(ac)) cat("area_comparison_all unavailable\n") else {
  keep <- intersect(c("country", "outcome", "approach", "mae_pp", "pearson_r",
                      "r_max", "r_share", "signal", "national_err_pp"), names(ac))
  cat(sprintf("rows: %d", nrow(ac)))
  if ("signal" %in% names(ac))
    cat(sprintf("   |   cells with NO detectable signal: %d of %d\n",
                sum(!ac$signal, na.rm = TRUE), nrow(ac))) else cat("\n")
  print(as.data.frame(head(ac[, keep, drop = FALSE], 12)), row.names = FALSE)
}

# ---- how many districts / countries the deck should claim -------------------
hdr("4. Coverage")
for (c in c("gambia", "ghana", "sierraleone", "malawi", "tanzania")) {
  g <- rd(paste0("gee_admin2_", c))
  s <- rd(paste0("svy_admin2_", c, "_child_vitA"))
  cat(sprintf("  %-12s districts=%-5s surveyed=%-5s gee_cols=%s\n", c,
              if (is.null(g)) "?" else nrow(g),
              if (is.null(s)) "?" else nrow(s),
              if (is.null(g)) "?" else sum(grepl("^gee_", names(g)))))
}

hdr("Slides to update in docs/mn_prediction_slides.qmd")
cat("  * Summary   - replace the AUC claim and the LOCO sentence\n")
cat("  * Results   - AUC tables and forest plots re-render from the store\n")
cat("  * anything quoting a bare r_max: state the estimator AND estimand\n")
cat("    (FINDINGS.md section 31 - the vitamin A ceiling spans 0.57-0.89\n")
cat("     depending on which is used)\n")
