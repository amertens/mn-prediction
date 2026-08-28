# =============================================================================
# scripts/build_transport_bakeoff.R
#
# Consolidate every AREA-LEVEL leave-one-country-out transport result into one
# table, so the model families can be read side by side.
#
# WHY A CONSOLIDATED TABLE NEEDS CARE
# -----------------------------------
# The source files were produced by separate experiments and do NOT share a
# reliability ceiling: loco_headline r_max 0.305, loco_bakeoff and loco_spatial
# 0.361, spatial_plus_soil_rescored 0.205. A raw Pearson r is therefore not
# comparable across files. Two things are carried so it can be read honestly:
#
#   r_max     the ceiling for that experiment's cells, i.e. the correlation
#             sampling noise alone permits
#   r_share   pearson / r_max, the fraction of the attainable correlation the
#             method achieves, which IS comparable across experiments
#
# An r_share above 1 means the method correlates better than sampling noise
# should allow on those cells, which indicates the number is not honest
# out-of-sample rather than that the method is exceptional.
#
# ANCHORING
# ---------
# Three of the four files carry an `anchor` column. `train_mean` is the honest
# setting: the map is placed using the training countries' mean. `oracle_national`
# uses the held-out country's own national value, which is unavailable in
# deployment. Both are kept and labelled; the headline should be read on
# train_mean.
#
# Run:
#   Rscript scripts/build_transport_bakeoff.R
#
# Writes results/tables/corrected/loco_transport_bakeoff.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})

SRC_DIR <- here("results", "tables", "frozen_2026-08", "sandbox_parsimony_out")
FILES <- c("loco_headline", "loco_bakeoff", "loco_spatial",
           "spatial_plus_soil_rescored")

rows <- list()
for (f in FILES) {
  p <- file.path(SRC_DIR, paste0(f, ".csv"))
  if (!file.exists(p)) { warning("missing ", p, call. = FALSE); next }
  d <- readr::read_csv(p, show_col_types = FALSE)
  if (!"anchor" %in% names(d)) d$anchor <- "not_applicable"
  rows[[f]] <- d |>
    dplyr::group_by(.data$variant, .data$anchor) |>
    dplyr::summarise(
      n_cells     = dplyr::n(),
      pearson_r   = round(mean(.data$pearson, na.rm = TRUE), 4),
      spearman_r  = round(mean(.data$spearman, na.rm = TRUE), 4),
      rmse_pp     = round(mean(.data$rmse_pp, na.rm = TRUE), 2),
      abs_level_bias_pp = round(mean(abs(.data$level_bias_pp), na.rm = TRUE), 2),
      r_max       = round(mean(.data$r_max, na.rm = TRUE), 3),
      .groups = "drop") |>
    dplyr::mutate(
      r_share = round(.data$pearson_r / .data$r_max, 3),
      source_file = paste0("sandbox_parsimony_out/", f, ".csv"))
}

if (!length(rows)) stop("No source files found.", call. = FALSE)
out <- dplyr::bind_rows(rows) |>
  dplyr::arrange(.data$source_file, .data$anchor, dplyr::desc(.data$pearson_r))

dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(out, here("results", "tables", "corrected", "loco_transport_bakeoff.csv"))

cat(sprintf("%d variant-by-anchor rows across %d source files\n",
            nrow(out), dplyr::n_distinct(out$source_file)))

# Printed one experiment at a time. Ranking across files would be misleading:
# the ceilings differ (0.205 to 0.361) on nominally the same 16 cells, so
# neither the raw correlation nor r_share is comparable between files.
for (sf in unique(out$source_file)) {
  cat(sprintf("\n=== %s, honest anchor ===\n", sf))
  print(as.data.frame(
    out |> dplyr::filter(.data$source_file == sf,
                         .data$anchor %in% c("train_mean", "not_applicable")) |>
      dplyr::select(variant, n_cells, pearson_r, r_max, r_share,
                    rmse_pp, abs_level_bias_pp) |>
      dplyr::arrange(dplyr::desc(.data$pearson_r))), row.names = FALSE)
}
cat("\nCeilings differ between files, so compare WITHIN an experiment, not across.\n")

cat("\nwrote results/tables/corrected/loco_transport_bakeoff.csv\n")
