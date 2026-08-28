# =============================================================================
# scripts/run_anchored_transport.R
#
# WS7. Anchored transport, new-country prediction intervals, and exceedance
# probabilities. See R/corrected/p13_anchored_transport.R.
#
# Run:
#   Rscript scripts/run_anchored_transport.R
#
# Writes
#   results/tables/corrected/anchored_transport_loco.csv
#   results/tables/corrected/new_country_intervals.csv
#   results/tables/corrected/exceedance_probabilities.csv
#   dashboard/data/exceedance_probabilities.rds   (data only; no module changes)
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(tidyr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

params    <- get_pipeline_params()
PROFILE   <- params$analysis_profile
configs   <- get_country_configs()
countries <- names(configs)
OUTCOMES  <- if (PROFILE == "smoke") "child_iron" else
  c("child_vitA", "women_vitA", "child_iron", "women_iron")

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

anchors_path <- here("metadata", "anchors", "national_anchors.csv")
if (!file.exists(anchors_path))
  stop("Run scripts/build_anchors.R first.", call. = FALSE)
anchors <- readr::read_csv(anchors_path, show_col_types = FALSE)

gee <- list()
for (cn in countries) {
  g <- read_target(paste0("gee_admin2_", tolower(cn)))
  if (!is.null(g) && ncol(g) > 2) gee[[cn]] <- g
}

cat(sprintf("\n[ws7] profile=%s outcomes=%s\n", PROFILE, paste(OUTCOMES, collapse = ", ")))

# ── 1. Anchored LOCO ────────────────────────────────────────────────────────
loco_rows <- list()
for (otag in OUTCOMES) {
  svy <- list()
  for (cn in countries) {
    s <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), otag))
    if (!is.null(s) && nrow(s) > 0) svy[[cn]] <- s
  }
  if (length(svy) < 2) { cat(sprintf("[skip] %s: <2 countries\n", otag)); next }
  p <- build_area_loco_dataset(svy, gee[names(svy)])
  d <- p$pooled_data
  if (is.null(d) || !nrow(d)) next
  if (!"n_svy" %in% names(d)) d$n_svy <- 1
  pool <- list(data = d, predictors = p$common_gee_vars,
               countries = p$country_names, outcome = otag)
  r <- tryCatch(run_anchored_loco(pool, otag, anchors),
                error = function(e) { message("  ERROR ", otag, ": ",
                                              conditionMessage(e)); NULL })
  if (is.null(r)) next
  loco_rows[[otag]] <- r
  cat(sprintf("[ok]   %-12s %d rows\n", otag, nrow(r)))
}

if (!length(loco_rows)) stop("No anchored-LOCO results.", call. = FALSE)
loco <- dplyr::bind_rows(loco_rows)
dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
suffix <- if (PROFILE == "smoke") "_SMOKE" else ""
readr::write_csv(loco, here("results", "tables", "corrected",
                            sprintf("anchored_transport_loco%s.csv", suffix)))

cat("\n=== anchored transport, by anchor source ===\n")
print(as.data.frame(loco |> dplyr::group_by(anchor_source) |>
  dplyr::summarise(n_folds = dplyr::n(),
                   mean_abs_nat_bias_pp = round(mean(abs(nat_bias_pp)), 3),
                   mean_mae_pp = round(mean(mae_pp), 3),
                   mean_rmse_pp = round(mean(rmse_pp), 3),
                   mean_pearson = round(mean(pearson_r, na.rm = TRUE), 4),
                   mean_spearman = round(mean(spearman_r, na.rm = TRUE), 4),
                   .groups = "drop")), row.names = FALSE)
cat(sprintf("\nranking preserved by the shift in every fold: %s\n",
            all(loco$ranking_preserved)))

# ── 2. New-country prediction intervals ─────────────────────────────────────
# Uses the same area-level log-biomarker response as the WS3 decomposition, so
# the interval describes the quantity WS3 showed dominates the transport error.
merged_all <- setNames(lapply(countries, function(c) read_target(paste0("merged_", tolower(c)))),
                       countries)
int_rows <- list()
for (otag in OUTCOMES) {
  frames <- list()
  for (cn in countries) {
    cc <- configs[[cn]]; oc <- cc$outcomes[[otag]]
    if (is.null(oc) || is.null(merged_all[[cn]])) next
    bio <- dist_continuous_biomarker(merged_all[[cn]], cc, oc)
    if (is.null(bio)) next
    m <- merged_all[[cn]]; val <- bio$value
    if (!is.null(cc$child_flag) && cc$child_flag %in% names(m)) {
      k <- !is.na(m[[cc$child_flag]]) & m[[cc$child_flag]] == oc$child_flag_val
      m <- m[k, , drop = FALSE]; val <- val[k]
    }
    w <- suppressWarnings(as.numeric(m[[cc$weight_col]])); w[!is.finite(w) | w <= 0] <- 1
    a <- as.character(m[[cc$admin2_col]])
    ok <- is.finite(val) & val > 0 & !is.na(a)
    if (sum(ok) < 30) next
    frames[[cn]] <- data.frame(country = cc$country, Admin2 = a[ok],
                               ly = log(val[ok]), w = w[ok],
                               stringsAsFactors = FALSE) |>
      dplyr::group_by(.data$country, .data$Admin2) |>
      dplyr::summarise(area_log_biomarker = stats::weighted.mean(.data$ly, .data$w),
                       .groups = "drop") |> as.data.frame()
  }
  if (length(frames) < 3) next
  ci <- new_country_interval(dplyr::bind_rows(frames))
  if (is.null(ci)) next
  ci$outcome <- otag
  int_rows[[otag]] <- ci
}
if (length(int_rows)) {
  ints <- dplyr::bind_rows(int_rows)
  readr::write_csv(ints, here("results", "tables", "corrected",
                              sprintf("new_country_intervals%s.csv", suffix)))
  cat("\n=== new-country prediction interval for the national LEVEL (log biomarker) ===\n")
  print(as.data.frame(ints |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))) |>
    dplyr::select(outcome, estimator, n_countries, sd_between,
                  lo_ratio_scale, centre_ratio_scale, hi_ratio_scale, boundary_fit)),
    row.names = FALSE)
  cat("\nThe *_ratio_scale columns are exp() of the log-scale interval, so they are in\n",
      "the marker's own units: the range a NEW country's geometric-mean biomarker\n",
      "level would be expected to fall in (ferritin ug/L for iron, RBP umol/L for\n",
      "vitamin A).\n", sep = "")
}

# ── 3. Exceedance probabilities ─────────────────────────────────────────────
thr_path <- here("metadata", "who_severity_thresholds.csv")
if (file.exists(thr_path)) {
  thr <- readr::read_csv(thr_path, show_col_types = FALSE)
  # Every threshold in the config carries a source citation and a
  # verified_against_source flag. Nothing here was taken from memory without a
  # citation, and nothing has been checked against the source PDF in this
  # session, so the flag is carried onto every output row.
  ex_rows <- list()
  for (otag in intersect(OUTCOMES, unique(thr$outcome))) {
    for (cn in countries) {
      s <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), otag))
      if (is.null(s) || !nrow(s)) next
      if (!all(c("svy_prev", "svy_prev_se") %in% names(s))) next
      tt <- thr[thr$outcome == otag, , drop = FALSE]
      for (i in seq_len(nrow(tt))) {
        ex_rows[[length(ex_rows) + 1L]] <- data.frame(
          country = configs[[cn]]$country, outcome = otag,
          Admin2 = as.character(s$Admin2),
          prevalence = s$svy_prev, se = s$svy_prev_se,
          threshold = tt$threshold[i], severity_label = tt$severity_label[i],
          p_exceeds = exceedance_probability(s$svy_prev, s$svy_prev_se, tt$threshold[i]),
          source_document = tt$source_document[i],
          verified_against_source = tt$verified_against_source[i],
          basis = "design-based Admin-2 estimate and its standard error",
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(ex_rows)) {
    ex <- dplyr::bind_rows(ex_rows)
    readr::write_csv(ex, here("results", "tables", "corrected",
                              sprintf("exceedance_probabilities%s.csv", suffix)))
    dir.create(here("dashboard", "data"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(ex, here("dashboard", "data",
                     sprintf("exceedance_probabilities%s.rds", suffix)))
    cat(sprintf("\n=== exceedance: %d rows across %d country-outcome-threshold combinations ===\n",
                nrow(ex), dplyr::n_distinct(paste(ex$country, ex$outcome, ex$threshold))))
    print(as.data.frame(ex |> dplyr::group_by(outcome, threshold, severity_label) |>
      dplyr::summarise(n_districts = dplyr::n(),
                       mean_p_exceeds = round(mean(p_exceeds, na.rm = TRUE), 3),
                       districts_p_gt_0.8 = sum(p_exceeds > 0.8, na.rm = TRUE),
                       .groups = "drop")), row.names = FALSE)
    cat("\nEvery threshold row carries verified_against_source = FALSE: the values are\n",
        "transcribed from the cited documents but were not checked against them here.\n", sep = "")
  }
}

cat("\nwrote results/tables/corrected/anchored_transport_loco", suffix, ".csv and companions\n", sep = "")
