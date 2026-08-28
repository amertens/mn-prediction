# =============================================================================
# scripts/run_level_decomposition.R
#
# WS3, deliverable 3. Decompose each LOCO national level bias, on the continuous
# biomarker scale, into a pure country intercept and a covariate-predicted
# component. See R/corrected/p11_level_decomposition.R for the algebra.
#
# The response is the survey-weighted Admin-2 mean of the LOG harmonized
# biomarker:
#   iron       BRINDA CRP+AGP adjusted ferritin  (brinda_adjust_ferritin)
#   vitamin A  BRINDA CRP+AGP adjusted RBP       (brinda_adjust_rbp)
# Both are the uniform definitions, so the decomposition is not measuring
# adjustment heterogeneity a second time.
#
# A weighted mean is used rather than a full survey design because only point
# estimates enter the decomposition; no standard error is reported from it.
#
# Run:
#   Rscript scripts/run_level_decomposition.R
#
# Writes results/tables/corrected/level_bias_decomposition.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

configs   <- get_country_configs()
countries <- names(configs)
OUTCOMES  <- c("child_iron", "women_iron", "child_vitA", "women_vitA")

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

#' Harmonized continuous biomarker for one country-outcome, on the merged data.
#' Returns NULL when the raw markers are unavailable, never a silent fallback.
.harmonized_biomarker <- function(merged, cc, oc) {
  pop <- if (grepl("^child", oc$tag)) "child" else "women"
  if (grepl("iron", oc$tag)) {
    spec <- brinda_ferritin_cols(cc$country)[[pop]]
    if (is.null(spec)) return(NULL)
    need <- unname(unlist(spec[c("fer", "crp", "agp")]))
    if (!all(need[!is.na(need)] %in% colnames(merged))) return(NULL)
    list(value = brinda_adjust_ferritin(merged[[spec$fer]], merged[[spec$crp]],
                                        if (is.null(spec$agp)) NULL else merged[[spec$agp]]),
         marker = "adjusted ferritin (ug/L)")
  } else {
    spec <- brinda_rbp_cols(cc$country)[[pop]]
    if (is.null(spec)) return(NULL)
    need <- unname(unlist(spec[c("rbp", "crp", "agp")]))
    if (!all(need[!is.na(need)] %in% colnames(merged))) return(NULL)
    list(value = brinda_adjust_rbp(merged[[spec$rbp]], merged[[spec$crp]],
                                   if (is.null(spec$agp)) NULL else merged[[spec$agp]]),
         marker = "adjusted RBP (umol/L)")
  }
}

#' Survey-weighted Admin-2 mean of the log biomarker, shaped like the tables
#' build_area_loco_dataset() consumes.
.area_log_mean <- function(merged, cc, oc, value) {
  d <- merged
  pop_col <- cc$child_flag
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    keep  <- !is.na(d[[pop_col]]) & d[[pop_col]] == oc$child_flag_val
    d     <- d[keep, , drop = FALSE]
    value <- value[keep]
  }
  y <- suppressWarnings(as.numeric(value))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]]))
  w[!is.finite(w) | w <= 0] <- 1
  a <- as.character(d[[cc$admin2_col]])
  ok <- is.finite(y) & y > 0 & !is.na(a)
  if (sum(ok) < 30) return(NULL)
  ly <- log(y[ok]); w <- w[ok]; a <- a[ok]
  # Admin1 is carried because gee_admin2_* is now keyed on the (Admin1, Admin2)
  # PAIR and therefore has duplicated NAMES where two districts share one
  # (Malawi: 243 rows, 239 distinct names). A name-only join against it
  # multiplies rows; before the join-key migration it could not, because that
  # table was name-deduplicated. See R/admin2_key_hygiene.R.
  a1 <- if (!is.null(cc$admin1_col) && cc$admin1_col %in% colnames(d))
    as.character(d[[cc$admin1_col]])[ok] else NA_character_
  agg <- data.frame(Admin1 = a1, Admin2 = a, ly = ly, w = w,
                    stringsAsFactors = FALSE) |>
    dplyr::group_by(.data$Admin1, .data$Admin2) |>
    dplyr::summarise(svy_prev = stats::weighted.mean(.data$ly, .data$w),
                     n_svy = dplyr::n(), .groups = "drop") |>
    as.data.frame()
  agg[is.finite(agg$svy_prev), , drop = FALSE]
}

gee <- list()
for (cn in countries) {
  g <- read_target(paste0("gee_admin2_", tolower(cn)))
  if (!is.null(g) && ncol(g) > 2) gee[[cn]] <- g
}

merged_all <- setNames(lapply(countries, function(c) read_target(paste0("merged_", tolower(c)))),
                       countries)

rows <- list()
for (otag in OUTCOMES) {
  svy <- list(); markers <- character(0)
  for (cn in countries) {
    cc <- configs[[cn]]; oc <- cc$outcomes[[otag]]
    if (is.null(oc) || is.null(merged_all[[cn]])) next
    hb <- .harmonized_biomarker(merged_all[[cn]], cc, oc)
    if (is.null(hb)) { cat(sprintf("[skip] %s %s: raw markers unavailable\n", cn, otag)); next }
    a <- .area_log_mean(merged_all[[cn]], cc, oc, hb$value)
    if (is.null(a)) { cat(sprintf("[skip] %s %s: too few usable rows\n", cn, otag)); next }
    svy[[cn]] <- a
    markers <- c(markers, hb$marker)
  }
  if (length(svy) < 2) { cat(sprintf("[skip] %s: fewer than 2 countries\n", otag)); next }

  pool <- build_level_decomp_dataset(svy, gee[names(svy)])
  if (is.null(pool)) { cat(sprintf("[skip] %s: no pooled area data\n", otag)); next }
  # build_area_loco_dataset() names the response svy_prev; here it is the
  # area mean of the LOG biomarker, not a prevalence.
  pool$data$area_log_biomarker <- pool$data$svy_prev

  r <- decompose_level_bias(pool, otag)
  if (is.null(r)) { cat(sprintf("[skip] %s: decomposition produced nothing\n", otag)); next }
  r$marker <- unique(markers)[1]
  rows[[otag]] <- r
  cat(sprintf("[decomp] %s: %d folds\n", otag, nrow(r)))
}

if (!length(rows)) stop("No decomposition rows produced.", call. = FALSE)
out <- dplyr::bind_rows(rows)
dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(out, here("results", "tables", "corrected", "level_bias_decomposition.csv"))

cat("\n=== level-bias decomposition (log biomarker scale) ===\n")
print(as.data.frame(out |> dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))) |>
  dplyr::select(outcome, held_out, level_ratio_train_over_holdout,
                total_bias_log, location_offset_log, pattern_term_log,
                location_share, pearson_z)), row.names = FALSE)

cat("\n=== how much of the level bias is a pure intercept ===\n")
print(as.data.frame(out |> dplyr::group_by(outcome) |>
  dplyr::summarise(n_folds = dplyr::n(),
                   median_abs_location = round(stats::median(abs_location_log), 4),
                   median_abs_pattern  = round(stats::median(abs_pattern_log), 4),
                   median_location_share_bounded =
                     round(stats::median(location_share_bounded, na.rm = TRUE), 3),
                   folds_location_larger = sum(abs_location_log > abs_pattern_log),
                   folds_pattern_opposes = sum(pattern_opposes_location, na.rm = TRUE),
                   min_level_ratio = round(min(level_ratio_train_over_holdout), 3),
                   max_level_ratio = round(max(level_ratio_train_over_holdout), 3),
                   .groups = "drop")), row.names = FALSE)

cat("\nwrote results/tables/corrected/level_bias_decomposition.csv\n")
