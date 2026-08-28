# =============================================================================
# scripts/run_uniform_brinda.R
#
# WS3, deliverable 1. Enumerate what inflammation adjustment each
# country-outcome actually receives, then re-derive the iron outcomes under one
# uniform protocol and measure how much of the national LEVEL bias moves.
#
# Why iron and not vitamin A: vitamin A was already harmonized in DC-H2.
# apply_brinda_vita_binary() overwrites every country's VAD binary with one
# BRINDA CRP+AGP definition, and brinda_country_method() returns
# "BRINDA CRP+AGP" for all four active countries, so a uniform-adjustment
# sensitivity for vitamin A is a no-op. Iron is the outcome where
# UNIFORM_TRANSPORT_TAGS applies a uniform CUTOFF to four differently ADJUSTED
# ferritin columns.
#
# The uniform scheme is applied by pointing compute_svy_admin2() at a
# newly-derived adjusted-ferritin column with the WHO cutoff on the original
# scale, so the survey design, weighting and domain handling are byte-identical
# between the two schemes and only the biomarker adjustment differs.
#
# Run:
#   Rscript scripts/run_uniform_brinda.R
#
# Writes
#   metadata/adjustment_inventory.csv
#   results/tables/corrected/uniform_brinda_prevalence.csv
#   results/tables/corrected/uniform_brinda_loco_level_bias.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

params    <- get_pipeline_params()
configs   <- get_country_configs()
countries <- names(configs)
IRON      <- c("child_iron", "women_iron")

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

# ── 1. Adjustment inventory ─────────────────────────────────────────────────
inv <- adjustment_inventory(configs)
dir.create(here("metadata"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(inv, here("metadata", "adjustment_inventory.csv"))
cat("\n=== adjustment inventory ===\n")
print(as.data.frame(inv[, c("country", "outcome", "configured_continuous",
                            "adjustment_harmonized_today", "uniform_method")]),
      row.names = FALSE)

# ── 2. Re-derive iron under the uniform protocol ────────────────────────────
# The adjustment is computed on the MERGED data, not on outcome_data.
# build_outcome_dataset() treats other outcomes' biomarker columns as leakage
# and drops them, so the raw ferritin column survives into outcome_data only for
# Gambia. Deriving from merged and rebuilding the outcome dataset keeps every
# country on the same footing and keeps the population filter, predictor
# pruning and column bookkeeping identical to production.

#' Return outcome_data rebuilt with a uniform-BRINDA adjusted ferritin column,
#' and the modified outcome config that points compute_svy_admin2() at it.
.uniform_iron <- function(merged, cc, oc) {
  pop  <- if (grepl("^child", oc$tag)) "child" else "women"
  spec <- brinda_ferritin_cols(cc$country)[[pop]]
  if (is.null(spec)) return(NULL)
  need <- unname(unlist(spec[intersect(names(spec), c("fer", "crp", "agp"))]))
  if (!all(need %in% colnames(merged))) {
    warning(sprintf("[uniform_brinda] %s %s: missing %s in merged data", cc$country,
                    oc$tag, paste(setdiff(need, colnames(merged)), collapse = ", ")))
    return(NULL)
  }
  newcol <- "ws3_ferritin_brinda"
  merged[[newcol]] <- brinda_adjust_ferritin(
    merged[[spec$fer]], merged[[spec$crp]],
    if (is.null(spec$agp)) NULL else merged[[spec$agp]])

  oc2 <- oc
  oc2$continuous   <- newcol
  oc2$cutoff       <- if (pop == "child") 12 else 15
  oc2$cutoff_scale <- "original"
  oc2$cutoff_dir   <- "less"
  # Setting oc2$continuous BEFORE the rebuild is what makes build_outcome_dataset()
  # keep the new column: it retains oc$continuous by name (R/data_prep.R:590) and
  # has no other reason to keep a column this script invented.
  #
  # Both arms are rebuilt here rather than read from the stored outcome_data_*
  # targets. The stored targets were built by an older code version and no longer
  # agree with the current population filter on row count, and using one stored
  # arm against one freshly built arm would compare two code versions rather than
  # two adjustments.
  od <- build_outcome_dataset(merged, cc, oc2)
  d  <- od$data

  # compute_svy_admin2() derives the binary internally but does not write it
  # back, and national_design_based() reads oc$binary off the data. Write the
  # derived binary in so the Admin-2 and national figures describe the same
  # outcome.
  d[[oc$binary]] <- apply_threshold(d[[newcol]], oc2$cutoff, "less")
  od$data <- d
  list(od = od, oc = oc2, raw_col = spec$fer, crp_col = spec$crp, agp_col = spec$agp,
       merged_pop = if (!is.null(cc$child_flag) && cc$child_flag %in% colnames(merged))
         merged[!is.na(merged[[cc$child_flag]]) &
                merged[[cc$child_flag]] == oc$child_flag_val, , drop = FALSE]
       else merged)
}

#' Apply the PRODUCTION uniform-cutoff rule to the configured arm, for the same
#' reason: without this the configured national figure would be computed from
#' the survey agency's own binary while its Admin-2 figures come from the
#' uniform cutoff, and the two arms would not be comparable.
.configured_uniform <- function(od, oc) {
  d <- od$data
  if (!is.null(oc$continuous) && oc$continuous %in% colnames(d) && !is.null(oc$cutoff))
    d[[oc$binary]] <- apply_threshold(suppressWarnings(as.numeric(d[[oc$continuous]])),
                                      oc$cutoff, oc$cutoff_dir %||% "less")
  od$data <- d
  od
}

prev_rows <- list()
svy_by_scheme <- list(configured = list(), uniform_brinda = list())

for (cn in countries) {
  cc <- configs[[cn]]
  for (otag in IRON) {
    oc <- cc$outcomes[[otag]]
    if (is.null(oc)) next
    merged <- read_target(paste0("merged_", tolower(cn)))
    if (is.null(merged)) { cat(sprintf("[skip] no merged data for %s\n", cn)); next }

    u <- .uniform_iron(merged, cc, oc)
    if (is.null(u)) next

    cat(sprintf("\n--- %s / %s ---\n", cn, otag))
    od_cfg  <- .configured_uniform(build_outcome_dataset(merged, cc, oc), oc)
    svy_cfg <- tryCatch(compute_svy_admin2(od_cfg, cc, oc),   error = function(e) NULL)
    svy_uni <- tryCatch(compute_svy_admin2(u$od,   cc, u$oc), error = function(e) NULL)
    nat_cfg <- tryCatch(national_design_based(od_cfg, cc, oc),   error = function(e) NULL)
    nat_uni <- tryCatch(national_design_based(u$od,   cc, u$oc), error = function(e) NULL)
    if (is.null(svy_cfg) || is.null(svy_uni)) next

    .nat <- function(x) if (is.null(x)) NA_real_ else
      as.numeric(if (is.list(x) && !is.null(x$prev)) x$prev else x[[1]])

    # Marker medians come from the population-filtered MERGED rows, because the
    # raw markers are dropped as leakage from the outcome dataset for three of
    # the four countries.
    mp <- u$merged_pop
    d  <- u$od$data
    prev_rows[[length(prev_rows) + 1L]] <- data.frame(
      country = cc$country, outcome = otag,
      raw_ferritin_col = u$raw_col,
      raw_ferritin_median = stats::median(suppressWarnings(as.numeric(mp[[u$raw_col]])), na.rm = TRUE),
      crp_median = stats::median(suppressWarnings(as.numeric(mp[[u$crp_col]])), na.rm = TRUE),
      agp_median = if (is.null(u$agp_col)) NA_real_ else
        stats::median(suppressWarnings(as.numeric(mp[[u$agp_col]])), na.rm = TRUE),
      adj_ferritin_median = stats::median(d[["ws3_ferritin_brinda"]], na.rm = TRUE),
      configured_continuous = oc$continuous,
      n_rows_configured = nrow(od_cfg$data), n_rows_uniform = nrow(d),
      n_admin2_configured = nrow(svy_cfg), n_admin2_uniform = nrow(svy_uni),
      admin2_mean_prev_configured = mean(svy_cfg$svy_prev, na.rm = TRUE),
      admin2_mean_prev_uniform    = mean(svy_uni$svy_prev, na.rm = TRUE),
      national_design_configured  = .nat(nat_cfg),
      national_design_uniform     = .nat(nat_uni),
      stringsAsFactors = FALSE)
    svy_by_scheme$configured[[paste(cn, otag, sep = "|")]]     <- svy_cfg
    svy_by_scheme$uniform_brinda[[paste(cn, otag, sep = "|")]] <- svy_uni
  }
}

prev <- dplyr::bind_rows(prev_rows)
dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(prev, here("results", "tables", "corrected", "uniform_brinda_prevalence.csv"))
cat("\n=== prevalence under the two schemes ===\n")
print(as.data.frame(prev[, c("country", "outcome", "raw_ferritin_median", "crp_median",
                             "adj_ferritin_median", "admin2_mean_prev_configured",
                             "admin2_mean_prev_uniform")]), row.names = FALSE)

# ── 3. Area-level LOCO under both schemes ───────────────────────────────────
gee <- list()
for (cn in countries) {
  g <- read_target(paste0("gee_admin2_", tolower(cn)))
  if (!is.null(g) && ncol(g) > 2) gee[[cn]] <- g
}

loco_rows <- list()
for (otag in IRON) {
  for (scheme in c("configured", "uniform_brinda")) {
    svy <- list()
    for (cn in countries) {
      k <- paste(cn, otag, sep = "|")
      s <- svy_by_scheme[[scheme]][[k]]
      if (!is.null(s) && nrow(s) > 0) svy[[cn]] <- s
    }
    if (length(svy) < 2) next
    p <- build_area_loco_dataset(svy, gee[names(svy)])
    d <- p$pooled_data
    if (is.null(d) || !nrow(d)) next
    if (!"n_svy" %in% names(d)) d$n_svy <- 1
    pool <- list(data = d, predictors = p$common_gee_vars,
                 countries = p$country_names, outcome = otag)
    res <- tryCatch(run_area_transport_loco(pool, AREA_TRANSPORT_RECIPE),
                    error = function(e) { message("  fit failed: ", conditionMessage(e)); NULL })
    if (is.null(res) || is.null(res$metrics)) next
    m <- res$metrics
    m$scheme  <- scheme
    m$outcome <- otag
    loco_rows[[paste(otag, scheme)]] <- m
  }
}

if (length(loco_rows)) {
  loco <- dplyr::bind_rows(loco_rows)
  readr::write_csv(loco, here("results", "tables", "corrected",
                              "uniform_brinda_loco_level_bias.csv"))
  cat("\n=== LOCO national level bias, by scheme ===\n")
  print(as.data.frame(loco |>
    dplyr::group_by(outcome, scheme) |>
    dplyr::summarise(n_folds = dplyr::n(),
                     mean_abs_nat_bias_pp = round(mean(abs(nat_bias_pp), na.rm = TRUE), 3),
                     mean_nat_bias_pp = round(mean(nat_bias_pp, na.rm = TRUE), 3),
                     mean_rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 3),
                     mean_pearson_r = round(mean(pearson_r, na.rm = TRUE), 4),
                     .groups = "drop")), row.names = FALSE)
} else {
  cat("\nNo LOCO metrics produced.\n")
}

cat("\nwrote:\n  metadata/adjustment_inventory.csv",
    "\n  results/tables/corrected/uniform_brinda_prevalence.csv",
    "\n  results/tables/corrected/uniform_brinda_loco_level_bias.csv\n")
