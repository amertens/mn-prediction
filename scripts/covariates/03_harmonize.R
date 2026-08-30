# =============================================================================
# scripts/covariates/03_harmonize.R                                 [STAGE 3]
#
# Turn the per-country raw predictor tables (stage 2) into ONE harmonised
# cross-country admin-2 dataset in a single canonical vocabulary.
#
# What happens here, in order:
#   1. map every raw column onto a canonical name  (harmonization_rules.csv)
#   2. collapse the several raw columns that mean the same thing in one country
#      (year-stamped duplicates, month climatologies)
#   3. convert values into the canonical unit                (unit_conversions.csv)
#   4. drop canonical variables that cannot be reconciled          (exclusions.csv)
#   5. keep the intersection present in EVERY country, and record what fell out
#
# Outputs (data/covariates/harmonized/):
#   predictors_admin2_harmonized.csv   the analysis-ready pooled table
#   column_map.csv                     raw name -> canonical, per country, with
#                                      the reason for every keep and every drop
#   coverage.csv                       canonical variable x country availability
#
#   Rscript scripts/covariates/03_harmonize.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
source(here("R", "covariates", "canonicalize.R"))

IN_DIR  <- here("data", "covariates", "country")
OUT_DIR <- here("data", "covariates", "harmonized")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(IN_DIR, pattern = "_predictors_admin2_raw\\.csv$", full.names = TRUE)
if (!length(files)) stop("no stage-2 outputs in ", IN_DIR, " -- run 02_build_country_predictors.R first")

# WHICH COUNTRIES DEFINE THE SHARED SET
# -------------------------------------
# The intersection is taken over the countries the pipeline actually ANALYSES,
# not over every country that happens to have a stage-2 table. Tanzania has been
# ingested but is not in get_country_configs() yet, and because it has no
# dhs_custom file its presence in the intersection cost 94 predictors -- every
# DHS indicator -- for an analysis it does not take part in. A country that is
# not modelled must not be able to shrink the covariate set of one that is.
#
# Non-active countries are still harmonised and still get their own per-country
# canonical table; they are simply not allowed a veto over the shared vocabulary.
# Set HARMONIZE_COUNTRIES to override (comma-separated), e.g. when preparing a
# run that will include a newly activated country.
active <- Sys.getenv("HARMONIZE_COUNTRIES", "")
active <- if (nzchar(active)) trimws(strsplit(active, ",")[[1]]) else
  tryCatch({ targets::tar_source(here("R", "config.R"))
             vapply(get_country_configs(), function(x) x$country, character(1)) },
           error = function(e) character(0))
active <- cov_country_key(active)
if (!length(active)) {
  message("[harmonize] could not resolve active countries; using every stage-2 table")
} else {
  message("[harmonize] shared set defined by active countries: ",
          paste(active, collapse = ", "))
}

rules <- cov_load_rules(); units <- cov_load_units(); excl <- cov_load_exclusions()
maps <- list(); harm <- list()

for (f in files) {
  country <- sub("_predictors_admin2_raw\\.csv$", "", basename(f))
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>% as.data.frame()
  res <- cov_harmonize_country(d, country, rules = rules, units = units, exclusions = excl)
  res$data$country <- country
  harm[[country]] <- res$data
  maps[[country]] <- res$map
  n_can <- ncol(res$data) - sum(c("Admin1", "Admin2", "country") %in% names(res$data))
  message(sprintf("%-12s %4d raw -> %4d canonical  (%d dropped, %d unmatched, %d excluded)",
                  country, nrow(res$map), n_can,
                  sum(res$map$action == "drop"), sum(res$map$action == "unmatched"),
                  sum(res$map$excluded)))
}

map_all <- dplyr::bind_rows(maps)
readr::write_csv(map_all, file.path(OUT_DIR, "column_map.csv"))

# ── Coverage: which canonical variables exist in which countries ────────────
canon_by_country <- lapply(harm, function(x)
  setdiff(names(x), c("Admin1", "Admin2", "country")))
all_canon <- sort(unique(unlist(canon_by_country)))
cov_tab <- data.frame(canonical = all_canon, stringsAsFactors = FALSE)
for (cn in names(canon_by_country))
  cov_tab[[cn]] <- all_canon %in% canon_by_country[[cn]]
cov_tab$n_countries <- rowSums(cov_tab[, names(canon_by_country), drop = FALSE])

# `in_all` means "present in every ACTIVE country" (see the note above). The
# per-country columns still show non-active countries, so a reader can see
# exactly what activating one would cost.
active_present <- intersect(cov_country_key(names(canon_by_country)), active)
if (length(active_present)) {
  key_of <- setNames(cov_country_key(names(canon_by_country)), names(canon_by_country))
  act_cols <- names(key_of)[key_of %in% active_present]
  cov_tab$n_active <- rowSums(cov_tab[, act_cols, drop = FALSE])
  cov_tab$in_all <- cov_tab$n_active == length(act_cols)
  message(sprintf("[harmonize] intersection over %d active countries (%s)",
                  length(act_cols), paste(act_cols, collapse = ", ")))
  skipped <- setdiff(names(canon_by_country), act_cols)
  if (length(skipped)) {
    lost <- sum(cov_tab$in_all &
                  rowSums(cov_tab[, skipped, drop = FALSE]) < length(skipped))
    message(sprintf("[harmonize] %s not in the intersection; including %s would cost %d predictors",
                    paste(skipped, collapse = ", "), paste(skipped, collapse = "/"), lost))
  }
} else {
  cov_tab$n_active <- cov_tab$n_countries
  cov_tab$in_all <- cov_tab$n_countries == length(canon_by_country)
}
readr::write_csv(cov_tab, file.path(OUT_DIR, "coverage.csv"))

shared <- cov_tab$canonical[cov_tab$in_all]
message(sprintf("\ncanonical variables: %d total | %d shared by all %d countries",
                length(all_canon), length(shared), length(harm)))

# ── Pooled harmonised table: the shared vocabulary, one row per admin-2 ──────
pooled_src <- if (length(active_present))
  harm[cov_country_key(names(harm)) %in% active_present] else harm
pooled <- dplyr::bind_rows(lapply(pooled_src, function(x) {
  keep <- c(intersect(c("Admin1", "Admin2", "country"), names(x)), shared)
  x[, keep, drop = FALSE]
}))
pooled <- pooled[, c("country", intersect(c("Admin1", "Admin2"), names(pooled)), shared)]
readr::write_csv(pooled, file.path(OUT_DIR, "predictors_admin2_harmonized.csv"))
message(sprintf("-> predictors_admin2_harmonized.csv (%d areas x %d predictors)",
                nrow(pooled), length(shared)))

# Country-specific canonical variables are not lost -- they stay in the
# per-country stage-2 tables and are listed here so the exclusion is visible.
lost <- cov_tab %>% filter(!in_all) %>% arrange(desc(n_countries))
if (nrow(lost))
  message(sprintf("%d canonical variables are NOT in all countries (see coverage.csv); top: %s",
                  nrow(lost), paste(utils::head(lost$canonical, 6), collapse = ", ")))
message("\nStage 3 complete.")
