# =============================================================================
# scripts/covariates/02_build_country_predictors.R                  [STAGE 2]
#
# Assemble ONE full admin-2 predictor table per country from every registered
# source, in each source's NATIVE vocabulary. No harmonisation happens here --
# stage 3 does that. Keeping the two apart means the raw table stays a faithful
# record of what each provider actually delivered, and any harmonisation
# decision can be re-run or revised without re-downloading anything.
#
# Sources come from metadata/covariates/admin2_sources.csv, so adding one is a
# row in a CSV plus (if it needs downloading) a stage-1 script.
#
# WHY THIS EXISTS
# ---------------
# The area-level model previously saw ONLY the Earth-Engine block: the DHS,
# MICS, SoilGrids, MapSPAM and ESPEN admin-2 aggregates were merged into the
# INDIVIDUAL-level data and never reached the primary area-level estimator.
# This stage puts every admin-2 source in one place so that is a modelling
# choice rather than an accident of plumbing.
#
# Output: data/covariates/country/<Country>_predictors_admin2_raw.csv
#         data/covariates/country/<Country>_source_manifest.csv
#
#   Rscript scripts/covariates/02_build_country_predictors.R          # all
#   Rscript scripts/covariates/02_build_country_predictors.R Ghana    # one
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
# tar_source() rather than a couple of explicit source() calls: admin2_analysis.R
# depends on the covariates helpers (harmonized_vocab_active), and sourcing a
# hand-picked subset made the legacy extraction fail with "could not find
# function" -- swallowed by a tryCatch and reported only as "gee failed".
targets::tar_source(here("R"))

COUNTRIES <- list(
  Gambia      = list(key = "gambia",      code = "GMB", dhs_name = "Gambia",       dhs_year = 2019L),
  Ghana       = list(key = "ghana",       code = "GHA", dhs_name = "Ghana",        dhs_year = 2014L),
  SierraLeone = list(key = "sierraleone", code = "SLE", dhs_name = "Sierra Leone", dhs_year = 2013L),
  Malawi      = list(key = "malawi",      code = "MWI", dhs_name = "Malawi",       dhs_year = 2015L),
  Tanzania    = list(key = "tanzania",    code = "TZA", dhs_name = "Tanzania",     dhs_year = 2010L)
)

TARGETS_STORE <- here("_targets_full")
OUT_DIR <- here("data", "covariates", "country")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sources <- readr::read_csv(here("metadata", "covariates", "admin2_sources.csv"),
                           show_col_types = FALSE) %>% filter(enabled)

expand_path <- function(tmpl, cc, country) {
  tmpl <- gsub("\\{key\\}", cc$key, tmpl)
  tmpl <- gsub("\\{code\\}", cc$code, tmpl)
  tmpl <- gsub("\\{country\\}", country, tmpl)
  tmpl <- gsub("\\{dhs_name\\}", cc$dhs_name, tmpl)
  gsub("\\{dhs_year\\}", cc$dhs_year, tmpl)
}

load_source <- function(src, cc, country) {
  path <- expand_path(src$path_template, cc, country)
  if (identical(src$kind, "csv")) {
    f <- here(path)
    if (!file.exists(f)) return(list(data = NULL, status = "missing", detail = basename(path)))
    d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
    return(list(data = as.data.frame(d), status = "ok", detail = basename(path)))
  }
  if (identical(src$kind, "rds")) {
    f <- here(path)
    if (!file.exists(f)) return(list(data = NULL, status = "missing", detail = basename(path)))
    return(list(data = as.data.frame(readRDS(f)), status = "ok", detail = basename(path)))
  }
  if (identical(src$kind, "targets_or_live")) {
    # IMMUTABLE SOURCE, ON PURPOSE.
    #
    # This originally read gee_admin2_<country> from the targets store. That is
    # a CIRCULAR input: the pipeline writes those targets, so once a harmonised
    # run has executed, the store holds CANONICAL names and stage 2 would
    # re-harmonise its own output. It happened on 2026-08-29 -- the per-country
    # raw tables silently shrank from ~660 columns to ~400 because the "raw"
    # Earth-observation input had already been harmonised.
    #
    # The raw layer must therefore come from something the pipeline never
    # rewrites: a one-time legacy extraction cached to disk. It costs ~15 min
    # per country once and is then read in milliseconds forever.
    cache <- here("data", "covariates", "raw",
                  sprintf("%s_gee_admin2_legacy.csv", country))
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    if (file.exists(cache)) {
      d <- suppressMessages(readr::read_csv(cache, show_col_types = FALSE))
      return(list(data = as.data.frame(d), status = "ok (cached)",
                  detail = basename(cache)))
    }

    # Accept the store only if it is genuinely in the LEGACY vocabulary.
    tgt <- sub("^.*::", "", path)
    d <- tryCatch(targets::tar_read_raw(tgt, store = TARGETS_STORE), error = function(e) NULL)
    if (!is.null(d) && any(grepl("^gee_", names(d)))) {
      readr::write_csv(as.data.frame(d), cache)
      return(list(data = as.data.frame(d), status = "ok (store->cached)",
                  detail = paste0("targets:", tgt)))
    }
    if (!is.null(d))
      message("    store copy of ", tgt, " is already harmonised - ignoring it ",
              "and extracting the legacy vocabulary from source")

    message("    extracting legacy covariates for ", country, " (slow, once)")
    old <- Sys.getenv("COVARIATE_VOCAB"); Sys.setenv(COVARIATE_VOCAB = "legacy")
    on.exit(Sys.setenv(COVARIATE_VOCAB = old), add = TRUE)
    ccfg <- get_country_configs()[[country]]
    d <- tryCatch(extract_gee_admin2(ccfg),
                  error = function(e) { message("    extraction ERROR: ",
                                                conditionMessage(e)); NULL })
    if (is.null(d) || !any(grepl("^gee_", names(d))))
      return(list(data = NULL, status = "failed", detail = tgt))
    readr::write_csv(as.data.frame(d), cache)
    return(list(data = as.data.frame(d), status = "ok (live->cached)", detail = basename(cache)))
  }
  list(data = NULL, status = "unknown kind", detail = src$kind)
}

args <- commandArgs(trailingOnly = TRUE)
sel  <- if (length(args)) args else names(COUNTRIES)

for (country in sel) {
  cc <- COUNTRIES[[country]]; if (is.null(cc)) { message("unknown: ", country); next }
  message("\n=== ", country, " ===")
  base <- NULL; manifest <- list()

  for (i in seq_len(nrow(sources))) {
    src <- sources[i, ]
    res <- load_source(src, cc, country)
    n_cols <- if (is.null(res$data)) 0L else
      ncol(res$data) - sum(c("Admin1", "Admin2") %in% names(res$data))
    manifest[[length(manifest) + 1L]] <- data.frame(
      country = country, source_id = src$source_id, status = res$status,
      file = res$detail, n_areas = if (is.null(res$data)) NA_integer_ else nrow(res$data),
      n_columns = n_cols, stringsAsFactors = FALSE)
    message(sprintf("  %-11s %-10s %4s cols  %s", src$source_id, res$status,
                    n_cols, res$detail))
    if (is.null(res$data) || !n_cols) next

    d <- res$data
    if (!"Admin2" %in% names(d)) next
    d$Admin2 <- trimws(as.character(d$Admin2))
    if ("Admin1" %in% names(d)) d$Admin1 <- trimws(as.character(d$Admin1))
    # Prefix collisions across sources would silently overwrite; make them loud.
    d <- tryCatch(clean_admin2_keys(d, what = paste(country, src$source_id)),
                  error = function(e) d)

    if (is.null(base)) { base <- d; next }
    # Join on the (Admin1, Admin2) PAIR whenever both sides carry Admin1 --
    # GADM admin-2 names are not unique within a country (Malawi has four
    # same-named district pairs), and a name-only join silently fans out rows.
    by <- if (all(c("Admin1", "Admin2") %in% names(base)) &&
              all(c("Admin1", "Admin2") %in% names(d))) c("Admin1", "Admin2") else "Admin2"
    dup <- intersect(setdiff(names(d), by), names(base))
    if (length(dup)) {
      message("    ", length(dup), " duplicate column name(s) skipped: ",
              paste(utils::head(dup, 4), collapse = ", "))
      d <- d[, setdiff(names(d), dup), drop = FALSE]
    }
    n_before <- nrow(base)
    base <- dplyr::left_join(base, d, by = by)
    if (nrow(base) != n_before)
      warning(sprintf("[%s] %s join changed row count %d -> %d",
                      country, src$source_id, n_before, nrow(base)))
  }

  if (is.null(base)) { message("  no sources resolved; skipping"); next }
  out <- file.path(OUT_DIR, sprintf("%s_predictors_admin2_raw.csv", country))
  readr::write_csv(base, out)
  readr::write_csv(dplyr::bind_rows(manifest),
                   file.path(OUT_DIR, sprintf("%s_source_manifest.csv", country)))
  message(sprintf("  -> %s (%d areas x %d predictors)", basename(out), nrow(base),
                  ncol(base) - sum(c("Admin1", "Admin2") %in% names(base))))
}
message("\nStage 2 complete.")
