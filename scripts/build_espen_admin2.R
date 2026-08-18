#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_espen_admin2.R
#
# Build subnational soil-transmitted-helminth (STH) and schistosomiasis (SCH)
# prevalence predictors aggregated to GADM Admin-2, writing
# data/ESPEN/<Country>_espen_admin2.csv with pipeline-ready `espen_` columns.
#
# WHY: hookworm causes chronic intestinal blood loss -> iron-deficiency anaemia
# (schistosomiasis likewise). It is currently absent from the predictor set and
# is the most mechanistically-relevant missing covariate for the IRON outcomes.
# Subnational (IU ~ Admin-2) -> helps within-country AND transfer.
# See docs/ra_new_data_sources_specs.md §1.
#
# -----------------------------------------------------------------------------
# ACCESS (read this — the API route is currently closed):
#
#   The documented ESPEN Data API now REQUIRES a key, and WHO/AFRO has PAUSED
#   issuing new keys (confirmed 2026-07 via the ESPENAPI package README). So the
#   reliable route today is the portal's public, no-key "Download data" page:
#
#     PRIMARY (manual download, no key) — do this once per country:
#       1. Go to  https://espen.afro.who.int/maps-data/data-query-tools/download-data
#       2. Select: Country = <country> ; Disease = STH ; Level = IU ;
#          From/To = full available year range ; click "Download Data".
#       3. Repeat with Disease = SCH.
#       4. Save both CSVs into  data/ESPEN/raw/  (any filename; this script
#          auto-detects country + disease from the file contents/name).
#
#     OPTIONAL (only if you already hold a valid key): set ESPEN_API_KEY and
#       install the ESPENAPI package; the script will pull directly instead.
#       Request access (if reopened): email  ntd.espen@who.int  or see
#       https://espen.afro.who.int/espen-api-platform .
#
# Written defensively: it PRINTS each source's columns and the field mapping it
# chose — confirm them on first run. Output columns/paths are kept compatible
# with sandbox_fe/21_espen_iron_test.R (baseline vs +helminth).
# =============================================================================
suppressMessages({library(here); library(dplyr); library(stringdist)})

COUNTRIES <- c("Gambia","Ghana","SierraLeone","Malawi","Tanzania")
ISO2 <- c(Gambia="GM", Ghana="GH", SierraLeone="SL", Malawi="MW", Tanzania="TZ")
ISO3 <- c(Gambia="GMB",Ghana="GHA",SierraLeone="SLE",Malawi="MWI",Tanzania="TZA")
DISEASES <- c("sth","sch")
LEVEL    <- "iu"                    # implementation unit ~ Admin-2 district
RAW_DIR  <- here::here("data","ESPEN","raw")
OUT_DIR  <- here::here("data","ESPEN")
dir.create(RAW_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)) COUNTRIES <- args
bad <- setdiff(COUNTRIES, names(ISO2))
if (length(bad)) stop("Unknown country: ", paste(bad, collapse=", "))

KEY <- Sys.getenv("ESPEN_API_KEY")
USE_API <- nzchar(KEY) && requireNamespace("ESPENAPI", quietly=TRUE)

# ---- field pickers (VERIFY on first run) -----------------------------------
pick_prev <- function(df) {
  cand <- grep("preval|^prev$|_prev|Prev_", names(df), ignore.case=TRUE, value=TRUE)
  cand <- cand[vapply(cand, function(c){
    v <- suppressWarnings(as.numeric(as.character(df[[c]]))); mean(is.finite(v)) > 0.5
  }, logical(1))]
  if (!length(cand)) return(NULL)
  # prefer an "any-species" / 1+ prevalence field if present
  pref <- grep("1plus|any|infection|_prev$|prevalence", cand, ignore.case=TRUE, value=TRUE)
  col  <- if (length(pref)) pref[1] else cand[1]
  cat("      prevalence field ->", col, "\n")
  v <- suppressWarnings(as.numeric(as.character(df[[col]])))
  if (max(v, na.rm=TRUE) > 1.5) v <- v/100   # percent -> proportion
  v
}
pick_admin2 <- function(df) {
  cand <- grep("IU_?NAME|IUs_NAME|ADMIN2|Admin2|IUName|district|^IU$", names(df),
               ignore.case=TRUE, value=TRUE)
  if (!length(cand)) cand <- grep("name", names(df), ignore.case=TRUE, value=TRUE)
  cat("      admin2/IU field  ->", cand[1], "\n")
  as.character(df[[cand[1]]])
}

match_to_gadm <- function(esp_names, gadm_names, thr=0.15) {
  if (!length(gadm_names)) return(esp_names)
  m <- stringdist::stringdistmatrix(tolower(esp_names), tolower(gadm_names), method="jw")
  best <- apply(m, 1, which.min); dist <- apply(m, 1, min)
  out <- gadm_names[best]; out[dist > thr] <- NA
  out
}

# ---- source readers --------------------------------------------------------
read_api <- function(cn, dz) {
  tryCatch(ESPENAPI::ESPEN_API_data(iso2=ISO2[[cn]], disease=dz, level=LEVEL),
           error=function(e){message("  ESPENAPI failed (",cn,"/",dz,"): ",conditionMessage(e)); NULL})
}
# match a downloaded raw CSV to (country, disease) by filename or iso content
read_raw <- function(cn, dz) {
  files <- list.files(RAW_DIR, pattern="\\.csv$", full.names=TRUE, ignore.case=TRUE)
  if (!length(files)) return(NULL)
  score <- function(f) {
    nm <- tolower(basename(f))
    s <- 0
    if (grepl(tolower(dz), nm)) s <- s + 2
    if (grepl(tolower(cn), nm) || grepl(tolower(ISO3[[cn]]), nm) ||
        grepl(paste0("\\b", tolower(ISO2[[cn]]), "\\b"), nm)) s <- s + 2
    s
  }
  cand <- files[order(-vapply(files, score, numeric(1)))]
  for (f in cand) {
    df <- tryCatch(utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE),
                   error=function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    # confirm country + disease from contents when columns allow
    cc <- names(df)
    iso_col <- cc[grep("iso|country|ADMIN0", cc, ignore.case=TRUE)][1]
    ok_country <- is.na(iso_col) ||
      any(grepl(paste(c(cn, ISO2[[cn]], ISO3[[cn]]), collapse="|"),
                df[[iso_col]], ignore.case=TRUE))
    dz_col <- cc[grep("disease|type", cc, ignore.case=TRUE)][1]
    ok_dz <- is.na(dz_col) || any(grepl(dz, df[[dz_col]], ignore.case=TRUE)) ||
             grepl(dz, tolower(basename(f)))
    if (ok_country && ok_dz) { cat("    using raw file:", basename(f), "\n"); return(df) }
  }
  NULL
}

source(here::here("R","data_prep.R"))  # load_gadm_cached()

for (cn in COUNTRIES) {
  cat("\n#####", cn, "#####\n")
  gadm <- tryCatch(load_gadm_cached(ISO3[[cn]], level=2), error=function(e) NULL)
  gadm_a2 <- if (!is.null(gadm)) unique(as.character(sf::st_drop_geometry(gadm)$NAME_2)) else character()
  per_disease <- list()
  for (dz in DISEASES) {
    cat("  ", toupper(dz), ":\n", sep="")
    df <- if (USE_API) read_api(cn, dz) else read_raw(cn, dz)
    if (is.null(df) && !USE_API) df <- read_api(cn, dz)   # last resort if key exists
    if (is.null(df) || !nrow(df)) { cat("    no data (download from portal into data/ESPEN/raw/)\n"); next }
    cat(sprintf("    %d rows; columns -> %s\n", nrow(df),
                paste(head(names(df),12), collapse=", ")))
    prev <- pick_prev(df); a2 <- pick_admin2(df)
    if (is.null(prev) || is.null(a2)) { cat("    !! verify field mapping\n"); next }
    sub <- data.frame(iu=a2, prev=prev, stringsAsFactors=FALSE)
    sub <- sub[is.finite(sub$prev), ]
    if (!nrow(sub)) { cat("    no finite prevalence values\n"); next }
    agg <- aggregate(prev ~ iu, sub, mean)        # collapse surveys/years per IU
    agg$Admin2 <- match_to_gadm(agg$iu, gadm_a2)
    n_match <- sum(!is.na(agg$Admin2))
    cat(sprintf("    IU->GADM matched %d/%d\n", n_match, nrow(agg)))
    agg <- agg[!is.na(agg$Admin2), c("Admin2","prev")]
    if (!nrow(agg)) next
    # Collapse multiple IUs mapping to the same GADM Admin2 -> one row per Admin2.
    # ESPEN carries no Admin1, so the pipeline join is on Admin2 alone; a unique
    # key here prevents row expansion for countries with repeated Admin2 names.
    agg <- aggregate(prev ~ Admin2, agg, mean)
    names(agg)[2] <- paste0("espen_", dz, "_prev")
    per_disease[[dz]] <- agg
  }
  if (!length(per_disease)) { cat("  nothing written for", cn, "\n"); next }
  out <- Reduce(function(a,b) merge(a,b,by="Admin2",all=TRUE), per_disease)
  if ("espen_sth_prev" %in% names(out)) out$helminth_sth_any_prev <- out$espen_sth_prev
  fn <- file.path(OUT_DIR, sprintf("%s_espen_admin2.csv", cn))
  write.csv(out, fn, row.names=FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows, %d espen_ cols)\n",
              basename(fn), nrow(out), sum(grepl("^espen_|^helminth_", names(out)))))
}
cat("\nDONE. Then: Rscript sandbox_fe/21_espen_iron_test.R  (baseline vs +helminth)\n")
cat("Config domain ESPEN=list(prefix=\"espen_\") and the guarded merge join are already wired.\n")
