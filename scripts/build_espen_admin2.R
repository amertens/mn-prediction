#!/usr/bin/env Rscript
# =============================================================================
# scripts/build_espen_admin2.R
#
# Pull subnational soil-transmitted-helminth (STH) and schistosomiasis (SCH)
# prevalence from the WHO/AFRO ESPEN portal and aggregate to GADM Admin-2,
# producing data/ESPEN/<Country>_espen_admin2.csv with pipeline-ready columns.
#
# WHY: hookworm is a leading cause of iron-deficiency anaemia and is currently
# absent from the predictor set. STH/SCH are subnational (help within-country
# AND transfer). See docs/ra_new_data_sources_specs.md.
#
# ACCESS: the ESPEN Data API now requires an API key (free):
#   1. Request a key at https://espen.afro.who.int  (or run ESPENAPI::espen_key_setup())
#   2. Set it once:   Sys.setenv(ESPEN_API_KEY = "....")  or add to ~/.Renviron
#   Endpoint (from the ESPENAPI package):
#     https://espenjapapi.afro.who.int/api/data?country=<C>&disease=<d>&level=<l>&...&api_key=<KEY>
#   diseases: sth, sch ; levels: iu (implementation unit ≈ district) or sitelevel
#
# This script prefers the ESPENAPI R package if installed, else calls the API
# directly with httr. It is written defensively: on first run, INSPECT the
# returned column names (printed) and confirm the prevalence field mapping below.
# =============================================================================
suppressMessages({library(here); library(dplyr); library(httr); library(jsonlite)
  library(stringdist)})

COUNTRIES <- c("Gambia","Ghana","SierraLeone","Malawi")
ISO2      <- c(Gambia="GM", Ghana="GH", SierraLeone="SL", Malawi="MW")
DISEASES  <- c("sth","sch")
LEVEL     <- "iu"                 # implementation unit ≈ Admin-2 district
OUT_DIR   <- here::here("data","ESPEN"); dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
KEY       <- Sys.getenv("ESPEN_API_KEY")
if (!nzchar(KEY)) stop("Set ESPEN_API_KEY (get a free key at https://espen.afro.who.int).")

# ---- low-level fetch (one country x disease) -------------------------------
espen_fetch <- function(iso2, disease, level=LEVEL) {
  if (requireNamespace("ESPENAPI", quietly=TRUE)) {
    return(tryCatch(ESPENAPI::ESPEN_API_data(iso2=iso2, disease=disease, level=level,
                                             api_key=KEY, df=TRUE),
                    error=function(e){message("  ESPENAPI failed: ",conditionMessage(e)); NULL}))
  }
  url <- sprintf("https://espenjapapi.afro.who.int/api/data?iso2=%s&disease=%s&level=%s&api_key=%s",
                 iso2, disease, level, KEY)
  r <- tryCatch(httr::GET(url, httr::timeout(60)), error=function(e) e)
  if (inherits(r,"error") || httr::http_error(r)) { message("  HTTP fail ",iso2,"/",disease); return(NULL) }
  txt <- httr::content(r, "text", encoding="UTF-8")
  out <- tryCatch(jsonlite::fromJSON(txt), error=function(e) NULL)
  if (is.list(out) && !is.data.frame(out)) out <- out[[which(vapply(out, is.data.frame, logical(1)))[1]]]
  out
}

# ---- map an ESPEN prevalence-like field; VERIFY on first run ---------------
# ESPEN IU tables typically carry an endemicity class and a prevalence value.
# Common field names: "Prevalence", "prev", "PrevalenceValue", "Endemicity".
pick_prev <- function(df) {
  cand <- grep("preval|^prev$|_prev", names(df), ignore.case=TRUE, value=TRUE)
  cand <- cand[vapply(cand, function(c) is.numeric(df[[c]]) || all(grepl("^[0-9.]*$", na.omit(df[[c]]))), logical(1))]
  if (!length(cand)) return(NULL)
  as.numeric(df[[cand[1]]])
}
pick_admin2 <- function(df) {
  cand <- grep("IU_?NAME|ADMIN2|IUs_NAME|district|Admin2|IUName", names(df), ignore.case=TRUE, value=TRUE)
  if (!length(cand)) cand <- grep("name", names(df), ignore.case=TRUE, value=TRUE)
  as.character(df[[cand[1]]])
}

# ---- fuzzy-match ESPEN IU names -> GADM Admin-2 (mirror project convention) -
match_to_gadm <- function(esp_names, gadm_names, thr=0.15) {
  m <- stringdist::stringdistmatrix(tolower(esp_names), tolower(gadm_names), method="jw")
  best <- apply(m, 1, which.min); dist <- apply(m, 1, min)
  out <- gadm_names[best]; out[dist > thr] <- NA
  out
}

source(here::here("R","data_prep.R"))  # load_gadm_cached()

for (cn in COUNTRIES) {
  cat("\n#####", cn, "#####\n")
  gadm <- tryCatch(load_gadm_cached(c(Gambia="GMB",Ghana="GHA",SierraLeone="SLE",Malawi="MWI")[[cn]], level=2),
                   error=function(e) NULL)
  gadm_a2 <- if (!is.null(gadm)) unique(as.character(gadm$NAME_2)) else character()
  per_disease <- list()
  for (dz in DISEASES) {
    df <- espen_fetch(ISO2[[cn]], dz)
    if (is.null(df) || !nrow(df)) { cat("  no data:", dz, "\n"); next }
    cat(sprintf("  %s: %d rows; columns -> %s\n", dz, nrow(df), paste(head(names(df),12),collapse=", ")))
    prev <- pick_prev(df); a2 <- pick_admin2(df)
    if (is.null(prev) || is.null(a2)) { cat("  !! verify field mapping for", dz, "\n"); next }
    sub <- data.frame(iu=a2, prev=prev, stringsAsFactors=FALSE)
    sub <- sub[is.finite(sub$prev), ]
    # collapse multiple surveys/years per IU to the mean prevalence
    agg <- aggregate(prev ~ iu, sub, mean)
    agg$Admin2 <- if (length(gadm_a2)) match_to_gadm(agg$iu, gadm_a2) else agg$iu
    agg <- agg[!is.na(agg$Admin2), c("Admin2","prev")]
    names(agg)[2] <- paste0("espen_", dz, "_prev")
    per_disease[[dz]] <- agg
  }
  if (!length(per_disease)) { cat("  nothing written for", cn, "\n"); next }
  out <- Reduce(function(a,b) merge(a,b,by="Admin2",all=TRUE), per_disease)
  # convenience: any-STH (max of sth) already in sth; add helminth_ alias for sth
  if ("espen_sth_prev" %in% names(out)) out$helminth_sth_any_prev <- out$espen_sth_prev
  fn <- file.path(OUT_DIR, sprintf("%s_espen_admin2.csv", cn))
  write.csv(out, fn, row.names=FALSE)
  cat(sprintf("  wrote %s (%d Admin-2 rows, %d cols)\n", basename(fn), nrow(out), ncol(out)-1))
}
cat("\nDONE. Then: Rscript sandbox_fe/21_espen_iron_test.R  (baseline vs +helminth)\n")
cat("To wire into the pipeline: add domain ESPEN=list(prefix=\"espen_\") in R/config.R and\n",
    "join data/ESPEN/<Country>_espen_admin2.csv on Admin2 in each src/<Country>/2_*_data_merge.R.\n")
