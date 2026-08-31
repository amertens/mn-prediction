# Shared setup. Sourcing only the files the guards need keeps the suite fast:
# targets::tar_source("R") pulls in the modelling stack and takes tens of
# seconds, which is enough to stop anyone running the tests.
suppressPackageStartupMessages(library(here))

for (f in c("lint_admin2_joins.R", "unit_counts.R")) {
  source(here::here("R", f))
}

# is_biomarker_column() lives in R/data_prep.R alongside the whole preparation
# stack. Sourcing that file needs its dependencies, so the guard predicate is
# read out of it in isolation instead.
.mnp_load_biomarker_guard <- function() {
  if (exists("is_biomarker_column", envir = globalenv(), inherits = FALSE)) return(TRUE)
  src <- readLines(here::here("R", "data_prep.R"), warn = FALSE)
  # BIOMARKER_DERIVED_TOKENS is the vector the predicate interpolates.
  tok <- grep("^BIOMARKER_DERIVED_TOKENS", src)
  fun <- grep("^is_biomarker_column <- function", src)
  if (!length(tok) || !length(fun)) return(FALSE)
  end_tok <- tok[1]
  while (end_tok < length(src) &&
         !grepl("\\)\\s*$", src[end_tok])) end_tok <- end_tok + 1L
  end_fun <- fun[1]
  depth <- 0L; started <- FALSE
  for (i in fun[1]:length(src)) {
    depth <- depth + lengths(regmatches(src[i], gregexpr("\\{", src[i]))) -
                     lengths(regmatches(src[i], gregexpr("\\}", src[i])))
    if (grepl("\\{", src[i])) started <- TRUE
    if (started && depth <= 0L) { end_fun <- i; break }
  }
  eval(parse(text = paste(c(src[tok[1]:end_tok], src[fun[1]:end_fun]), collapse = "\n")),
       envir = globalenv())
  exists("is_biomarker_column", envir = globalenv(), inherits = FALSE)
}

#' Read an artefact the guard script writes, or skip the test that needs it.
mnp_artefact <- function(rel) {
  p <- here::here("results", "tables", rel)
  if (!file.exists(p))
    testthat::skip(paste0(rel, " absent; run scripts/accuracy_impact/ws7a_guards.R"))
  utils::read.csv(p, stringsAsFactors = FALSE)
}
