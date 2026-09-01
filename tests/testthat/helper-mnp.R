# Shared setup. Sourcing only the files the guards need keeps the suite fast:
# targets::tar_source("R") pulls in the modelling stack and takes tens of
# seconds, which is enough to stop anyone running the tests.
suppressPackageStartupMessages(library(here))

for (f in c("lint_admin2_joins.R", "unit_counts.R",
            "area_reliability.R", "reliability_empirical.R",
            "calibration_gate.R", "eb_district_estimator.R",
            "downscale_admin1.R", "eb_stack_target.R",
            "loco_outlyingness.R")) {
  source(here::here("R", f))
}

# is_biomarker_column() lives in R/data_prep.R alongside the whole preparation
# stack. Sourcing that file needs its dependencies, so the guard predicate is
# read out of it in isolation instead.
.mnp_load_biomarker_guard <- function() {
  if (exists("is_biomarker_column", envir = globalenv(), inherits = FALSE)) return(TRUE)
  src <- readLines(here::here("R", "data_prep.R"), warn = FALSE)
  # BIOMARKER_DERIVED_TOKENS is the vector the predicates interpolate.
  tok <- grep("^BIOMARKER_DERIVED_TOKENS", src)
  if (!length(tok)) return(FALSE)
  end_tok <- tok[1]
  while (end_tok < length(src) &&
         !grepl("\\)\\s*$", src[end_tok])) end_tok <- end_tok + 1L

  # Read one top-level function out of the source by brace balance.
  .extract <- function(name) {
    fun <- grep(paste0("^", name, " <- function"), src)
    if (!length(fun)) return(NULL)
    depth <- 0L; started <- FALSE; end_fun <- fun[1]
    for (i in fun[1]:length(src)) {
      depth <- depth + lengths(regmatches(src[i], gregexpr("\\{", src[i]))) -
                       lengths(regmatches(src[i], gregexpr("\\}", src[i])))
      if (grepl("\\{", src[i])) started <- TRUE
      if (started && depth <= 0L) { end_fun <- i; break }
    }
    src[fun[1]:end_fun]
  }

  # dhs_measurement_class() must come first: both predicates below call it.
  # WS-C4 added it, and a suite that loaded the callers without it reported
  # "could not find function" rather than a leak, which is the wrong failure.
  need <- c("dhs_measurement_class", "is_biomarker_column",
            "biomarker_column_class", "allowed_under_arm")
  code <- src[tok[1]:end_tok]
  for (nm in need) {
    got <- .extract(nm)
    if (is.null(got)) return(FALSE)
    code <- c(code, got)
  }
  eval(parse(text = paste(code, collapse = "
")), envir = globalenv())
  all(vapply(need, exists, logical(1), envir = globalenv(), inherits = FALSE))
}

#' Read an artefact the guard script writes, or skip the test that needs it.
mnp_artefact <- function(rel) {
  p <- here::here("results", "tables", rel)
  if (!file.exists(p))
    testthat::skip(paste0(rel, " absent; run scripts/accuracy_impact/ws7a_guards.R"))
  utils::read.csv(p, stringsAsFactors = FALSE)
}
