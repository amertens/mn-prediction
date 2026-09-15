# =============================================================================
# R/mns_dhs_overlap.R
#
# Same-survey overlap between a micronutrient survey and the DHS round that
# supplies its proxy predictors.
#
# THE RULE (2026-09-15). A predictor is leakage when it is measured on the same
# individuals as the outcome, i.e. when it comes from the same survey instance.
# A predictor from an independent source is not leakage whatever it measures:
# a DHS round's anaemia prevalence, UNICEF's supplementation coverage or an
# IHME anaemia surface are external information about a district, exactly what
# a country without a micronutrient survey would have. Name-based exclusion of
# "anaemia", "iron" or "vitamin A" columns from external sources was therefore
# withdrawn (LK-02); the `gw_` guard in R/data_prep.R, which covers the
# biomarker survey's own columns, is untouched.
#
# The one place the rule bites is Malawi. The MNS 2015-16 was a subsample of the
# MDHS 2015-16 -- 105 of its 850 clusters -- and 3,097 of 3,099 MNS rows link to
# a DHS person record (docs/findings/WSC4_MALAWI_BLOCK.md). Every DHS-derived
# Admin-2 aggregate for Malawi must be built from the OTHER 745 clusters, or the
# outcome individuals sit inside their own district's predictors. This applies
# to every DHS indicator, not only the blood-draw ones: shared cluster noise
# leaks too.
#
# The second place it bites (found 2026-09-15): The Gambia. The GMNS 2018 was
# fielded inside the MICS6 2018 sample - its 3,209 respondents sit in 70 of the
# 390 MICS clusters and carry the MICS household number - so every
# MICS-derived predictor for The Gambia must be built from the other 320
# clusters. (The Gambia's DHS 2019-20 is a separate sample and is unaffected.)
#
# metadata/mns_dhs_overlap_clusters.csv is the single source of truth (one row
# per country x programme x cluster to exclude; `programme` is DHS or MICS, so
# a MICS cluster number is never applied to a DHS recode), written by
# scripts/covariates/build_mns_overlap_clusters.R. Countries absent from the
# file for a programme have no overlap and nothing is filtered.
# =============================================================================

.mns_overlap_path <- function() here::here("metadata", "mns_dhs_overlap_clusters.csv")

.mns_key <- function(x) tolower(gsub("[^a-z]", "", tolower(as.character(x))))

#' DHS cluster numbers that the country's micronutrient survey re-sampled.
#'
#' @param country Country name in any of the project's spellings ("Malawi",
#'   "SierraLeone", "Sierra Leone").
#' @return Sorted integer vector; `integer(0)` when there is no overlap.
#' @param programme "DHS" (default) or "MICS": which survey programme's cluster
#'   numbers are wanted. Cluster numbers are only meaningful within a programme.
mns_overlap_clusters <- function(country, programme = c("DHS", "MICS"), path = .mns_overlap_path()) {
  programme <- match.arg(programme)
  if (!file.exists(path))
    stop("overlap table not found: ", path,
         " (run scripts/covariates/build_mns_overlap_clusters.R)")
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"programme" %in% names(m)) m$programme <- "DHS"
  m <- m[.mns_key(m$country) == .mns_key(country) & toupper(m$programme) == programme, , drop = FALSE]
  sort(unique(as.integer(m$dhs_cluster)))
}

#' Drop the rows of a DHS recode (or surveyPrev indicator table) that belong to
#' the micronutrient survey's own clusters.
#'
#' A no-op for countries without overlap and for frames that do not carry the
#' cluster column, so it can be applied uniformly at every DHS ingest point.
#'
#' @param df A data frame holding one row per respondent, household or cluster.
#' @param country Country name.
#' @param cluster_col Name of the DHS cluster column (`v001`, `hv001`,
#'   `cluster`, `DHSCLUST`, ...).
#' @param quiet Suppress the one-line report.
drop_mns_overlap <- function(df, country, cluster_col, quiet = FALSE, programme = "DHS") {
  if (is.null(df) || !cluster_col %in% names(df)) return(df)
  cl <- mns_overlap_clusters(country, programme = programme)
  if (!length(cl)) return(df)
  v <- df[[cluster_col]]
  v <- if (is.factor(v)) as.character(v) else unclass(v)   # haven_labelled -> plain
  v <- suppressWarnings(as.integer(as.numeric(v)))
  hit <- !is.na(v) & v %in% cl
  if (!quiet)
    message(sprintf("[mns-overlap] %s (%s): dropped %d of %d rows in %d cluster(s) (%s) sampled by the micronutrient survey",
                    country, programme, sum(hit), nrow(df), length(unique(v[hit])), cluster_col))
  df[!hit, , drop = FALSE]
}
