# =============================================================================
# R/malawi_outcomes.R
#
# Binary deficiency indicators for the Malawi MNS 2015-16 outcomes that the
# survey file carries only as continuous biomarkers. Used by load_merged_data()
# (R/data_prep.R) the way folate_def and b12_def are derived there, and by
# scripts/protocol_v2/61_malawi_selenium_iodine.R.
#
#   sel_def  plasma selenium < 84.6 ug/L (1.07 umol/L, optimal GPX3 activity,
#            Thomson 2004; the threshold of the MNS analyses, Phiri et al. 2019)
#   iod_def  urinary iodine < 100 ug/L (WHO insufficient intake, non-pregnant
#            women). WHO classifies a population by its MEDIAN UIC; the
#            individual share below 100 is the convention IO-01 used.
# The cut-offs live in R/config.R (get_country_configs()$Malawi$outcomes).
# =============================================================================

#' Derive a less-than binary indicator from a continuous column, leaving an
#' existing indicator column untouched and NA where the biomarker is NA.
derive_malawi_binary <- function(d, cont_col, bin_col, cutoff) {
  if (bin_col %in% names(d) || !cont_col %in% names(d)) return(d)
  x <- suppressWarnings(as.numeric(unclass(d[[cont_col]])))
  x[is.finite(x) & x <= 0] <- NA_real_          # sentinel codes (-100, 0) are not measurements
  b <- ifelse(x < cutoff, 1L, 0L); b[is.na(x)] <- NA_integer_
  d[[bin_col]] <- b
  d
}
