# =============================================================================
# scripts/build_assay_lineage.R
#
# WS3, deliverable 2. Generate metadata/assay_lineage.csv: one row per
# survey x biomarker x population, recording the specimen, instrument or kit,
# laboratory, adjustment applied, cut-point and source document.
#
# The table is GENERATED, not hand-written, so it cannot drift from the code it
# describes. Everything that can be read off R/config.R and
# R/brinda_adjustment.R is read off them. Everything else is written as UNKNOWN.
#
# Nothing here is inferred. A cell is UNKNOWN unless a file in this repository
# states it. Instrument, kit and laboratory are not recorded anywhere in the
# repository for any of the four surveys, so they are UNKNOWN throughout, and
# the per-cell source column says where a reader would have to go to fill them
# in. Populating them requires reading the survey reports (for example
# data/IPD/Gambia/Gambia_Child data_Codebook.pdf), which is a data-entry task
# with a human in the loop, not something to guess at.
#
# Run:
#   Rscript scripts/build_assay_lineage.R
#
# Writes metadata/assay_lineage.csv
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
targets::tar_source(here("R"))

UNKNOWN <- "UNKNOWN"

# Which biomarker underlies each outcome family, and in what units the cut-point
# is expressed. Read from the configs where possible; the biomarker identity is
# a property of the outcome name and is declared here once.
BIOMARKER <- list(
  vitA   = list(marker = "retinol-binding protein (RBP)", units = "umol/L"),
  iron   = list(marker = "serum ferritin",                units = "ug/L"),
  folate = list(marker = "serum or plasma folate",        units = "nmol/L"),
  b12    = list(marker = "serum vitamin B12",             units = "pmol/L"),
  zinc   = list(marker = "serum or plasma zinc",          units = "ug/dL")
)

.family <- function(tag) {
  for (f in names(BIOMARKER)) if (grepl(f, tag, ignore.case = TRUE)) return(f)
  NA_character_
}

configs <- get_country_configs()
rows <- list()

for (cn in names(configs)) {
  cc <- configs[[cn]]
  for (on in names(cc$outcomes)) {
    oc  <- cc$outcomes[[on]]
    fam <- .family(oc$tag)
    pop <- if (grepl("^child", oc$tag)) "child" else "women"

    # Adjustment: read from the code that actually applies it.
    if (identical(fam, "vitA")) {
      adjustment <- brinda_country_method(cc$country)
      adj_src    <- "R/brinda_adjustment.R: brinda_country_method(); applied by apply_brinda_vita_binary() and compute_svy_admin2()"
      spec       <- brinda_rbp_cols(cc$country)[[pop]]
      raw_col    <- spec$rbp %||% UNKNOWN
    } else if (identical(fam, "iron")) {
      # The configured iron outcome is NOT harmonized; the configured continuous
      # column carries whatever the survey agency applied. WS3 adds a uniform
      # alternative but does not make it the default.
      adjustment <- paste0("survey-provided, not harmonized (configured column ",
                           oc$continuous, "); uniform alternative available: ",
                           brinda_iron_method(cc$country, pop))
      adj_src    <- "R/config.R outcome definition; uniform alternative in R/brinda_adjustment.R: brinda_adjust_ferritin()"
      spec       <- brinda_ferritin_cols(cc$country)[[pop]]
      raw_col    <- spec$fer %||% UNKNOWN
    } else {
      adjustment <- "none defined in this pipeline"
      adj_src    <- "R/config.R outcome definition (no inflammation adjustment configured)"
      raw_col    <- UNKNOWN
    }

    rows[[length(rows) + 1L]] <- data.frame(
      survey            = sprintf("%s %d", cc$country, cc$survey_year),
      country           = cc$country,
      survey_year       = cc$survey_year,
      outcome_tag       = oc$tag,
      population        = oc$population,
      biomarker         = if (is.na(fam)) UNKNOWN else BIOMARKER[[fam]]$marker,

      # Specimen: venous vs capillary vs dried blood spot is a real
      # comparability question and the repository does not record it.
      specimen          = UNKNOWN,
      source_specimen   = "not recorded in this repository; see the survey report or codebook under data/IPD/<country>/",

      instrument_or_kit = UNKNOWN,
      source_instrument = "not recorded in this repository; see the survey report or codebook under data/IPD/<country>/",

      laboratory        = UNKNOWN,
      source_laboratory = "not recorded in this repository; see the survey report or codebook under data/IPD/<country>/",

      adjustment_applied = adjustment,
      source_adjustment  = adj_src,

      cutoff            = oc$cutoff %||% NA_real_,
      cutoff_scale      = oc$cutoff_scale %||% UNKNOWN,
      cutoff_units      = if (is.na(fam)) UNKNOWN else BIOMARKER[[fam]]$units,
      source_cutoff     = "R/config.R outcome definition",

      raw_column        = raw_col,
      configured_continuous_column = oc$continuous %||% UNKNOWN,
      configured_binary_column     = oc$binary %||% UNKNOWN,
      uniform_cutoff_applied = oc$tag %in% UNIFORM_TRANSPORT_TAGS,
      source_columns    = "R/config.R and R/brinda_adjustment.R column maps",

      microdata_file    = cc$data_path %||% UNKNOWN,
      generated_by      = "scripts/build_assay_lineage.R",
      stringsAsFactors  = FALSE
    )
  }
}

lin <- dplyr::bind_rows(rows)
lin$microdata_file <- sub(".*[/\\\\]data[/\\\\]", "data/", lin$microdata_file)

dir.create(here("metadata"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(lin, here("metadata", "assay_lineage.csv"))

cat(sprintf("assay lineage: %d rows across %d surveys\n",
            nrow(lin), dplyr::n_distinct(lin$survey)))
unk <- vapply(lin[, c("specimen", "instrument_or_kit", "laboratory")],
              function(x) sum(x == UNKNOWN), integer(1))
cat("UNKNOWN cells:\n"); for (k in names(unk)) cat(sprintf("  %-18s %d of %d\n", k, unk[[k]], nrow(lin)))
cat("\nwrote metadata/assay_lineage.csv\n")
