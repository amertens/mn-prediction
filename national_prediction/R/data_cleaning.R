# =============================================================================
# National Prediction Pipeline — Data Cleaning
# =============================================================================

#' Load national-level proxy predictor dataset
#'
#' Reads the combined DHS/IHME/etc. predictor dataset and fixes column naming.
#' @param path Path to combined_dataset.dta
#' @return data.frame of national-level predictors with country + year keys
load_predictor_data <- function(path) {
  data <- haven::read_dta(path)
  colnames(data)[1] <- "country"
  data
}


#' Clean BRINDA prevalence data and merge with predictors
#'
#' Loads BRINDA prevalence RDS, row-binds country datasets, parses prevalence
#' strings (strips CI text), and merges with national predictor data.
#' @param brinda_path Path to clean_brinda_data.rds
#' @param predictor_data data.frame from load_predictor_data()
#' @return Merged data.frame ready for modeling
clean_brinda_data <- function(brinda_path, predictor_data) {
  brinda <- readRDS(brinda_path)
  brinda_df <- data.table::rbindlist(brinda, fill = TRUE)

  # Parse prevalence columns: extract numeric value before " ("
  prev_cols <- c("prev_ferr", "prev_stfr", "prev_rol", "prev_rbp", "prev_zinc")
  for (col in prev_cols) {
    if (col %in% colnames(brinda_df)) {
      brinda_df[[col]] <- as.numeric(
        stringr::str_split(brinda_df[[col]], " \\(", simplify = TRUE)[, 1]
      )
    }
  }

  # Merge with predictors on country + year
  d <- dplyr::left_join(brinda_df, predictor_data, by = c("country", "year"))
  as.data.frame(d)
}


#' Clean VMNIS indicator data and merge with predictors
#'
#' Loads WHO VMNIS long-format data, renames key columns, selects relevant
#' fields, and merges with national predictor data.
#' @param vmnis_path Path to VMNISIndicator_long_format.dta
#' @param predictor_data data.frame from load_predictor_data()
#' @return Merged data.frame ready for indicator-level filtering
clean_vmnis_data <- function(vmnis_path, predictor_data) {
  vmnis_full <- haven::read_dta(vmnis_path)

  vmnis <- vmnis_full %>%
    dplyr::rename(country = Country, year = Beginyear) %>%
    dplyr::select(
      country, year, Indicator, Representativeness, Population,
      Mean, Geomean, Median,
      dplyr::starts_with("Prevalence"),
      Dataadjustedfor, Surveymethodology, Methodofanalysis,
      Season, Agefrom, Ageto
    )

  # Merge with predictors on country + year
  d <- dplyr::left_join(vmnis, predictor_data, by = c("country", "year"))
  as.data.frame(d)
}


#' Filter VMNIS data to a specific indicator
#'
#' @param vmnis_merged Merged VMNIS data from clean_vmnis_data()
#' @param indicator VMNIS indicator string (e.g., "Zinc (plasma or serum)")
#' @return Filtered data.frame
filter_vmnis_indicator <- function(vmnis_merged, indicator) {
  vmnis_merged %>%
    dplyr::filter(Indicator == indicator)
}


#' Build predictor variable list for BRINDA models
#'
#' @param d Merged BRINDA data
#' @return Character vector of predictor column names
get_brinda_xvars <- function(d) {
  # Predictors start after the BRINDA-specific columns (first 8 cols) + year + age
  predictor_start <- 9
  c("year", "age", colnames(d)[predictor_start:ncol(d)])
}


#' Build predictor variable list for VMNIS models
#'
#' @param d Merged VMNIS data
#' @return Character vector of predictor column names
get_vmnis_xvars <- function(d) {
  # VMNIS metadata columns end around column 31; predictors start after
  # Include age range and metadata as potential predictors
  predictor_start <- which(colnames(d) == "Ageto") + 1
  # Find where the predictor columns from combined_dataset begin
  # These come after the VMNIS-specific columns
  predictor_cols <- colnames(d)[predictor_start:ncol(d)]
  c("year", "Agefrom", "Ageto", "Dataadjustedfor", "Season", predictor_cols)
}


#' Prepare VMNIS data for heatmap visualization
#'
#' Filters to national, both urban/rural surveys in LMICs, collapses population
#' groups and indicator categories, and summarizes the latest survey per triad.
#' @param vmnis_path Path to VMNISIndicator_long_format.dta
#' @param lmic_countries Character vector of LMIC country names
#' @return data.frame with Country, Population, mn_group, mean, n_surveys
clean_vmnis_for_heatmap <- function(vmnis_path, lmic_countries) {
  d <- haven::read_dta(vmnis_path)

  # Map indicators to broad micronutrient groups
  mn_map <- function(x) {
    dplyr::case_when(
      x %in% c("Ferritin", "Serum transferrin receptor",
               "Iron (serum or plasma)", "Total iron binding capacity (TIBC)",
               "Transferrin saturation", "Mean corpuscular volume (MCV)") ~ "Iron",
      x %in% c("Retinol (plasma or serum)", "Beta-carotene")              ~ "Vitamin A",
      x == "Zinc (plasma or serum)"                                       ~ "Zinc",
      x %in% c("Goitre by ultrasound", "Iodine")                          ~ "Iodine",
      x %in% c("Folate (plasma or serum)", "Folate (red blood cell)")     ~ "Folate",
      x == "Vitamin B12"                                                  ~ "Vitamin B12",
      x == "25-Hydroxyvitamin D"                                          ~ "Vitamin D",
      x == "Vitamin E (alpha-tocopherol)"                                 ~ "Vitamin E",
      x == "Vitamin C"                                                    ~ "Vitamin C",
      x == "Haemoglobin"                                                  ~ "Haemoglobin",
      TRUE                                                                ~ NA_character_
    )
  }

  plotdf <- d %>%
    dplyr::filter(
      Representativeness == "national",
      Areacovered == "both urban and rural",
      Country %in% lmic_countries
    ) %>%
    dplyr::mutate(
      Population = dplyr::case_when(
        Population %in% c(
          "Infants (1-11 months)", "Children (12-23 months)",
          "Neonates (0-28 days)", "Young infants (0-5 months)",
          "Children (24-59 months)"
        ) ~ "Children <5",
        TRUE ~ Population
      ),
      mn_group = mn_map(Indicator)
    ) %>%
    tidyr::drop_na(mn_group)

  # Keep most recent survey per Country x Population x micronutrient
  survey_summary <- plotdf %>%
    dplyr::group_by(Indicator) %>%
    dplyr::mutate(Mean = scale(Mean)) %>%
    dplyr::group_by(Country, Population, mn_group) %>%
    dplyr::filter(Date == max(Date)) %>%
    dplyr::summarise(
      n_surveys = dplyr::n(),
      mean      = mean(Mean, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::filter(!is.na(mean))

  survey_summary
}
