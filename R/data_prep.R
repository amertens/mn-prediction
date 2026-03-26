# =============================================================================
# R/data_prep.R
#
# Functions for loading and constructing analysis-ready datasets.
# Wrapped as pure functions so targets can track inputs/outputs.
# =============================================================================

#' Load the merged country dataset from disk and harmonize columns
#'
#' Handles cross-country differences:
#' - Derives gw_child_flag if missing (Ghana) from child biomarker presence
#' - Derives binary VAD columns if missing (Sierra Leone) from continuous RBP
#' - Malawi uses raw biomarker names (rbp, fer, vitA_def) and a text
#'   `population` column — no harmonization needed, just load.
#'
#' @param data_path Path to the .rds file
#' @return data.frame with harmonized columns
load_merged_data <- function(data_path) {
  stopifnot(file.exists(data_path))
  d <- readRDS(data_path)
  cat(sprintf("[data_prep] Loaded %s: %d rows x %d cols\n",
              basename(data_path), nrow(d), ncol(d)))

  # ── Detect Malawi by its column naming pattern ─────────────────────────
  # Malawi has a text `population` column and no gw_-prefixed biomarkers.
  # The config handles population filtering via child_flag_val text matching,
  # so no harmonization is needed here.
  is_malawi <- "population" %in% colnames(d) &&
    is.character(d[["population"]]) &&
    any(d[["population"]] == "preschool children", na.rm = TRUE)

  if (is_malawi) {
    cat("  Malawi MNS data detected — using raw biomarker columns\n")
    return(d)
  }

  # ── Derive gw_child_flag if missing ────────────────────────────────────
  # Ghana doesn't have this column. Derive from child-specific biomarkers:
  # rows with non-NA child RBP (gw_cRBPAdjThurn or gw_cRBPAdj) are children.
  if (!"gw_child_flag" %in% colnames(d)) {
    child_rbp <- NULL
    for (col in c("gw_cRBPAdjThurn", "gw_cRBPAdjBrinda", "gw_cRBPAdj",
                  "gw_cRBP", "gw_childid")) {
      if (col %in% colnames(d) && any(!is.na(d[[col]]))) {
        child_rbp <- col
        break
      }
    }
    if (!is.null(child_rbp)) {
      d$gw_child_flag <- ifelse(!is.na(d[[child_rbp]]), 1L, 0L)
      cat(sprintf("  Derived gw_child_flag from %s: %d children, %d women\n",
                  child_rbp, sum(d$gw_child_flag == 1), sum(d$gw_child_flag == 0)))
    } else {
      warning("Cannot derive gw_child_flag — no child biomarker column found")
    }
  }

  # ── Derive binary VAD columns if missing ────────────────────────────────
  # Sierra Leone has continuous RBP but no pre-computed binary VAD column.
  # Derive as RBP < 0.70 µmol/L (standard WHO threshold).
  derive_binary <- function(d, cont_col, bin_col, threshold = 0.70) {
    if (cont_col %in% colnames(d) && !(bin_col %in% colnames(d))) {
      d[[bin_col]] <- ifelse(d[[cont_col]] < threshold, 1L, 0L)
      d[[bin_col]][is.na(d[[cont_col]])] <- NA_integer_
      n_def <- sum(d[[bin_col]] == 1, na.rm = TRUE)
      n_tot <- sum(!is.na(d[[bin_col]]))
      cat(sprintf("  Derived %s from %s < %.2f: %d/%d deficient\n",
                  bin_col, cont_col, threshold, n_def, n_tot))
    }
    d
  }

  d <- derive_binary(d, "gw_cRBPAdj", "gw_cVAD", 0.70)
  d <- derive_binary(d, "gw_wRBPAdj", "gw_wVAD", 0.70)

  d
}


#' Build an outcome-specific dataset (one population x one micronutrient)
#'
#' Filters to the correct population, selects predictors, removes leakage
#' columns, and drops rows with missing outcome.
#'
#' @param merged_data The full merged dataset
#' @param cc Country config (from get_country_configs())
#' @param oc Outcome config (one element of cc$outcomes)
#' @return list with components: data, Xvars_full, Xvars, domain_vars
build_outcome_dataset <- function(merged_data, cc, oc) {

  d <- merged_data

  # Filter to population
  pop_col <- cc$child_flag
  if (!is.null(pop_col) && pop_col %in% colnames(d)) {
    d <- d[d[[pop_col]] == oc$child_flag_val, ]
  }

  # Require non-missing outcome
  if (!is.null(oc$continuous) && oc$continuous %in% colnames(d))
    d <- d[!is.na(d[[oc$continuous]]), ]
  if (!is.null(oc$binary) && oc$binary %in% colnames(d))
    d <- d[!is.na(d[[oc$binary]]), ]

  # Recode binary outcome if coded as factor-like 1/2 instead of 0/1.
  # Some surveys (e.g., Sierra Leone MNS) code deficiency as 1=yes, 2=no.
  # Convert to standard 0/1 where 1=deficient, 0=not deficient.
  if (!is.null(oc$binary) && oc$binary %in% colnames(d)) {
    yb <- d[[oc$binary]]
    yb_vals <- sort(unique(as.numeric(yb[!is.na(yb)])))
    if (length(yb_vals) == 2 && all(yb_vals == c(1, 2))) {
      cat(sprintf("  Recoding %s from {1,2} to {1,0} (1=deficient, 2=not → 0)\n",
                  oc$binary))
      d[[oc$binary]] <- ifelse(d[[oc$binary]] == 1, 1, 0)
    }
  }

  cat(sprintf("[build_outcome] %s — %s: %d observations\n",
              cc$country, oc$tag, nrow(d)))

  # Identify predictor columns by domain prefix
  all_cols <- colnames(d)
  domain_vars <- list()
  Xvars_all <- character()


  for (nm in names(cc$domains)) {
    dom <- cc$domains[[nm]]
    prefix_cols <- all_cols[startsWith(all_cols, dom$prefix)]
    if (!is.null(dom$extra)) prefix_cols <- c(prefix_cols, intersect(dom$extra, all_cols))

    # Remove leakage columns from GW domain
    if (nm == "GW" && !is.null(cc$gw_exclude_patterns)) {
      leak <- Reduce("|", lapply(cc$gw_exclude_patterns,
                                  function(p) grepl(p, prefix_cols, ignore.case = TRUE)))
      prefix_cols <- prefix_cols[!leak]

      # Also exclude the outcome columns themselves
      exclude_exact <- c(oc$continuous, oc$binary, cc$admin1_col, cc$admin2_col,
                         cc$cluster_id, cc$weight_col, cc$strata_col, cc$child_flag)
      prefix_cols <- setdiff(prefix_cols, exclude_exact)
    }

    domain_vars[[nm]] <- prefix_cols
    Xvars_all <- c(Xvars_all, prefix_cols)
  }

  Xvars_full <- unique(Xvars_all)
  Xvars_no_gw <- unique(setdiff(Xvars_all, domain_vars[["GW"]]))

  # Keep only columns we need
  keep_cols <- unique(c(
    oc$continuous, oc$binary,
    cc$cluster_id, cc$admin1_col, cc$admin2_col,
    cc$weight_col, cc$strata_col, cc$psu_col,
    Xvars_full
  ))
  keep_cols <- intersect(keep_cols, colnames(d))
  d <- d[, keep_cols, drop = FALSE]

  list(
    data        = d,
    Xvars_full  = Xvars_full,
    Xvars       = Xvars_no_gw,
    domain_vars = domain_vars
  )
}
