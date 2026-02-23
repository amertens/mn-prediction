# =============================================================================
# src/analysis/01_load_and_construct.R
#
# Loads the Gambia merged dataset, drops geometry, constructs the predictor
# variable lists, removes outcome-leaking gw_ variables, and builds the four
# analysis datasets (child/women × VitA/iron).
#
# Inputs  : cfg (from src/analysis/config.R, must be in environment)
# Outputs : (assigned into the calling environment via list2env)
#   gw_data_list  - named list with df_child_VitA, df_mom_VitA,
#                   df_child_iron, df_mom_iron
#   Xvars_full    - full predictor variable vector (incl. gw_ non-outcome vars)
#   Xvars         - predictor vector without gw_ domain vars (external data only)
#   domain_vars   - named list of domain-specific predictor columns
# =============================================================================

cat("\n[01] Loading and constructing analysis datasets...\n")

# ---- 1. Load ----------------------------------------------------------------
df <- readRDS(here::here(cfg$data_path))
df <- sf::st_drop_geometry(df)

# Unique row identifier (matches existing script convention)
df$dataid <- paste0("gambia", seq_len(nrow(df)))

# ---- 2. Identify predictor columns by domain prefix -------------------------
domain_vars <- lapply(cfg$domains, function(dom) {
  vars <- colnames(df)[grepl(dom$prefix, colnames(df), fixed = TRUE)]
  if (!is.null(dom$extra)) {
    extra_present <- dom$extra[dom$extra %in% colnames(df)]
    vars <- unique(c(vars, extra_present))
  }
  vars
})
names(domain_vars) <- names(cfg$domains)

# ---- 3. Remove outcome-leaking gw_ variables --------------------------------
# Exclude any gw_ column matching the leakage patterns in cfg$gw_exclude_patterns
# This prevents direct biomarker measures from appearing as predictors.
leakage_pattern <- paste(cfg$gw_exclude_patterns, collapse = "|")
domain_vars[["GW"]] <- domain_vars[["GW"]][
  !grepl(leakage_pattern, domain_vars[["GW"]])
]

cat(sprintf(
  "  gw_ predictors after leakage removal: %d\n",
  length(domain_vars[["GW"]])
))

# ---- 4. Build Xvars vectors -------------------------------------------------
# Xvars_full : identifiers + all domain predictors (incl. cleaned gw_ vars)
# Xvars      : identifiers + external data only (no gw_ predictors)
# Both include dataid, Admin1, gw_month as meta columns.

meta_cols   <- c("dataid", cfg$admin1, cfg$month)

Xvars_full <- unique(c(
  meta_cols,
  domain_vars[["GW"]],
  domain_vars[["DHS"]],
  domain_vars[["MICS"]],
  domain_vars[["IHME"]],
  domain_vars[["LSMS"]],
  domain_vars[["MAP"]],
  domain_vars[["WFP"]],
  domain_vars[["FLUNET"]],
  domain_vars[["GEE"]]
))

Xvars <- unique(c(
  meta_cols,
  domain_vars[["DHS"]],
  domain_vars[["MICS"]],
  domain_vars[["IHME"]],
  domain_vars[["LSMS"]],
  domain_vars[["MAP"]],
  domain_vars[["WFP"]],
  domain_vars[["FLUNET"]],
  domain_vars[["GEE"]]
))

# Keep only columns that actually exist in df
Xvars_full <- Xvars_full[Xvars_full %in% colnames(df)]
Xvars      <- Xvars[Xvars %in% colnames(df)]

cat(sprintf("  Xvars_full length: %d\n", length(Xvars_full)))
cat(sprintf("  Xvars (no gw_)  : %d\n", length(Xvars)))

# ---- 5. Build per-outcome analysis datasets ---------------------------------
# Each dataset includes:
#   - The continuous outcome column
#   - The binary deficiency column (if present in df)
#   - gw_cnum (cluster id)
#   - All Xvars_full columns
# Rows are filtered to the correct population and to non-missing continuous outcome.

build_dataset <- function(outcome_cfg, df, Xvars_full, cluster_id) {
  pop_flag   <- outcome_cfg$child_flag_val
  cont_col   <- outcome_cfg$continuous
  bin_col    <- outcome_cfg$binary
  flag_col   <- cfg$child_flag

  # Columns to always keep in the dataset
  keep_cols <- unique(c(cont_col, cluster_id, Xvars_full))

  # Also keep binary outcome column if it exists (for reference / direct binary modeling)
  if (!is.null(bin_col) && bin_col %in% colnames(df)) {
    keep_cols <- unique(c(keep_cols, bin_col))
  }
  # Keep existing cols only
  keep_cols <- keep_cols[keep_cols %in% colnames(df)]

  d <- df %>%
    dplyr::filter(.data[[flag_col]] == pop_flag) %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(.data[[cont_col]]))

  cat(sprintf(
    "  %-25s : %d rows, %d cols\n",
    outcome_cfg$tag, nrow(d), ncol(d)
  ))
  d
}

cat("  Building per-outcome datasets:\n")
gw_data_list <- list(
  child_vitA = build_dataset(cfg$outcomes$child_vitA, df, Xvars_full, cfg$cluster_id),
  women_vitA = build_dataset(cfg$outcomes$women_vitA, df, Xvars_full, cfg$cluster_id),
  child_iron = build_dataset(cfg$outcomes$child_iron, df, Xvars_full, cfg$cluster_id),
  women_iron = build_dataset(cfg$outcomes$women_iron, df, Xvars_full, cfg$cluster_id)
)

# Summary
cat(sprintf("\n  Vitamin A (children): n = %d\n", nrow(gw_data_list$child_vitA)))
cat(sprintf("  Vitamin A (women)   : n = %d\n", nrow(gw_data_list$women_vitA)))
cat(sprintf("  Iron (children)     : n = %d\n", nrow(gw_data_list$child_iron)))
cat(sprintf("  Iron (women)        : n = %d\n", nrow(gw_data_list$women_iron)))

cat("[01] Done.\n\n")
