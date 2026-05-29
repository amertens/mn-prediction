# =============================================================================
# sandbox/18_continuous_biomarker_audit.R
#
# Two questions:
#   1. Is the LOCO failure for women_vitA / women_b12 / women_folate a
#      RARE-OUTCOME problem (truly low prevalence everywhere) or a
#      BETWEEN-COUNTRY THRESHOLD INCONSISTENCY problem (different
#      countries using different cutoffs)?
#   2. Do the underlying CONTINUOUS biomarker concentrations transport
#      better than the binary thresholded prevalence? If yes, modelling
#      the continuous outcome could be the fix.
#
# Audits: for each outcome, look at admin-2-level mean biomarker
# concentration across countries, the within-country distribution, and
# the gap between country means.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))

# ── Look for continuous biomarker columns in merged_<country> ──────────────
countries_lower <- c("gambia","ghana","sierraleone","malawi")
country_labels  <- c(Gambia = "gambia", Ghana = "ghana",
                      SierraLeone = "sierraleone", Malawi = "malawi")
merged_list <- setNames(lapply(countries_lower, function(c)
  tryCatch(targets::tar_read_raw(paste0("merged_", c), store = STORE),
            error = function(e) NULL)),
  names(country_labels))

cat("===== merged_ frames =====\n")
for (n in names(merged_list)) {
  m <- merged_list[[n]]
  if (is.null(m)) { cat(sprintf("  %s: NULL\n", n)); next }
  cat(sprintf("  %s: %d rows x %d cols\n", n, nrow(m), ncol(m)))
}

# ── Search for biomarker concentration columns ─────────────────────────────
biomarker_patterns <- list(
  vitA   = c("rbp", "retinol", "vita", "vit_a"),
  iron   = c("ferritin", "hemoglobin", "haemoglobin", "stfr", "transferrin"),
  zinc   = c("zinc"),
  b12    = c("b12", "cobalamin", "_b12_"),
  folate = c("folate", "rbcf", "_fol_")
)

cat("\n===== candidate biomarker columns per country =====\n")
biomarker_cols <- list()
for (n in names(merged_list)) {
  m <- merged_list[[n]]
  if (is.null(m)) next
  cat(sprintf("\n--- %s ---\n", n))
  per_biomarker <- list()
  for (bm in names(biomarker_patterns)) {
    pats <- biomarker_patterns[[bm]]
    hits <- character()
    for (p in pats) {
      hits <- c(hits, grep(p, colnames(m), ignore.case = TRUE, value = TRUE))
    }
    hits <- unique(hits)
    # Filter to numeric-ish columns: drop ones that are clearly indicator vars
    # (gw_*_<digit> style).
    hits <- hits[!grepl("^gw_.*_\\d+$", hits)]
    if (length(hits) > 0) {
      cat(sprintf("  %s (%d): %s\n", bm, length(hits),
                  paste(head(hits, 6), collapse=", ")))
      per_biomarker[[bm]] <- hits
    }
  }
  biomarker_cols[[n]] <- per_biomarker
}

# ── For each biomarker, identify a CONSISTENT column across countries ─────
cat("\n===== common biomarker columns across all 4 countries =====\n")
common_per_bm <- list()
for (bm in names(biomarker_patterns)) {
  per_country <- lapply(biomarker_cols, function(x) x[[bm]])
  common <- Reduce(intersect, per_country)
  if (length(common) > 0) {
    cat(sprintf("\n%s — common across countries (%d):\n", bm, length(common)))
    print(head(common, 10))
    common_per_bm[[bm]] <- common
  } else {
    cat(sprintf("\n%s — NO common columns across all 4 countries\n", bm))
  }
}

# ── For the most likely "primary" column per biomarker, compute admin-2 ────
#    mean and audit cross-country distribution.

# Hand-pick candidates after seeing what's there. We'll check after running
# the audit; the priority columns are typically things like:
#   gw_ferritin_ngml, gw_rbp_umol, gw_zinc_ugdl, gw_b12_pgml, gw_folate_ngml
#   (column names may vary; we'll print the actual ones found).
cat("\n===== ALL columns containing biomarker-like substrings =====\n")
for (n in names(merged_list)) {
  m <- merged_list[[n]]
  if (is.null(m)) next
  bio_search <- c("ferritin", "retinol", "rbp", "haemoglobin", "hemoglobin",
                   "_zinc", "_b12", "folate", "cobalamin", "stfr",
                   "transferrin", "hgb", "_hb_", "biomarker", "concentration",
                   "umol", "ngml", "ngdl", "ugdl", "pgml")
  hits <- character()
  for (p in bio_search) {
    hits <- c(hits, grep(p, colnames(m), ignore.case = TRUE, value = TRUE))
  }
  hits <- unique(hits)
  if (length(hits) > 0) {
    cat(sprintf("\n%s: %d biomarker-like columns\n", n, length(hits)))
    print(head(sort(hits), 30))
  }
}

cat("\nDONE — see output above to identify which continuous biomarker\n")
cat("columns can be used as transport targets for the rare-outcome\n")
cat("question. The follow-up script will build admin-2 means and run LOCO.\n")
