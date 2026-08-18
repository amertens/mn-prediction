# =============================================================================
# scripts/check_covariate_hygiene.R
#
# Diagnostic ONLY -- reports covariate-hygiene problems in the GEE Admin-2
# predictor set. Changes nothing and depends on no toggle.
#
# Reports, per country:
#   * how many columns are cross-band _annual_* summaries
#   * which raster families emit them, and whether their bands are commensurable
#     (per the declarations in R/gee_band_semantics.R)
#   * exact-duplicate column groups (a value-based check, run here as a
#     diagnostic only -- the pipeline's pruning is name-based on purpose, see
#     R/gee_band_semantics.R)
#   * any family emitting summaries that is not yet classified
#
# Run (reads the {targets} store, so gee_admin2_* must exist):
#   Rscript scripts/check_covariate_hygiene.R
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("R", "config.R"))
source(here("R", "gee_band_semantics.R"))

countries <- names(get_country_configs())
sum_rx <- "_annual_(mean|sd|min|max|range)$"

exact_duplicate_groups <- function(df, cols) {
  num <- cols[vapply(df[cols], is.numeric, logical(1))]
  if (!length(num)) return(list())
  key <- vapply(df[num], function(x) paste0(format(round(x, 10)), collapse = "|"),
                character(1))
  tab <- table(key)
  keys <- names(tab)[tab > 1]
  lapply(keys, function(k) names(key)[key == k])
}

all_unclassified <- character(0)

for (cn in countries) {
  x <- tryCatch(targets::tar_read_raw(paste0("gee_admin2_", tolower(cn))),
                error = function(e) NULL)
  if (is.null(x)) { message(sprintf("%-14s -- not in store", cn)); next }
  v <- setdiff(names(x), c("Admin1", "Admin2"))
  if (!length(v)) { message(sprintf("%-14s -- no gee_ columns", cn)); next }

  summaries <- grep(sum_rx, v, value = TRUE)
  bad       <- gee_noncommensurable_summaries(v)
  yeardups  <- gee_static_year_duplicates(v)
  dupgroups <- exact_duplicate_groups(x, v)
  n_dupcols <- sum(lengths(dupgroups)) - length(dupgroups)
  unclass   <- gee_unclassified_families(v)
  all_unclassified <- union(all_unclassified, unclass)

  cat(sprintf("\n=== %s ===\n", cn))
  cat(sprintf("  gee_ columns                       %5d\n", length(v)))
  cat(sprintf("  cross-band _annual_* summaries     %5d  (%.0f%%)\n",
              length(summaries), 100 * length(summaries) / length(v)))
  cat(sprintf("    of which non-commensurable       %5d\n", length(bad)))
  cat(sprintf("  static-layer per-year duplicates   %5d\n", length(yeardups)))
  cat(sprintf("  columns that exactly copy another  %5d  (%d groups)\n",
              n_dupcols, length(dupgroups)))
  cat(sprintf("  would be dropped by hygiene        %5d\n",
              length(unique(c(bad, yeardups)))))
  if (length(unclass))
    cat(sprintf("  UNCLASSIFIED families (kept)       %s\n",
                paste(unclass, collapse = ", ")))

  if (length(dupgroups)) {
    cat("  example duplicate groups:\n")
    for (g in head(dupgroups, 3)) cat("    {", paste(g, collapse = ", "), "}\n")
  }
}

cat("\n=== declared band semantics ===\n")
for (pat in names(GEE_BAND_SEMANTICS))
  cat(sprintf("  %-34s %s\n", pat, GEE_BAND_SEMANTICS[[pat]]))

if (length(all_unclassified)) {
  cat("\nFamilies emitting cross-band summaries with NO declaration.\n")
  cat("They keep their summaries; classify them in R/gee_band_semantics.R:\n")
  for (f in sort(all_unclassified)) cat("  ", f, "\n")
} else {
  cat("\nEvery family emitting cross-band summaries is classified.\n")
}
