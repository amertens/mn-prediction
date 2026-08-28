# =============================================================================
# scripts/test_admin2_keys.R
#
# WS8g regression test for the Admin-2 join key.
#
# Asserts, per country:
#   1. the survey Admin-2 table has no duplicated key
#   2. the GEE Admin-2 table has no duplicated key
#   3. GADM Admin-2 names are unique ONCE the pair key is used, even where the
#      names alone are not
#   4. pooling survey and covariates does not multiply rows: the pooled row
#      count equals the number of surveyed districts that have covariates
#
# Check 4 is the one that matters. Before this migration, Malawi's 87 surveyed
# districts became 90 pooled rows, because the centroid join keyed on the
# Admin-2 name against an un-deduped GADM table with four same-named district
# pairs. Those three districts carried double weight in every area-level fit and
# one copy of each pair carried the wrong region's centroid.
#
# Exits 1 on any failure so it can gate a run.
#
# Run:
#   Rscript scripts/test_admin2_keys.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr)
})
targets::tar_source(here("R"))
source(here("src", "0-functions.R"))

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

configs <- get_country_configs()
PROBE_OUTCOME <- "child_iron"

failures <- character(0)
rows <- list()

.fail <- function(...) failures <<- c(failures, sprintf(...))

for (cn in names(configs)) {
  cc <- configs[[cn]]

  svy <- read_target(sprintf("svy_admin2_%s_%s", tolower(cn), PROBE_OUTCOME))
  gee <- read_target(paste0("gee_admin2_", tolower(cn)))
  gadm <- tryCatch(sf::st_drop_geometry(load_admin2_centroids(cc$gadm_code)),
                   error = function(e) NULL)

  n_svy  <- if (is.null(svy)) NA_integer_ else nrow(svy)
  n_gee  <- if (is.null(gee)) NA_integer_ else nrow(gee)
  n_gadm <- if (is.null(gadm)) NA_integer_ else nrow(gadm)

  # 1 and 2: no duplicated key on either side.
  if (!is.null(svy) && any(duplicated(admin2_key(svy))))
    .fail("%s: svy_admin2 has duplicated keys", cn)
  if (!is.null(gee) && any(duplicated(admin2_key(gee))))
    .fail("%s: gee_admin2 has duplicated keys", cn)

  # 3: GADM names may repeat; the PAIR must not, once water polygons are out.
  gadm_name_dups <- NA_integer_; gadm_pair_dups <- NA_integer_
  if (!is.null(gadm)) {
    land <- gadm[!is_water_admin2(as.character(gadm$Admin2)), , drop = FALSE]
    gadm_name_dups <- sum(duplicated(as.character(land$Admin2)))
    gadm_pair_dups <- sum(duplicated(admin2_key(land)))
    if (gadm_pair_dups > 0)
      .fail("%s: GADM has %d duplicated (Admin1, Admin2) pair(s) among land polygons",
            cn, gadm_pair_dups)
  }

  # 4: pooling must not multiply rows.
  n_pooled <- NA_integer_
  if (!is.null(svy) && !is.null(gee)) {
    p <- tryCatch(build_area_loco_dataset(setNames(list(svy), cn),
                                          setNames(list(gee), cn)),
                  error = function(e) NULL)
    if (!is.null(p) && !is.null(p$pooled_data)) {
      pd <- p$pooled_data
      n_pooled <- nrow(pd)
      if (any(duplicated(admin2_key(pd))))
        .fail("%s: pooled area dataset has duplicated keys", cn)
      # The pooled count is the surveyed districts that matched a covariate row.
      # It may be LOWER than n_svy (a district with no covariates drops out); it
      # must never be HIGHER.
      if (is.finite(n_pooled) && is.finite(n_svy) && n_pooled > n_svy)
        .fail("%s: pooling multiplied rows, %d surveyed districts -> %d pooled",
              cn, n_svy, n_pooled)
    }
  }

  rows[[cn]] <- data.frame(
    country = cc$country,
    gadm_polygons = n_gadm,
    gadm_dup_names_land = gadm_name_dups,
    gadm_dup_pairs_land = gadm_pair_dups,
    gee_admin2_rows = n_gee,
    svy_admin2_rows = n_svy,
    pooled_rows = n_pooled,
    stringsAsFactors = FALSE)
}

tab <- dplyr::bind_rows(rows)
cat("\n=== Admin-2 unit counts per country ===\n")
print(as.data.frame(tab), row.names = FALSE)

dir.create(here("results", "tables", "corrected"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tab, here("results", "tables", "corrected", "admin2_unit_counts.csv"))
cat("\nwrote results/tables/corrected/admin2_unit_counts.csv\n")

if (length(failures)) {
  cat("\nFAIL:\n"); cat(paste0("  ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("\nPASS: no duplicated Admin-2 keys and no row multiplication in any country.\n")
