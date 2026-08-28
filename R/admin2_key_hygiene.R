# =============================================================================
# R/admin2_key_hygiene.R
#
# Every join in this pipeline is keyed on the Admin-2 NAME. Two things break
# that, and both reach the published maps.
#
# 1. GADM SHIPS INLAND WATER BODIES AS ADMIN-2 POLYGONS.
#    Malawi has 4, Tanzania 5. They are not populated prediction units, and the
#    zonal statistics taken over them are NDVI, soil chemistry and land surface
#    temperature of open water.
#
#    They are also not confined to the covariate side. `Lake Malawi` appears in
#    the Malawi survey tables with n_svy 21-33, and `Lake Victoria` appears in
#    the Tanzania tables with n_svy 60-110 -- one of that country's larger
#    "districts". Respondents have been geocoded onto a lake. Those rows were
#    being used to FIT the area-level models and were being mapped.
#
# 2. ADMIN-2 NAMES REPEAT.
#    Partly because the lake polygons are multi-part, partly for real: Malawi
#    has TA Lundu, TA Malemia, TA Ngabu and TA Pemba in more than one region.
#    An unguarded name join then multiplies rows. extract_gee_admin2() already
#    collapsed duplicates, but only on the local-raster path -- a country
#    without a raster directory returns early through the legacy-parity CSV and
#    skips the dedup entirely, which is why gee_admin2_tanzania carried 186 rows
#    for 176 districts and its inner_join with the survey produced 167 rows for
#    163 surveyed districts.
#
# These helpers are applied at every point where an Admin-2 keyed table is
# produced, so the guarantee holds regardless of which code path a country
# takes.
# =============================================================================

# =============================================================================
# 2026-08-28: the join key is now Admin1 + Admin2, not Admin2 alone.
#
# WHAT WAS WRONG. Malawi's GADM Admin-2 layer has 256 polygons under 243
# distinct names but 256 distinct (Admin1, Admin2) pairs: 13 water polygons plus
# four genuine same-name districts in different regions (TA Lundu in Blantyre and
# Chikwawa, TA Ngabu in Chikwawa and Nsanje, TA Pemba in Dedza and Salima, TA
# Malemia in Nsanje and Zomba). Three of the four are surveyed.
#
# Two distinct failures followed, both measured before this change:
#
#   a) dedupe_admin2_key() AVERAGED each same-name pair into one covariate row,
#      so respondents in TA Lundu, Chikwawa were joined to covariates that were
#      half from TA Lundu, Blantyre.
#   b) build_area_loco_dataset()'s centroid left_join(by = "Admin2") MULTIPLIED
#      rows against the un-deduped GADM table. Malawi's 87 surveyed districts
#      became 90 pooled rows. Those three districts carried double weight in
#      every area-level fit, and one copy of each pair carried the wrong region's
#      centroid, which corrupts the spatial comparators specifically.
#
# WHY A PAIR KEY IS SAFE HERE. Every one of Malawi's 89 surveyed
# (Admin1, Admin2) pairs matches a GADM pair exactly, and Gambia, Ghana and
# Sierra Leone have no name collisions at all (37, 260 and 14 polygons, all
# uniquely named). So the pair key changes nothing for three countries and
# separates three districts in the fourth.
#
# FALLBACK. Every helper below degrades to the name-only behaviour when Admin1
# is missing on either side, so a table that never carried Admin1 is unaffected.
# =============================================================================

#' Names GADM uses for inland water polygons.
#' Deliberately anchored so a real district called e.g. "Lakeside" is not caught.
ADMIN2_WATER_PATTERN <- "^Lake |^Water|^Lac |^Sea$|Reservoir"

#' The canonical Admin-2 join key.
#'
#' Admin1 + Admin2 where Admin1 is available, Admin2 alone otherwise. The
#' separator is chosen not to occur in an administrative name.
#'
#' @param d data.frame with Admin2 and optionally Admin1
#' @return character vector of keys
admin2_key <- function(d) {
  if (is.null(d) || !"Admin2" %in% names(d)) return(character(0))
  a2 <- as.character(d$Admin2)
  if (!"Admin1" %in% names(d)) return(a2)
  a1 <- as.character(d$Admin1)
  ifelse(is.na(a1) | !nzchar(a1), a2, paste(a1, a2, sep = "||"))
}

#' Can two tables be joined on the pair key?
#'
#' TRUE only when BOTH carry a usable Admin1. Used to decide per join rather
#' than globally, so a partially migrated pipeline stays correct.
can_pair_join <- function(x, y) {
  ok <- function(d) !is.null(d) && "Admin1" %in% names(d) &&
    any(!is.na(d$Admin1) & nzchar(as.character(d$Admin1)))
  ok(x) && ok(y)
}

#' The join-by vector for an Admin-2 join: the pair when both sides support it.
admin2_join_by <- function(x, y) if (can_pair_join(x, y)) c("Admin1", "Admin2") else "Admin2"

#' Rows whose Admin-2 name denotes a water body rather than a district.
is_water_admin2 <- function(x) grepl(ADMIN2_WATER_PATTERN, x, ignore.case = TRUE)

#' Drop water-body rows from any Admin-2 keyed table.
#'
#' @param d data.frame with an Admin2 column
#' @param what label for the log line
drop_water_admin2 <- function(d, what = "table") {
  if (is.null(d) || !"Admin2" %in% names(d)) return(d)
  w <- is_water_admin2(as.character(d$Admin2))
  if (any(w, na.rm = TRUE)) {
    cat(sprintf("[admin2_hygiene] %s: dropped %d water-body row(s): %s\n",
                what, sum(w, na.rm = TRUE),
                paste(unique(d$Admin2[w]), collapse = ", ")))
    d <- d[!w, , drop = FALSE]
  }
  d
}

#' Collapse repeated Admin-2 names to one row.
#'
#' Numeric columns are averaged, everything else takes the first value. This
#' matches what the survey side already does -- svy_prev is computed per Admin-2
#' NAME -- so collapsing keeps the two sides on the same key.
#'
#' @param d data.frame with an Admin2 column
#' @param what label for the log line
dedupe_admin2_key <- function(d, what = "table") {
  if (is.null(d) || !"Admin2" %in% names(d)) return(d)

  # Collapse on the PAIR when Admin1 is available. Two same-named districts in
  # different regions are different districts and averaging them was the bug
  # this migration fixes; what remains to collapse is genuine multi-part
  # geometry, where the same district appears more than once within one region.
  key <- admin2_key(d)
  if (!any(duplicated(key))) return(d)

  by_pair <- "Admin1" %in% names(d)
  dupk <- unique(key[duplicated(key)])
  cat(sprintf("[admin2_hygiene] %s: collapsing %d duplicated %s: %s\n",
              what, length(dupk),
              if (by_pair) "(Admin1, Admin2) key(s)" else "Admin-2 name(s)",
              paste(utils::head(dupk, 6), collapse = ", ")))

  d$.a2key <- key
  grp <- ".a2key"
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  oth <- setdiff(names(d), c(grp, num))
  out <- d |>
    dplyr::group_by(.data[[grp]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(oth), dplyr::first),
      dplyr::across(dplyr::all_of(num), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop") |>
    as.data.frame()
  out[[grp]] <- NULL
  d$.a2key <- NULL
  out[, intersect(names(d), names(out)), drop = FALSE]
}

#' Both fixes at once: the standard treatment for any Admin-2 keyed table.
clean_admin2_keys <- function(d, what = "table") {
  dedupe_admin2_key(drop_water_admin2(d, what), what)
}

#' Reassign survey rows that landed in a lake to the nearest LAND district.
#'
#' WHY SNAP RATHER THAN DROP. The respondents assigned to a water polygon are
#' not people living on a lake and they are not a coding error -- they are real
#' lakeside communities whose cluster GPS was displaced into the water by the
#' DHS/MICS confidentiality procedure (up to 2 km urban, 5 km rural, 10 km for a
#' random 1% of rural clusters). Measured distances to the nearest land district:
#'
#'   Malawi   0.00 km -> TA Musisya      (shoreline)
#'   Malawi   3.84 km -> TA Timbiri      (within the 5 km rural limit)
#'   Malawi   7.95 km -> TA Timbiri      (within the 10 km 1%-of-rural limit)
#'   Tanzania 0.36 km -> Chato           (shoreline)
#'
#' Every one falls inside the documented displacement bounds, so dropping them
#' discards real observations -- 87 individuals in Malawi and 175 in Tanzania,
#' the latter one of that survey's larger clusters. Snapping to the nearest land
#' polygon is what the displacement procedure implies.
#'
#' Water polygons are still DROPPED from the covariate/prediction side by
#' drop_water_admin2(): a lake is not a populated prediction unit, and its zonal
#' statistics are NDVI and soil chemistry of open water.
#'
#' @param d data.frame with Admin2 and lon/lat columns (individual or cluster level)
#' @param gadm_code ISO3 code for load_admin2_centroids()
#' @param max_km refuse to move a point further than this; beyond the 10 km
#'   displacement limit the assignment is not explained by displacement and the
#'   row is left alone (and so will be dropped downstream) rather than moved to
#'   a district it may have nothing to do with
#' @return `d` with Admin2 (and Admin1, if present) reassigned for water rows
snap_water_to_land <- function(d, gadm_code, max_km = 12, what = "table") {
  if (is.null(d) || !all(c("Admin2", "lon", "lat") %in% names(d))) return(d)
  w <- is_water_admin2(as.character(d$Admin2)) &
    is.finite(suppressWarnings(as.numeric(d$lon)))
  if (!any(w)) return(d)
  if (!requireNamespace("sf", quietly = TRUE)) {
    warning("sf unavailable; leaving water-assigned rows in place")
    return(d)
  }
  polys <- tryCatch(load_admin2_centroids(gadm_code), error = function(e) NULL)
  if (is.null(polys)) {
    warning(sprintf("[admin2_hygiene] %s: GADM unavailable, cannot snap", what))
    return(d)
  }
  land <- polys[!is_water_admin2(as.character(polys$Admin2)), , drop = FALSE]
  if (!nrow(land)) return(d)

  # one lookup per distinct coordinate, not per row
  key <- paste(round(as.numeric(d$lon[w]), 6), round(as.numeric(d$lat[w]), 6))
  uk <- !duplicated(key)
  pts <- sf::st_as_sf(
    data.frame(lon = as.numeric(d$lon[w])[uk], lat = as.numeric(d$lat[w])[uk]),
    coords = c("lon", "lat"), crs = 4326)
  dm <- sf::st_distance(pts, sf::st_transform(land, 4326))
  near <- apply(dm, 1, which.min)
  km <- apply(dm, 1, min) / 1000
  ok <- km <= max_km

  map_a2 <- setNames(as.character(land$Admin2)[near], key[uk])
  map_a1 <- setNames(as.character(land$Admin1)[near], key[uk])
  map_ok <- setNames(ok, key[uk])

  moved <- map_ok[key]
  moved[is.na(moved)] <- FALSE
  idx <- which(w)[moved]
  if (length(idx)) {
    d$Admin2[idx] <- map_a2[key[moved]]
    if ("Admin1" %in% names(d)) d$Admin1[idx] <- map_a1[key[moved]]
  }
  cat(sprintf(paste0("[admin2_hygiene] %s: snapped %d row(s) at %d displaced ",
                     "cluster location(s) from water to the nearest land ",
                     "district (max %.2f km); %d row(s) beyond %g km left in place\n"),
              what, length(idx), sum(ok), if (any(ok)) max(km[ok]) else 0,
              sum(w) - length(idx), max_km))
  d
}

#' Guard for a join about to be keyed on Admin2.
#'
#' Duplicated keys on either side silently multiply rows, and the resulting
#' table looks perfectly normal. Fail loudly instead.
assert_unique_admin2 <- function(d, what = "table") {
  if (is.null(d) || !"Admin2" %in% names(d)) return(invisible(d))
  key <- admin2_key(d)
  if (any(duplicated(key)))
    stop(sprintf("%s has %d duplicated Admin-2 key(s) (%s) - run clean_admin2_keys() first",
                 what, sum(duplicated(key)),
                 paste(utils::head(unique(key[duplicated(key)]), 5), collapse = ", ")))
  invisible(d)
}

#' Report survey rows that a pair-keyed join dropped but a name join would keep.
#'
#' A pair join is only correct when both sides agree on Admin1. They can
#' disagree while a store is PARTIALLY rebuilt: the name-only dedupe kept
#' whichever region came first, so Malawi's stored covariate table labels
#' TA Lundu as Blantyre and TA Malemia as Nsanje, while the survey observed both
#' in Chikwawa and Zomba respectively. Once the survey side gains Admin1 and the
#' covariate side has not yet been rebuilt, those districts stop matching.
#'
#' Dropping them is the safe failure (better an absent district than one wearing
#' another region's covariates), but it must not be silent, because a quietly
#' shorter table looks entirely normal. This says so and names the districts.
#'
#' @param x the left (survey) table
#' @param y the right (covariate) table
#' @param by the join key actually used
#' @param what label for the message
#' @return invisible character vector of the affected Admin-2 names
report_pair_join_losses <- function(x, y, by, what = "join") {
  if (!identical(by, c("Admin1", "Admin2"))) return(invisible(character(0)))
  if (is.null(x) || is.null(y)) return(invisible(character(0)))
  xk <- admin2_key(x); yk <- admin2_key(y)
  lost <- !(xk %in% yk)
  if (!any(lost)) return(invisible(character(0)))
  # Of the lost rows, which would have matched on the NAME alone? Those are the
  # Admin1 disagreements, i.e. the stale-table case rather than a genuinely
  # absent district.
  nm <- as.character(x$Admin2)[lost]
  rescuable <- nm[nm %in% as.character(y$Admin2)]
  if (length(rescuable))
    warning(sprintf(paste0("[admin2_hygiene] %s: %d district(s) dropped by the ",
                           "(Admin1, Admin2) join that WOULD match on the name ",
                           "alone: %s. The two sides disagree on Admin1, which ",
                           "usually means one of them predates the join-key ",
                           "migration. Rebuild the covariate targets."),
                    what, length(rescuable),
                    paste(utils::head(rescuable, 6), collapse = ", ")), call. = FALSE)
  invisible(rescuable)
}

#' Warn, rather than stop, when a join has multiplied rows.
#'
#' Used at join sites where failing hard would take down a whole pipeline run
#' for a defect that is worth surfacing but not worth aborting on. Row
#' multiplication from a duplicated key is silent and the resulting table looks
#' entirely normal, which is how Malawi's 87 surveyed districts became 90 pooled
#' rows unnoticed.
#'
#' @param before row count before the join
#' @param after data.frame after the join
#' @param what label for the message
#' @return `after`, unchanged
warn_if_join_multiplied <- function(before, after, what = "join") {
  if (is.null(after)) return(after)
  if (nrow(after) > before) {
    key <- admin2_key(after)
    dupk <- unique(key[duplicated(key)])
    warning(sprintf(paste0("[admin2_hygiene] %s multiplied rows: %d -> %d. ",
                           "Duplicated key(s): %s"),
                    what, before, nrow(after),
                    paste(utils::head(dupk, 5), collapse = ", ")), call. = FALSE)
  }
  after
}
