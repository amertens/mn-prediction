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

#' Names GADM uses for inland water polygons.
#' Deliberately anchored so a real district called e.g. "Lakeside" is not caught.
ADMIN2_WATER_PATTERN <- "^Lake |^Water|^Lac |^Sea$|Reservoir"

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
  if (!any(duplicated(d$Admin2))) return(d)
  dup <- unique(d$Admin2[duplicated(d$Admin2)])
  cat(sprintf("[admin2_hygiene] %s: collapsing %d duplicated Admin-2 name(s): %s\n",
              what, length(dup), paste(utils::head(dup, 6), collapse = ", ")))
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  oth <- setdiff(names(d), c("Admin2", num))
  out <- d |>
    dplyr::group_by(Admin2) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(oth), dplyr::first),
      dplyr::across(dplyr::all_of(num), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop") |>
    as.data.frame()
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
  if (any(duplicated(d$Admin2)))
    stop(sprintf("%s has %d duplicated Admin-2 key(s) (%s) - run clean_admin2_keys() first",
                 what, sum(duplicated(d$Admin2)),
                 paste(utils::head(unique(d$Admin2[duplicated(d$Admin2)]), 5),
                       collapse = ", ")))
  invisible(d)
}
