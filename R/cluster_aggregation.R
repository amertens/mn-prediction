# =============================================================================
# R/cluster_aggregation.R
#
# Aggregate biomarker outcomes to SURVEY-CLUSTER GPS buffers — a geographic
# unit finer than admin-2 — following Bondi-Kelly et al. (2022), who aggregate
# individuals to satellite pixels. Each survey cluster is one GPS point with a
# handful of sampled individuals; we give it a survey-weighted prevalence label
# and re-extract GEE covariates on a small buffer around its location.
#
# WHY CLUSTERS, NOT 25 m PIXELS: like the paper (where only ~0.02% of pixels
# contained any individuals), almost every fine grid cell would be unlabeled —
# only cells containing a survey cluster get an outcome. So the survey cluster
# IS the natural finest labeled unit; a buffer around it sets the spatial scale
# over which covariates are averaged (and absorbs minor GPS imprecision).
#
# The output is a DROP-IN for the admin-2 comparators: it carries svy_prev,
# n_svy, lon, lat, Admin1, and gee_* columns named identically to the admin-2
# extracts, so run_area_benchmarks_within / _loro / _loco and every fit_predict_*
# work on it unchanged.
#
# DATA NOTE: per-cluster GPS lives in a separate CSV per country (e.g.
# data/IPD/Gambia/Gambia_GMS_GPS_cleaned.csv), NOT in the merged RDS (the merge
# step drops coordinates into geometry and then drops geometry). load_cluster_gps()
# rejoins it, auto-detecting the cluster-key / lon / lat columns.
# =============================================================================

#' Load and join per-cluster GPS coordinates for a country.
#'
#' Auto-detects the lon/lat columns and the cluster-key column (the one whose
#' values best overlap the merged data's cluster ids), so it works across the
#' countries' differing GPS file layouts without hard-coding column names.
#'
#' @param cc Country config (uses data_path's directory + cluster_id)
#' @param merged_cluster_ids Optional vector of cluster ids from the merged data,
#'   used to pick the join key by maximum overlap.
#' @return data.frame(cluster_id, lon, lat) or NULL if no usable GPS file.
load_cluster_gps <- function(cc, merged_cluster_ids = NULL) {
  ddir <- dirname(cc$data_path)
  cand <- list.files(ddir, pattern = "gps.*\\.csv$", full.names = TRUE,
                     ignore.case = TRUE)
  if (length(cand) == 0) {
    cat(sprintf("[cluster_gps] %s: no GPS csv found in %s\n", cc$country, ddir))
    return(NULL)
  }
  g <- utils::read.csv(cand[1], stringsAsFactors = FALSE, check.names = TRUE)
  nm <- names(g)
  # Prefer explicitly-named longitude/latitude columns; only fall back to the
  # ambiguous lon/lng/x (resp. lat/y) forms if no full name exists. (A bare "X"
  # is often a row-index column, so it must NOT win over "longitude".)
  lon <- nm[grepl("longitude", nm, ignore.case = TRUE)][1]
  if (is.na(lon)) lon <- nm[grepl("^lon$|^lng$|^x_?coord", nm, ignore.case = TRUE)][1]
  if (is.na(lon)) lon <- nm[grepl("^x$", nm, ignore.case = TRUE)][1]
  lat <- nm[grepl("latitude", nm, ignore.case = TRUE)][1]
  if (is.na(lat)) lat <- nm[grepl("^lat$|^y_?coord", nm, ignore.case = TRUE)][1]
  if (is.na(lat)) lat <- nm[grepl("^y$", nm, ignore.case = TRUE)][1]
  if (is.na(lon) || is.na(lat)) {
    cat(sprintf("[cluster_gps] %s: could not find lon/lat columns (%s)\n",
                cc$country, paste(nm, collapse = ", ")))
    return(NULL)
  }
  # Choose the cluster key: column whose values best overlap merged ids.
  key <- NA_character_
  if (!is.null(merged_cluster_ids)) {
    best <- 0L
    for (k in nm) {
      ov <- length(intersect(as.character(merged_cluster_ids),
                             as.character(g[[k]])))
      if (ov > best) { best <- ov; key <- k }
    }
    if (best == 0) key <- NA_character_
  }
  if (is.na(key)) key <- nm[grepl("cluster", nm, ignore.case = TRUE)][1]
  if (is.na(key)) {
    cat(sprintf("[cluster_gps] %s: could not identify a cluster-key column\n", cc$country))
    return(NULL)
  }
  out <- data.frame(cluster_id = as.character(g[[key]]),
                    lon = suppressWarnings(as.numeric(g[[lon]])),
                    lat = suppressWarnings(as.numeric(g[[lat]])),
                    stringsAsFactors = FALSE)
  out <- out[is.finite(out$lon) & is.finite(out$lat) & !is.na(out$cluster_id) &
             abs(out$lon) <= 180 & abs(out$lat) <= 90, ]
  out <- out[!duplicated(out$cluster_id), ]
  cat(sprintf("[cluster_gps] %s: %d cluster coords (key=%s, lon=%s, lat=%s)\n",
              cc$country, nrow(out), key, lon, lat))
  out
}


#' Aggregate individual biomarker outcomes to survey-cluster prevalences.
#'
#' Mirrors compute_svy_admin2()'s outcome handling (population filter + uniform
#' WHO-cutoff override for transport tags), but bins to the survey cluster and
#' joins GPS coordinates. Cluster prevalence is the survey-weighted mean of the
#' binary outcome among that cluster's sampled individuals.
#'
#' @return data.frame(cluster_id, lon, lat, Admin1, Admin2, svy_prev, n_svy)
aggregate_outcome_to_clusters <- function(cc, oc, min_per_cluster = 1L) {
  if (!file.exists(cc$data_path)) {
    warning("[cluster] merged data not found: ", cc$data_path); return(NULL)
  }
  d <- readRDS(cc$data_path)
  if (is.list(d) && !is.data.frame(d) && is.data.frame(d$data)) d <- d$data
  cid <- cc$cluster_id; bin <- oc$binary
  if (is.null(cid) || !cid %in% names(d) || !bin %in% names(d)) {
    warning("[cluster] missing cluster id or outcome column"); return(NULL)
  }

  # Population filter (children vs women).
  if (!is.null(cc$child_flag) && cc$child_flag %in% names(d) &&
      !is.null(oc$child_flag_val)) {
    d <- d[!is.na(d[[cc$child_flag]]) & d[[cc$child_flag]] == oc$child_flag_val,
           , drop = FALSE]
  }

  # Uniform cross-country outcome override (match compute_svy_admin2()).
  if (!is.null(oc$tag) && exists("UNIFORM_TRANSPORT_TAGS") &&
      oc$tag %in% UNIFORM_TRANSPORT_TAGS &&
      !is.null(oc$continuous) && oc$continuous %in% names(d) &&
      !is.null(oc$cutoff)) {
    d[[bin]] <- apply_threshold(suppressWarnings(as.numeric(d[[oc$continuous]])),
                                oc$cutoff, oc$cutoff_dir %||% "less")
  }

  wcol <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d))
            cc$weight_col else NA_character_
  a1 <- cc$admin1_col; a2 <- cc$admin2_col
  d <- d[!is.na(d[[bin]]) & !is.na(d[[cid]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$.bin <- as.numeric(d[[bin]])
  d$.w   <- if (!is.na(wcol)) suppressWarnings(as.numeric(d[[wcol]])) else 1
  d$.w[!is.finite(d$.w)] <- 1

  parts <- split(d, as.character(d[[cid]]))
  rows <- lapply(names(parts), function(k) {
    s <- parts[[k]]
    if (nrow(s) < min_per_cluster) return(NULL)
    data.frame(
      cluster_id = k,
      Admin1 = if (!is.null(a1) && a1 %in% names(s))
                 as.character(s[[a1]][which(!is.na(s[[a1]]))[1]]) else NA_character_,
      Admin2 = if (!is.null(a2) && a2 %in% names(s))
                 as.character(s[[a2]][which(!is.na(s[[a2]]))[1]]) else NA_character_,
      svy_prev = sum(s$.bin * s$.w, na.rm = TRUE) / sum(s$.w, na.rm = TRUE),
      n_svy = nrow(s), stringsAsFactors = FALSE)
  })
  cl <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(cl) || nrow(cl) == 0) return(NULL)

  # Join GPS.
  gps <- load_cluster_gps(cc, merged_cluster_ids = unique(cl$cluster_id))
  if (is.null(gps)) { warning("[cluster] no GPS — cannot place clusters"); return(NULL) }
  cl <- merge(cl, gps, by = "cluster_id", all.x = TRUE)
  cl <- cl[is.finite(cl$lon) & is.finite(cl$lat), , drop = FALSE]
  cat(sprintf("[cluster] %s/%s: %d clusters with GPS (median n/cluster = %d)\n",
              cc$country, oc$tag, nrow(cl), stats::median(cl$n_svy)))
  cl[, c("cluster_id", "lon", "lat", "Admin1", "Admin2", "svy_prev", "n_svy")]
}


#' Canonical gee_ variable name from a raster filename (mirrors the admin-2
#' extractor so cluster columns are named identically and stay drop-in).
.cluster_gee_varname <- function(tif) {
  base <- tools::file_path_sans_ext(basename(tif))
  for (cname in c("Gambia", "Ghana", "Sierra_Leone", "Sierra Leone",
                  "Malawi", "Cote_dIvoire", "Cote d'Ivoire")) {
    base <- gsub(paste0("_?", gsub("[' ]", ".", cname), "_?"), "_", base)
  }
  v <- paste0("gee_", tolower(gsub("[^A-Za-z0-9]+", "_", base)))
  sub("^_|_$", "", gsub("_+", "_", v))
}

#' Re-extract GEE covariates on buffers around each cluster.
#'
#' Builds a `buffer_km` buffer around each cluster GPS point and zonal-means the
#' rasters onto it, producing gee_* columns named identically to the admin-2
#' extracts (so the result is a drop-in for the comparators).
#'
#' Crash-safety: the admin-2 extractor loops exact_extract once PER BAND, which
#' segfaults on the project's huge multi-band rasters (e.g. Atmosphere = 628
#' daily bands -> 628 calls). Here, any raster with more than `max_bands_individual`
#' layers is first COLLAPSED across bands with terra (one mean/sd layer) and then
#' extracted once — same annual-summary semantics, no 628-call blow-up.
extract_gee_cluster_buffers <- function(cc, clusters, buffer_km = 2,
                                        raster_dir = NULL,
                                        max_bands_individual = 24L) {
  if (is.null(raster_dir)) raster_dir <- .resolve_raster_dir(cc$raster_dir)
  if (is.null(raster_dir)) { warning("[cluster] no raster dir"); return(clusters) }
  pts  <- sf::st_as_sf(clusters, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  bufs <- sf::st_buffer(pts, dist = buffer_km * 1000)   # s2 geodesic buffer (metres)
  tifs <- sort(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
  cat(sprintf("[cluster_gee] %s: %d rasters x %d buffers (%gkm)\n",
              cc$country, length(tifs), nrow(bufs), buffer_km))

  xx <- function(layer) tryCatch(
    exactextractr::exact_extract(layer, bufs, "mean", progress = FALSE),
    error = function(e) rep(NA_real_, nrow(bufs)))

  cols <- list()
  for (tif in tifs) {
    base <- .cluster_gee_varname(tif)
    r <- tryCatch(terra::rast(tif), error = function(e) NULL)
    if (is.null(r)) next
    nl <- terra::nlyr(r)
    if (nl == 1) {
      cols[[base]] <- xx(r)
    } else if (nl <= max_bands_individual) {
      lyr_names <- tryCatch(names(r), error = function(e) NULL)
      mat <- vapply(seq_len(nl), function(j) xx(r[[j]]), numeric(nrow(bufs)))
      for (j in seq_len(nl)) {
        nmj <- if (!is.null(lyr_names) && nzchar(lyr_names[j]))
                 tolower(gsub("[^A-Za-z0-9]+", "_", lyr_names[j]))
               else if (nl == 12) tolower(month.abb[j]) else sprintf("b%02d", j)
        cols[[sub("_$", "", gsub("_+", "_", paste0(base, "_", nmj)))]] <- mat[, j]
      }
      cols[[paste0(base, "_annual_mean")]] <- rowMeans(mat, na.rm = TRUE)
      cols[[paste0(base, "_annual_sd")]]   <- apply(mat, 1, sd,  na.rm = TRUE)
      cols[[paste0(base, "_annual_min")]]  <- apply(mat, 1, min, na.rm = TRUE)
      cols[[paste0(base, "_annual_max")]]  <- apply(mat, 1, max, na.rm = TRUE)
    } else {
      # Big stack (e.g. 628 daily bands): collapsing ALL bands with terra::app
      # exhausts memory on larger countries. Instead SUBSAMPLE ~12 evenly-spaced
      # bands and summarise those — a fine annual approximation, memory-bounded.
      idx <- unique(round(seq(1, nl, length.out = min(nl, 12L))))
      mat <- tryCatch(
        vapply(idx, function(j) xx(r[[j]]), numeric(nrow(bufs))),
        error = function(e) NULL)
      if (!is.null(mat) && is.matrix(mat)) {
        cols[[paste0(base, "_annual_mean")]] <- rowMeans(mat, na.rm = TRUE)
        cols[[paste0(base, "_annual_sd")]]   <- apply(mat, 1, sd, na.rm = TRUE)
      }
    }
  }
  if (length(cols) == 0) return(clusters)
  gee_df <- as.data.frame(cols, check.names = FALSE)
  cat(sprintf("[cluster_gee] %s: %d gee_ covariates extracted\n",
              cc$country, ncol(gee_df)))
  cbind(clusters, gee_df)
}


#' Build a cluster-resolution modeling dataset (outcome + buffered covariates).
#'
#' @return drop-in data.frame for the area-level comparators: cluster_id, lon,
#'   lat, Admin1, Admin2, svy_prev, n_svy, gee_*; plus a `country` column.
build_cluster_dataset <- function(cc, oc, buffer_km = 2) {
  cl <- aggregate_outcome_to_clusters(cc, oc)
  if (is.null(cl) || nrow(cl) < 5) return(NULL)
  out <- extract_gee_cluster_buffers(cc, cl, buffer_km = buffer_km)
  out$country <- cc$country
  out
}
