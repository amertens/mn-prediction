# =============================================================================
# src/Gambia/gambia_DHS_aggregation_admin2.R
#
# Compute DHS indicator prevalences at Admin-2 level for Gambia using the
# surveyPrev package.  Uses three estimation approaches:
#   1. Direct (design-weighted) estimates
#   2. Fay-Herriot (area-level) smoothed estimates (BYM2)
#   3. Cluster-level model (BYM2, unstratified)
#
# Approach follows the surveyPrev vignette (Dong et al., 2024) and mirrors
# the existing admin-1 workflow in gambia_DHS_aggregation_admin1.R, but
# uses raw microdata + GPS coordinates to assign clusters to Admin-2 areas.
#
# Requirements:
#   - rdhs configured with valid DHS credentials (rdhs::set_rdhs_config())
#   - INLA installed (install.packages("INLA", repos = ...))
#   - Internet access for DHS data download and GADM boundaries
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(sf)
  library(geodata)
  library(surveyPrev)
  library(SUMMER)
  library(rdhs)
  library(ggplot2)
  library(patchwork)
})

# Optional: load INLA if available (needed for Fay-Herriot and cluster models)
has_inla <- requireNamespace("INLA", quietly = TRUE)
if (has_inla) library(INLA) else message("INLA not installed; skipping smoothed models.")

# =============================================================================
# 0. Configuration
# =============================================================================

COUNTRY      <- "Gambia"
COUNTRY_ISO2 <- "GM"
COUNTRY_ISO3 <- "GMB"
SURVEY_YEAR  <- 2019       # Gambia 2019-20 DHS

# Output directory
out_dir <- here("data", "DHS", "clean")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Indicators to estimate (surveyPrev built-in IDs)
# Use surveyPrevIndicators to see the full list; we attempt all 22
data(surveyPrevIndicators)
INDICATOR_IDS <- surveyPrevIndicators$ID
cat(sprintf("Will attempt %d built-in surveyPrev indicators.\n", length(INDICATOR_IDS)))

# =============================================================================
# 1. Download DHS survey data
# =============================================================================
# getDHSdata() downloads the appropriate recode file for each indicator.
# We cache all unique recode types to avoid re-downloading.

cat("\n[01] Downloading DHS survey data for", COUNTRY, SURVEY_YEAR, "...\n")

# We need to download the data for each indicator; getDHSdata handles
# selecting the right recode file (IR, KR, PR, HR, etc.)
# Cache downloaded data to avoid repeated API calls
dhs_data_cache <- list()
indicator_data <- list()  # stores getDHSindicator() output per indicator
indicator_log  <- tibble(indicator = character(), ok = logical(), message = character())

for (ind_id in INDICATOR_IDS) {
  cat(sprintf("  [%s] ", ind_id))

  res <- tryCatch({
    # Step 1: Download raw DHS data for this indicator
    cache_key <- ind_id
    dhsData <- getDHSdata(
      country   = COUNTRY,
      indicator = ind_id,
      year      = SURVEY_YEAR
    )

    # Step 2: Process into standardized binary indicator
    dat <- getDHSindicator(dhsData, indicator = ind_id)

    n_valid <- sum(!is.na(dat$value))
    cat(sprintf("OK (n = %d)\n", n_valid))
    list(ok = TRUE, data = dat, msg = "")
  }, error = function(e) {
    cat(sprintf("SKIP (%s)\n", e$message))
    list(ok = FALSE, data = NULL, msg = e$message)
  })

  indicator_log <- bind_rows(indicator_log, tibble(
    indicator = ind_id, ok = res$ok, message = res$msg
  ))

  if (res$ok) indicator_data[[ind_id]] <- res$data
}

cat(sprintf("\n  Successfully processed: %d / %d indicators\n",
            sum(indicator_log$ok), nrow(indicator_log)))

if (length(indicator_data) == 0) stop("No indicators could be processed. Check DHS credentials and data availability.")

# =============================================================================
# 2. Download GPS data and GADM boundaries
# =============================================================================

cat("\n[02] Downloading GPS coordinates and GADM boundaries...\n")

# GPS cluster locations
geo <- getDHSgeo(country = COUNTRY, year = SURVEY_YEAR)
cat(sprintf("  GPS data: %d clusters\n", nrow(geo)))

# GADM Admin boundaries
poly.adm1 <- geodata::gadm(country = COUNTRY_ISO3, level = 1, path = tempdir())
poly.adm1 <- sf::st_as_sf(poly.adm1)
cat(sprintf("  Admin1 polygons: %d\n", nrow(poly.adm1)))

poly.adm2 <- geodata::gadm(country = COUNTRY_ISO3, level = 2, path = tempdir())
poly.adm2 <- sf::st_as_sf(poly.adm2)
cat(sprintf("  Admin2 polygons: %d\n", nrow(poly.adm2)))

# =============================================================================
# 3. Assign clusters to Admin areas
# =============================================================================

cat("\n[03] Assigning clusters to Admin1/Admin2 areas...\n")

cluster.info <- clusterInfo(
  geo       = geo,
  poly.adm1 = poly.adm1,
  poly.adm2 = poly.adm2,
  by.adm1   = "NAME_1",
  by.adm2   = "NAME_2"
)

n_valid_clusters <- nrow(cluster.info$data)
n_wrong <- length(cluster.info$wrong.points)
cat(sprintf("  Valid clusters: %d | Invalid GPS: %d\n", n_valid_clusters, n_wrong))

if (n_wrong > 0) {
  cat("  Clusters with invalid GPS (excluded):", paste(cluster.info$wrong.points, collapse = ", "), "\n")
}

# Admin info objects (without population — added later if needed)
admin.info1 <- adminInfo(
  poly.adm      = poly.adm1,
  admin         = 1,
  by.adm        = "NAME_1"
)

admin.info2 <- adminInfo(
  poly.adm      = poly.adm2,
  admin         = 2,
  by.adm        = "NAME_2",
  by.adm.upper  = "NAME_1"
)

cat(sprintf("  Admin1 areas: %d | Admin2 areas: %d\n",
            nrow(admin.info1$data), nrow(admin.info2$data)))

# =============================================================================
# 4. Compute Admin-2 direct estimates for all indicators
# =============================================================================

cat("\n[04] Computing Admin-2 direct estimates...\n")

direct_results <- list()
direct_log <- tibble(indicator = character(), ok = logical(), n_admin2 = integer(), message = character())

for (ind_id in names(indicator_data)) {
  dat <- indicator_data[[ind_id]]

  res <- tryCatch({
    est <- directEST(
      data         = dat,
      cluster.info = cluster.info,
      admin        = 2,
      aggregation  = FALSE
    )

    n_est <- sum(!is.na(est$res.admin2$direct.est))
    cat(sprintf("  [%s] %d Admin2 areas with estimates\n", ind_id, n_est))
    list(ok = TRUE, result = est, n = n_est, msg = "")
  }, error = function(e) {
    cat(sprintf("  [%s] FAILED: %s\n", ind_id, e$message))
    list(ok = FALSE, result = NULL, n = 0L, msg = e$message)
  })

  direct_log <- bind_rows(direct_log, tibble(
    indicator = ind_id, ok = res$ok, n_admin2 = as.integer(res$n), message = res$msg
  ))

  if (res$ok) direct_results[[ind_id]] <- res$result
}

cat(sprintf("\n  Direct estimates computed for %d / %d indicators\n",
            sum(direct_log$ok), nrow(direct_log)))

# Also compute Admin-1 and national direct estimates for reference
cat("\n  Computing Admin-1 and national direct estimates for reference...\n")
direct_admin1_results <- list()
direct_national_results <- list()

for (ind_id in names(indicator_data)) {
  dat <- indicator_data[[ind_id]]

  # Admin 1
  tryCatch({
    est1 <- directEST(data = dat, cluster.info = cluster.info, admin = 1)
    direct_admin1_results[[ind_id]] <- est1
  }, error = function(e) NULL)

  # National
  tryCatch({
    est0 <- directEST(data = dat, cluster.info = cluster.info, admin = 0, strata = "all")
    direct_national_results[[ind_id]] <- est0
  }, error = function(e) NULL)
}

# =============================================================================
# 5. Smoothed estimates (Fay-Herriot BYM2) at Admin-2
# =============================================================================

fh_results <- list()

if (has_inla) {
  cat("\n[05] Computing Fay-Herriot (BYM2) smoothed estimates at Admin-2...\n")

  fh_log <- tibble(indicator = character(), ok = logical(), message = character())

  for (ind_id in names(direct_results)) {
    dat <- indicator_data[[ind_id]]

    res <- tryCatch({
      # Identify and exclude regions with near-zero design variance
      dir_est <- direct_results[[ind_id]]$res.admin2
      bad_admin2 <- dir_est$admin2.name.full[
        !is.na(dir_est$direct.var) & dir_est$direct.var < 1e-30
      ]
      bad_clusters <- cluster.info$data$cluster[
        cluster.info$data$admin2.name.full %in% bad_admin2
      ]

      dat_clean <- if (length(bad_clusters) > 0) {
        subset(dat, !cluster %in% bad_clusters)
      } else {
        dat
      }

      est <- fhModel(
        data         = dat_clean,
        cluster.info = cluster.info,
        admin.info   = admin.info2,
        admin        = 2,
        model        = "bym2",
        aggregation  = FALSE
      )

      n_est <- sum(!is.na(est$res.admin2$mean))
      cat(sprintf("  [%s] FH: %d Admin2 areas (removed %d bad clusters)\n",
                  ind_id, n_est, length(bad_clusters)))
      list(ok = TRUE, result = est, msg = "")
    }, error = function(e) {
      cat(sprintf("  [%s] FH FAILED: %s\n", ind_id, e$message))
      list(ok = FALSE, result = NULL, msg = e$message)
    })

    fh_log <- bind_rows(fh_log, tibble(
      indicator = ind_id, ok = res$ok, message = res$msg
    ))

    if (res$ok) fh_results[[ind_id]] <- res$result
  }

  cat(sprintf("\n  Fay-Herriot estimates computed for %d / %d indicators\n",
              sum(fh_log$ok), nrow(fh_log)))
} else {
  cat("\n[05] Skipping Fay-Herriot models (INLA not installed).\n")
}

# =============================================================================
# 6. Cluster-level model (BYM2, unstratified) at Admin-2
# =============================================================================

cluster_results <- list()

if (has_inla) {
  cat("\n[06] Computing cluster-level (BYM2, unstratified) estimates at Admin-2...\n")

  cl_log <- tibble(indicator = character(), ok = logical(), message = character())

  for (ind_id in names(indicator_data)) {
    dat <- indicator_data[[ind_id]]

    res <- tryCatch({
      est <- clusterModel(
        data           = dat,
        cluster.info   = cluster.info,
        admin.info     = admin.info2,
        stratification = FALSE,
        model          = "bym2",
        admin          = 2,
        aggregation    = FALSE,
        CI             = 0.95
      )

      n_est <- sum(!is.na(est$res.admin2$mean))
      cat(sprintf("  [%s] Cluster: %d Admin2 areas\n", ind_id, n_est))
      list(ok = TRUE, result = est, msg = "")
    }, error = function(e) {
      cat(sprintf("  [%s] Cluster FAILED: %s\n", ind_id, e$message))
      list(ok = FALSE, result = NULL, msg = e$message)
    })

    cl_log <- bind_rows(cl_log, tibble(
      indicator = ind_id, ok = res$ok, message = res$msg
    ))

    if (res$ok) cluster_results[[ind_id]] <- res$result
  }

  cat(sprintf("\n  Cluster-level estimates computed for %d / %d indicators\n",
              sum(cl_log$ok), nrow(cl_log)))
} else {
  cat("\n[06] Skipping cluster-level models (INLA not installed).\n")
}

# =============================================================================
# 7. Compile results into tidy data frames
# =============================================================================

cat("\n[07] Compiling results...\n")

# --- 7a. Direct estimates wide table ---
# Pivot: one row per Admin2 area, one column per indicator
compile_direct_wide <- function(results, admin_level = 2) {
  rows <- list()
  for (ind_id in names(results)) {
    if (admin_level == 2) {
      df <- results[[ind_id]]$res.admin2
      df$indicator <- ind_id
    } else {
      df <- results[[ind_id]]$res.admin1
      df$indicator <- ind_id
    }
    rows[[ind_id]] <- df
  }
  if (length(rows) == 0) return(tibble())
  bind_rows(rows)
}

# Long-format table with all direct estimates
direct_long <- compile_direct_wide(direct_results, admin_level = 2)
cat(sprintf("  Direct (long): %d rows × %d indicators\n",
            nrow(direct_long), length(unique(direct_long$indicator))))

# Wide-format: one row per admin2, columns = indicator prevalences
if (nrow(direct_long) > 0) {
  direct_wide <- direct_long %>%
    select(admin2.name.full, indicator, direct.est, direct.se) %>%
    tidyr::pivot_wider(
      names_from  = indicator,
      values_from = c(direct.est, direct.se),
      names_glue  = "{indicator}_{.value}"
    )
  cat(sprintf("  Direct (wide): %d Admin2 areas × %d columns\n",
              nrow(direct_wide), ncol(direct_wide)))
} else {
  direct_wide <- tibble()
}

# --- 7b. Smoothed estimates long table ---
if (length(fh_results) > 0) {
  fh_rows <- list()
  for (ind_id in names(fh_results)) {
    df <- fh_results[[ind_id]]$res.admin2
    df$indicator <- ind_id
    fh_rows[[ind_id]] <- df
  }
  fh_long <- bind_rows(fh_rows)
  cat(sprintf("  Fay-Herriot (long): %d rows × %d indicators\n",
              nrow(fh_long), length(unique(fh_long$indicator))))
} else {
  fh_long <- tibble()
}

if (length(cluster_results) > 0) {
  cl_rows <- list()
  for (ind_id in names(cluster_results)) {
    df <- cluster_results[[ind_id]]$res.admin2
    df$indicator <- ind_id
    cl_rows[[ind_id]] <- df
  }
  cl_long <- bind_rows(cl_rows)
  cat(sprintf("  Cluster-level (long): %d rows × %d indicators\n",
              nrow(cl_long), length(unique(cl_long$indicator))))
} else {
  cl_long <- tibble()
}

# =============================================================================
# 8. Visualization: maps for a sample indicator
# =============================================================================

cat("\n[08] Generating sample maps...\n")

# Create admin2.name.full in polygon for mapping
poly.adm2$admin2.name.full <- paste0(poly.adm2$NAME_1, "_", poly.adm2$NAME_2)

# Pick the first indicator with valid direct estimates for demo maps
map_ind <- if (length(direct_results) > 0) names(direct_results)[1] else NULL

if (!is.null(map_ind)) {
  cat(sprintf("  Mapping indicator: %s\n", map_ind))

  # Collect available estimates for this indicator
  map_data <- list()

  # Direct
  out_dir_est <- direct_results[[map_ind]]$res.admin2 %>%
    select(admin2.name.full, mean = direct.est, cv) %>%
    mutate(model = "Direct")
  map_data[["Direct"]] <- out_dir_est

  # Fay-Herriot
  if (map_ind %in% names(fh_results)) {
    out_fh <- fh_results[[map_ind]]$res.admin2 %>%
      select(admin2.name.full, mean, cv) %>%
      mutate(model = "Fay-Herriot (BYM2)")
    map_data[["FH"]] <- out_fh
  }

  # Cluster
  if (map_ind %in% names(cluster_results)) {
    out_cl <- cluster_results[[map_ind]]$res.admin2 %>%
      select(admin2.name.full, mean, cv) %>%
      mutate(model = "Cluster (BYM2)")
    map_data[["Cluster"]] <- out_cl
  }

  map_combined <- bind_rows(map_data)

  # Get indicator label
  ind_label <- surveyPrevIndicators$Description[surveyPrevIndicators$ID == map_ind]
  if (length(ind_label) == 0) ind_label <- map_ind

  # Mean map
  tryCatch({
    g_mean <- SUMMER::mapPlot(
      data      = map_combined,
      geo       = poly.adm2,
      by.data   = "admin2.name.full",
      by.geo    = "admin2.name.full",
      is.long   = TRUE,
      variable  = "model",
      value     = "mean",
      legend.label = "Prevalence"
    ) + ggtitle(paste0(ind_label, "\n(", COUNTRY, " ", SURVEY_YEAR, ")"))

    fig_path <- file.path(here("results", "figures"),
                          paste0(COUNTRY_ISO2, "_admin2_", map_ind, "_comparison.png"))
    dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(fig_path, g_mean, width = 14, height = 6, dpi = 150)
    cat(sprintf("  Saved: %s\n", basename(fig_path)))
  }, error = function(e) {
    cat(sprintf("  Map plotting failed: %s\n", e$message))
  })
}

# =============================================================================
# 9. Save outputs
# =============================================================================

cat("\n[09] Saving outputs...\n")

# Direct estimates (long format — most useful for downstream analysis)
if (nrow(direct_long) > 0) {
  f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin2_direct.rds"))
  saveRDS(direct_long, file = f)
  cat(sprintf("  Saved: %s\n", basename(f)))

  f_csv <- sub("\\.rds$", ".csv", f)
  write.csv(direct_long, file = f_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", basename(f_csv)))
}

# Direct estimates (wide — one row per admin2)
if (nrow(direct_wide) > 0) {
  f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin2_direct_wide.rds"))
  saveRDS(direct_wide, file = f)
  cat(sprintf("  Saved: %s\n", basename(f)))
}

# Fay-Herriot smoothed (long)
if (nrow(fh_long) > 0) {
  f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin2_fh_bym2.rds"))
  saveRDS(fh_long, file = f)
  cat(sprintf("  Saved: %s\n", basename(f)))
}

# Cluster-level (long)
if (nrow(cl_long) > 0) {
  f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin2_cluster_bym2.rds"))
  saveRDS(cl_long, file = f)
  cat(sprintf("  Saved: %s\n", basename(f)))
}

# Admin-1 direct (for comparison / aggregation checks)
direct_admin1_long <- compile_direct_wide(direct_admin1_results, admin_level = 1)
if (nrow(direct_admin1_long) > 0) {
  f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin1_direct.rds"))
  saveRDS(direct_admin1_long, file = f)
  cat(sprintf("  Saved: %s\n", basename(f)))
}

# Cluster info and admin info (useful for downstream modeling)
f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_cluster_admin_info.rds"))
saveRDS(list(
  cluster.info = cluster.info,
  admin.info1  = admin.info1,
  admin.info2  = admin.info2,
  poly.adm1    = poly.adm1,
  poly.adm2    = poly.adm2
), file = f)
cat(sprintf("  Saved: %s\n", basename(f)))

# Processing log
f <- file.path(out_dir, paste0(COUNTRY, "_", SURVEY_YEAR, "_dhs_admin2_log.rds"))
saveRDS(list(
  indicator_log = indicator_log,
  direct_log    = direct_log,
  fh_log        = if (exists("fh_log")) fh_log else NULL,
  cl_log        = if (exists("cl_log")) cl_log else NULL
), file = f)
cat(sprintf("  Saved: %s\n", basename(f)))

# =============================================================================
# 10. Summary
# =============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("  Country: %s | Survey: %d\n", COUNTRY, SURVEY_YEAR))
cat(sprintf("  Indicators attempted: %d | Processed: %d\n",
            length(INDICATOR_IDS), length(indicator_data)))
cat(sprintf("  Admin2 areas: %d | Clusters: %d (invalid GPS: %d)\n",
            nrow(admin.info2$data), n_valid_clusters, n_wrong))
cat(sprintf("  Direct estimates: %d indicators\n", length(direct_results)))
cat(sprintf("  Fay-Herriot (BYM2): %d indicators\n", length(fh_results)))
cat(sprintf("  Cluster-level (BYM2): %d indicators\n", length(cluster_results)))
cat("\nDone.\n")
