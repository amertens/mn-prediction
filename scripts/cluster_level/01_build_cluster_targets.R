# =============================================================================
# scripts/cluster_level/01_build_cluster_targets.R   [CL-01]
#
# CLUSTER-LEVEL OUTCOME TABLE: the parallel analysis's first artefact.
#
# The district and regional estimations are untouched (protocol v2,
# scripts/protocol_v2/). This directory builds the same targets one rung
# down, at the survey cluster, where the outcome is measured and where a
# covariate can be extracted at the place the respondents live rather than
# averaged over a district most of them do not resemble. Design and the
# reasons are in docs/findings/CLUSTER_LEVEL_DESIGN_2026-09.md.
#
# For every country x outcome x cluster:
#   n_raw, n_kish     respondents with the outcome; Kish effective n from weights
#   y_prev            survey-weighted prevalence of the UNIFORM outcome
#                     (resolve_uniform_outcome, as targets_v2)
#   y_level, sd_level survey-weighted mean and SD of the negated log biomarker
#   n_eff, n_eff_cont n_raw / national design effect (deff_v2.csv; the cluster
#                     IS the PSU, so the clustering part of deff is 1 here and
#                     only the weighting part applies: n_eff = n_kish)
#   Admin1, Admin2    the cluster's district (modal, if respondents disagree)
#   lat, lon          cluster GPS (data/IPD/<country>/*_GPS_cleaned.csv)
#   urban             GPS-file urban/rural flag where the survey carries one
#   date_med          fieldwork date (FW-01, fieldwork_windows_cluster.csv)
#
#   Rscript scripts/cluster_level/01_build_cluster_targets.R
# -> results/tables/cluster_level/targets_cluster.csv
# -> results/tables/cluster_level/cluster_gps_match.csv (join diagnostics)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
STORE <- "_targets_full"; OUTDIR <- "results/tables/cluster_level"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
GPS <- list(Gambia = list(f = "data/IPD/Gambia/Gambia_GMS_GPS_cleaned.csv", id = "MICS_Cluster_Number", urban = NULL),
            Ghana = list(f = "data/IPD/Ghana/Ghana_GMS_GPS_cleaned.csv", id = "cnum", urban = "urban/ rural"),
            Malawi = list(f = "data/IPD/Malawi/Malawi_GMS_GPS_cleaned.csv", id = "gw_cnum", urban = NULL),
            SierraLeone = list(f = "data/IPD/Sierra Leone/Sierra Leone_GMS_GPS_cleaned.csv", id = "cnum", urban = NULL))
FW <- tryCatch(read.csv("results/tables/protocol_v2/fieldwork_windows_cluster.csv", stringsAsFactors = FALSE), error = function(e) NULL)
cfgs <- get_country_configs()
modal <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] }
rows <- list(); diag <- list()
for (lc in names(COUNTRIES)) { cn <- COUNTRIES[[lc]]; cc <- cfgs[[cn]]
  g <- GPS[[cn]]; gp <- read.csv(g$f, stringsAsFactors = FALSE, check.names = FALSE)
  gp <- data.frame(cluster = as.character(gp[[g$id]]), lat = as.numeric(gp$latitude), lon = as.numeric(gp$longitude),
                   urban = if (!is.null(g$urban) && g$urban %in% names(gp)) as.integer(tolower(trimws(gp[[g$urban]])) == "urban") else NA_integer_, stringsAsFactors = FALSE)
  gp <- gp[is.finite(gp$lat) & is.finite(gp$lon), ]; gp <- gp[!duplicated(gp$cluster), ]
  for (on in names(cc$outcomes)) { oc <- cc$outcomes[[on]]
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL); if (is.null(od)) next
    d <- od$data; if (!cc$cluster_id %in% names(d)) next
    w <- .v2_num(d[[cc$weight_col]]); w[!is.finite(w) | w <= 0] <- NA
    yb <- tryCatch(resolve_uniform_outcome(d, cc, oc, label = "[cl]"), error = function(e) NULL)
    ybin <- if (!is.null(yb)) .v2_num(yb) else .v2_num(d[[oc$binary]])
    ycont <- rep(NA_real_, nrow(d))
    if (!is.null(oc$continuous) && oc$continuous %in% names(d)) { v <- .v2_num(d[[oc$continuous]])
      t <- if (identical(oc$cutoff_scale, "log")) v else { v[!is.finite(v) | v <= 0] <- NA; log(v) }; ycont <- -t }
    dd <- data.frame(cluster = as.character(d[[cc$cluster_id]]), Admin1 = as.character(d$Admin1), Admin2 = as.character(d$Admin2), ybin = ybin, ycont = ycont, w = w, stringsAsFactors = FALSE)
    agg <- dd |> filter(!is.na(cluster)) |> group_by(cluster) |> summarise(
      Admin1 = modal(Admin1), Admin2 = modal(Admin2), n_raw = sum(is.finite(ybin)), n_raw_cont = sum(is.finite(ycont)),
      n_kish = kish_n_v2(w[is.finite(ybin)]), n_kish_cont = kish_n_v2(w[is.finite(ycont)]),
      y_prev = .v2_wmean(ybin, w), y_level = .v2_wmean(ycont, w), sd_level = sqrt(.v2_wvar(ycont, w)), .groups = "drop")
    agg <- agg[agg$n_raw > 0 | agg$n_raw_cont > 0, ]
    # the cluster is the PSU, so only the weighting component of the design effect applies
    agg$n_eff <- pmax(1, agg$n_kish); agg$n_eff_cont <- pmax(1, agg$n_kish_cont)
    agg <- left_join(agg, gp, by = "cluster")
    if (!is.null(FW)) { f <- FW[FW$country == cn, c("cluster", "date_med", "month_med")]; f$cluster <- as.character(f$cluster); agg <- left_join(agg, f, by = "cluster") }
    agg$country <- cn; agg$outcome <- on
    rows[[paste(cn, on)]] <- agg
    diag[[paste(cn, on)]] <- data.frame(country = cn, outcome = on, clusters_in_data = nrow(agg), with_gps = sum(is.finite(agg$lat)), gps_file_clusters = nrow(gp),
                                        gps_unmatched = length(setdiff(gp$cluster, agg$cluster)), median_n = stats::median(agg$n_raw), stringsAsFactors = FALSE)
    cat(sprintf("  %-12s %-13s clusters %3d | with GPS %3d | median n %2d | prevalence %.3f\n", cn, on, nrow(agg), sum(is.finite(agg$lat)), stats::median(agg$n_raw), .v2_wmean(agg$y_prev, agg$n_raw)))
  }
}
TC <- bind_rows(rows); DG <- bind_rows(diag)
TC <- TC[, c("country", "outcome", "cluster", "Admin1", "Admin2", "lat", "lon", "urban", "date_med", "month_med", "n_raw", "n_raw_cont", "n_kish", "n_eff", "n_eff_cont", "y_prev", "y_level", "sd_level")]
write.csv(TC, file.path(OUTDIR, "targets_cluster.csv"), row.names = FALSE); write.csv(DG, file.path(OUTDIR, "cluster_gps_match.csv"), row.names = FALSE)
cat("\n===== CL-01: cluster targets =====\n"); print(as.data.frame(DG), row.names = FALSE)
cat(sprintf("\nrows %d | clusters with GPS %d of %d | countries %s\n", nrow(TC), sum(is.finite(TC$lat[!duplicated(paste(TC$country, TC$cluster))])), sum(!duplicated(paste(TC$country, TC$cluster))), paste(unique(TC$country), collapse = ", ")))
cat("\nDONE\n")
