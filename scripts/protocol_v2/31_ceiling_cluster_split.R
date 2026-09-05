# =============================================================================
# scripts/protocol_v2/31_ceiling_cluster_split.R   [CE-01]
#
# HOW MUCH OF THE RELIABILITY CEILING IS CLUSTER EFFECT?
#
# The empirical ceiling of record (split-half within each unit, Spearman-Brown
# corrected) treats agreement between two halves of a unit's respondents as
# evidence of real district geography. But 57-85% of Admin-2 units in Gambia,
# Ghana and Malawi contain a single survey cluster, so the two halves share
# everything a cluster shares -- village, field team, day, assay batch -- and
# that shared component is counted as signal. The ceiling and the "headroom"
# are inflated by it wherever units are single-cluster.
#
# This script sizes the inflation on the units where it CAN be sized: those
# with at least two clusters. On the same units it computes
#   within   halves formed by splitting RESPONDENTS at random within the unit
#   cluster  halves formed by splitting CLUSTERS at random within the unit
# Both give a half-half correlation across units, Spearman-Brown corrected to
# the ceiling for the full unit. The cluster split removes what clusters
# share; the gap between the two is the cluster-effect share of the ceiling.
# Sierra Leone (median 4 clusters per district) is the clean case.
#
# Outcome definition is the project's uniform one (resolve_uniform_outcome).
#
#   Rscript scripts/protocol_v2/31_ceiling_cluster_split.R
# -> results/tables/protocol_v2/ceiling_cluster_split.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
OUTDIR <- "results/tables/protocol_v2"; STORE <- "_targets_full"; NSPLIT <- 200L; set.seed(20260903L)
cfg <- get_country_configs()

sb <- function(r) {
  # NA in (e.g. the cluster split left < 8 units with two usable halves) -> NA out;
  # the first run died here on a Malawi cell with `if (NA)`
  if (!is.finite(r)) return(NA_real_)
  r <- max(min(r, 0.999), -0.999); rel <- 2 * r / (1 + r); if (rel <= 0) 0 else sqrt(rel) }
half_corr <- function(y, w, unit, grp) {
  # grp: per-row half label (1/2) within unit; returns Pearson r of unit prevalences across halves
  a <- tapply(seq_along(y), unit, function(i) {
    i1 <- i[grp[i] == 1]; i2 <- i[grp[i] == 2]
    if (length(i1) < 2 || length(i2) < 2) return(c(NA, NA))
    c(stats::weighted.mean(y[i1], w[i1]), stats::weighted.mean(y[i2], w[i2])) })
  m <- do.call(rbind, a); m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) < 8) return(NA_real_)
  suppressWarnings(stats::cor(m[, 1], m[, 2]))
}

rows <- list()
for (cn in names(cfg)) { cc <- cfg[[cn]]; lc <- tolower(cn); clu <- cc$cluster_id; a2 <- cc$admin2_col %||% "Admin2"
  for (on in names(cc$outcomes)) { oc <- cc$outcomes[[on]]
    od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL)
    if (is.null(od)) { cat(sprintf("  %-12s %-12s skip: not in store
", cn, on)); next }
    d <- od$data; if (!all(c(clu, a2) %in% names(d))) { cat(sprintf("  %-12s %-12s skip: missing %s/%s
", cn, on, clu, a2)); next }
    # plain call; the first version wrapped this in capture.output() and lost
    # the result, so every cell was skipped without a message
    y <- tryCatch(resolve_uniform_outcome(d, cc, oc), error = function(e) { cat("  resolver error:", cn, on, conditionMessage(e), "
"); NULL })
    if (is.null(y) || all(is.na(y))) { cat(sprintf("  %-12s %-12s skip: resolver returned nothing
", cn, on)); next }
    y <- as.numeric(y); w <- if (!is.null(cc$weight_col) && cc$weight_col %in% names(d)) as.numeric(d[[cc$weight_col]]) else rep(1, nrow(d))
    w[!is.finite(w) | w <= 0] <- 1
    # CE_MALAWI_DISTRICT=1 (MW-01): analyse Malawi at its DISTRICT rung (our
    # "Admin1", 27 units of ~3 clusters) instead of Traditional Authorities,
    # which are mostly single clusters. Other countries unchanged.
    unit_col <- if (cn == "Malawi" && Sys.getenv("CE_MALAWI_DISTRICT", "0") == "1") (cc$admin1_col %||% "Admin1") else a2
    unit <- as.character(d[[unit_col]]); cl <- as.character(d[[clu]])
    ok <- is.finite(y) & !is.na(unit) & !is.na(cl); y <- y[ok]; w <- w[ok]; unit <- unit[ok]; cl <- cl[ok]
    ncl <- tapply(cl, unit, function(z) length(unique(z))); multi <- names(ncl)[ncl >= 2]
    keep <- unit %in% multi; if (sum(keep) < 40 || length(multi) < 8) { cat(sprintf("  %-12s %-12s skip: %d multi-cluster units\n", cn, on, length(multi))); next }
    y <- y[keep]; w <- w[keep]; unit <- unit[keep]; cl <- cl[keep]
    rw <- rc <- numeric(NSPLIT)
    for (s in seq_len(NSPLIT)) {
      # within: random halves of respondents inside each unit
      gw <- integer(length(y)); for (u in unique(unit)) { i <- which(unit == u); gw[i] <- sample(rep(1:2, length.out = length(i))) }
      # cluster: random halves of CLUSTERS inside each unit
      gc <- integer(length(y)); for (u in unique(unit)) { i <- which(unit == u); cls <- unique(cl[i]); lab <- sample(rep(1:2, length.out = length(cls)))
        gc[i] <- lab[match(cl[i], cls)] }
      rw[s] <- half_corr(y, w, unit, gw); rc[s] <- half_corr(y, w, unit, gc)
    }
    rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, n_units_multi = length(multi), n_resp = length(y),
      r_within = median(rw, na.rm = TRUE), r_cluster = median(rc, na.rm = TRUE),
      ceiling_within = sb(median(rw, na.rm = TRUE)), ceiling_cluster = sb(median(rc, na.rm = TRUE)), stringsAsFactors = FALSE)
    cat(sprintf("  %-12s %-12s multi-cluster units %2d | ceiling within %.2f | cluster %.2f\n", cn, on, length(multi),
                sb(median(rw, na.rm = TRUE)), sb(median(rc, na.rm = TRUE))))
  } }
R <- bind_rows(rows); if (!nrow(R)) stop("no cells produced a ceiling -- see skip messages above")
R$inflation <- round(R$ceiling_within - R$ceiling_cluster, 3)
write.csv(R, file.path(OUTDIR, if (Sys.getenv("CE_MALAWI_DISTRICT", "0") == "1") "ceiling_cluster_split_malawi_district.csv" else "ceiling_cluster_split.csv"), row.names = FALSE)
cat("\n===== CE-01: ceiling on multi-cluster units, within-split vs cluster-split =====\n")
print(as.data.frame(R |> mutate(across(c(r_within, r_cluster, ceiling_within, ceiling_cluster), ~ round(.x, 3)))), row.names = FALSE)
cat("\nby country (mean over outcomes):\n")
print(aggregate(cbind(ceiling_within, ceiling_cluster, inflation) ~ country, data = R, FUN = function(z) round(mean(z, na.rm = TRUE), 3)), row.names = FALSE)
cat(sprintf("\noverall: within %.3f | cluster %.3f | cluster-effect share of the within ceiling %.0f%%\n",
            mean(R$ceiling_within, na.rm = TRUE), mean(R$ceiling_cluster, na.rm = TRUE),
            100 * (1 - mean(R$ceiling_cluster, na.rm = TRUE) / mean(R$ceiling_within, na.rm = TRUE))))
cat("\nDONE\n")
