# =============================================================================
# scripts/covariates/build_mns_overlap_clusters.R
#
# Write metadata/mns_dhs_overlap_clusters.csv: the DHS clusters that a
# country's micronutrient survey re-sampled, one row per country x cluster.
#
# Malawi qualifies for DHS and, since 2026-09-15, The Gambia for MICS (the
# GMNS 2018 was fielded inside the MICS6 2018 sample). Malawi: the MNS 2015-16 was a subsample of the MDHS 2015-16
# (105 of 850 clusters; README: "use the cluster, household and line number
# (MCLUSTER, MNUMBER, M01) to merge with DHS data"). `gw_cnum` in the merged
# dataset is MCLUSTER. The Gambia (GMNS 2018 vs DHS 2019-20), Ghana (GMS 2017
# vs DHS 2014) and Sierra Leone (SLMS 2013 vs DHS 2013) surveys drew their own
# samples, so they contribute no rows.
#
# The list is checked against the DHS GPS file of the same round before it is
# written: every MNS cluster must be a cluster of MWGE7AFL, otherwise the
# merged dataset and the recodes are not the same round and nothing downstream
# should trust the exclusion.
#
#   Rscript -e "source('scripts/covariates/build_mns_overlap_clusters.R')"
# -> metadata/mns_dhs_overlap_clusters.csv
# =============================================================================
suppressPackageStartupMessages(library(here))

OUT   <- here("metadata", "mns_dhs_overlap_clusters.csv")
CACHE <- "C:/Users/andre/AppData/Local/andre/rdhs/Cache/datasets"

mns <- readRDS(here("data", "IPD", "Malawi", "Malawi_merged_dataset.rds"))
stopifnot("gw_cnum" %in% names(mns))
cl <- sort(unique(as.integer(unclass(mns$gw_cnum))))
cl <- cl[!is.na(cl)]
cat(sprintf("[mns-overlap] Malawi merged dataset: %d rows, %d distinct MNS clusters\n",
            nrow(mns), length(cl)))

ge_path <- file.path(CACHE, "MWGE7AFL.rds")
if (file.exists(ge_path)) {
  ge <- readRDS(ge_path)
  gecl <- suppressWarnings(as.integer(as.numeric(unclass(ge$DHSCLUST))))
  miss <- setdiff(cl, gecl)
  cat(sprintf("[mns-overlap] DHS GPS file MWGE7AFL: %d clusters; MNS clusters not in it: %d\n",
              length(unique(gecl)), length(miss)))
  if (length(miss)) stop("MNS clusters absent from MWGE7AFL: ", paste(miss, collapse = ", "))
} else {
  cat("[mns-overlap] WARNING: MWGE7AFL not in the rdhs cache; GPS-round check skipped\n")
}

out <- data.frame(
  country     = "Malawi", programme = "DHS",
  dhs_survey  = "MW2015DHS (MWIR7ADT/MWKR7ADT/MWPR7ADT/MWHR7ADT/MWBR7ADT; GPS MWGE7AFL)",
  dhs_cluster = cl,
  source      = "unique(gw_cnum) of data/IPD/Malawi/Malawi_merged_dataset.rds; MNS 2015-16 sampled 105 of the 850 MDHS 2015-16 clusters (WSC4: 3,097 of 3,099 rows link to a DHS person)",
  stringsAsFactors = FALSE)

# The Gambia (2026-09-15): the GMNS 2018 was fielded inside the MICS6 2018
# sample. gw_MICS_Cluster_Number in the merged dataset is the MICS HH1; the
# respondents sit in 70 of the 390 MICS clusters and carry the MICS household
# number. Checked against hh.sav before writing.
gm <- readRDS(here("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"))
gm <- sf::st_drop_geometry(gm)
stopifnot("gw_MICS_Cluster_Number" %in% names(gm))
gcl <- sort(unique(suppressWarnings(as.integer(as.numeric(unclass(gm$gw_MICS_Cluster_Number)))))); gcl <- gcl[!is.na(gcl)]
mics_hh <- here("data", "MICS", "Gambia 2018", "hh.sav")
if (file.exists(mics_hh)) {
  h1 <- unique(as.integer(haven::read_sav(mics_hh, col_select = "HH1")$HH1))
  miss <- setdiff(gcl, h1)
  cat(sprintf("[mns-overlap] Gambia MICS6 2018 hh.sav: %d clusters; GMNS clusters not in it: %d\n", length(h1), length(miss)))
  if (length(miss)) stop("GMNS MICS cluster numbers absent from hh.sav: ", paste(miss, collapse = ", "))
} else cat("[mns-overlap] WARNING: data/MICS/Gambia 2018/hh.sav absent; MICS cluster check skipped\n")
cat(sprintf("[mns-overlap] Gambia merged dataset: %d rows, %d distinct MICS clusters\n", nrow(gm), length(gcl)))
out <- rbind(out, data.frame(
  country     = "Gambia", programme = "MICS",
  dhs_survey  = "GMB 2018 MICS6 (hh/ch/wm.sav HH1; GPS GambiaMICS2018GPS)",
  dhs_cluster = gcl,
  source      = "unique(gw_MICS_Cluster_Number) of data/IPD/Gambia/Gambia_merged_dataset.rds; GMNS 2018 was fielded in 70 of the 390 MICS6 2018 clusters and carries the MICS household number",
  stringsAsFactors = FALSE))
utils::write.csv(out, OUT, row.names = FALSE)
cat(sprintf("[mns-overlap] wrote %d rows -> %s\n", nrow(out), OUT))
