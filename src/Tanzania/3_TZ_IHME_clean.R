# =============================================================================
# src/Tanzania/3_TZ_IHME_clean.R
#
# Build Tanzania (TDHS 2009-10) IHME admin-2 and admin-1 predictor tables from
# the GLOBAL Local Burden of Disease (LBD) geospatial CSVs already in the repo
# (data/IHME/<family>/...). Tanzania is present in those global files, so NO new
# download is needed — this is a filter + reshape, run analogous to the other
# four countries.
#
# Outputs (consumed by src/Tanzania/2_GW_Tanzania_data_merge.R):
#   data/IHME/tanzania_2010_merged_IHME_data.csv          (admin-2, ihme_)
#   data/IHME/tanzania_2010_merged_IHME_admin1_data.csv   (admin-1, ihme_adm1_)
#
# Year: IHME geospatial models run 2000-2017/2018/2019. We extract YEAR = 2010
# (the survey year) to match how the existing countries pin to their survey year
# (Sierra Leone 2012, Ghana 2017, Malawi 2015, Gambia 2018).
#
# After running, check the PARITY line printed for the admin-2 build: it reports
# how many Tanzania ihme_ columns match Ghana's existing column names. High
# overlap (minus the known Malaria/education admin-2 gap) confirms Tanzania will
# align in the pooled / LOCO common-vocabulary set.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({ library(here) })

source(here("src", "IHME", "build_ihme_admin2.R"))
source(here("src", "IHME", "build_ihme_admin1.R"))

COUNTRY <- "Tanzania"   # confirm this matches ADM0_NAME in the LBD files
YEAR    <- 2010L

# ── Admin-2 ────────────────────────────────────────────────────────────────
build_ihme_admin2(
  country_name  = COUNTRY,
  year_filter   = YEAR,
  out_path      = here("data", "IHME", "tanzania_2010_merged_IHME_data.csv"),
  reference_csv = here("data", "IHME", "ghana_2017_merged_IHME_data.csv")
)

# ── Admin-1 ────────────────────────────────────────────────────────────────
build_ihme_admin1(
  country_name = COUNTRY,
  year_filter  = YEAR,
  out_path     = here("data", "IHME", "tanzania_2010_merged_IHME_admin1_data.csv")
)

cat("\nDone. Tanzania IHME admin-1/2 tables written to data/IHME/.\n")
