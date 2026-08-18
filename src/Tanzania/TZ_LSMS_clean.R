# =============================================================================
# src/Tanzania/TZ_LSMS_clean.R
#
# SCAFFOLD — build admin-1 (region) LSMS predictors for Tanzania from the
# Tanzania National Panel Survey (NPS / LSMS-ISA), analogous to
# src/Ghana/LSMS_clean.R. Output: data/LSMS/Tanzania_LSMS_clean.RDS with
# `lsms_` columns keyed on region (Admin1), joined in the merge on Admin1.
#
# LSMS is OPTIONAL and admin-1 (coarse) — only Ghana and Gambia use it among the
# existing countries. Prioritise after GEE/soil/IHME.
#
# -----------------------------------------------------------------------------
# BLOCKING: requires a registered download (cannot be auto-fetched).
#   Source : World Bank Microdata Library — Tanzania National Panel Survey
#            https://microdata.worldbank.org/index.php/catalog/lsms  (search NPS)
#   Wave   : NPS 2010-2011 (Wave 2) — best match for the TDHS 2009-10 survey year
#            (Wave 1 2008-09 is an acceptable alternative).
#   Access : free account + accept the data-use terms; download the Stata (.dta)
#            package. Place the files under data/LSMS/Tanzania_NPS_2010/.
#
# The NPS file/variable layout differs from Ghana's GLSS "g7aggregates" files,
# so the two placeholders below (CONSUMPTION FILE + REGION/WEIGHT vars) MUST be
# confirmed against the downloaded data before running. The aggregation logic
# (survey-weighted region means of numeric welfare vars) is the same as Ghana.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(janitor); library(survey); library(here)
})

nps_dir <- here("data", "LSMS", "Tanzania_NPS_2010")
if (!dir.exists(nps_dir))
  stop("Download the Tanzania NPS 2010-11 (World Bank Microdata) into ", nps_dir,
       " first — see the header of this script.")

# --- CONFIRM THESE against the downloaded NPS files -------------------------
# The NPS consumption/welfare aggregate file (household-level, one row per hh,
# with a region code and household weight). Typical NPS names include a
# consumption aggregate like "ConsumptionNPS2010.dta" or a welfare file with
# `hhweight`/`y2_weight` and a region variable `region` / `hh_a01_1`.
HH_FILE     <- file.path(nps_dir, "ConsumptionNPS2010.dta")  # <-- confirm filename
REGION_VAR  <- "region"      # <-- confirm admin-1 region code/name variable
WEIGHT_VAR  <- "hhweight"    # <-- confirm household weight variable
# Optional map of NPS numeric region codes -> region names (if REGION_VAR is a
# code). Leave NULL if REGION_VAR is already a name. Tanzania mainland + Zanzibar
# regions (~26 in 2010). Fill from the NPS codebook if needed.
REGION_LABELS <- NULL
# ---------------------------------------------------------------------------

if (!file.exists(HH_FILE))
  stop("Consumption/welfare file not found: ", HH_FILE,
       " — set HH_FILE to the actual NPS aggregate filename.")

hh <- read_dta(HH_FILE) |> clean_names()
stopifnot(REGION_VAR %in% names(hh), WEIGHT_VAR %in% names(hh))

# Numeric welfare/consumption variables to aggregate (exclude ids/weights/geo).
id_like <- c("hhid", "y2_hhid", "clusterid", "clust", "ea_id", REGION_VAR,
             WEIGHT_VAR, "district", "ward", "urban", "rural", "domain")
num_vars <- hh |> select(where(is.numeric)) |> select(-any_of(id_like)) |> names()

# Survey-weighted region means (same approach as src/Ghana/LSMS_clean.R).
des <- svydesign(ids = ~1, weights = as.formula(paste0("~", WEIGHT_VAR)), data = hh)
region_means <- svyby(
  as.formula(paste("~", paste(num_vars, collapse = "+"))),
  as.formula(paste0("~", REGION_VAR)),
  des, svymean, na.rm = TRUE, vartype = NULL)

# Drop near-zero-variance columns.
nzv <- caret::nzv(region_means)
if (length(nzv)) region_means <- region_means[, -nzv, drop = FALSE]

# Attach a clean Admin1 name column for the merge join.
if (!is.null(REGION_LABELS)) {
  region_means$Admin1 <- REGION_LABELS[as.character(region_means[[REGION_VAR]])]
} else {
  region_means$Admin1 <- trimws(as.character(region_means[[REGION_VAR]]))
}

# Prefix lsms_ (keep Admin1 unprefixed for the join key).
val_cols <- setdiff(names(region_means), "Admin1")
names(region_means)[match(val_cols, names(region_means))] <- paste0("lsms_", val_cols)

saveRDS(region_means, here("data", "LSMS", "Tanzania_LSMS_clean.RDS"))
cat(sprintf("Saved data/LSMS/Tanzania_LSMS_clean.RDS (%d regions x %d lsms_ cols)\n",
            nrow(region_means), sum(grepl("^lsms_", names(region_means)))))
