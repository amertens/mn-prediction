# =============================================================================
# src/Tanzania/1_GW_Tanzania_data_clean.R
#
# Tanzania (TDHS 2009-10) vitamin A outcome construction -> Tanzania_GMS_cleaned.rds
# (consumed by 2_GW_Tanzania_data_merge.R). VALIDATED 2026-07-01.
#
# Pipeline: OB vitamin A (RBP) + PR blood-id crosswalk + GE cluster GPS.
#
# TWO issues found during validation and fixed here:
#   1. haven AND foreign SEGFAULT on the large PR recode (TZPR63FL.DTA). The
#      blood-id -> cluster/weight/demographics crosswalk is therefore built in
#      Python (pandas) by src/Tanzania/0_extract_pr_crosswalk.py -> pr_crosswalk.csv;
#      this script reads that CSV. (VITA reads fine in R with encoding="latin1"
#      — the default encoding also segfaults.)
#   2. RBP in the OB file is in mg/L (median ~25), NOT umol/L. The WHO VAD cutoff
#      is <0.70 umol/L, so RBP is converted mg/L -> umol/L (divide by the RBP4
#      molar mass ~21.2 kDa) before thresholding. Direct <0.70 on mg/L gave an
#      implausible 5% child VAD; the conversion gives ~25% (women ~8%), matching
#      expectations. CONFIRM the 21.2 factor / cutoff against TZ61BIOMARKER.DOC.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(haven); library(here); library(sf)
})

tz_dir <- here("data", "IPD", "Tanzania 2010")
RBP_MW <- 21.2   # RBP4 ~21.2 kDa: divisor to convert RBP mg/L -> umol/L
THRESH <- 0.70   # WHO VAD cutoff, umol/L

# ── 1. OB vitamin A (encoding="latin1" avoids a haven segfault on this file) ──
vita <- read_dta(file.path(tz_dir, "TZOB61DT", "TZ61VITA.DTA"), encoding = "latin1")
ob_na <- function(x) { x <- suppressWarnings(as.numeric(x))
                       ifelse(!is.na(x) & x >= 999.99, NA_real_, x) }
vita <- vita %>% transmute(
  blood_id  = toupper(trimws(as.character(lbar))),
  rbp_umol  = ob_na(rbpadcrp) / RBP_MW,   # CRP-adjusted RBP, mg/L -> umol/L
  rbp_raw_umol = ob_na(rbpres) / RBP_MW,  # unadjusted, for reference
  crp       = ob_na(crpres),
  stfr      = ob_na(stfrres)               # iron biomarker (not an outcome here)
)
cat(sprintf("OB VITA: %d samples, %d non-missing RBP\n",
            nrow(vita), sum(!is.na(vita$rbp_umol))))

# ── 2. PR blood-id crosswalk (built by 0_extract_pr_crosswalk.py) ────────────
xw_path <- file.path(tz_dir, "pr_crosswalk.csv")
if (!file.exists(xw_path))
  stop("Missing pr_crosswalk.csv — run: python src/Tanzania/0_extract_pr_crosswalk.py")
xw <- read.csv(xw_path, stringsAsFactors = FALSE)
xw$blood_id <- toupper(trimws(xw$blood_id))

# ── 3. Join biomarkers onto the crosswalk ────────────────────────────────────
df <- dplyr::inner_join(xw, vita, by = "blood_id")
cat(sprintf("OB<->PR match: %d of %d OB samples (%.1f%%)\n",
            nrow(df), nrow(vita), 100 * nrow(df) / nrow(vita)))
if (nrow(df) < 0.8 * nrow(vita))
  warning("Low OB<->PR match rate — check blood-id formatting.")

# ── 4. Derive vitamin A outcomes (umol/L; VAD < 0.70) ────────────────────────
df <- df %>% mutate(
  cRBPAdjCRP = ifelse(child_flag == 1L, rbp_umol, NA_real_),
  wRBPAdjCRP = ifelse(child_flag == 0L, rbp_umol, NA_real_),
  cVAD       = ifelse(child_flag == 1L, as.integer(rbp_umol < THRESH), NA_integer_),
  wVAD       = ifelse(child_flag == 0L, as.integer(rbp_umol < THRESH), NA_integer_)
)
w <- df$svy_weight
wprev <- function(y) { ok <- !is.na(y); 100 * sum(y[ok] * w[ok]) / sum(w[ok]) }
cat(sprintf("VAD (weighted) — children: %.1f%%   women: %.1f%%\n",
            wprev(df$cVAD), wprev(df$wVAD)))

# ── 5. Cluster GPS from the GE shapefile ─────────────────────────────────────
ge  <- sf::st_read(file.path(tz_dir, "TZGE61FL", "TZGE61FL.shp"), quiet = TRUE)
gps <- ge %>% sf::st_drop_geometry() %>%
  transmute(cnum = as.integer(DHSCLUST),
            latitude = as.numeric(LATNUM), longitude = as.numeric(LONGNUM)) %>%
  mutate(across(c(latitude, longitude), ~ ifelse(. == 0, NA_real_, .)))
df <- left_join(df, gps, by = "cnum")
cat(sprintf("GPS join: %d rows missing coordinates\n", sum(is.na(df$latitude))))

# ── 6. gw_-prefix all survey columns (except lat/long) and save ──────────────
# Keep rbp_raw_umol: brinda_rbp_cols("Tanzania") points the harmonized VitA
# outcome at the RAW RBP + CRP (R/brinda_adjustment.R). It is excluded from the
# predictor set by the config's gw_exclude_patterns ("rbp"), so it cannot leak.
keep_geo <- c("latitude", "longitude")
colnames(df)[!(colnames(df) %in% keep_geo)] <-
  paste0("gw_", colnames(df)[!(colnames(df) %in% keep_geo)])

expected <- c("gw_cnum", "gw_child_flag", "gw_svy_weight", "gw_strata",
              "gw_cRBPAdjCRP", "gw_wRBPAdjCRP", "gw_cVAD", "gw_wVAD",
              "gw_rbp_raw_umol", "gw_crp")
miss <- setdiff(expected, colnames(df))
if (length(miss)) warning("Missing config-expected columns: ", paste(miss, collapse = ", "))

out <- file.path(tz_dir, "Tanzania_GMS_cleaned.rds")
saveRDS(df, out)
cat(sprintf("\nSaved %s  (%d rows, %d cols)\n", out, nrow(df), ncol(df)))
