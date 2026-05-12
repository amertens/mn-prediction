# =============================================================================
# 02_gbd_placeholder.R
#
# Build placeholder GBD/IHME comparison data for the Burden of Disease tab.
# This is INTENTIONALLY filled with reasonable-looking but non-authoritative
# values. The dashboard surfaces a clear note that an RA should replace these
# with actual GBD Results Tool exports.
#
# To replace with real data, an RA should:
#   1. Visit https://vizhub.healthdata.org/gbd-results/
#   2. Download "Prevalence" for each of:
#       - Vitamin A deficiency
#       - Iron-deficiency anaemia
#       - (B12 and folate are not standalone GBD outcomes — leave NA)
#       - Zinc deficiency (not standalone in GBD — leave NA)
#     For ages: 6–59 months (children); 15–49 years (women)
#     For countries: GMB, GHA, SLE, MWI, CIV
#     For years: survey year and 2023
#   3. Save the CSV to dashboard/data-raw/gbd_estimates.csv
#   4. Re-run this script to convert to RDS
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

DASHBOARD_DATA <- here::here("dashboard", "data")

countries <- c("Gambia", "Ghana", "Sierra Leone", "Malawi")

# Placeholder GBD prevalences. These are illustrative values consistent with
# published Lancet / WHO regional estimates but are NOT from the GBD Results
# Tool. The dashboard surfaces this caveat prominently.
gbd_placeholder <- tibble::tribble(
  ~country,        ~outcome,        ~year, ~gbd_prev, ~gbd_lo, ~gbd_hi, ~status,
  "Gambia",        "child_vitA",     2018,    0.180,   0.140,   0.220, "placeholder",
  "Gambia",        "women_iron",     2018,    0.250,   0.200,   0.300, "placeholder",
  "Gambia",        "child_iron",     2018,    0.350,   0.290,   0.410, "placeholder",
  "Gambia",        "women_vitA",     2018,    0.040,   0.020,   0.060, "placeholder",
  "Ghana",         "child_vitA",     2017,    0.210,   0.170,   0.250, "placeholder",
  "Ghana",         "women_iron",     2017,    0.110,   0.080,   0.140, "placeholder",
  "Ghana",         "child_iron",     2017,    0.180,   0.140,   0.220, "placeholder",
  "Ghana",         "women_vitA",     2017,    0.020,   0.010,   0.040, "placeholder",
  "Sierra Leone",  "child_vitA",     2013,    0.220,   0.180,   0.270, "placeholder",
  "Sierra Leone",  "child_iron",     2013,    0.080,   0.040,   0.130, "placeholder",
  "Sierra Leone",  "women_vitA",     2013,    0.030,   0.010,   0.050, "placeholder",
  "Malawi",        "child_vitA",     2016,    0.210,   0.170,   0.260, "placeholder",
  "Malawi",        "women_iron",     2016,    0.090,   0.060,   0.130, "placeholder",
  "Malawi",        "child_iron",     2016,    0.140,   0.100,   0.180, "placeholder",
  "Malawi",        "women_vitA",     2016,    0.050,   0.030,   0.080, "placeholder",
  "Malawi",        "child_zinc",     2016,    0.620,   0.520,   0.720, "placeholder",
  # 2023 projections — flat assumption (replace with actual GBD trends)
  "Gambia",        "child_vitA",     2023,    0.150,   0.110,   0.190, "placeholder",
  "Gambia",        "women_iron",     2023,    0.230,   0.180,   0.280, "placeholder",
  "Gambia",        "child_iron",     2023,    0.330,   0.270,   0.390, "placeholder",
  "Ghana",         "child_vitA",     2023,    0.180,   0.140,   0.220, "placeholder",
  "Ghana",         "women_iron",     2023,    0.100,   0.070,   0.130, "placeholder",
  "Sierra Leone",  "child_vitA",     2023,    0.190,   0.150,   0.240, "placeholder",
  "Malawi",        "child_vitA",     2023,    0.180,   0.140,   0.230, "placeholder",
  "Malawi",        "child_zinc",     2023,    0.580,   0.480,   0.690, "placeholder"
)

# If the RA has populated dashboard/data-raw/gbd_estimates.csv, prefer that
real_gbd_path <- here::here("dashboard", "data-raw", "gbd_estimates.csv")

if (file.exists(real_gbd_path)) {
  cat("Loading real GBD data from gbd_estimates.csv\n")
  gbd <- read.csv(real_gbd_path, stringsAsFactors = FALSE)
  gbd$status <- "GBD_actual"
} else {
  cat("Using PLACEHOLDER GBD values — flag this in the dashboard.\n")
  gbd <- gbd_placeholder
}

# Metadata note for dashboard to surface
gbd_meta <- list(
  using_placeholder = !file.exists(real_gbd_path),
  source_url = "https://vizhub.healthdata.org/gbd-results/",
  ra_instructions = paste(
    "RA TASK: Download GBD prevalence estimates from the IHME GBD Results Tool",
    "(URL above). Save as dashboard/data-raw/gbd_estimates.csv with columns:",
    "country, outcome, year, gbd_prev, gbd_lo, gbd_hi.",
    "Then re-run dashboard/data-raw/02_gbd_placeholder.R to refresh.",
    "Outcome codes match dashboard internal codes (child_vitA, women_iron, etc.)."
  )
)

saveRDS(list(data = gbd, meta = gbd_meta),
        file.path(DASHBOARD_DATA, "gbd_estimates.rds"))

cat(sprintf("\nSaved %d GBD comparison rows (status: %s)\n",
            nrow(gbd),
            if (gbd_meta$using_placeholder) "PLACEHOLDER" else "actual GBD data"))
