# Research Assistant: Data Acquisition Instructions

## Micronutrient Prediction Project — External Predictor Data

**Last updated:** March 2026

This document provides step-by-step instructions for obtaining API credentials and downloading datasets needed to expand the predictor set for the micronutrient deficiency prediction pipeline. Each section is independent — complete them in any order.

**Countries:** Gambia, Ghana, Sierra Leone, Malawi, Côte d'Ivoire (out-of-sample)
**Key survey years:** Gambia 2018, Ghana 2017, Sierra Leone 2012, Malawi 2016

---

## PART A: API Credentials (set up accounts)

### A1. NASA Earthdata (for VIIRS Nighttime Lights)

**What it provides:** Satellite-derived nighttime light radiance — a proxy for economic activity, urbanization, and infrastructure at ~500m resolution.

**Steps:**
1. Go to https://urs.earthdata.nasa.gov/users/new
2. Create a free account (any institutional affiliation works)
3. After account creation, go to https://urs.earthdata.nasa.gov/profile
4. Click "Generate Token" under "Bearer Token"
5. Copy the token string

**Deliverable:** Send me the bearer token string. I'll set it as:
```
Sys.setenv(NASA_EARTHDATA_TOKEN = "your_token_here")
```

**Time estimate:** 5 minutes

---

### A2. ACLED — Armed Conflict Location & Event Data

**What it provides:** Geolocated conflict events (battles, protests, violence against civilians) with dates, fatalities, and admin-2 location. Relevant because conflict disrupts food systems and health services.

**Steps:**
1. Go to https://developer.acleddata.com/
2. Click "Register" — select "Researcher/Academic" as user type
3. Fill in the registration form:
   - Organization: [your university]
   - Purpose: "Academic research on nutrition and health determinants in West Africa"
   - Data use: "Integrating conflict indicators as predictors in a machine learning model for micronutrient deficiency prevalence estimation"
4. You'll receive an email with your **access key**
5. Your registered **email address** is also needed

**Deliverable:** Send me:
- Your registered email address
- Your ACLED access key

I'll set them as:
```
Sys.setenv(ACLED_EMAIL_ADDRESS = "your.email@university.edu")
Sys.setenv(ACLED_ACCESS_KEY = "your_access_key")
```

**Time estimate:** 10 minutes (approval may take up to 24 hours)

---

### A3. IPC API — Food Security Phase Classifications

**What it provides:** Integrated Food Security Phase Classification (IPC) and Cadre Harmonisé data at admin-2 level — population counts in each food crisis phase (1=minimal through 5=famine).

**Steps:**
1. Go to https://docs.api.ipcinfo.org
2. Look for "Get API Key" or "Register" — the IPC API may require filling out a request form
3. If there's no self-service registration:
   - Email: ipc@ipcinfo.org
   - Subject: "API access request for academic research"
   - Body: "We are researchers at [institution] working on predicting micronutrient deficiency prevalence in West Africa using machine learning. We would like API access to download IPC/CH phase classification data for Gambia, Ghana, Sierra Leone, and Malawi at the subnational level. The data will be used as predictor variables in our models. Thank you."
4. You should receive an API key

**Deliverable:** Send me the API key. I'll set it as:
```
Sys.setenv(IPC_API_KEY = "your_key")
```

**Time estimate:** 10 minutes to request; 1-5 business days for approval

**Fallback:** If the API key is not obtainable, we already have Cadre Harmonisé data from the Excel file (`data/raw/CadreHarmonise/`). The API would give us cleaner, more complete data but is not strictly necessary.

---

### A4. DHS Program Account (already have — verify spatial data access)

**What it provides:** Modeled raster surfaces of DHS indicators (stunting, anemia, ITN use, etc.) at 5×5 km resolution. Much better spatial granularity than our current Admin-1 aggregated DHS variables.

**Steps:**
1. Log in at https://dhsprogram.com/ (we already have an account)
2. Go to https://spatialdata.dhsprogram.com/modeled-surfaces/
3. Verify you can access the download page
4. If prompted, you may need to add a "project" that includes spatial data access:
   - Go to https://dhsprogram.com/data/new-user-registration.cfm
   - Add a project requesting "Geographic Datasets" for Gambia, Ghana, Sierra Leone, and Malawi

**Deliverable:** Confirm that you can access the modeled surfaces download page for all 4 countries. The actual downloads are in Part B below.

**Time estimate:** 5 minutes (if account already has access); up to 48 hours if a new data request needs approval

---

## PART B: Direct Downloads (manual file downloads)

### B1. DHS Modeled Surfaces (5×5 km rasters)

**What it provides:** Spatially interpolated surfaces of key health indicators at 5km resolution. These are modeled from DHS microdata using Bayesian geostatistics. Each raster package includes a mean estimate surface and an uncertainty surface.

**Priority: HIGH — this directly addresses our biggest data limitation (DHS indicators currently at Admin-1 only)**

**Steps:**
1. Go to https://spatialdata.dhsprogram.com/modeled-surfaces/
2. For EACH of the following countries, download ALL available indicator packages:

   **Gambia:**
   - Select "Gambia" from the country dropdown
   - Check all available indicators (likely: stunting, wasting, anemia in children, anemia in women, ITN use, vaccination, improved water, improved sanitation, education)
   - Click "Download" — you'll get a .zip file per indicator

   **Ghana:**
   - Select "Ghana" — download all indicators
   - Ghana likely has the most indicators available (multiple DHS rounds: 2008, 2014, 2022)
   - Download the most recent available round for each indicator

   **Sierra Leone:**
   - Select "Sierra Leone" — download all indicators

   **Malawi:**
   - Select "Malawi" — download all indicators
   - Malawi 2015-16 DHS should have modeled surfaces; check for 2024 DHS too

3. Unzip each downloaded package
4. Place the raster files (.tif) in the following directory structure:
```
data/DHS_modeled_surfaces/
├── Gambia/
│   ├── GMGMS2018_stunting_mean.tif
│   ├── GMGMS2018_stunting_uncertainty.tif
│   ├── GMGMS2018_anemia_children_mean.tif
│   └── ...
├── Ghana/
│   ├── GHGMS2017_stunting_mean.tif
│   └── ...
├── SierraLeone/
│   └── ...
└── Malawi/
    └── ...
```

**Note:** The exact filenames will vary — just maintain the country subfolder structure. Include BOTH mean and uncertainty .tif files.

**Deliverable:** The filled `data/DHS_modeled_surfaces/` directory with all downloaded rasters.

**What to record:** For each country, note which indicators were available and which DHS round they came from. Create a simple text file `data/DHS_modeled_surfaces/README.txt` listing what you downloaded.

**Time estimate:** 30-60 minutes (depends on number of indicators per country)

---

### B2. HarvestStat Africa — Subnational Crop Statistics

**What it provides:** Subnational crop production, harvested area, and yields for 33 Sub-Saharan African countries (1980-2022). Based on FEWS NET agricultural statistics. Relevant because crop yields proxy for food availability and dietary quality.

**Steps:**
1. Go to https://datadryad.org/dataset/doi:10.5061/dryad.vq83bk42w
2. Download the main CSV file (the harmonized subnational crop statistics)
3. Also download the GeoPackage file (administrative boundaries used by FEWS NET)
4. Place files in:
```
data/external_cache/harveststat_africa.csv
data/external_cache/harveststat_boundaries.gpkg    (if available)
```

**Deliverable:** The CSV file placed in the correct location.

**Time estimate:** 5 minutes

---

### B3. Global Data Lab — Subnational Human Development Index

**What it provides:** Subnational HDI (education, health, income indices) at Admin-1 level for 161 countries, 1990-2021. Provides composite development indicators that capture multiple dimensions of wellbeing.

**Steps:**
1. Go to https://globaldatalab.org/shdi/
2. Click on "Table" or "Download" (look for download options)
3. Download the full SHDI dataset as CSV
   - If prompted for parameters: select "All countries", "All years", "All indicators"
   - The indicators needed are: SHDI, Health Index, Education Index, Income Index, Expected Years of Schooling, Mean Years of Schooling, Life Expectancy, GNI per capita
4. Place the file in:
```
data/external_cache/gdl_shdi_full.csv
```

**Alternative download path:**
- Go to https://globaldatalab.org/shdi/archive/
- Download the latest version of the full dataset

**Deliverable:** The CSV file in the correct location.

**Time estimate:** 5-10 minutes

---

### B4. Côte d'Ivoire GADM Boundaries (for out-of-sample prediction)

**What it provides:** Administrative boundary shapefiles for Côte d'Ivoire, needed for out-of-sample prediction mapping.

**Steps:**
This should download automatically when we run the pipeline, but in case of network issues:
1. Go to https://gadm.org/download_country.html
2. Select "Côte d'Ivoire"
3. Download the "R (sp)" level 1 and level 2 files
4. Place in `data/gadm/gadm/`

**Deliverable:** Should already be cached. Only needed if internet access is limited during pipeline runs.

**Time estimate:** 2 minutes

---

## PART C: R Package Installation

Before running the data extraction, the following R packages need to be installed. Run in R:

```r
# CRAN packages (straightforward)
install.packages(c(
  "ripc",           # IPC/CH API client
  "acled.api",      # ACLED conflict data
  "chirps",         # CHIRPS rainfall (optional — we also download directly)
  "malariaAtlas",   # Malaria Atlas Project rasters
  "blackmarbler",   # VIIRS nighttime lights
  "rdhs"            # DHS data API
))

# GitHub packages (if not on CRAN)
# WorldPop download tools
devtools::install_github("wpgp/wpgpDownloadR")

# Check installations
for (pkg in c("ripc", "acled.api", "chirps", "malariaAtlas",
              "blackmarbler", "rdhs", "wpgpDownloadR")) {
  cat(sprintf("  %s: %s\n", pkg,
      if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"))
}
```

**Deliverable:** All packages installed without errors.

**Time estimate:** 10-15 minutes

---

## PART D: Verification Checklist

Once all credentials and downloads are in place, run this verification script in R from the project root:

```r
# Set all credentials
Sys.setenv(NASA_EARTHDATA_TOKEN = "your_token")
Sys.setenv(ACLED_EMAIL_ADDRESS = "your_email")
Sys.setenv(ACLED_ACCESS_KEY = "your_key")
Sys.setenv(IPC_API_KEY = "your_key")  # if obtained

# Verify credentials
cat("=== Credential Check ===\n")
cat("NASA token:", ifelse(nchar(Sys.getenv("NASA_EARTHDATA_TOKEN")) > 0, "SET", "MISSING"), "\n")
cat("ACLED email:", ifelse(nchar(Sys.getenv("ACLED_EMAIL_ADDRESS")) > 0, "SET", "MISSING"), "\n")
cat("ACLED key:", ifelse(nchar(Sys.getenv("ACLED_ACCESS_KEY")) > 0, "SET", "MISSING"), "\n")
cat("IPC key:", ifelse(nchar(Sys.getenv("IPC_API_KEY")) > 0, "SET", "MISSING (optional)"), "\n")

# Verify downloads
cat("\n=== File Check ===\n")
files_to_check <- c(
  "data/external_cache/harveststat_africa.csv",
  "data/external_cache/gdl_shdi_full.csv",
  "data/DHS_modeled_surfaces/Gambia/",
  "data/DHS_modeled_surfaces/Ghana/",
  "data/DHS_modeled_surfaces/SierraLeone/",
  "data/DHS_modeled_surfaces/Malawi/"
)
for (f in files_to_check) {
  exists <- file.exists(here::here(f)) || dir.exists(here::here(f))
  cat(sprintf("  %s: %s\n", f, ifelse(exists, "FOUND", "MISSING")))
}

# Verify packages
cat("\n=== Package Check ===\n")
for (pkg in c("ripc", "acled.api", "chirps", "malariaAtlas",
              "blackmarbler", "rdhs", "wpgpDownloadR")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %s: %s\n", pkg, ifelse(ok, "OK", "MISSING")))
}

# Quick test: try one extraction
cat("\n=== Quick Test ===\n")
source(here::here("R", "external_data.R"))
source(here::here("R", "config.R"))
cc <- get_country_configs()[["Ghana"]]
# Just test CHIRPS (no credentials needed)
admin2_sf <- sf::st_as_sf(geodata::gadm("GHA", level = 2,
                           path = here::here("data", "gadm")))
admin2_sf$Admin2 <- admin2_sf$NAME_2
cat("Testing CHIRPS extraction for Ghana (may take 2-3 min)...\n")
chirps_test <- tryCatch(
  extract_chirps(admin2_sf, "GHA", 2017, here::here("data", "external_cache")),
  error = function(e) { cat("  FAILED:", e$message, "\n"); NULL }
)
if (!is.null(chirps_test)) cat("  CHIRPS test: OK,", ncol(chirps_test) - 1, "variables\n")
```

---

## Summary: Priority Order

| Priority | Task | Time | Blocking? |
|----------|------|------|-----------|
| 1 | **B1: DHS Modeled Surfaces** | 30-60 min | Yes — highest impact new data |
| 2 | **A1: NASA Earthdata token** | 5 min | Yes — needed for nightlights |
| 3 | **A2: ACLED registration** | 10 min + wait | No — pipeline runs without it |
| 4 | **B2: HarvestStat CSV** | 5 min | No — pipeline runs without it |
| 5 | **B3: GDL CSV** | 5-10 min | No — pipeline runs without it |
| 6 | **C: R packages** | 10-15 min | Yes — needed before any extraction |
| 7 | **A3: IPC API key** | 10 min + wait | No — we have CH data already |
| 8 | **A4: DHS spatial access** | 5 min | Only if B1 requires new project |

**Start with C (packages), then B1 (DHS surfaces) and A1 (NASA), as these are highest impact and fastest.**

---

## Questions?

Contact Andrew Mertens for any questions about:
- Where to place downloaded files
- Which indicators to prioritize if not all are available
- Whether a specific dataset is relevant
