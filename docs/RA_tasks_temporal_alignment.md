# RA Tasks: Temporal Alignment of External Predictors

## Context

The micronutrient prediction pipeline links biomarker survey data to external predictors (satellite imagery, DHS indicators, climate data, etc.). For the predictors to be valid proxies, they should reflect conditions at or shortly before the survey fieldwork period. This document lists tasks to improve temporal alignment, ordered by priority.

### Survey Fieldwork Dates (reference)

| Country | Survey | Fieldwork | Pipeline Year |
|---------|--------|-----------|---------------|
| Gambia | GMNS | 13 Mar -- 4 May 2018 | 2018 |
| Ghana | GMS | 27 Apr -- 9 Jun 2017 | 2017 |
| Sierra Leone | SLMS | 11 Nov -- 2 Dec 2013 | 2013 |
| Malawi | MNS | 8 Dec 2015 -- 16 Feb 2016 | 2016 |

---

## Task 1: Generate Ghana DHS 2017 Admin-2 Indicators (HIGH PRIORITY)

**Problem:** The pipeline currently uses Ghana DHS 2014 indicators for the 2017 micronutrient survey --- a 3-year gap. Ghana DHS 2017 exists and the processing script is already configured for it.

**What to do:**

1. Ensure `rdhs` is configured with valid DHS Program credentials:
   ```r
   rdhs::set_rdhs_config(email = "your_email@example.com",
                          project = "your_project_name",
                          config_path = "~/.rdhs.json")
   ```

2. Run the DHS Admin-2 aggregation script (it will download Ghana 2017 microdata and produce smoothed Admin-2 estimates):
   ```r
   # In src/DHS/DHS_admin2_aggregation.R, line 56:
   # Set COUNTRIES_TO_RUN to only process Ghana 2017:
   COUNTRIES_TO_RUN <- c("Ghana")
   # Then source the script:
   source("src/DHS/DHS_admin2_aggregation.R")
   ```

3. Run the custom indicators script for the same year:
   ```r
   source("src/DHS/DHS_custom_admin2_indicators.R")
   ```

4. Verify output files were created in `data/DHS/clean/`:
   - `Ghana_2017_dhs_admin2_wide.rds`
   - `Ghana_2017_dhs_custom_admin2_wide.rds`
   - `Ghana_2017_dhs_admin2_fh_bym2.rds`
   - `Ghana_2017_dhs_admin1_direct.rds`

5. After these files exist, update the pipeline config:
   - In `R/config.R` line ~120: change `dhs_year = 2014L` to `dhs_year = 2017L`
   - In `R/config.R` line ~199: can optionally remove DHS2014 and DHS2016 from Ghana domains, keeping only `DHS2017 = list(prefix = "dhs2017_")`

**Requirements:** R with `rdhs`, `surveyPrev`, `SUMMER`, `INLA` packages installed. DHS Program account with approved data access for Ghana.

**Expected impact:** Eliminates a 3-year temporal gap for ~110 DHS indicators in Ghana models.

---

## Task 2: Export NDVI Rasters for All Countries (HIGH PRIORITY)

**Problem:** NDVI (Normalized Difference Vegetation Index) rasters currently exist only for Ghana (2010--2015) and are missing for Gambia, Sierra Leone, and Malawi entirely. NDVI is a strong proxy for agricultural productivity and food availability.

**What to do:**

Export annual NDVI composites from Google Earth Engine for each country at the survey year and the year prior. Use MODIS MOD13A1 (500m, 16-day) or MOD13A2 (1km, 16-day).

**GEE export parameters per country:**

| Country | Years to export | GEE collection | Filename pattern |
|---------|----------------|----------------|-----------------|
| Gambia | 2017, 2018 | MODIS/006/MOD13A1 | `NDVI_Gambia_{year}.tif` |
| Ghana | 2016, 2017 | MODIS/006/MOD13A1 | `NDVI_Ghana_{year}.tif` |
| Sierra Leone | 2012, 2013 | MODIS/006/MOD13A1 | `NDVI_Sierra_Leone_{year}.tif` |
| Malawi | 2015, 2016 | MODIS/006/MOD13A1 | `NDVI_Malawi_{year}.tif` |

**GEE script template:**
```javascript
var country = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  .filter(ee.Filter.eq('country_na', 'Gambia'));  // change per country
var year = 2018;  // change per country-year

var ndvi = ee.ImageCollection("MODIS/006/MOD13A1")
  .filterDate(year + '-01-01', year + '-12-31')
  .filterBounds(country)
  .select('NDVI')
  .mean()
  .clip(country);

Export.image.toDrive({
  image: ndvi,
  description: 'NDVI_Gambia_' + year,
  folder: 'GEE_exports',
  region: country.geometry(),
  scale: 500,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
```

**Save files to:** `data/{Country}_GEE_rasters/NDVI_{Country}_{year}.tif`

The pipeline will automatically pick them up --- `R/admin2_analysis.R` reads all `.tif` files from the raster directory.

---

## Task 3: Export Updated LST Night Rasters for Gambia and Ghana (MEDIUM PRIORITY)

**Problem:** Land Surface Temperature (nighttime) rasters exist only for 2014--2015 across all countries. For Gambia (2018 survey) this is a 3-year gap; for Ghana (2017) a 2-year gap.

**What to do:**

Export MODIS LST nighttime annual composites from GEE:

| Country | Years to export | Current years available |
|---------|----------------|----------------------|
| Gambia | 2017, 2018 | 2014, 2015 |
| Ghana | 2016, 2017 | 2014, 2015 |

Sierra Leone (2013) and Malawi (2016) already have adequate coverage with 2012--2015 rasters.

**GEE collection:** `MODIS/006/MOD11A2` (1km, 8-day)
**Band:** `LST_Night_1km`
**Aggregation:** Annual mean

**GEE script template:**
```javascript
var country = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  .filter(ee.Filter.eq('country_na', 'Gambia'));
var year = 2018;

var lst = ee.ImageCollection("MODIS/006/MOD11A2")
  .filterDate(year + '-01-01', year + '-12-31')
  .filterBounds(country)
  .select('LST_Night_1km')
  .mean()
  .multiply(0.02).subtract(273.15)  // Scale to Celsius
  .clip(country);

Export.image.toDrive({
  image: lst,
  description: 'LST_Night_' + year + '_Annual_Mean_Gambia',
  folder: 'GEE_exports',
  region: country.geometry(),
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
```

**Save files to:** `data/{Country}_GEE_rasters/LST_Night_{year}_Annual_Mean_{Country}.tif`

---

## Task 4: Export FLDAS 2016 Raster for Malawi (MEDIUM PRIORITY)

**Problem:** Malawi FLDAS rasters cover 2014--2015, but the survey spans Dec 2015 -- Feb 2016. The 2016 FLDAS would better capture conditions during the majority of fieldwork.

**What to do:**

Export FLDAS 2016 annual mean from GEE:

**GEE collection:** `NASA/FLDAS/NOAH01/C/GL/M/V001`
**Variables:** All bands (soil moisture, evapotranspiration, etc.)
**Aggregation:** Annual mean across 12 monthly composites

**Save to:** `data/Malawi_GEE_rasters/FLDAS_2016_Annual_Mean_Malawi.tif`

---

## Task 5: Investigate ESA Land Cover Alternatives to WorldCereal (LOW PRIORITY)

**Problem:** ESA WorldCereal 2021 is used for all countries but dates from 3--8 years after the surveys. The ESA CCI Land Cover product provides annual maps from 1992--2020 at 300m resolution.

**What to do:**

1. Check ESA CCI Land Cover availability on GEE: `ESA/WorldCover/v100` (2020 only, 10m) or `users/...` community uploads
2. If annual products are available for 2013--2018, export country-specific rasters at the survey year
3. Alternative: use MODIS Land Cover (`MODIS/006/MCD12Q1`) which has annual global coverage 2001--present at 500m

**Save to:** Replace or supplement `ESA_WorldCereal_2021_{Country}.tif` in each country's raster directory

---

## Task 6: Verify GPW Demographic 2010 Necessity (LOW PRIORITY)

**Problem:** `GPW_Demographic_2010_*.tif` provides 2010 census-based population demographics for all countries. This is 3--8 years behind the surveys. WorldPop (already in the pipeline) provides annually-updated population estimates.

**What to do:**

1. Check what variables GPW Demographic provides beyond total population (age structure? sex ratio? urban/rural?)
2. If it only provides total population, it's redundant with WorldPop and can be dropped
3. If it provides unique demographic breakdowns (age/sex), check whether UN World Population Prospects or WorldPop age-sex structures have more recent versions
4. Document findings in a short memo

---

## Task 7: Improve WFP Food Price Admin Matching for Ghana and Malawi (LOW PRIORITY)

**Problem:** The newly integrated WFP food price data matched well for Gambia (11/37 Admin-2 units) and Sierra Leone (13/14), but poorly for Ghana (0/260) and Malawi (3/326). This is because Ghana and Malawi WFP data uses admin1-level market locations that don't match the GADM Admin-2 district names.

**What to do:**

1. Open `data/food_price/wfp_food_prices_gha.csv` and check the `admin1` and `admin2` columns
2. Compare the WFP admin names to GADM Admin-1 names for Ghana and Malawi
3. Create a manual crosswalk CSV mapping WFP market names to GADM Admin-1 or Admin-2 names
4. Alternatively, modify the extraction to broadcast Admin-1 level prices to all Admin-2 districts within that region (same approach used for GDL HDI data)

---

## Summary Priority Matrix

| Task | Priority | Effort | Impact on Models |
|------|----------|--------|-----------------|
| 1. Ghana DHS 2017 | **HIGH** | 2--3 hours | Eliminates 3-year gap for ~110 indicators |
| 2. NDVI for all countries | **HIGH** | 1--2 hours (GEE export) | Adds vegetation proxy to 3 countries |
| 3. LST Night 2017/2018 | MEDIUM | 1 hour (GEE export) | Closes 2--3 year gap for 2 countries |
| 4. FLDAS 2016 for Malawi | MEDIUM | 30 min (GEE export) | Closes 1-year gap for 1 country |
| 5. ESA Land Cover alternatives | LOW | 2 hours (research + export) | Replaces 3--8 year stale cropland map |
| 6. GPW Demographic review | LOW | 1 hour (research) | May drop redundant predictor |
| 7. WFP price matching fix | LOW | 2 hours | Improves food price coverage for 2 countries |
