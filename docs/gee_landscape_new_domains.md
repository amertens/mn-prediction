# GEE Catalog Landscape: Additional Predictor Domains for Micronutrient-Deficiency Prediction

**Scope:** Sub-Saharan Africa (Gambia, Ghana, Sierra Leone, Malawi). Outcomes: Admin-2 prevalence of iron, vitamin A, folate, B12, and zinc deficiency.

**Purpose:** Identify Google Earth Engine (GEE) public-catalog datasets covering predictor domains plausibly linked to micronutrient status that are **NOT** already in the current pipeline.

**Already in the pipeline (excluded here):** CHIRPS rainfall, temperature, NDVI/vegetation, land cover, VIIRS nighttime lights, WorldPop population density, Malaria Atlas Project, SoilGrids soil properties, generic cropland/agriculture.

---

## Prioritized table

Priority reflects (a) mechanistic plausibility for micronutrient status, (b) coverage/quality over the four study countries, and (c) likely added information beyond what current domains already capture. Resolution and temporal coverage are as of the latest catalog versions in mid-2026; verify in-catalog before ingest.

| # | Domain | Dataset / name | GEE asset ID | Resolution | Temporal coverage | Rationale (link to micronutrient status) |
|---|--------|----------------|--------------|------------|-------------------|------------------------------------------|
| 1 | Market access / remoteness | Accessibility to Cities 2015 (Weiss et al. 2018, Oxford MAP) | `Oxford/MAP/accessibility_to_cities_2015_v1_0` | ~1 km | 2015 (single epoch) | Travel time to nearest city proxies access to diverse markets, fortified foods, animal-source foods, and health services — a strong correlate of dietary quality and deficiency in rural SSA. |
| 2 | Drought / climatic water balance | TerraClimate (U. Idaho) | `IDAHO_EPSCOR/TERRACLIMATE` | ~4 km | 1958–present, monthly | PDSI, climate water deficit, actual ET, soil moisture, runoff, vapor-pressure deficit. Water stress drives crop failure and dietary-diversity loss; multi-band drought signal richer than rainfall alone. |
| 3 | Soil moisture (root zone) | SMAP L4 Global 3-hourly surface & root-zone soil moisture | `NASA/SMAP/SPL4SMGP/008` | 9 km | 2015-03–present, 3-hourly | Root-zone moisture governs crop growth, rainfed yields, and pasture productivity — upstream determinants of food availability and micronutrient intake. Complements rainfall (CHIRPS) with actual water available to plants. |
| 4 | Evapotranspiration / productivity | MODIS Terra Net ET 8-Day (gap-filled) | `MODIS/061/MOD16A2GF` (or `MODIS/061/MOD16A2`) | 500 m | 2001–present (GF), 8-day | ET reflects vegetation water use and primary productivity; low/declining ET flags agricultural stress affecting local food production and diet quality. |
| 5 | Surface water | JRC Global Surface Water (occurrence, seasonality, recurrence) | `JRC/GSW1_4/GlobalSurfaceWater` | 30 m | 1984–2021 | Permanent/seasonal water presence proxies irrigation potential, fishing access (a key dietary source of zinc, iron, B12, and vitamin A), and dry-season water security. |
| 6 | Air pollution (PM2.5) | Global satellite-derived PM2.5 (ACAG / van Donkelaar, V6.GL.02) | `projects/sat-io/open-datasets/GLOBAL-SATELLITE-PM25/ANNUAL` (and `/MONTHLY`) | ~1 km (0.01°) | 2000–2022, annual & monthly | Chronic PM2.5 exposure (incl. household biomass smoke) is linked to systemic inflammation that elevates ferritin/hepcidin and suppresses iron absorption — confounds iron-status interpretation and may track deficiency risk. |
| 7 | Livestock density | Gridded Livestock of the World (GLW3/GLW4) — cattle, goats, sheep, poultry, pigs | FAO/Harvard Dataverse rasters (ingest as assets; not a native `ee` ID). GEE-native annual variant: AGLW (community). | ~10 km (5 arc-min) | GLW3=2010, GLW4=2015/2020 | Livestock density proxies availability of animal-source foods (meat, milk, eggs, offal) — the principal dietary sources of bioavailable iron, zinc, B12, and preformed vitamin A. |
| 8 | Crop system / irrigation | GFSAD1000 Cropland Extent & irrigated-vs-rainfed mask | `USGS/GFSAD1000_V1` | 1 km | nominal 2010 | Distinguishes **irrigated vs. rainfed** cropping (beyond generic cropland) — irrigation buffers seasonal food gaps and shapes the crop mix (cereals vs. pulses/vegetables) that determines micronutrient intake. |
| 9 | Terrain / topography | SRTM 30 m DEM (derive slope, TPI, ruggedness, elevation) | `USGS/SRTMGL1_003` (or `MERIT/DEM/v1_0_3` for hydro-conditioned) | 30 m (90 m MERIT) | 2000 (SRTM) | Elevation/slope shape agro-ecology, malaria transmission, soil retention, and accessibility — joint upstream drivers of diet and infection-mediated micronutrient status. Cheap, stable covariates. |
| 10 | Fire / biomass burning | FireCCI51 Burned Area (or MCD64A1) | `ESA/CCI/FireCCI/5_1` (or `MODIS/061/MCD64A1`) | 250 m (500 m MCD64A1) | 2001–2020, monthly | Burned area proxies smoke exposure (inflammation/iron axis), shifting cultivation intensity, and biomass-fuel reliance — links to both diet-source disruption and PM2.5 confounding of iron markers. |
| 11 | Air pollution / aerosols (recent) | Sentinel-5P TROPOMI: UV Aerosol Index, NO2, SO2, CO | `COPERNICUS/S5P/OFFL/L3_AER_AI`, `.../L3_NO2`, `.../L3_SO2`, `.../L3_CO` | ~1 km (3.5×5.5 km native) | 2018–present, ~daily | Higher-temporal-resolution pollution/aerosol signal for the survey period (post-2018); complements the historical ACAG PM2.5 climatology for inflammation-related iron confounding. |
| 12 | Drought indices (standardized) | Derived SPI/SPEI from TerraClimate; or GRIDMET DROUGHT (CONUS-only, reference) | derive from `IDAHO_EPSCOR/TERRACLIMATE` (`pr`, `pet`) | ~4 km | 1958–present, monthly | Standardized precipitation-evapotranspiration anomalies capture episodic shocks (lean-season severity) that depress dietary diversity and elevate deficiency risk; not available as a ready-made global GEE band, must be computed. |

### Secondary / lower-priority candidates

| Domain | Dataset | GEE asset ID | Resolution | Coverage | Rationale |
|--------|---------|--------------|------------|----------|-----------|
| Enhanced soil moisture | NASA-USDA Enhanced SMAP Global Soil Moisture | `NASA_USDA/HSL/SMAP10KM_soil_moisture` | ~10 km | 2015–2022 | Surface/subsurface SM + anomalies; alternative/complement to SMAP L4 with explicit anomaly bands. |
| Yield / production (gridded) | MapSPAM (crop area, yield, production, 40+ crops) | ingest as raster (no native `ee` ID) | ~10 km | 2010 (v2010) | Direct crop-specific production/yield maps (pulses, vegetables, cereals) — the closest proxy to local micronutrient-rich food supply. |
| Land-cover dynamics | ESA WorldCover / ESA-CCI Land Cover | `ESA/WorldCover/v200`; `ESA/CCI/...` | 10 m / 300 m | 2020–2021 / 1992–present | Finer or time-series land cover than current layer; captures cropland/wetland/forest mosaics tied to wild-food and fishery access. |
| Crop yield (field-scale) | Note: no global ready-made yield band in GEE core catalog | — | — | — | Yield typically modeled, not a primitive dataset; flagged for awareness. |

---

## Notes and caveats

- **Mechanistic clustering.** The candidates fall into three causal pathways to micronutrient status: (i) **food availability/diet quality** — market access (#1), surface water/fishing (#5), livestock/ASF (#7), crop system (#8), yield (MapSPAM); (ii) **agro-climatic stress on production** — drought/water balance (#2, #12), soil moisture (#3), ET (#4), fire (#10); (iii) **infection/inflammation confounding of biomarkers** — PM2.5 (#6), aerosols (#11), and partly fire (#10). Pathway (iii) is especially relevant given the memory note that cross-survey ferritin *levels* are confounded; inflammation proxies may help explain or adjust ferritin-based iron estimates.
- **Effective-n constraint.** Per project memory, proxy predictors are Admin-2-level with effective n = number of areas (14–87). Adding many high-resolution layers risks overfitting; prioritize a small set of mechanistically distinct, Admin-2-aggregated summaries (mean/SD per area) rather than dumping all bands.
- **Single-epoch vs. time-varying.** Accessibility (2015), GFSAD (2010), GLW (2010–2020), SRTM (2000) are effectively static — fine as structural covariates but cannot capture survey-year dynamics. TerraClimate, SMAP, MOD16, S5P, FireCCI are time-varying and should be windowed to each survey's fieldwork period.
- **Ingest-required datasets.** GLW (global) and MapSPAM are not first-class `ee.ImageCollection` IDs in the core catalog; they are distributed via FAO/Harvard Dataverse / mapspam.info and must be uploaded as Earth Engine assets (or accessed via community mirrors). The AGLW community implementation provides a GEE-native annual livestock series (1961–2021, ~5 km).
- **Resolution mismatch.** Several datasets (SMAP 9–10 km, TerraClimate 4 km, GLW 10 km, accessibility 1 km) are coarser than Admin-2 units in small countries (e.g., The Gambia). Aggregate with care; very coarse layers may contribute little within-country variance.

---

## Sources

- SMAP L4 (`NASA/SMAP/SPL4SMGP/008`): https://developers.google.com/earth-engine/datasets/catalog/NASA_SMAP_SPL4SMGP_008
- SMAP intro / NASA-USDA enhanced: https://developers.google.com/earth-engine/tutorials/community/smap-soil-moisture ; tag: https://developers.google.com/earth-engine/datasets/tags/soil-moisture
- MODIS MOD16A2 ET: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD16A2 ; gap-filled: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD16A2GF
- JRC Global Surface Water v1.4: https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_GlobalSurfaceWater
- Accessibility to Cities 2015 (Weiss et al. 2018, Nature): https://developers.google.com/earth-engine/datasets/catalog/Oxford_MAP_accessibility_to_cities_2015_v1_0 ; https://www.nature.com/articles/nature25181
- ACAG / van Donkelaar PM2.5 (community catalog): https://gee-community-catalog.org/projects/global_pm25/ ; ACAG: https://sites.wustl.edu/acag/surface-pm2-5/
- Gridded Livestock of the World (GLW): https://www.fao.org/land-water/land/land-governance/land-resources-planning-toolbox/category/details/en/c/1236449/ ; https://www.nature.com/articles/sdata2018227 ; community: https://gee-community-catalog.org/projects/gridded_livestock/
- GFSAD1000 cropland extent (`USGS/GFSAD1000_V1`): https://developers.google.com/earth-engine/datasets/catalog/USGS_GFSAD1000_V1
- SRTM 30 m (`USGS/SRTMGL1_003`): https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003 ; MERIT DEM: https://developers.google.com/earth-engine/datasets/catalog/MERIT_DEM_v1_0_3
- FireCCI51 (`ESA/CCI/FireCCI/5_1`): https://developers.google.com/earth-engine/datasets/catalog/ESA_CCI_FireCCI_5_1 ; MCD64A1: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD64A1
- Sentinel-5P TROPOMI (NO2/SO2/AER_AI): https://developers.google.com/earth-engine/datasets/catalog/sentinel-5p/
- TerraClimate (`IDAHO_EPSCOR/TERRACLIMATE`): https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE ; SPEI-from-TerraClimate method: https://benny.istan.to/blog/20211114-terraclimate-data-and-standardized-precipitation-evapotranspiration-index-spei
- GRIDMET DROUGHT (US reference): https://developers.google.com/earth-engine/datasets/catalog/GRIDMET_DROUGHT
- MapSPAM: https://www.mapspam.info/data/ ; https://essd.copernicus.org/articles/12/3545/2020/
