# =============================================================================
# src/GEE/gee_layer_manifest.R
#
# Canonical Earth Engine layer manifest for the proxy-covariate extraction,
# reconstructed from data/GEE/GEE_export_metadata.xlsx (the analyst's record of
# exactly which EE assets each raster family came from) and cross-checked
# against the raster filenames in data/<Country>_GEE_rasters/.
#
# This is the single source of truth for the auto-extraction in
# src/GEE/extract_gee_ee_api.R, so a new country (Tanzania) pulls the SAME
# layers as Gambia/Ghana/Sierra Leone/Malawi.
#
# Each entry:
#   family   : the raster-family name (matches the legacy .tif / gee_a2_ stem)
#   asset    : EE asset id
#   kind     : "ImageCollection" or "Image"
#   bands    : character vector to .select(), or NULL for all bands
#   temporal : "filter_year"  -> filterDate(year-01-01, year-12-31) then reduce
#              "static"        -> single Image / time-invariant
#              "nearest_epoch" -> ImageCollection of epochs; pick image nearest year
#   reducer  : "mean" (default) — server-side ee.Reducer
#   scale    : nominal export scale (m), from the metadata sheet
#   avail    : [min,max] year the asset supports (for nearest-year fallback)
#   scale_factor : optional multiplier applied after reduction (e.g. MODIS NDVI)
#   note     : provenance / caveats
#
# NON-EE families (handled outside this manifest — already global statics in the
# repo, NOT pulled from the API): SPAM crop (data/SPAM/, Harvard Dataverse),
# Koppen-Geiger (data/Koppen_geiger_tif/), Agro-Ecological Zones
# (data/"Agro-Ecological Zones..."/), FEWS NET "NDVI Mean Anomaly" (USGS),
# GLW4 livestock (data/GEE/GLW4-2020...tif, FAO). These flow into gee_a2_ via
# the existing local-raster path in build_gee_admin2.R, unchanged.
# =============================================================================

GEE_LAYER_MANIFEST <- list(

  # ── Vegetation / greenness ──────────────────────────────────────────────
  list(family = "NDVI", asset = "NOAA/CDR/AVHRR/NDVI/V5", kind = "ImageCollection",
       bands = "NDVI", temporal = "filter_year", yoy = TRUE, scale = 5000, avail = c(1981, 2013),
       scale_factor = 1e-4, note = "NOAA CDR AVHRR NDVI v5; daily->annual mean."),
  list(family = "DailyEVI", asset = "MODIS/061/MOD13Q1", kind = "ImageCollection",
       bands = "EVI", temporal = "filter_year", yoy = TRUE, scale = 250, avail = c(2000, 2024),
       scale_factor = 1e-4,
       note = "MOD13Q1 16-day EVI composite (~23 imgs/yr). Swapped from MOD09GA daily EVI (365 imgs) for speed; scale 0.0001."),
  list(family = "LAI8days", asset = "MODIS/061/MOD15A2H", kind = "ImageCollection",
       bands = c("Lai_500m"), temporal = "filter_year", yoy = TRUE, scale = 500, avail = c(2000, 2026),
       note = "MOD15A2H 8-day LAI; annual mean."),
  list(family = "Productivity", asset = "MODIS/061/MOD17A3HGF", kind = "ImageCollection",
       bands = "Npp", temporal = "filter_year", yoy = TRUE, scale = 500, avail = c(2001, 2025),
       note = "Annual NPP (gap-filled)."),
  list(family = "WAPOR", asset = "FAO/WAPOR/2/L1_NPP_D", kind = "ImageCollection",
       bands = NULL, temporal = "filter_year", yoy = TRUE, scale = 250, avail = c(2009, 2023),
       note = "WaPOR dekadal NPP; annual mean."),

  # ── Climate / precipitation / temperature ───────────────────────────────
  list(family = "TRMM", asset = "TRMM/3B43V7", kind = "ImageCollection",
       bands = "precipitation", temporal = "filter_year", yoy = TRUE, scale = 5000, avail = c(1998, 2019),
       note = "TRMM 3B43 monthly precip; mean of 12 months."),
  list(family = "TerraClimate", asset = "IDAHO_EPSCOR/TERRACLIMATE", kind = "ImageCollection",
       bands = NULL, temporal = "filter_year", yoy = TRUE, scale = 4638, avail = c(1958, 2025),
       note = "Monthly climate/water-balance; annual mean of all bands."),
  list(family = "FLDAS", asset = "NASA/FLDAS/NOAH01/C/GL/M/V001", kind = "ImageCollection",
       bands = NULL, temporal = "filter_year", yoy = TRUE, scale = 11132, avail = c(1982, 2026),
       note = "FLDAS monthly land-surface; annual mean of all *_tavg/_inst bands."),
  list(family = "LST_Night", asset = "Oxford/MAP/LST_Night_5km_Monthly", kind = "ImageCollection",
       bands = NULL, temporal = "filter_year", yoy = TRUE, scale = 5000, avail = c(2001, 2015),
       note = "DEPRECATED Oxford MAP gap-filled night LST; monthly->annual mean. Asset may be retired; if unavailable substitute MODIS/061/MOD11A2 LST_Night_1km."),
  list(family = "Atmosphere", asset = "MODIS/061/MOD08_M3", kind = "ImageCollection",
       bands = c("Cloud_Fraction_Mean_Mean", "Atmospheric_Water_Vapor_Mean_Mean"),
       temporal = "filter_year", yoy = TRUE, scale = 111320, avail = c(2000, 2026),
       note = "MOD08_M3 monthly; restricted to cloud fraction + water vapor (all-bands was 628 cols of noise). Skips gracefully if a band name is off."),
  list(family = "AerosolOptical", asset = "MODIS/061/MOD08_M3",
       bands = "Aerosol_Optical_Depth_Land_Ocean_Mean_Mean", kind = "ImageCollection",
       temporal = "filter_year", yoy = TRUE, scale = 111320, avail = c(2000, 2026),
       note = "MOD08_M3 monthly AOD (12 imgs/yr). Swapped from MCD19A2 daily (365 imgs) for speed. If band name errors it skips; Atmosphere (MOD08_M3 all-bands) also carries AOD."),

  # ── Population / settlement / built environment ─────────────────────────
  list(family = "GHSPOP", asset = "JRC/GHSL/P2023A/GHS_POP", kind = "nearest_epoch",
       bands = "population_count", temporal = "nearest_epoch", scale = 100, avail = c(1975, 2030),
       note = "GHSL population, 5-yr epochs; pick epoch nearest survey year."),
  list(family = "GHSBUILTS", asset = "JRC/GHSL/P2023A/GHS_BUILT_S", kind = "nearest_epoch",
       bands = c("built_surface", "built_surface_nres"), temporal = "nearest_epoch",
       scale = 100, avail = c(1975, 2030), note = "GHSL built-up surface, 5-yr epochs."),
  list(family = "GHSL_SMOD", asset = "JRC/GHSL/P2023A/GHS_SMOD", kind = "nearest_epoch",
       bands = "smod_code", temporal = "nearest_epoch", scale = 1000, avail = c(1975, 2030),
       note = "GHSL Degree of Urbanisation (settlement model), 5-yr epochs."),
  list(family = "PopDensity", asset = "CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density",
       kind = "nearest_epoch", bands = "unwpp-adjusted_population_density",
       temporal = "nearest_epoch", scale = 1000, avail = c(2000, 2020),
       note = "GPW v4.11 UN-adjusted pop density, 5-yr epochs."),
  list(family = "GPW_Demographic", asset = "CIESIN/GPWv411/GPW_Basic_Demographic_Characteristics",
       kind = "ImageCollection", bands = NULL, temporal = "static", to_bands = TRUE,
       scale = 1000, avail = c(2010, 2010),
       note = "GPW v4.11 basic demographics — 2010 round. It's an ImageCollection of age-sex groups; toBands() stacks them into one multiband image."),
  list(family = "WSF", asset = "DLR/WSF/WSF2015/v1", kind = "Image",
       bands = NULL, temporal = "static", scale = 10, avail = c(2015, 2015),
       note = "World Settlement Footprint 2015 (single epoch)."),
  list(family = "CCNL", asset = "BNU/FGS/CCNL/v1", kind = "nearest_epoch",
       bands = "b1", temporal = "nearest_epoch", scale = 1000, avail = c(1992, 2013),
       note = "Consistent corrected nighttime lights (DMSP-OLS), annual; pick nearest <=2013."),

  # ── Accessibility / human modification / elevation ──────────────────────
  list(family = "Accessibility",
       asset = "Oxford/MAP/accessibility_to_healthcare_2019",
       kind = "Image", bands = "accessibility", temporal = "static", scale = 1000,
       avail = c(2019, 2019), note = "MAP accessibility to healthcare 2019 (single epoch)."),
  list(family = "GlobalHumanModification", asset = "CSP/HM/GlobalHumanModification",
       kind = "ImageCollection", bands = "gHM", temporal = "static", scale = 1000,
       avail = c(2016, 2016), note = "CSP gHM 2016 — single-image ImageCollection; mean() returns it."),
  list(family = "Elevation", asset = "USGS/SRTMGL1_003", kind = "Image",
       bands = "elevation", temporal = "static", scale = 30, avail = c(2000, 2000),
       note = "SRTM 30m DEM (static)."),

  # ── Land cover / cropland ───────────────────────────────────────────────
  list(family = "LandCoverType", asset = "MODIS/061/MCD12Q1", kind = "nearest_epoch",
       bands = "LC_Type1", temporal = "nearest_epoch", scale = 500, avail = c(2001, 2024),
       note = "MCD12Q1 yearly land-cover; pick nearest year."),
  list(family = "LandCoverLayers", asset = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
       kind = "nearest_epoch", bands = "discrete_classification", temporal = "nearest_epoch",
       scale = 100, avail = c(2015, 2019),
       note = "Copernicus CGLS-LC100 C3; 2015 is nearest for pre-2015 surveys."),
  list(family = "ESA_WorldCereal", asset = "ESA/WorldCereal/2021/MARKERS/v100",
       kind = "ImageCollection", bands = "classification", temporal = "static",
       scale = 10, avail = c(2021, 2021),
       note = "Active cropland markers, 2021 only (single epoch for all countries)."),
  list(family = "GPW_Grasslands",
       asset = "projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c",
       kind = "nearest_epoch", bands = NULL, temporal = "nearest_epoch", scale = 30,
       avail = c(2000, 2020), note = "Global Pasture Watch grassland class; nearest 5-yr."),

  # ── iSDAsoil micronutrients / chemistry (Africa, 2001-2017 synthesis) ────
  # NOTE: these are iSDAsoil (gee_a2_Soil*), DISTINCT from the SoilGrids/ISRIC
  # `soil_` domain built by scripts/build_soilgrids_admin2.R. Both are kept.
  list(family = "SoilCEC",        asset = "ISDASOIL/Africa/v1/cation_exchange_capacity", kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil CEC 0-20cm."),
  list(family = "SoilCalcium",    asset = "ISDASOIL/Africa/v1/calcium_extractable",      kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable Ca."),
  list(family = "SoilIron",       asset = "ISDASOIL/Africa/v1/iron_extractable",         kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable Fe."),
  list(family = "SoilMagnesium",  asset = "ISDASOIL/Africa/v1/magnesium_extractable",    kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable Mg."),
  list(family = "SoilPhosphorus", asset = "ISDASOIL/Africa/v1/phosphorus_extractable",   kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable P."),
  list(family = "SoilPotassium",  asset = "ISDASOIL/Africa/v1/potassium_extractable",     kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable K."),
  list(family = "SoilSulfur",     asset = "ISDASOIL/Africa/v1/sulphur_extractable",       kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable S."),
  list(family = "SoilZinc",       asset = "ISDASOIL/Africa/v1/zinc_extractable",          kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable Zn."),
  list(family = "SoilTotalCarbon",asset = "ISDASOIL/Africa/v1/carbon_total",              kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil total C."),
  list(family = "SoilNitrogen",   asset = "ISDASOIL/Africa/v1/nitrogen_total",            kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil total N."),
  list(family = "SoilAluminium",  asset = "ISDASOIL/Africa/v1/aluminium_extractable",     kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil extractable Al."),
  list(family = "SoilpH",         asset = "ISDASOIL/Africa/v1/ph",                        kind = "Image", bands = c("mean_0_20", "mean_20_50"), temporal = "static", scale = 30, avail = c(2001, 2017), note = "iSDAsoil pH 0-20cm — governs Zn/Fe crop bioavailability (soil->crop->human micronutrient link; Broadley/Joy GeoNutrition)."),

  # ── High-value additions (literature-driven; new causal pathways) ─────────
  # Animal-source-food / fisheries pathway (B12, heme iron, preformed vit A, Zn):
  list(family = "SurfaceWater", asset = "JRC/GSW1_4/GlobalSurfaceWater", kind = "Image",
       bands = c("occurrence", "seasonality", "recurrence"), temporal = "static",
       scale = 30, avail = c(1984, 2021),
       note = "JRC Global Surface Water — district water occurrence/seasonality; proxy for inland fisheries, irrigation, flood-recession & dairy value chains. New hydrology pathway."),
  # Market access -> dietary diversity (complements accessibility_to_healthcare):
  list(family = "AccessibilityCities", asset = "Oxford/MAP/accessibility_to_cities_2015_v1_0",
       kind = "Image", bands = "accessibility", temporal = "static", scale = 1000,
       avail = c(2015, 2015),
       note = "Travel time to nearest city (Weiss et al. 2018, Nature). Market access -> ability to buy diverse nutrient-dense foods."),
  # Agricultural seasonality (seasonal food availability -> seasonal deficiency):
  list(family = "Phenology", asset = "MODIS/061/MCD12Q2", kind = "ImageCollection",
       bands = c("Greenup_1", "Senescence_1", "Dormancy_1"), temporal = "nearest_epoch",
       scale = 500, avail = c(2001, 2023),
       note = "MODIS land-cover dynamics — growing-season timing/length; you have NDVI/EVI levels but no phenology.")
)

#' Resolve the nearest available year for a layer given a target survey year.
#' @return integer year within the layer's avail window closest to target.
gee_nearest_year <- function(layer, target_year) {
  lo <- layer$avail[1]; hi <- layer$avail[2]
  max(lo, min(hi, target_year))
}

#' Families NOT pulled from EE (already global statics in the repo).
GEE_NON_EE_FAMILIES <- c("SPAM", "Koppen_geiger", "Agro_Ecological_Zones",
                         "NDVI Mean Anomaly", "GLW4")
