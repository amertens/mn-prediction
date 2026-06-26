#!/usr/bin/env python3
"""Export fishing-effort and market-access rasters from Google Earth Engine.

Motivated by Bondi-Kelly et al. 2022 (AAAI, "Predicting Micronutrient Deficiency
with Publicly Available Satellite Data"), whose strongest predictors were
distance-to-market and aquatic-food access. We hold no aquatic-food or
market-access layer; travel-time-to-healthcare is the closest proxy.

This project does NOT call Earth Engine live -- the R pipeline auto-discovers
per-country GeoTIFFs dropped into `data/<Country>_GEE_rasters/` (see
`.append_gee_zonal_cols()` in R/admin2_analysis.R). This script produces those
rasters. After the EE export tasks finish, download/move each GeoTIFF into the
matching country raster folder. The file name's country token must be one of
{Gambia, Ghana, Sierra_Leone, Malawi, Cote_dIvoire} so the R name-stripper turns
`GFWFishingHours_Ghana.tif` into the column `gee_gfwfishinghours`.

CAVEAT (coastal vs landlocked): fishing effort is meaningful for the coastal
countries (Gambia, Ghana, Sierra Leone, Cote d'Ivoire) but Malawi is landlocked,
so its raster is ~all-zero. We `unmask(0)` because inland marine-fishing access
genuinely is 0 (the correct value, not missing). The all-zero Malawi column is
then dropped by `.screen_vars()` from the cross-country LOCO intersection -- so
fishing effort helps the *within-country* area-level SuperLearner for coastal
countries but, correctly, does not enter the LOCO transportability models. This
is option #1 in the plan; it is honest and requires no extra handling.

Usage:
    earthengine authenticate          # one-time
    python scripts/export_gee_fishing_market.py
The exports land in your Google Drive folder `mn_prediction_gee`; download them
into data/<Country>_GEE_rasters/.
"""

import ee

ee.Initialize()

DRIVE_FOLDER = "mn_prediction_gee"

# Survey years match the four micronutrient surveys (+ CIV external validation).
COUNTRIES = {
    "Gambia":       {"iso": "GMB", "year": 2018},
    "Ghana":        {"iso": "GHA", "year": 2017},
    "Sierra_Leone": {"iso": "SLE", "year": 2013},
    "Malawi":       {"iso": "MWI", "year": 2016},   # landlocked -> fishing ~ 0
    "Cote_dIvoire": {"iso": "CIV", "year": 2017},
}

# FAO GAUL level-0 national boundaries, joined by ISO3 code.
GAUL0 = ee.FeatureCollection("FAO/GAUL_SIMPLIFIED_500m/2015/level0")


def country_geom(iso3):
    """National geometry for an ISO3 code (falls back to ADM0 name match)."""
    fc = GAUL0.filter(ee.Filter.eq("ISO3_CODE", iso3))
    return fc.geometry()


def export(image, name, geom, scale=1000):
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=name,
        folder=DRIVE_FOLDER,
        fileNamePrefix=name,
        region=geom,
        scale=scale,
        maxPixels=int(1e10),
    )
    task.start()
    print(f"  started export: {name} (scale={scale} m)")


def main():
    for name, c in COUNTRIES.items():
        geom = country_geom(c["iso"])
        yr = c["year"]
        print(f"{name} ({c['iso']}, {yr}):")

        # --- Fishing effort: Global Fishing Watch daily apparent fishing hours ---
        # GFW/GFF/V2/fishing_hours is a daily ImageCollection; sum to an annual
        # total per pixel, then unmask(0) so inland pixels read 0 rather than NA.
        fishing = (
            ee.ImageCollection("GFW/GFF/V2/fishing_hours")
            .filterDate(f"{yr}-01-01", f"{yr}-12-31")
            .select("fishing_hours")
            .sum()
            .unmask(0)
            .clip(geom)
            .rename("fishing_hours")
        )
        export(fishing, f"GFWFishingHours_{name}", geom)

        # --- Market access: Weiss et al. 2018 travel time to cities (= markets) ---
        access = (
            ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
            .select("accessibility")
            .clip(geom)
            .rename("market_access_tt")
        )
        export(access, f"MarketAccessTT_{name}", geom)

        # --- Optional stronger market signal: GHSL built-up surface (2020) ---
        # We already extract a GHSL built-up layer (GHSBUILTS_*); keep this finer
        # variant only if it adds signal beyond that, else delete the .tif.
        built = (
            ee.Image("JRC/GHSL/P2023A/GHS_BUILT_S/2020")
            .select("built_surface")
            .clip(geom)
            .rename("built_surface")
        )
        export(built, f"BuiltUpSurface_{name}", geom)

    print(
        "\nAll export tasks started. Monitor at https://code.earthengine.google.com/tasks\n"
        f"When complete, download from Drive folder '{DRIVE_FOLDER}' into\n"
        "data/<Country>_GEE_rasters/ (e.g. data/Ghana_GEE rasters/GFWFishingHours_Ghana.tif)."
    )


if __name__ == "__main__":
    main()
