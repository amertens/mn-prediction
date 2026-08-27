"""
sandbox_parsimony/python/extract_tanzania_gee.py

Tanzania is the only country without data/<Country>_GEE_rasters/, which is why
the pooled covariate intersection loses accessibility, LST, TerraClimate,
productivity, EVI, LAI and land cover for EVERY country (FINDINGS.md section 4).

This pulls those families for Tanzania's Admin-2 polygons directly through the
Earth Engine API -- no Drive export, no task queue. Values come back inline via
reduceRegions().getInfo(), chunked so the payload stays under EE's response
limit.

AUTH: uses the refresh token already stored in ~/.config/earthengine/credentials
      (account amertens@berkeley.edu). ee.Initialize() picks it up with no
      interactive step.

YEAR: Tanzania's survey_year is 2010 (R/config.R), so annual composites are
      built for 2010 wherever the product covers it. Copernicus land cover only
      starts in 2015, so its earliest year is used and the column is named
      accordingly -- a temporal mismatch that must not be silently hidden.

OUTPUT: sandbox_parsimony/out/tanzania_gee_admin2.csv, one row per Admin-2
        polygon, columns named to match the gee_* vocabulary the other four
        countries already use.
"""
import json
import os
import sys
import time

import ee

OUT = "sandbox_parsimony/out/tanzania_gee_admin2.csv"
GEOJSON = "sandbox_parsimony/out/_tz_admin2.geojson"
YEAR = 2010
CHUNK = 20          # polygons per getInfo call
SCALE = 1000        # metres; matches the 1 km products

ee.Initialize()
print("EE initialised", flush=True)

if not os.path.exists(GEOJSON):
    sys.exit(f"missing {GEOJSON} -- run the R exporter first")

with open(GEOJSON) as fh:
    gj = json.load(fh)
feats = gj["features"]
print(f"{len(feats)} Admin-2 polygons", flush=True)


def annual(collection, bands, reducer="mean", year=YEAR):
    """Annual composite of one collection, renamed with a gee_ prefix."""
    ic = (ee.ImageCollection(collection)
          .filterDate(f"{year}-01-01", f"{year+1}-01-01")
          .select(bands))
    img = ic.mean() if reducer == "mean" else ic.reduce(getattr(ee.Reducer, reducer)())
    return img


def build_stack():
    """One multi-band image holding every covariate we want."""
    parts = []

    # --- travel time to cities (static, 2015) -----------------------------
    parts.append(ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
                 .select(["accessibility"])
                 .rename(["gee_accessibility"]))

    # --- MODIS land surface temperature, day and night --------------------
    lst = annual("MODIS/061/MOD11A2", ["LST_Day_1km", "LST_Night_1km"])
    # scale factor 0.02, kelvin
    parts.append(lst.multiply(0.02)
                 .rename(["gee_lst_day_annual_mean", "gee_lst_night_annual_mean"]))

    # --- MODIS vegetation indices ----------------------------------------
    veg = annual("MODIS/061/MOD13A2", ["NDVI", "EVI"]).multiply(0.0001)
    parts.append(veg.rename(["gee_ndvi_modis_annual_mean",
                             "gee_dailyevi_annual_mean"]))

    # --- MODIS leaf area index -------------------------------------------
    lai = annual("MODIS/061/MCD15A3H", ["Lai"]).multiply(0.1)
    parts.append(lai.rename(["gee_lai8days_annual_mean"]))

    # --- MODIS productivity (annual GPP / NPP) ----------------------------
    prod = (ee.ImageCollection("MODIS/061/MOD17A3HGF")
            .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01")
            .select(["Gpp", "Npp"]).mean().multiply(0.0001))
    parts.append(prod.rename(["gee_productivity_gpp", "gee_productivity_npp"]))

    # --- TerraClimate: the full climate/water suite -----------------------
    tc_bands = ["aet", "def", "pdsi", "pet", "pr", "ro", "soil", "srad",
                "tmmn", "tmmx", "vap", "vpd", "vs"]
    tc = annual("IDAHO_EPSCOR/TERRACLIMATE", tc_bands)
    parts.append(tc.rename([f"gee_terraclimate_{b}" for b in tc_bands]))

    # --- Copernicus land cover fractions ----------------------------------
    # NB 2015 is the earliest year this product exists; the survey is 2010.
    # The column names carry _2015 so the mismatch stays visible.
    lc = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2015")
    lc_bands = ["bare-coverfraction", "crops-coverfraction",
                "grass-coverfraction", "shrub-coverfraction",
                "tree-coverfraction", "urban-coverfraction",
                "water-permanent-coverfraction"]
    parts.append(lc.select(lc_bands).rename(
        ["gee_landcoverlayers_" + b.replace("-", "_") + "_2015" for b in lc_bands]))

    return ee.Image.cat(parts)


stack = build_stack()
band_names = stack.bandNames().getInfo()
print(f"{len(band_names)} bands: {band_names}", flush=True)

rows = []
for i in range(0, len(feats), CHUNK):
    chunk = feats[i:i + CHUNK]
    fc = ee.FeatureCollection([
        ee.Feature(ee.Geometry(f["geometry"]), {"Admin2": f["properties"]["Admin2"]})
        for f in chunk])
    for attempt in range(3):
        try:
            res = stack.reduceRegions(
                collection=fc, reducer=ee.Reducer.mean(), scale=SCALE
            ).getInfo()
            break
        except Exception as exc:                      # transient EE errors
            if attempt == 2:
                print(f"  chunk {i} FAILED: {str(exc)[:200]}", flush=True)
                res = None
            else:
                time.sleep(5 * (attempt + 1))
    if res is None:
        continue
    for f in res["features"]:
        p = f["properties"]
        rows.append({"Admin2": p.get("Admin2"),
                     **{b: p.get(b) for b in band_names}})
    print(f"  {min(i + CHUNK, len(feats))}/{len(feats)} districts", flush=True)

if not rows:
    sys.exit("no rows extracted")

import csv
cols = ["Admin2"] + band_names
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"\nwrote {OUT}: {len(rows)} rows x {len(cols)} cols", flush=True)
