"""
sandbox_parsimony/python/fetch_tanzania_rasters.py

Fill data/Tanzania_GEE_rasters/ from the Earth Engine API.

WHY THIS AND NOT A ZONAL EXTRACT. Tanzania is the only country without a local
raster directory, and because the pooled model intersects covariates by NAME,
its absence deletes accessibility, LST, TerraClimate, productivity, EVI, LAI and
land cover for EVERY country (FINDINGS.md section 4). Writing GeoTIFFs into the
directory the pipeline already looks in fixes the cause: extract_gee_admin2()
and extract_area_covariates() then treat Tanzania exactly like the other four,
using the same .append_gee_zonal_cols() naming, with no new code path and no
legacy-parity CSV fallback.

AUTH: the refresh token in ~/.config/earthengine/credentials
      (amertens@berkeley.edu). ee.Initialize() finds it; no interactive step.

YEAR: Tanzania's survey_year is 2010 (R/config.R). Annual composites use 2010
      wherever the product covers it. Copernicus land cover begins in 2015, so
      its filename carries 2015 -- the temporal mismatch stays visible rather
      than being buried in a generic name.

Each request is capped at a few bands so the GeoTIFF stays inside Earth
Engine's direct-download size limit.
"""
import io
import os
import sys
import time
import zipfile

import ee
import requests

OUT_DIR = "data/Tanzania_GEE_rasters"
YEAR = 2010
SCALE = 1000

ee.Initialize()
os.makedirs(OUT_DIR, exist_ok=True)

# Tanzania bounding box, from the GADM polygons
REGION = ee.Geometry.Rectangle([29.30, -11.80, 40.50, -0.95], proj="EPSG:4326",
                               geodesic=False)


def annual_mean(cid, bands, scale_factor=1.0, year=YEAR):
    ic = (ee.ImageCollection(cid)
          .filterDate(f"{year}-01-01", f"{year + 1}-01-01")
          .select(bands))
    img = ic.mean()
    return img.multiply(scale_factor) if scale_factor != 1.0 else img


# Each entry becomes one .tif named to match the other countries' vocabulary:
# .append_gee_zonal_cols() lowercases the filename and prefixes gee_, so
# "Accessibility_Tanzania_2015.tif" -> gee_accessibility, matching the existing
# Accessibility_<Country>_2019.tif files.
JOBS = [
    ("Accessibility_Tanzania_2015",
     lambda: ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
     .select(["accessibility"])),

    ("LST_Day_Monthly_Tanzania_2010",
     lambda: annual_mean("MODIS/061/MOD11A2", ["LST_Day_1km"], 0.02)),

    ("LST_Night_Annual_Mean_Tanzania_2010",
     lambda: annual_mean("MODIS/061/MOD11A2", ["LST_Night_1km"], 0.02)),

    ("NDVI_Tanzania_2010",
     lambda: annual_mean("MODIS/061/MOD13A2", ["NDVI"], 0.0001)),

    ("DailyEVI_Tanzania_2010",
     lambda: annual_mean("MODIS/061/MOD13A2", ["EVI"], 0.0001)),

    ("LAI8Days_Tanzania_2010",
     lambda: annual_mean("MODIS/061/MCD15A3H", ["Lai"], 0.1)),

    ("Productivity_Tanzania_2010",
     lambda: (ee.ImageCollection("MODIS/061/MOD17A3HGF")
              .filterDate(f"{YEAR}-01-01", f"{YEAR + 1}-01-01")
              .select(["Gpp", "Npp"]).mean().multiply(0.0001))),

    # TerraClimate carries 13 distinct variables; split so each download stays
    # inside the size limit.
    ("TerraClimate_Tanzania_2010_a",
     lambda: annual_mean("IDAHO_EPSCOR/TERRACLIMATE",
                         ["aet", "def", "pdsi", "pet"])),
    ("TerraClimate_Tanzania_2010_b",
     lambda: annual_mean("IDAHO_EPSCOR/TERRACLIMATE",
                         ["pr", "ro", "soil", "srad"])),
    ("TerraClimate_Tanzania_2010_c",
     lambda: annual_mean("IDAHO_EPSCOR/TERRACLIMATE",
                         ["tmmn", "tmmx", "vap", "vpd", "vs"])),

    ("LandCoverLayers_Tanzania_2015_a",
     lambda: ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2015")
     .select(["bare-coverfraction", "crops-coverfraction",
              "grass-coverfraction"])),
    ("LandCoverLayers_Tanzania_2015_b",
     lambda: ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2015")
     .select(["shrub-coverfraction", "tree-coverfraction",
              "urban-coverfraction"])),
]


def fetch(name, builder):
    dest = os.path.join(OUT_DIR, name + ".tif")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"  cached  {name}", flush=True)
        return True
    for attempt in range(3):
        try:
            img = builder().toFloat()
            url = img.getDownloadURL({
                "region": REGION, "scale": SCALE, "crs": "EPSG:4326",
                "format": "GEO_TIFF", "filePerBand": False})
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            payload = r.content
            # EE returns either a bare GeoTIFF or a zip of them
            if payload[:2] == b"PK":
                with zipfile.ZipFile(io.BytesIO(payload)) as z:
                    inner = [n for n in z.namelist() if n.endswith(".tif")]
                    if not inner:
                        raise RuntimeError("zip contained no tif")
                    payload = z.read(inner[0])
            with open(dest, "wb") as fh:
                fh.write(payload)
            print(f"  OK      {name}  ({len(payload) / 1e6:.1f} MB)", flush=True)
            return True
        except Exception as exc:
            msg = str(exc)[:180]
            if attempt == 2:
                print(f"  FAILED  {name}: {msg}", flush=True)
                return False
            time.sleep(5 * (attempt + 1))
    return False


ok = 0
for name, builder in JOBS:
    if fetch(name, builder):
        ok += 1

print(f"\n{ok}/{len(JOBS)} rasters in {OUT_DIR}", flush=True)
if ok < len(JOBS):
    sys.exit(1)
