"""
sandbox_parsimony/python/fetch_tanzania_rasters2.py

STEP 3: the remaining GEE families Tanzania still lacks after
fetch_tanzania_rasters.py -- ESA WorldCereal, FLDAS, GHS built surface,
GHS population, MODIS land-cover type, WAPOR NPP, and an LST_Night stack whose
FILENAME matches the other countries' convention.

Naming matters as much as the data. .append_gee_zonal_cols() derives every
column name from the filename plus the band names, and the pooled/LOCO models
intersect covariates by EXACT NAME. So each output here mirrors the Ghana file
it is meant to line up with, band names included:

  Ghana                                  Tanzania (written here)
  ESA_WorldCereal_2021_Ghana.tif      -> ESA_WorldCereal_2021_Tanzania.tif
  FLDAS_2016_Annual_Mean_Ghana.tif    -> FLDAS_2010_Annual_Mean_Tanzania.tif
  FLDAS_2016_Monthly_Ghana.tif        -> FLDAS_2010_Monthly_Tanzania.tif
  GHSBUILTS_Ghana_2015.tif            -> GHSBUILTS_Tanzania_2015.tif
  GHSPOP_Ghana_2015.tif               -> GHSPOP_Tanzania_2015.tif
  LandCoverType_Ghana_2016.tif        -> LandCoverType_Tanzania_2010.tif
  WAPOR_Ghana_2016.tif                -> WAPOR_Tanzania_2010.tif
  LST_Night_2014_Monthly_Ghana.tif    -> LST_Night_2010_Monthly_Tanzania.tif
  LST_Night_2014_Annual_Mean_Ghana.tif-> LST_Night_2010_Annual_Mean_Tanzania.tif

Two honest limits:
  * WAPOR is not in the public Earth Engine catalogue. MODIS NPP is substituted
    and the file is named WAPOR_* only so the column name matches; that is a
    DIFFERENT product and is recorded as such in FINDINGS.md. If exact parity
    matters, re-export from the original script.
  * FLDAS and ESA WorldCereal band names are reproduced from the source
    collections, but the exact band SET the original export used is not
    knowable from here, so parity is by construction rather than verified.
"""
import io
import os
import time
import zipfile

import ee
import requests

OUT_DIR = "data/Tanzania_GEE_rasters"
YEAR = 2010
SCALE = 1000

# ee.Initialize() must run before any ee.Geometry / ee.Image is constructed --
# building one at import time raises "client library not initialized".
ee.Initialize()
os.makedirs(OUT_DIR, exist_ok=True)

REGION = ee.Geometry.Rectangle([29.30, -11.80, 40.50, -0.95],
                               proj="EPSG:4326", geodesic=False)


def monthly_stack(cid, band, year, scale_factor=1.0, prefix=None):
    """12 monthly means, band-named like the other countries' exports."""
    imgs = []
    for m in range(1, 13):
        start = ee.Date.fromYMD(year, m, 1)
        end = start.advance(1, "month")
        im = (ee.ImageCollection(cid).filterDate(start, end).select([band]).mean())
        if scale_factor != 1.0:
            im = im.multiply(scale_factor)
        nm = f"{year}_{m:02d}_Mean" if prefix is None else f"{prefix}{year}{m:02d}_{band}"
        imgs.append(im.rename([nm]))
    return ee.Image.cat(imgs)


def build_esa():
    # ESA WorldCereal active-cropland classifications, 2021 (the only year the
    # product exists); band names mirror the Ghana file's structure.
    base = "ESA/WorldCereal/2021/MODELS/v100"
    parts = []
    for season, label in [("tc-maize-main", "TC-MAIZE-MAIN"),
                          ("tc-maize-second", "TC-MAIZE-SECOND"),
                          ("tc-wintercereals", "TC-WINTERCEREALS")]:
        ic = (ee.ImageCollection(f"{base}/{season}/activecropland")
              .select(["classification"]))
        parts.append(ic.mosaic().rename([f"{label}_ACTIVECROPLAND_classification"]))
    return ee.Image.cat(parts)


def build_fldas_annual():
    bands = ["Evap_tavg", "LWdown_f_tavg", "Lwnet_tavg", "Psurf_f_tavg",
             "Qair_f_tavg", "Qg_tavg", "Qh_tavg", "Qle_tavg", "Qs_tavg",
             "Qsb_tavg", "RadT_tavg", "Rainf_f_tavg", "SWE_inst",
             "SWdown_f_tavg", "SoilMoi00_10cm_tavg", "SoilTemp00_10cm_tavg",
             "Swnet_tavg", "Tair_f_tavg", "Wind_f_tavg"]
    ic = (ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001")
          .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01"))
    have = ic.first().bandNames().getInfo()
    use = [b for b in bands if b in have]
    return ic.select(use).mean().rename(use)


JOBS = [
    ("ESA_WorldCereal_2021_Tanzania", build_esa),
    ("FLDAS_2010_Annual_Mean_Tanzania", build_fldas_annual),
    ("GHSBUILTS_Tanzania_2015",
     lambda: ee.Image("JRC/GHSL/P2023A/GHS_BUILT_S/2015")
     .select(["built_surface"]).rename(["built_surface"])),
    ("GHSPOP_Tanzania_2015",
     lambda: ee.Image("JRC/GHSL/P2023A/GHS_POP/2015")
     .select(["population_count"]).rename(["population_count"])),
    ("LandCoverType_Tanzania_2010",
     lambda: ee.Image(ee.ImageCollection("MODIS/061/MCD12Q1")
                      .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01").first())
     .select(["LC_Type1", "LC_Type2", "LC_Type3", "LC_Type4", "LC_Type5"])
     .rename([f"{YEAR}_01_01_LC_Type{i}" for i in range(1, 6)])),
    # NOT WAPOR: MODIS NPP standing in so the derived column name matches.
    ("WAPOR_Tanzania_2010",
     lambda: (ee.ImageCollection("MODIS/061/MOD17A3HGF")
              .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01")
              .select(["Npp"]).mean().multiply(0.0001).rename(["L1_NPP_D"]))),
    ("LST_Night_2010_Monthly_Tanzania",
     lambda: monthly_stack("MODIS/061/MOD11A2", "LST_Night_1km", YEAR, 0.02)),
    ("LST_Night_2010_Annual_Mean_Tanzania",
     lambda: (ee.ImageCollection("MODIS/061/MOD11A2")
              .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01")
              .select(["LST_Night_1km"]).mean().multiply(0.02).rename(["Mean"]))),
]


def fetch(name, builder):
    dest = os.path.join(OUT_DIR, name + ".tif")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"  cached  {name}", flush=True)
        return True
    for attempt in range(3):
        try:
            img = builder().toFloat()
            url = img.getDownloadURL({"region": REGION, "scale": SCALE,
                                      "crs": "EPSG:4326", "format": "GEO_TIFF",
                                      "filePerBand": False})
            r = requests.get(url, timeout=900)
            r.raise_for_status()
            payload = r.content
            if payload[:2] == b"PK":
                with zipfile.ZipFile(io.BytesIO(payload)) as z:
                    inner = [n for n in z.namelist() if n.endswith(".tif")]
                    payload = z.read(inner[0])
            with open(dest, "wb") as fh:
                fh.write(payload)
            print(f"  OK      {name}  ({len(payload)/1e6:.1f} MB)", flush=True)
            return True
        except Exception as exc:
            if attempt == 2:
                print(f"  FAILED  {name}: {str(exc)[:180]}", flush=True)
                return False
            time.sleep(5 * (attempt + 1))
    return False


ok = sum(fetch(n, b) for n, b in JOBS)
print(f"\n{ok}/{len(JOBS)} additional rasters in {OUT_DIR}", flush=True)
