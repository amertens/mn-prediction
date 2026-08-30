"""
sandbox_parsimony/python/fetch_tanzania_rasters3.py

STEP 3, pass 2: the three files pass 1 could not deliver, plus a fix.

  ESA WorldCereal  -- the seasons live in .../MARKERS/v100 as a FLAT collection
                      filtered on the `season` property, not as separate
                      per-season collections. Pass 1 guessed the latter path
                      and every season came back MISSING.
  FLDAS annual     -- 143 MB in one request against a 50 MB cap. Downloaded in
                      band chunks and merged locally with terra.
  LST monthly      -- 90 MB, same treatment (2 chunks of 6 months).
  GHSPOP           -- pass 1 wrote a layer with mean -10.96, i.e. the nodata
                      fill was being averaged as if it were population. Masked
                      to non-negative here and rewritten.

Chunked files are merged by merge_tanzania_chunks.R so the final band set and
filename still match the other countries' exports.
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

ee.Initialize()
os.makedirs(OUT_DIR, exist_ok=True)
REGION = ee.Geometry.Rectangle([29.30, -11.80, 40.50, -0.95],
                               proj="EPSG:4326", geodesic=False)


def download(img, dest, scale=SCALE):
    for attempt in range(3):
        try:
            url = img.toFloat().getDownloadURL(
                {"region": REGION, "scale": scale, "crs": "EPSG:4326",
                 "format": "GEO_TIFF", "filePerBand": False})
            r = requests.get(url, timeout=900)
            r.raise_for_status()
            payload = r.content
            if payload[:2] == b"PK":
                with zipfile.ZipFile(io.BytesIO(payload)) as z:
                    payload = z.read([n for n in z.namelist() if n.endswith(".tif")][0])
            with open(dest, "wb") as fh:
                fh.write(payload)
            print(f"  OK      {os.path.basename(dest)} ({len(payload)/1e6:.1f} MB)", flush=True)
            return True
        except Exception as exc:
            if attempt == 2:
                print(f"  FAILED  {os.path.basename(dest)}: {str(exc)[:170]}", flush=True)
                return False
            time.sleep(5 * (attempt + 1))
    return False


# --- 1. ESA WorldCereal: flat collection, filter on the season property -----
MARKERS = ee.ImageCollection("ESA/WorldCereal/2021/MARKERS/v100")
seasons = MARKERS.aggregate_array("season").distinct().getInfo()
print("ESA seasons available:", seasons, flush=True)
parts, names = [], []
for s in seasons:
    sub = MARKERS.filter(ee.Filter.eq("season", s)).select(["classification"])
    if sub.size().getInfo() == 0:
        continue
    parts.append(sub.mosaic())
    names.append(f"{s.upper()}_ACTIVECROPLAND_classification")
if parts:
    download(ee.Image.cat(parts).rename(names),
             os.path.join(OUT_DIR, "ESA_WorldCereal_2021_Tanzania.tif"))

# --- 2. GHSPOP with the nodata fill masked ---------------------------------
pop = ee.Image("JRC/GHSL/P2023A/GHS_POP/2015").select(["population_count"])
download(pop.updateMask(pop.gte(0)).rename(["population_count"]),
         os.path.join(OUT_DIR, "GHSPOP_Tanzania_2015.tif"))

# --- 3. FLDAS annual mean, in band chunks ----------------------------------
FL = (ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001")
      .filterDate(f"{YEAR}-01-01", f"{YEAR+1}-01-01"))
fl_bands = FL.first().bandNames().getInfo()
print(f"FLDAS bands: {len(fl_bands)}", flush=True)
CH = 5
for i in range(0, len(fl_bands), CH):
    grp = fl_bands[i:i + CH]
    download(FL.select(grp).mean().rename(grp),
             os.path.join(OUT_DIR, f"_chunk_FLDAS_{YEAR}_Annual_Mean_Tanzania_{i//CH}.tif"))

# --- 4. LST night monthly, in two chunks of six ----------------------------
def month_img(m):
    start = ee.Date.fromYMD(YEAR, m, 1)
    return (ee.ImageCollection("MODIS/061/MOD11A2")
            .filterDate(start, start.advance(1, "month"))
            .select(["LST_Night_1km"]).mean().multiply(0.02)
            .rename([f"{YEAR}_{m:02d}_Mean"]))

for k, months in enumerate([range(1, 7), range(7, 13)]):
    download(ee.Image.cat([month_img(m) for m in months]),
             os.path.join(OUT_DIR, f"_chunk_LST_Night_{YEAR}_Monthly_Tanzania_{k}.tif"))

print("\ndone; run merge_tanzania_chunks.R to assemble the chunked files",
      flush=True)
