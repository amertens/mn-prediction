"""
explore/scripts/10_gee_cluster_monthly.py   [probe TM-01, phase 2 - GATED]

MONTHLY COVARIATE STACKS AT SURVEY CLUSTERS, 36 MONTHS BEFORE EACH BLOOD DRAW

Phase 1 (09_temporal_phase1.R) can only test seasonal PHASE, because the
monthly rasters on disk cover parts of 2014-15 and the four surveys run
2013-2018. The hypothesis phase 1 cannot reach is the one that matters:

    ferritin and retinol integrate MONTHS of intake, so the quality of the
    PREVIOUS GROWING SEASON (lags 9-15 months) should predict status better
    than the long-run climatology of the pixel, and better than the 3-month
    window the cluster table already carries.

DESIGN. One reduceRegions per layer over a multiband monthly image, not one
call per month: a monthly ImageCollection over the full 2010-2019 span is
turned into a single ~120-band image with toBands(), reduced once over all
323 cluster buffers, and the lag alignment is done afterwards in pandas
against each cluster's own fieldwork date. That is 5 Earth Engine calls
total rather than 5 x 120.

Buffers follow the cluster track's existing convention: 2 km urban / 5 km
rural, urban taken from the GHSL SMOD flag already in predictors_cluster.csv.
DHS displaces cluster GPS by up to 2 km urban / 5 km rural, so a finer radius
than this would be measuring the displacement, not the place.

  python explore/scripts/10_gee_cluster_monthly.py
  EXP_GEE_LAYERS=chirps,ndvi   restrict to a subset
-> explore/out/10_gee_cluster_monthly.csv   (long: cluster x layer x lag)
"""
import csv
import os
import sys
import time
from collections import defaultdict

import ee

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
PROJECT = "mn-prediction-420517"
TARGETS = ROOT + "results/tables/cluster_level/targets_cluster.csv"
PREDS = ROOT + "data/covariates/cluster/predictors_cluster.csv"
OUT = ROOT + "explore/out/10_gee_cluster_monthly.csv"

Y0, Y1 = 2010, 2019          # span covering every survey window minus 36 months
N_LAG = 36                   # months of history kept per cluster
URBAN_KM, RURAL_KM = 2.0, 5.0

LAYERS = {
    # name: (collection, band, scale_m, reducer_over_month)
    "chirps":  ("UCSB-CHG/CHIRPS/DAILY", "precipitation", 5566, "sum"),
    "ndvi":    ("MODIS/061/MOD13Q1", "NDVI", 250, "mean"),
    "lst":     ("MODIS/061/MOD11A2", "LST_Day_1km", 1000, "mean"),
    "soilmoi": ("NASA/FLDAS/NOAH01/C/GL/M/V001", "SoilMoi00_10cm_tavg", 11132, "mean"),
    "evap":    ("NASA/FLDAS/NOAH01/C/GL/M/V001", "Evap_tavg", 11132, "mean"),
}
want = os.environ.get("EXP_GEE_LAYERS", "")
if want:
    LAYERS = {k: v for k, v in LAYERS.items() if k in want.split(",")}

ee.Initialize(project=PROJECT)


def getinfo_retry(obj, what):
    for attempt in range(5):
        try:
            return obj.getInfo()
        except Exception as e:
            print(f"   retry {attempt + 1} ({what}): {str(e)[:160]}", flush=True)
            time.sleep(20 * (attempt + 1))
    raise RuntimeError("Earth Engine call failed five times: " + what)


# ── cluster buffers ─────────────────────────────────────────────────────────
urban_of = {}
with open(PREDS, newline="", encoding="utf-8-sig") as fh:
    for r in csv.DictReader(fh):
        u = (r.get("urban") or r.get("smod") or "").strip()
        try:
            urban_of[(r["country"], r["cluster"])] = float(u) >= 1
        except ValueError:
            urban_of[(r["country"], r["cluster"])] = False

seen, feats, meta = set(), [], []
with open(TARGETS, newline="", encoding="utf-8-sig") as fh:
    for r in csv.DictReader(fh):
        key = (r["country"], r["cluster"])
        if key in seen:
            continue
        try:
            lat, lon = float(r["lat"]), float(r["lon"])
        except (ValueError, KeyError):
            continue
        if not (lat or lon):
            continue
        seen.add(key)
        urb = urban_of.get(key, False)
        rad = (URBAN_KM if urb else RURAL_KM) * 1000.0
        cid = f"{r['country']}__{r['cluster']}"
        feats.append(ee.Feature(ee.Geometry.Point([lon, lat]).buffer(rad), {"cid": cid}))
        meta.append({"cid": cid, "country": r["country"], "cluster": r["cluster"],
                     "urban": int(urb), "buffer_km": rad / 1000.0,
                     "date_med": r.get("date_med", "")})

FC = ee.FeatureCollection(feats)
print(f"{len(feats)} cluster buffers "
      f"({sum(m['urban'] for m in meta)} urban at {URBAN_KM} km, rest at {RURAL_KM} km)")


def monthly_image(coll_id, band, how):
    """A single image whose bands are one value per calendar month, Y0..Y1."""
    coll = ee.ImageCollection(coll_id).select(band)

    def one(ym):
        ym = ee.Number(ym)
        y = ym.divide(12).floor().add(Y0)
        m = ym.mod(12).add(1)
        start = ee.Date.fromYMD(y, m, 1)
        end = start.advance(1, "month")
        sub = coll.filterDate(start, end)
        img = ee.Algorithms.If(how == "sum", sub.sum(), sub.mean())
        return (ee.Image(img)
                .rename(ee.String("m_").cat(y.format("%04d")).cat("_").cat(m.format("%02d")))
                .set("ym", ym))

    n = (Y1 - Y0 + 1) * 12
    return ee.ImageCollection(ee.List.sequence(0, n - 1).map(one)).toBands()


rows = defaultdict(dict)
for name, (coll_id, band, scale, how) in LAYERS.items():
    print(f"-> {name}: {coll_id} [{band}] at {scale} m, monthly {how}", flush=True)
    img = monthly_image(coll_id, band, how)
    fc = img.reduceRegions(collection=FC, reducer=ee.Reducer.mean(), scale=scale)
    res = getinfo_retry(fc, name)
    for f in res["features"]:
        p = f["properties"]
        cid = p["cid"]
        for k, v in p.items():
            if k == "cid" or v is None:
                continue
            # band names arrive as "<index>_m_YYYY_MM"; keep the YYYY_MM tail
            parts = k.split("m_")
            if len(parts) < 2:
                continue
            rows[cid][(name, parts[-1])] = v
    print(f"   got {len(res['features'])} features", flush=True)

# ── align each cluster to its own fieldwork date and write long ─────────────
os.makedirs(os.path.dirname(OUT), exist_ok=True)
n_written = 0
with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["country", "cluster", "urban", "buffer_km", "date_med",
                "layer", "lag_months", "year_month", "value"])
    for m in meta:
        d = m["date_med"]
        if len(d) < 7:
            continue
        fy, fm = int(d[:4]), int(d[5:7])
        base = fy * 12 + (fm - 1)
        for (layer, ym), val in rows.get(m["cid"], {}).items():
            yy, mm = int(ym[:4]), int(ym[5:7])
            lag = base - (yy * 12 + (mm - 1))
            if 0 <= lag < N_LAG:
                w.writerow([m["country"], m["cluster"], m["urban"], m["buffer_km"],
                            d, layer, lag, ym, val])
                n_written += 1

print(f"wrote {n_written:,} rows -> {OUT}")
print("NOTE: lag 0 is the fieldwork month itself; lags 9-15 are the previous "
      "growing season, which is the hypothesis this extraction exists to test.")
