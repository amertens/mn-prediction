"""
scripts/covariates/gee_water_coast_distance.py

Distance to surface water and to the coast, from Earth Engine, for the Admin-2
spine and the survey clusters.

  wdist_perm_km   distance (km) to the nearest pixel that is water in >= 50% of
                  the 1984-2021 JRC Global Surface Water observations
  wdist_any_km    distance (km) to water present in >= 10% of observations
                  (seasonal rivers, floodplains, reservoirs)
  wdist_coast_km  distance (km) to the ocean (USDOS LSIB 2017 land polygons)

Why: fish and aquatic foods (iron, zinc, B12, vitamin A), and the classic
iodine geography (inland and mountainous soils are iodine-poor). Both are
missing from the vocabulary; the nearest existing columns are the water
land-cover fractions, which say how much water a district contains, not how
far its people are from it.

Method: binary mask -> fastDistanceTransform in Web Mercator at 500 m (water)
or 2 km (coast), multiplied by cos(latitude) to undo the Mercator stretch.
Saturates at 1024 pixels (512 km for water, 2048 km for coast), which no
district in the four countries reaches. Admin-2: mean and minimum over the
polygon; clusters: mean over the 2 km / 5 km buffer.

Run with the reticulate venv python (earthengine-api 0.1.370 pinned):
  ".virtualenvs/r-reticulate/Scripts/python.exe" scripts/covariates/gee_water_coast_distance.py
-> data/covariates/harmonized/predictors_admin2_water_distance.csv
   data/covariates/cluster/predictors_cluster_water_distance.csv
"""
import csv, json, math, os, sys, time
import ee

PROJECT = "mn-prediction-420517"
ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
GEOM = ROOT + "data/external_cache/gee_geoms/"
OUT_A2 = ROOT + "data/covariates/harmonized/predictors_admin2_water_distance.csv"
OT = os.environ.get("CL_OUT_TAG", "")
OUT_CL = ROOT + "data/covariates/cluster/predictors_cluster_water_distance" + OT + ".csv"
BATCH = 25

ee.Initialize(project=PROJECT)


def dist_km(mask, scale):
    """km from every pixel to the nearest mask == 1 pixel, computed on a Web Mercator grid at `scale` m."""
    m = mask.reproject(crs="EPSG:3857", scale=scale)
    d = m.fastDistanceTransform(1024, "pixels", "squared_euclidean").sqrt().multiply(scale).divide(1000.0)
    coslat = ee.Image.pixelLonLat().select("latitude").multiply(math.pi / 180.0).cos()
    return d.multiply(coslat)


occ = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").unmask(0)
water_perm = occ.gte(50)
water_any = occ.gte(10)
land = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
ocean = ee.Image.constant(1).paint(land, 0)          # 1 = not land
# One image per layer, each reduced at its own scale. Reducing a composite of
# differently projected distance transforms at one scale silently recomputed
# the coast transform on the reduce grid and returned distances of thousands
# of km for inland districts (found 2026-09-07); per-layer reduction matches
# point checks (Tamale 415 km, Lilongwe 555 km, Banjul 0-2 km).
LAYERS = {
    "wdist_perm_km":  (dist_km(water_perm, 500).rename("wdist_perm_km"), 500),
    "wdist_any_km":   (dist_km(water_any, 500).rename("wdist_any_km"), 500),
    "wdist_coast_km": (dist_km(ocean, 2000).rename("wdist_coast_km"), 2000),
}
BANDS = list(LAYERS)


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)["features"]


def reduce_batch(feats, reducer, scale_unused, keep):
    """Reduce every layer at its own scale; merge the per-layer properties by feature order."""
    fc = ee.FeatureCollection(feats)
    merged = None
    for band, (img, scale) in LAYERS.items():
        res = img.reduceRegions(collection=fc, reducer=reducer, scale=scale, tileScale=4)
        out = None
        for attempt in range(4):
            try:
                out = res.getInfo()["features"]; break
            except Exception as e:                   # transient EE errors: back off and retry
                print("   retry", attempt + 1, str(e)[:120], flush=True)
                time.sleep(15 * (attempt + 1))
        if out is None:
            raise RuntimeError("Earth Engine reduce failed four times")
        # a single-band image reduced with mean/min yields properties named "mean"/"min" (no band prefix): re-key them
        for ft in out:
            p = ft["properties"]
            for k in ("mean", "min"):
                if k in p:
                    p[band + "_" + k if COMBINED else band] = p.pop(k)
        if merged is None:
            merged = out
        else:
            for a, b in zip(merged, out):
                a["properties"].update({k: v for k, v in b["properties"].items() if k.startswith("wdist_")})
    return merged


COMBINED = True   # True: mean+min reducer (Admin-2); False: single mean (clusters)


def run(geojson, geom_fn, id_keys, reducer, scale, out_path, expected_cols):
    global COMBINED
    COMBINED = len(expected_cols) > len(LAYERS)
    feats_in = load(geojson)
    rows = []
    t0 = time.time()
    for i in range(0, len(feats_in), BATCH):
        chunk = feats_in[i:i + BATCH]
        feats = [ee.Feature(geom_fn(ft), {k: ft["properties"].get(k) for k in id_keys}) for ft in chunk]
        out = reduce_batch(feats, reducer, scale, id_keys)
        for ft in out:
            p = ft["properties"]
            rows.append({k: p.get(k) for k in id_keys + expected_cols})
        print(f"   {min(i + BATCH, len(feats_in))}/{len(feats_in)} done, {time.time() - t0:.0f} s", flush=True)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=id_keys + expected_cols)
        w.writeheader(); w.writerows(rows)
    print("wrote", out_path, len(rows), "rows", flush=True)
    return rows


if __name__ == "__main__":
    # Admin-2: mean and min of each distance over the polygon
    red_a2 = ee.Reducer.mean().combine(reducer2=ee.Reducer.min(), sharedInputs=True)
    cols_a2 = [b + "_mean" for b in BANDS] + [b + "_min" for b in BANDS]
    if not os.environ.get("CL_SKIP_ADMIN2"):
        print("Admin-2 polygons", flush=True)
        run(GEOM + "admin2_simplified.geojson", lambda ft: ee.Geometry(ft["geometry"]),
            ["country", "Admin1", "Admin2"], red_a2, 500, OUT_A2, cols_a2)
    # clusters: mean over the buffer (2 km urban / 5 km rural)
    print("clusters", flush=True)
    run(GEOM + "clusters" + OT + ".geojson",
        lambda ft: ee.Geometry.Point(ft["geometry"]["coordinates"]).buffer(1000.0 * float(ft["properties"].get("radius_km", 5))),
        ["country", "cluster", "radius_km"], ee.Reducer.mean(), 500, OUT_CL, BANDS)
    print("DONE", flush=True)
