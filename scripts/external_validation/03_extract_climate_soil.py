"""
scripts/external_validation/03_extract_climate_soil.py               [XV-01]

CLIMATE AND SOIL FOR THE FOUR PANEL COUNTRIES AND THE FOUR NEW ONES,
FROM ONE EXTRACTOR

The pre-registered parsimonious candidate is the two-domain climate + soil
index (P2). Both domains are global/continental rasters with no survey
linkage, which is why a new country needs no data-sharing agreement -- and why
the new countries must be in Africa: the soil half is iSDAsoil, Africa only.

WHY THIS RE-EXTRACTS THE PANEL COUNTRIES TOO. The existing
predictors_admin2_shared.csv soil/climate columns come from an older analyst
export whose scaling and back-transforms are not fully recorded. The pooled
model matches covariates BY NAME and treats an absent or differently-scaled
name as if it were comparable (gee_legacy_name_vocabulary; WS-G found a
country silently zeroing a whole block). Rather than try to match that export
column by column and hope, both sides of the comparison are extracted here by
the same code, so training and held-out values are commensurate by
construction. The cost is that this is the index re-derived on a common
vocabulary, not the literal pre-registered column set; the column count is
reported so the two can be told apart.

Area-weighted only. Population-weighting the PREDICTORS costs 0.05-0.15
Spearman (AG-01), so the `_pw` half of script 54 is deliberately not computed.

  python scripts/external_validation/03_extract_climate_soil.py [blocks]
      blocks: clim,isda  (default both)
-> data/external_validation/gee_clim_admin2.csv
   data/external_validation/gee_isda_admin2.csv
Needs the miniconda interpreter (earthengine-api): C:/Users/andre/miniconda3/python.exe
"""
import csv
import json
import os
import sys
import time
import ee

PROJECT = "mn-prediction-420517"
ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
GEOJSONS = [
    ROOT + "data/external_cache/gee_geoms/admin2_simplified.geojson",
    ROOT + "data/external_validation/admin2_new_countries.geojson",
    ROOT + "data/external_validation/admin2_offcontinent.geojson",
]
OUTDIR = ROOT + "data/external_validation/"
BATCH = int(os.environ.get("XV_BATCH", "60"))
CLIM_Y0, CLIM_Y1 = 1991, 2020
LST_Y0, LST_Y1 = 2003, 2020

ee.Initialize(project=PROJECT)


def survey_years():
    out = {}
    with open(ROOT + "metadata/survey_years.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[r["country"]] = int(r["survey_year"])
    return out


SY = survey_years()

TC_SCALE = {"pr": 1.0, "tmmx": 0.1, "tmmn": 0.1, "pet": 0.1, "def": 0.1,
            "aet": 0.1, "soil": 0.1, "vpd": 0.01, "srad": 0.1}


def getinfo_retry(obj, what):
    for attempt in range(5):
        try:
            return obj.getInfo()
        except Exception as e:
            print("   retry %d (%s): %s" % (attempt + 1, what, str(e)[:140]), flush=True)
            time.sleep(20 * (attempt + 1))
    raise RuntimeError("Earth Engine call failed five times: " + what)


def monthly_climatology(col, band, factor, y0, y1, prefix):
    imgs = []
    for m in range(1, 13):
        sub = (col.filter(ee.Filter.calendarRange(y0, y1, "year"))
                  .filter(ee.Filter.calendarRange(m, m, "month")).select(band))
        imgs.append(sub.mean().multiply(factor).rename("%s_m%02d" % (prefix, m)))
    return ee.Image.cat(imgs)


def climatology_image(country):
    """Identical to scripts/protocol_v2/54_extract_climatology_terrain_soil.py's
    climatology block, minus the fieldwork-window term (fieldwork months are not
    reported in VMNIS for Zambia or Sudan, so a window would be invented)."""
    tc = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    pr = monthly_climatology(tc, "pr", 1.0, CLIM_Y0, CLIM_Y1, "pr")
    tmax = monthly_climatology(tc, "tmmx", 0.1, CLIM_Y0, CLIM_Y1, "tmax")
    ann = []
    for b, nm in (("tmmn", "tmin_ann"), ("pet", "pet_ann"), ("def", "def_ann"),
                  ("aet", "aet_ann"), ("soil", "soilm_ann"), ("vpd", "vpd_ann"),
                  ("srad", "srad_ann")):
        sub = tc.filter(ee.Filter.calendarRange(CLIM_Y0, CLIM_Y1, "year")).select(b)
        ann.append(sub.mean().multiply(TC_SCALE[b]).rename(nm))
    yearly = ee.ImageCollection([
        tc.filter(ee.Filter.calendarRange(y, y, "year")).select("pr").sum().rename("pr_year")
        for y in range(CLIM_Y0, CLIM_Y1 + 1)])
    parts = [pr, tmax] + ann + [
        yearly.mean().rename("pr_ann_mean"),
        yearly.reduce(ee.Reducer.stdDev()).rename("pr_ann_sd")]
    sy = SY[country]
    parts.append(tc.filter(ee.Filter.calendarRange(sy, sy, "year"))
                 .select("pr").sum().rename("pr_sy"))
    parts.append(tc.filter(ee.Filter.calendarRange(sy, sy, "year"))
                 .select("tmmx").mean().multiply(0.1).rename("tmax_sy"))
    lst = ee.ImageCollection("MODIS/061/MOD11A2")
    parts.append(monthly_climatology(lst, "LST_Day_1km", 0.02, LST_Y0, LST_Y1, "lstd")
                 .subtract(273.15))
    parts.append(monthly_climatology(lst, "LST_Night_1km", 0.02, LST_Y0, LST_Y1, "lstn")
                 .subtract(273.15))
    return ee.Image.cat(parts).toFloat().resample("bilinear")


ISDA = {
    "ph": ("ph", "div10"), "clay": ("clay_content", "pct"),
    "sand": ("sand_content", "pct"), "silt": ("silt_content", "pct"),
    "bd": ("bulk_density", "div100"), "oc": ("carbon_organic", "explog"),
    "ntot": ("nitrogen_total", "explog100"),
    "cec": ("cation_exchange_capacity", "explog"),
    "zn": ("zinc_extractable", "explog"), "fe": ("iron_extractable", "explog"),
    "ca": ("calcium_extractable", "explog"),
    "mg": ("magnesium_extractable", "explog"),
    "k": ("potassium_extractable", "explog"),
    "p": ("phosphorus_extractable", "explog"),
    "s": ("sulphur_extractable", "explog"),
    "al": ("aluminium_extractable", "explog"),
}


def back(img, how):
    if how == "pct":
        return img
    if how == "div10":
        return img.divide(10)
    if how == "div100":
        return img.divide(100)
    if how == "explog":
        return img.divide(10).exp().subtract(1)
    if how == "explog100":
        return img.divide(100).exp().subtract(1)
    raise ValueError(how)


def isda_image(country):
    """Both the mean and the stdev band, both depths -- the stdev bands carry
    within-district soil heterogeneity and are part of the Soil domain."""
    parts = []
    for short, (suffix, how) in ISDA.items():
        im = ee.Image("ISDASOIL/Africa/v1/" + suffix)
        for stat in ("mean", "stdev"):
            for depth in ("0_20", "20_50"):
                parts.append(back(im.select(stat + "_" + depth), how)
                             .rename("%s_%s_%s" % (short, stat, depth)))
    return ee.Image.cat(parts).toFloat()


def zonal(img, feats, scale, tile=4):
    bands = img.bandNames().getInfo()
    fc = ee.FeatureCollection(feats)
    out = getinfo_retry(
        img.reduceRegions(collection=fc, reducer=ee.Reducer.mean(),
                          scale=scale, tileScale=tile),
        "area mean")["features"]
    rows = []
    for ft in out:
        p = ft["properties"]
        row = {"country": p["country"], "Admin1": p["Admin1"], "Admin2": p["Admin2"]}
        for b in bands:
            row[b] = p.get(b)
        rows.append(row)
    return rows


# ── the substituted soil block: global SoilGrids v2.0 (ISRIC) ───────────────
# iSDAsoil is Africa-only, so the off-continent arm cannot use it. SoilGrids is
# global but carries only physical/chemical properties -- it has NO
# plant-available micronutrients (Zn, Fe, Ca, Mg, P, K, S), which are the
# mechanistically interesting half of iSDA for a micronutrient outcome. The
# substitution is therefore NOT neutral, which is why this block is extracted
# for the African countries too: running both soil blocks on Africa prices the
# substitution, so an off-continent result can be read against it instead of
# being confounded with it.
# Integer -> natural units per the SoilGrids v2.0 conversion table.
SGRID = {"phh2o": 10.0, "soc": 10.0, "nitrogen": 100.0, "cec": 10.0,
         "clay": 10.0, "sand": 10.0, "silt": 10.0, "bdod": 100.0,
         "cfvo": 10.0, "ocd": 10.0}
SGRID_DEPTHS = ("0-5cm", "5-15cm", "15-30cm", "30-60cm")


def sgrid_image(country):
    parts = []
    for prop, div in SGRID.items():
        im = ee.Image("projects/soilgrids-isric/" + prop + "_mean")
        for d in SGRID_DEPTHS:
            parts.append(im.select("%s_%s_mean" % (prop, d)).divide(div)
                         .rename("sg_%s_%s" % (prop, d.replace("-", "_"))))
    return ee.Image.cat(parts).toFloat()


def load_feats():
    feats = {}
    for path in GEOJSONS:
        if not os.path.exists(path):
            print("  (missing, skipped) " + path, flush=True)
            continue
        gj = json.load(open(path, encoding="utf-8"))
        for f in gj["features"]:
            feats.setdefault(f["properties"]["country"], []).append(f)
    return feats


def write(path, rows):
    keys = set()
    for r in rows:
        keys.update(r.keys())
    keys -= {"country", "Admin1", "Admin2"}
    cols = ["country", "Admin1", "Admin2"] + sorted(keys)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        for r in rows:
            out = {}
            for k in cols:
                v = r.get(k)
                out[k] = "" if v in (None, "") else (round(v, 6) if isinstance(v, float) else v)
            wr.writerow(out)


def run_block(name, image_fn, scale, feats_by_country):
    out_path = OUTDIR + "gee_" + name + "_admin2.csv"
    done = set()
    rows = []
    if os.path.exists(out_path):   # resume: an interrupted run keeps its work
        with open(out_path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                rows.append(r)
                done.add((r["country"], r["Admin1"], r["Admin2"]))
        print("  resuming %s: %d polygons already extracted" % (name, len(done)), flush=True)
    t0 = time.time()
    for country, feats in feats_by_country.items():
        todo = [f for f in feats
                if (country, f["properties"]["Admin1"],
                    f["properties"]["Admin2"]) not in done]
        if not todo:
            print("  %s: complete" % country, flush=True)
            continue
        img = image_fn(country)
        print("  %s: %d polygons to do (%d bands)"
              % (country, len(todo), img.bandNames().size().getInfo()), flush=True)
        for i in range(0, len(todo), BATCH):
            chunk = todo[i:i + BATCH]
            ef = []
            for ft in chunk:
                props = {k: ft["properties"][k] for k in ("country", "Admin1", "Admin2")}
                ef.append(ee.Feature(ee.Geometry(ft["geometry"]), props))
            rows.extend(zonal(img, ef, scale))
            print("    %d/%d %s, %.0fs"
                  % (min(i + BATCH, len(todo)), len(todo), country, time.time() - t0),
                  flush=True)
            write(out_path, rows)
    write(out_path, rows)
    print("wrote %s: %d rows" % (out_path, len(rows)), flush=True)


def main(blocks):
    feats = load_feats()
    total = sum(len(v) for v in feats.values())
    print("%d polygons over %d countries; blocks: %s"
          % (total, len(feats), ", ".join(blocks)), flush=True)
    if "clim" in blocks:
        print("[climatology]", flush=True)
        run_block("clim", climatology_image, 1000, feats)
    if "isda" in blocks:   # Africa only; the off-continent countries are skipped
        print("[isda]", flush=True)
        afr = {k: v for k, v in feats.items() if k not in ("Pakistan", "India")}
        run_block("isda", isda_image, 250, afr)
    if "sgrid" in blocks:
        print("[soilgrids]", flush=True)
        run_block("sgrid", sgrid_image, 250, feats)
    print("DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(",") if len(sys.argv) > 1 else ["clim", "isda", "sgrid"])
