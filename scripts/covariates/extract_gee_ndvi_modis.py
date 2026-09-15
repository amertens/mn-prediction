"""
scripts/covariates/extract_gee_ndvi_modis.py   [NV-01, 2026-09-15]

MODIS NDVI AT THE SURVEY, FOR EVERY ADMIN-2 POLYGON

Replaces the eight NOAA AVHRR CDR v5 year columns (ndvi_y2011..ndvi_y2018):
that asset is declared available only to 2013, its post-2018 values collapse
country-dependently (exclusions.csv), and a column per calendar year is a
lag of up to 8 years from the survey. MOD13Q1 (250 m, 16-day, 2000-present)
covers every survey year, and the same collection already supplies evi_t0.

Six columns, all from MOD13Q1 NDVI (scale 0.0001) with SummaryQA <= 1
(good or marginal; snow/ice and cloud masked):
  ndvi_modis_t0        mean NDVI over the survey calendar year
  ndvi_modis_win       mean NDVI over the 12 months ending in the fieldwork
                       median month (metadata/survey_years.csv)
  ndvi_modis_peak_t0   maximum 16-day NDVI in the survey year (peak greenness)
  ndvi_modis_amp_t0    peak minus minimum in the survey year (seasonal amplitude)
  ndvi_modis_clim      mean of the 2001-2020 annual means (long-run greenness)
  ndvi_modis_anom_t0   (survey-year mean - climatology mean) / interannual s.d.
                       of the 2001-2020 annual means, per pixel: a z-score of
                       the survey year against its own history
Each is reduced two ways, area-weighted (_aw) and WorldPop-weighted at the
survey year (_pw, as in script 54); the append step takes the _aw columns
under the plain names and leaves both in the block file.

  python scripts/covariates/extract_gee_ndvi_modis.py
-> data/covariates/harmonized/gee_ndvi_modis_admin2.csv
   data/covariates/harmonized/predictors_admin2_ndvi_modis.csv (+ _metadata.csv)
"""
import csv, json, os, sys, time
import ee

sys.path.insert(0, "C:/Users/andre/OneDrive/Documents/mn-prediction/scripts/protocol_v2")
from survey_years import SURVEY_YEAR, survey_years_table

PROJECT = "mn-prediction-420517"
ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
GEOJSON = ROOT + "data/external_cache/gee_geoms/admin2_simplified.geojson"
OUTDIR = ROOT + "data/covariates/harmonized/"
ISO = {"Gambia": "GMB", "Ghana": "GHA", "Malawi": "MWI", "SierraLeone": "SLE"}
BATCH = int(os.environ.get("FE_BATCH", "40"))
CLIM_Y0, CLIM_Y1 = 2001, 2020
SCALE = 250
VARS = ["ndvi_modis_t0", "ndvi_modis_win", "ndvi_modis_peak_t0", "ndvi_modis_amp_t0", "ndvi_modis_clim", "ndvi_modis_anom_t0"]

ee.Initialize(project=PROJECT)


def fieldwork_median(country):
    for r in survey_years_table():
        if r["country"] == country and r.get("respondent_median_date"):
            return int(r["respondent_median_date"][:4]), int(r["respondent_median_date"][5:7])
    return SURVEY_YEAR[country], 6


def getinfo_retry(obj, what):
    for attempt in range(5):
        try:
            return obj.getInfo()
        except Exception as e:
            print(f"   retry {attempt + 1} ({what}): {str(e)[:140]}", flush=True)
            time.sleep(20 * (attempt + 1))
    raise RuntimeError("Earth Engine call failed five times: " + what)


def pop_weight(country, year):
    iso = ISO[country]
    col = ee.ImageCollection("WorldPop/GP/100m/pop").filter(ee.Filter.eq("country", iso))
    years = sorted(set(col.aggregate_array("year").getInfo()))
    pick = min(years, key=lambda y: abs(y - year))
    pop = ee.Image(col.filter(ee.Filter.eq("year", pick)).first()).select("population")
    w = (pop.unmask(0).reduceResolution(reducer=ee.Reducer.sum().unweighted(), maxPixels=1024)
         .reproject(crs="EPSG:4326", scale=250)).rename("w")
    return w, pick


def masked_ndvi():
    """MOD13Q1 NDVI in natural units, good + marginal quality only."""
    def prep(img):
        qa = img.select("SummaryQA")
        return img.select("NDVI").multiply(0.0001).updateMask(qa.lte(1)).rename("ndvi").copyProperties(img, ["system:time_start"])
    return ee.ImageCollection("MODIS/061/MOD13Q1").map(prep)


def ndvi_image(country):
    sy = SURVEY_YEAR[country]
    my, mm = fieldwork_median(country)
    col = masked_ndvi()
    year = col.filter(ee.Filter.calendarRange(sy, sy, "year"))
    t0 = year.mean().rename("ndvi_modis_t0")
    peak = year.max().rename("ndvi_modis_peak_t0")
    amp = year.max().subtract(year.min()).rename("ndvi_modis_amp_t0")
    # 12 months ending in the fieldwork median month
    end = ee.Date.fromYMD(my, mm, 1).advance(1, "month")
    win = col.filterDate(end.advance(-12, "month"), end).mean().rename("ndvi_modis_win")
    annual = ee.ImageCollection([col.filter(ee.Filter.calendarRange(y, y, "year")).mean().rename("a") for y in range(CLIM_Y0, CLIM_Y1 + 1)])
    clim_mean = annual.mean().rename("ndvi_modis_clim")
    clim_sd = annual.reduce(ee.Reducer.stdDev()).rename("sd")
    anom = t0.subtract(clim_mean).divide(clim_sd).rename("ndvi_modis_anom_t0")
    img = ee.Image.cat([t0, win, peak, amp, clim_mean, anom])
    return img, (sy, my, mm)


def zonal(img, feats, w):
    bands = VARS
    fc = ee.FeatureCollection(feats)
    aw = img.reduceRegions(collection=fc, reducer=ee.Reducer.mean(), scale=SCALE, tileScale=4)
    rows = []
    for ft in getinfo_retry(aw, "area mean")["features"]:
        p = ft["properties"]
        rows.append({"country": p["country"], "Admin1": p["Admin1"], "Admin2": p["Admin2"], **{b + "_aw": p.get(b) for b in bands}})
    num = img.multiply(w).rename([b + "_num" for b in bands])
    den = ee.Image.cat([w.updateMask(img.select(b).mask()).rename(b + "_den") for b in bands])
    pw = num.addBands(den).reduceRegions(collection=fc, reducer=ee.Reducer.sum(), scale=250, tileScale=4)
    for r, ft in zip(rows, getinfo_retry(pw, "population-weighted sums")["features"]):
        p = ft["properties"]
        for b in bands:
            n, d = p.get(b + "_num"), p.get(b + "_den")
            r[b + "_pw"] = (n / d) if (n is not None and d not in (None, 0)) else None
    return rows


def main():
    gj = json.load(open(GEOJSON, encoding="utf-8"))
    feats_by_country = {}
    for ft in gj["features"]:
        feats_by_country.setdefault(ft["properties"]["country"], []).append(ft)
    rows = []
    t_start = time.time()
    for country in ("Gambia", "Ghana", "Malawi", "SierraLeone"):
        feats = feats_by_country[country]
        img, (sy, my, mm) = ndvi_image(country)
        w, pick = pop_weight(country, sy)
        print(f"  {country}: {len(feats)} polygons, survey year {sy}, window ends {my}-{mm:02d}, WorldPop {pick}", flush=True)
        for i in range(0, len(feats), BATCH):
            chunk = feats[i:i + BATCH]
            efeats = [ee.Feature(ee.Geometry(ft["geometry"]), {k: ft["properties"][k] for k in ("country", "Admin1", "Admin2")}) for ft in chunk]
            rows.extend(zonal(img, efeats, w))
            print(f"    {min(i + BATCH, len(feats))}/{len(feats)} done, {time.time() - t_start:.0f} s", flush=True)
    cols = ["country", "Admin1", "Admin2"] + [v + s for v in VARS for s in ("_aw", "_pw")]
    raw = OUTDIR + "gee_ndvi_modis_admin2.csv"
    with open(raw, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols); wr.writeheader()
        for r in rows:
            wr.writerow({k: ("" if r.get(k) is None else (round(r[k], 6) if isinstance(r[k], float) else r[k])) for k in cols})
    print(f"wrote {raw}: {len(rows)} rows", flush=True)

    # the block for the shared set: area-weighted values under the plain names
    blk = OUTDIR + "predictors_admin2_ndvi_modis.csv"
    with open(blk, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=["country", "Admin1", "Admin2"] + VARS); wr.writeheader()
        for r in rows:
            wr.writerow({"country": r["country"], "Admin1": r["Admin1"], "Admin2": r["Admin2"],
                         **{v: ("" if r.get(v + "_aw") is None else round(r[v + "_aw"], 6)) for v in VARS}})
    desc = {
        "ndvi_modis_t0": "Mean MOD13Q1 NDVI over the survey calendar year",
        "ndvi_modis_win": "Mean MOD13Q1 NDVI over the 12 months ending in the fieldwork median month",
        "ndvi_modis_peak_t0": "Maximum 16-day MOD13Q1 NDVI in the survey year (peak greenness)",
        "ndvi_modis_amp_t0": "Peak minus minimum 16-day NDVI in the survey year (seasonal amplitude)",
        "ndvi_modis_clim": "Mean of the 2001-2020 annual mean NDVI (long-run greenness)",
        "ndvi_modis_anom_t0": "Survey-year mean NDVI minus the 2001-2020 climatology, divided by the interannual s.d. of the annual means (per-pixel z-score, then area-weighted)",
    }
    countries = sorted(set(r["country"] for r in rows))
    with open(OUTDIR + "predictors_admin2_ndvi_modis_metadata.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=["column", "source", "domain", "subnational", "assumption", "n_countries", "countries", "completeness"]); wr.writeheader()
        for v in VARS:
            fin = [r for r in rows if r.get(v + "_aw") is not None]
            cs = sorted(set(r["country"] for r in fin))
            wr.writerow({"column": v, "source": "GEE", "domain": "Ecosystem productivity/greenness", "subnational": "TRUE",
                         "assumption": desc[v] + ". MODIS/061/MOD13Q1 at 250 m, SummaryQA <= 1, area-weighted zonal mean over the GADM Admin-2 polygon (scripts/covariates/extract_gee_ndvi_modis.py; the WorldPop-weighted twin is in gee_ndvi_modis_admin2.csv). Replaces the AVHRR CDR v5 year columns ndvi_y2011-2018 (NV-01).",
                         "n_countries": len(cs), "countries": ";".join(cs), "completeness": round(len(fin) / len(rows), 3)})
    print("wrote the ndvi_modis block and metadata", flush=True)


if __name__ == "__main__":
    main()
