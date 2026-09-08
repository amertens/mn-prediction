"""
scripts/protocol_v2/54_extract_climatology_terrain_soil.py   [FE-01]

CLIMATOLOGY, TERRAIN AND iSDA SOIL PROPERTIES FOR EVERY ADMIN-2 POLYGON,
AREA-WEIGHTED AND POPULATION-WEIGHTED

Three raw blocks for the climate / soil feature-engineering test (sandbox log
FE-01). Every variable comes out twice: `<v>_aw` is the ordinary area-weighted
zonal mean and `<v>_pw` the WorldPop-weighted mean at the country's survey
year (AU-01 finding 10: respondents live where people are, not where the
district's area is). Everything is computed server-side; the per-polygon
reductions run in batches with retries.

  climatology (TerraClimate 1991-2020 monthly, 4 km; MODIS LST 2003-2020, 1 km)
      pr_m01..pr_m12   monthly precipitation climatology (mm)
      tmax_m01..m12    monthly mean daily maximum temperature (C)
      tmin_ann, pet_ann, def_ann, aet_ann, soilm_ann, vpd_ann, srad_ann   annual means
      pr_ann_mean, pr_ann_sd          mean and s.d. of the 30 annual totals
      pr_sy, tmax_sy                  survey-year annual precipitation / mean tmax
      pr_win                          precipitation over the 12 months ending in the fieldwork median month
      lstd_m01..m12, lstn_m01..m12    monthly day / night LST climatology (C)
  terrain (MERIT Hydro 93 m; Geomorpho90m 90 m)
      elv, hnd, upa (km2), slope, tri, tpi, roughness, vrm, cti, elev_stdev
  isda (iSDAsoil Africa v1, 30 m, reduced at 250 m; back-transformed to natural units)
      <property>_0_20, <property>_20_50 for ph, clay, sand, silt, bd, oc, ntot, cec,
      zn, fe, ca, mg, k, p, s, al

  python scripts/protocol_v2/54_extract_climatology_terrain_soil.py [blocks]   (blocks: clim,terrain,isda; default all)
-> data/covariates/harmonized/gee_climatology_admin2.csv
   data/covariates/harmonized/gee_terrain_admin2.csv
   data/covariates/harmonized/gee_isda_admin2.csv
Script 55 turns these into the engineered add-on blocks.
"""
import csv, json, math, os, sys, time
import ee
from survey_years import SURVEY_YEAR, survey_years_table

PROJECT = "mn-prediction-420517"
ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
GEOJSON = ROOT + "data/external_cache/gee_geoms/admin2_simplified.geojson"
OUTDIR = ROOT + "data/covariates/harmonized/"
ISO = {"Gambia": "GMB", "Ghana": "GHA", "Malawi": "MWI", "SierraLeone": "SLE"}
BATCH = int(os.environ.get("FE_BATCH", "40"))
CLIM_Y0, CLIM_Y1 = 1991, 2020
LST_Y0, LST_Y1 = 2003, 2020

ee.Initialize(project=PROJECT)


# ── helpers ──────────────────────────────────────────────────────────────────
def fieldwork_month(country):
    for r in survey_years_table():
        if r["country"] == country and r.get("respondent_median_date"):
            return int(r["respondent_median_date"][5:7])
    return 6


def pop_weight(country, year):
    """WorldPop count at the survey year, aggregated to a 250 m grid by summing the 100 m cells."""
    iso = ISO[country]
    col = ee.ImageCollection("WorldPop/GP/100m/pop").filter(ee.Filter.eq("country", iso))
    years = sorted(set(col.aggregate_array("year").getInfo()))
    pick = min(years, key=lambda y: abs(y - year))
    pop = ee.Image(col.filter(ee.Filter.eq("year", pick)).first()).select("population")
    w = (pop.unmask(0).reduceResolution(reducer=ee.Reducer.sum().unweighted(), maxPixels=1024)
         .reproject(crs="EPSG:4326", scale=250)).rename("w")
    return w, pick


def getinfo_retry(obj, what):
    for attempt in range(5):
        try:
            return obj.getInfo()
        except Exception as e:
            print(f"   retry {attempt + 1} ({what}): {str(e)[:140]}", flush=True)
            time.sleep(20 * (attempt + 1))
    raise RuntimeError("Earth Engine call failed five times: " + what)


def zonal(img, feats, scale, w=None, tile=4):
    """Area-weighted mean of every band, and (if w is given) the population-weighted mean.

    The weighted mean is sum(v * w) / sum(w over v's mask) per band, so a band's own
    missing pixels do not enter its denominator."""
    bands = img.bandNames().getInfo()
    fc = ee.FeatureCollection(feats)
    aw = img.reduceRegions(collection=fc, reducer=ee.Reducer.mean(), scale=scale, tileScale=tile)
    out_aw = getinfo_retry(aw, "area mean")["features"]
    rows = []
    for ft in out_aw:
        p = ft["properties"]
        rows.append({"country": p["country"], "Admin1": p["Admin1"], "Admin2": p["Admin2"],
                     **{b + "_aw": p.get(b) for b in bands}})
    if w is not None:
        num = img.multiply(w).rename([b + "_num" for b in bands])
        den = ee.Image.cat([w.updateMask(img.select(b).mask()).rename(b + "_den") for b in bands])
        st = num.addBands(den)
        pw = st.reduceRegions(collection=fc, reducer=ee.Reducer.sum(), scale=250, tileScale=tile)
        out_pw = getinfo_retry(pw, "population-weighted sums")["features"]
        for r, ft in zip(rows, out_pw):
            p = ft["properties"]
            for b in bands:
                n, d = p.get(b + "_num"), p.get(b + "_den")
                r[b + "_pw"] = (n / d) if (n is not None and d not in (None, 0)) else None
    return rows


def run_block(name, image_fn, scale, feats_by_country):
    """image_fn(country) -> (image, needs_pop_weight). Writes one CSV over all countries."""
    out_path = OUTDIR + f"gee_{name}_admin2.csv"
    rows = []
    t0 = time.time()
    for country, feats in feats_by_country.items():
        img = image_fn(country)
        w, pick = pop_weight(country, SURVEY_YEAR[country])
        print(f"  {country}: {len(feats)} polygons, WorldPop {pick}, {img.bandNames().size().getInfo()} bands", flush=True)
        for i in range(0, len(feats), BATCH):
            chunk = feats[i:i + BATCH]
            efeats = [ee.Feature(ee.Geometry(ft["geometry"]), {k: ft["properties"][k] for k in ("country", "Admin1", "Admin2")}) for ft in chunk]
            rows.extend(zonal(img, efeats, scale, w))
            print(f"    {min(i + BATCH, len(feats))}/{len(feats)} done, {time.time() - t0:.0f} s", flush=True)
    cols = ["country", "Admin1", "Admin2"] + sorted({k for r in rows for k in r if k not in ("country", "Admin1", "Admin2")})
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols); wr.writeheader()
        for r in rows:
            wr.writerow({k: ("" if r.get(k) is None else (round(r[k], 6) if isinstance(r[k], float) else r[k])) for k in cols})
    print(f"wrote {out_path}: {len(rows)} rows x {len(cols) - 3} columns", flush=True)


# ── block 1: climatology ─────────────────────────────────────────────────────
TC_SCALE = {"pr": 1.0, "tmmx": 0.1, "tmmn": 0.1, "pet": 0.1, "def": 0.1, "aet": 0.1, "soil": 0.1, "vpd": 0.01, "srad": 0.1}


def monthly_climatology(col, band, factor, y0, y1, prefix):
    """12 bands: the mean of `band` over the years y0..y1 for each calendar month."""
    imgs = []
    for m in range(1, 13):
        sub = col.filter(ee.Filter.calendarRange(y0, y1, "year")).filter(ee.Filter.calendarRange(m, m, "month")).select(band)
        imgs.append(sub.mean().multiply(factor).rename(f"{prefix}_m{m:02d}"))
    return ee.Image.cat(imgs)


def climatology_image(country):
    tc = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    pr = monthly_climatology(tc, "pr", 1.0, CLIM_Y0, CLIM_Y1, "pr")
    tmax = monthly_climatology(tc, "tmmx", 0.1, CLIM_Y0, CLIM_Y1, "tmax")
    ann = []
    for b, nm in (("tmmn", "tmin_ann"), ("pet", "pet_ann"), ("def", "def_ann"), ("aet", "aet_ann"), ("soil", "soilm_ann"), ("vpd", "vpd_ann"), ("srad", "srad_ann")):
        sub = tc.filter(ee.Filter.calendarRange(CLIM_Y0, CLIM_Y1, "year")).select(b)
        ann.append(sub.mean().multiply(TC_SCALE[b]).rename(nm))
    # annual totals of precipitation, one image per year, then their mean and s.d.
    yearly = ee.ImageCollection([tc.filter(ee.Filter.calendarRange(y, y, "year")).select("pr").sum().rename("pr_year") for y in range(CLIM_Y0, CLIM_Y1 + 1)])
    pr_ann_mean = yearly.mean().rename("pr_ann_mean"); pr_ann_sd = yearly.reduce(ee.Reducer.stdDev()).rename("pr_ann_sd")
    sy = SURVEY_YEAR[country]
    pr_sy = tc.filter(ee.Filter.calendarRange(sy, sy, "year")).select("pr").sum().rename("pr_sy")
    tmax_sy = tc.filter(ee.Filter.calendarRange(sy, sy, "year")).select("tmmx").mean().multiply(0.1).rename("tmax_sy")
    fm = fieldwork_month(country)
    end = ee.Date.fromYMD(sy, fm, 1).advance(1, "month"); start = end.advance(-12, "month")
    pr_win = tc.filterDate(start, end).select("pr").sum().rename("pr_win")
    lst = ee.ImageCollection("MODIS/061/MOD11A2")
    lstd = monthly_climatology(lst, "LST_Day_1km", 0.02, LST_Y0, LST_Y1, "lstd").subtract(273.15)
    lstn = monthly_climatology(lst, "LST_Night_1km", 0.02, LST_Y0, LST_Y1, "lstn").subtract(273.15)
    img = ee.Image.cat([pr, tmax] + ann + [pr_ann_mean, pr_ann_sd, pr_sy, tmax_sy, pr_win, lstd, lstn]).toFloat()
    return img.resample("bilinear")


# ── block 2: terrain ─────────────────────────────────────────────────────────
def terrain_image(country):
    mh = ee.Image("MERIT/Hydro/v1_0_1")
    parts = [mh.select("elv").rename("elv"), mh.select("hnd").rename("hnd"), mh.select("upa").rename("upa")]
    for nm in ("slope", "tri", "tpi", "roughness", "vrm", "cti", "elev-stdev"):
        ic = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/" + nm)
        parts.append(ic.mosaic().rename(nm.replace("-", "_")))
    return ee.Image.cat(parts).toFloat()


# ── block 3: iSDA soils ──────────────────────────────────────────────────────
ISDA = {  # short name: (dataset suffix, back-transform)
    "ph": ("ph", "div10"), "clay": ("clay_content", "pct"), "sand": ("sand_content", "pct"), "silt": ("silt_content", "pct"),   # texture is stored as plain %, the three sum to 100 (checked 2026-09-08)
    "bd": ("bulk_density", "div100"), "oc": ("carbon_organic", "explog"), "ntot": ("nitrogen_total", "explog100"),
    "cec": ("cation_exchange_capacity", "explog"), "zn": ("zinc_extractable", "explog"), "fe": ("iron_extractable", "explog"),
    "ca": ("calcium_extractable", "explog"), "mg": ("magnesium_extractable", "explog"), "k": ("potassium_extractable", "explog"),
    "p": ("phosphorus_extractable", "explog"), "s": ("sulphur_extractable", "explog"), "al": ("aluminium_extractable", "explog"),
}


def back(img, how):
    if how == "pct": return img
    if how == "div10": return img.divide(10)
    if how == "div100": return img.divide(100)
    if how == "explog": return img.divide(10).exp().subtract(1)
    if how == "explog100": return img.divide(100).exp().subtract(1)
    raise ValueError(how)


def isda_image(country):
    parts = []
    for short, (suffix, how) in ISDA.items():
        im = ee.Image("ISDASOIL/Africa/v1/" + suffix)
        for depth in ("0_20", "20_50"):
            parts.append(back(im.select("mean_" + depth), how).rename(f"{short}_{depth}"))
    return ee.Image.cat(parts).toFloat()


# ── main ─────────────────────────────────────────────────────────────────────
def main(blocks):
    gj = json.load(open(GEOJSON, encoding="utf-8"))
    feats_by_country = {c: [f for f in gj["features"] if f["properties"]["country"] == c] for c in ISO}
    print(f"{sum(len(v) for v in feats_by_country.values())} polygons; blocks: {', '.join(blocks)}", flush=True)
    if "clim" in blocks:
        print("[climatology]", flush=True); run_block("climatology", climatology_image, 1000, feats_by_country)
    if "terrain" in blocks:
        print("[terrain]", flush=True); run_block("terrain", terrain_image, 93, feats_by_country)
    if "isda" in blocks:
        print("[isda]", flush=True); run_block("isda", isda_image, 250, feats_by_country)
    print("DONE", flush=True)


if __name__ == "__main__":
    blocks = sys.argv[1].split(",") if len(sys.argv) > 1 else ["clim", "terrain", "isda"]
    main(blocks)
