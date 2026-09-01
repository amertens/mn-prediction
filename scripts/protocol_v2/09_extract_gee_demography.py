# =============================================================================
# scripts/protocol_v2/09_extract_gee_demography.py
#
# Earth Engine extraction of DEMOGRAPHIC COMPOSITION and URBANISATION at
# Admin-2, for the four analysis countries.
#
# WHY THESE TWO, AND WHY THEY WERE MISSING
# ----------------------------------------
# The vocabulary has 428 predictors and not one describes WHO LIVES THERE. Age
# and sex composition is the most direct area-level determinant of nutritional
# REQUIREMENT: a district with a high share of under-fives and women of
# reproductive age has a higher per-capita requirement for iron, vitamin A and
# folate than one with the same food supply and an older population. WorldPop
# age-sex has been used in this project only for post-stratification weights,
# never as a predictor. Degree of urbanisation (GHS-SMOD) is the standard
# settlement classifier and is more interpretable than the night-lights and
# built-fraction proxies currently standing in for it.
#
# WHY THIS IS A PYTHON SCRIPT AND NOT rgee
# ----------------------------------------
# rgee's ee_Initialize() rejects the stored credential with an "expired" check
# that the Python API itself does not agree with: ee.Initialize() authenticates
# and serves requests from the same credential file. Rather than re-authorise
# interactively (which needs a browser), this talks to Earth Engine directly.
# Note that ee.Initialize(project=...) fails for this account - the default
# project is the one that works.
#
# METHOD
#   - polygons: GADM Admin-2 for the four countries, simplified to 0.005 deg,
#     exported by the caller to GeoJSON (554 features).
#   - reduceRegions server-side, sum over each polygon, at 1 km. WorldPop is
#     native 100 m; summing counts at 1 km is exact enough for SHARES, which is
#     all that is derived here, and keeps the request small.
#   - batched, because a single request with 554 detailed polygons exceeds the
#     payload limit. Batch failures are reported, never silently skipped.
#
# DERIVED COLUMNS (shares, so they are comparable across countries)
#   wpop_share_under5        M_0+M_1+F_0+F_1 over total
#   wpop_share_wra           F_15..F_45 over total   (women of reproductive age)
#   wpop_share_over60        M_60..M_80 + F_60..F_80 over total
#   wpop_dependency_ratio    (under-15 + over-60) / working age
#   wpop_sex_ratio_wra       men 15-49 over women 15-49
#   wpop_log_density         log1p(total population per polygon)  [scale-bearing]
#   ghsl_smod_mean           mean GHS-SMOD class, 10 rural .. 30 urban centre
#
#   python scripts/protocol_v2/09_extract_gee_demography.py <in.geojson> <out.csv>
# =============================================================================
import json, sys, time
import ee

SURVEY_YEAR = {"Gambia": 2020, "Ghana": 2017, "Malawi": 2015, "SierraLeone": 2013}
ISO = {"Gambia": "GMB", "Ghana": "GHA", "Malawi": "MWI", "SierraLeone": "SLE"}
BATCH = 20
SCALE = 1000


def band_sum(img, names):
    """Sum a list of bands into one image, tolerating absent bands."""
    have = img.bandNames()
    keep = [n for n in names]
    return img.select(keep).reduce(ee.Reducer.sum())


def main(geojson_path, out_csv):
    ee.Initialize()
    gj = json.load(open(geojson_path, encoding="utf-8"))
    feats = gj["features"]
    print(f"{len(feats)} polygons", flush=True)

    # WorldPop age-sex bands, 5-year groups
    M = ["M_0", "M_1"] + [f"M_{a}" for a in range(5, 85, 5)]
    F = ["F_0", "F_1"] + [f"F_{a}" for a in range(5, 85, 5)]
    under5 = ["M_0", "M_1", "F_0", "F_1"]
    wra = [f"F_{a}" for a in (15, 20, 25, 30, 35, 40, 45)]
    mra = [f"M_{a}" for a in (15, 20, 25, 30, 35, 40, 45)]
    over60 = [f"M_{a}" for a in range(60, 85, 5)] + [f"F_{a}" for a in range(60, 85, 5)]
    under15 = under5 + ["M_5", "M_10", "F_5", "F_10"]
    working = [f"M_{a}" for a in range(15, 60, 5)] + [f"F_{a}" for a in range(15, 60, 5)]

    smod = ee.ImageCollection("JRC/GHSL/P2023A/GHS_SMOD")

    rows = {}
    for country, iso in ISO.items():
        sub = [f for f in feats if f["properties"]["country"] == country]
        if not sub:
            continue
        yr = SURVEY_YEAR[country]
        wp = ee.ImageCollection("WorldPop/GP/100m/pop_age_sex").filter(
            ee.Filter.eq("country", iso))
        n = wp.size().getInfo()
        if n == 0:
            print(f"  {country}: no WorldPop image", flush=True)
            continue
        # nearest available year
        years = sorted(wp.aggregate_array("year").getInfo())
        pick = min(years, key=lambda y: abs(y - yr))
        img = ee.Image(wp.filter(ee.Filter.eq("year", pick)).first())

        # GHS-SMOD epoch nearest the survey
        sm = ee.Image(smod.filter(
            ee.Filter.calendarRange(1990 + 5 * round((pick - 1990) / 5),
                                    1990 + 5 * round((pick - 1990) / 5),
                                    "year")).first())
        try:
            sm.bandNames().getInfo()
        except Exception:
            sm = ee.Image(smod.sort("system:time_start", False).first())

        stack = (img.select(["population"]).rename("pop_total")
                 .addBands(band_sum(img, under5).rename("pop_under5"))
                 .addBands(band_sum(img, wra).rename("pop_wra"))
                 .addBands(band_sum(img, mra).rename("pop_mra"))
                 .addBands(band_sum(img, over60).rename("pop_over60"))
                 .addBands(band_sum(img, under15).rename("pop_under15"))
                 .addBands(band_sum(img, working).rename("pop_working")))

        print(f"  {country}: WorldPop {pick}, {len(sub)} polygons", flush=True)
        for i in range(0, len(sub), BATCH):
            chunk = sub[i:i + BATCH]
            fc = ee.FeatureCollection([
                ee.Feature(ee.Geometry(f["geometry"]),
                           {"k": f"{f['properties']['country']}||"
                                 f"{f['properties']['Admin1']}||"
                                 f"{f['properties']['Admin2']}"})
                for f in chunk])
            for attempt in range(3):
                try:
                    a = stack.reduceRegions(fc, ee.Reducer.sum(), SCALE).getInfo()
                    b = sm.reduceRegions(fc, ee.Reducer.mean(), SCALE).getInfo()
                    bm = {f["properties"]["k"]: f["properties"].get("mean")
                          for f in b["features"]}
                    for f in a["features"]:
                        p = f["properties"]
                        rows[p["k"]] = {**p, "smod": bm.get(p["k"])}
                    break
                except Exception as e:
                    if attempt == 2:
                        print(f"    batch {i}-{i+len(chunk)} FAILED: "
                              f"{str(e)[:110]}", flush=True)
                    else:
                        time.sleep(3)
        print(f"    {len([k for k in rows if k.startswith(country)])} done",
              flush=True)

    # write derived shares
    import csv
    cols = ["country", "Admin1", "Admin2", "wpop_share_under5", "wpop_share_wra",
            "wpop_share_over60", "wpop_dependency_ratio", "wpop_sex_ratio_wra",
            "wpop_log_density", "ghsl_smod_mean"]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        import math
        for k, p in sorted(rows.items()):
            c, a1, a2 = k.split("||")
            tot = p.get("pop_total") or 0
            def sh(x):
                v = p.get(x)
                return round(v / tot, 6) if (v is not None and tot and tot > 0) else ""
            work = p.get("pop_working") or 0
            dep = ((p.get("pop_under15") or 0) + (p.get("pop_over60") or 0)) / work \
                if work else ""
            wrav, mrav = p.get("pop_wra") or 0, p.get("pop_mra") or 0
            sr = round(mrav / wrav, 6) if wrav else ""
            w.writerow([c, a1, a2, sh("pop_under5"), sh("pop_wra"),
                        sh("pop_over60"),
                        round(dep, 6) if dep != "" else "",
                        sr,
                        round(math.log1p(tot), 6) if tot else "",
                        round(p["smod"], 4) if p.get("smod") is not None else ""])
    print(f"wrote {out_csv}: {len(rows)} rows", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
