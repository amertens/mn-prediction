# =============================================================================
# scripts/protocol_v2/11_extract_gee_rwi_density.py
#
# Three additions, and the answer to "is there anything better time-matched
# than WorldPop's single 2020 vintage?"
#
# 1. RELATIVE WEALTH INDEX (Meta / Data for Good, via the Earth Engine
#    community catalog as a POINT FeatureCollection at 2.4 km). The current SES
#    block is nine DHS aggregates from a prior survey round; RWI is a
#    high-resolution wealth surface validated against DHS wealth measurements,
#    and it is the highest-resolution SES signal publicly available for these
#    countries. Extracted as the mean and SD of member points per Admin-2, plus
#    the point count, so a district resting on three points is distinguishable
#    from one resting on three hundred.
#
# 2. YEAR-MATCHED POPULATION DENSITY. This IS the better-time-matched source.
#    WorldPop/GP/100m/pop carries every year from 2000 to 2020, so density is
#    taken at each country's own survey year (Ghana 2017, Malawi 2015, Sierra
#    Leone 2013, Gambia 2020) rather than at a single global vintage. It
#    supersedes the 2020-only wpop_log_density built in step 09.
#
# 3. GPW 2010 AGE-SEX, as the second bracket. Checked directly against the
#    catalogue: WorldPop's age-sex collection exists ONLY for 2020, and GPW
#    v4.11 basic demographic characteristics exist ONLY for 2010. There is no
#    annual age-sex product. Those two vintages BRACKET the 2013-2021 survey
#    window, so both are carried and the metadata says which is which; age and
#    sex composition moves slowly, so a bracket is an honest representation of
#    what is knowable rather than a false precision.
#
#   python scripts/protocol_v2/11_extract_gee_rwi_density.py <in.geojson> <out.csv>
# =============================================================================
import csv, json, math, sys, time
import ee

SURVEY_YEAR = {"Gambia": 2018, "Ghana": 2017, "Malawi": 2015, "SierraLeone": 2013}   # Gambia fieldwork Jan-Apr 2018 (was 2020)
ISO = {"Gambia": "GMB", "Ghana": "GHA", "Malawi": "MWI", "SierraLeone": "SLE"}
BATCH = 20
SCALE = 100        # WorldPop 100 m has a MEAN pyramid: summing at 1 km returned 1/100 of the count (found 2026-09-07); GPW shares are ratios and unaffected


def main(geojson_path, out_csv):
    ee.Initialize()
    feats = json.load(open(geojson_path, encoding="utf-8"))["features"]
    print(f"{len(feats)} polygons", flush=True)

    rwi = ee.FeatureCollection(
        "projects/sat-io/open-datasets/facebook/relative_wealth_index")
    gpw = ee.ImageCollection("CIESIN/GPWv411/GPW_Basic_Demographic_Characteristics")

    # GPW 2010 age-sex: band ids encode age range, sex and year
    def gpw_band(code):
        img = gpw.filter(ee.Filter.stringContains("system:index", code)).first()
        return ee.Image(img).rename(code)

    # both totals (bt) for under-5, under-15, and female 15-49
    g_u5 = gpw_band("a000_004bt_2010")
    g_u15 = gpw_band("a000_014bt_2010")
    g_f1549 = gpw_band("a015_049ft_2010")
    g_tot = gpw_band("atotpopbt_2010")

    rows = {}
    for country, iso in ISO.items():
        sub = [f for f in feats if f["properties"]["country"] == country]
        if not sub:
            continue
        yr = SURVEY_YEAR[country]
        pop = ee.ImageCollection("WorldPop/GP/100m/pop").filter(
            ee.Filter.And(ee.Filter.eq("country", iso), ee.Filter.eq("year", yr)))
        if pop.size().getInfo() == 0:
            years = sorted(set(ee.ImageCollection("WorldPop/GP/100m/pop")
                               .filter(ee.Filter.eq("country", iso))
                               .aggregate_array("year").getInfo()))
            yr = min(years, key=lambda y: abs(y - SURVEY_YEAR[country]))
            pop = ee.ImageCollection("WorldPop/GP/100m/pop").filter(
                ee.Filter.And(ee.Filter.eq("country", iso),
                              ee.Filter.eq("year", yr)))
        popimg = ee.Image(pop.first()).rename("pop_year")
        print(f"  {country}: WorldPop density {yr}, {len(sub)} polygons", flush=True)

        stack = (popimg
                 .addBands(ee.Image.pixelArea().rename("area_m2"))
                 .addBands(g_u5).addBands(g_u15)
                 .addBands(g_f1549).addBands(g_tot))

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
                    # RWI: mean/sd/count of the point features inside each polygon
                    def rwi_stats(feat):
                        pts = rwi.filterBounds(feat.geometry())
                        return feat.set({
                            "rwi_mean": pts.aggregate_mean("rwi"),
                            "rwi_sd": pts.aggregate_total_sd("rwi"),
                            "rwi_n": pts.size(),
                        })
                    b = fc.map(rwi_stats).getInfo()
                    bm = {f["properties"]["k"]: f["properties"]
                          for f in b["features"]}
                    for f in a["features"]:
                        p = f["properties"]
                        rows[p["k"]] = {**p, **bm.get(p["k"], {}), "_yr": yr}
                    break
                except Exception as e:
                    if attempt == 2:
                        print(f"    batch {i}-{i+len(chunk)} FAILED: "
                              f"{str(e)[:110]}", flush=True)
                    else:
                        time.sleep(4)
        print(f"    {len([k for k in rows if k.startswith(country)])} done",
              flush=True)

    cols = ["country", "Admin1", "Admin2", "rwi_mean", "rwi_sd", "rwi_n_points",
            "wpop_log_density_survey_year", "wpop_density_year",
            "gpw2010_share_under5", "gpw2010_share_under15",
            "gpw2010_share_female_1549"]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for k, p in sorted(rows.items()):
            c, a1, a2 = k.split("||")
            tot = p.get("atotpopbt_2010") or 0

            def sh(b):
                v = p.get(b)
                return round(v / tot, 6) if (v is not None and tot) else ""
            # population per km2: the polygon count divided by its area (AU-01 finding 3; was the raw count)
            dens = (p.get("pop_year") / (p.get("area_m2") / 1e6)) if (p.get("pop_year") is not None and p.get("area_m2")) else None
            w.writerow([
                c, a1, a2,
                round(p["rwi_mean"], 5) if p.get("rwi_mean") is not None else "",
                round(p["rwi_sd"], 5) if p.get("rwi_sd") is not None else "",
                p.get("rwi_n", ""),
                round(math.log1p(dens), 6) if dens else "",
                p.get("_yr", ""),
                sh("a000_004bt_2010"), sh("a000_014bt_2010"),
                sh("a015_049ft_2010"),
            ])
    print(f"wrote {out_csv}: {len(rows)} rows", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
