"""scripts/policy_deck/03c_civ_climatology_gee.py   [CV-01]

TerraClimate / MODIS LST climatology for Cote d'Ivoire's 33 Admin-2 polygons,
through the same functions as scripts/protocol_v2/54_extract_climatology_terrain_soil.py
(imported from that file), so the columns match gee_climatology_admin2.csv exactly.
CIV has no biomarker survey: the "survey year" the anomaly columns need is
CIV_REF_YEAR (default 2016). The script's ISO / SURVEY_YEAR tables gain a
CoteDIvoire entry for this run only, and its output directory is redirected.

    "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe" scripts/policy_deck/03c_civ_climatology_gee.py
-> data/external_cache/gee_geoms/civ/gee_climatology_admin2.csv
"""
import importlib.util
import json
import os
import sys

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
sys.path.insert(0, ROOT + "scripts/protocol_v2")
spec = importlib.util.spec_from_file_location("s54", ROOT + "scripts/protocol_v2/54_extract_climatology_terrain_soil.py")
s54 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s54)

REF_YEAR = int(os.environ.get("CIV_REF_YEAR", "2016"))
s54.ISO["CoteDIvoire"] = "CIV"
s54.SURVEY_YEAR["CoteDIvoire"] = REF_YEAR
OUT = ROOT + "data/external_cache/gee_geoms/civ/"
os.makedirs(OUT, exist_ok=True)
s54.OUTDIR = OUT

gj = json.load(open(ROOT + "data/external_cache/gee_geoms/civ_admin2_simplified.geojson", encoding="utf-8"))
feats = gj["features"]
for f in feats:
    f["properties"] = {k: f["properties"][k] for k in ("country", "Admin1", "Admin2")}
print(f"{len(feats)} CIV polygons, reference year {REF_YEAR}", flush=True)
s54.run_block("climatology", s54.climatology_image, 1000, {"CoteDIvoire": feats})
print("DONE", flush=True)
