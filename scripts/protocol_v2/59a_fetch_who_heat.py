"""
scripts/protocol_v2/59a_fetch_who_heat.py   [AB-01, 2026-09-15]

Download the WHO Health Inequality Data Repository (HEAT) datasets that carry
subnational-region estimates from household surveys, and keep the rows for the
four countries. No key, no login; the files are the bulk exports linked from
https://www.who.int/data/inequality-monitor/data (the datasafe Azure endpoint).

Consumed by scripts/protocol_v2/59_build_addback_sources.R (section G), which
takes the MICS rows: the only MICS source on disk for Sierra Leone (2010, 2017)
and Malawi (2014, district level). The DHS rows duplicate the project's own
DHS-derived Admin-2 block and are kept in the CSV for reference only.

    python scripts/protocol_v2/59a_fetch_who_heat.py            # skip files already present
    python scripts/protocol_v2/59a_fetch_who_heat.py --force    # re-download
-> data/WHO_HEAT/<dataset>.xlsx
-> data/WHO_HEAT/heat_subnational_4countries.csv
"""
import os, sys, urllib.request
import pandas as pd

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
OUT = os.path.join(ROOT, "data", "WHO_HEAT")
DATASETS = ["rep_dhs_unicef_rmnch", "rep_imm", "rep_dhs_unicef_malaria"]
URL = "https://datasafe-h5afbhf4gwctabaa.z01.azurefd.net/api/Download/TOP/{ds}/data"
ISO3 = ["GMB", "GHA", "MWI", "SLE"]
force = "--force" in sys.argv

os.makedirs(OUT, exist_ok=True)
frames = []
for ds in DATASETS:
    f = os.path.join(OUT, f"{ds}.xlsx")
    if force or not os.path.exists(f):
        print(f"[heat] downloading {ds} ...", flush=True)
        urllib.request.urlretrieve(URL.format(ds=ds), f)
    d = pd.read_excel(f, sheet_name="Data")
    d = d[d.iso3.isin(ISO3) & (d.dimension == "Subnational region")]
    print(f"[heat] {ds}: {len(d)} subnational rows for the four countries", flush=True)
    frames.append(d)
D = pd.concat(frames, ignore_index=True)
D.to_csv(os.path.join(OUT, "heat_subnational_4countries.csv"), index=False)
print(D.groupby(["iso3", "source", "date"]).agg(n_ind=("indicator_abbr", "nunique"),
                                              n_reg=("subgroup", "nunique")).to_string())
print(f"-> {os.path.join(OUT, 'heat_subnational_4countries.csv')} ({len(D)} rows)")
