# =============================================================================
# src/Tanzania/0_extract_pr_crosswalk.py
#
# Extract the blood-sample -> cluster/weight/demographics crosswalk from the
# Tanzania 2010 PR (household-member) recode. Written in Python because R's
# haven AND foreign both SEGFAULT on TZPR63FL.DTA (a large DHS household recode);
# pandas.read_stata reads it cleanly. Run once before 1_GW_Tanzania_data_clean.R.
#
#   <venv>/python.exe src/Tanzania/0_extract_pr_crosswalk.py
#
# Output: data/IPD/Tanzania 2010/pr_crosswalk.csv  (one row per blood sample)
#   blood_id, cnum, svy_weight, strata, region_code, age, sex, child_flag
# Blood ids: women = SH615 (child_flag 0), children = SH516 (child_flag 1).
# =============================================================================
import os, sys
import pandas as pd

D = "data/IPD/Tanzania 2010"
PR = os.path.join(D, "TZPR63DT", "TZPR63FL.DTA")
OUT = os.path.join(D, "pr_crosswalk.csv")

cols = ["hv001", "hv005", "hv022", "hv024", "hv104", "hv105", "sh615", "sh516"]
pr = pd.read_stata(PR, columns=cols, convert_categoricals=False)

def build(idcol, flag):
    x = pr[[idcol, "hv001", "hv005", "hv022", "hv024", "hv104", "hv105"]].copy()
    x["blood_id"] = x[idcol].astype(str).str.strip().str.upper()
    x = x[(x["blood_id"] != "") & (x["blood_id"].str.lower() != "nan")]
    x["child_flag"] = flag
    return x.drop(columns=[idcol])

xw = pd.concat([build("sh615", 0), build("sh516", 1)], ignore_index=True)
xw = xw.rename(columns={"hv001": "cnum", "hv022": "strata",
                        "hv024": "region_code", "hv104": "sex", "hv105": "age"})
xw["svy_weight"] = xw["hv005"] / 1e6   # DHS weight scaling
xw = xw[["blood_id", "cnum", "svy_weight", "strata", "region_code", "age", "sex", "child_flag"]]

# Drop duplicate blood ids (keep first) so the downstream join is 1:1.
xw = xw.drop_duplicates(subset="blood_id", keep="first")

xw.to_csv(OUT, index=False)
print(f"Wrote {OUT}: {len(xw)} rows "
      f"({(xw.child_flag==0).sum()} women, {(xw.child_flag==1).sum()} children), "
      f"{xw.cnum.nunique()} clusters")
