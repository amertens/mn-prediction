"""
explore/scripts/26_mimi_malawi_block.py   [HC-06, task 1]

Malawi district-level micronutrient inadequacy from the MIMI/WFP paper's
supplementary tables into an Admin-2 predictor block.

Source: Tang et al., "The risk of dietary multiple micronutrient inadequacies is
widespread and geographically varied in Malawi", BMC Nutrition 2026,
doi 10.1186/s40795-026-01369-2, CC BY 4.0. IHS5 (April 2019 - April 2020).
  Additional file 3 (Table 2): district prevalence, vitamins
  Additional file 4 (Table 3): district prevalence, minerals Ca Fe Se Zn

The project's Malawi spine is Admin1 = district, Admin2 = Traditional Authority,
so the paper's districts join on Admin1 and broadcast to the TAs beneath. That
is a coarser rung than the analysis unit, exactly as the existing Ghana mimi_
columns are - recorded in the metadata so it is never mistaken for a TA-level
measurement.

  python explore/scripts/26_mimi_malawi_block.py
-> explore/out/26_mimi_malawi_admin2.csv (+ _metadata.csv, _unmatched.csv)
"""
import csv
import os
import re
import zipfile
from xml.etree import ElementTree as ET

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
SRC = ROOT + "explore/data/mimi_malawi/"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def doc_table(path):
    x = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    for t in x.iter(W + "tbl"):
        rows = []
        for tr in t.iter(W + "tr"):
            rows.append([
                "".join(n.text or "" for n in tc.iter(W + "t")).strip()
                for tc in tr.findall(W + "tc")])
        return rows
    return []


def pct(cell):
    """'50.8 (42.5-59.1)' -> (50.8, 42.5, 59.1); handles the en-dash mojibake."""
    if not cell:
        return (None, None, None)
    c = cell.replace("\u2013", "-").replace("\ufffd", "-").replace("�", "-")
    m = re.match(r"\s*([\d.]+)\s*\(\s*([\d.]+)\s*-\s*([\d.]+)\s*\)", c)
    if m:
        return tuple(float(g) for g in m.groups())
    m = re.match(r"\s*([\d.]+)\s*$", c)
    return (float(m.group(1)), None, None) if m else (None, None, None)


REGION_ROWS = {"northern region", "central region", "southern region"}
# paper column label -> our short nutrient tag
KEEP = {"vitamin a (rae)": "vita", "vitamin b9": "folate", "vitamin b12": "b12",
        "fe": "iron", "zn": "zinc", "se": "selenium", "ca": "calcium",
        "vitamin c": "vitc", "vitamin e": "vite", "vitamin b2": "b2",
        "vitamin b3": "b3", "vitamin b6": "b6"}

rows_out = {}
for fn in ("40795_2026_1369_MOESM3_ESM.docx", "40795_2026_1369_MOESM4_ESM.docx"):
    tab = doc_table(SRC + fn)
    if not tab:
        print("no table in", fn)
        continue
    hdr = [h.strip() for h in tab[0]]
    print(f"{fn}: {len(tab)} rows | columns: {hdr}")
    for r in tab[1:]:
        if not r or not r[0]:
            continue
        name = r[0].strip()
        if name.lower() in REGION_ROWS or name.lower().startswith("district"):
            continue
        rec = rows_out.setdefault(name, {"district": name})
        if len(r) > 1 and r[1].strip().isdigit():
            rec["n_households"] = int(r[1])
        for j, h in enumerate(hdr):
            tag = KEEP.get(h.strip().lower())
            if tag is None or j >= len(r):
                continue
            v, lo, hi = pct(r[j])
            if v is not None:
                rec[f"mimiMW_{tag}_inadequate_pct"] = v
                rec[f"mimiMW_{tag}_lo"] = lo
                rec[f"mimiMW_{tag}_hi"] = hi

D = list(rows_out.values())
print(f"\nparsed {len(D)} districts")

# ── join to the project's Malawi spine on Admin1 ────────────────────────────
spine = [r for r in csv.DictReader(open(ROOT + "metadata/admin2_spine.csv",
                                        encoding="utf-8-sig"))
         if r["country"] == "Malawi"]
print(f"spine Malawi rows: {len(spine)} | distinct Admin1: "
      f"{len({r['Admin1'] for r in spine})}")


def key(s):
    return re.sub(r"[^a-z]", "", (s or "").lower())


# The paper splits Blantyre, Lilongwe and Zomba into city and non-city, and
# reports Mzuzu city separately (it sits inside Mzimba); the spine has one unit
# for each. Merge the parts by household-weighted mean, which is what a district
# figure would have been. "Tchisi" is the paper's spelling of Ntchisi.
ALIAS = {"tchisi": "ntchisi"}
MERGE = {"blantyre": ["blantyre city", "blantyre non-city"],
         "lilongwe": ["lilongwe city", "lilongwe non-city"],
         "zomba":    ["zomba city", "zomba non-city"],
         "mzimba":   ["mzimba", "mzuzu city"]}

byk = {}
for d in D:
    byk[ALIAS.get(key(d["district"]), key(d["district"]))] = d

for tgt, parts in MERGE.items():
    got = [byk.get(key(p)) or next((x for x in D if key(x["district"]) == key(p)), None)
           for p in parts]
    got = [g for g in got if g]
    if len(got) < 2:
        continue
    wts = [float(g.get("n_households") or 0) for g in got]
    if sum(wts) <= 0:
        wts = [1.0] * len(got)
    merged = {"district": tgt, "n_households": sum(
        float(g.get("n_households") or 0) for g in got)}
    for c in {k for g in got for k in g if k.startswith("mimiMW_")}:
        num = [(float(g[c]), w) for g, w in zip(got, wts)
               if g.get(c) is not None]
        if num:
            merged[c] = sum(v * w for v, w in num) / sum(w for _, w in num)
    byk[tgt] = merged
    print(f"  merged {parts} -> {tgt} "
          f"(n_hh {merged['n_households']:.0f})")
val_cols = sorted({k for d in D for k in d
                   if k.startswith("mimiMW_")} | {"n_households"})

out, miss_spine, matched = [], set(), set()
for r in spine:
    d = byk.get(key(r["Admin1"]))
    row = {"country": "Malawi", "Admin1": r["Admin1"], "Admin2": r["Admin2"]}
    if d:
        matched.add(key(r["Admin1"]))
        for c in val_cols:
            row[c] = d.get(c)
    else:
        miss_spine.add(r["Admin1"])
        for c in val_cols:
            row[c] = None
    out.append(row)

unmatched_paper = [d["district"] for d in D if key(d["district"]) not in matched]
os.makedirs(ROOT + "explore/out", exist_ok=True)
with open(ROOT + "explore/out/26_mimi_malawi_admin2.csv", "w", newline="",
          encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["country", "Admin1", "Admin2"] + val_cols)
    w.writeheader()
    w.writerows(out)

cov = sum(1 for r in out if r.get("mimiMW_iron_inadequate_pct") is not None)
print(f"\nwrote {len(out)} spine rows | {cov} with values "
      f"({100*cov/len(out):.0f}%) | {len(matched)} of {len(D)} paper districts matched")
print("spine Admin1 with no paper district:", sorted(miss_spine)[:12])
print("paper districts not in the spine   :", sorted(unmatched_paper)[:12])

with open(ROOT + "explore/out/26_mimi_malawi_unmatched.csv", "w", newline="",
          encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["side", "name"])
    for n in sorted(miss_spine):
        w.writerow(["spine_admin1", n])
    for n in sorted(unmatched_paper):
        w.writerow(["paper_district", n])

meta = [{"column": c,
         "domain": "Dietary inadequacy (MODELLED SURFACE, HCES)",
         "source": "Tang et al. 2026 BMC Nutrition (MIMI/WFP), Malawi IHS5 "
                   "2019/20, Additional files 3-4, CC BY 4.0",
         "subnational": "TRUE",
         "resolution": "district (spine Admin1), broadcast to Traditional "
                       "Authorities beneath - NOT a TA-level measurement",
         "tier": "open"} for c in val_cols]
with open(ROOT + "explore/out/26_mimi_malawi_admin2_metadata.csv", "w",
          newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(meta[0]))
    w.writeheader()
    w.writerows(meta)
print(f"\n{len(val_cols)} columns: {', '.join(val_cols[:8])} ...")
