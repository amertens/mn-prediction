"""
explore/scripts/25_glss7_food_stream.py   [HC-06, step 1]

Stream GLSS7 section 9b (3.2 GB, household x item x six visits) down to a
compact household x item table, so the R side can classify items with the
project's OWN classify_item() and hh_from_items() and stay pooled-comparable
with Malawi, The Gambia and Sierra Leone.

Section 9b is the frequently-purchased-items diary: for each visit N there is
s9bqNa (amount spent), s9bqNb (quantity acquired) and s9bqNc (unit). A household
is taken to have CONSUMED an item if any visit records a positive amount or
quantity, matching the `consumed` flag the other three countries use.

Also writes the freqcd item-code -> label map, which is what the classifier
needs. Item labels cover food and non-food (the diary includes soap and pens),
so the classifier's "other" bucket does the selecting, as it does elsewhere.

  python explore/scripts/25_glss7_food_stream.py
-> explore/out/25_glss7_hh_item.csv       hid, code, consumed, value, qty
   explore/out/25_glss7_item_labels.csv   code, item
"""
import os
import pandas as pd

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
SRC = ROOT + "data/LSMS/GHA_2017/g7sec9b.dta"
OUT_HH = ROOT + "explore/out/25_glss7_hh_item.csv"
OUT_LAB = ROOT + "explore/out/25_glss7_item_labels.csv"
CHUNK = 250_000

VAL = [f"s9bq{i}a" for i in range(1, 7)]   # amount spent at visit i
QTY = [f"s9bq{i}b" for i in range(1, 7)]   # quantity acquired at visit i

# ── item labels ─────────────────────────────────────────────────────────────
rdr0 = pd.io.stata.StataReader(SRC)
vls = rdr0.value_labels()
lab = None
for k, v in vls.items():
    if k.lower().startswith("freqcd") or "freq" in k.lower():
        lab = v
        break
if lab is None:                      # fall back to the largest label set
    lab = max(vls.values(), key=len)
pd.DataFrame({"code": list(lab.keys()), "item": list(lab.values())}).to_csv(
    OUT_LAB, index=False)
print(f"item labels: {len(lab)} -> {OUT_LAB}", flush=True)

# ── stream ──────────────────────────────────────────────────────────────────
keep = ["hid", "clust", "freqcd"] + VAL + QTY
parts, n_in, n_out = [], 0, 0
rdr = pd.read_stata(SRC, chunksize=CHUNK, convert_categoricals=False,
                    columns=keep)
for i, ch in enumerate(rdr):
    n_in += len(ch)
    v = ch[VAL].fillna(0).sum(axis=1)
    q = ch[QTY].fillna(0).sum(axis=1)
    m = (v > 0) | (q > 0)
    if m.any():
        out = pd.DataFrame({
            "hid": ch.loc[m, "hid"].astype(str),
            "clust": ch.loc[m, "clust"],
            "code": ch.loc[m, "freqcd"].astype("Int64"),
            "value": v[m].values, "qty": q[m].values})
        parts.append(out)
        n_out += len(out)
    if i % 20 == 0:
        print(f"  chunk {i:4d}  read {n_in:>10,}  kept {n_out:>9,}", flush=True)

D = pd.concat(parts, ignore_index=True)
D["consumed"] = True
D.to_csv(OUT_HH, index=False)
print(f"\nread {n_in:,} rows, kept {n_out:,} with positive value or quantity")
print(f"households {D.hid.nunique():,} | items {D.code.nunique()} | -> {OUT_HH}")
print(f"output size {os.path.getsize(OUT_HH)/1e6:.1f} MB")
