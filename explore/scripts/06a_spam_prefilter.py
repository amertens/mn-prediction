"""
explore/scripts/06a_spam_prefilter.py   [probe MX-01, step 1]

Stream the 286 MB global MapSPAM production grid and keep only the four study
countries and the 42 crop columns, so the R step never has to hold the global
grid in memory. Flat memory: one row at a time.

  python explore/scripts/06a_spam_prefilter.py
-> explore/out/06a_spam_prod_4countries.csv
"""
import csv
import os
import sys

ROOT = r"C:\Users\andre\OneDrive\Documents\mn-prediction"
SRC = os.path.join(ROOT, "data", "MapSPAM", "raw", "spam2010V2r0_global_P_TA.csv")
DST = os.path.join(ROOT, "explore", "out", "06a_spam_prod_4countries.csv")
KEEP_ISO = {"GMB", "GHA", "MWI", "SLE"}

CROPS = ["whea", "rice", "maiz", "barl", "pmil", "smil", "sorg", "ocer",
         "pota", "swpo", "yams", "cass", "orts",
         "bean", "chic", "cowp", "pige", "lent", "opul", "soyb",
         "grou", "cnut", "oilp", "sunf", "rape", "sesa", "ooil",
         "sugc", "sugb", "cott", "ofib", "acof", "rcof", "coco", "teas", "toba",
         "bana", "plnt", "trof", "temf", "vege", "rest"]

csv.field_size_limit(1 << 24)

with open(SRC, newline="", encoding="utf-8", errors="replace") as fh:
    rdr = csv.reader(fh)
    hdr = next(rdr)
    lc = [h.strip().lower() for h in hdr]

    def find(name):
        try:
            return lc.index(name)
        except ValueError:
            sys.exit(f"column {name!r} not found; header starts {hdr[:8]}")

    i_iso, i_x, i_y = find("iso3"), find("x"), find("y")
    # crop columns are "<code>_<tech>", e.g. whea_a
    bare = [c[:-2] if len(c) > 2 and c[-2] == "_" else c for c in lc]
    crop_idx, crop_name = [], []
    for c in CROPS:
        if c in bare:
            crop_idx.append(bare.index(c))
            crop_name.append(c)
    missing = [c for c in CROPS if c not in crop_name]
    print(f"crop columns found: {len(crop_name)} of {len(CROPS)}"
          + (f"; missing {missing}" if missing else ""))

    os.makedirs(os.path.dirname(DST), exist_ok=True)
    n_in = n_out = 0
    with open(DST, "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out)
        w.writerow(["iso3", "x", "y"] + crop_name)
        for row in rdr:
            n_in += 1
            iso = row[i_iso].strip().upper()
            if iso not in KEEP_ISO:
                continue
            n_out += 1
            w.writerow([iso, row[i_x], row[i_y]] + [row[j] for j in crop_idx])

print(f"scanned {n_in:,} rows, kept {n_out:,} for {sorted(KEEP_ISO)}")
print("->", DST)
