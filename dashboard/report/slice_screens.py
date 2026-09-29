"""Slice tall dashboard screenshots into page-sized JPEGs for the Word export.

    python dashboard/report/slice_screens.py <capture_dir>

Reads manifest.json in the capture folder and writes slices/<name>_<k>.jpg
plus slices.json ({png name: [[jpg, width, height], ...]}). A slice is at most
1.3 times as tall as it is wide, so it fits a portrait page at full width.
"""
import json
import os
import sys

from PIL import Image

cap = sys.argv[1]
man = json.load(open(os.path.join(cap, "manifest.json"), encoding="utf-8"))
os.makedirs(os.path.join(cap, "slices"), exist_ok=True)
out = {}
for page in man["pages"]:
    for png in page["shots"]:
        im = Image.open(os.path.join(cap, png)).convert("RGB")
        w, h = im.size
        max_h = int(w * 1.3)
        parts, top, k = [], 0, 0
        while top < h:
            bottom = min(h, top + max_h)
            if h - bottom < w * 0.25:  # don't leave a sliver on its own
                bottom = h
            piece = im.crop((0, top, w, bottom))
            name = f"{os.path.splitext(png)[0]}_{k}.jpg"
            piece.save(os.path.join(cap, "slices", name), quality=82, optimize=True)
            parts.append([name, w, bottom - top])
            top, k = bottom, k + 1
        out[png] = parts
json.dump(out, open(os.path.join(cap, "slices.json"), "w"), indent=1)
print(f"sliced {sum(len(v) for v in out.values())} pieces from {len(out)} screenshots")
