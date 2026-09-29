"""Copy the best slides of the full talk into the rendered v4 MNF15 deck, each one next
to the v4 slide it most resembles (so the two can be compared and combined), the rest
into the appendix, applying Andrew's comments of 27 September on the way.

    bash scripts/render_deck.sh docs/slides/MNF15-talk-2026-09-v4.qmd --forum --no-label --text 18,12 \
         --concept docs/slides/MNF15-talk-2026-09-v4.concept.yaml
    python scripts/policy_deck/32_merge_full_talk_slides.py

Source: docs/slides/MN-proxy-full-talk-best slides.pptx (Andrew's selection, with his
comments; read only). Both decks use the Micronutrient Forum reference template, so a
slide is copied shape by shape onto the same layout, with its pictures re-registered and
its speaker notes kept; comments are not copied. Every copied slide carries a small grey
"full talk" tag, bottom left, so it can be told apart from its v4 neighbour.

Edits made to copied slides (Andrew's comments):
  7, 8   titles become descriptions ("questions don't work for this and the next slide title")
  17     notes answer his question; v4's "How was it tested?" is the corrected map
  18     figure redrawn with a readable legend and points (results/figures/mnf15_v4/v4_ft_prev_error.png)
  21     each combination gets its share of district pairs and its worst-fifth reach
  22     Ghana women's iron replaced by Malawi women's B12, the best-ranked combination
  31     an explanation of what the rank intervals show
  41     an updated version is a main-body conclusion slide in the qmd; the original goes to the appendix
Slides 29 onward go to the appendix, as his comment on slide 29 asks.

V5 (27 Sep): Andrew's MNF15-talk-2026-09-v5.pptx deleted full-talk slides 4, 7, 8, 9 and the
v4 predictors grid, and edited two copied slides by hand (Goals: an icon removed; the
domains slide re-laid out as "575 proxy variables from 29 conceptual domains" with his
image). Those two are copied from V5 itself ("v5:N" entries in PLAN).
"""
import copy
import io
import os
import sys

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.opc.constants import RELATIONSHIP_TYPE as RT
from pptx.util import Inches, Pt

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
SRC = os.path.join(ROOT, "docs/slides/MN-proxy-full-talk-best slides.pptx")
DST = os.path.join(ROOT, "docs/slides/MNF15-talk-2026-09-v4.pptx")
V5 = os.path.join(ROOT, "docs/slides/MNF15-talk-2026-09-v5.pptx")   # Andrew's V5 (27 Sep, read only): slides he re-laid out by hand
FIG = os.path.join(ROOT, "results/figures/mnf15_v4")

# (full-talk slide number, v4 slide title it follows, or None to follow the previous copy)
# A title prefixed with "=" means "after the slide with exactly this title".
PLAN = [
    ("v5:3", "Introduction and study objectives"),                 # Goals, as edited in V5 (third icon removed)
    ("v5:7", "Which of these causes can public data measure?"),   # V5's re-laid-out domains slide, with Andrew's image
    (6, "Which surveys does the model predict to?"),
    (12, "How does the model work, and how was it tested?"),
    (15, "How the model is built and tested, step by step"),
    (16, "How was it tested?"), (17, None),
    (3, "Why not a more complex machine-learning model?"), (13, None), (14, None),
    (22, "What about subnational prevalence estimates?"),
    (20, "How well does it rank districts, nutrient by nutrient?"), (24, None),
    (25, "Accuracy by nutrient, against the best possible"),
    (18, "How do proxy models compare with the survey's own information?"), (19, None),
    (23, "Which deficiencies can it map?"),
    (21, "What does this mean for policy? The vitamin B12 model"),
    (27, "What have we learned, and what comes next?"),
    # appendix
    (42, "Is it just geography? Comparison with DHS-style geostatistics"), (28, None),
    (39, "What is the cheapest way to get district numbers?"),
    (38, "How can it help plan the next survey?"),
    (44, "How does it compare with other methods, including machine learning?"), (45, None),
    (35, "What can we learn about missing predictor domains?"), (36, None), (37, None),
    (32, "What does the model rely on, for each women's outcome?"), (43, None),
    (33, "=Malawi: why the weights are not causes"),
    (34, "Malawi: why the weights are not causes"),
    (5, "Which outcomes and cut-offs does the analysis use?"), (10, None),
    (31, "How firmly is each Côte d'Ivoire district placed?"),
    (26, "Does targeting reach more deficient people?"), (29, None), (30, None),
    (40, "Accuracy pooled over all nutrients, in district pairs"),
    (41, "Person-level models: a coin toss on proxies alone"),
]
NEW_TITLES = {7: "The land, climate and agriculture layers", 8: "The health, food and household layers",
              22: "Malawi: survey, held-out prediction, every area ranked, and uncertainty"}
EXTRA_NOTES = {
    17: ("\n\nANSWER TO ANDREW'S COMMENT (27 Sep): the map is wrong for in-fill. Most of Ghana's districts have no survey data "
         "and should be grey, not training. And in-fill folds are NOT grouped in space: five random parts of the surveyed "
         "districts, ten draws; the region hold-out is the spatially grouped test. v4's 'How was it tested?' slide, just before, "
         "is the corrected version (grey = no survey data, the benchmark's own first fold)."),
    20: ("\n\nUPDATED (27 Sep): the numbers on this slide are current (benchmarks_v2_cells.csv, 19 Sep). The per-nutrient "
         "version in district pairs with 95% intervals is v4's 'How well does it rank districts, nutrient by nutrient?', just before."),
    18: ("\n\nREDRAWN (27 Sep) for Andrew's comment: a legend in words, the national number as a bar, larger points "
         "(scripts/policy_deck/31_mnf15_v4_january_updates.R). Means over 18 cells: model 10.3 points, regional average 11.8, "
         "one national number 11.5."),
    21: ("\n\nADDED (27 Sep): share of district pairs in the survey's order and the share of deficient people in the model's worst "
         "fifth (v3_percell_pairs.csv, nce_targeting_metrics.csv). Note The Gambia's child vitamin A and women's iron: order is "
         "strong, but the worst fifth by rate does not hold more deficient people than the regional averages' pick."),
    22: ("\n\nREPLACED (27 Sep) for Andrew's comment: Ghana women's iron (ranking 0.40) swapped for Malawi women's B12, the "
         "best-ranked combination (73 of 100 pairs on the prevalence scale, each area hidden). Panel 3 shows the deployed order; "
         "panel 4 the calibrated chance of 20% or more, which is low almost everywhere because the calibrated prevalence is "
         "shrunk toward the national 11 per cent. Tables: viz/oof_women_b12.csv, viz/deploy_malawi_women_b12_cal.csv "
         "(scripts/policy_deck/10_viz_tables.R A,B2 with VZ_* overrides)."),
    31: ("\n\nEXPLANATION ADDED (27 Sep) for Andrew's comment: see the text box on the slide."),
    41: ("\n\nThe original full-talk slide. The updated version is the main-body slide 'What have we learned, and what comes next?'."),
}
SLIDE21_ADD = {  # after the paragraph starting with this text, add this paragraph (same style)
    "Malawi: women's B12": "75 of 100 district pairs in order; its worst fifth holds 40% of deficient women (regional averages 27%)",
    "The Gambia: women's vitamin A": "74 of 100 pairs in order; worst fifth holds 25% of deficient women (regional averages 15%)",
    "The Gambia: child vitamin A": "74 of 100 pairs in order; worst fifth holds 14% of deficient children (regional averages 16%)",
    "The Gambia: women's iron": "74 of 100 pairs in order; worst fifth holds 14% of deficient women (regional averages 18%)",
}
SLIDE31_TEXT = [
    ("What this shows", True),
    ("Left: we refit the model 100 times on resampled training districts and take each district's middle 90% of ranks. "
     "If that range were a true 90% interval, the survey's rank would fall inside it 90% of the time (dashed line). "
     "It does so only 15 to 70% of the time, 38% on average: the ranges show how stable the ranking is, not how right it is.", False),
    ("Right: to contain the survey's rank 90% of the time, an interval would have to span 55% of the list either side, "
     "almost as wide as a random ordering (68%). So we report the order and calibrated prevalence bands, not rank intervals.", False),
]


def title_of(slide):
    t = slide.shapes.title
    return t.text_frame.text.strip().replace("’", "'").replace("‘", "'") if t is not None else ""   # pandoc curls apostrophes


def copy_slide(src_slide, dst):
    layout = next(l for l in dst.slide_layouts if l.name == src_slide.slide_layout.name)
    new = dst.slides.add_slide(layout)
    for ph in list(new.placeholders):
        ph._element.getparent().remove(ph._element)
    rid = {}
    for r_id, rel in src_slide.part.rels.items():
        if rel.reltype == RT.IMAGE:
            _, rid[r_id] = new.part.get_or_add_image_part(io.BytesIO(rel.target_part.blob))
        elif rel.reltype == RT.HYPERLINK and rel.is_external:
            rid[r_id] = new.part.relate_to(rel.target_ref, RT.HYPERLINK, is_external=True)
    tree = new.shapes._spTree
    for shp in src_slide.shapes:
        el = copy.deepcopy(shp._element)
        for node in el.iter():
            for a in list(node.attrib):
                if a.split("}")[-1] in ("embed", "link", "id") and node.attrib[a] in rid:
                    node.attrib[a] = rid[node.attrib[a]]
        tree.append(el)
    if src_slide.has_notes_slide:
        new.notes_slide.notes_text_frame.text = src_slide.notes_slide.notes_text_frame.text
    return new


def move_after(dst, new_slide, after_idx):
    lst = dst.slides._sldIdLst
    el = next(e for e in lst if e.rId == dst.part.relate_to(new_slide.part, RT.SLIDE))
    lst.remove(el)
    lst.insert(after_idx + 1, el)


def index_of(dst, title):
    for i, s in enumerate(dst.slides):
        if title_of(s) == title:
            return i
    sys.exit(f"v4 slide not found: {title!r}")


def set_title(slide, text):
    tf = slide.shapes.title.text_frame
    p = tf.paragraphs[0]
    if p.runs:
        p.runs[0].text = text
        for r in p.runs[1:]:
            r.text = ""
    else:
        tf.text = text


def tag(slide):
    tb = slide.shapes.add_textbox(Inches(0.25), Inches(7.08), Inches(1.6), Inches(0.3))
    p = tb.text_frame.paragraphs[0]
    r = p.add_run(); r.text = "full talk"; r.font.size = Pt(10); r.font.color.rgb = RGBColor(0x99, 0x99, 0x99); r.font.italic = True


def pictures(slide):
    return [s for s in slide.shapes if s.shape_type == 13]


def replace_picture(slide, png):
    pic = max(pictures(slide), key=lambda p: p.width * p.height)
    left, top, w, h = pic.left, pic.top, pic.width, pic.height
    pic._element.getparent().remove(pic._element)
    from PIL import Image
    iw, ih = Image.open(png).size
    scale = min(w / iw, h / ih)
    nw, nh = int(iw * scale), int(ih * scale)
    slide.shapes.add_picture(png, left + (w - nw) // 2, top + (h - nh) // 2, nw, nh)


def add_after_paragraph(slide, lead, text):
    """Slide 21: each combination is a heading box with a metrics box just below it, in the same column.
    Append one line to that metrics box, styled like its first line."""
    A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
    boxes = [sh for sh in slide.shapes if sh.has_text_frame]
    head = next((sh for sh in boxes if sh.text_frame.text.strip().startswith(lead)), None)
    if head is None:
        return False
    below = [sh for sh in boxes if abs(sh.left - head.left) < Inches(0.25) and sh.top > head.top]
    if not below:
        return False
    m = min(below, key=lambda sh: sh.top)
    first = m.text_frame.paragraphs[0]._p
    new_p = copy.deepcopy(first)
    runs = new_p.findall(A + "r")
    for k, r in enumerate(runs):
        r.find(A + "t").text = text if k == 0 else ""
    m.text_frame.paragraphs[-1]._p.addnext(new_p)
    return True


def main():
    src = Presentation(SRC)
    dst = Presentation(DST)
    src_slides = list(src.slides)
    v5_slides = list(Presentation(V5).slides)
    last_idx = None
    for num, after in PLAN:
        if isinstance(num, str):   # "v5:N": copy slide N of V5 as Andrew edited it (already tagged, no further edits)
            new = copy_slide(v5_slides[int(num.split(":")[1]) - 1], dst)
            idx = index_of(dst, after)
            move_after(dst, new, idx)
            last_idx = idx + 1
            print(f"  {num:>6} -> after {title_of(dst.slides[idx])[:60]!r}")
            continue
        s = src_slides[num - 1]
        new = copy_slide(s, dst)
        if num in NEW_TITLES:
            set_title(new, NEW_TITLES[num])
        if num == 18:
            replace_picture(new, os.path.join(FIG, "v4_ft_prev_error.png"))
        if num == 22:
            replace_picture(new, os.path.join(FIG, "v4_ft_malawi_b12_four.png"))
        if num == 21:
            for lead, text in SLIDE21_ADD.items():
                if not add_after_paragraph(new, lead, text):
                    print("  slide 21: could not place", lead)
        if num == 31:
            pic = max(pictures(new), key=lambda p: p.width * p.height)
            ratio = pic.height / pic.width
            pic.width = Inches(8.6); pic.height = int(pic.width * ratio); pic.left = Inches(4.5)
            tb = new.shapes.add_textbox(Inches(0.45), Inches(1.75), Inches(3.9), Inches(5.0))
            tb.text_frame.word_wrap = True
            for k, (txt, bold) in enumerate(SLIDE31_TEXT):
                p = tb.text_frame.paragraphs[0] if k == 0 else tb.text_frame.add_paragraph()
                r = p.add_run(); r.text = txt; r.font.size = Pt(14 if bold else 12.5); r.font.bold = bold
                r.font.color.rgb = RGBColor(0x0F, 0x7B, 0x8A) if bold else RGBColor(0x33, 0x33, 0x33)
                p.space_after = Pt(6)
        if num in EXTRA_NOTES:
            nt = new.notes_slide.notes_text_frame
            nt.text = nt.text + EXTRA_NOTES[num]
        tag(new)
        if after is None:
            idx = last_idx
        elif after.startswith("="):
            idx = index_of(dst, after[1:]) - 1   # before that slide
        else:
            idx = index_of(dst, after)
        move_after(dst, new, idx)
        last_idx = idx + 1
        print(f"  full-talk {num:2d} -> after {title_of(dst.slides[idx])[:60]!r}")
    dst.save(DST)
    print(f"merged {len(PLAN)} full-talk slides into {DST}: {len(dst.slides)} slides")
    # readable speaker notes (Andrew, 27 Sep: "the speaker notes are tiny and I can't increase them")
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import pptx_notes_size
    pptx_notes_size.main(DST, 18)


if __name__ == "__main__":
    main()
