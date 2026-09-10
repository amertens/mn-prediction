"""
scripts/pptx_forum_postprocess.py  <deck.pptx> [more.pptx ...]

Run after pandoc has rendered a deck on the Micronutrient Forum template
(MNF2026-reference.pptx). Two adjustments pandoc cannot make itself:

1. Figure slides. Pandoc fits a picture inside the body placeholder and never
   lets it grow past 4.75 in high, so on the 13.33 x 7.5 in Forum slide a wide
   figure fills about two-thirds of the width. On every "Title and Content"
   slide that holds one picture and no body text, the picture is scaled to fill
   the area below the title (12.2 x 5.25 in), keeping its aspect ratio.

2. The label the Forum guidelines ask for on slides with unpublished data
   ("UNPUBLISHED DATA - DO NOT COPY OR DISTRIBUTE"). A small grey line is added
   at the foot of every slide except the title slide and section headers.

Edits in place.
"""
import sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE_TYPE

BOX_X, BOX_Y, BOX_W, BOX_H = Inches(0.55), Inches(1.75), Inches(12.23), Inches(5.1)
# two-column slides: picture left, text right
COL_PIC = (Inches(0.45), Inches(1.8), Inches(7.4))          # x, y, width
COL_TXT = (Inches(8.05), Inches(1.8), Inches(4.85), Inches(5.0))  # x, y, width, height
LABEL = "Unpublished data. Do not copy or distribute."
SKIP_LAYOUTS = {"Title Slide", "Section Header"}


def body_text_shapes(slide):
    out = []
    for sh in slide.shapes:
        if not sh.has_text_frame or not sh.text_frame.text.strip():
            continue
        if sh.is_placeholder and sh.placeholder_format.type is not None:
            t = str(sh.placeholder_format.type)
            if "TITLE" in t or "SLIDE_NUMBER" in t or "FOOTER" in t or "DATE" in t:
                continue
        out.append(sh)
    return out


def enlarge_pictures(slide):
    pics = [sh for sh in slide.shapes if sh.shape_type == MSO_SHAPE_TYPE.PICTURE]
    texts = body_text_shapes(slide)
    layout = slide.slide_layout.name
    if layout == "Title and Content" and len(pics) == 1 and not texts:
        pic = pics[0]
        scale = min(BOX_W / pic.width, BOX_H / pic.height)
        w, h = int(pic.width * scale), int(pic.height * scale)
        pic.left, pic.top, pic.width, pic.height = int(BOX_X + (BOX_W - w) / 2), int(BOX_Y), w, h
        return 1
    if layout == "Two Content" and len(pics) == 1 and len(texts) == 1:
        # pandoc puts the picture in the left column at the column width; widen it
        # and move the text column right
        pic, txt = pics[0], texts[0]
        x, y, w = COL_PIC
        h = int(pic.height * (w / pic.width))
        pic.left, pic.top, pic.width, pic.height = int(x), int(y), int(w), h
        txt.left, txt.top, txt.width, txt.height = [int(v) for v in COL_TXT]
        return 1
    return 0


def add_label(slide):
    if slide.slide_layout.name in SKIP_LAYOUTS:
        return 0
    tb = slide.shapes.add_textbox(Inches(0.92), Inches(6.98), Inches(7.5), Inches(0.35))
    tf = tb.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    r = p.add_run()
    r.text = LABEL
    r.font.size = Pt(11)
    r.font.color.rgb = RGBColor(0x80, 0x80, 0x80)
    return 1


def main(paths):
    for path in paths:
        prs = Presentation(path)
        n_pic = sum(enlarge_pictures(s) for s in prs.slides)
        n_lab = sum(add_label(s) for s in prs.slides)
        prs.save(path)
        print("%s: enlarged %d figure(s), labelled %d slide(s)" % (path, n_pic, n_lab))


if __name__ == "__main__":
    main(sys.argv[1:])
