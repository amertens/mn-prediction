"""Sonja Hess's conceptual framework of iron deficiency, with a coverage badge on
each of HER boxes, for the v2 MNF15 talk (slide 3, a two-slide build).

    python scripts/policy_deck/27_mnf15_v2_framework.py

Base figure: docs/slides/img/old_s06_framework_1.png, the same framework as the
slide Sonja sent on 23 September (Hess et al., Ann N Y Acad Sci 2023, adapted)
without the text cursor that is baked into her screenshot after "Ecology".
Box positions were measured from its border lines (27 September). The coverage
call for each box follows the 19 September coverage slide where that slide had
the box; the boxes that slide dropped are marked DRAFT below and need Sonja's
sign-off (obesity, family planning, knowledge, built environment, other
micronutrient deficiencies, genetic risk).

Writes results/figures/mnf15_v2/v2_framework_plain.png (slide 3a) and
v2_framework_coverage.png (slide 3b). The legend and the list of training
surveys are editable text on the slide (docs/slides/MNF15-talk-2026-09-v2.concept.yaml).
"""
import os

from PIL import Image, ImageDraw

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
BASE = os.path.join(ROOT, "docs/slides/img/old_s06_framework_1.png")
OUT = os.path.join(ROOT, "results/figures/mnf15_v2")
os.makedirs(OUT, exist_ok=True)

GREEN, AMBER, GREY, TEAL = (46, 125, 50), (224, 161, 0), (150, 150, 150), (15, 123, 138)

# (x0, y0, x1, y1) in the 927 x 696 base image, and the coverage call
BOXES = {
    "Physiological vulnerability":           ((152, 8, 588, 77), AMBER),    # age-sex structure and fertility only
    "Genetic risk factors":                  ((599, 8, 814, 77), GREEN),    # Malaria Atlas sickle-cell / G6PD surfaces (DRAFT)
    "Ecology, climate, geography, ...":      ((152, 111, 814, 156), GREEN),
    "Poverty":                               ((153, 179, 278, 258), GREEN),
    "Low educational attainment":            ((287, 179, 412, 258), GREEN),
    "Cultural norms and behaviors":          ((420, 179, 545, 258), GREY),
    "Inadequate built environment":          ((554, 179, 679, 258), GREEN),  # built-up area, night lights, housing (DRAFT)
    "Health policies":                       ((690, 179, 815, 258), AMBER),  # national fortification / supplementation only
    "Food insecurity":                       ((111, 295, 230, 398), GREEN),
    "Inadequate maternal and child care":    ((236, 295, 355, 398), GREEN),
    "Inadequate family planning":            ((361, 295, 480, 398), GREEN),  # DHS fertility and reproductive health (DRAFT)
    "Limited access to services":            ((486, 295, 611, 398), GREEN),
    "Inadequate health/nutrition knowledge": ((616, 295, 741, 398), AMBER),  # schooling yes, nutrition knowledge no (DRAFT)
    "Inadequate WASH":                       ((747, 295, 872, 398), GREEN),
    "Inadequate intake, absorption":         ((99, 443, 304, 500), AMBER),   # diet yes, absorption no
    "Infectious / chronic disease":          ((366, 443, 520, 500), GREEN),
    "Obesity":                               ((565, 443, 679, 500), GREEN),  # IHME overweight, DHS adult BMI (DRAFT)
    "Gynecologic and obstetric":             ((724, 444, 872, 501), AMBER),
    "Other micronutrient deficiencies":      ((142, 535, 284, 580), AMBER),  # modelled anaemia, national fortification (DRAFT)
    "Inherited red blood cell disorders":    ((304, 535, 446, 580), GREEN),
    "Inflammation":                          ((463, 535, 645, 580), AMBER),  # infection only, no CRP / AGP
    "Blood loss":                            ((656, 535, 838, 580), AMBER),  # helminths only
}
OUTCOME = (415, 615, 597, 684)
S = 2          # upscale factor
PAD = 22       # room above the top row for its badges, and at the right edge (base pixels)


def build(with_badges):
    base = Image.open(BASE).convert("RGB")
    bw, bh = base.size
    img = Image.new("RGB", ((bw + PAD) * S, (bh + PAD) * S), "white")
    img.paste(base.resize((bw * S, bh * S), Image.LANCZOS), (0, PAD * S))
    if with_badges:
        d = ImageDraw.Draw(img)
        r = 10 * S
        for (x0, y0, x1, y1), col in BOXES.values():
            cx, cy = (x1 - 2) * S, (y0 + PAD - 4) * S          # on the top-right corner, nudged up off the box text
            d.ellipse([cx - r - 5, cy - r - 5, cx + r + 5, cy + r + 5], fill="white")
            d.ellipse([cx - r, cy - r, cx + r, cy + r], fill=col)
        x0, y0, x1, y1 = OUTCOME   # the outcome the surveys measured
        d.rectangle([x0 * S - 5, (y0 + PAD) * S - 5, x1 * S + 5, (y1 + PAD) * S + 5], outline=TEAL, width=8)
    return img


for flag, name in [(False, "v2_framework_plain.png"), (True, "v2_framework_coverage.png")]:
    im = build(flag)
    im.save(os.path.join(OUT, name), dpi=(300, 300))
    print("wrote", name, im.size)
