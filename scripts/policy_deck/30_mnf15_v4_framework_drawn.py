"""Sonja Hess's conceptual framework of iron deficiency, REDRAWN in code so each
box can be filled by whether public district-level data reach it (v4 MNF15 talk,
docs/slides/MNF15-talk-2026-09-v4.qmd, the coverage slide).

    python scripts/policy_deck/30_mnf15_v4_framework_drawn.py

Layout and wording follow the slide Sonja sent on 23 September (Hess et al.,
Ann N Y Acad Sci 2023, adapted); box positions are the ones measured from that
image for scripts/policy_deck/27_mnf15_v2_framework.py (927 x 696 pixels), so
the redrawn figure keeps her geometry. The v2/v3 badges become box fills, with
the "yes" call split in two, as Andrew asked on 27 September:
  AREA  available, and area-level by nature (climate, ecology, built-up land):
        the data describe the place itself
  AGG   available as a district average standing in for individuals (poverty,
        schooling, water, infection): ecological, so an indirect proxy for what
        happens to any one person
  PART  partly: national or modelled surfaces only, or only part of the box
  NONE  not publicly available
The calls that changed from v3: inherited red blood cell disorders move from
yes to partly (the only layers are the Malaria Atlas 2012 modelled
allele-frequency surfaces for HbS, HbC and G6PD); genetic risk factors move to
none (28 September: those three layers belong to the red cell box, and nothing
public covers other genetic risk, so colouring both boxes counted them twice); blood loss stays
partly (ESPEN hookworm and schistosomiasis prevalence only, nothing on
menstrual or obstetric loss). Calls still marked DRAFT for Sonja: built
environment, family planning, obesity, knowledge, other micronutrient
deficiencies.

28 September review, every box against the layers on disk. DHS counts as
public (the headline model holds it out, but that is a modelling choice, not
a gap in the data), so family planning, obesity (women's BMI), obstetric
history and pregnancy keep their fills. Changed: health / nutrition knowledge
to none (the only layers are schooling and literacy, already counted under
educational attainment); two boxes only part of which is measured are drawn
in two fills (SPLITS): fundamental drivers (climate, ecology, geography
measured everywhere; politics is ACLED conflict events only, economy district
wealth averages, inequity one layer, the within-district spread of Meta's
relative wealth index) and infectious / chronic disease (infection well
covered; chronic disease nothing beyond HIV).

Writes results/figures/mnf15_v4/v4_framework_coverage.png. The legend is
editable text on the slide (docs/slides/MNF15-talk-2026-09-v4.concept.yaml),
so its colours are printed below for the YAML.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
OUT = os.path.join(ROOT, "results/figures/mnf15_v4")
os.makedirs(OUT, exist_ok=True)

FILL = {"AREA": "#2E7D32", "AGG": "#B9E0B0", "PART": "#F6D27A", "NONE": "#E3E3E3", "OUT": "#0F7B8A"}
TEXT = {"AREA": "white", "AGG": "#1A1A1A", "PART": "#1A1A1A", "NONE": "#666666", "OUT": "white"}

# (label, (x0, y0, x1, y1) in the 927 x 696 base image, call)
BOXES = [
    ("Physiological vulnerability of women,\nchildren and adolescents\u00b9", (152, 8, 593, 77), "PART"),  # age structure; pregnancy (DHS)
    ("Genetic risk factors", (593, 8, 814, 77), "NONE"),                        # HbS/HbC/G6PD counted under inherited red cell disorders
    ("Poverty", (153, 179, 278, 258), "AGG"),
    ("Low educational\nattainment", (287, 179, 412, 258), "AGG"),
    ("Cultural norms\nand behaviors", (420, 179, 545, 258), "NONE"),
    ("Inadequate built\nenvironment", (554, 179, 679, 258), "AREA"),          # built-up area, night lights (DRAFT)
    ("Health policies", (690, 179, 815, 258), "PART"),                         # national fortification only
    ("Food insecurity\n(quality and\nquantity)", (111, 295, 230, 398), "AGG"),
    ("Inadequate\nmaternal and\nchild care", (236, 295, 355, 398), "AGG"),
    ("Inadequate\nfamily planning\u00b2", (361, 295, 480, 398), "AGG"),       # DHS fertility (DRAFT)
    ("Limited access to\nhealth / nutrition\nservices and\ninterventions", (486, 295, 611, 398), "AGG"),
    ("Inadequate\nhealth / nutrition\nknowledge and\neducation", (616, 295, 741, 398), "NONE"),   # schooling only, already under education (DRAFT)
    ("Inadequate access\nto water, sanitation\nand hygiene", (747, 295, 872, 398), "AGG"),
    ("Inadequate nutrient intake,\nabsorption, utilization", (99, 443, 304, 500), "PART"),        # diet yes, absorption no
    ("Obesity", (565, 443, 679, 500), "AGG"),                                  # IHME overweight, DHS BMI (DRAFT)
    ("Gynecologic and\nobstetric conditions", (724, 444, 872, 501), "PART"),
    ("Other micronutrient\ndeficiencies", (142, 535, 284, 580), "PART"),      # modelled anaemia, fortification (DRAFT)
    ("Inherited red blood\ncell disorders\u00b3", (304, 535, 446, 580), "PART"),   # modelled allele frequencies only
    ("Inflammation", (463, 535, 645, 580), "PART"),                            # infection only, no CRP / AGP
    ("Blood loss", (656, 535, 838, 580), "PART"),                              # hookworm, schistosomiasis only
    ("Iron deficiency", (415, 615, 597, 684), "OUT"),
]
# Boxes public data reach only in part, drawn in two fills:
# ((first label, second label), box, (first call, second call), split, "v" left|right at x or "h" top/bottom at y)
SPLITS = [
    (("Ecology, climate, geography", "Politics, economy, inequity"), (152, 111, 814, 156), ("AREA", "PART"), 483, "v"),
    (("Infectious disease", "Chronic disease"), (366, 443, 520, 500), ("AGG", "PART"), 471.5, "h"),   # chronic: HIV only
]
ROWS = [("Intrinsic\nrisk factors", 42, 90), ("Fundamental\ndrivers", 133, 168), ("Underlying\nrisk factors", 218, 278),
        ("Intermediate\nrisk factors", 346, 420), ("Direct\ncauses", 510, None)]

W, H = 927, 700
fig = plt.figure(figsize=(W / 100, H / 100), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")

B = {}
for lab, (x0, y0, x1, y1), call in BOXES:
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=FILL[call], edgecolor="#4D4D4D", linewidth=0.8, zorder=2))
    size = 9.6 if call == "OUT" else (8.2 if lab.count("\n") >= 3 else 8.6)
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, lab, ha="center", va="center", fontsize=size, color=TEXT[call],
            fontweight="bold" if call in ("OUT",) else "normal", linespacing=1.15, zorder=3)
    B[lab.split("\n")[0]] = (x0, y0, x1, y1)

for labs, (x0, y0, x1, y1), calls, s, how in SPLITS:
    parts = [(x0, y0, s, y1), (s, y0, x1, y1)] if how == "v" else [(x0, y0, x1, s), (x0, s, x1, y1)]
    for lab, (a0, b0, a1, b1), call in zip(labs, parts, calls):
        ax.add_patch(Rectangle((a0, b0), a1 - a0, b1 - b0, facecolor=FILL[call], edgecolor="none", zorder=2))
        ax.text((a0 + a1) / 2, (b0 + b1) / 2, lab, ha="center", va="center", fontsize=8.6, color=TEXT[call], zorder=3)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="none", edgecolor="#4D4D4D", linewidth=0.8, zorder=3))
    ax.plot([s, s] if how == "v" else [x0, x1], [y0, y1] if how == "v" else [s, s], color="#4D4D4D", linewidth=0.5, zorder=3)

for lab, yc, sep in ROWS:   # row labels and the dotted separators of the original
    ax.text(6, yc, lab, ha="left", va="center", fontsize=8.4, style="italic", color="#333333")
    if sep is not None:
        ax.plot([6, 100], [sep, sep], linestyle=(0, (1, 2)), color="#777777", linewidth=0.8)


def arrow(p, q, both=False):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>" if both else "-|>", mutation_scale=9, color="#333333",
                                 linewidth=0.9, shrinkA=0, shrinkB=0, zorder=1))


def line(xs, ys):
    ax.plot(xs, ys, color="#333333", linewidth=0.9, zorder=1)


arrow((483, 156), (483, 179))                     # fundamental -> underlying
arrow((483, 258), (483, 295))                     # underlying -> intermediate
line([483, 483], [398, 420]); line([201, 798], [420, 420])   # intermediate -> the direct-cause bus
for x, top in [(201, 443), (443, 443), (622, 443), (798, 444)]:
    arrow((x, 420), (x, top))
arrow((304, 471), (366, 471), both=True)          # intake <-> infection
arrow((213, 500), (213, 535), both=True)          # intake <-> other micronutrient deficiencies
line([120, 120], [500, 650]); arrow((120, 650), (415, 650))   # intake -> iron deficiency
arrow((443, 500), (520, 535))                     # infection -> inflammation
arrow((600, 500), (575, 535))                     # obesity -> inflammation
arrow((470, 500), (720, 535))                     # infection -> blood loss (hookworm, schistosomiasis)
arrow((798, 501), (760, 535))                     # gynecologic and obstetric -> blood loss
arrow((213, 580), (430, 628))                     # other micronutrient deficiencies -> iron
arrow((375, 580), (440, 615))                     # inherited red cell disorders -> iron
arrow((545, 580), (530, 615), both=True)          # inflammation <-> iron
arrow((747, 580), (597, 640))                     # blood loss -> iron

fig.savefig(os.path.join(OUT, "v4_framework_coverage.png"), dpi=300, facecolor="white")
print("wrote v4_framework_coverage.png")
print("legend colours for the YAML:", {k: v.lstrip("#") for k, v in FILL.items()})
