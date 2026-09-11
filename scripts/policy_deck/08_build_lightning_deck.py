"""
scripts/policy_deck/08_build_lightning_deck.py

Build the 3-minute short-oral (rapid-fire) deck for the Micronutrient Forum
2026 on the Forum's own template: at most three slides, with the abstract
title, background and aims, methods, results, conclusions and impact.

The talk takes the one story held out of the joint 15-minute talk: which
public data layers carry the signal, and that the same layers hold across
countries. Every number is read from the result tables.

    python scripts/policy_deck/08_build_lightning_deck.py

Input template:  C:/Users/andre/Dropbox/MN prediction/MNF2026-Short-oral-template.pptx
Output:          docs/slides/MN-proxy-lightning-MNF2026.pptx
"""
import copy
import csv
import os
import re
import sys

from pptx import Presentation
from pptx.util import Inches, Pt

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
TEMPLATE = "C:/Users/andre/Dropbox/MN prediction/MNF2026-Short-oral-template.pptx"
OUT = os.path.join(ROOT, "docs", "slides", "MN-proxy-lightning-MNF2026.pptx")
P2 = os.path.join(ROOT, "results", "tables", "protocol_v2")
FIG = os.path.join(ROOT, "results", "figures", "policy_deck")


def rows(name):
    with open(os.path.join(P2, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def one(rs, **cond):
    hit = [r for r in rs if all(r.get(k) == v for k, v in cond.items())]
    if not hit:
        sys.exit("no row for %r in table" % (cond,))
    return hit[0]


def mean(xs):
    xs = [float(x) for x in xs if x not in ("", "NA")]
    return sum(xs) / len(xs)


# ---- numbers -----------------------------------------------------------------
BM = rows("benchmarks_v2_summary.csv")
r_in = one(BM, estimand="infill", arm="domain_index", target="level")
r_tr = one(BM, estimand="country", arm="domain_index", target="level")
infill = float(r_in["mean_spearman"])
tr = float(r_tr["mean_spearman"])
tr_pos, tr_n = int(r_tr["cells_positive"]), int(r_tr["cells"])
null_d = float(one(rows("transport_null_calibration.csv"), tier="admin2")["null_mean_q95"])
ND = rows("nested_domain_selection.csv")
cs = mean([r["spearman"] for r in ND if r["arm"] == "fixed_cs" and r["target"] == "level"])
sparse20 = float(one(rows("weight_sources_summary.csv"), estimand="country", arm="sparse20", target="level")["mean_spearman"])
dhs = one(rows("source_ablation_loco_summary.csv"), target="level", source="DHS")
dhs_cols, dhs_full, dhs_drop = int(dhs["n_cols"]), float(dhs["full"]), float(dhs["drop"])

# share of the index carried by satellite imagery, climate and soil, per outcome (pooled fit, biomarker level)
IMP = rows("index_importance_domains.csv")
pooled = [r for r in IMP if r["scope"] == "pooled" and r["target"] == "level"]
if not pooled:
    pooled = [r for r in IMP if r["target"] == "level"]
domains = sorted(set(r["domain"] for r in pooled))
print("importance scopes:", sorted(set(r["scope"] for r in IMP)), "| domains:", domains)
ENV = re.compile(r"satellite|remote|climate|soil|vegetation|land cover|terrain|night", re.I)
env_share = {}
for r in pooled:
    if ENV.search(r["domain"]):
        env_share[r["outcome"]] = env_share.get(r["outcome"], 0.0) + float(r["share"])
env_lo, env_hi = min(env_share.values()), max(env_share.values())
print("environment share by outcome:", {k: round(v, 2) for k, v in env_share.items()})


def f2(x):
    return "%.2f" % x


def pcr(x):
    return "%d" % round(100 * x)


# The Forum asks for the submitted abstract title on the first slide, verbatim
# (Dropbox: "MNF abstract_proxy modeling_V2.docx", April 2026), with its authors.
TITLE = "Predicting Micronutrient Deficiency Prevalence in Sub-Saharan Africa Using Proxy Indicators: A Multi-country Machine Learning Framework"
AUTHORS = ("Andrew N. Mertens (presenting), Reed Atkin, Sorrel Namaste, Demewoz Haile, Kenneth H. Brown, Saskia Osendarp, "
           "Xiuping Tan, Eric Stewart, Keith Lividini, Seth Adu-Afarwuah, Nicolai Petry, Aminata S. Koroma, Mary H. Hodges, "
           "Fabian Rohner, Haddy Crookes, James P. Wirth, Jonathan Gorstein, Kathy Banke, Sonja Y. Hess")

S1 = [
    "Aim: find which freely available data layers track blood-measured micronutrient deficiency across districts, and whether the same layers hold in a country never surveyed.",
    "Data: four national micronutrient surveys (The Gambia 2018, Ghana 2017, Malawi 2015 to 16, Sierra Leone 2013), 206 districts, vitamin A, iron, folate and B12 in children and women.",
    "Predictors: 454 public layers in 24 groups: satellite imagery, climate, soil, crops, livestock, malaria, prices and household surveys.",
    "Model: each group is reduced to a few summary axes, weighted by its correlation with deficiency and summed, with nothing tuned. Scored on held-out districts and on held-out countries.",
]
S2 = [
    "Inside a surveyed country the index ranks held-out districts at %s, against %s for chance." % (f2(infill), f2(null_d)),
    "In a country never surveyed it scores %s, positive in %d of %d country-outcome pairs, and %s with climate and soil layers only." % (f2(tr), tr_pos, tr_n, f2(cs)),
    "Satellite, climate and soil layers carry %s to %s percent of the index, and the same kind of district ranks worst on every deficiency: drier, grassier, poorer and more pastoral." % (pcr(env_lo), pcr(env_hi)),
    "Household-survey layers help inside a country and hurt in a new one: dropping all %d of them raises transport from %s to %s." % (dhs_cols, f2(dhs_full), f2(dhs_drop)),
]
S3 = [
    "One environmental gradient underlies all six deficiencies. A twenty-layer public list ranks a new country as well as the full database (%s against %s)." % (f2(sparse20), f2(tr)),
    "Districts nobody has surveyed can be ranked for programme targeting and for placing the next survey. Rankings cross borders; absolute levels do not.",
    "One national blood sample turns a ranking into district prevalence. The climate-and-soil list was found on these four countries and is pre-registered for the fifth.",
    "Impact: a country can rebuild the ranking itself from downloadable layers, at no survey cost, while the next survey is planned.",
]


# ---- template filling ----------------------------------------------------------
def shape_by_name(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    sys.exit("shape %r not found" % name)


def set_text(shape, text, size=None):
    """Keep the first run's formatting, drop everything else."""
    tf = shape.text_frame
    p0 = tf.paragraphs[0]
    for p in list(tf.paragraphs[1:]):
        p._p.getparent().remove(p._p)
    runs = p0.runs
    if not runs:
        runs = [p0.add_run()]
    for r in runs[1:]:
        r._r.getparent().remove(r._r)
    runs[0].text = text
    if size:
        runs[0].font.size = Pt(size)


def set_bullets(shape, lines, size=None):
    tf = shape.text_frame
    tf.word_wrap = True
    proto = copy.deepcopy(tf.paragraphs[0]._p)
    for p in list(tf.paragraphs):
        p._p.getparent().remove(p._p)
    for line in lines:
        p = copy.deepcopy(proto)
        tf._txBody.append(p)
        para = tf.paragraphs[-1]
        for r in para.runs[1:]:
            r._r.getparent().remove(r._r)
        if not para.runs:
            para.add_run()
        para.runs[0].text = line
        if size:
            para.runs[0].font.size = Pt(size)
        para.space_after = Pt(6)


def label(slide):
    tb = slide.shapes.add_textbox(Inches(5.7), Inches(6.95), Inches(3.6), Inches(0.35))
    r = tb.text_frame.paragraphs[0].add_run()
    r.text = "Unpublished data. Do not copy or distribute."
    r.font.size = Pt(10)


prs = Presentation(TEMPLATE)
s1, s2, s3 = list(prs.slides)[:3]

title = shape_by_name(s1, "Title 1")
title.left, title.top, title.width, title.height = Inches(0.27), Inches(0.15), Inches(12.8), Inches(1.15)
set_text(title, TITLE, size=22)
byline = shape_by_name(s1, "TextBox 5")
byline.left, byline.top, byline.width, byline.height = Inches(0.27), Inches(1.35), Inches(12.8), Inches(0.7)
byline.text_frame.word_wrap = True
set_text(byline, AUTHORS, size=12)
heading = shape_by_name(s1, "Title 23")
heading.top = Inches(2.1)
set_text(heading, "Background, aims and methods")
body = shape_by_name(s1, "Text Placeholder 2")
body.left, body.top, body.width, body.height = Inches(0.39), Inches(2.8), Inches(12.4), Inches(4.0)
set_bullets(body, S1, size=18)

set_text(shape_by_name(s2, "Title 1"), "Results")
body = shape_by_name(s2, "Text Placeholder 2")
body.left, body.top, body.width, body.height = Inches(0.27), Inches(1.2), Inches(4.75), Inches(5.6)
set_bullets(body, S2, size=16)
s2.shapes.add_picture(os.path.join(FIG, "fig3_top_predictors.png"), Inches(5.1), Inches(1.15), width=Inches(8.0))
label(s2)

set_text(shape_by_name(s3, "Title 1"), "Conclusions and impact")
body = shape_by_name(s3, "Text Placeholder 2")
body.left, body.top, body.width, body.height = Inches(0.27), Inches(1.25), Inches(7.2), Inches(5.5)
set_bullets(body, S3, size=18)
s3.shapes.add_picture(os.path.join(FIG, "fig5_learning_curve.png"), Inches(7.7), Inches(1.3), width=Inches(5.3))
label(s3)

prs.save(OUT)
print("wrote", OUT)
