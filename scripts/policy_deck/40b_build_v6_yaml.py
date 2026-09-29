"""Build docs/slides/MNF15-talk-2026-09-v6.concept.yaml from the v5 YAML (see script 40
for the v6 changes): the outline worded as the talk now runs, the variable-importance
slide as figure plus findings, the B12 error card on the post-fix figure, and one
closing slide that mirrors the outline.

    python scripts/policy_deck/40b_build_v6_yaml.py
"""
import re

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
Y5 = ROOT + "docs/slides/MNF15-talk-2026-09-v5.concept.yaml"
Y6 = ROOT + "docs/slides/MNF15-talk-2026-09-v6.concept.yaml"

y = open(Y5, encoding="utf-8").read()
head = y[:y.index("slides:")]
entries = {}
for m in re.finditer(r"^- title: ([^\n]+)\n(.*?)(?=^- title: |^# ----|\Z)", y, re.S | re.M):
    entries[m.group(1).strip()] = "- title: " + m.group(1).strip() + "\n" + m.group(2).rstrip() + "\n"


def get(title, subs=()):
    e = entries[title]
    for a, b in subs:
        assert a in e, (title, a[:60])
        e = e.replace(a, b)
    return e


QR = '''  image: {file: img/qr_dashboard.png, width: 1.9}
  texts:
  - {x: 10.25, y: 3.95, w: 2.45, h: 1.0, paras: [{text: "Every district, five countries:", size: 11, color: "555555"}, {text: "amertens.shinyapps.io/", size: 11, color: "555555"}, {text: "micronutrient-burden", size: 11, color: "555555"}]}
'''

out = [head.replace("MNF15-talk-2026-09-v5", "MNF15-talk-2026-09-v6").replace("version 5", "version 6"), "slides:\n\n"]
out.append(get("title slide"))
out.append(get("Introduction and study objectives"))
out.append('''- title: Talk outline
  layout: icon_rows
  heading_size: 20
  text_size: 14
  items:
  - icon: vials
    heading: "Can public data tell us how common a deficiency is, nationally or by district?"
  - icon: ranking-star
    heading: "Can they rank a country's districts from most to least deficient?"
  - icon: earth-africa
    heading: "Does that work in a country that has never had a survey?"
  - icon: satellite
    heading: "Which data matter, and for which nutrients?"
  - icon: map-location-dot
    heading: "How could programmes and survey planners use this?"
''')
_unused_outline = (get("Talk outline", subs=[
    ('heading: A simple model does best\n    caption: "There are few places to learn from, so the simplest model wins"',
     'heading: A simple model does as well as complex ones\n    caption: "There are few places to learn from"'),
    ('heading: It ranks districts no survey reached\n    caption: "Tested on districts, regions and countries it never saw"',
     'heading: It ranks districts, even where no survey reached\n    caption: "Tested on districts, regions and countries it never saw; levels still need a survey"')]))
out.append(get("Conceptual framework of iron deficiency"))
out.append(get("Which of these causes can public data measure?"))
out.append(get("Why not a more complex machine-learning model?"))
out.append('''- title: What variables drive model predictions
  layout: freeform
  texts:
  - x: 0.45
    y: 1.8
    w: 3.75
    h: 5.2
    space_before: 9
    paras:
    - {text: "Climate and soil carry the cross-border ranking", size: 16, bold: true, color: "0F7B8A"}
    - {text: "Drop either and the ranking of a country held out gets worse, in all four countries. Climate alone does about as well as all 575 layers.", size: 13.5}
    - {text: "The same families lead almost every model", size: 16, bold: true, color: "0F7B8A"}
    - {text: "Climate, soil and a satellite imagery summary are in the top five of all 22 cross-country models and 21 of 24 single-country models.", size: 13.5}
    - {text: "Only B12 has a map of its own", size: 16, bold: true, color: "0F7B8A"}
    - {text: "B12 follows meat or fish in children's diets, and wealth. The iron and vitamin A maps mark largely the same poorer districts. Markers of where, not causes.", size: 13.5}
  images:
  - {file: ../../results/figures/mnf15_v5/v4_domains_child_iron.png, x: 4.35, y: 1.85, w: 8.75}
''')
out.append('''- title: Which deficiencies can it map?
  layout: freeform
  images:
  - {file: ../../results/figures/mnf15_v5/v2_nutrients.png, x: 1.3, y: 1.55, w: 10.7}
  texts:
  - {x: 1.3, y: 6.15, w: 10.7, h: 0.6, paras: [{text: "Only the B12 map is specific to its nutrient: the iron and vitamin A maps mark largely the same poorer districts.", size: 15, bold: true, color: "0F7B8A"}]}
''')
out.append(get("What does this mean for policy? The vitamin B12 model", subs=[
    ('heading: About 6 points\n    lines:\n    - "Typical error in district prevalence, each district hidden"\n    - "90% ranges about 21 to 27 points either side"',
     'heading: About 9 points\n    lines:\n    - "Typical error in district prevalence, each district hidden"\n    - "About as close as the survey\'s own regional figures"')]))
out[-1] = out[-1].replace('"65 to 71 with the country\'s own survey removed"', '"65 to 71 with the country\'s own survey removed (average B12 level)"').replace('"The fifth of districts it ranks worst hold 39 to 46% of deficient women; a random fifth holds 20%"', '"Inside a surveyed country, the fifth of districts it ranks worst hold 39 to 46% of deficient women; a random fifth holds 20%"')
assert 'average B12 level' in out[-1] and 'Inside a surveyed country' in out[-1]
out.append(get("Next steps"))
out.append('''- title: Conclusions
  layout: icon_rows
  heading_size: 16
  text_size: 13
''' + QR + '''  items:
  - icon: vials
    heading: "How common is a deficiency? Not from public data alone"
    caption: "Without a survey, national estimates missed by up to 38 points; district levels need the survey's national figure."
  - icon: ranking-star
    heading: "Ranking districts? Yes, modestly"
    caption: "66 of 100 district pairs in the survey's order, about as well as the survey's own regional figures (65)."
  - icon: earth-africa
    heading: "With no survey? Mostly"
    caption: "62 of 100 pairs, and the right direction in 37 of 41 checks in six more countries; not for every nutrient or country."
  - icon: satellite
    heading: "Which data? Climate and soil above all"
    caption: "They carry the cross-border ranking in all four countries. Only B12 has a map of its own; the iron and vitamin A maps mark the same poorer districts."
  - icon: map-location-dot
    heading: "For programmes: match the ranking to the budget"
    caption: "Per person, rank by the model (31% of the deficient reached with a fifth of the people; random 20%). Per district, start with the most populous. It does not replace a survey."
''')
out.append("\n# ---- appendix\n")
out.append(get("How does the model work, and how was it tested?"))
out.append(get("How can it help plan the next survey?", subs=[
    ("About 40 blood draws per nutrient anchor the ranking; district figures then typically land within about 10 points.", "About 40 blood draws per nutrient give the national level; district figures then land within about 10 points, about as close as the national figure alone."),
    ("Ours reaches 0.4 to 0.6.", "Ours reaches 0.4 to 0.6; for zinc, nothing.")]))
out.append(get("What can we learn about the consistently most important predictor domains?", subs=[
    ('heading: Climate, soil and satellite land cover\n      caption: "In the top five of all 22 cross-country models and every single-country model; with greenness, about half the weight."',
     'heading: Climate and soil carry the cross-border ranking\n      caption: "Dropping either makes the held-out ranking worse in all four countries. With a satellite imagery summary they are in the top five of all 22 cross-country models and 21 of 24 single-country models."'),
    ('heading: Each alone does almost as well as everything\n      caption: "Climate or soil on its own ranks a country held out about as well as all 575 layers together."',
     'heading: Climate alone does about as well as everything\n      caption: "Climate on its own ranks a country held out about as well as all 575 layers (better in 3 of 4 countries); soil alone is less consistent."')]))
out.append(get("What can we learn about unused predictor domains?"))
import glob as _glob, os as _os
_txt = "\n".join(e.rstrip() + "\n" for e in out)
for _f in sorted(_os.path.basename(x) for x in _glob.glob(ROOT + "results/figures/mnf15_v6/*.png")):
    _txt = _txt.replace("mnf15_v5/" + _f, "mnf15_v6/" + _f)   # figures redrawn on the post-fix screen
open(Y6, "w", encoding="utf-8").write(_txt)
print("wrote", Y6, len(re.findall(r"^- title:", open(Y6, encoding="utf-8").read(), re.M)), "entries")
