"""Build docs/slides/MNF15-talk-2026-09-v5.concept.yaml from the v4 YAML: the entries of
the slides v5 keeps, with Andrew's text edits from v4-ANM, Sonja's request for more
specific objectives, the post-fix numbers, and the variable-importance findings.

    python scripts/policy_deck/36b_build_v5_yaml.py
"""
import re

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
Y4 = ROOT + "docs/slides/MNF15-talk-2026-09-v4.concept.yaml"
Y5 = ROOT + "docs/slides/MNF15-talk-2026-09-v5.concept.yaml"
B = "\u2022"

y = open(Y4, encoding="utf-8").read()
head = y[:y.index("slides:")]
entries = {}
for m in re.finditer(r"^- title: ([^\n]+)\n(.*?)(?=^- title: |^# ----|\Z)", y, re.S | re.M):
    entries[m.group(1).strip()] = "- title: " + m.group(1).strip() + "\n" + m.group(2).rstrip() + "\n"


def get(title, new_title=None, subs=()):
    e = entries[title]
    if new_title:
        e = e.replace("- title: " + title, "- title: " + new_title, 1)
    for a, b in subs:
        assert a in e, (title, a[:60])
        e = e.replace(a, b)
    return e


out = [head.replace("MNF15-talk-2026-09-v4.qmd", "MNF15-talk-2026-09-v5.qmd").replace("version 4", "version 5"), "slides:\n\n"]
out.append(get("title slide"))

out.append(f'''- title: Introduction and study objectives
  layout: freeform
  texts:
  - x: 0.55
    y: 1.6
    w: 5.4
    h: 5.4
    space_before: 7
    paras:
    - {{text: "{B} Information on the prevalence of vitamin and mineral deficiencies (VMD) is limited", size: 16}}
    - {{text: "{B} In an inter-institutional, multi-disciplinary workshop we brainstormed potential approaches, and relevant data domains and data sources", size: 16}}
    - {{text: "{B} Developed conceptual frameworks to guide machine learning", size: 16}}
    - {{text: "Objectives", size: 17, bold: true, color: "0F7B8A"}}
    - {{text: "{B} To explore whether machine learning can help fill the data gap in VMD prevalence using aggregated proxy variables:", size: 16}}
    - {{text: "   - estimate deficiency prevalence, nationally and by district, validated against biomarker surveys", size: 15}}
    - {{text: "   - identify the districts at greatest risk: a ranking programmes can act on", size: 15}}
  - x: 6.1
    y: 5.78
    w: 6.95
    h: 0.7
    paras:
    - {{text: "Source: WHO Vitamin and Mineral Nutrition Information System (VMNIS), Micronutrients Database, February 2025. Nationally representative surveys measuring ferritin, retinol or RBP, zinc, folate or vitamin B12; surveys not yet deposited are not shown.", size: 10, color: "808080"}}
  images:
  - {{file: ../../results/figures/mnf15_v4/v4_vmnis_coverage.png, x: 6.1, y: 1.5, w: 6.95}}
''')

out.append(get("What we will tell you", "Talk outline", subs=[
    ("heading: Free public data carry real signal about deficiency\n    caption: \"Climate, soil and land cover, measured everywhere, every year\"",
     "heading: Aggregate public data carry real signal about deficiency\n    caption: \"Top predictors: climate, soil and land cover\"")]))
out.append(get("Conceptual framework of iron deficiency"))
out.append(get("Which of these causes can public data measure?", subs=[
    ('"Available, intrinsically ecological: describes the place (climate, land, built-up area)"', '"Available, intrinsically ecological (climate, land, built-up area)"'),
    ('"Available as a district average standing in for individuals (poverty, schooling, infection)"', '"Individual data linkable to VMD as a district aggregate (poverty, schooling, infection)"'),
    ('"Partly: national or modelled only, or part of the box"', '"Partly: national or modelled only, or only some surrogate variables"'),
    ('"The outcome the surveys measured"', '"Survey outcome"')]))
out.append(get("Why not a more complex machine-learning model?", subs=[
    ('{text: "Cross-validation: every score comes from districts the model never saw, the same hidden districts for every method.", size: 16}',
     '{text: "We compared 14 machine-learning methods on the same held-out districts, including an ensemble (SuperLearner) of twelve.", size: 16}'),
    ('{text: "With so few places, methods that tune themselves chase noise: the untuned Domain-PC index matched or beat them all, including an ensemble of twelve.", size: 16, bold: true, color: "0F7B8A"}',
     '{text: "A simple model did as well as complex ones: none was more than slightly better on average across nutrients.", size: 16, bold: true, color: "0F7B8A"}'),
    ("../../results/figures/mnf15_v4/v4_ml_comparison.png", "../../results/figures/mnf15_v5/v4_ml_comparison.png")]))

out.append('''- title: What can we learn about the consistently most important predictor domains?
  layout: two_panel
  heading_size: 15
  text_size: 12.5
  panels:
  - icon: earth-africa
    title: In every model, every nutrient
    items:
    - icon: cloud-sun-rain
      heading: Climate, soil and satellite land cover
      caption: "In the top five of all 22 cross-country models and every single-country model; with greenness, about half the weight."
    - icon: mound
      heading: Each alone does almost as well as everything
      caption: "Climate or soil on its own ranks a country held out about as well as all 575 layers together."
    - icon: wheat-awn
      heading: The same markers, the same direction
      caption: "Grassland, ruminants and cereal cropping mark more deficiency; plant productivity, market access and child overweight mark less."
  - icon: magnifying-glass-chart
    title: Markers for particular nutrients
    items:
    - icon: drumstick-bite
      heading: B12
      caption: "Young children eating meat or fish, and household wealth, mark better B12 status."
    - icon: droplet
      heading: Iron and folate
      caption: "Modelled anaemia marks more deficiency; the farming system matters for iron."
    - icon: child
      heading: Vitamin A
      caption: "Child wasting and underweight, and a young population, mark more deficiency."
    - icon: mosquito
      heading: Malaria, the other way round
      caption: "Marks less iron and vitamin A deficiency, most likely because inflammation shifts the biomarkers."
''')
out.append(get("What does this mean for policy? The vitamin B12 model", subs=[
    ('"The fifth of districts it ranks worst hold 40 to 46% of deficient women; a random fifth holds 20%"',
     '"The fifth of districts it ranks worst hold 39 to 46% of deficient women; a random fifth holds 20%"')]))
out.append(get("Conclusions"))
out.append(get("What have we learned, and what comes next?", subs=[
    ('"Inside a surveyed country, 67 of 100 district pairs in the survey\'s order; its regional averages, 62."',
     '"Inside a surveyed country, 66 of 100 district pairs in the survey\'s order; its regional averages, 61."')]))
out.append(get("Next steps"))
out.append("\n# ---- appendix\n")
out.append(get("How does the model work, and how was it tested?"))
out.append(get("How can it help plan the next survey?", subs=[
    ("../../results/figures/mnf15_v4/v3_ceiling_by_nutrient.png", "../../results/figures/mnf15_v5/v3_ceiling_by_nutrient.png")]))
out.append(get("What can we learn about missing predictor domains?", "What can we learn about unused predictor domains?"))
open(Y5, "w", encoding="utf-8").write("\n".join(e.rstrip() + "\n" for e in out))
print("wrote", Y5, len(re.findall(r"^- title:", open(Y5, encoding="utf-8").read(), re.M)), "entries")
