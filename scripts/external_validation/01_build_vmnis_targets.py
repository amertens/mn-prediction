"""
scripts/external_validation/01_build_vmnis_targets.py            [XV-01]

HELD-OUT ADMIN-1 TARGETS FOR FOUR NEW COUNTRIES, FROM WHO VMNIS

The transport claim (P1, PREREGISTRATION_NEW_COUNTRIES_2026-09.md) has never
been scored against a country outside the four-country panel. The one existing
external check -- scripts/policy_deck/24_civ_2007_survey_check.R, women's B12
over nine Cote d'Ivoire eco-regions, rank agreement 0.95 -- was a one-off on a
single outcome with a compass crosswalk.

VMNIS deposits sub-national rows (6,803 at "1st administration level", 3,665
"regional (within country)") that the project's national-only pull
(scripts/pull_vmnis_validation.R, Representativeness == "national") discards.
Those rows are a held-out label set for exactly the estimand the project
claims: WITHIN-country ranking, which is invariant to the cross-survey level
offset that defeats level transport, and to assay/cut-off differences, which
are constant inside one survey.

Four countries clear the bar (>= 6 units, populated prevalence, our
populations, Africa so the pre-registered climate+soil index runs unmodified):

  Zambia   2023   9 provinces   (1st admin level)
  Ethiopia 2015  11 regions     (1st admin level)
  Sudan    2018  15 states      (1st admin level)
  Nigeria  2021   6 zones       (regional within country; NFCMS 2021 was
                                 deposited at zone, not state, resolution)

Boundary vintage is reconciled to each survey's own framework, not to GADM's:
Zambia's Muchinga (created 2011) is returned to Northern except Chama
(Eastern); Sudan's Central Darfur, East Darfur and West Kurdufan (created
2012-13) are returned to their pre-split parents; Abyei is dropped as disputed.
Every one of those is recorded in the crosswalk with a `note`.

  python scripts/external_validation/01_build_vmnis_targets.py
-> data/external_validation/vmnis_admin1_targets.csv
   data/external_validation/gadm_to_vmnis_crosswalk.csv
"""
import os, sys, zipfile, warnings
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
ZIP = os.path.join(ROOT, "data/RA_2026-09/VMNIS.zip")
OUT = os.path.join(ROOT, "data/external_validation")
os.makedirs(OUT, exist_ok=True)

# survey year comes from metadata/survey_years.csv (CLAUDE.md: ONE place)
SY = pd.read_csv(os.path.join(ROOT, "metadata/survey_years.csv"))
SY = SY.set_index("country")

TARGETS = {
    # Africa arm (XV-01): iSDAsoil exists, so the pre-registered climate + soil
    # index runs unmodified
    "Zambia":   dict(iso="ZMB", rep="1st administration level", arm="africa"),
    "Ethiopia": dict(iso="ETH", rep="1st administration level", arm="africa"),
    "Sudan":    dict(iso="SDN", rep="1st administration level", arm="africa"),
    "Nigeria":  dict(iso="NGA", rep="regional (within country)", arm="africa"),
    # off-continent arm (XV-02): iSDAsoil is Africa-only, so the soil half is
    # substituted with global SoilGrids and the comparison is run BOTH ways on
    # the Africa arm to price the substitution
    "Pakistan": dict(iso="PAK", rep="1st administration level", arm="offcontinent"),
    # India's CNNS sampled preschool/school-age children and adolescents, no
    # women, so only the child outcomes can be scored there
    "India":    dict(iso="IND", rep="1st administration level", arm="offcontinent"),
}

# VMNIS indicator x population -> the panel outcome it can be scored against
OUTCOME = {
    ("Retinol (plasma or serum)", "child"): "child_vitA",
    ("Retinol binding protein",   "child"): "child_vitA",
    ("Retinol (plasma or serum)", "women"): "women_vitA",
    ("Retinol binding protein",   "women"): "women_vitA",
    ("Ferritin",                  "child"): "child_iron",
    ("Ferritin",                  "women"): "women_iron",
    ("Vitamin B12",               "women"): "women_b12",
    ("Folate (plasma or serum)",  "women"): "women_folate",
    ("Folate (red blood cell)",   "women"): "women_folate",
    ("Zinc (plasma or serum)",    "child"): "child_zinc",
    ("Zinc (plasma or serum)",    "women"): "women_zinc",
}
# case variants of one unit inside a single deposit
ALIAS = {"North west": "North West"}
POPGROUP = {"Preschool-age children": "child",
            "Non-pregnant women (NPW)": "women",
            "Women of reproductive age": "women"}

INDICATORS = sorted({k[0] for k in OUTCOME})


def load_vmnis():
    """Read the per-indicator Export sheets. Prevalence lives in an
    indicator-specific column ('Depleted iron stores prevalence',
    'Inadequacy prevalence', ...), which is why a generic
    'Prevalenceofdeficiency' pull finds nothing sub-national."""
    z = zipfile.ZipFile(ZIP)
    frames = []
    for info in z.infolist():
        base = os.path.basename(info.filename)
        if not base.startswith("VMNISIndicator_"):
            continue
        ind = base.replace("VMNISIndicator_", "").rsplit("_", 1)[0]
        if ind not in INDICATORS:
            continue
        with z.open(info) as fh:
            d = pd.read_excel(fh, sheet_name="Export")
        pcols = [c for c in d.columns
                 if "prevalence" in c.lower() and "cut-off" not in c.lower()
                 and not any(k in c.lower() for k in
                             ("marginal", "severe", "elevated", "excess"))]
        d["prev"] = d[pcols].bfill(axis=1).iloc[:, 0] if pcols else pd.NA
        d["level"] = d[["Mean", "Geo mean", "Median"]].bfill(axis=1).iloc[:, 0]
        d["indicator"] = ind
        frames.append(d[["Country", "Begin year", "Representativeness.",
                         "Representativeness name", "Population", "Area covered",
                         "Sample size", "indicator", "prev", "level"]])
    return pd.concat(frames, ignore_index=True)


def main():
    V = load_vmnis()
    rows = []
    for country, spec in TARGETS.items():
        year = int(SY.loc[country, "survey_year"])
        d = V[(V.Country == country) & (V["Begin year"] == year)
              & (V["Representativeness."] == spec["rep"])].copy()
        # Keep VMNIS's own spelling -- the crosswalk in script 02 is written to
        # match it. Only collapse whitespace and fix the one case variant in the
        # Nigeria deposit ('North West' / 'North west'); blanket .title() would
        # silently turn Ethiopia's 'SNNP' into 'Snnp' and drop that region.
        d["unit"] = (d["Representativeness name"].astype(str)
                     .str.strip().str.replace(r"\s+", " ", regex=True)
                     .replace(ALIAS))
        d["popgroup"] = d["Population"].map(POPGROUP)
        d = d[d.popgroup.notna()]
        d["outcome"] = [OUTCOME.get((i, p)) for i, p in zip(d.indicator, d.popgroup)]
        d = d[d.outcome.notna()]
        # urban/rural splits would duplicate a unit; keep the combined row
        both = d[d["Area covered"] == "both urban and rural"]
        d = both if len(both) else d
        # one row per unit x outcome: prefer a row carrying prevalence, then median
        d = (d.sort_values("prev", na_position="last")
               .groupby(["unit", "outcome", "indicator"], as_index=False)
               .agg(prev=("prev", "median"), level=("level", "median"),
                    n=("Sample size", "median")))
        d["country"], d["iso3"], d["survey_year"] = country, spec["iso"], year
        d["arm"] = spec["arm"]
        rows.append(d)
    T = pd.concat(rows, ignore_index=True)
    # where two assays serve one outcome (retinol vs RBP, serum vs RBC folate),
    # keep the one with the most units covered, then the most prevalences
    pick = (T.groupby(["country", "outcome", "indicator"])
              .agg(units=("unit", "nunique"), prevs=("prev", lambda s: s.notna().sum()))
              .reset_index()
              .sort_values(["country", "outcome", "units", "prevs"], ascending=False)
              .drop_duplicates(["country", "outcome"]))
    T = T.merge(pick[["country", "outcome", "indicator"]],
                on=["country", "outcome", "indicator"])
    T = T[["country", "iso3", "arm", "survey_year", "unit", "outcome", "indicator",
           "prev", "level", "n"]].sort_values(["country", "outcome", "unit"])
    T.to_csv(os.path.join(OUT, "vmnis_admin1_targets.csv"), index=False)
    print(f"wrote vmnis_admin1_targets.csv: {len(T)} rows")
    print(T.groupby(["country", "outcome"])
           .agg(units=("unit", "nunique"), prev=("prev", lambda s: s.notna().sum()),
                lvl=("level", lambda s: s.notna().sum())).to_string())


if __name__ == "__main__":
    main()
