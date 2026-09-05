# External check against Tang et al. (Nature Food 2026)

**Date:** 2026-09-02
**Source:** Tang K, Agu U, Mbodj S, et al. *Diversified food fortification
portfolios can enhance coverage of micronutrient-vulnerable populations in West
Africa.* Nature Food (2026). doi:10.1038/s43016-026-01412-2. Open access, CC BY.
WFP Nutrition and Food Quality Service (the MIMI team).
Code: <https://github.com/MIMI-wfp/West-Africa-LSFF-Analysis>

**Scripts:** `sandbox_lsff/01_compare_tang_adm1.R`, `02_prespecified_and_spatial.R`
**Outputs:** `results/tables/tang_lsff_adm1_agreement.csv`,
`results/tables/tang_lsff_prespecified.csv`
**Reference data:** `data/WFP_LSFF_2026/`

## What was compared

Two entirely independent measurements of the same geography:

| | ours | theirs |
|---|---|---|
| quantity | biomarker deficiency prevalence | dietary inadequacy |
| data | DHS/MICS biomarkers + GEE covariates | HCES apparent consumption |
| model | area-level SAE / elastic net | AFE allocation vs Harmonized Average Requirement |

Levels are not comparable and were not compared. The question is **rank
agreement**: do the two approaches order the same areas the same way?

Geography is a lucky break for Côte d'Ivoire: the paper's ADM1 (HDX COD-AB,
33 regions) is the same tier as our **Admin2** (GADM level 2, 33 regions), so
the join is unit-for-unit at our own working resolution. Ghana required
aggregating our 260 Admin2 up to the pre-2018 **10** regions that GLSS 2016/17
uses (crosswalk in `data/WFP_LSFF_2026/`), population-weighted.

## Result

Pre-specified family: our prediction vs the *same nutrient's* % vulnerable,
10 tests, Benjamini-Hochberg. `partial` = both variables rank-residualised on
centroid latitude, because this project already knows that most of its skill
is a north-south gradient.

| country | n | outcome | rho | p | q | rho (partial lat) | p |
|---|---|---|---|---|---|---|---|
| Côte d'Ivoire | 33 | women_vitA | **0.495** | 0.004 | **0.043** | **0.479** | 0.005 |
| Côte d'Ivoire | 33 | child_vitA | 0.395 | 0.024 | 0.118 | 0.323 | 0.067 |
| Ghana | 10 | women_b12 | 0.539 | 0.114 | 0.381 | -0.546 | 0.092 |
| Ghana | 10 | women_folate | -0.418 | 0.232 | 0.390 | -0.368 | 0.293 |
| Ghana | 10 | child_iron | 0.421 | 0.232 | 0.390 | -0.477 | 0.167 |
| Ghana | 10 | women_vitA | 0.415 | 0.234 | 0.390 | -0.026 | 0.945 |
| Côte d'Ivoire | 33 | child_iron | 0.178 | 0.322 | 0.460 | 0.220 | 0.220 |
| Ghana | 10 | women_iron | -0.274 | 0.446 | 0.557 | -0.655 | 0.045 |
| Ghana | 10 | child_vitA | 0.207 | 0.562 | 0.625 | -0.282 | 0.419 |
| Côte d'Ivoire | 33 | women_iron | 0.053 | 0.772 | 0.772 | 0.116 | 0.519 |

### 1. Côte d'Ivoire vitamin A corroborates, and it is not latitude

This is the finding worth keeping. Côte d'Ivoire has **no biomarker survey** in
this project — its Admin-2 estimates are pure out-of-sample transport from
models fitted on the other countries (`oos_civ_*`). Until now there was nothing
to check that ranking against.

An independent, published, differently-sourced estimate agrees at rho ~ 0.4-0.5
across 33 districts, for both women (q = 0.043) and children (same direction,
weaker). It survives latitude adjustment almost unchanged (0.495 -> 0.479),
and the two sources are *not* agreeing by both tracking north-south: our
prediction correlates 0.73 with latitude, theirs only 0.24.

This does not validate the *levels*, which remain untransportable. It supports
the claim the project actually makes — that the transported product is a usable
**priority ordering**.

### 2. Iron corroborates nowhere

rho = 0.05-0.18 in Côte d'Ivoire; in Ghana the sign is unstable and flips under
latitude adjustment. Consistent with what the project already knows about
ferritin-defined iron (cross-survey level offset, adjustment heterogeneity).
Do not present iron transport as externally supported.

### 3. Ghana is uninformative, and instructively so

n = 10 gives essentially no power, and both sources are dominated by the
north-south gradient (each correlates 0.67-0.88 with latitude). Every Ghana
association collapses or reverses once latitude is removed. This is a clean
illustration of the deck's own "what geography alone already gives you" point,
and a warning against reading the raw Ghana numbers.

An earlier cut of this analysis compared against the paper's **overall MPI**
(the five-nutrient average) rather than the matched nutrient, and produced
"Ghana child_iron rho = 0.73, p = 0.02". That is a latitude artifact — the
partial correlation is -0.48. Do not quote it.

## Caveats

- Their MPI is *modelled* dietary inadequacy, not measured intake, and not
  deficiency. Divergence is expected on biological grounds (absorption,
  inflammation, infection) and is not by itself evidence either model is wrong.
- Ghana GLSS is 2016/17 and Côte d'Ivoire EHCVM is 2021/22; our biomarker
  surveys are different vintages again.
- 10 tests here; 28 in the wider exploratory cut, where **nothing** survived
  BH (min q = 0.138). The family definition does real work in this result and
  should be stated whenever the 0.043 is quoted.
