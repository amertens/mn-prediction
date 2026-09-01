# WS-C3. The iron-anaemia bridge

The one covariate model in this project with a physiological mechanism rather
than a correlational one, and the transport test the session was built around.

`results/tables/iron_bridge.csv`, 8 cells (child and women iron, four
countries), region-blocked LORO within country and leave-one-country-out across.

## The pre-specified set

Nine covariates, one per mechanistic concept, named rather than filtered:

| Concept | Column | Source |
|:---|:---|:---|
| Child haemoglobin | `np_child_hb_gdl` | prior DHS round, WS-F2 |
| Women haemoglobin | `np_women_hb_gdl` | prior DHS round, WS-F2 |
| Child anaemia | `np_child_anaemia_any` | prior DHS round, WS-F2 |
| Women anaemia | `np_women_anaemia_any` | prior DHS round, WS-F2 |
| Iron in pregnancy | `np_iron_pregnancy` | prior DHS round, WS-F2 |
| Malaria parasite rate | `map2_explorer_2020_global_pfpr` | MAP, via the individual merge |
| Sickle HbS allele | `map2_blood_disorders_..._hbs_al` | MAP |
| HbC allele | `map2_blood_disorders_..._hbc_allele_frequency` | MAP |
| Bednet coverage | `dhs_hh_itn_any` | harmonised 294 |

All nine resolve in all four countries.

**A first version of this script arrived at 31 covariates** by accepting every
MAP layer that passed a filter, including Duffy negativity, G6PD deficiency,
temperature suitability and non-malarial fever. A 31-column model with a top-8
screen is exploratory, not the mechanism test the workstream specifies. The set
is now named. That first version also accepted a MAP column only where it was
constant within a district, which kept 22 columns for three countries and
**one** for Malawi, whose layers are extracted at the cluster; the district mean
is taken instead, and the number of columns that genuinely vary within a
district is counted rather than assumed.

## The result: the bridge does not beat the 294, and on average does not work

Pooled over the 8 cells:

| Evaluation | Arm | Mean r | Median r | MAE pp | Predictors |
|:---|:---|---:|---:|---:|---:|
| LORO | **bridge** | **0.027** | -0.017 | 10.70 | 9 |
| LORO | flat regional mean (jackknife) | 0.181 | 0.192 | 10.19 | 0 |
| LORO | full 294 | **0.364** | 0.372 | 9.84 | ~294 |
| LOCO | **bridge** | **-0.022** | -0.064 | 28.05 | 9 |
| LOCO | full 294 | **0.210** | 0.257 | 26.55 | 294 |

The mechanistic set is worse than the 294 both within and across countries, and
worse than a covariate-free jackknifed regional mean within country.

## The heterogeneity the pooled means hide, and it matters

Eight cells is too few for a mean to summarise. The two countries whose
crosswalk is exact disagree with each other.

**Gambia (crosswalk 30 of 30 exact):**

| Evaluation | Outcome | Bridge (9) | Flat regional (jackknife) | Full 294 |
|:---|:---|---:|---:|---:|
| LORO | women_iron | **0.509** | 0.441 | 0.748 |
| LORO | child_iron | 0.179 | 0.119 | 0.600 |
| LOCO | women_iron | **0.459** | not applicable | 0.467 |
| LOCO | child_iron | -0.093 | not applicable | 0.557 |

For Gambia women's iron the bridge reaches **0.459 across countries with nine
covariates against 0.467 for two hundred and ninety-four**. That is the
transport property the workstream was looking for, and it appears in exactly one
of the eight cells.

**Sierra Leone (crosswalk 14 of 14 exact):** the bridge is negative in all four
of its cells (-0.397, -0.092 under LORO; -0.549, -0.147 under LOCO) while the
full 294 is positive under LORO (0.402, 0.341).

Two countries with a clean crosswalk, one where the bridge works for women's
iron and one where it fails everywhere. No conclusion about the mechanism
survives that.

**One place the bridge is clearly better: level transport.** Sierra Leone
child_iron under LOCO has bridge MAE **4.94 pp** against **66.56 pp** for the
294. A nine-covariate model does not extrapolate wildly; a 294-covariate one
does. The pooled MAE hides this because Gambia runs the other way (49.63 against
40.94).

## The limitation that weakens the test, stated because it is large

The bridge covariates come from the prior DHS round, whose district system is
not the analytic one. Measured:

| Country | DHS districts | Analytic districts | Matched | DHS districts discarded |
|:---|---:|---:|---:|---:|
| Gambia | 35 | 30 | 30 | 5 |
| Ghana | 219 | 75 | **69** | **150** |
| Sierra Leone | 14 | 14 | 14 | 0 |
| Malawi | 219 | 87 | 87 | **132** |

For Ghana and Malawi roughly two thirds of the DHS districts are discarded, and
each analytic district takes the value of the single same-named DHS district
rather than a population-weighted aggregate of the DHS districts inside it. That
injects measurement error into the bridge covariates for the two largest
countries, and it is a plausible partial explanation for the null there. It does
not explain Sierra Leone, whose crosswalk is exact and whose bridge is negative.

Building a proper crosswalk needs a DHS-district to analytic-district
correspondence that does not exist on disk. It is recorded as the fix that would
make this test fair rather than attempted here.

## What this does not contradict

The four-cell result that field haemoglobin adds **+0.394** to Gambia women's
iron, reaching **0.848**, used **concurrent** haemoglobin measured on the same
respondents in the same survey. This uses **prior-round district-aggregated**
haemoglobin from a different survey and a different district system. They are
different predictors answering different questions, and the second is the
deployable one. The concurrent result stands; it is simply not available to a
country that has not yet run the survey.

## Reviewer pass, statistical

**Shared-noise check, as guardrail 9 requires.** For each arm, the fraction of
the evaluation target's sampling noise the prediction can see:

- **bridge**: zero. Every covariate comes from a different survey (prior DHS
  round) or from a raster. No respondent contributes to both the predictor and
  the target.
- **full 294**: zero for the same reason, with one caveat below.
- **flat regional mean (jackknife)**: zero by construction, since a district's
  regional estimate excludes its own respondents. The un-jackknifed variant is
  not used here.

The caveat on the 294: 97 of its columns are prior-round DHS aggregates, which
share the district system but not the respondents. Shared noise is zero;
shared *district-level* confounding is not, and that is a different concern
which the region-partialed family of WS-A addresses.

**Identical cells and folds.** All three arms run inside one loop per cell on
the same rows, the same Admin-1 blocking and the same scoring. The flat regional
arm reports slightly fewer units (12 against 14, 28 against 30) because a
district in a single-district region has no other district to average, and those
are dropped rather than filled.

**Seeds.** 20260921. The ridge is deterministic.

**Multiplicity.** Eight cells and three arms is 24 numbers; no correction is
applied and none is claimed. The Gambia women's iron result is reported as one
cell out of eight, not as a finding.

## Reviewer pass, reproducibility

**Targets graph.** No new target.

**Joins.** Survey to covariates and covariates to nutrition-proximal both go
through `admin2_join_by()`, so Malawi is keyed on the pair. Verified: no
duplicate `(country, Admin1, Admin2)` rows exist in
`nutrition_proximal.csv`, and no cell's unit count rose after the join.

**Stamps.** Reads `nutrition_proximal.csv` (WS-F2, provenance recorded),
`predictors_admin2_harmonized.csv` and the `_targets_full` store.

**Additive.** `iron_bridge.csv` is new; nothing is overwritten.

**Runtime.** About four minutes for 8 cells and three arms under both
evaluations.
