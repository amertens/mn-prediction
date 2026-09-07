# Two readings of the same results — 4 September 2026 (revision f)

Supersedes `TWO_READINGS_2026-09e.md`, which is preserved beside it. Revision e
folded in the sandbox work of 2–3 September. This revision adds the work of
4 September (`docs/findings/SANDBOX_LOG_2026-09.md`, scripts `protocol_v2/34`–`42`):
a variance-components reliability ceiling, an urbanicity test, an
anchor-and-rank survey design, fieldwork dates for every cluster, three add-on
feature blocks scored under one harness, a zinc collection test — and one bug
whose correction changes a headline number. Both abstracts describe the same
evidence; each is written so that a reader who checks the sources finds
nothing overstated.

**The correction first.** The population file spells Sierra Leone with a
space; every earlier result that joined population (burden capture, the
comparator-fairness test, the consistent-rung baseline, the admin-1 transport
and aggregation runs) silently dropped the country. The "12 of 12 regional
cells at 0.50–0.56" of revisions d and e was therefore four outcomes in three
countries. Re-run on all 22 cells the regional transport is **mean Spearman
0.32, 17 of 22 positive** (0.43 and 12 of 12 in the two countries with eight
or more regions). The district-level, in-fill and burden numbers were not
affected: the in-fill folds cannot be cut for Sierra Leone's 14 districts, so
its absence changed nothing there.

*Numbers refreshed twice on 4 September and once on 7 September. At 08:17, after the AlphaEarth
embedding was moved out of the agriculture domain and Gambia's food-price
window was corrected to its 2018 fieldwork, the regional full-index transport
moved from 0.29 to 0.31 and the district full-index transport from 0.25 to
0.26. At 23:25, after Gambia's missing IHME block was restored (15 of 25 empty
columns; sandbox log IH-01), the regional full-index transport settled at
0.30, the penalised district fit moved from 0.24 (21 of 22) to 0.27 (19 of 22)
and the regional permutation null's 95th percentile from 0.17 to 0.16. On 7
September three domains were added (livestock density, water and coast
proximity, helminth burden) and the IHME block was rebuilt from its 5 km
surfaces (sandbox log NL-01; 473 columns, 24 domains): the district
zero-tuning index moved from 0.26 to 0.27 (prevalence 0.16 to 0.18), the
penalised fit to 0.27 (20 of 22), the regional full index to 0.32 and its
null's 95th percentile to 0.15; the climate + soil, in-fill and burden
figures did not move.*

Every number comes from `results/tables/protocol_v2/` under one protocol:
replicated folds, three separately scored estimands, information-matched
baselines, precision-weighted scoring, a 451-predictor vocabulary reduced to
domain principal components, permutation nulls for transport, and — new in
this revision — a variance-components ceiling and a common add-on harness.

---

## Reading A: the optimistic abstract

**Remotely sensed climate and soil rank districts for micronutrient targeting
across borders, a national estimate is enough to put levels on the ranking,
and the model already sits at the reliability ceiling in a third of cells**

*Background.* Sub-national targeting of micronutrient programmes is limited by
the cost of biomarker surveys. Earlier negative findings — including our own —
rested on evaluation designs whose power and information symmetry had not been
quantified. Correcting them left a positive result whose mechanism, rival
explanations and practical use had not been tested. This revision tests them.

*Methods.* Across 24 country–outcome cells in four national biomarker surveys
(146 districts on a common administrative rung) we scored a 451-predictor
vocabulary, reduced to principal components within 18 domains, under three
estimands with information-matched baselines. Random folds were replicated;
districts were weighted by effective sample size. Transport counts were
calibrated against a country-block permutation null. The reliability ceiling
was re-estimated from a respondent-in-cluster-in-district-in-region variance
model. An urbanicity composite was partialled out of the transported index.
A survey design that measures only a national prevalence and lets the
transported covariates order the districts was costed against a district
survey of equal size. Fieldwork dates were recovered for every cluster, and
three new feature blocks — food prices and night-time temperature matched to
the fieldwork window, and a 64-dimension satellite embedding — were scored
under one harness.

*Results.* **Covariates beat the honest survey baseline within a country**
(0.398 vs 0.320 on the biomarker level, 0.286 vs 0.193 on prevalence, 15 of
18 cells). **Rankings transport to a country never seen in training**:
district-level mean Spearman 0.25–0.31 with the full vocabulary and 0.37 with
climate and soil alone (22 of 22 positive); regional-level 0.32 over all 22
cells and 0.45 with climate and soil (20 of 22), against a permutation null
whose 95th percentile is 0.15. **It is not an urban–rural map**: partialling
out night lights, population density, built surface and travel time moves
the climate + soil transport by less than 0.02, and urbanicity alone does not
transport (−0.07). **A national anchor is enough to put levels on the
ranking**: five percent of a survey's sample — fifty to seventy respondents —
plus the transported ranking gives a district error of 10 percentage points,
against 20 for a district survey of the same size and 12 for a regional one;
the district survey needs about 40 percent of the full sample, roughly 400
respondents, before it does better. **The model is at the honest ceiling in
five of 24 cells.** With cluster effects removed, the attainable correlation
averages 0.44 on prevalence; the model reaches or exceeds it in Gambia and
Ghana for child vitamin A and in Gambia for women's iron and vitamin A.
**Nothing we added helps, which is itself a result**: fieldwork-matched food
prices and night-time temperature, and the satellite embedding as a domain of
its own, each change in-fill and transport by less than 0.01. The vocabulary
is saturated; the two remotely sensed domains carry it, and the associations
they rest on replicate in sign across all four countries.

*Conclusion.* The deliverable is a **ranked district priority list** for any
country, and a **level-calibrated map for any country willing to measure one
national prevalence**. Seven pre-registered predictions for the next two
countries will test whether the two-domain index transports at least as well
as the full one; the current evidence says it transports better.

---

## Reading B: the pessimistic abstract

**Geography does the work; the covariate contribution is a ranking whose
regional strength was overstated by a join bug, the ceiling shows little room
left, and no new data source moves anything**

*Background.* Correcting an over-stated negative does not establish a
positive, and one of the positive numbers turned out to rest on a dropped
country. This reading states what survives when each result is scored against
the cheapest alternative, the fairest baseline, its own null, and the data it
was not selected on.

*Results.* **The regional headline was a subset.** "Twelve of twelve at
0.50–0.56" was three countries and four outcomes; on all 22 cells the regional
transport is 0.32, 17 positive, the same figure the null calibration had
reported all along. **Within a surveyed country geography is sufficient**: a
covariate-free smoother scores 0.391 to the full model's 0.398 and covariates
add nothing to its residuals. **The reliability ceiling was inflated by cluster
effects and the model has less headroom than claimed, not more**: with the
cluster removed the attainable correlation falls from 0.54 to 0.44 on
prevalence, one in five of it was cluster effect, six of 24 cells have a
ceiling below 0.30, and child zinc has no district geography at all — its
between-district and between-TA variance is zero, and what a split-half
called geography was the cluster. **The anchored design wins on level error
only.** It loses to a district survey of any size on burden captured (0.15
vs 0.29 in the worst-ranked fifth) because burden is concentrated in populous
districts, which a survey measures precisely, and it loses on ordering once
the survey exceeds a tenth of its sample. **The model does not substitute for
survey sample** (+3.36 vs +3.33 pp degradation from full to 15 percent). **The
climate + soil result was found on the data it is scored on**; honest nested
selection recovers +0.03. **No new data source adds anything**: food prices
matched to the fieldwork month, night-time temperature matched to the
fieldwork month, and a satellite embedding are all within ±0.01 of the base
vocabulary, and the embedding was already in the vocabulary — filed under
agriculture, where it made up 64 of 93 columns and so was two thirds of what
the domain ablation called "agriculture". **Fieldwork timing is a nuisance,
not a predictor**: the survey-month temperature columns are among the most
strongly and consistently associated with deficiency in the scan (four of four
countries for child vitamin A) yet improve no fold, because a covariate
defined by when the survey happened cannot rank districts that have not been
surveyed.

*Conclusion.* The honest product is a **ranked list, mean correlation about
0.3 at either tier, for countries without a survey**, plus a level calibration
that needs one national estimate. Within a surveyed country a smoother
delivers most of what the covariates do. The vocabulary is saturated, the
headroom is small, and every claim about levels, burden or survey savings
should be read as an upper bound or withdrawn.

---

## What moved since revision e

| Revision e claim | Status in f |
|:---|:---|
| Regional transport 0.50–0.56, 12 of 12 | **Corrected (BUG-01).** Three countries only. All 22 cells: 0.32 mean, 17 positive; climate + soil 0.45 (20 of 22); 0.43 and 12 of 12 in countries with ≥ 8 regions. |
| DA-03: climate + soil "a wash" at Admin-1 | **Reversed.** Same bug. On 22 cells climate + soil 0.45 vs full 0.32 on the level, 0.38 vs 0.33 on prevalence. |
| District transport 0.28–0.31 (full), 0.368 (climate + soil) | **Stands.** No population join involved. |
| Burden capture 24 / 19 / 20 percent; in-fill 15 of 18; CF-01 | **Stand.** Sierra Leone's 14 districts cannot be folded, so its absence changed nothing. |
| Reliability ceiling "an upper bound; 16–27% cluster effect on multi-cluster units" | **Quantified on all units (VC-01).** Honest ceiling 0.44 vs 0.54; 21% cluster share; model at ceiling in 5 of 24 cells; child zinc ceiling 0.27 / 0.00. |
| "Not an urbanicity map" was untested | **Tested (UR-01).** Partial correlations within 0.02 of raw; urbanicity alone −0.07. |
| Levels do not transport; own anchor halves error | **Made a design (AR-01).** 5% national sample + ranking: 10 pp vs 20 (district survey) vs 12 (regional); crossover at 40% of sample. Loses on burden at every size. |
| Zinc "cannot be distinguished from a collection artefact" | **Sharpened (ZN-02, VC-01).** No district variance exists; afternoon draw −3.5%, date +0.7%/day, neither explains district differences because there are none. |
| — | **New (FW-01):** fieldwork dates for all 328 clusters; Sierra Leone reconstructed from date of birth + age. |
| — | **New (AD-food, AD-climate, AD-alpha):** three add-on blocks, none earns its place. Script 07 used 2021 for Gambia's survey year (fieldwork was Jan–Apr 2018). |
| — | **New:** `CLUSTER_LEVEL_DESIGN_2026-09.md` — cluster linkage exists already; 328 clusters, a 1.6–2.2× gain in units, not "several hundred per country". |

## Paper outline

Unchanged in structure from revision e. Section 4 (prediction) replaces the
regional headline with the 22-cell figure and adds the anchor-and-rank
design; Section 5 (protocol) replaces the split-half ceiling with the
variance-components one and records BUG-01 as a worked example of why a
dropped country has to be impossible, not merely unlikely.
