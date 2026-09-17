# Two readings of the same results — 17 September 2026 (revision g)

Revision g rewrites both abstracts on the RR-12 results (575-column set, headline
tiers open + public survey microdata, no DHS) and the analyses added between
15 and 17 September: the survey-report reconciliation (iron and zinc binaries,
RBP rule), the exact pair count (PC-01), the per-cell performance tables
(CB-01), the five-domain transport candidate (DA-03), the extended Cote
d'Ivoire database and rank uncertainty for both candidates (CV-01), and the
calibration check on the resampling intervals (VZ-01). Every number below is
read from results/tables/protocol_v2 or results/tables/policy_deck through the
full-talk deck's setup chunk (docs/slides/MN-proxy-full-talk-2026-09.qmd), so
the two abstracts and the deck cannot disagree. Revision f (4 September) is
superseded; its "what moved" section is at the foot of this file.

---

## Reading A: the optimistic abstract

**Public environmental data rank districts by micronutrient deficiency inside
and across countries; where the survey resolves districts the model reaches
0.64, a national estimate puts a level on the ranking, and the product has
been delivered for a country with no survey at all**

*Background.* Sub-national targeting of micronutrient programmes is limited
by the cost and rarity of biomarker surveys: about one per country per decade,
resolving regions, not districts. We asked whether public data alone — climate,
soil, satellite imagery, modelled health surfaces, crop and livestock grids,
public household-survey microdata — can order a country's districts by
deficiency, inside a surveyed country and in one never surveyed.

*Methods.* Four national biomarker surveys (The Gambia 2018, Ghana 2017,
Malawi 2015-16, Sierra Leone 2013; 206 surveyed districts, 323 GPS clusters)
were reconciled to one outcome definition per nutrient against each survey
report (iron, vitamin A, folate, B12, zinc; 24 country-outcome cells) and linked
to 575 district-level predictors from 28 open or public sources in 29 domains,
of which 383 enter the headline models. Predictors are rank-normalised within
country and reduced to domain principal components learned on training
districts only; the estimator is a zero-tuning index that weights each
component by its training-fold rank correlation, so nothing is tuned on
held-out data. Three estimands with information-matched baselines: in-fill
(5-fold by district, ten draws), region extrapolation, and leave-one-country-out
transport scored against a country-block permutation null. Pairwise ordering
accuracy is counted exactly; burden capture is population-weighted; the
reliability ceiling comes from a respondent-in-cluster-in-district variance
model; rank uncertainty from 400 stratified resamples of the training set,
and its calibration from the same procedure with each surveyed country held
out.

*Results.* **Inside a surveyed country the index ranks districts at Spearman
0.40 on the biomarker level** (prevalence 0.28) against 0.31 for the survey's
own jackknifed regional average (better in 13 of 18 cells) and 0.39 for a
covariate-free spatial smoother; it orders 60% of district pairs the way the
survey does, against 56% for the regional average and 50% by chance. **The
average hides where it works: six of eighteen in-country cells clear 0.5 and
average 0.64** — Malawi women's B12 (0.70), The Gambia women's vitamin A
(0.70), child vitamin A (0.69) and women's iron (0.68), Ghana women's B12
(0.57) and child iron (0.53) — and in those cells the model sits at 85-90% of
its own attainable ceiling (0.62-0.78). **Accuracy follows the survey's
resolution more than the nutrient**: The Gambia, whose 30 districts each hold
several clusters, scores 0.41-0.70 on all four outcomes; the single-cluster
districts of Ghana and Malawi average 0.34. **Rankings transport to a country
never seen in training**: 0.30 at district level with the full vocabulary
(18 of 22 cells positive; null 0.08), 0.38 with climate and soil alone
(21 of 22) and 0.45 at regional level; a five-domain set adding the anaemia
surfaces, agriculture and infection reaches 0.41 (better than climate and soil
in 15 of 22, all 22 positive) and enters the pre-registration as a second
candidate. It is not an urban-rural map: partialling out an urbanicity
composite moves the transported index by 0.01, and urbanicity alone does not
transport (-0.07). **Levels follow from one national number**: a 5% national
sample plus the transported ranking gives a median district error of 9.3
percentage points, which a district-representative survey matches only at 25%
of the full sample. **Uncertainty is on every map and has been checked**: 90%
rank intervals from resampling the training set are 5 places wide out of 33
for Cote d'Ivoire, and their coverage under leave-one-country-out (38%)
identifies exactly the transport error the correlation reports, so they are
reported as stability, not as truth. **The product exists for a country with
no survey**: Cote d'Ivoire's 33 districts are ranked for six outcomes on a
212-column open database; the two pre-registered candidates agree at Spearman
0.88-0.97 with 83% of the worst fifth in common, the northern savanna belt is
the answer under either, and the one external map (WFP/MIMI vitamin A) ranks
the same belt worst. **The vocabulary is saturated and its signal is
replicated**: climate and soil are the only load-bearing domains across
borders, the environmental domains carry 39-54% of the index, eleven added
blocks each move it by less than 0.01, and the district-level associations
the index rests on (drier, grassier, poorer, more pastoral districts rank
worse) hold in sign in every country.

*Conclusion.* The deliverable is a ranked district priority list for any
country from public data, a level-calibrated planning map for any country
willing to measure one national prevalence, a live dashboard carrying both
with their uncertainty, and a pre-registered pair of transport indices for the
fifth survey. The evidence supports using the ranking to decide which
districts to reach first and where the next survey's clusters go; the case for
each survey added is that it improves every other country's map (0.20 with
one training survey, 0.30 with three, not yet flattened).

---

## Reading B: the pessimistic abstract

**A ranking of about 0.3 in a new country, no better than the survey's
regional averages on burden or level inside a surveyed one, with intervals
that measure stability rather than truth; the strong cells are where the
survey is good, not where the model is**

*Background.* The corrected protocol removed the leakage and fold artefacts
of the January analyses and left a positive result. This reading scores each
surviving claim against the cheapest alternative, its own null, the survey's
noise, and the data it was not selected on.

*Results.* **Inside a surveyed country the covariates add almost nothing to
geography**: the index scores 0.40, a covariate-free smoother 0.39, and
covariates fitted to the smoother's residuals 0.39. **The gain over the
survey's own regional average is real but small**: 0.40 vs 0.31 on the
ranking, 60% vs 56% of district pairs, and on the two things a programme
acts on it disappears — the worst-ranked fifth holds 22% of the deficient
people against 20% for the regional average (a perfect map 48%), and the
district prevalence error of the index is no better than the regional
average's even in its four best cells (10.6 vs 8.8 points for Malawi B12,
10.3 vs 8.1 for Gambian child vitamin A), because a rank-rescaled score is
not a calibrated level. On WHO vitamin A severity bands the regional average
is slightly better (66% vs 64% exact, 94% vs 91% within one band). **Two
thirds of the in-country cells fail**: twelve of eighteen average 0.28,
Malawi's vitamin A cells at under 1% prevalence cannot be ranked by anything,
zinc has no district signal (-0.07), and the six strong cells are the cells
where the survey resolves districts (The Gambia's multi-cluster districts),
so "0.64" describes the survey design as much as the model. **The model does
not work uniformly across the settlement gradient**: within the most rural
third of districts the held-out ranking is 0.31 in Ghana and -0.35 in Malawi,
where deficiency is highest, so in Malawi the index partly separates
lakeshore towns from uplands rather than ordering rural districts. **What it
learns can be a surrogate**: in its best cell the modelled anaemia and malaria
surfaces carry negative weight for B12 because they mark the fish-eating
lakeshore, an ecological association with no mechanism and no guarantee of
transport. **Transport is 0.30 with the vocabulary the protocol specifies**,
positive in 18 of 22 cells against a null of 0.08, 0.28 at regional level;
the 0.38 climate-and-soil figure and the 0.41 five-domain figure were both
chosen on the 22 cells they are scored on, and honest nested selection inside
the training fold gives 0.31, no better than the full index. Person-level
classifiers on the same proxies are a coin toss (AUC 0.50). **The intervals
on the maps are not intervals for the truth**: the 90% rank intervals from
resampling the training set cover the survey's rank 16-70% of the time when
a surveyed country is held out, mean 38%; they describe how the ranking
would move under a different training draw and nothing more, and no
calibrated interval yet exists. **Cote d'Ivoire is untested**: the two
candidates agree with each other (0.88-0.97) and give nearly the same map for
every outcome (0.93 between outcomes), which says the transported product is
a map of environmental deprivation rather than of any nutrient, and nothing
inside the country can say whether it is right. **The anchored survey design
wins on level error only**: at 5% of a sample it loses to a district survey
of any size on burden captured, and the model does not substitute for sample.
**Nothing added helps**: eleven candidate blocks, three alternative weightings,
sparse composites and a cluster-level fit are all within 0.01-0.02 of the
index, and the one thing that would raise the ceiling — more clusters per
district — is a survey decision, not a modelling one.

*Conclusion.* The honest product is a ranked list, mean Spearman about 0.3 in
a country without a survey and 0.4 inside one, whose value lies only where
the survey's regional averages do not exist; a level requires a national
estimate and remains a planning figure; the uncertainty shown is stability,
not calibration. Claims about burden captured, survey savings or nutrient-
specific maps in unsurveyed countries should be read as upper bounds or
withdrawn until a fifth survey scores the pre-registered candidates.

---

## What moved since revision f (4 September)

- Predictor set 451 -> 575 columns (RR-11/12), 383 in the headline tiers, no
  DHS; outcomes re-reconciled to the survey reports (iron/zinc binaries were
  the IDA columns, RBP < 0.70 rule for vitamin A).
- In-fill 0.398 -> 0.399 (level), transport 0.31 -> 0.30 full, 0.37 -> 0.38
  climate + soil; regional transport 0.31 -> 0.28 full, 0.46 -> 0.45 c+s.
- New: exact pair count (60 / 56 / 55%); per-cell tables and the strong-cell
  finding (6 of 18, mean 0.64; Gambia 0.41-0.70); five-domain candidate 0.41
  (DA-03); Cote d'Ivoire on a 212-column database with both candidates and
  rank uncertainty for six outcomes; coverage of the resampling intervals 38%
  (VZ-01); urbanicity strata (Malawi rural third -0.35); surrogate-marker
  reading of the Malawi B12 weights; anchored design 10 -> 9.3 points at 5%,
  matched by a district survey at 25% (was 40%).
- Withdrawn from the optimistic reading: "at the ceiling in three of 24
  cells" (the per-cell ceiling now reads 0.62-0.78 in the strong cells and the
  model is at 85-90% of it, not above it); "nothing we added helps" is kept
  but the five-domain set is now the one exception, labelled post hoc.
