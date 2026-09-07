# Sandbox log — methodological experiments, September 2026

Running record of experiments run OUTSIDE the `{targets}` pipeline, in
`scripts/protocol_v2/`, with results in `results/tables/protocol_v2/`. Each
entry states the question, the design, the result, and what (if anything) it
changes in the project's claims. Nothing here modifies `_targets.R`. The one
production change made on request — `fit_predict_sl_prescreened()` now using
the weighted, blocked helper in `R/area_superlearner.R` — is recorded in
[PROTOCOL_V2.md](PROTOCOL_V2.md) and in the memory note
`superlearner_package_choice.md`.

Companion documents: [PROTOCOL_V2.md](PROTOCOL_V2.md) (the corrected
evaluation layer), [TWO_READINGS_2026-09d.md](TWO_READINGS_2026-09d.md)
(claim status), [SLIDE_OUTLINE_2026-09.md](SLIDE_OUTLINE_2026-09.md).

---

## Standing conclusions this log has reached so far

1. **Evaluation design moved the answers more than any model.** Replicated
   folds, information-symmetric baselines, and power checks reversed several
   published negatives (see PROTOCOL_V2.md).
2. **The signal is regional, and it is agro-ecological.** Transport of
   rankings to an unseen country: 12/12 at the first sub-national tier (mean
   Spearman ~0.50), 19/22 at the district rung (0.31). Within a surveyed
   country, geography alone captures nearly all district-level accuracy;
   covariates add nothing on top of a spatial smoother (nested test, 4
   estimand/target combos, all null). The replicated covariate associations
   run *opposite* to individual-level nutrition for legumes and cattle — a
   rural-subsistence axis, not diet.
3. **Capacity is a liability at n = 14–87.** A zero-tuning domain index beats
   every tuned learner and a SuperLearner that contains it. The SuperLearner's
   MSE meta-learner selects shrinkage (a constant learner gets the largest
   weight); a rank-aligned meta-learner recovers the index's performance but
   does not exceed it.
4. **Each added training country buys ~0.05 of transported rank accuracy**,
   monotone across all four specifications tested (1→2→3 countries).
5. **Survey size buys regional precision, not district resolution — and the
   model does not substitute for sample.** District MAE from survey alone is
   flat in sample size (10.6 → 9.2 pp from 5% to 100%); regional MAE falls
   steeply. Fitted symmetrically on the same reduced survey, the model's
   accuracy degrades in lockstep with the survey's (+3.4 vs +3.3 pp from 100%
   to 15% of sample); what it adds is a roughly constant ~0.6 pp and a much
   better district ordering at any given sample size (G4-02).
6. **Transport is carried by remotely sensed climate and soil.** Dropping any
   other domain costs nothing or helps; a fixed climate+soil index transports
   better than the full index (0.368 vs 0.252 on level, positive in 22/22) —
   but it was chosen on these four countries, so it is a pre-registered
   prediction for the next ones, not a result. Honest nested selection gains
   only +0.03.
7. **The model's advantage over the survey's own regional averages is in
   ranking (large, robust), not in burden captured (small, inconsistent).**
   Part of the "regional averages ≈ chance" result is the jackknife
   over-correcting at 3–5 districts per region.
8. **Most "district" estimates are single survey clusters** (Gambia 57%,
   Ghana 83%, Malawi 85% of units; Sierra Leone 0%). The split-half ceiling
   therefore counts cluster-level effects as geography wherever units are
   single-cluster: on multi-cluster units 16–27% of the within ceiling is
   cluster effect, and at Malawi's district rung only child iron retains
   reliable geography once clusters are split. Malawi's outsized headroom was
   mostly artefact; quoted headroom is an upper bound.

---

## Entries

### 2026-09-03 · SL-01 — domain index inside a SuperLearner (script 19)
**Q.** Does a SuperLearner select the domain index if it is in the library?
**Design.** 24 cells × 5 reps × 5 district folds; library {mean, enet, rf,
domain_index}; survey weights; classic `SuperLearner` (mlr3superlearner cannot
take a custom learner, has no weights, and ignores `group=` for regression).
**Result.** Discrete SL picks the index in **1%** of 480 fits; NNLS weight
0.22 on the index, **0.42 on the constant**. OOF Spearman: index 0.285 > rf
0.247 > SL-NNLS 0.231 > SL-discrete 0.173. SL beats the index it contains in
2 of 18 cells. Weighted/blocked SuperLearner beats mlr3superlearner discrete
in 11/18 (+0.037).
**Why.** CV squared error at this n is minimised by shrinking toward the mean;
the decision metric is a ranking. Loss and decision disagree.
**Changes.** Nothing in the NCE (it cites the index directly). Paper: the
parsimony claim strengthens — the index beats an ensemble that includes it.

### 2026-09-03 · SL-02 — rank-aligned meta-learner (script 20)
**Q.** Does aligning the meta-learner loss with ranking fix SL-01?
**Design.** Same base fits and folds; four meta-learners derived from the same
`Z`: MSE-discrete, MSE-NNLS, rank-discrete (argmax CV Spearman), rank-NNLS
(NNLS on ranks). Production object `method.asl_rank` validated: reproduces the
derived coefficients exactly.
**Result.** Index selected in **64%** of fits (was 1%); ensemble weight on the
index 0.50 (constant falls to 0.19). Spearman: index 0.285 ≥ rank-NNLS 0.269 ≈
rank-discrete 0.266 > rf 0.247 > MSE-NNLS 0.234 > MSE-discrete 0.181.
Top-quintile overlap: rank-NNLS 0.412, index 0.410, MSE-discrete 0.347.
Rank beats MSE in 16/18 (discrete) and 14/18 (NNLS) cells; rank beats the
index alone in only 4–6 of 18 (median −0.015).
**Changes.** Helper exposes `meta = "rank"`; default left at `"mse"` pending
SL-03 and the production-metric check (SL-04).

### 2026-09-03 · SL-03 — population-aware meta-learners (script 21)
**Q.** Can a meta-learner that knows district population beat the index?
**Design.** As SL-02 plus `wrank` (population-weighted Spearman / weighted
NNLS on ranks) and `burden` (share of prevalence × population captured in the
top fifth; Nelder-Mead over softmax weights with restarts). Evaluated on
Spearman, top-quintile overlap, burden capture with test-row population, and
MAE. Factories `make_method_asl_wrank()` / `make_method_asl_burden()` validated
against derived coefficients.
**Result.** Factories reproduce the derived coefficients exactly. Index
selected by the discrete rule: mse 0.7%, rank 64.7%, wrank 60.9%, burden
32.7%. Out of fold (cell medians; capture random = 0.20):

| arm | Spearman | top-quintile | capture | MAE pp |
|---|---|---|---|---|
| **domain_index alone** | **0.285** | 0.410 | **0.242** | 13.97 |
| burden_discrete | 0.204 | 0.379 | 0.238 | 12.48 |
| burden_ens | 0.207 | 0.370 | 0.237 | 11.89 |
| rank_nnls | 0.269 | **0.412** | 0.228 | 11.84 |
| wrank_nnls | 0.259 | 0.401 | 0.219 | 11.67 |
| mse_discrete (old production) | 0.181 | 0.347 | 0.195 | 11.22 |

burden_ens beats the index on capture in 13 of 18 cells but by a hair
(median +0.010) and loses Spearman in 15 of 18 (median −0.058): optimising
the top fifth directly buys a marginal, inconsistent capture gain at a large
cost in overall ordering, and the mean capture still favours the index.
wrank does not improve on rank. **No meta-learner beats the index alone.**
**Verdict (with SL-04).** The rank-family losses fix the SuperLearner
(vs MSE: +0.08 Spearman, 15 of 18 cells) but cannot beat the best single
ranker at n = 14–87, and on the production benchmark — no index in the
library, level metrics — rank does not win. **Default stays `meta = "mse"`;
`rank`, `wrank`, `burden` remain available.** The paper sentence: four
meta-learner losses tested; an ensemble at best matches the simple index.

### 2026-09-03 · QA-01 — sign-convention regression test (script 24) — PASS
12 assertions against the p1/p4 scan outputs and the dashboard bundle
(legumes/child vitamin A +, cattle/child iron +, root crops/child iron −,
wasting/women's vitamin A +; p4 negation present in code and header). Exit
status is non-zero on failure so it can gate a rebuild.

### 2026-09-03 · CF-01 — comparator fairness (script 22)
**Q.** How much of "regional averages are no better than chance" is the
jackknife over-correcting at 3–5 surveyed districts per region?
**Design.** 18 in-fill cells; four regional-mean comparators on identical
data — `full` (includes the district's own respondents; leaks), `jk` (the
protocol-v2 baseline), `jk_shrunk` (jackknife shrunk toward the national mean
by k/(k+1), k = other districts in the region), `split` (mean of a random half
of the other districts, 20 splits) — against the model's saved cell results.
**Result (mean over cells; burden capture, random = 0.20):**

| arm | Spearman | capture |
|---|---|---|
| full (with self) | 0.562 | **0.300** |
| model (domain_index) | 0.286 | 0.240 |
| split | 0.199 | 0.229 |
| jk | 0.190 | 0.213 |
| jk_shrunk | 0.204 | 0.210 |

Per cell on capture: model beats jk in **8 of 18** (median −0.017), beats
split in 8/18, beats full in 6/18. Shrinkage does not repair the jackknife
(jk_shrunk vs jk: 5/18, median 0.000) — shrinking toward the national mean
removes ranking information rather than restoring the district's own signal.
By region size: with ≤3.5 districts per region the model captures 0.272 vs jk
0.213; with 3.5–5 it is 0.214 vs 0.213 — a tie.
**Reading.** The model's ranking advantage over the survey baseline is robust
(Spearman 0.286 vs 0.19–0.20; 15 of 18 cells). Its BURDEN-capture advantage
is not: pooled 0.240 vs 0.213–0.229 for any comparator that never sees the
district, and no better than a coin flip per cell, with the pooled margin
concentrated in the countries whose regions are tiniest — exactly where the
jackknife penalises the baseline most. A regional average that includes the
district (what a surveyed district's region average really is) captures 0.300
and beats the model.
**Changes.** The NCE sentence "24% vs 19%" is correct as computed but should
not be read as a robust programmatic margin; the defensible claims are (i)
ranking, 15 of 18, and (ii) the settings where no survey baseline exists at
all (transport). Recommend the burden sentence be scoped to "districts a survey
has not reached" and not repeated as a headline. Recorded for the NCE thread.

### 2026-09-03 · SL-04 — rank meta-learner on the production benchmark — NO WIN
`fit_predict_sl_prescreened()` now reads `ASL_META` (default "mse"); the switch
was verified to change the pick and weights on a fixed split (women's vitamin
A, Malawi held out: mse picks lasso, rank picks enet). Full production driver
rerun with `ASL_META=rank`, 22 LOCO holdouts, same folds:

| metric | mse | rank | rank better in |
|---|---|---|---|
| Pearson r | 0.034 | 0.059 | 4 of 20 |
| Spearman | 0.039 | 0.044 | 5 of 20 |
| MAE (pp) | 14.96 | 16.14 | 2 of 22 |
| |bias| (pp) | 1.02 | 3.05 | 5 of 22 |

Median change 0.000 on every metric: in most holdouts both losses pick the
same learner (the production library has no domain index, and a tree
ensemble dominates both criteria). Where they differ, rank is worse on level.
**Verdict: does not win on the production benchmark; default stays `mse`.**
The sandbox result (SL-02) stands for the ranking use case with the index in
the library; the two are not in conflict — rank loss helps when there is a
strong ranker to select, and there is none in the production library.

### 2026-09-03 · DA-01 — domain ablation under LOCO (script 23)
**Q.** Which of the 18 domains present carry cross-country transport?
**Design.** Zero-tuning domain index on identical LOCO folds; drop-one and
only-one ablations; 22 cells per target. First launch crashed because the
number of PCs per domain is chosen on training rows and differs by fold;
fixed by intersecting columns per fold.
**Result (mean LOCO Spearman).** Full index: level 0.252, prevalence 0.151.

| domain | Δ when dropped (level / prev) | alone (level / prev) |
|---|---|---|
| **Soil characteristics** | +0.041 / +0.044 | **0.318 / 0.249** |
| **Climate and weather** | +0.034 / +0.024 | **0.351 / 0.214** |
| Agricultural production, land use | +0.014 / +0.012 | 0.216 / 0.102 |
| Malaria incidence and treatment | +0.010 / +0.006 | 0.129 / 0.078 |
| all other 14 domains | ≤ +0.001, ten are **negative** | ≤ 0.19 |

Removing household assets (−0.018), built environment (−0.014), healthcare
access, fertility or child morbidity (−0.010 each) *improves* transport.
**Reading.** Climate alone and soil alone each transport better than the
full index. The survey-derived socioeconomic, health, WASH and dietary
domains are dead weight or harmful across borders — consistent with their
meaning and scale being country-specific even after within-country rank
normalisation, while remotely sensed soil and climate mean the same thing
everywhere. This is the agro-ecology axis again, now as a transport result.
**Caveat.** "Only" sets were identified on all-country LOCO, so choosing them
and then reporting their LOCO score is optimistic. DA-02 (queued) selects
domains inside the training countries only.
**Changes.** "What to collect" for transport: the remotely sensed layers.
Candidate simplification of the production index to climate + soil (+
agriculture), pending DA-02.


### 2026-09-03 · DA-01 — domain ablation under LOCO (script 23) — RUNNING
First launch crashed ("subscript out of bounds"): the number of PCs per domain
is chosen to 80% variance on the TRAINING rows, so D's column set differs by
held-out country and a column list from one fold does not index another.
Fixed by intersecting per fold. Relaunched.

### 2026-09-03 · DA-02 — nested domain selection under LOCO (script 25)
**Q.** Does a few-domain index transport better when the domains are chosen
honestly (inside the training countries only)?
**Design.** For each held-out country, inner LOCO over the three training
countries scores domains and greedily adds them (cap 5); the chosen set is fit
on all three and scored once on the held-out country. Comparators: the full
18-domain index; fixed climate+soil and climate+soil+agriculture (pre-specified
from DA-01, therefore optimistic).
**Result (LOCO Spearman, 22 cells per target):**

| set | level | cells + | prev | cells + |
|---|---|---|---|---|
| **fixed climate+soil** (2 domains) | **0.368** | **22/22** | **0.268** | 21/22 |
| fixed climate+soil+agri (3) | 0.349 | 21/22 | 0.246 | 20/22 |
| nested (honest, mean 3.3 domains) | 0.280 | 19/22 | 0.095 | 15/22 |
| full (18 domains) | 0.252 | 17/22 | 0.151 | 15/22 |

Nested vs full: level better in 13/22 (+0.028), prevalence worse (8/22,
−0.040). Domains chosen by inner LOCO are scattered (soil 57%, climate 39%,
malaria 34%, greenness 23%, WASH 23%…): with three training countries the
inner selection is too noisy to reliably find the two-domain set.
**Reading.** The climate+soil index is a strong, stable hypothesis — positive
transport in every one of 22 combinations on the biomarker level — but it was
identified on these four countries, so 0.368 is an upper bound, not a
validated number. Honest selection recovers a modest gain on level and none
on prevalence. **The right use is as a pre-registered prediction for the next
countries (Ethiopia, Pakistan): a two-domain remotely-sensed index should
transport rankings at least as well as the full index.**
**Changes.** Adds a concrete, testable claim for the NCE's "add countries"
activity; nothing in the current claims changes.

### 2026-09-03 · HR-01 — where the headroom is
Achieved in-fill Spearman (domain_index, biomarker level) against the
empirical split-half ceiling ("within" scheme; `headroom_by_cell.csv`, 21
cells — Sierra Leone has no empirical ceiling in that scheme). Mean gap by
country: **Gambia 0.12** (achieved 0.63 vs ceiling 0.76), **Ghana 0.11**
(0.46 vs 0.56), **Malawi 0.53** (achieved 0.12 vs ceiling 0.65). Malawi — the
Traditional-Authority rung, median 9 biomarker measurements per unit — is
where nearly all remaining headroom sits, and its ceiling is *not* low. The
zinc outcomes are the extreme: ceilings ~0.8, achieved negative — a reliable
target the proxies do not touch at all. Women's B12 sits above its ceiling,
which means that ceiling is uncertain, not that the model is superhuman.
**Reading.** "Two-thirds of the attainable ceiling" is a pooled average of
one country at the ceiling and one far below it. The NCE headroom sentence is
still fair, but the optimism it implies is concentrated in Malawi (and in
zinc), not spread evenly.

### 2026-09-03 · R6-02 — the survey baseline on the consistent district rung (script 27)
**Q.** Run 6 put all four countries on the same rung but could not compute
the regional survey baseline for Malawi (GADM has no region tier for it). With
a verified district→region lookup (27/27 mapped; Northern 5, Central 9,
Southern 13), what does the model-vs-survey comparison look like on 146
genuinely comparable units?
**Result (in-fill, prevalence, 18 cells, 10 replicates):**

| arm | Spearman | burden capture |
|---|---|---|
| domain_index | **0.311** | 0.204 |
| spatial_plus_domain | 0.300 | 0.234 |
| spatial | 0.297 | **0.237** |
| region_mean_jk (survey baseline) | 0.140 | 0.188 |
| national mean | — | 0.140 |

Model vs survey baseline: Spearman better in 12/18 (index, median +0.141) and
15/18 (spatial+domain, +0.106); burden capture better in 11/18 (+0.015) and
12/18 (+0.033).
**Reading.** Same shape as the mixed rung and CF-01: the ranking advantage
over the survey's own regional averages is large and robust; the burden
advantage is small. Two new details. With Malawi on its three real regions
the regional baseline ranks worse than on the mixed rung (0.140 vs 0.193),
because three regional means carry almost no ordering information over 27
districts — so on this rung the baseline is weaker for a *structural* reason,
not a jackknife artefact. And the spatial arms capture more burden than the
index here (0.237 vs 0.204) while ranking slightly worse — burden weights
populous districts, and the smoother places them better.
**Changes.** None to claims. The consistent-rung numbers can now be quoted
in full (ranking 0.31 vs 0.14; capture 0.20–0.24 vs 0.19) if the NCE moves to
the same-rung framing.

### 2026-09-03 · DA-03 — the climate+soil index at the first sub-national tier (script 28)
**Q.** Does the district-rung finding (climate+soil beats the full index) hold
one tier up, where the headline 12/12 result lives?
**Result (Admin-1 LOCO, 12 cells per target; all sets 12/12 positive):**

| set | level mean / median | prev mean / median | capture (prev) |
|---|---|---|---|
| full (18 domains) | **0.556** / 0.541 | **0.498** / 0.518 | 0.188 |
| climate+soil | 0.508 / **0.587** | 0.490 / **0.600** | **0.233** |
| climate+soil+agri | 0.527 / 0.556 | 0.488 / 0.596 | 0.208 |

Per cell climate+soil beats full in only 4–5 of 12 (median −0.04 / −0.02).
Excluding the two thin countries (n < 8 units) the ordering is the same.
**Reading.** At the regional tier the simplification is a wash: means favour
the full index, medians and burden capture favour the two-domain set, and
every set transports in every cell. The dead-weight domains hurt at the
district rung, where units are small and noisy, and stop mattering once
units are large enough. So the pre-registered prediction for new countries
should read "a climate+soil index transports **at least as well** as the
full index," not "better."

### 2026-09-03 · ZN-01 — zinc: a reliable target the proxies do not touch (script 29)
**Q.** Why are the zinc outcomes the only ones with a high empirical ceiling
(~0.80) and no in-fill skill at all? Malawi only (87 Traditional
Authorities), so this is a diagnosis, not a replicated claim.
**Result.**
1. The target is coherent: child zinc and women's zinc prevalence, two
   independent samples of the same places, correlate at **ρ = 0.55**; ceilings
   0.80 (0.73–0.86) and 0.79 (0.68–0.87).
2. It is orthogonal to the axis the model learns: |ρ| ≤ 0.30 with every other
   Malawi outcome, and *negative* with iron and folate.
3. No indicator survives a family-wise permutation test for either zinc
   outcome (top |r| 0.30–0.33, p_fwer ≥ 0.38, 428 predictors); the same scan
   finds one survivor for child iron as a positive control.
4. Soil zinc columns: r 0.09–0.15 and in the direction of *more* deficiency
   with more soil zinc — the wrong sign and negligible.
**Reading.** Zinc has a real, shared spatial pattern that none of the 451
proxies track — a genuine proxy gap rather than a noisy target. Two
explanations are worth separating before it is called a biological gap.
(a) Zinc status may genuinely follow a pathway the proxies miss (dietary
phytate:zinc, crop varieties). (b) **Serum zinc is unusually sensitive to
collection conditions** — time since last meal, time of day, haemolysis,
tube type — and those vary by field team and schedule, which are themselves
spatially clustered. A TA-level pattern driven by *who collected the sample
and when* would look exactly like this: highly reliable across two
populations drawn by the same teams on the same days, and unrelated to
anything on the ground. ZN-02 (queued) checks that.

### 2026-09-03 · G4-02 — does the model let a country field a smaller survey? (script 26)
**Q.** Script 17's model arm had seen the full survey while the survey arm
saw a subsample. Here BOTH see only the smaller survey (incremental sampling
noise added to the district estimates; the model refitted on the noisy
outcome, 5-fold out of fold), both scored against the full survey.
**Result (median district MAE, pp; 18 cells × 8 replicates):**

| share of survey | survey (regional mean) | spatial + covariates | index alone | survey level + index pattern |
|---|---|---|---|---|
| 100% | 9.26 | 9.38 | 10.66 | 12.07 |
| 60% | 9.86 | 11.24 | 13.29 | 13.67 |
| 25% | 11.31 | 11.20 | 13.20 | 17.70 |
| 15% | 12.59 | 12.74 | 17.32 | 20.87 |
| degradation 100% → 15% | **+3.33** | **+3.36** | +6.66 | +8.80 |

Head-to-head on identical draws: spatial+covariates has lower MAE than the
survey arm in 99 of 144 cell-fractions (median −0.6 pp); the index alone in
49 of 144 (+2.1 pp); the anchored pattern in 5 of 108. On ranking the model
is clearly better at every fraction ≥ 25% (Spearman ~0.25 vs ~0.10–0.14) and
collapses with the survey at 15%.
**Reading.** **The model does not substitute for survey sample.** Its
accuracy falls in lockstep with the survey's as the sample shrinks (the
spatial + covariate arm) or faster (the index). What it adds is roughly
constant at any sample size: ~0.6 pp lower district error and a much better
district *ordering* than the regional mean. So the honest sentence is "at a
given survey size the model gives better district rankings," not "the model
lets you field a smaller survey."
**Changes.** The NCE sentence on survey size was already softened to a
question; the evidence now answers it in the negative. Recommend the draft
say the extension will quantify what the model adds *at* a given survey size
and the survey-size economics of regional vs district estimates — and drop
the cost-savings-from-a-smaller-sample framing.

### 2026-09-03 · ZN-02 (lite) and CS-01 — single-cluster units and what the ceiling measures
**Q.** Can the zinc "reliability" be a collection artefact? The store has no
interview date/team/time columns, so the direct test needs raw MNS field
files. Proxy: the reliability file's two split schemes — respondents split
at random *within* each unit vs whole *clusters* split — should diverge if
unit-level agreement comes from shared collection conditions.
**Result.** For Malawi the cluster-split ceiling is 0.000 for child zinc and
child vitamin A, but that is degeneracy, not evidence: **74 of Malawi's 87
Traditional Authorities contain exactly one survey cluster.** The check
generalises and is the real finding:

| country | units | clusters per unit (median, max) | single-cluster units |
|---|---|---|---|
| Gambia | 30 | 1 (max 13) | 17 (57%) |
| Ghana | 75 | 1 (max 3) | 62 (83%) |
| Malawi | 87 | 1 (max 3) | 74 (85%) |
| Sierra Leone | 14 | 4 (max 8) | 0 (0%) |

**Reading.** In three of four countries a "district prevalence" is, for most
districts, one enumeration area's dozen respondents. Two consequences.
(1) The empirical split-half ceiling — computed by splitting respondents
*within* a unit — then measures within-cluster agreement, which includes
everything a cluster shares (village, team, day, assay batch) as if it were
district geography. **The ceiling, and therefore the headroom, is inflated by
cluster-level effects wherever units are single-cluster**, which is most of
Ghana and Malawi. The "two-thirds of the attainable ceiling" sentence is
correct as computed but the remaining third is partly cluster noise no proxy
can or should predict. (2) Zinc's high ceiling with no proxy signal is
exactly what a cluster-level effect would look like; it cannot be separated
from real geography at the TA rung. Women's zinc holds a 0.71 ceiling under
the (partially degenerate) cluster scheme, so real geography is likely
present too. Verdict: undetermined at this rung; needs the field metadata.
**Changes.** Headroom should be re-estimated on multi-cluster units (CE-01,
queued). Reinforces Run 6: Malawi belongs at its district rung (27 units,
~3 clusters each), not TAs. Adds to the memory note on survey resolution.

### 2026-09-03 · LV-01 — transporting levels via anchors (already established in the corrected layer)
`R/corrected/p13_anchored_transport.R` / `results/tables/corrected/anchored_transport_loco.csv`
(35 holdouts): anchoring a transported map to the held-out country's **own**
national survey value cuts MAE from 13.7 to 8.9 pp; anchoring to an
**external** published estimate makes it far worse (44.7 pp) because of year
and definition mismatch. Rankings are untouched by construction. So the level
problem is solvable only with a number the target country itself supplies —
the argument for a small in-country survey — and not with published
external anchors as they stand. Nothing to re-run.

### 2026-09-03 · TC-02 — training-country curve for the climate+soil index (script 30)
**Q.** Does a two-domain remotely-sensed index also improve as training
countries are added, or is it saturated at three? Same subsets as script 15.
**Result (mean LOCO Spearman; % of fits with positive transport):**

| training countries | full index, level | climate+soil, level | full, prev | climate+soil, prev |
|---|---|---|---|---|
| 1 | 0.159 (68%) | **0.318 (98%)** | 0.076 (57%) | **0.231 (89%)** |
| 2 | 0.209 (76%) | 0.361 (100%) | 0.119 (67%) | 0.261 (91%) |
| 3 | 0.265 (81%) | 0.390 (100%) | 0.196 (75%) | 0.291 (94%) |
| slope per country | +0.052 | +0.038 | +0.056 | +0.030 |

**Reading.** The two-domain index starts far higher — trained on a *single*
country it transports positively to another in 98% of fits — and still gains
with each added country, more slowly than the full index because it has
less country-specific noise to average out. Both NCE arguments hold together:
collect the remotely sensed layers, and add countries. The usual caveat: the
domain set was identified on these four countries, so the curve's *level* is
optimistic; its *shape* (near-universal positivity, monotone gain) is the
thing to pre-register.

### 2026-09-03 · CE-01 — how much of the reliability ceiling is cluster effect? (script 31)
**Design.** On units with ≥ 2 survey clusters (13–14 per country), the
split-half ceiling computed two ways on the same units: respondents split at
random *within* the unit (the method of record) vs whole *clusters* split.
The gap is the share of "reliability" that clusters share and geography does
not. 200 splits, Spearman–Brown corrected, uniform outcome definition.
**Result (ceiling within → cluster):**

| | child vit A | women vit A | child iron | women iron | mean inflation |
|---|---|---|---|---|---|
| Gambia | 0.83 → 0.76 | 0.84 → 0.81 | 0.83 → 0.73 | **0.75 → 0.16** | +0.20 |
| Ghana | 0.00 → 0.00 | 0.59 → 0.82 | 0.90 → 0.90 | 0.77 → 0.82 | −0.07 |
| Sierra Leone | 0.00 → 0.00 | 0.41 → 0.30 | **0.26 → 0.00** | **0.53 → 0.00** | +0.23 |
| Malawi | 0.00 → 0.00 | (NA) | 0.69 → 0.75 | 0.77 → 0.81 | −0.03 |

Pooled: within 0.546 → cluster 0.456; **~16% of the within ceiling is
cluster effect** on these units. Heterogeneous: Gambia and Sierra Leone
inflate (iron collapses to near zero under the cluster split); Ghana and
Malawi do not. On multi-cluster units **child vitamin A has no reliable
geography at all** in Ghana, Sierra Leone or Malawi (r_within ≤ 0) — its
file-of-record ceiling there came from single-cluster units.
**Caveats.** 13–14 units per country makes every cell noisy (Ghana women's
vitamin A moves +0.22 the "wrong" way). Folate and B12 skipped: the uniform
resolver returns nothing for them in Ghana and Sierra Leone. The pooled 16%
is the defensible figure; the iron pattern is suggestive, not established.
**Changes.** Headroom re-read against the cluster-split ceiling (HR-02, in
the run log): for Gambia women's iron the model (0.66) now sits *above* the
cluster ceiling (0.16) — which means either that ceiling is too noisy on 13
units to be a bound, or the model is partly predicting cluster-level
structure. Both readings say the same thing for reporting: at 13 units per
country the cluster-split ceiling is a diagnostic, not a quotable bound; the
pooled 16% is the number to carry. The NCE's "two-thirds of the attainable
ceiling" stands as computed but its implied headroom should be described as
an upper bound. P7 in the pre-registration revised to the pooled pattern.

### 2026-09-03 · PS-01 — pre-registered predictions for Ethiopia and Pakistan
Written to [PREREGISTRATION_NEW_COUNTRIES_2026-09.md](PREREGISTRATION_NEW_COUNTRIES_2026-09.md):
seven predictions (rankings transport / levels do not; climate+soil ≥ full
index; each added country buys +0.02–0.06; survey-derived domains do not help
transport; burden margin small; model does not substitute for sample; cluster
split lowers the ceiling), each with script, metric and threshold fixed in
advance, plus what would falsify the central claim.

### 2026-09-03 · MW-01 — Malawi at its district rung: does the headroom survive? (script 31, `CE_MALAWI_DISTRICT=1`)
**Q.** HR-01 put nearly all remaining headroom in Malawi (achieved 0.12 vs
ceiling 0.65), computed on 87 Traditional Authorities of which 74 are single
clusters. At Malawi's district rung (26 multi-cluster districts, ~3 clusters
each), what does the ceiling look like when clusters are split?
**Result (ceiling within → cluster):** child vitamin A **0.56 → 0.00**;
women's vitamin A 0.00 → 0.00; child iron **0.76 → 0.64**; women's iron
0.06 → 0.00. Mean 0.343 → 0.161; **27% of the within ceiling is cluster
effect** pooled over the four countries with Malawi at this rung.
**Reading.** At the rung Malawi should be analysed on, only child iron has
reliable district geography. Child vitamin A's apparent district ceiling is
entirely what clusters share; women's vitamin A and iron have none. Malawi's
outsized "headroom" in HR-01 was therefore mostly single-cluster artefact,
not unexploited signal. This is the strongest single reason to treat the
NCE's headroom sentence as an upper bound and not to promise large Malawi
gains. (Zinc, folate and B12 could not be scored: the uniform resolver
returns nothing for them — see note below.)

### 2026-09-03 · AG-01 — is the regional 12/12 sensitive to aggregation weights? (script 32)
**Result (LOCO Spearman, domain_index, 12 cells per target):**

| aggregation | level | prev | positive cells |
|---|---|---|---|
| outcome by effective n, predictors simple mean (record) | 0.556 | 0.498 | 12/12, 12/12 |
| outcome by population, predictors simple mean | 0.551 | 0.503 | 12/12, 12/12 |
| outcome AND predictors by population | 0.402 | 0.444 | 12/12, 12/12 |

**Reading.** The positivity result is robust to every scheme. The *level* is
insensitive to how the outcome is weighted (Δ ≤ 0.005) but drops 0.05–0.15
when predictors are population-weighted — a region's environment is better
represented by its area than by where its people live. Protocol note: keep
predictors area-averaged. Added to the pre-registration protocol.

### 2026-09-03 · RS-01 — are the zinc, folate and B12 targets uniformly defined? (resolved: yes, effectively)
**What looked wrong.** `resolve_uniform_outcome()` returns NULL for zinc,
folate and B12, and `01_build_targets_v2.R` falls back to each survey's own
binary for them, so `targets_v2.csv` mixes harmonised (vitamin A, iron) and
survey-defined (zinc, folate, B12) targets.
**Why.** Not missing cutoffs — every config carries one (folate < 10 nmol/L,
B12 < 148 pmol/L, zinc < 65/66 µg/dL). The resolver is gated on
`UNIFORM_TRANSPORT_TAGS <- c(child_iron, women_iron, child_vitA, women_vitA)`
(`R/admin2_analysis.R:20`): uniform re-derivation was deliberately limited to
the transport outcomes.
**Does it matter?** Checked directly — survey binary vs the configured cutoff
applied to the continuous column:

| cell | survey binary | cutoff-derived | agreement |
|---|---|---|---|
| Ghana women's folate (< 10 nmol/L) | 54.7% | 54.7% | 100% |
| Ghana women's B12 (< 148 pmol/L) | 8.5% | 8.2% | 99.8% |
| Malawi women's folate | 20.8% | 20.8% | 100% |
| Malawi women's B12 | 2.7% | 2.7% | 100% |
| Malawi child zinc (< 65 µg/dL flat) | 59.1% | 67.7% | 91% |
| Malawi women's zinc (< 66 µg/dL flat) | 61.8% | 74.3% | 88% |

Folate and B12 are already uniform in practice. Zinc differs because the
survey's `zinc_def` applies IZiNCG cutoffs stratified by age, sex, time of
day and fasting, which a flat cutoff over-calls — the survey definition is
the *better* one, and zinc is Malawi-only so cross-country uniformity does not
arise. **No change to the targets; HV-01 withdrawn from the queue.** One
useful by-product: the zinc definition depends on time of collection, so the
raw MNS files must record it — the ZN-02 artefact test is feasible with them.

### 2026-09-03 · NC-01 — how surprising is "N of N cells positive" under no transport? (script 33)
**Design.** Keep predictors, folds, PC orientation and predictions fixed;
permute each held-out country's outcome across its units — the *same*
permutation applied to all of that country's outcomes, so the within-country
correlation between outcomes (the reason cells are not independent) is
preserved while any covariate–outcome relation is destroyed. 500 replicates;
biomarker-level target; domain_index; all outcomes with ≥ 3 countries (22
cells per tier, thin countries included).
**Result:**

| tier | observed positive | observed mean ρ | null positive (median / 95th pct) | null mean ρ 95th pct | p(positive) | p(mean ρ) |
|---|---|---|---|---|---|---|
| regional | 17 / 22 | 0.292 | 11 / 15 | 0.159 | 0.010 | **0.000** |
| district | 17 / 22 | 0.252 | 11 / 16 | 0.080 | 0.022 | **0.000** |

**Reading.** The transport signal is unambiguous on the *mean* rank
correlation: not one of 500 null replicates reached the observed value at
either tier. But the **count of positive cells has a fat null** — with
correlated outcomes inside a country, "no transport" still yields 15–16 of 22
positive cells 5% of the time, so 17/22 is significant only at p ≈ 0.01–0.02,
and a "12 of 12" on a smaller cell set would be less surprising than it
sounds. The impressive-looking statistic is the weaker one.
**Changes.** Report transport as **mean Spearman against its permutation
null** (0.29 vs a 95th percentile of 0.16 regional; 0.25 vs 0.08 district),
and give cell counts as a secondary descriptive. Added to the pre-registration
(P1 now carries a null-calibrated criterion) and to the NCE implications note.

### Queue (in priority order)
- **ZN-02** (blocked on data) Zinc collection-artefact check needs the raw MNS
  field metadata (interview date, team, time of draw); the store has none.

- **ZN-02** Zinc collection-artefact check. The store carries no interview
  date/team/time columns for the Malawi zinc data (1,397 columns, none are
  collection metadata), so the direct test needs the raw MNS field files —
  a data request. Proxy check from the existing reliability schemes below.
- **SL-04** Production-metric check of the rank default: rerun
  `scripts/run_sl_prescreened_main.R` with `meta = "rank"` and compare MAE /
  Pearson / Spearman against the MSE run, before flipping the default.
- **CF-01** Comparator fairness: how much of "regional averages ≈ chance" is
  the jackknife over-correcting at 3–5 units per region? Compare jackknife vs
  shrunk (empirical-Bayes-style) vs split-half regional means, on Spearman and
  burden capture. Directly affects NCE wording.
- **DA-01** Domain ablation under leave-one-country-out: which of the 20
  domains carry transport? Informs what to collect and the RA annotation
  priorities.
- **QA-01** Regression test for the signal-scan sign convention (legumes,
  cattle must be positive) so the 2026-09-02 inversion cannot recur silently.

---

## State at the end of 2026-09-03

Open items in the queue above are the only ones left; everything else ran and
is recorded. Nothing in `_targets.R` was touched. Files that changed today
outside `scripts/protocol_v2/` and `results/tables/protocol_v2/`:
`R/area_superlearner.R` (new), `R/benchmark_models.R`
(`fit_predict_sl_prescreened` → weighted/blocked helper; `ASL_META` switch),
`scripts/run_sl_prescreened_main.R` (loads code via `tar_source`),
`results/tables/sl_prescreened_main.csv` (regenerated, MSE meta-learner),
`archive/` unchanged. Documents: this log, `PROTOCOL_V2.md` (addendum),
`TWO_READINGS_2026-09d.md` (four status rows; superseded by `TWO_READINGS_2026-09e.md`, both abstracts rewritten), `PREREGISTRATION_NEW_COUNTRIES_2026-09.md`
(new), `NCE_IMPLICATIONS_2026-09-03.md` (new). Memory notes updated:
`superlearner_package_choice`, `transport_domains_climate_soil`,
`survey_weighting_issues`, `signal_scan_sign_convention`.

Where to pick up: (1) decide the NCE wording changes in
`NCE_IMPLICATIONS_2026-09-03.md`; (2) the spawned task to migrate the three
`area_level_comparison.R` SuperLearner sites; (3) when Ethiopia/Pakistan data
arrive, score the seven pre-registered predictions; (4) request the raw MNS
field metadata (interview date, team, time of draw) to run ZN-02 properly.


---

# 2026-09-04

## BUG-01 · Sierra Leone silently dropped from every population-joined result

`dashboard/data/admin2_population.rds` spells the country `"Sierra Leone"`;
`targets_v2.csv` and the predictors spell it `"SierraLeone"`. Every script that
joined population and then filtered on `is.finite(pop)` lost the country with no
message: 12 (burden capture, the NCE 24/19/20 numbers), 16 (admin-1 transport,
the "12 of 12" regional headline: 4 outcomes x 3 countries, not 4), 21 (SL-03),
22 (CF-01), 27 (R6-02: 132 units, not 146), 28, 32 (AG-01). Found when AR-01
returned 12 cells. All seven scripts now normalise `POP$country` after loading
and were re-run on 2026-09-04 (results below as they complete). NC-01 (no
population join) was never affected, which is why its regional mean of 0.29
over 22 cells never matched the 0.50-0.56 over "12".

## VC-01 · Variance-components ceiling (script 34)

Respondent in cluster in district in region, REML, all units. **Validation:**
with the cluster counted as geography the VC ceiling reproduces the published
split-half `r_max_emp` (0.605 vs 0.612 over 21 prevalence cells, r = 0.79).
**Honest ceiling** (cluster removed), prevalence, Admin-2 rung: mean 0.435
against 0.543 with cluster-as-geography; cluster share 21%; 6 of 24 cells below
0.30. By country (prev / level): Gambia 0.70 / 0.69, Ghana 0.54 / 0.66, Malawi
0.40 / 0.52, Sierra Leone 0.21 / 0.45. **Child zinc** has essentially no
district geography: honest ceiling 0.27 on prevalence and 0.00 on the level
against 0.80 / 0.78 with cluster counted (Malawi district rung: 0.46 / 0.16).
**Headroom** (achieved in-fill Spearman / honest ceiling) reaches or exceeds 1
in Gambia child vitamin A, Gambia women's iron and Ghana child vitamin A on
prevalence, and Gambia child and women's vitamin A on the level; elsewhere
0.4-0.7. Fourteen of 24 prevalence fits are singular (a region or district
component at zero). Replaces the "two-thirds of attainable accuracy" sentence:
on the honest ceiling the model is at the ceiling in a third of cells and the
remaining headroom is smaller than the split-half implied.
-> `variance_components_ceiling.csv`

## UR-01 · Urbanicity conditioning (script 36)

Composite from night lights, population density, built surface, urban
land-cover fraction and travel time to healthcare (negated). LOCO, climate+soil
index: district level 0.368 raw, 0.355 partial given urbanicity, 0.367
residualised with training-country coefficients; regional level 0.452 / 0.414 /
0.413; prevalence 0.268 / 0.253 / 0.266. **Urbanicity alone does not transport**
(mean -0.07, 8 of 22 positive) and the climate+soil index is barely urban
(median Spearman with the composite +0.13; the full index is anti-urban at
-0.30 to -0.58). **The transported signal is not an urban-rural map.**
-> `urbanicity_conditioning.csv`

## FW-01 · Fieldwork windows (script 37)

Per-cluster dates recovered for all four surveys: Gambia 24 Jan-19 Apr 2018
(70 clusters, Stata day counts), Ghana 15 Apr-15 Jun 2017 (90; 82 in May),
Malawi 8 Dec 2015-15 Feb 2016 (105; from the clean file, 34% undated
respondents), Sierra Leone reconstructed as date of birth + age in days,
3 Nov-4 Dec 2013 (60 clusters, 370 of 486 children). Three surveys sit in the
dry/post-harvest season, Malawi in the lean season.
-> `fieldwork_windows_cluster.csv`, `fieldwork_windows_admin2.csv`

## FP-02 · Time-matched food prices (script 38)

23 columns from the WFP market series matched to each district's fieldwork
months: seasonal amplitude, ALPS-style anomaly z, seasonal position, spatial
relative price, USD/kg staple level, nutritious-to-staple relative prices,
months to the staple price peak; the same columns at cluster GPS points.
Coverage is complete over surveyed districts for Gambia and Sierra Leone;
Ghana's window has no pulses/animal/vegetable observations and Malawi's no
animal/vegetable ones. Found while building it: script 07 used 2021 as
Gambia's survey year (fieldwork was Jan-Apr 2018), so the nine `fprice_*`
columns in the vocabulary are three years late for Gambia. Scored by the
add-on harness (script 39, AD-food) - result below when complete.

## ZN-02 · Zinc collection artefact, preliminary (script 40, children)

Malawi clean file, time of draw and interview date. Afternoon draw lowers
serum zinc 3.5% (t = -2.0); date +0.7%/day (t = 2.2). Variance components of
log zinc: district 0, TA 0, cluster 0.013, residual 0.070 - **there is no
district geography in child zinc to explain**; what the split-half called
geography is cluster effect, and the share of afternoon draws does not track
district prevalence (Spearman -0.05). Women and school-age children pending
the rerun.

## Design note

`CLUSTER_LEVEL_DESIGN_2026-09.md`: how cluster linkage would work (the buffer
path already exists in `data/GEE/DHS GEE merge.Rmd`), the field's conventions
(DHS 2/10 km buffers; LBD/MAP 5 km annual + synoptic-with-lags; Fourier
seasonality), a ~100-column recipe, and the correction that clusters number
328, a 1.6-2.2x gain over district units, not "several hundred per country".

## BUG-01 · outcomes of the re-runs

Unchanged (Sierra Leone's 14 districts cannot be folded for in-fill, so its
absence never entered): R6-02 (index 0.311 / capture 0.204 vs jackknifed
regional mean 0.140 / 0.188), the NCE burden numbers (24 / 19 / 20 percent),
CF-01 (8 of 18). **Changed:** admin-1 transport (script 16) is now 22 cells,
domain_index mean 0.292 (median 0.40), 17 positive on the level and 0.318, 18
positive on prevalence (0.310 and 0.316 after the 08:17 re-run, RR-01 below); 0.431 and 12 of 12 in the two countries with eight or
more regions (Ghana, Malawi). The former "12 of 12 at 0.50-0.56" was Gambia,
Ghana and Malawi x four outcomes. AG-01 (script 32) on 22 cells: n_eff
aggregation 0.292 / 0.318, population 0.289 / 0.280, both-sides population
0.224 / 0.221 - n_eff still best, and "12 of 12" now reads 17 of 22. DA-03
(script 28) reverses: climate + soil 0.452 (20 of 22) vs full 0.292 on the
level, 0.378 vs 0.318 on prevalence - the two-domain index is better at the
regional tier as well. `TWO_READINGS_2026-09f.md`, the pre-registration and
`NCE_IMPLICATIONS_2026-09-03.md` item 4 carry the corrections.

## AR-01 · Anchor-and-rank survey design (script 35), 22 cells

A national anchor from a fraction f of the survey's effective sample plus the
transported ranking (BLP shrinkage: rho_train x sd_train x z), against a
district survey and a regional survey of the same size, scored on the full
survey's district prevalences (consistent rung, Malawi at 27 districts).
Median district MAE at f = 0.05 (50-70 respondents): anchor + ranking 10.4 pp,
regional survey 12.3, district survey 19.6; at f = 0.25: 8.7 / 8.0 / 10.0; the
district survey overtakes the anchored design at a median f of 0.40 (about
394 respondents), in every cell somewhere on the grid. Regional anchors
(A2) overtake the national one at f ~ 0.15-0.25. Anchor + ranking beats the
district survey on MAE in 22 of 22 cells at f = 0.05 and 17-20 of 22 at
f = 0.25, and the regional survey in 19-21 of 22. **It loses on burden
captured at every size** (0.15 vs 0.29 for the district survey, 0.24 for the
regional one) because burden sits in populous districts that a survey
measures precisely; and its Spearman (0.25-0.37) is passed by the district
survey at f ~ 0.05-0.10. The comparison is conservative for the anchored
design (the district survey converges to the truth by construction).
-> `anchor_and_rank.csv`, `_summary.csv`, `_crossover.csv`

## AD-food, AD-climate, AD-alpha · Add-on blocks under one harness (script 39)

Same folds, zero-tuning index, base vocabulary vs base + block vs block alone
vs (block replacing a named set); LOCO at both tiers with and without the
climate + soil restriction.
- **Time-matched food prices (23 cols, script 38):** in-fill +0.005 / +0.004
  (10 and 9 of 18 cells); LOCO +0.002 / -0.003 (admin2), -0.006 / -0.009
  (admin1); block alone -0.08 to -0.12. Replacing the nine `fprice_*` columns:
  +0.005 / +0.006. Scan: the strongest replicated signs are *negative* -
  higher local staple price and higher survey-month anomaly go with less
  deficiency (women's iron meta z -3.0, 3 of 4 countries) - the urban gradient,
  not food access. **No gain.**
- **Time-matched climate (14 cols, script 41; LST monthly from the shared
  table, CHIRPS monthly zonal means for the survey year):** in-fill
  +0.001 / -0.002; LOCO +0.002 / +0.005 (admin2); climate + soil with the
  block -0.015 / -0.003. Scan: night LST in the fieldwork months tracks more
  child vitamin A deficiency in 4 of 4 countries (meta z 3.3) and rainfall
  amplitude tracks less women's vitamin A (4 of 4, -3.9) - real associations
  that add nothing once the annual columns are in. **No gain.**
- **AlphaEarth embedding as its own domain (64 cols, script 42):** already in
  the vocabulary as `aef_A00..63`, filed under "Agricultural production, land
  use" (64 of 93 columns), so DA-01's "agriculture" was mostly the embedding.
  As its own domain: in-fill 0.000 / -0.003; LOCO replace +0.009 / +0.006
  (admin2), +0.019 / -0.002 (admin1); embedding alone 0.20 at district level
  (18 of 22), below climate + soil's 0.37; climate + soil + embedding worse
  than climate + soil (-0.03 / -0.04). Strongly replicated single-dimension
  associations (90 of 384 column x outcome pairs agree in sign across all
  countries; sat_A03 with women's vitamin A meta z 6.7). **No gain**; fix the
  metadata label.
-> `addon_{food_tm,climate_tm,alphaearth}_{scan,infill,loco}.csv`

## ZN-02 · Zinc collection artefact, complete (script 40)

Afternoon draw lowers serum zinc in every group: children -3.5% (t -2.0),
school-age -7.3% (t -3.3), women -5.5% (t -3.0); calendar date +0.7-1.3% per
day where dated (children, school-age). Variance components of log zinc,
district | TA | cluster | residual: children 0 | 0 | 0.013 | 0.070; school-age
0.003 | 0.006 | 0.009 | 0.078; women 0.002 | 0.004 | 0.011 | 0.059. Adjusting
for time of draw and date leaves the district component unchanged. So the
time-of-draw effect is real at the individual level but **is not what
separates districts - almost nothing does**; the between-district share of
zinc variance is 0-3%, the cluster share 10-16%. The women's district
correlation between share of afternoon draws and prevalence is -0.44 (the
ecological sign is the reverse of the individual one), a team/area confound
rather than an explanation. Team identifiers exist for Ghana (10) and Sierra
Leone (8) but not Malawi.
-> `zinc_collection_artefact.csv`

## State at the end of 2026-09-04

New scripts 34-42 in `scripts/protocol_v2/` (all parse and ran); new tables
in `results/tables/protocol_v2/`; new feature files in
`data/covariates/harmonized/` (`predictors_admin2_food_tm.csv`,
`predictors_cluster_food_tm.csv`, `predictors_admin2_climate_tm.csv`,
`predictors_admin2_alphaearth.csv`, each with metadata) - none merged into
`predictors_admin2_shared.csv`. Seven earlier scripts patched for BUG-01 and
re-run. Documents: `TWO_READINGS_2026-09f.md` (current), banner on e,
`NCE_IMPLICATIONS_2026-09-03.md` item 4 corrected and 4b added,
pre-registration amended, `CLUSTER_LEVEL_DESIGN_2026-09.md` (new). Memory:
`population_file_sierra_leone_spelling` (new), `two_readings_current`,
`transport_domains_climate_soil`, `protocol_v2`, `survey_weighting_issues`
updated. `_targets.R` untouched; nothing committed.

Decisions for the user: (1) fix the AlphaEarth domain label in
`predictors_admin2_shared_metadata.csv` (changes the agriculture PCs and so
every domain-PC result); (2) fix script 07's Gambia survey year (2021 -> 2018)
and rebuild the nine `fprice_*` columns; (3) whether the manuscript body
adopts revision f's regional figure; (4) the cluster-level build in
`CLUSTER_LEVEL_DESIGN_2026-09.md`.

*All four were taken up the same day at the user's direction; see RR-01 and
the cluster-level entries below.*

## RR-01 · Full re-run after the AlphaEarth relabel and the Gambia food-year fix

`scripts/covariates/build_shared_predictor_set.R` now files `aef_*` under
"Satellite embedding" (source "AlphaEarth (GEE)") and the live metadata was
relabelled to match; script 07 takes Gambia's survey year as 2018 and the
nine `fprice_*` columns were rebuilt (window 2016-2020). Then the four
in-fill shards, the merge/LOCO step and fourteen downstream scripts were
re-run (rerun_all.ps1, 08:03-08:17). What moved, on the full 22/24 cells:

| Figure | Before | After |
|---|---|---|
| Transport, district, zero-tuning index (level / prev) | 0.252 / 0.151 | 0.261 / 0.159 |
| Transport, district, penalised domain fit (level) | 0.281 (21 of 22) | 0.244 (21 of 22) |
| Transport, regional, full index (level / prev) | 0.292 / 0.318 | 0.310 / 0.316 |
| NC-01 null 95th percentile, regional / district | 0.159 / 0.080 | 0.166 / 0.079 |
| Climate + soil, district / regional (level) | 0.368 / 0.452 | 0.368 / 0.452 |
| In-fill index vs jackknifed regional mean (level; prev) | 0.398 vs 0.320; 0.286 vs 0.193 | 0.388 vs 0.320; 0.285 vs 0.193 |
| Burden captured, index / jk / spatial | 0.240 / 0.190 / 0.229 | 0.241 / 0.190 / 0.229 |
| DA-01 load-bearing domains (level) | soil +0.041, climate +0.034 | soil +0.036, climate +0.033, agriculture +0.016, satellite embedding +0.012 |
| TC-01 slope per training country (level / prev) | +0.052 / +0.056 | +0.052 / +0.060 |
| DA-02 nested selection vs full (level, median) | +0.028 | +0.018 |

The satellite embedding, once its own domain, is mildly load-bearing on the
level target and harmful on prevalence (-0.012); "agriculture" without it is
+0.016 / +0.021. Nothing else changed by more than 0.01. The manuscript,
revision f, the NCE note and the memory notes now carry 0.31 for the regional
figure and quote the zero-tuning index (0.26, 17 of 22) rather than the
penalised fit (0.24, 21 of 22) for the district figure.

## CL-01 / CL-02 / CL-03 · Cluster-level analysis, first pass (scripts/cluster_level/)

A parallel track, not a replacement: same three estimands, same arms, same
conventions, fitted at the survey cluster. `01` builds the cluster outcome
table (2,002 rows; 323 clusters, all with GPS; median 8-19 respondents per
cluster; the cluster is the PSU so n_eff is the Kish n). `02` extracts 117
covariates at 2 km (urban, GHSL SMOD >= 21) / 5 km (rural) buffers from
rasters already on disk - no Earth Engine call - with climatology, amplitude,
peak month and fieldwork-window columns for the dynamic layers; 113 of 117
columns exist in all four countries. Three things had to be fixed to get
there: cropping a 30 m country raster before extraction (150 s for eight
buffers; now no crop for country rasters and area weights only for coarse
ones - ~2-4 min per country), Sierra Leone's rasters being named both
"Sierra_Leone" and "Sierra Leone" (44 columns lost), and seven of its
space-named soil GeoTIFFs crashing GDAL silently part-way through a read
(clean single-band LZW copies written in separate processes, `*_clean.tif`).
`03` fits and scores everything twice: across clusters, and after aggregating
cluster predictions to districts (n_eff-weighted), which is the number
comparable with protocol v2.

**Result: fitting at the cluster does not beat fitting at the district
anywhere in this pass, and is worse for transport.** Zero-tuning index,
aggregated to districts, against the district-fitted protocol-v2 value on
the same cells:

| Estimand | Target | Cluster-fitted (agg.) | District-fitted | Cluster better |
|---|---|---|---|---|
| In-fill | level | 0.357 | 0.388 | 10 of 18 |
| In-fill | prev | 0.218 | 0.285 | 8 of 18 |
| Region | level | 0.355 | 0.366 | 10 of 18 |
| Region | prev | 0.201 | 0.268 | 9 of 18 |
| Transport, Admin-2, climate+soil | level | 0.244 (20 of 22) | 0.368 (22 of 22) | 9 of 22 (vs full 0.261) |
| Transport, Admin-2, climate+soil | prev | 0.166 (18 of 22) | 0.268 (21 of 22) | 10 of 22 (vs full 0.159) |
| Transport, Admin-1, climate+soil | level | 0.249 (17 of 22) | 0.452 (20 of 22) | - |
| Transport, Admin-1, climate+soil | prev | 0.303 (20 of 22) | 0.378 (20 of 22) | 9 of 22 (vs full 0.313) |

At the cluster itself the in-fill Spearman is 0.31 (level) / 0.20 (prev):
the cluster outcome is the noisier unit. The spatial smoother at the cluster
(0.31 aggregated) is also below its district-fitted value (0.39). The one
gain is the fieldwork-window block, which is only usable in-country: +0.018
on the level target (0.375 vs 0.357, better in 14 of 24 cells), nothing on
prevalence - consistent with the district-level AD-climate result.

Why, probably: (1) a cluster of 8-19 respondents is a much noisier outcome
than a district, and the index weights domains by an UNWEIGHTED Spearman
with that outcome; (2) the cluster vocabulary is the remotely sensed layers
on disk (iSDA soil elements, TerraClimate, FLDAS, CHIRPS, LST), not the
district vocabulary's 50 climate + 45 soil columns (SoilGrids among them),
so "climate + soil" is not the same block; (3) one radius, no shrinkage, no
tuning. What to try before concluding anything: n_eff-weighted Spearman in
the index; empirical-Bayes shrinkage of cluster outcomes toward their
district before fitting; the district vocabulary's exact layers at the
cluster; a 10 km rural radius; and prediction on a grid rather than at
district-mean covariates. None of these were run today.
-> `results/tables/cluster_level/` (`targets_cluster.csv`,
`benchmarks_cluster_{raw,cells,loco}[_cs].csv`, `cluster_vs_district.csv`),
`data/covariates/cluster/predictors_cluster.csv` (+ metadata).

## DA-04 · Transport ablation by DATA SOURCE (script 43)

The January 2026 Ghana deck asked which data sources contribute and answered
with a change-in-MSE importance from one country's prescreened SuperLearner.
Redone under the protocol: drop every column from one source (or keep only
that source), rebuild the domain PCs, score the zero-tuning index under
leave-one-country-out, 22 cells. Loss when dropped (biomarker level / prev):
SoilGrids-iSDA +0.036 / +0.033, Earth Engine +0.031 / +0.008, Malaria Atlas
+0.018 / +0.008, MapSPAM +0.013 / +0.012, AlphaEarth +0.012 / -0.012, WFP
prices +0.004 / +0.001, Koppen/AEZ ~0. **Dropping the 138 DHS/MICS columns
IMPROVES transport: 0.261 -> 0.307 on the level and 0.159 -> 0.226 on
prevalence.** Alone, soil (0.318) and Earth Engine (0.269) each beat the full
vocabulary (0.261); DHS alone 0.144. Consistent with DA-01 by domain: the
survey-derived indicators are dead weight or worse across borders, the
remotely sensed layers carry it.
-> `source_ablation_loco.csv`, `source_ablation_loco_summary.csv`

## Deck · Ghana presentation updated (docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd)

The January 2026 Ghana deck (26 slides, Dropbox) rebuilt as a Quarto
document rendering to PowerPoint on the same template (reference-doc), 31
slides, every number read from `results/tables/` at render time (protocol
v2 as re-run 08:17). Photos and the Hess et al. framework figure reused from
the old file (`docs/slides/img/`). New live figures: Ghana district maps of
survey prevalence beside held-out index predictions, and the training-country
curve. Rendering lessons: pandoc's pptx writer puts any table or figure on a
slide of its own and moves following text to an untitled slide, so tables
carry their interpretation in speaker notes; the template's body font needs
one-line titles and three short bullets per slide. Render with the RStudio
Quarto binary (`quarto render ... --to pptx`); visual QA via LibreOffice PDF
and PyMuPDF. The deck's "no longer in the pipeline" slides list what the
January deck had that the pipeline does not: individual-level SuperLearner
metrics (sensitivity layer), the national VMNIS/BRINDA track (separate
sub-pipeline, April 2026), iodine, risk-category accuracy, LSMS and FluNet
predictors, source-level importance (re-added as DA-04).

## Add-backs from the January deck (scripts 44-47), at the user's direction

**RC-01 · Risk-category accuracy (script 44).** WHO 2009 vitamin A bands
(< 2 / 2-10 / 10-20 / >= 20 percent) on out-of-fold in-fill predictions
averaged over ten draws, and on anchored transport (held-out country's own
national prevalence + transported ranking, AR-01 design A1). Vitamin A, 8
cells, share of districts in the exact class / within one class: index 0.53
/ 0.89, smoother 0.53 / 0.91, jackknifed regional mean 0.52 / 0.89; at
Admin-1 0.61-0.67 / 0.96-0.98; anchored transport 0.50 / 0.87 (district) and
0.55-0.61 / 0.97-0.98 (region). Ghana child vitamin A in-fill: index 0.40 /
0.67 exact / within-one at district (January deck: 0.44 / 0.75 full model,
0.31 / 0.49 proxies), 0.50 at region. For iron and the others only the 20
percent threshold has a WHO source: two-class accuracy 0.84-0.85 in-fill,
0.75-0.79 under anchored transport, 0.80-0.91 at Admin-1. The class table
adds nothing the rank correlations do not already say, but it is the
programme-facing number and it now exists under the protocol.
-> `risk_category_accuracy.csv`, `_summary.csv`

**AD-lsms · LSMS and FluNet (script 45 + harness 39).** LSMS exists for
Ghana only (GLSS7, a svyby object of 10 old-region means); 134 numeric
region means broadcast to Ghana's districts through the 16-to-10 crosswalk
(260 of 260 mapped). FluNet reports for Ghana and Sierra Leone only and is a
national weekly count: three survey-year indicators built, constant within
country. In-fill, Ghana, biomarker level, base -> base + LSMS: child iron
0.588 -> 0.572, child vitamin A 0.302 -> 0.282, women's B12 0.584 -> 0.600,
folate 0.404 -> 0.380, women's iron 0.394 -> 0.380, women's vitamin A 0.495
-> 0.502; prevalence: B12 0.239 -> 0.291, the rest unchanged or lower. LSMS
alone 0.33 vs base 0.39. FluNet contributes nothing by construction and
cannot be scored across borders. **Neither earns a place**; both would need
to be collected sub-nationally in the next countries to be worth anything.
-> `predictors_admin2_lsms_flunet.csv`, `addon_lsms_flunet_{scan,infill,loco}.csv`

**IO-01 · Iodine within country (script 47).** The Gambia: 1,285 women with
continuous UIC, median 143 ug/L, 32.5 percent below 100; 8 percent of
households with iodised salt. Sierra Leone: 516 non-lactating women in six
WHO categories (labels assumed from the WHO scheme; the codebook label set
could not be extracted), 18.6 percent below 100; 81 percent of households
with adequately iodised salt. District targets (27 Gambia, 14 Sierra Leone).
In-fill, Gambia, UIC < 100: index 0.36, jackknifed regional mean 0.33,
smoother 0.31; median log UIC: index 0.49, regional mean 0.35, smoother 0.27.
Leave-one-region-out: Gambia 0.40-0.51 (smoother best on prevalence, index
best on the level), Sierra Leone index +0.27 with every other arm negative
on 14 districts. Descriptive: district share of iodised salt against UIC
< 100 is -0.48 (Gambia) and -0.77 (Sierra Leone); salt ppm -0.16 and -0.66.
Iodine behaves like the other outcomes in-country; it simply has no
cross-border test.
-> `iodine_targets.csv`, `iodine_in_country.csv`

**IL-01 · Individual-level survey-only vs proxies vs both (script 46).**
Person-level SuperLearner (mean / elastic net / ranger, survey weights,
inner folds blocked by district), outer 5-fold by district, survey columns
capped at 80 by coverage, run in four country shards (Gambia, Ghana, Sierra
Leone complete: 15 cells; Malawi's respondent columns carry no gw_ prefix,
the shard failed and was relaunched with a prefix-free candidate rule). Mean
over 15 cells: AUC 0.50 / 0.50 / 0.51 and Brier skill 0.007 / 0.006 / 0.012
for survey-only / proxies-only / both; PR gain 1.1-1.2. The only cells with
real person-level skill: Ghana child iron (survey-only AUC 0.76, Brier skill
0.14; both 0.77 / 0.15; proxies 0.68 / 0.07), Ghana child vitamin A (both
0.63 / 0.04), Gambia women's iron (proxies 0.61 / 0.04). Sierra Leone: the
discrete SuperLearner picks the mean learner in most folds, i.e. no skill.
Survey-only beats proxies in 6 of 15 cells; both beats survey-only in 9 of
15. **The January figure (Brier skill 20-60 percent, B12 AUC 0.69, PR gain
3.1x) does not survive district-blocked folds and the removal of
biomarker-adjacent survey variables**; what it measured was mostly
within-district similarity and the correlates of the biomarker itself.
Kept as a sensitivity analysis, as before.
*Leakage caught and removed:* Malawi's respondent columns carry no `gw_`
prefix, and the first prefix-free candidate rule let the biomarkers
themselves through (`sf_reg`, `fol_nmol`), giving AUC 0.99 for iron and
folate. That shard was discarded, the script now excludes every configured
outcome variable by name and the `map2_` aggregates, and Malawi was re-run;
its cells enter the table only if the printed column sample is clean.
*Malawi, final:* the store's Malawi data holds only the outcome and district
aggregates, so its cells were built from `clean_malawi_mn_data.RDS`
(questionnaire items m01-m126, hunger score, fortification exposure; the
survey's own outcome definitions; no folate or B12 there). Column sample
clean. Child iron: survey-only AUC 0.78 / Brier skill 0.10, proxies 0.54 /
0.01; women's iron 0.75 / 0.10 vs 0.52 / -0.01; vitamin A nothing. **Final,
19 cells:** AUC 0.51 / 0.50 / 0.53 and Brier skill 0.012 / 0.003 / 0.019 for
survey-only / proxies-only / both; survey-only beats proxies in 8 of 19,
both beats survey-only in 12 of 19. Person-level skill exists for iron in
Ghana and Malawi (AUC 0.75-0.78, from the questionnaire) and nowhere else.
-> `individual_level_models.csv` (+ per-country shards)

## CL-04 · Weighted index x shrunken cluster outcome (script 04)

The two repairs the first pass called for, crossed: domain weights from a
Kish-n-weighted Spearman, and empirical-Bayes shrinkage of each cluster's
outcome toward its district (region for single-cluster districts) on the
modelling scale, with between-cluster variance by the method of moments from
the training fold only. Scored, as before, against the unshrunk district
outcome. **The shrinkage factor tells the story before the scores do:** the
median B is 0.2 on the biomarker level and 0.0 on prevalence - on prevalence
the between-cluster variance within districts is no larger than sampling
noise, so the "shrunk" outcome IS the district (or region) mean, and that
variant is a district-level fit with cluster-level covariates.

Zero-tuning index, aggregated to districts, mean over cells:

| Estimand / set | Target | Raw, unweighted (CL-03) | Best repair | District-fitted |
|---|---|---|---|---|
| In-fill, transportable | level | 0.357 | 0.363 (weighted; 14 of 24 cells) | 0.388 |
| In-fill, transportable | prev | 0.218 | 0.222 (shrunk; 12 of 24) | 0.285 |
| Region, transportable | level | 0.355 | 0.362 (shrunk + weighted; 13 of 24) | 0.366 |
| Region, transportable | prev | 0.201 | 0.252 (shrunk + weighted; 16 of 24) | 0.268 |
| Transport Admin-2, climate + soil | level | 0.244 | 0.269 (shrunk + weighted; 19 of 22) | 0.368 |
| Transport Admin-2, climate + soil | prev | 0.166 | 0.166 (none helps) | 0.268 |
| Transport Admin-1, climate + soil | level | 0.249 | 0.308 (shrunk + weighted; 18 of 22) | 0.452 |
| Transport Admin-1, climate + soil | prev | 0.303 | 0.364 (shrunk; 19 of 22) | 0.378 |

The repairs recover a third to a half of the gap where the outcome is
noisiest (region extrapolation on prevalence, transport at Admin-1) and
almost nothing in-fill; cluster-fitted models still trail district-fitted
ones on every row. Since the prevalence "shrunk" variant carries the
district outcome on every cluster, its failure to match the district fit
says the 2-5 km covariates are not better predictors than district means on
this vocabulary - the resolution argument does not pay off here. What is
left to try is the vocabulary itself (the district set's exact layers,
SoilGrids included), a 10 km rural radius, and grid prediction; on this
evidence the cluster track is a sensitivity analysis, not a replacement.
-> `benchmarks_cluster_ws_{raw,cells,loco}.csv`

## IH-01 · Gambia's missing IHME block, and the re-run after fixing it (RR-02)

A coverage audit after the round was committed (finite share by country x
domain, surveyed districts only) found three holes the column counts hide:

- **Gambia: 25 of 34 IHME/GFDx columns all-NA** (finite share 0.26). Cause:
  `08_build_extra_sources.R` picked ONE year per country, nearest the survey,
  then filtered every indicator to that year. The IHME series end in different
  years (child growth failure, education, WASH, ORS stop at 2017; anaemia,
  EBF, MCV run to 2019); with Gambia's year wrongly set to 2021 the pick was
  2019 and 18 of Gambia's 29 indicators vanished. Fixed: nearest year per
  indicator, and Gambia 2018 (the same year bug as script 07). Gambia is now
  0.70 finite (15 columns recovered; 41 of 41 IHME names match). The 10
  still-empty columns are source gaps: GFDx has no Gambia rows, and IHME does
  not model circumcision or the s* sanitation indicators for Gambia.
- **Ghana: 26% NA on IHME columns** among surveyed districts, from name
  matching against the 2019 district splits (185 of 293 IHME names exact, 13
  unmatched after fuzzy matching). Not fixed here; the fix is zonal means of
  the IHME GeoTIFF surfaces already on disk (`data/IHME/*/GeoTIFF`), as the
  cluster track does.
- **Cluster track: WorldPop as four country-specific columns** (the band name
  is the file name), all-NA outside their own country and so excluded from the
  shared vocabulary. Coalesced into one `worldpop` column (323 of 323 clusters
  finite); `02` patched so the next extraction does this itself.

Expected, not bugs: Ghana has no sibling-history block (its survey is not a
DHS); the time-matched food and climate features exist only for dated,
surveyed districts, by construction; WFP carries no pulses or animal series
for Ghana and no animal or vegetable series for Malawi;
`dhs_w_barrier_transport` is empty in all four countries; ESPEN (helminths)
was never populated - `data/ESPEN` is empty.

Everything downstream of the shared set was then re-run (rerun_all.ps1,
22:45-23:25; scripts 43, 44, 47; the food and climate add-ons; the four
individual-level shards; cluster 03 and 04). What moved, RR-01 -> RR-02:

| Figure | RR-01 | RR-02 |
|---|---|---|
| Transport, district, zero-tuning index (level / prev) | 0.261 / 0.159 | 0.261 / 0.164 |
| Transport, district, penalised domain fit (level) | 0.244 (21 of 22) | 0.273 (19 of 22) |
| Transport, regional, full index (level / prev) | 0.310 / 0.316 (17 / 17 positive) | 0.300 / 0.312 (17 / 18) |
| Regional, the two countries with >= 8 regions (level) | 0.43, 12 of 12 | 0.42, 12 of 12 |
| NC-01 null 95th percentile, regional / district | 0.166 / 0.079 | 0.161 / 0.080 |
| Climate + soil, district / regional (level) | 0.368 / 0.452 | 0.368 / 0.452 |
| In-fill index vs jackknifed regional mean (level; prev) | 0.388 vs 0.320; 0.285 vs 0.193 | 0.389 vs 0.320; 0.286 vs 0.193 |
| Region estimand, index (level / prev) | 0.366 / 0.268 | 0.369 / 0.272 |
| Burden captured, index / jk / spatial | 0.241 / 0.190 / 0.229 | 0.241 / 0.190 / 0.229 |
| DA-01 load-bearing domains (level) | soil +0.036, climate +0.033, agriculture +0.016, embedding +0.012 | soil +0.037, climate +0.028, agriculture +0.015, malaria +0.013, embedding +0.010 |
| TC-01 slope per training country (level / prev) | +0.052 / +0.060 | +0.052 / +0.061 |
| DA-02 nested selection vs full (level, median) | +0.018 | +0.022 |
| AR-01 at f = 0.05: A1 / regional survey / district survey MAE (pp) | 10.4 / 12.5 / 19.6; crossover f 0.40 | 10.4 / 12.5 / 19.6; crossover f 0.40 (394 respondents) |
| IL-01 AUC survey / proxies / both; Brier skill | 0.51 / 0.50 / 0.53; 0.012 / 0.003 / 0.019 | 0.51 / 0.51 / 0.53; 0.012 / 0.003 / 0.018 |
| Cluster track, index transport at Admin-2 (prev / level) | 0.139 / 0.196 | 0.135 / 0.194 |

Only Gambia's cells and the leave-one-country-out folds that train on Gambia
moved; the largest single-cell changes are Gambia women's vitamin A at the
regional tier (0.89 -> 0.77 on six units) and Malawi women's vitamin A
(0.39 -> 0.27). The penalised fit gained 0.03 on the mean and lost two
positive cells. The food-price and climate add-ons stay within +/-0.005 of
"no gain" on every comparison (they had been scored at 07:20 against the
pre-relabel set; they now match the current set). The AlphaEarth add-on was
not re-run: `aef_*` has been its own domain since the relabel, so "plus"
duplicates the base and "replace" equals it; DA-04's "only AlphaEarth" (0.20)
is the remaining informative arm. The individual-level result moved only in
Gambia (women's vitamin A AUC 0.46 -> 0.67, child iron 0.53 -> 0.58); the
conclusion - skill only for iron, none for the vitamins - is unchanged.
Cluster in-fill and region results are identical to three decimals; the
coalesced WorldPop column adds nothing.

Documents updated to 0.30 / 0.27 (19 of 22) / null 0.16: TWO_READINGS_2026-09f,
NCE_IMPLICATIONS_2026-09-03, manuscript_mcn.qmd, the deck (re-rendered) and
the memory notes. The RR-01 entry above is left as the record of that run.
-> `scripts/protocol_v2/08_build_extra_sources.R`,
`scripts/cluster_level/02_extract_cluster_covariates.R`; all tables under
`results/tables/protocol_v2/` and `results/tables/cluster_level/` refreshed.

## NL-01 · Four new layers: livestock, water and coast distance, helminths, IHME surfaces (scripts 48, 49, covariates/gee_water_coast_distance.py; RR-03)

Requested 2026-09-07: the Ghana IHME spatial join, the ESPEN download, and
livestock density and distance-to-water layers from Earth Engine. All four
are built, validated, scored as add-ons, and now in the shared set (451 ->
473 columns, 21 -> 24 domains); the cluster table carries the same layers
(`scripts/cluster_level/05_merge_new_layers.R`).

- **IHME surfaces by zonal mean (script 48).** Every IHME indicator with a
  5 km GeoTIFF on disk (`data/IHME/*/GeoTIFF`) rebuilt as a WorldPop-weighted
  zonal mean over the spine polygons at the nearest year to each survey.
  Against the name-joined rollup on matched districts: Spearman 0.59-0.98 for
  growth failure, anaemia, EBF, HIV, MCV1 and the WASH surfaces (Malawi's
  TA-level values against a district value broadcast, 0.36-0.55). The
  diarrhoea rate surfaces do not reproduce `ihme_incidence` / `ihme_deaths`
  (rho about 0) and stay tabular, as do the six indicators without a surface.
  Script 08 swaps in 17 columns where rho >= 0.5. Ghana IHME missingness
  among surveyed districts 26% -> 0; Gambia's Infant/child domain 33% -> 12%.
- **Livestock density (script 48).** GLW4 2015 dasymetric 5-arc-minute
  rasters (Harvard Dataverse, CC-BY 4.0) for cattle, sheep, goats, pigs and
  chickens: head per km2, tropical livestock units per km2 and per person
  (WorldPop), ruminant share of TLU. Eight columns, complete everywhere.
- **Water and coast distance (`scripts/covariates/gee_water_coast_distance.py`).**
  JRC Global Surface Water occurrence (>= 50% = permanent, >= 10% = any) and
  LSIB land polygons; fastDistanceTransform in Web Mercator corrected by
  cos(latitude); polygon mean and minimum, and the buffer mean at clusters.
  Two pitfalls cost an hour: layers reprojected at different scales must be
  reduced separately (a composite returned thousands of km for inland
  districts), and a single-band reduceRegions names its outputs `mean`/`min`.
  Validated: Accra 2 km, Cape Coast 9, Kumasi 176, Tamale 415, Bolgatanga
  565; minimum permanent-water distance is 0 in 99% of districts with
  permanent-water land cover; coast distance vs elevation Spearman 0.81.
- **Helminths (script 49).** ESPEN's portal export needs no key
  (`/api/download-data/{ISO2}/{sth|sch}/iu/{from}/{to}`, found in the site's
  route table; the keyed API stays closed): eight implementation-unit files,
  2014-2025. No continuous prevalence, only the programme's endemicity class
  and MDA delivery / coverage per IU-year, so the features are the class
  midpoint at the survey year and at baseline, the share of 2014-18 with MDA,
  and mean coverage. IU = ADM2 in all four (Malawi = district, broadcast);
  aliases for post-split districts (Accra Metropolitan Area -> Accra, Wa,
  Ho, Basse / Jimara / Tumana -> Fulladu East, Karene -> Bombali, Falaba ->
  Koinadugu). Coverage 95-100%. Scan against the iron outcomes: STH class
  rho +0.07 to +0.09, MDA share -0.12 (children) / -0.19 (women): the
  expected signs, weak.

**Add-on tests (script 39, base = the 451-column set):** every block within
+/-0.005 of the base on plus-vs-base at every tier. Alone: livestock matches
the full vocabulary at district transport on the level target (median
+0.025, 12 of 22 better); the IHME raster replacement is +0.003 / +0.004 at
district level (13-14 of 22 better); water and coast alone -0.10 to -0.37;
helminths alone -0.19 to -0.45. -> `addon_{espen,livestock,ihme_raster,water_distance}_*.csv`

**RR-03, everything downstream re-run on the 473-column set** (rerun_all.ps1
12:08-12:53, scripts 43/44/47, four individual-level shards, cluster 03-05).
What moved, RR-02 -> RR-03:

| Figure | RR-02 | RR-03 |
|---|---|---|
| Transport, district, zero-tuning index (level / prev) | 0.261 / 0.164 (17 / 15 positive) | 0.268 / 0.184 (17 / 16) |
| Transport, district, penalised domain fit (level) | 0.273 (19 of 22) | 0.280 (21 of 22) |
| Transport, regional, full index (level / prev) | 0.300 / 0.312 (17 / 18) | 0.321 / 0.326 (17 / 17) |
| Regional, the two countries with >= 8 regions (level) | 0.42, 12 of 12 | 0.43, 12 of 12 |
| NC-01 null 95th percentile, regional / district | 0.161 / 0.080 | 0.148 / 0.081 |
| Climate + soil, district / regional (level) | 0.368 / 0.452 | 0.368 / 0.452 |
| In-fill index vs jackknifed regional mean (level; prev) | 0.389 vs 0.320; 0.286 vs 0.193 | 0.392 vs 0.320; 0.291 vs 0.193 |
| Region estimand, index (level / prev) | 0.369 / 0.272 | 0.373 / 0.275 |
| Burden captured, index / jk / spatial | 0.241 / 0.190 / 0.229 | 0.243 / 0.190 / 0.229 |
| DA-01 load-bearing domains (level, drop delta) | soil +0.037, climate +0.028, agriculture +0.015, malaria +0.013, embedding +0.010 | soil +0.021, agriculture +0.010, climate +0.010, livestock +0.008 |
| DA-04 dropping the 138 DHS columns (level / prev) | 0.261 -> 0.307 / 0.164 -> 0.227 | 0.268 -> 0.326 / 0.184 -> 0.239 |
| TC-01 slope per training country (level / prev) | +0.052 / +0.061 | +0.051 / +0.063 |
| DA-02 nested selection vs full (level, median) | +0.022 | +0.014 |
| AR-01 at f = 0.05: A1 / regional survey / district survey MAE (pp) | 10.4 / 12.5 / 19.6; crossover f 0.40 | 10.4 / 12.5 / 19.6; crossover f 0.40 |
| RC-01 vitamin A bands, in-fill exact / within one | 0.53 / 0.89 | 0.52 / 0.90 |
| UR-01 partial-vs-raw (climate + soil), district / regional | +0.006 / +0.003 | -0.006 / -0.016; full index vs urbanicity -0.35 / -0.31 |
| IL-01 AUC survey / proxies / both; Brier skill | 0.51 / 0.51 / 0.53; 0.012 / 0.003 / 0.018 | 0.51 / 0.51 / 0.53; 0.012 / 0.001 / 0.017 |
| IO-01 Gambia in-fill index, UIC < 100 / median UIC | 0.37 / 0.49 | 0.35 / 0.48 |
| Cluster track, index transport at Admin-2 (level / prev) | 0.194 / 0.135 | 0.241 / 0.173 |
| Cluster track, in-fill index, transportable set (level / prev) | 0.312 / 0.201 | 0.319 / 0.206 |

Reading: the three domains and the raster IHME block lift the full index by
0.007 (district level), 0.020 (district prevalence) and 0.02 (regional
level), and the penalised fit regains 21 of 22 positive cells; the
climate + soil index is unchanged by construction and still leads at both
tiers (0.37 / 0.45), so the pre-registered P2 stands. The per-domain drop
deltas shrink as domains multiply (each is a smaller share of the PCs); soil
stays first, livestock enters fourth on the level. The cluster track gains
most (transport 0.19 -> 0.24 on the level) because the new layers are a
larger share of its 117-column vocabulary. Nothing else moved by more than
0.02. Documents updated to 0.27 / 0.28 (21 of 22) / 0.32 / null 0.15:
TWO_READINGS_2026-09f, NCE_IMPLICATIONS_2026-09-03, the pre-registration
note (amendment dated 2026-09-07), manuscript_mcn.qmd, the deck
(re-rendered) and the memory notes.

Not done, for the record: the GLW4 2020 release exists only on the FAO
catalog (2015 used); ESPEN site-level survey prevalence is not exposed on the
export route; the remaining Gambia gaps are source gaps (GFDx has no Gambia
rows; IHME does not model circumcision or the s* sanitation indicators
there).
