# Micronutrient prediction: methods, findings, and interpretation

**Prepared for independent methods review.**
Date: 2026-08-31. Repository: `mn-prediction`, branch `covariate-harmonisation-and-honest-baselines`.
Relevant commits from the session this documents: `fe30fa8`, `59f4510`, `9eb8bb6`, `b15d723`.

---

## 0. How to read this document, and what to attack

This is a working document written to be **falsified**, not admired. The project's
headline results are negative — proxy covariates do not predict district-level
micronutrient deficiency — and negative results are exactly the kind that a
methods bug can manufacture. Three separate leakage or join defects were found
during the session this documents, each of which had already produced a number
that went into a document or a dashboard. That is the base rate a reviewer
should assume for the rest.

The highest-value review targets, in order:

1. **§5 — the individual-level anchor.** This is the single load-bearing result:
   it is what licenses the claim that the ceiling is the *target*, not the
   *predictors*. It was wrong twice before it was right (r = 0.973, then
   r = 0.986, both leakage). If it is still leaking, or if the questionnaire arm
   is unfairly handicapped, the project's central interpretation collapses.
2. **§7 — the assay guard** (`is_biomarker_column()`). The predicate that makes
   §5 meaningful. It is a regex over column names with no ground-truth labels
   available. Its failure mode is silent.
3. **§4 — the anchoring result.** The only strongly *positive* finding. It is
   also the one most at risk of being a tautology: an anchor uses the survey's
   own regional totals, and the evaluation scores against the survey's own
   district totals. §4.4 argues it is not circular; that argument deserves
   scrutiny.
4. **§9 — the noise ceiling.** `r_share > 1` in a large fraction of cells,
   meaning measured skill exceeds the estimated ceiling. Either the ceiling is
   biased low or something is optimistic. This is unresolved.

Numbers below are labelled **[measured this session]** or **[carried forward]**.
Carried-forward numbers were produced earlier in the project and were *not*
re-verified during this session; treat them as second-hand.

---

## 1. Study design and data

### 1.1 Surveys

Four sub-Saharan African national micronutrient/DHS surveys with individual-level
biomarker measurements, harmonised to a common schema:

| Country | Admin-2 units (analytic) | Survey year used for national covariates |
|---|---:|---|
| The Gambia | 30 | 2021 |
| Ghana | 75 | 2017 |
| Malawi | 87 | 2015 |
| Sierra Leone | 14 | 2013 |

Outcomes are binary deficiency indicators derived from biomarkers with
inflammation adjustment (BRINDA/Thurnham) where applicable, and IZiNCG
time-of-day cutoffs for zinc. The core four outcomes present in all countries are
`child_iron`, `child_vitA`, `women_iron`, `women_vitA`; Ghana, Malawi and Sierra
Leone add subsets of `child_zinc`, `women_zinc`, `women_folate`, `women_b12`.
This yields **24 country × outcome cells** for analyses using all outcomes, and
**16 cells** for analyses restricted to the four shared outcomes.

### 1.2 Predictors

294 harmonised area-level (Admin-2) proxy predictors, spanning: earth
observation (GEE — NDVI, climate, soil, night lights, accessibility), malaria
(MAP), disease burden (IHME), food prices (WFP), crop mix (MapSPAM), helminths
(ESPEN), and DHS-derived district aggregates from *other* surveys.

**Primary models exclude all `gw_` survey variables** — the concurrent
questionnaire items from the biomarker survey itself — because a model that
requires a concurrent survey cannot be deployed where no survey exists. `gw_`
variables appear only in the anchor arm of §5, and only with the blood draw
removed.

### 1.3 Critical structural fact

**Proxy predictors are Admin-2-resolved, so the effective sample size for
learning the predictor→outcome map is the number of AREAS, not individuals.**
Depending on the analysis this is 14–87 units per country. [carried forward]

Worse, the areas are themselves coarse relative to the survey design: in Ghana,
**62 of 75 districts contain exactly one survey cluster** [carried forward]. A
model fitted at the district cannot see within-district variation in the
predictors because there generally is none.

---

## 2. Estimation and evaluation

### 2.1 Learners

- Individual-level: SuperLearner via `mlr3superlearner`, NNLS-stacked
  (`.awsl_stack` / `.awsl_predict` / `.awsl_screen` in `R/area_weighted_sl.R`).
- Area-level: the same stack fitted directly on areas, plus comparators —
  training-mean null, covariate-free spatial smoother, penalised regression
  (ridge/lasso), random forest, Fay–Herriot, BYM2, CORAL.

### 2.2 Cross-validation

- **LORO** (leave-one-region-out): within-country. Every Admin-2 unit is
  predicted by a model that never saw its own Admin-1 region. Region blocks are
  used rather than random folds because spatial autocorrelation makes random
  CV optimistic.
- **LOCO** (leave-one-country-out): cross-country transport. Every survey from
  the held-out country is removed.

### 2.3 Selection optimism — the reason folds are strict

Two otherwise-identical pipelines, both holding out whole regions, differing
only in whether the covariate prescreen saw the held-out regions, differ by a
mean of **+0.182 in rank correlation, positive in all 20 measurable cells**,
reaching +0.51. [carried forward]

**Implication for the reviewer:** any comparison in this project between a
number produced here and a number from the literature must account for whether
the literature number screened on all data. Most do.

### 2.4 The noise ceiling

For a district with a binomial-sampled prevalence, observed between-district
variance contains sampling noise. Subtracting it gives the reliability
`lambda`, and `r_max = sqrt(lambda)` bounds any predictor's achievable
correlation. Median `r_max` across the 24 cells is **0.098**, and 4 of 24 cells
have no detectable signal above noise. [carried forward]

**Known problem, unresolved:** see §9.

---

## 3. Established results carried into this session

These frame everything below. None were re-verified this session.

| Result | Value |
|---|---|
| Individual-level prediction aggregated to Admin-2 | median r 0.516, r_share 0.92 |
| Area-level SuperLearner vs national-mean null | MAE 9.31–9.84 vs 9.55 pp — does **not** beat null |
| Covariate-free spatial smoother | r 0.304–0.393, beats every covariate arm |
| Covariates beat geography at Admin-2 | 6 of 24 cells |
| Admin-1 covariates vs Admin-1 spatial | 0.437 vs 0.209 |
| WHO category accuracy, Admin-2 | 33% vs 32% null — no better than null |
| WHO category accuracy, Admin-1 | 56% vs 35% null |
| National prevalence | 24/24 inside survey 95% CI, mean abs error 0.96 pp |
| Predictors surviving FDR control | **0 of 294, in all 24 cells** |
| Penalised regression retained predictors | median 0; 13 of 18 cells retain none |
| SuperLearner beats best of 21 comparators | 0 of 16 LOCO holdouts |
| Positive control (EO → district education) | r 0.48–0.71 — linkage is sound |

The positive control matters: the pipeline **can** predict a district-level
quantity from earth observation when that quantity is well measured. Whatever is
wrong at Admin-2 is not the geospatial linkage.

---

## 4. FINDING: where a district estimate's *level* comes from

**[measured this session]** — `scripts/covariates/08_admin1_arms.R` →
`results/tables/admin1_arms.csv`

### 4.1 Design

A district model produces two separable things: a **pattern** (which districts
are worse) and a **level** (how high prevalence sits). Five arms were scored on
identical cells under identical LORO folds:

| Arm | Cells | Mean r | Median r | MAE (pp) | Mean \|bias\| (pp) |
|:---|---:|---:|---:|---:|---:|
| Admin-1 anchor (hard) | 24 | **0.413** | 0.405 | 8.85 | **1.59** |
| Admin-1 anchor (shrunk) | 24 | 0.318 | 0.325 | 9.36 | 2.76 |
| Fit at Admin-1, extrapolate to Admin-2 | 22 | 0.230 | 0.318 | 8.78 | 3.18 |
| National anchor | 24 | 0.170 | 0.184 | 11.98 | 5.85 |
| No anchor | 24 | 0.164 | 0.192 | 10.71 | 3.24 |

The hard Admin-1 anchor beats no anchor in **20 of 24 cells**.

Each anchor rescales predictions so their (population-weighted) mean within each
region equals the design-based survey estimate for that region. The rescaling is
a **one-parameter shift on the logit scale**, hence monotone, hence the district
*ranking* is identical in every arm by construction. Only the level moves.

### 4.2 The national anchor is worse than useless

Anchoring to a single national total adds +0.006 to r and makes mean absolute
bias **worse** (5.85 vs 3.24 pp). A single number displaces every district by
the same amount, including districts that were already correct.

This is the mechanism that recurs in §6.

### 4.3 Hard beats shrunk

Shrinking regional targets toward the national value — the standard hedge
against treating a thinly sampled region as an exact anchor — recovers only
about half the gain (0.318 vs 0.413). On these data the hard anchor is
preferable despite the theoretical risk.

### 4.4 Is this circular? — the argument, for the reviewer to attack

**The concern:** the anchor uses the survey's regional totals; the evaluation
scores against the survey's district totals. The anchor is therefore
constructed from an aggregate of the very thing being predicted.

**The argument that it is not circular:**

- The anchor supplies the region's *mean*. The evaluation asks how well
  individual *districts* are predicted. Within a region, the anchor supplies no
  information about which districts are above or below.
- The shift is monotone, so **rank correlation within a region is mathematically
  unchanged**. The r gain comes entirely from correcting between-region level
  offsets — regions were being systematically mis-levelled relative to each
  other, and the anchor fixes that.
- The counterfactual arm is present: the national anchor uses the same kind of
  information at coarser resolution and gains nothing. If the gain were pure
  information leakage from "seeing the truth", the national anchor would also
  gain.

**What remains genuinely arguable:** the r gain measures how much the *model's*
between-region level errors were suppressing an otherwise-real within-region
signal. A reviewer could reasonably hold that this makes r the wrong summary
statistic here, and that only MAE/bias should be reported. Note that on MAE the
hard anchor gains far less (10.71 → 8.85 pp) than on r (0.164 → 0.413).

### 4.5 Practical consequence

The survey supplies the level; the model supplies the within-region pattern. A
district map for a country **with** a survey may be read as prevalence. A
transported map for a country **without** one has no anchor available and should
be read as a ranking only. §6 tests whether a predicted level can substitute. It
cannot.

---

## 5. FINDING: the individual-level anchor — the ceiling is the target, not the predictors

**[measured this session]** — `scripts/covariates/18_individual_anchor.R` →
`results/tables/individual_anchor.csv`

### 5.1 The question

A district r of 0.15 is uninterpretable alone. Against a ceiling of 0.13 it is
most of what is attainable; against a questionnaire model reaching 0.60 it is a
poor showing. The comparison that settles it replaces the geospatial proxies
with the **household questionnaire administered to these same respondents** —
wealth, education, diet, WASH, recent illness, supplementation — while removing
**every blood-assay column**, so the arm answers: *what would a household survey
without a blood draw buy?*

Same learner, same rows, same LORO folds. ~1,345 proxy predictors vs ~1,548 with
the questionnaire added.

### 5.2 Result

| Unit | Arm | Cells | Mean r | Median r | MAE (pp) | Mean predictors |
|:---|:---|---:|---:|---:|---:|---:|
| District | proxy | 16 | 0.154 | 0.126 | 8.33 | 1345 |
| District | **questionnaire** | 16 | **0.228** | 0.264 | 7.92 | 1548 |
| Cluster | proxy | 16 | 0.146 | 0.126 | 12.79 | 1345 |
| Cluster | questionnaire | 16 | 0.229 | 0.212 | 9.27 | 1548 |

- Questionnaire better in **10 of 16** cells; mean gain **+0.075**.
- Clears r = 0.4 in **3 of 16** cells.
- In **Malawi it adds nothing**: gains of 0.000, 0.000, 0.002, 0.004 across its
  four outcomes.
- **Max r anywhere across all 64 rows: 0.544** — the leak check (see §7).

Per-cell district results, sorted by gain:

| Country | Outcome | r proxy | r quest | gain |
|:---|:---|---:|---:|---:|
| Ghana | women_iron | −0.104 | 0.530 | +0.634 |
| Sierra Leone | women_vitA | −0.426 | −0.216 | +0.210 |
| Sierra Leone | women_iron | 0.166 | 0.355 | +0.189 |
| Gambia | child_vitA | 0.136 | 0.304 | +0.168 |
| Gambia | child_iron | 0.107 | 0.269 | +0.162 |
| Ghana | women_vitA | 0.028 | 0.052 | +0.024 |
| Ghana | child_iron | 0.398 | 0.418 | +0.020 |
| Gambia | women_iron | 0.454 | 0.466 | +0.012 |
| Malawi | women_vitA | 0.060 | 0.064 | +0.004 |
| Malawi | child_vitA | 0.079 | 0.081 | +0.002 |
| Malawi | child_iron | 0.366 | 0.366 | 0.000 |
| Malawi | women_iron | 0.117 | 0.117 | 0.000 |
| Sierra Leone | child_vitA | 0.103 | 0.098 | −0.005 |
| Ghana | child_vitA | 0.288 | 0.260 | −0.028 |
| Gambia | women_vitA | 0.399 | 0.319 | −0.080 |
| Sierra Leone | child_iron | 0.286 | 0.171 | −0.115 |

### 5.3 Interpretation — why this is the most important result

A null result of this kind invites three explanations, all of which this
excludes:

- **"The geospatial data is badly linked."** The questionnaire arm contains no
  geospatial data at all and also fails.
- **"The inflammation adjustment is wrong."** Both arms use the identical
  outcome. An outcome-definition error would move both.
- **"You are overfitting the 294 proxies."** The questionnaire arm has *more*
  predictors and the same folds, and does only marginally better.

The signal is weak at this resolution **for any predictor set that can be
constructed**, and the geospatial proxies already recover most of the little
that exists. This is consistent with §3's positive control (EO predicts district
*education* at r 0.48–0.71): the pipeline works; the biomarker target at Admin-2
does not carry recoverable district-level signal.

### 5.4 A prediction this falsified

The script was written to test a specific mechanism: since effective n is the
number of *clusters* rather than individuals, linking at the survey cluster
should help, and should help **most where clusters most outnumber districts** —
Sierra Leone, 58 clusters vs 14 districts.

| Country | Districts | Clusters | Mean gain | Better |
|:---|---:|---:|---:|:---|
| Gambia | 30 | 68 | **+0.017** | 6/8 |
| Ghana | 72 | 86 | −0.003 | 3/8 |
| Malawi | 83 | 97 | −0.002 | 2/8 |
| Sierra Leone | 14 | 58 | **−0.025** | 3/8 |

No gain anywhere, and **Sierra Leone is the worst, not the best** — the exact
opposite of the prediction. Additional units do not help when each unit's
outcome rests on a handful of respondents: the units become noisier in the same
proportion as they become more numerous.

### 5.5 Weaknesses of this analysis a reviewer should press on

- **Lighter learner.** The anchor uses a single NNLS stack with a top-40
  prescreen, not the production 12-learner stack, because 16 cells × 4 arms
  would take hours. The claim rests on the *gap* between arms under a consistent
  learner, but a richer learner could plausibly exploit the questionnaire more.
  **This is untested and is the most likely way the conclusion is wrong.**
- **Aggregation.** Individual out-of-fold predictions are averaged to the unit
  and units with n < 5 are dropped. Different thresholds were not swept.
- **The questionnaire is not a *good* questionnaire for this purpose.** DHS
  items were not designed to predict micronutrient status. A purpose-built
  dietary-intake instrument might do better; this result does not bound that.
- **Median vs mean.** District median r rises 0.126 → 0.264 (a doubling) while
  the mean rises 0.154 → 0.228. The distribution is skewed by Ghana women_iron
  (+0.634). Conclusions should be robust to dropping that cell — this was not
  formally checked.

---

## 6. FINDING: a predicted national level cannot substitute for a survey

**[measured this session]** — `scripts/covariates/19_national_composition.R` →
`results/tables/national_composition.csv`, `national_vmnis_loco.csv`,
`national_vmnis_ceiling.csv`

### 6.1 Motivation

§4 says the level must come from a survey. A country without one needs it from
somewhere else. The natural candidate: fit a model where the unit of observation
is the *country-survey*, using the WHO Vitamin and Mineral Nutrition Information
System (VMNIS) joined to World Bank national indicators.

### 6.2 The VMNIS model itself is competent

Leave-one-country-out; every survey from the target country removed at every
date. Surveys are never dropped for missing covariates (a complete-case rule
over ~1,600 candidate columns would discard nearly all of them); instead
missingness indicators, kNN imputation with training-only donors, and the
prescreen all happen **inside the fold**.

| Panel | Model | Surveys | Countries | Mean prev (pp) | MAE (pp) | Pearson | Spearman |
|:---|:---|---:|---:|---:|---:|---:|---:|
| Vitamin A / preschool | rf | 108 | 69 | 24.9 | 11.75 | **0.655** | 0.625 |
| Vitamin A / preschool | ridge | 108 | 69 | 24.9 | 12.88 | 0.562 | 0.611 |
| Vitamin A / preschool | null | 108 | 69 | 24.9 | 16.23 | −0.637 | −0.836 |
| Folate / NPW | rf | 32 | 28 | 28.5 | 16.71 | 0.529 | 0.511 |
| Zinc / preschool | ridge | 41 | 30 | 31.7 | 14.32 | 0.409 | 0.361 |
| Vitamin B12 / NPW | ridge | 23 | 21 | 16.4 | 9.38 | −0.861 | −0.882 |
| Vitamin A / NPW | rf | 29 | 22 | 11.7 | 9.07 | −0.133 | −0.020 |

**Reviewer note on the null arm's correlation.** The null's Pearson r is
strongly negative (−0.637, −0.847, −0.881). This is an artefact, not a finding:
the null predicts the training-country mean, which varies only slightly across
LOCO folds, so its correlation is dominated by trivial fold-to-fold variation
and is **not interpretable**. Only the null's MAE should be read. The same
caution applies to the B12 and Vitamin A/NPW panels, where the fitted models are
also near-degenerate.

### 6.3 The ceiling, and a correction to the record

`national_noise_ceiling()` decomposes
`logit(prev) = country + method + year trend + residual + sampling`,
computing the sampling term directly from VMNIS's recorded sample sizes
(design effect 1.5) rather than inferring error from a grouping choice.

| Panel | Surveys | Countries | sd country | sd method | sd resid | sd sampling | r_max (report) | r_max (standardised) | best r | r_share |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Vitamin A / preschool | 528 | 70 | 1.411 | 0.564 | 0.000 | 0.816 | **0.818** | 0.866 | 0.655 | **0.80** |
| Zinc / preschool | 219 | 30 | 1.130 | 0.468 | 0.254 | 0.207 | 0.892 | 0.960 | 0.409 | 0.46 |
| Folate / NPW | 184 | 28 | 1.764 | 0.950 | 0.442 | 0.331 | 0.849 | 0.954 | 0.529 | 0.62 |
| Vitamin A / NPW | 93 | 22 | 1.232 | 1.996 | 0.000 | 0.409 | 0.517 | 0.949 | −0.008 | −0.02 |

Two ceilings had previously been quoted for the vitamin A panel — 0.66 and
0.505. **They differ only in whether repeat surveys were grouped by country-year
or by country-year-*method***, an arbitrary choice that decides what counts as
error. Estimating the components directly supersedes both.

**This corrects a claim made earlier in this project** that the VMNIS model was
"98% saturated" (r 0.655 against a ceiling of 0.66). It is at **80%** of a
ceiling of 0.818. There is real headroom.

**Reviewer flags on this table:**
- `sd_resid` hits the **boundary at 0.000** in two panels. A variance component
  pinned at zero means the `lmer` fit is degenerate there and `r_max` is not
  trustworthy for those rows.
- For Vitamin A / NPW, **method variance exceeds country variance** (1.996 vs
  1.232). That panel is measuring survey methodology more than country.

### 6.4 The composition fails, and the oracle proves it cannot be fixed

Each transported district map was shifted onto a national level and scored
against the survey's own district prevalences.

| Arm | Cells | MAE (pp) | Mean \|bias\| (pp) | Better than no level |
|:---|---:|---:|---:|:---|
| Transported, no level | 8 | **5.81** | 3.35 | — |
| \+ national null level | 8 | 8.22 | 4.69 | 1/8 |
| \+ VMNIS predicted level | 8 | 12.70 | 10.03 | 1/8 |
| \+ **true** national level (oracle) | 8 | 5.58 | **0.73** | 4/8 |

Predicted vs true national levels:

| Outcome | Country | True (pp) | VMNIS (pp) | Null (pp) | VMNIS err | Null err |
|:---|:---|---:|---:|---:|---:|---:|
| child_vitA | Gambia | 19.78 | 20.27 | 17.19 | 0.49 | 2.59 |
| child_vitA | Ghana | 14.93 | 20.77 | 17.45 | 5.84 | 2.53 |
| child_vitA | Sierra Leone | 12.19 | **40.90** | 17.45 | **28.72** | 5.27 |
| child_vitA | Malawi | 9.18 | **31.21** | 17.13 | **22.03** | 7.95 |
| women_vitA | Gambia | 2.50 | 6.44 | 6.51 | 3.94 | 4.01 |
| women_vitA | Ghana | 1.51 | 6.18 | 6.51 | 4.66 | 5.00 |
| women_vitA | Sierra Leone | 1.49 | 7.58 | 6.51 | 6.09 | 5.02 |
| women_vitA | Malawi | 1.32 | 7.24 | 6.51 | 5.93 | 5.20 |

**Two compounding failures:**

1. **The predicted level is on a different scale.** The no-covariate null beats
   the covariate model in 6 of 8 cells. This is not a weak model but a model of
   a *different quantity*: VMNIS records what national surveys report under
   their own assays and inflammation adjustments. This is the same cross-survey
   biomarker level offset already documented as defeating LOCO transport.

2. **Decisively: even a perfect national level buys 0.23 pp** (5.81 → 5.58),
   better in only 4 of 8 cells. **The ceiling on this entire approach is a
   quarter of a percentage point**, so no further investment in national
   modelling can raise it. This is §4.2's mechanism from the other side.

### 6.5 Scope limitation — stated because it materially weakens the claim

VMNIS carries **Folate, Vitamin A, Vitamin B12, Vitamin D, Zinc — and no iron or
ferritin panel**. The LOCO transport predictions cover iron and vitamin A. The
intersection is **vitamin A only**, so this is 2 outcomes × 4 countries = 8
cells. Vitamin A / preschool is also VMNIS's *strongest* panel.

**This is a best case, not an average.** The oracle argument (point 2) is
outcome-agnostic and does generalise; the specific VMNIS-level failure (point 1)
is measured on 8 cells and should not be over-read.

Also note: `pearson`/`spearman` are undefined in some cells because the
transported map is flat across a country's districts (Gambia women_vitA, and
Sierra Leone after a large shift saturates the logit). Correlations are averaged
over the cells that have one, with the count reported.

---

## 7. DEFECT FOUND AND FIXED: outcome leakage via un-guarded assay columns

### 7.1 History — three iterations, two of which produced published-looking numbers

| Attempt | Guard | Result |
|---|---|---|
| 1 | `Xvars_full` as-is (drops only the outcome's own columns) | **r = 0.973** |
| 2 | Local regex anchored `(^|_)` on analyte tokens | **r = 0.986** |
| 3 | Shared `is_biomarker_column()`, three classes | max r 0.544 |

Attempt 2 failed because the analyte tokens sit **after** a population prefix:
`(^|_)FER` does not match `gw_wFER`. Attempt 1 was a genuine production bug —
`Xvars_full` is documented as "everything the survey collected except the blood
draw" but carried `gw_cFER`, `gw_cSTFR`, `gw_cCRP`, `gw_cAGP` (Gambia 15
columns, Ghana 14, Sierra Leone 10, Malawi 0).

### 7.2 The decisive step: measure, don't guess

Rather than iterate on regexes, every surviving column was ranked by |r| against
the outcome. The top of the list was:

| Column | \|r\| | What it is |
|---|---:|---|
| `gw_cAnemiaYN` | 0.800 | Hb-derived anaemia status |
| `gw_cAnemiaCat` | 0.776 | Hb-derived anaemia category |
| `gw_bs2` | 0.686 | blood-sample field |
| `gw_bis` | 0.668 | blood-sample field |
| `gw_cID_NoAdj` | 0.571 | **the same iron-deficiency indicator, unadjusted** |

**None of these names an analyte.** No pattern over analyte names could ever
have caught them. Analyte names were never the whole leak.

### 7.3 The guard as it now stands

`is_biomarker_column()` in `R/data_prep.R` — one definition, shared by
`build_outcome_dataset()` and by every downstream script.

Five pattern classes:

1. **Derived tokens** (case-insensitive): `FerAdj`, `RBPAdj`, `Retinol`,
   `VitADef`, `VitAInsuff`, `VitADefic`, `FeDef`, `FolDef`, `B12Def`, `ZincDef`
2. **Free analyte names** (case-insensitive, match anywhere — no English word
   contains them): `crp|agp|stfr|ferritin|transferrin|haemoglob|hemoglob`
3. **Case-sensitive analyte codes** (upper-case after a population prefix):
   `FER|RBP|TFR|Ferr`; plus bounded lower-case `_(fer|rbp)([^a-z]|$)`; plus
   case-sensitive `Hb`
4. **Derived status and unadjusted twins**: `_[cw]An[a]?emia`, `NoAdj`,
   `_[cw]ID($|_|s[TR])`, `_[cw]VAI`
5. **The blood-sample block**: `^gw_(bs[0-9]|bis$|rpb[0-9])`

### 7.4 Verification performed

Against **all 4,906 unique columns** across the four surveys' `Xvars_full`:
32 columns caught by class 4–5 (all draw-derived), 54 total removed per cell.

Deliberately **preserved** (verified individually):

- Identifiers: `gw_HHID`, `gw_MotherID`, `gw_SAMPLE_ID`, `gw_WomanID`,
  `gw_hID`, `gw_indivID`, `gw_momID`, `gw_old_indivID`,
  `gw_hHRostCaretakerID02–19`
- Knowledge item: `gw_wHeardAnemia`
- Beliefs-about-fortified-food block: `gw_wFFAnemia`, `gw_wFFGivesBlood`,
  `gw_wFFImpHealth`, `gw_wFFMort`, `gw_wFFSchool`, `gw_wFFWork`
- Exposures: `gw_wIodSalt`, `gw_wFeSuppl`, `gw_wFeSupplTime`,
  `gw_wSupplastPregFeFA`, `gw_wVitASuppl`, `gw_wFeRichFood`, `gw_wFvoPreferBuy`

### 7.5 Residual risk — please attack this

- **No ground truth.** The merged data has **no variable labels** (`attr(x,
  "label")` is NULL for every column tested), so the guard is a regex over
  opaque names. `gw_bis`, `gw_rpb1`, `gw_bs2` were classified as blood-sample
  fields **by name pattern and correlation**, not by documentation.
- **The 0.544 maximum is reassuring but not proof.** A leaked column with a
  genuinely moderate correlation would not stand out.
- **Over-blocking is possible.** `gw_wFFAnemia` is *assumed* to be a
  fortified-food belief item on the strength of its `gw_wFF*` siblings. If it is
  actually a measured status, it should be blocked. Conversely, blocking
  `NoAdj` wholesale may remove legitimate unadjusted *exposures* if any exist.
- **Only four surveys were checked.** A fifth country (Tanzania is planned) will
  have its own naming conventions and this guard will need re-verification.

**Suggested reviewer action:** obtain the original Stata `.dta` files with
labels intact and re-derive the guard from labels rather than names. That is the
correct fix; the current one is a workaround.

---

## 8. DEFECT FOUND AND FIXED: the tenth pair-key consumer

Malawi has **six Admin-2 names that occur in more than one Admin-1 region**
(`Lake Chilwa` ×3, `Lake Malawi` ×8, `TA Lundu`, `TA Malemia`, `TA Ngabu`,
`TA Pemba` — GADM level 2 has 256 polygons but only 243 unique names). Joining
Admin-2 tables on the district **name** therefore fans rows.

Nine consumers were migrated to a pair key `(Admin1, Admin2)` earlier in the
project. A tenth — `scripts/covariates/08_admin1_arms.R`, which produces the
anchoring table in §4 — was found **only by checking the deployed dashboard
against the app's own prediction table** (90 districts vs 87). It also contained
`key[!duplicated(key$Admin2), ]`, which silently *discards* the second district
of each duplicate-named pair.

Three shared library functions matched population weights on the name alone and
are now pair-aware, falling back to the name when a caller's table lacks Admin1:
`benchmark_admin2_table()`, `benchmark_admin2_to_admin1()` (both
`R/benchmark_area.R`), `aggregate_covariates_to_admin1()`
(`R/downscale_admin1.R`).

**Effect on the §4 result — conclusion unchanged, slightly strengthened:**

| Arm | r before → after | \|bias\| pp before → after |
|:---|:---|:---|
| No anchor | 0.170 → **0.164** | 3.60 → 3.24 |
| National anchor | 0.176 → **0.170** | 5.47 → 5.85 |
| Admin-1 anchor (shrunk) | 0.306 → **0.318** | 2.41 → 2.76 |
| Admin-1 anchor (hard) | 0.406 → **0.413** | 1.49 → 1.59 |
| Fit at Admin-1 | 0.244 → **0.230** | 3.42 → 3.18 |

Still 20 of 24 cells.

**Reviewer note:** this was the *tenth* instance of one defect class. A tenth
occurrence after nine fixes strongly suggests the codebase should not be relying
on developers remembering to call `admin2_join_by()`. A structural fix — a typed
key class, or a lint/test that fails on any `by = "Admin2"` — is warranted.
**No such test currently exists.** This is an open recommendation, not
implemented.

---

## 9. UNRESOLVED: `r_share > 1`

Mean `r_share` (measured r ÷ estimated noise ceiling) across the anchoring arms:

| Arm | mean r_share |
|:---|---:|
| Admin-1 anchor (hard) | **2.06** |
| Admin-1 anchor (shrunk) | 1.64 |
| Fit at Admin-1 | 1.77 |
| National anchor | 1.34 |
| No anchor | 1.35 |

**Measured skill routinely exceeds the estimated ceiling**, by a factor of two
for the best arm. It was already known that half of all cells show r ≥ r_max.
Either:

- the ceiling is **biased low** — plausible: `r_max = sqrt(lambda)` with lambda
  estimated by subtracting binomial sampling variance from observed
  between-district variance is itself very noisy at n ≈ 6–36 per district, and
  cannot go above 1 while r can fluctuate; or
- something in the evaluation is **optimistic** — the anchoring arms in
  particular use survey-derived regional totals, which is exactly where an
  optimism would enter (see §4.4).

**This is not resolved and should not be presented as if it were.** The
manuscript currently reports the ceiling with the caveat that half of cells
exceed it; it does not explain why. A reviewer able to distinguish these two
explanations would add a great deal.

---

## 10. Corrections to previously reported numbers

Made during this session. Any earlier document quoting the left column is stale.

| Claim | Was | Now | Why |
|:---|:---|:---|:---|
| VMNIS ceiling / saturation | r_max 0.66, "98% saturated" | r_max **0.818**, **80%** | old ceiling depended on an arbitrary survey-grouping choice (§6.3) |
| Anchoring gain | 0.170 → 0.406, bias 3.60 → 1.49 | **0.164 → 0.413**, bias **3.24 → 1.59** | pair-key fix (§8) |
| National anchoring | "buys almost nothing" | adds nothing to r **and makes bias worse** (3.24 → 5.85) | full arm table (§4.2) |
| Method standardisation | "MAE 12.46 → 10.81, a win" | MAE improves but **Pearson drops 0.592 → 0.514**; no-op for zinc and folate | it is a trade, and standardised prevalence is not the reportable quantity |
| Individual anchor | r 0.973, then 0.986 | **0.228** | leakage, twice (§7) |
| Cluster linkage | predicted to help most in Sierra Leone | **hurts most** in Sierra Leone | §5.4 |

---

## 11. The interpretation as it now stands

A single coherent account, with each claim's support:

1. **Admin-2 is below the resolution these surveys can support.** Median district
   contributes 6–36 biomarker measurements; median `r_max` 0.098. [carried
   forward, but see §9]

2. **This is a property of the target, not the predictors.** The household
   questionnaire administered to the same people barely does better than
   geospatial proxies (0.154 → 0.228), and the pipeline demonstrably *can*
   predict a well-measured district quantity (education, r 0.48–0.71). [§5]

3. **Geography carries most of what transportable signal exists.** A
   covariate-free spatial smoother matches the full 294-predictor set; no
   predictor survives FDR; penalised regression retains a median of zero.
   [carried forward]

4. **The level and the pattern are different problems with different answers.**
   The covariates carry pattern, not level. [§4, §6]

5. **The level must be resolved regionally to be worth anything.** Regional
   anchoring is worth more than any estimator choice tested (0.164 → 0.413);
   national anchoring is worth nothing and costs bias; an *oracle* national
   level buys 0.23 pp. [§4, §6.4]

6. **Therefore:** national estimates are the defensible deliverable; district
   maps are defensible as *rankings*, and as prevalences only where a survey
   supplies regional anchors; a transported map for a country with no survey is
   a ranking and must be labelled as one.

---

## 12. What a reviewer should run

```bash
Rscript scripts/covariates/08_admin1_arms.R          # -> admin1_arms.csv (§4)
Rscript scripts/covariates/18_individual_anchor.R    # -> individual_anchor.csv (§5)
Rscript scripts/covariates/19_national_composition.R # -> national_composition.csv (§6)
```

Requires the built targets store `_targets_full`. Set `R_USER` and `HOME` to the
OneDrive Documents directory so `~/.rdhs.json` resolves.

### Key implementation files

| Concern | File |
|:---|:---|
| Assay guard (§7) | `R/data_prep.R` — `is_biomarker_column()` |
| Pair-key helpers (§8) | `R/admin2_key_hygiene.R` — `admin2_join_by()`, `can_pair_join()` |
| Anchoring/benchmarking (§4) | `R/benchmark_area.R`, `R/downscale_admin1.R` |
| VMNIS model (§6) | `R/national_vmnis.R`, `R/national_covariates.R` |
| Area comparison arms | `R/area_level_comparison.R` |
| SuperLearner stack | `R/area_weighted_sl.R` |
| Noise ceiling | `admin2_reliability()`; `national_noise_ceiling()` in `R/national_vmnis.R` |

### Specific things to check

- [ ] Re-derive the assay guard from Stata variable labels; compare to the
      name-based guard. Any disagreement is a finding. (§7.5)
- [ ] Re-run §5 with the full 12-learner production stack on a subset of cells.
      If the questionnaire arm improves materially, §5's conclusion weakens.
      (§5.5)
- [ ] Drop Ghana women_iron (+0.634) from §5 and confirm the conclusion holds.
- [ ] Resolve §9 — is the ceiling biased low, or is the evaluation optimistic?
- [ ] Check whether §4's r gain survives a within-region-only correlation
      analysis, which would settle §4.4 definitively.
- [ ] Add a test that fails on `by = "Admin2"` anywhere in the codebase. (§8)
- [ ] Confirm the `lmer` boundary fits (`sd_resid = 0.000`) in §6.3 and either
      refit or mark those ceilings unusable.

---

## 13. Open questions

1. Would a purpose-built dietary/intake instrument beat both arms in §5? Nothing
   here bounds that — DHS items were not designed for this.
2. Is there a resolution *between* Admin-1 and Admin-2 at which covariates
   remain informative and the target is adequately measured? The Admin-1 result
   (0.437 vs 0.209 spatial) suggests the crossover is somewhere in that range.
3. How few, and how thinly sampled, can the anchor regions be before §4's gain
   disappears? This is the question that determines survey design, and it is not
   answered.
4. Does the level offset in §6.4 have a correctable structure — e.g. is it
   predictable from the assay and adjustment metadata VMNIS records? Outcome
   standardisation was tested and traded MAE against correlation (§10), so the
   answer is "not by that route".
