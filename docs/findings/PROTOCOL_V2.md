# Protocol v2 — the five fixes, implemented and measured

1 September 2026. Code: `R/protocol_v2.R`, `scripts/protocol_v2/`.
Outputs: `results/tables/protocol_v2/`.

This implements the five changes recommended after the 2026-09-01 methods audit
and signal-probe battery (`docs/findings/TWO_READINGS_2026-09c.md`), and
measures what each one bought. Every comparison below varies **one** thing on
identical cells, folds and learner.

## Why a parallel layer and not an edit to the DAG

Nothing here is in the `targets` graph. `R/protocol_v2.R` adds only new
functions, so no cached target is invalidated. The upstream fix for `n_svy`
semantics is a short addition to `compute_svy_admin2()` — carry `n_eff`
alongside `n_svy` — but that edit invalidates every `svy_admin2_*` target and
everything downstream, costing a multi-hour rebuild. It is deliberately **not**
made; this layer recomputes the corrected quantities from `_targets_full`
instead. Move `effective_n_v2()` upstream when a rebuild is acceptable.

## Verification before use

The rebuilt district outcomes reproduce the pipeline's own `svy_prev`
**exactly in all 24 cells** (r = 1.000, maximum difference 0.000 pp).

That check earned its place. The first build used the configured `oc$binary`
column and matched only 19 of 24 cells, diverging by up to **50 pp** on the
iron cells. Root cause: `compute_svy_admin2()` overwrites the configured binary
with one re-derived from the adjusted continuous biomarker under a uniform
cross-country cutoff (`resolve_uniform_outcome()` — BRINDA for vitamin A, a WHO
threshold otherwise). The pipeline is right and the first build was wrong; it
now calls the same resolver.

---

## FIX 1 — replicated folds and effective sample size

### 1a. How much can one fold draw move a result?

`results/tables/protocol_v2/fold_draw_risk.csv`, 10 draws per cell × arm.

| Quantity | Value |
|:---|---:|
| Median SD of a cell's Spearman across 10 draws | **0.056** |
| Median range (max − min) across 10 draws | **0.174** |
| Maximum range observed | **0.705** |
| Cell × arm combinations with range > 0.2 | **82 of 216** |

A single-draw report is typically good to about ±0.09 and can be off by ±0.35.
This is the mechanism behind the published `median r 0.058`, and it is why no
scheme with a random component is reported from one draw here. The region
scheme is exhaustive leave-one-Admin-1-out and therefore carries no draw luck
at all — the preferred within-country scheme for that reason.

### 1b. Effective sample size, and precision-weighted scoring

Design effects are estimated once per country × outcome at the **national**
level, because a district-level design effect is not estimable: Ghana 64 of 75,
Malawi 74 of 87 and Gambia 17 of 30 districts contain a single PSU, so
`svy_prev_se` is degenerate (< 1e-10) for them.

**Measured design effects are 1.0–5.6, median ≈ 2.4** — not the 1.5 the project
assumes. Effective sample size is therefore 1.3–3.0× smaller than the `n_svy`
in use:

| Country | median `n_svy` | median `n_eff` | ratio | districts with `n_eff` < 5 |
|:---|---:|---:|---:|---:|
| Gambia | 23 | 9.5 | 2.4 | 18% |
| Ghana | 10 | 7.2 | 1.4 | 33% |
| Malawi | 9 | 3.0 | 3.0 | **77%** |
| Sierra Leone | 40.5 | 31.3 | 1.3 | 0% |

Scoring districts by `n_eff` rather than equally moves mean absolute error from
**12.37 to 11.64 pp**, lower in **210 of 234** cell × arm combinations. The
unweighted evaluation was charging models for districts that carry almost no
information.

### 1c. Weight diagnostic (measured, not silently changed)

Two flagged weight problems, quantified rather than assumed:

- **Ghana's configured `gw_sWeight` takes 3 distinct values, matching
  `gw_strata` exactly.** It is a stratum constant, not a survey weight, while
  `gw_PSU_weight` / `gw_PSUStrat_weight` carry 90 distinct values (one per
  cluster). Switching moves district prevalence by 0.26 pp on average but up to
  **4.8 pp** for a single district, and up to 22 pp using `gw_cwt`.
- **Gambia's unused `gw_c_blood_weight` is identical to the configured
  weight** for the child biomarker cells (r = 1.000, difference 0.000). That
  standing concern is resolved: it is a non-issue.
- Sierra Leone `child_vitA` is the one cell where an alternative weight
  diverges materially (r = 0.63, mean difference 3.8 pp).

Changing a weight changes every published estimate, so this layer reports the
size of the decision and leaves it to the PI.

---

## FIX 2 — model the biomarker, not the dichotomised indicator

Continuous target minus binary prevalence, same arm, same folds, same cells:

| Estimand | continuous | binary | gain | continuous better in |
|:---|---:|---:|---:|---:|
| A. in-fill | 0.259 | 0.150 | **+0.109** | 84 of 126 |
| B. region extrapolation | 0.142 | 0.029 | **+0.113** | 75 of 108 |
| C. country transport | 0.161 | 0.060 | **+0.102** | 28 of 43 |

The gain is the same size in all three estimands, which is what a real property
of the target looks like rather than a protocol artifact.

---

## FIX 3 — rank-normalise within country before pooling

`scripts/protocol_v2/04_fix3_normalisation.R`. Leave-one-country-out transport,
identical cells and learner, varying only how predictors are put on a common
scale. Mean Spearman over 22 held-out country × outcome cells:

**A correction to an earlier version of this table.** The first run set
`standardize = FALSE` in `cv.glmnet` for every scheme, which is right when the
caller has already scaled the columns but punishes the `raw_pooled` arm for its
units rather than judging its information. Re-run with glmnet standardising
internally, the elastic-net rows change substantially and the conclusion
narrows. The corrected table, mean Spearman over 22 held-out country × outcome
cells:

| Scaling | level / index | level / enet | prev / index | prev / enet |
|:---|---:|---:|---:|---:|
| `raw_pooled` | **0.014** | 0.212 | **0.008** | 0.083 |
| `z_pooled` | 0.135 | **0.246** | 0.045 | **0.118** |
| `rank_within` | **0.174** | 0.231 | **0.103** | 0.091 |

Paired across all cells and arms: `rank_within` minus `raw_pooled` is
**+0.082 mean gain, better in 49 of 80**.

The honest reading is narrower than "raw pooling destroys transport":

- For a **penalised regression**, which standardises internally, the scaling
  scheme barely matters (0.212 against 0.246 on the level target). The
  estimator absorbs the offsets on its own.
- For **any composite built by averaging** — the domain index, and every
  index-like score — raw pooling is catastrophic: **0.014 against 0.174**. An
  unweighted mean of columns with wildly different units is simply the
  largest-unit column wearing a domain's name.

So the audit's "22 percent of the pooled vocabulary carries >10× offsets,
uncentered" finding costs little if the only consumer is a standardising
learner, and costs almost everything if any averaging happens first. Since the
best-transporting arm in this project *is* an average of averages, the fix is
load-bearing here — but the general claim must be stated with its scope.

### A companion defect found while building this

Building domain scores **per country** — the obvious implementation — lets each
country learn its own PC1 sign orientation, so a domain score means the
opposite thing in two countries and transport collapses. That implementation
scored **−0.059**; building the scores once on the pooled matrix, with the
orientation learned from the training countries only, scores **+0.151** on the
same data. Any composite or index built separately per country carries this
hazard, and it is invisible unless a cross-country comparison is run.

---

## FIX 4 — domain scores instead of raw columns

| Comparison | estimand / target | A | B | gain | A better in |
|:---|:---|---:|---:|---:|---:|
| domain enet vs raw enet | in-fill / prev | 0.161 | 0.103 | +0.058 | 12 of 18 |
| domain enet vs raw enet | region / prev | 0.028 | −0.085 | +0.113 | 12 of 18 |
| domain enet vs raw enet | in-fill / level | 0.272 | 0.308 | −0.036 | 8 of 18 |
| domain enet vs raw enet | region / level | 0.159 | 0.194 | −0.035 | 11 of 18 |
| **zero-tuning index vs raw enet** | in-fill / prev | 0.285 | 0.103 | **+0.182** | **17 of 18** |
| **zero-tuning index vs raw enet** | region / prev | 0.271 | −0.085 | **+0.357** | **16 of 18** |
| **zero-tuning index vs raw enet** | in-fill / level | 0.361 | 0.308 | +0.053 | 13 of 18 |
| **zero-tuning index vs raw enet** | region / level | 0.346 | 0.194 | +0.151 | 14 of 18 |

The result is sharper than "fewer predictors is better". A **penalised fit** on
18 domain scores beats one on 373 raw columns on the binary target and loses
slightly on the continuous one. What wins consistently — in 16 of 18 and 17 of
18 cells, by up to +0.357 — is the **zero-tuning domain index**: weight each
domain by its training-fold association and add them up. Nothing is selected,
so there is no selection to overfit. At these sample sizes the tuning is the
problem, not the predictor count.

---

## FIX 5 — three estimands, information-matched baselines

The corrected leaderboard, median Spearman over 24 cells:

| Estimand / target | spatial + domain | spatial | domain index | raw enet |
|:---|---:|---:|---:|---:|
| A. in-fill / level | **0.457** | 0.453 | 0.337 | 0.352 |
| A. in-fill / prev | 0.326 | **0.337** | 0.257 | 0.116 |
| B. region / level | **0.433** | 0.433 | 0.302 | 0.327 |
| B. region / prev | 0.252 | **0.258** | 0.235 | −0.053 |
| C. transport / level | — | — | 0.151 (enet 0.176) | — |
| C. transport / prev | — | — | 0.078 (enet 0.040) | — |

Against the published `within_country_consolidated` medians of **0.058**
(prevalence) and **0.206** (level), the corrected protocol reports **0.235 to
0.337** and **0.302 to 0.457**. The prevalence figure rises roughly fourfold to
sixfold, and none of that comes from a better model — it comes from replicated
folds, the corrected target, and an honest baseline.

Three comparisons the old leaderboard could not make:

- **Covariates beat the honest covariate-free survey baseline.** Against the
  *jackknifed* regional mean — the version that does not read the held-out
  district's own respondents — covariates score **0.323 vs 0.255**, better in
  **24 of 36**. The withdrawn version of that baseline (r 0.516) is not an arm
  in any estimand here.
- **Covariates and the spatial smoother are tied when folds match**: 0.323 vs
  0.327, covariates ahead in 21 of 36. The published claim that the smoother
  beat every covariate arm compared a leave-one-district-out smoother with
  leave-one-region-out covariate arms.
- **The same arm scores very differently under different estimands** (domain
  enet: 0.216 in-fill vs 0.094 region extrapolation). Conflating them, as the
  single old leaderboard did, is not a small approximation.

---

## Follow-up: domain representation and whether LOCO transport works

`scripts/protocol_v2/05_domain_representation.R`. Seven ways of representing
each domain, all learned inside the fold (PCA rotations from training countries
only), scored on leave-one-country-out ranking over 22 cells.

| Representation | predictors | index, level | index, prev | enet, level |
|:---|---:|---:|---:|---:|
| `mean1` sign-aligned mean (the current score) | 17 | 0.151 | 0.078 | 0.193 |
| `sup1` training-correlation-weighted mean | 16 | 0.196 | 0.100 | 0.185 |
| `pc1` first PC per domain | 16 | 0.218 | 0.123 | 0.230 |
| `pc12` first two PCs | 32 | 0.244 | 0.140 | 0.232 |
| `pc123` first three PCs | 48 | 0.235 | 0.123 | 0.250 |
| `pcvar` PCs to 80% of domain variance | 100 | **0.255** | **0.142** | 0.221 |
| `raw` no domain structure | 358 | — | — | **0.256** |

**The single sign-aligned mean is the worst representation in every column.**
It is what the domain index has been using, and replacing it with several
principal components per domain raises transport by +0.06 to +0.10. `pcvar` is
best for the untuned index and gives the best top-quartile capture (0.455 on
level, 0.331 on prevalence); for the elastic net, the full raw matrix is
marginally best but `pc123` matches it with 48 predictors instead of 358.

Supervision does not help: `sup1` is no better than `pc1`, so the gain comes
from representing more of each domain's variance, not from pointing it at the
outcome.

**Does LOCO transport work now?** For ranking, yes. On the biomarker level it
reaches mean Spearman **0.25, positive in 18 or 19 of 22 cells**, with
top-quartile capture of **0.44 to 0.46 against a chance level of 0.25** — about
a 1.8× targeting lift into a country the model has never seen. On prevalence it
is weaker, 0.14 with capture 0.33. This is a working **ranking** model and is
not a working **level** model: outcomes are standardised within country before
pooling, so no prevalence is claimed and none should be reported.

## What this does not change

Within a surveyed country, **geography still does most of the work**: the
spatial smoother alone matches or beats every covariate arm, and adding
covariates to it buys almost nothing (0.457 vs 0.453 on the best target). The
covariates' distinct contribution is in the two places geography cannot go —
beating a survey-derived baseline that has to be jackknifed to be honest, and
transporting to a country with no survey at all, where no smoother can be
fitted. Transport at Admin-2 remains modest (0.15 to 0.18) and is a **ranking**
claim only: outcomes are within-country standardised before pooling because
biomarker levels carry large cross-survey offsets, so no level is claimed.

Sierra Leone (14 areas) remains uninformative, the malaria domain remains null,
and women's vitamin A at Admin-2 remains a near-degenerate target — exactly
zero in 47 to 87 percent of districts, with 76 of Malawi's 87 districts sitting
on the logit clamp.

---

---

## Reading the leaderboard: what each arm is, in plain terms

Written for a collaborator who is not going to read the code. Every arm below
is scored on exactly the same districts, the same folds and the same outcome,
so the only thing that differs between two rows is the sentence describing it.

### The question each estimand asks

**A. In-fill.** *We surveyed a country but not every district. Can we fill in
the gaps?* Districts are held out at random, so a held-out district usually has
surveyed neighbours in its own region. This is the easiest of the three and it
is the one most published work implicitly reports.

**B. Region extrapolation.** *Can we predict a region we did not visit at all?*
Whole Admin-1 regions are held out, so the model must reach somewhere it has no
nearby information about. Harder, and the honest test for a survey that skipped
regions.

**C. Transport.** *Can a model built in three countries rank districts in a
fourth it has never seen?* This is the deployment case for a country with no
survey. Only RANKING is scored here: biomarker levels differ so much between
surveys - different assays, different adjustments - that a transported
prevalence is not a quantity this design can validate. A ranked priority list
is.

### The arms, from least to most information

**null_train_mean** - "give every district the average of the districts we
trained on." It uses no covariates and no geography. It exists to answer *would
a model have been better than assuming everything is average?* Its correlation
is always about -0.5 and that number is meaningless: a near-constant prediction
has no ranking, so its r is dominated by trivial fold-to-fold wobble. **Read
its error (MAE), never its correlation.**

**region_mean_jk** - "give every district the average of the OTHER surveyed
districts in its own region." This is the covariate-free comparator that
actually matters: it is what a survey statistician would do without any model,
and it is the honest version of a baseline this project previously reported at
r 0.516. That earlier version let each district see its own survey response,
which no covariate model can do; jackknifing removes that and the number drops
to 0.076 nationally. Available only under in-fill, because it needs other
surveyed districts inside the held-out unit's own region.

**spatial** - "fit a smooth surface over the map and read off the value." A
generalised additive model on district centroid coordinates. No covariates at
all - just the observation that neighbouring districts resemble each other.
This is the arm to beat, because if geography alone does the job then the whole
covariate programme buys nothing. Note it CANNOT be used for transport: a
smoother needs observed outcomes nearby, and an unsurveyed country has none.

**domain_index** - the deliberately simple covariate model. Each of the 18
conceptual domains (soil, climate, agriculture, water and sanitation, and so
on) is summarised into a few numbers; each summary is weighted by how strongly
it tracks the outcome in the training data; the weighted values are added up.
**Nothing is selected and nothing is tuned**, so there is no opportunity to
overfit. It is the arm that has held up best across every version of this
analysis.

**domain_enet** - the same domain summaries, but handed to a penalised
regression (elastic net) that decides which to keep and how much to shrink
them. More flexible than the index, and at these sample sizes that flexibility
is usually a liability rather than an asset.

**raw_enet** - the same penalised regression given ALL 451 raw predictor
columns instead of the domain summaries. This is closest to what the project
used to do. It is the clearest demonstration of the small-sample problem: with
14 to 87 districts and 451 columns, the model spends its effort choosing
between predictors rather than learning from them, and under region
extrapolation on the prevalence target it goes NEGATIVE.

**spatial_plus_domain** - fit the smooth spatial surface first, then let the
covariates explain what is left over. This is the sensible combination for a
country that HAS a survey, and it is usually the best or joint-best arm in
estimands A and B.

### The two targets

**level** is the district mean of the biomarker concentration itself (log
ferritin, log RBP, and so on). **prev** is the share of people below the
deficiency cut-off. The level target carries roughly 0.10 more correlation
everywhere, because dichotomising throws away information - a district where
everyone sits just above the cut-off looks identical to one where everyone sits
far above it.

### The columns

**mean r / median r** - Spearman rank correlation between predicted and
observed district values, averaged over cells. Rank rather than Pearson because
the deliverable is a priority ordering.

**positive** - in how many of the country-outcome cells the correlation was
above zero. A mean of 0.28 built from 21 positive cells is a very different
claim from the same mean built from 12 positive and 12 negative.

**topk (top-quartile capture)** - of the districts that truly are in the worst
quarter, what share does the model also place in its worst quarter. **Chance is
0.25.** This is the most decision-relevant number in the table: 0.44 means that
targeting the model's worst quarter reaches about 1.75 times as much of the
real burden as picking districts at random.

## How to re-run

```bash
Rscript scripts/protocol_v2/01_build_targets_v2.R          # targets, deff, weights
Rscript scripts/protocol_v2/02_run_benchmarks_v2.R malawi  # one country shard
Rscript scripts/protocol_v2/02b_merge_and_loco.R           # merge + transport
Rscript scripts/protocol_v2/03_effect_of_each_fix.R        # isolate each fix
Rscript scripts/protocol_v2/04_fix3_normalisation.R        # fix 3 in isolation
Rscript scripts/protocol_v2/05_domain_representation.R     # PCs per domain, LOCO
Rscript scripts/protocol_v2/06_build_variable_sheet.R      # RA annotation sheet
```

Step 02 takes a country name as its only argument and writes a shard, so the
four countries run in parallel; with no argument it runs everything in one
process. `PROFILE=smoke` runs Ghana at 3 replications. `V2_REPS` overrides the
replication count (default 10). Seeds are fixed; re-running reproduces the
tables.
