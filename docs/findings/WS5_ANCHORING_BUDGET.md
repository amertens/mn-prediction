# WS5. The survey-design experiment

Open question 3 of `docs/SESSION_FINDINGS_FOR_REVIEW.md`: how few, and how thinly
sampled, can the anchor regions be before the Section 4 gain disappears?

WS2 changed the shape of the question. The Section 4 gain does not survive a
jackknife, and a covariate-free arm that assigns each district its region's
design-based survey estimate outperforms the anchored covariate model. The useful
design question is therefore **how thinly a region can be sampled before its own
mean stops being a good estimate**, which is what this workstream measures. The
anchor is a jackknife throughout, so a district never contributes to the number
used to predict it.

`results/tables/anchoring_design_curve.csv`, 25 replicates per setting, seed
20260908, scored against the full survey's district and region estimates.

## A correction made before the result was read

The first summary pooled on the **absolute number** of anchored regions. The four
countries hold 4, 6, 16 and 28 Admin-1 regions, so "2 regions anchored" was
mostly Sierra Leone and Gambia while "12 regions anchored" was Ghana and Malawi.
The pooled curve showed error **rising** with more regions anchored, which is a
country effect and not a budget effect. The summary now works within country and
on the region share.

## The answer

### Thinning clusters barely moves the Admin-2 map

Admin-2 mean absolute error, all regions anchored, by share of clusters retained:

| Country | 25% | 50% | 75% | 100% |
|:---|---:|---:|---:|---:|
| Gambia | 14.91 | 14.29 | 13.89 | 13.66 |
| Ghana | 12.23 | 11.79 | 11.57 | 11.21 |
| Malawi | 14.18 | 14.13 | 13.73 | 13.44 |
| Sierra Leone | 6.46 | 5.75 | 5.39 | 5.09 |

Keeping a quarter of the clusters instead of all of them costs a mean of
**1.054 pp** and a median of **0.825 pp** in Admin-2 MAE, worse in 22 of 24
cells. Three quarters of the survey buys about one percentage point of district
accuracy.

### Dropping whole regions barely moves it either

Admin-2 MAE, all clusters retained, by share of regions anchored:

| Country | ~25-33% | ~50% | ~67-75% | 100% |
|:---|---:|---:|---:|---:|
| Gambia | 15.20 | 14.96 | 14.68 | 13.66 |
| Ghana | 12.30 | 11.78 | 11.46 | 11.21 |
| Malawi | 13.22 | 13.21 | 13.27 | 13.44 |
| Sierra Leone | 5.09 | 5.09 | 5.02 | 5.09 |

Malawi and Sierra Leone are flat to within noise; Malawi is marginally **worse**
at full coverage. Ghana and Gambia gain about 1.1 and 1.5 pp.

### The Admin-1 estimate is where the budget actually bites

Admin-1 mean absolute error, all regions anchored:

| Share of clusters retained | Admin-1 MAE (pp) |
|---:|---:|
| 25% | 5.998 |
| 50% | 3.667 |
| 75% | 2.514 |
| 100% | 0.297 |

**A factor of twenty between the thinnest and the fullest design.** The regional
estimate, which WS2 identifies as the quantity that carries the information, is
highly sensitive to how many clusters a region keeps. The district map is not.

## Acceptance: the smallest budget within a stated tolerance

At **Admin-2**, with a tolerance of 1 pp against the full-anchor result, keeping
**25 percent of clusters in every region** qualifies: the mean penalty is 1.054
pp, with a median of 0.825 pp and a 10th-to-90th-percentile spread across
replicates that overlaps the full-anchor result in every country.

At **Admin-1**, with the same 1 pp tolerance, **no setting below full retention
qualifies.** The nearest is 75 percent of clusters at 2.514 pp, still more than
eight times the full-anchor 0.297 pp.

## What a survey planner should take from this

The two findings point the same way and they are the practical output of this
workstream. Spend the budget on **clusters within the regions you sample**, not
on covering more regions, because the regional estimate is what degrades and it
degrades steeply. And do not expect a bigger survey to buy a better district map:
Admin-2 error moves by about one percentage point across a fourfold change in
sample size, because that error is not budget-limited. It is limited by the fact
that a regional estimate is being applied to districts that differ within the
region, which is the same conclusion WS2 reaches from the other direction.

## Reviewer pass, statistical

**Leakage.** The jackknife is applied inside every replicate: a district's
regional mean excludes its own respondents. Without it this experiment would
reproduce the circularity WS2 found.

**Sampling unit.** Whole clusters are dropped, not individuals, because a survey
planner buys clusters. Thinning within every cluster would price a design that
cannot be fielded.

**Replication and seeds.** 25 replicates per setting, seed 20260908, with
100 to 600 settings contributing to each summary row depending on how many
countries reach that region share.

**Identical scoring.** Every replicate is scored against the same full-survey
district estimates, so all budgets are compared on one yardstick.

**A limit.** The yardstick is the full survey, not the truth. A budget that
reproduces the full survey exactly would score zero error while still carrying
the full survey's own sampling error, which WS1 measures as substantial at
Admin-2. The Admin-2 numbers here should be read as distance from the survey, not
distance from the population.

**A second limit.** The region-share grid is coarse and lands on different exact
shares per country, so the pooled curve remains partly country-weighted even
after the correction. The within-country tables above are the ones to read.

## Reviewer pass, reproducibility

**Targets graph.** No new target.

**Additive results.** `anchoring_design_curve.csv` and
`anchoring_design_summary.csv` are new files; nothing is overwritten.

**Stamp targets.** Reads only `_targets_full` targets. No untracked input.

**Paths and determinism.** `here()` throughout; seeded; re-running reproduces.

**Reuse.** The workstream specification suggested reusing the replicate structure
in `scripts/run_subsample.R`. That machinery samples for a different purpose,
scoring model arms under subsampling rather than scoring an anchor budget, and
its stratum vocabulary (`has_clusters`, `zero_clusters`) does not express a
region-share grid. A fresh loop with the same replicate discipline was written
instead, and the choice is recorded here rather than left implicit.

**Smoke profile.** `PROFILE=smoke` runs Ghana with 5 replicates. Full run over 24
cells and 36 settings took about ten minutes.
