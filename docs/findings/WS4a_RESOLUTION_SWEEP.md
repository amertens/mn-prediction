# WS4a. The resolution sweep

Open question 2 of `docs/SESSION_FINDINGS_FOR_REVIEW.md`: is there a resolution
*between* Admin-1 and Admin-2 at which covariates remain informative and the
target is adequately measured? Section 13 suggests the crossover is somewhere in
that range, on the strength of the Admin-1 result (0.437 against 0.209 spatial).

`results/tables/resolution_sweep.csv`, 24 cells at four levels, 96 rows, seed
20260911.

## The answer: no, and the direction is monotone

**Skill peaks at Admin-1 and falls with every step toward Admin-2. There is no
intermediate optimum.**

Paired within cell, over the **14 cells present at all four levels**
(source: `results/tables/resolution_sweep.csv`, column `r_oof`):

| Level | Mean r | Median r | Median measurements per unit |
|:---|---:|---:|---:|
| Admin-1 | **0.315** | **0.228** | 33.0 |
| Admin-1 split into 2 | 0.166 | 0.203 | 19.5 |
| Admin-1 split into 3 | 0.093 | 0.112 | 16.0 |
| Admin-2 | 0.086 | 0.076 | 10.0 |

- **Admin-1 is the best level in 9 of 14 cells** and beats Admin-2 in **13 of
  14**.
- An intermediate level beats **both** Admin-1 and Admin-2 in only **4 of 14**
  cells, and never by a margin that survives the cell-to-cell scatter.
- Mean r falls by more than half from Admin-1 to the two-way split, and by
  roughly three quarters by Admin-2.

Across all available cells rather than the paired subset, the same ordering
holds on error: median MAE **5.32 pp** at Admin-1, **6.64** at the two-way split,
**7.32** at the three-way split, **9.66** at Admin-2.

## The ceiling falls too, which is why r_share is not the right summary

Median empirical `r_max` by level: **0.804**, **0.726**, **0.734**, **0.605**.
The ceiling degrades with resolution, as it must, since each unit is measured on
fewer people. Median `r_share` across levels is **0.270, 0.446, 0.260, 0.338**,
which is non-monotone and should not be read as a resolution signal: it is a
ratio of two quantities that both fall, estimated at each level on a different
number of units, and at Admin-1 the ceiling itself rests on as few as 16 units.

One cell shows what that does. Ghana `women_vitA` at Admin-1 returns `r_share`
**1.873**, from an r of 0.743 against an `r_max` of 0.323 measured on 16 units.
The raw correlation is the reliable half of that pair. **The raw-r comparison is
the finding; the r_share comparison is reported for completeness and is noisier.**

## The constraint this surfaced, which is arguably the more useful result

**Two of the four countries cannot support a within-country covariate model at
Admin-1 at all.**

| Country | Admin-1 units | Cells with a fitted model at Admin-1 |
|:---|---:|---:|
| Gambia | 6 | **0 of 4** |
| Ghana | 16 | 6 of 6 |
| Malawi | 27 | 8 of 8 |
| Sierra Leone | 4 | **0 of 6** |

Under leave-one-region-out blocking, Sierra Leone would train a 20-covariate
ridge on **three** units and Gambia on **five**. Both fall below the eight-unit
floor this script sets, so they return no estimate rather than a meaningless one.

This is a structural tension rather than a threshold artefact. The resolution at
which the survey measures the target well is the resolution at which there are
too few units to learn a covariate map. Where you can trust the estimate you
cannot fit the model, and where you can fit the model you cannot trust the
estimate. It is the same conclusion WS2 reaches from the other direction, where
a covariate-free regional mean outperforms every covariate arm.

## Measurements per unit: where the crossover sits

Expressed as the survey-design quantity, over the paired cells: skill is highest
at a median of **33 measurements per unit** and has fallen by half by **19.5**.
Across all cells the Admin-1 median is **54.25** and the Admin-2 median is
**13.0**.

The practical reading is that below roughly **20 measurements per unit** these
covariates carry very little, and the district level, at a median of 10 to 13,
is well inside that region. This is consistent with the direct measurement in
`docs/FITNESS_FOR_USE.md` section 2 and gives it a mechanism rather than only a
boundary.

The minimum column is a caution: some units hold as few as **1** measurement at
the finer levels, so the medians above conceal a long left tail.

## Method, and its choices

**The split rule.** k-means on district centroids within each region, k = 2 or 3,
25 restarts, seeded. This yields spatially **compact** parts. It does not
guarantee graph-theoretic **contiguity**: a region shaped around a bay can
produce a part whose members are close in Euclidean distance without sharing a
border. The parts are described as compact throughout and the workstream
specification's word "contiguous" is not claimed.

**Centroid imprecision, confined to Malawi.** `data/admin2_centroids.rds` is
keyed on country and district name with no Admin1. Malawi has **6** names
occurring in more than one region; their centroids are averaged before use, so
those districts share a coordinate and may be assigned to the same part when
they should not be. Affects Malawi only, and the count is carried in the output
as `n_shared_centroid_districts`.

**What is held constant.** The learner (`.ds_fit` ridge, top-20 in-fold
prescreen), the fold blocking (always on Admin1, so a held-out region is held out
whole at every level), the covariate set, and the aggregation rule
(population-weighted mean of district covariates). Only the unit changes.

**The ceiling is re-measured at every level** using the WS1a split-half
estimator, rather than assumed or carried from Admin-2. Its `min_units` floor is
lowered from 5 to 4 for this sweep so Sierra Leone's four regions can be
estimated at all; that change is a new argument with the old default, so WS1a's
published values are unaffected.

## Reviewer pass, statistical

**Paired comparison.** The headline uses the 14 cells present at all four levels,
so a level is never advantaged by being computable in easier cells. The
all-cells medians are reported beside it and agree in direction. The gap between
the two cell counts is itself informative: Admin-1 is computable in 14 cells and
Admin-2 in 24.

**Fold blocking is constant, and this matters.** Blocking on Admin1 at every
level means that at Admin-1 leave-one-region-out is leave-one-**unit**-out, which
is the thinnest possible training set. If anything this handicaps the Admin-1
arm, so the finding that Admin-1 nonetheless wins is conservative.

**Prescreen placement.** Inside the fold at every level.

**Seeds.** 20260911, covering both the k-means restarts and the split-half
resampling.

**What this does not test.** A different learner. The ridge with a top-20
prescreen is held constant so that resolution is the only thing varying, which
means the sweep bounds the resolution effect and not the learner effect. A
learner that exploits many units well might reorder the levels, though nothing in
WS3's protocol work suggests the learner is where the leverage is.

**A second limit.** At Admin-1 only Ghana and Malawi contribute, so the
Admin-1 column is two countries and the Admin-2 column is four. The paired
subset controls for this by construction, but the Admin-1 mean rests on 14 cells
from 2 countries.

## Reviewer pass, reproducibility

**Targets graph.** No new target; the script writes to `results/tables/` outside
the graph, following the `scripts/covariates/` convention.

**Additive results.** `resolution_sweep.csv` and
`results/figures/resolution_sweep.png` are new files. Nothing is overwritten;
the WS9 regression gate confirms the 19 frozen tables remain byte-identical.

**Shared-function change.** `split_half_reliability()` gains a `min_units`
argument defaulting to 5, the value WS1a used, so WS1a's committed output is
unchanged. Verified: the test suite passes 37 of 37 after the change.

**Joins.** District key built through `admin2_join_by()`, so Malawi is keyed on
the pair. The centroid table has no Admin1 and is therefore deduplicated by name
before use, which is the imprecision recorded above.

**Stamp targets.** Reads `data/admin2_centroids.rds` and
`data/covariates/harmonized/predictors_admin2_harmonized.csv`, both tracked, plus
the `_targets_full` store.

**Runtime.** About 50 minutes for 96 fits on a machine also running two of the
repository owner's jobs. `PROFILE=smoke` runs Ghana alone in about five minutes.
