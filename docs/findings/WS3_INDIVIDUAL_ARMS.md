# WS3. Section 5 under one protocol

## Scope actually delivered

The full 16-cell re-run was started and abandoned under the run-or-reframe rule.
It progressed at roughly six of 192 fits per half hour, because the machine was
running the repository owner's own R jobs (`06_tmle_full_cohort_effects.R` and
`21_gose_specifications.R`) throughout. What ran is the **four-cell protocol
parity subset** the workstream specifies for WS3a, spanning all four countries
and both nutrients: Ghana `child_iron`, Gambia `women_iron`, Malawi `child_vitA`,
Sierra Leone `women_iron`
(source: `results/tables/individual_arms_2026-09_PARITY.csv`).

WS3d, the permutation importance for Ghana `women_iron`, is **not yet computed**.
WS3b and WS3c are computed on the published table's own arithmetic, which is
stated as such below.

## WS3a. The two individual-level numbers differ by fold construction

Stage 0 established the mechanism by reading code. `aggregate_admin2_sl()`
aggregates `res$yhat_full`, which `src/analysis/sl_helpers.R` produces via
`origami::make_folds(cluster_ids = id_vec, V = folds)`. That is a
**cluster-blocked K-fold**, not the region-blocked leave-one-region-out Section 5
uses. WS3 measures what the difference is worth.

District-level mean `r` over the same four cells, same learner, same rows:

| Arm | Region-blocked LORO | Cluster-blocked K-fold | Gap |
|:---|---:|---:|---:|
| proxy | 0.274 | 0.594 | **+0.319** |
| questionnaire | 0.261 | 0.609 | **+0.348** |
| questionnaire + field Hb | 0.342 | 0.628 | **+0.286** |

**Fold construction alone is worth about +0.32.** That accounts for most of the
distance between Section 3's carried-forward 0.516 and Section 5's 0.154, and it
does so without invoking any difference in learner, arm or aggregation.

The analysis of record for the individual-level result is therefore the
**region-blocked** row. Section 3's number should not be quoted as an
out-of-sample individual-level result without stating that its folds are blocked
on clusters, which leaves neighbouring clusters inside the same region in the
training set, and that its preprocessing recipe is fitted on all rows before
folds are formed.

The remaining cell of the 2 by 2, the production 12-learner stack under both
protocols, is **not yet computed**. What is bounded here is the **protocol**
effect, not the learner effect.

## WS3e. What the questionnaire buys, and what field haemoglobin buys

Under the strict protocol, with the guard corrected and the arms nested so that
proxy is a strict subset of questionnaire:

| Country | Outcome | proxy | questionnaire | + field Hb | gain, questionnaire | gain, +Hb |
|:---|:---|---:|---:|---:|---:|---:|
| Gambia | women_iron | 0.454 | 0.466 | **0.848** | +0.012 | **+0.394** |
| Ghana | child_iron | 0.398 | 0.410 | 0.425 | +0.012 | +0.027 |
| Sierra Leone | women_iron | 0.166 | 0.089 | 0.016 | -0.077 | -0.150 |
| Malawi | child_vitA | 0.079 | 0.079 | 0.079 | 0.000 | 0.000 |
| **Mean** | | **0.274** | **0.261** | **0.342** | **-0.013** | **+0.068** |

**The questionnaire gain reverses sign.** The published mean gain is **+0.075**
over 16 cells; measured here with the leak removed and the arms nested it is
**-0.013** over four cells, better in 2 of 4. Two changes produce this:

1. **The leak.** The published questionnaire arm saw `gw_wm_whbc` (women's
   haemoglobin) and `gw_gchb` (child haemoglobin), plus thirteen further
   blood-derived columns found from Stata labels (WS7a, WS7b).
2. **The non-nesting.** The published questionnaire arm applied
   `is_biomarker_column()` to `Xvars_full`, which also stripped the MAP
   sickle-cell rasters and the DHS mean-haemoglobin Admin-2 aggregates that the
   proxy arm keeps, because the proxy arm uses `Xvars` unfiltered. The two arms
   differed in more than the questionnaire.

Section 5's conclusion, that a household questionnaire administered to the same
people barely beats geospatial proxies, is **strengthened, not weakened**: on
these four cells it does not beat them at all.

**Field haemoglobin is a different matter, and it is outcome-specific.** It adds
+0.394 for Gambia women's iron, taking that cell to 0.848, and +0.027 for Ghana
child iron. It adds nothing for Malawi child vitamin A and costs 0.150 for
Sierra Leone women's iron. That pattern is what physiology predicts:
haemoglobin is informative about iron status and uninformative about vitamin A.

Two cautions on the haemoglobin arm. Its mean absolute error at district level
is **19.96 pp** against 10.12 pp for the proxy arm, so it gains on correlation
while losing badly on level. And the whole positive mean rests on one cell.

## WS3b. Malawi has no questionnaire block

Malawi's `GW` domain contains **zero** columns: `Xvars_full` equals `Xvars` at
1222 columns, against 378 `gw_` columns for Ghana (measured from the
`outcome_data_*` targets). Three independent confirmations:

- In the published table its questionnaire arm has **fewer** predictors than its
  proxy arm, 1141 against 1147, the difference being the `sd > 0` filter after
  imputation.
- The WS7a leakage report finds Malawi's proxy and questionnaire maxima
  **identical to four decimal places in all eight outcomes**.
- Here, all three arms return **exactly 0.0789** for Malawi `child_vitA`.

The published gains of 0.000, 0.000, 0.002 and 0.004 across Malawi's four
outcomes are the signature of one arm scored twice, not a finding about Malawi.

**Scoping the ingestion, as instructed.** The data exists: 242 questionnaire
columns in `data/IPD/Malawi/Malawi_merged_dataset.rds`, coded `m01`, `m115a`,
`m220h`. **Zero carry a variable label.** The only local documentation is a
349-character README that names four population files and points technical
questions to `immpact@cdc.gov`. A name-based guard is blind to `m115a`, and this
branch has now found the eleventh through twenty-fourth instances of that defect
class in surveys whose names were at least semi-informative. Ingesting 242
opaque codes into an arm whose entire claim is "no blood draw" is not defensible
without a codebook. **Held as a to-do**, with the specific request being the
Malawi MNS questionnaire codebook from CDC IMMPaCt, or the DHS merge the README
describes via `MCLUSTER`/`MNUMBER`/`M01`, which would supply labelled
equivalents.

## WS3c. Dropping Ghana women_iron

On the published table's own arithmetic
(source: `frozen_2026-09/individual_anchor.csv`, `unit == "district"`):

| Subset | Cells | Mean gain | Median gain | Questionnaire better |
|:---|---:|---:|---:|:---|
| As published | 16 | +0.0748 | +0.0080 | 10 of 16 |
| Excluding Malawi | 12 | +0.0993 | +0.0220 | 8 of 12 |
| **Excluding Ghana women_iron** | 15 | **+0.0375** | **+0.0040** | 9 of 15 |
| Excluding both | 11 | +0.0506 | +0.0200 | 7 of 11 |

The published +0.075 and 10 of 16 reproduce exactly. Dropping the one cell
Section 5.5 names halves the mean gain and takes the median to +0.004. Section
5.5 says this "was not formally checked"; it now is, and the conclusion holds.

## WS3f. Cluster against district

Under the strict protocol, mean gain from scoring at the cluster rather than the
district:

| Country | Cells | Mean gain | Better |
|:---|---:|---:|:---|
| Gambia | 3 | +0.018 | 2 of 3 |
| Ghana | 3 | -0.008 | 0 of 3 |
| Malawi | 3 | -0.001 | 0 of 3 |
| Sierra Leone | 3 | **-0.034** | 0 of 3 |

**Section 5.4's falsification is confirmed and its ordering is reproduced.** The
published values are Gambia +0.017, Ghana -0.003, Malawi -0.002, Sierra Leone
-0.025. Sierra Leone, where clusters outnumber districts by the widest margin
and where the mechanism predicted the largest gain, is again the worst.

Note this is not the cluster-train, district-evaluate comparison the workstream
describes. Both variants here are trained on individuals and differ only in the
unit of aggregation, which is what the published script does. Training at the
cluster and scoring after aggregation to the district is **not yet computed**.

## Reviewer pass, statistical

**Identical cells, rows and folds.** All three arms and both protocols run inside
one loop per cell from one outcome vector and one fold vector, so an arm
difference is a predictor-set difference and a protocol difference is a fold
difference. The cluster-blocked folds are drawn once per arm from the same seed
stream, which means the two protocols are not paired on identical fold
assignments across arms; the protocol gap of +0.32 is far larger than any
plausible fold-draw variation, but the design is worth stating.

**Prescreen placement.** Inside the fold, via `.awsl_screen` called on training
rows only.

**Leakage.** The arms use the corrected guard, which now blocks 90 of the 4,906
columns against 77 before this branch. The questionnaire arm's predictor count
rose relative to the published run (Ghana 2086 against 2081) because the nesting
fix restores external columns the old filter stripped, while the guard fix
removes survey columns it should have stripped.

**Seeds.** 20260906, covering the cluster-fold assignment. The NNLS stack is
seeded internally at 20260829 by `.awsl_stack`.

**Aggregation.** Unweighted district means with `n >= 5`, matching
`scripts/covariates/18_individual_anchor.R` so the comparison against the
published table is like for like. WS1e established that survey weights are
constant within district in three of the four countries, so the weighted and
unweighted aggregates coincide there.

**The main limit.** Four cells. Every aggregate above is a mean over four
numbers and two of the four gains are zero by construction (Malawi) or driven by
one cell (Gambia). The protocol gap is the robust part, because it is measured
within each cell and has the same sign in all four.

## Reviewer pass, reproducibility

**Targets graph.** No new target. `biomarker_column_class()` and
`allowed_under_arm()` are added to `R/data_prep.R`; `is_biomarker_column()` is
extended. The WS7a `leakage_report` target depends on `outcome_data_*` and so
re-runs over any new country's columns.

**Store staleness.** The corrected guard is applied at read time, exactly as
`scripts/covariates/18_individual_anchor.R` already does and documents, so this
re-run does not require rebuilding `_targets_full`. WS7a measured that the guard
change adds no newly outdated targets: 806 of 845 were outdated before and after,
under the store's own build settings.

**Additive results.** `individual_arms_2026-09_PARITY.csv` and
`..._SMOKE.csv` are new files. `individual_anchor.csv` is untouched and the WS9
regression gate confirms it.

**Smoke profile.** `PROFILE=smoke` runs Ghana `child_iron` across three arms and
both protocols in about 25 minutes on a contended machine.
`WS3_CELLS=parity` runs the four-cell subset.
