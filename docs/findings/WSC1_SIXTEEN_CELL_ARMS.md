# WS-C1. The individual-level arms at sixteen cells

`results/tables/individual_arms_2026-09_16CELL.csv`. Three nested arms, two fold
protocols, two aggregation units, all sixteen shared-outcome cells.

## This corrects a claim I made from four cells

The previous report stated that the questionnaire gain **reverses sign** to
**-0.013** once the arms are decontaminated and nested. That was measured on the
four-cell parity subset. On all twelve computable cells it is **+0.036**.

| Basis | Cells | Mean questionnaire gain |
|:---|---:|---:|
| Published (contaminated, non-nested, includes Malawi) | 16 | +0.0748 |
| Published, excluding Malawi's four structurally-empty cells | 12 | +0.0993 |
| **Corrected, decontaminated and nested** | **12** | **+0.0358** |
| Corrected, four-cell parity subset (superseded) | 4 | -0.0130 |

The like-for-like comparison is the second and third rows, both on the same
twelve cells: decontamination and nesting remove **roughly two thirds** of the
published gain. **It does not reverse the sign.** The four-cell subset was not
representative and the sign reversal it showed does not survive.

Section 5's substantive conclusion is unchanged and slightly softened: a
household questionnaire administered to the same people adds a little, and much
less than the published figure.

## The three arms under the strict protocol

Per cell, district level, region-blocked folds:

| Country | Outcome | proxy | questionnaire | + field Hb | gain quest | gain Hb |
|:---|:---|---:|---:|---:|---:|---:|
| Gambia | child_iron | 0.107 | 0.269 | **0.865** | +0.161 | **+0.758** |
| Gambia | child_vitA | 0.136 | 0.304 | 0.304 | +0.168 | +0.168 |
| Gambia | women_iron | 0.454 | 0.466 | **0.848** | +0.012 | +0.394 |
| Gambia | women_vitA | 0.399 | 0.437 | 0.457 | +0.039 | +0.058 |
| Ghana | child_iron | 0.398 | 0.410 | 0.425 | +0.012 | +0.027 |
| Ghana | child_vitA | 0.289 | 0.244 | 0.260 | -0.045 | -0.029 |
| Ghana | women_iron | -0.104 | -0.012 | **0.593** | +0.092 | **+0.697** |
| Ghana | women_vitA | 0.028 | 0.052 | 0.105 | +0.024 | +0.077 |
| Sierra Leone | child_iron | 0.286 | 0.163 | 0.163 | -0.122 | -0.122 |
| Sierra Leone | child_vitA | 0.103 | 0.096 | 0.099 | -0.007 | -0.004 |
| Sierra Leone | women_iron | 0.166 | 0.089 | 0.016 | -0.077 | -0.150 |
| Sierra Leone | women_vitA | -0.426 | -0.253 | -0.211 | +0.173 | +0.214 |
| **Mean over 12** | | | | | **+0.036** | **+0.174** |

## Field haemoglobin is the substantial arm, and it is iron-specific

Mean gain **+0.174** over twelve cells, and the three largest gains are all iron
outcomes: Gambia child iron **+0.758** (reaching 0.865), Ghana women iron
**+0.697** (reaching 0.593 from a negative baseline), Gambia women iron
**+0.394** (reaching 0.848). Vitamin A gains are small or zero. That is the
pattern physiology predicts, since haemoglobin responds to iron status and not
to vitamin A status.

**It buys ranking and costs level.** Under the strict protocol its district MAE
is **11.07 pp** against **8.33 pp** for the proxy arm. Reported as a ranking
instrument it is the strongest arm in the project; reported as a prevalence it
is the worst of the three.

**Sierra Leone is the exception in both arms.** Its questionnaire and
haemoglobin gains are negative in three of four cells. With fourteen districts
and four regions, a leave-one-region-out fold trains on about ten districts, and
adding predictors there costs more than it buys.

## The protocol gap, at sixteen cells

Proxy arm, district level: **0.416** under cluster-blocked K-fold against
**0.154** under region-blocked folds, a gap of **+0.262**. The four-cell subset
gave +0.319. Both confirm that fold construction accounts for most of the
distance between Section 3's carried-forward individual-level number and
Section 5's, without invoking any difference in learner or arm.

## Malawi is `not_computed`, not scored twice

Malawi's four cells report `status = not_computed` for both questionnaire arms,
with the reason recorded on each row. Its `GW` domain holds zero columns, so a
questionnaire arm there would be the proxy arm under another label, which is
what produced the published gains of 0.000 to 0.004. **The script now detects an
empty survey domain and refuses to score the arm** rather than emitting a
duplicate.

WS-C4, the automatic DHS crosswalk that would let Malawi enter, is
**not_computed**; see the WS-C4 note.

## Reviewer pass, statistical

**Shared-noise check.** The proxy arm sees no survey column, so zero. The
questionnaire and haemoglobin arms see columns measured on **the same
respondents** whose outcomes form the target, so the shared-noise fraction there
is **not zero**: an individual's questionnaire response and their biomarker come
from the same person in the same visit. This is not sampling noise shared with
the yardstick in the WS-B1 sense, because the target is the district aggregate
and the predictor enters through a model fitted out of region. It is worth
stating plainly that these two arms are not deployable predictions but upper
bounds on what a survey that already collected the data could recover.

**Nesting.** `allowed_under_arm()` scopes the filter to the concurrent survey,
so proxy is a strict subset of questionnaire, which is a subset of
questionnaire-plus-haemoglobin. Predictor counts confirm it: 1345, 1687, 1693 on
average.

**Identical cells and folds.** All arms and protocols run from one outcome
vector and one fold vector per cell. The gain columns are paired within cell.

**Seeds.** 20260906 for the cluster-fold assignment; the NNLS stack seeds
internally.

**Twelve cells, not sixteen, for every arm comparison.** Malawi contributes to
the proxy column and to nothing else, so the proxy mean of 0.154 is over sixteen
cells and the arm means of 0.189 and 0.327 are over twelve. Only the paired
per-cell gains should be read across arms.

## Reviewer pass, reproducibility

**Additive.** `individual_arms_2026-09_16CELL.csv` is new;
`..._PARITY.csv` and the frozen `..._4CELL` snapshot are untouched.

**Status column.** Every row carries `status` and `reason`, following the
`ws9_session_summary.R` convention, so a consumer can distinguish a measured
zero from an uncomputed one.

**Runtime.** About five hours for 192 rows on a machine also running the
repository owner's jobs.
