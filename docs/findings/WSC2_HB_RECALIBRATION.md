# WS-C2. Recalibrating the field-haemoglobin arm

**Status: `stopped_partial`.** Six of roughly twelve cells completed before the
run was halted. The estimate remaining at that point was **two days, ten hours**,
against a machine also running the tournament and the repository owner's jobs.
The six completed cells are reported here; the run was not carried to
completion and no pooled figure is claimed from it.

## What was being tested

WS-C1 found field haemoglobin the strongest arm in the project for **ranking**
districts (mean gain +0.174, concentrated in iron outcomes) and the worst for
**level** (district MAE 11.07 pp against 8.33 for the proxy arm). The natural
repair is to recalibrate the haemoglobin arm's predictions against the survey's
own level before scoring, and ask whether the ranking gain survives.

## The six completed cells

| Country | Outcome | r before | r after | MAE before | MAE after |
|:---|:---|---:|---:|---:|---:|
| Gambia | child_iron | 0.865 | **0.794** | 6.45 | 7.45 |
| Gambia | child_vitA | 0.304 | **0.287** | 11.96 | 10.85 |
| Gambia | women_iron | 0.848 | **0.746** | 6.39 | 8.66 |
| Gambia | women_vitA | 0.457 | **0.106** | 3.20 | 3.65 |
| Ghana | child_iron | 0.425 | **0.313** | 12.83 | 13.77 |
| Ghana | child_vitA | 0.260 | **0.200** | 9.21 | 9.15 |

**Correlation falls in six of six cells.** Error rises in four of six. The
largest loss is Gambia women_vitA, from 0.457 to 0.106.

## Reading this honestly

The direction is consistent and the mechanism is plausible: recalibrating to the
national level removes between-district variation that was carrying real
ranking signal, so the arm gives up the thing it was good at in exchange for the
thing it was bad at. On these six cells that is a bad trade.

But **six cells is not twelve**, four of the six are one country, and Gambia is
the country where covariate signal is strongest throughout this project. A
pooled claim would be unsupportable. What can be said is narrower and still
useful:

> In every cell that completed, recalibrating the haemoglobin arm reduced its
> district-ranking correlation, and in two thirds of them it also increased
> error. There is no cell in which recalibration helped on both metrics.

The absence of a single cell where it helps on both is the part that does not
depend on the sample size.

## Why the run was stopped rather than resumed

The remaining cells are the large ones — Malawi at 87 districts and Sierra
Leone's six outcomes — and the per-cell cost is dominated by a SuperLearner
refit inside each fold. Two and a half days of exclusive machine time to convert
a consistent six-cell direction into a pooled twelve-cell mean was not a good
use of the budget while the estimator tournament, which decides what ships, was
competing for the same cores.

If it is resumed, the cheap version is to fix the learner to the ridge rather
than restacking per fold, which is what the tournament does and which costs
roughly an order of magnitude less.

## What it does not change

WS-C1's finding stands as reported: field haemoglobin buys ranking and costs
level, and it is an upper bound on what a survey that already drew blood could
recover rather than a deployable prediction. WS-C2 was asking whether that
trade could be improved, and on the evidence available the answer is no.
