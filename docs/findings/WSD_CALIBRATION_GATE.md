# WS-D. The calibration gate

`results/tables/calibration_gate_report.csv`, 24 cells, one row per
country-outcome.

## What the gate is for

Correlation is silent about level. A district model can rank districts
plausibly while placing the whole country tens of points from where the survey
puts it, and every metric this project reported before WS-D — Pearson r,
Spearman, even MAE against a noisy district yardstick — can look acceptable
while that is happening. The gate asks one question the correlations cannot:

> **Does the model's implied national prevalence agree with the survey's?**

It is a necessary condition, not a sufficient one. Passing it does not make a
cell's predictions good; failing it makes them unusable as prevalences whatever
their ranking looks like.

## The threshold

```
threshold = max(10 pp, 2 x CI half-width)
```

Two components, and both are needed. The 10-point floor keeps the gate from
passing a badly-wrong cell merely because its survey estimate is imprecise. The
CI-scaled term keeps it from failing a cell whose disagreement is inside what
sampling noise alone would produce.

## The result

**18 of 24 cells pass. 6 fail.**

| Country | Outcome | Predicted | Survey | Gap |
|:---|:---|---:|---:|---:|
| Sierra Leone | child_vitA | 89.6% | 12.0% | **+77.6 pp** |
| Ghana | women_folate | 29.2% | 54.6% | **−25.4 pp** |
| Malawi | child_iron | 25.6% | 10.3% | +15.3 pp |
| Gambia | child_iron | 50.6% | 37.5% | +13.0 pp |
| Sierra Leone | women_folate | 91.7% | 79.2% | +12.5 pp |
| Gambia | child_vitA | 6.7% | 17.3% | −10.5 pp |

Sierra Leone child_vitA is the one to quote. The model predicts **89.6 percent**
of children vitamin-A deficient where the survey measures **12.0 percent** — a
77-point error, in a cell that would have been reported without comment on the
strength of its correlation alone.

## What the failures have in common

Four of the six are Sierra Leone or Gambia cells, and every failure clears the
10 pp floor rather than the CI-scaled term -- these are large absolute errors,
not marginal ones. This is the same effective-n story that
runs through the project: 14 districts across 4 regions leaves a
leave-one-region-out fit trained on about ten points, and an unconstrained ridge
on ten points can land anywhere.

The failures are also **not** confined to cells with weak correlation. That is
the point of having the gate at all: the two quantities are close to
independent, and a project that reports only one of them will ship the other's
failures.

## How it is used downstream

Failing cells are **marked, not deleted**. Every row of
`results/deliverables/district_estimates.csv` carries `calibration_status`, so a
consumer can see why a cell is absent from a map rather than finding a silent
hole. Six cells carry `calibration_failed`; cells the gate did not assess carry
`not_assessed` rather than being defaulted to a pass.

This composes with, and is separate from, the reliability suppression at
`r_max < 0.30` that withholds five Admin-2 layers. A cell can fail either, both,
or neither, and the two reasons are recorded independently.

## Reviewer pass

**The gate is applied to the un-anchored ridge**, which is the arm whose level
is at risk. The shipped empirical Bayes blend is anchored to the survey's own
district estimates by construction, so its national total is close to the
survey's by design — running the gate on it would be close to tautological and
is not evidence of anything.

**Country-key matching is on a squashed key**, so `SierraLeone` and
`Sierra Leone` resolve to the same cell. An earlier version matched on the raw
string and silently assessed nothing for that country.

**Additive.** `calibration_gate_report.csv` is new. Six tests cover the gate,
including one that a cell whose CI is very wide still fails on a large absolute
gap, which the CI-scaled term alone would let through.
