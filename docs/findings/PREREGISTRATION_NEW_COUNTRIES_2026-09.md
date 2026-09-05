# Pre-registered predictions for the next countries (Ethiopia, Pakistan)

Written 3 September 2026, before any Ethiopian or Pakistani biomarker data has
been seen. The point of writing these down now is that several of the
project's strongest-looking results were *identified* on the four countries
already in hand, so quoting their numbers as findings would be selection on
the test. Stated as predictions in advance, they become tests instead. Each
prediction names the script that will score it, the metric, and the threshold
that counts as confirmation.

Source experiments: `docs/findings/SANDBOX_LOG_2026-09.md` (entries DA-01,
DA-02, DA-03, TC-01/02, G4-02, CF-01, CE-01) and `PROTOCOL_V2.md`.

## Protocol that applies to every prediction

- Outcome: uniform definition (`resolve_uniform_outcome()`), z-scored within
  country for transport; effective n from the measured design effect.
- Predictors: rank-normalised within country; domain PCs to 80% variance,
  oriented from the training countries only (`sign_rows = training rows`).
  When aggregating to regions, predictors are area-averaged (simple mean of
  districts), not population-weighted — population-weighting predictors cost
  0.05–0.15 Spearman on the current four (AG-01); outcomes may be weighted by
  effective n or population, which makes no difference.
- Estimand: leave-one-country-out, the new country held out; district rung
  is the country's own second sub-national tier (Malawi at its district rung,
  not Traditional Authorities); regional tier is the first sub-national tier.
- Metric: Spearman correlation between predicted and measured ordering,
  reported per outcome; a cell is "positive" if Spearman > 0.
- Scripts: `scripts/protocol_v2/02b_merge_and_loco.R`,
  `15_training_country_curve.R`, `16_admin1_transport.R`,
  `25_nested_domain_selection.R`, `30_training_curve_climate_soil.R`.

## Predictions

**P1 — Rankings transport; levels do not.** For a new country held out of
training, the full domain index ranks its districts with Spearman > 0 in at
least 75% of outcomes, and its regions with Spearman > 0 in every outcome
with ≥ 8 regions. Predicted national prevalence will differ from the survey's
by more than the survey's own confidence interval in most outcomes.
Primary criterion: the new country's mean Spearman across outcomes exceeds
the 95th percentile of a country-block permutation null (script 33) at each
tier; the cell counts above are secondary. *Basis:* mean ρ 0.29 (regional)
and 0.25 (district) against null 95th percentiles of 0.16 and 0.08 on the
current four (NC-01); counts of positive cells have a fat null when outcomes
within a country are correlated.

**P2 — A climate + soil index transports at least as well as the full
index.** On the same folds, the two-domain index (Climate and weather; Soil
characteristics) achieves district-rung Spearman within 0.03 of, or above,
the full 18-domain index, with positive transport in ≥ 90% of outcomes.
*Basis:* 0.368 vs 0.252 on level (22/22 positive) and 0.268 vs 0.151 on
prevalence on the current four — but that set was chosen on these countries,
which is exactly why this is a prediction and not a result. At the regional
tier the two are expected to be indistinguishable (DA-03: a wash).
*Amended 2026-09-04:* the DA-03 "wash" was computed on three countries
(BUG-01, `SANDBOX_LOG_2026-09.md`); on all 22 cells the climate + soil index
transports better than the full index at the regional tier too (0.45 vs 0.29
on the level, 0.38 vs 0.32 on prevalence). The prediction is therefore the
same at both tiers: climate + soil at least as good as the full index.

**P3 — Each added training country buys transported accuracy.** Adding the
new country to the training pool raises leave-one-country-out Spearman for
the *other* countries by between +0.02 and +0.06 per country for the full
index and by a smaller, positive amount (+0.01 to +0.04) for the climate+soil
index. *Basis:* monotone curves 1→2→3 countries, slopes +0.052/+0.056 (full)
and +0.038/+0.030 (climate+soil), across all four specifications.

**P4 — Survey-derived domains will not help transport.** Dropping household
assets, built environment, healthcare access, fertility and child morbidity
from the index will not reduce transported Spearman to the new country by
more than 0.01, and may raise it. *Basis:* all five have negative or zero
drop-one deltas on the current four (DA-01).

**P5 — The burden-capture margin over regional survey averages will be
small.** For districts a survey did not reach, the model's worst-ranked fifth
will capture 0.02–0.06 more of national burden than a jackknifed regional
mean, and will *not* beat a regional mean that includes the district's own
respondents. *Basis:* CF-01 — 0.240 vs 0.213 (jackknifed) vs 0.300 (with
self); model better than jackknife in only 8–11 of 18 cells.

**P6 — The model will not substitute for survey sample.** Fitted on a
reduced survey and scored against the full one, the model's district MAE
degrades within ±1 pp of the survey-only degradation across sample fractions
from 100% to 15%. *Basis:* G4-02, +3.36 vs +3.33 pp.

**P7 — The reliability ceiling will be lower on multi-cluster units.** Where
the new country's districts contain ≥ 2 survey clusters, the cluster-split
ceiling averaged over outcomes will be at least 0.05 below the within-split
ceiling (i.e. ≥ 10% of the within ceiling is cluster effect). No prediction
is made about which outcome shows it most: on the current four the iron
outcomes collapse in Gambia and Sierra Leone but not in Ghana or Malawi, and
child vitamin A shows no reliable geography on multi-cluster units in three
countries. *Basis:* CE-01 — pooled within 0.546 vs cluster 0.456.

## What would falsify the project's central claim

If P1 fails — district Spearman ≤ 0 in half or more of the new country's
outcomes — then the transport result on the current four countries was a
property of West/Southern Africa rather than of the proxies, and the
deliverable for unsurveyed countries has to be withdrawn. P2–P4 failing would
change the *recipe*, not the claim. P5–P7 are calibration of what the
project promises programmes, not of whether it works.

## Scoring

Run the scripts above with the new country added to `get_country_configs()`
and `targets_v2.csv` rebuilt by `01_build_targets_v2.R`; record outcomes in
`SANDBOX_LOG_2026-09.md` under a dated entry, one line per prediction:
confirmed / not confirmed / not testable (with the reason). Do not revise a
threshold after seeing the data.
