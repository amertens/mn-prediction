# XO-01 — Cross-outcome borrowing under transport (design)

Date: 2026-09-19. Sandbox ID XO-01. Status: approved, implementing.

## Question

Can a country that measured one biomarker (say iron) get a district ranking
for a biomarker it did not measure (say folate) by borrowing the
outcome-to-outcome relationship learned in countries that measured both?

## What the surveyed data already say (targets_v2.csv, district-level
Spearman of the continuous biomarker, within country)

| Pair | Gambia | Ghana | Malawi | S. Leone | pooled |
|---|---|---|---|---|---|
| child iron ↔ women iron | 0.74 | 0.38 | 0.31 | 0.44 | 0.41 |
| child vitA ↔ women vitA | 0.79 | 0.48 | 0.30 | 0.19 | 0.43 |
| child iron ↔ women folate | – | 0.28 | 0.08 | 0.16 | 0.17 |
| women iron ↔ women folate | – | 0.17 | 0.17 | −0.20 | 0.14 |
| women folate ↔ women B12 | – | +0.34 | −0.53 | −0.08 | −0.13 |

Same nutrient across populations is stable (~0.4); cross-nutrient is weak
and sign-unstable. The transported index already produces near enough one
map for every outcome (CIV rankings agree at 0.92 across outcomes).

## Pre-registered predictions

1. Cross-nutrient borrowing: mean delta vs base within ±0.05.
2. Same-nutrient, other population: delta +0.05 to +0.10.
3. Folate and B12 cells: no gain from any arm.

## Estimand

Transport only (estimand C, leave-one-country-out). The held-out country's
OTHER biomarkers are visible as predictors; its target outcome is not.

## Design

Script `scripts/protocol_v2/63_cross_outcome_borrowing.R`, standalone in
the style of `scripts/policy_deck/04_civ_climate_soil_prediction.R`. The
biomarker block is built in memory and never written to
`predictors_admin2_shared.csv`, so the headline vocabulary and the leakage
policy are untouched.

Parameters: `V2_DOMAIN_SET = cs | cs_top5 | all` (default `cs`, the deployed
transport candidate); `V2_PREDICTOR_TIERS` default `open,survey_public`.

Cells: target outcome Y in {child_iron, child_vitA, women_iron, women_vitA,
women_folate, women_b12} × scale in {level, prev}; held-out country = each
country with Y; needs ≥ 3 countries (zinc is Malawi-only and drops out).
The block for a cell is the set of outcomes measured in EVERY country of the
cell, minus Y (column-intersection rule, so the block is identical across
folds). Block columns are the other outcomes' `y_level`, rank-normalised
within country by `prep_predictors_v2()` like every other predictor, under a
new domain "Survey biomarkers (other outcomes)".

Arms (each through `domain_index`, with `domain_enet` alongside):

- `base` — covariates only.
- `same_nutrient` — base + the other population of the same nutrient
  (skipped where none exists: folate, B12).
- `other_nutrients` — base + every other-nutrient outcome, both populations.
- `all_other` — base + the whole block.

Scoring per held-out country with `score_v2()`; level metrics suppressed as
in `02b` (transport is a ranking claim). Deltas vs `base` per cell.

Noise bound: iron and folate are measured on the same women in the same
clusters, so the test overstates borrowable signal. Every row is tagged
`shares_population` (target and every block column in the same population)
and the summary reports cross-population-only rows separately. Median
clusters per district is 1 in Gambia, Ghana and Malawi, so a cluster split
is infeasible; CE-01 (`ceiling_cluster_split.csv`) bounds the cluster share.

Guards:

1. `base` on `cs` reproduces `civ_transport_guards.csv` ("4 countries only")
   to 3 dp.
2. No block column name is present in the covariate matrix.
3. Every block column has ≥ 12 finite values in every country of the cell.

Outputs:

- `results/tables/protocol_v2/cross_outcome_loco.csv` (cell × arm rows)
- `results/tables/protocol_v2/cross_outcome_summary.csv`
- `results/figures/protocol_v2/fig_cross_outcome_delta.png`
- XO-01 entry in `docs/findings/SANDBOX_LOG_2026-09.md`
