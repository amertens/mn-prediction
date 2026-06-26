# Simplified subset — all micronutrient outcomes, 4 countries

A small, self-contained slice of the micronutrient-prediction project so a
collaborator can run SuperLearner analyses without the full pipeline. **All
micronutrient-deficiency outcomes**, four countries, three aggregation levels,
fully documented columns, and a runnable example. Everything here is built from
the project's real pipeline outputs — no values are fabricated.

## Wide multi-outcome layout

One row per area; a block of four columns **per outcome** plus a shared set of
16 proxy predictors. The wide layout lets you run per-outcome models **and**
multi-outcome / sequential models (predict one deficiency using another as a
covariate, exploiting the covariance among deficiencies).

| File | Level | Rows | Cols |
|---|---|---|---|
| [`data/mn_cluster.csv`](data/mn_cluster.csv) | survey cluster (PSU) | 323 | 52 |
| [`data/mn_admin2.csv`](data/mn_admin2.csv) | district (Admin-2) | 206 | 52 |
| [`data/mn_admin1.csv`](data/mn_admin1.csv) | region (Admin-1) | 53 | 51 |

Rows per country (Admin-2): Gambia 30, Ghana 75, Malawi 87, Sierra Leone 14.

### Columns
- **Identifiers:** `country`, `admin1`, `admin2` (+ `cluster_id` in the cluster
  file; `n_clusters` in the Admin-2/Admin-1 files).
- **Per-outcome block** (one set per nutrient, `NA` where that outcome was not
  measured in a country):
  - `prev_<outcome>` — survey-weighted prevalence (the value to model)
  - `n_<outcome>` — sample size / denominator (good area weight)
  - `ndef_<outcome>` — number deficient / numerator (unweighted prev = ndef/n)
  - `var_<outcome>` — design-based sampling variance of the prevalence
    (`p(1-p)/effective_n`, `effective_n = n/deff`; deff = 1.5 at Admin-2/Admin-1,
    1 at cluster). Use as the known sampling-error variance in an area-level /
    Fay-Herriot model, or `1/var` as a precision weight.
- **The 16 shared proxy predictors** (same for every outcome).

### Outcomes present (and availability)
`women_iron`, `women_vitA`, `child_iron`, `child_vitA` — all four countries.
`women_b12`, `women_folate` — Ghana, Malawi, Sierra Leone (not Gambia).
`women_zinc`, `child_zinc` — Malawi only.
(Counts of Admin-2 areas with a value: iron/vitA 206, b12/folate 176, zinc 87.)

## The predictors (16) and why

Source: the project's permutation variable-importance for iron in women, then two
filters — (1) keep only predictors present in all four countries (this drops the
year-stamped DHS area indicators, leaving the geospatial proxies, which are
extracted identically everywhere); (2) drop redundancy (greedy |r| < 0.85). The
result is a compact, non-redundant, cross-country proxy set spanning rainfall
(TRMM), temperature (MODIS LST), vegetation (NDVI), land-use / human footprint,
and malaria (Malaria Atlas). They are kept as a **common predictor set for all
outcomes**. Every column (description, units, source, role) is in
[`data_dictionary.csv`](data_dictionary.csv) — no other reference needed.

## How to run the example

```r
install.packages(c("SuperLearner", "glmnet", "ranger"))   # once, if needed
```
```
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla "simplified subset/run_superlearner_example.R"
```
(Quote the path — the folder name contains a space.)

[`run_superlearner_example.R`](run_superlearner_example.R) loads one level,
fits a SuperLearner (mean / glm / elastic-net / random forest) weighted by
sample size, prints **cross-validated** RMSE/R² and ensemble weights, and then
runs a **sequential multi-outcome demo** (predict `women_vitA` from proxies,
then from proxies + `women_iron`) to illustrate borrowing strength across
deficiencies. Edit `LEVEL` / `OUTCOME` at the top to explore.

To rebuild/extend the datasets from the pipeline, run
[`build_simplified_subset.R`](build_simplified_subset.R) (needs the project's
`_targets_full` store; collaborators do not need this to use the CSVs).

## Caveats

- **Outcome definitions vary slightly by country** (e.g. iron uses
  country-specific inflammation adjustment — BRINDA/Thurnham/ferritin<15); the
  Gambia iron prevalence is notably higher, partly definitional. Treat
  cross-country *level* comparisons with care; *ranking* is more robust.
- **DHS predictors are deliberately excluded** (year-stamped, not shareable
  across countries), so a country's single strongest predictor may not appear.
- **Predictors are cluster-resolution geospatial buffers**, constant within a
  cluster and averaged up to districts/regions (mean over the area's clusters);
  not survey-derived, so no outcome leakage.
- **No population denominators** — `n_<outcome>` is the surveyed sample size, for
  weighting, not population counts.
- **Sparse areas:** Sierra Leone has only 14 districts / 4 regions; Ghana's 90
  clusters spread over 75 districts (~1 cluster/district). Area-level models will
  be noisy there.
- **Sequential transport caveat:** when conditioning one outcome on another for
  *cross-country* prediction, feed the cross-validated *prediction* of the
  upstream outcome, not its observed value (which is unavailable in a target
  country) — see [`ANDY_KIM_PROJECT_PLAN.md`](ANDY_KIM_PROJECT_PLAN.md).
