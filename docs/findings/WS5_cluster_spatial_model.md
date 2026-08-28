# WS5. Cluster-level spatial model

## Scope change, and why

The plan gated this workstream on a displacement-integrated Earth Engine
extraction of covariates at cluster points, estimated at 8 to 20 hours. That was
dropped before any Earth Engine time was spent, on three pieces of evidence
already in the repository:

1. `sandbox_parsimony/FINDINGS.md` section 12 measured the reliability ceiling by
   spatial level: Admin-1 `r_max` 0.59, Admin-2 0.31, cluster 0.22. Going finer
   lowers the achievable correlation, because units get noisier faster than they
   multiply. That is a sampling-noise fact no covariate improvement touches.
2. `results/cluster_vs_admin2_comparison.csv` has cluster-level correlation
   worse than Admin-2 in 6 of 10 comparable cells and MAE worse in 8 of 10.
3. The spatial smoother needs coordinates, not cluster covariates, and the
   coordinates already exist: 70 distinct cluster locations for Gambia, 90 for
   Ghana, 60 for Sierra Leone, 105 for Malawi.

The covariate arms therefore use Admin-2 covariates joined onto clusters. Those
are constant within district (FINDINGS.md section 12: Ghana has 75 distinct
covariate values across 2,450 rows, a between-district variance share of 1.00).
That limitation is recorded per row in `covariates_district_constant` rather
than left implicit, and it is precisely the limitation the dropped extraction
would have fixed. As the results below show, the covariate arms lose to the
covariate-free ones by a wide margin, so fixing them is not where the value is.

## Design

`R/cluster_mbg.R`. Leave-one-Admin-2-out: every cluster in a district is held
out together, the model is fitted on the rest, and the district is predicted.
That is the unsurveyed-district scenario, and the survey design makes it a
demanding one. The median district has **one** cluster in Gambia, Ghana and
Malawi; 62 of 75 Ghana districts and 76 of 89 Malawi districts have exactly one.
Only Sierra Leone has 2 to 8 clusters per district. Holding out a district
usually means holding out a single point and asking the smoother to extrapolate.

Five arms: `national_mean` (the null a district with no clusters would
otherwise get), `covariates_only`, `spatial_only` (`s(lon, lat)` alone),
`spatial_plus_covariates`, and `matern_spamm` (a spaMM Matérn field, wired as an
optional engine). The thin-plate basis is capped at one third of the training
clusters so the smoother cannot interpolate its own training data.

District predictions aggregate cluster predictions weighted by cluster biomarker
count. Population weighting needs WorldPop, which was not run; every row records
`aggregation = "cluster-count weighted (not population weighted)"`.

Scope: all 24 country-outcome cells, no skips
(source: `results/tables/corrected/cluster_mbg_within_country.csv`).

## Result

Restricted to the 21 cells whose reliability ceiling is distinguishable from
zero (source: same file, grouped by `arm`):

| Arm | Mean Pearson | Mean MAE (pp) | Mean RMSE (pp) | Mean `r_share` |
|---|---|---|---|---|
| `spatial_only` | 0.1357 | 10.087 | 12.663 | 0.319 |
| `matern_spamm` | 0.0062 | 9.947 | 12.591 | -0.099 |
| `covariates_only` | 0.1374 | 13.873 | 18.773 | 0.380 |
| `spatial_plus_covariates` | 0.1371 | 15.106 | 20.701 | 0.288 |
| `national_mean` | not interpretable | 10.599 | 13.019 | not interpretable |

The `national_mean` correlation is blanked because it is mechanical rather than
noisy: its prediction is the leave-one-out training mean, which is
anti-correlated with the held-out value by construction, since removing a
high-prevalence district drags the mean down. Its MAE and RMSE are meaningful
and are what the null is read on.

**Two things stand out.**

**Covariates hurt.** `covariates_only` has MAE 13.873 against the null's 10.599,
and adding covariates to the spatial model moves it from 10.087 to 15.106. The
covariate arms are worse than doing nothing. Their slightly higher correlation is
not a redeeming feature: a model can rank districts marginally better while
placing all of them at the wrong level, which is the same pattern WS3 and WS7
document.

**Spatial borrowing helps, slightly.** Against the null, per cell
(source: `results/tables/corrected/cluster_mbg_vs_null.csv`, and the per-cell
comparison the script prints):

- `spatial_only` is better in **13 of 24** cells, mean delta **-0.48 pp** MAE,
  median -0.108 pp.
- `matern_spamm` is better in 14 of 24, mean delta -0.60 pp.

The largest gains are Gambia women's iron (-4.217 pp), Malawi women's folate
(-3.746) and Ghana child iron (-3.663). The losses concentrate in Sierra Leone,
which has only 14 districts, so holding one out removes a large share of the
data and forces a long extrapolation, and in Malawi's zinc cells.

## The convergence worth noting

A mean gain of roughly half a percentage point of MAE over "give the district
the national average", in about half of cells, is close to the ~0.7 pp
zero-coverage gain that `docs/PROJECT_STATUS_2026-08_UPDATE.md` section 5
reports for the covariate-based area model in the archived subsample study.

Two different models, two different evaluation designs, and the same order of
magnitude. That figure is not yet reproducible from code in this repository (see
WS0), so this is a convergence of magnitudes rather than a confirmation. But it
does mean the narrow framing in that status note, a smarter guess for districts
you deliberately do not sample, survives being approached from a second
direction, and that a point-referenced spatial model does not change it.

## The Matérn arm is not recommendable

`matern_spamm` has the best mean MAE by a hair and essentially zero mean
correlation (0.0062), with 4 of 24 cells below -0.5 and a range of -0.747 to
0.579. The run produced **391** "nearly-singular correlation matrix" warnings
from spaMM. With 60 to 105 clusters and a leave-one-district-out loop, the
Matérn correlation matrix is frequently near-singular and the resulting
predictions are unstable. It is reported for completeness and should not be
used. `spatial_only` achieves the same error with a stable fit.

## Answer to the question the plan posed

Does spatial borrowing change the picture for unsurveyed districts? Marginally.
It beats the national mean by about half a percentage point of MAE in roughly
half of cells, it never beats it by much, and it loses in the country with the
fewest districts. It does not rescue the covariate story, and combining the two
is worse than either alone.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* Leave-one-Admin-2-out with every cluster in the district held out
  together. Covariate standardization uses training clusters only. The null arm's
  training mean is recomputed per fold.
- *Degenerate and mechanical correlations.* Two separate guards. Any arm whose
  predictions are near-constant relative to the truth has its correlation set to
  NA with `correlation_degenerate` recorded. The null arm is additionally
  flagged `correlation_is_mechanical` and blanked from correlation summaries,
  because its anti-correlation is a property of the CV scheme rather than of the
  model. Without this the null would appear to have a large negative correlation
  and would flatter every other arm.
- *Survey weights.* The cluster response is the survey-weighted cluster
  prevalence, and `mgcv` receives it with `weights = n`. This is stated in the
  code as a deliberate deviation from `cbind(k, n - k)`: count form cannot carry
  a non-integer survey weight, and with equal weights the two are identical.
- *Basis dimension.* Capped at one third of the training clusters, so the
  smoother cannot interpolate its own training data. Without the cap, mgcv's
  default of 30 basis functions on 60 clusters would.
- *Aggregation.* Cluster-count weighted, not population weighted, recorded on
  every row. This is a documented approximation, not a population-weighted
  national quantity.
- *Failure handling.* An arm that cannot be fitted returns NULL and is absent
  from the output; it never falls back to another arm's prediction. This matters
  after the WS3 experience with a silent fallback masking a failure.
- *Seeds.* Recorded per row.
- *Inference.* Twenty-four cells, counts and means only. No p-value; cells are
  not independent within country.

**Reproducibility reviewer.**

- No `targets` node added. The model is script-driven pending a decision on
  whether to adopt it; wiring it in would add a target that informs nothing
  downstream.
- `cluster_mbg_covariates()` is a function rather than a constant because
  `tar_source()` loads `R/` alphabetically and `R/cluster_mbg.R` precedes
  `R/corrected/p12_distributional.R`. A top-level reference would not resolve.
  Deferring to call time keeps one definition instead of two that can drift.
- The spaMM formula sets its environment to the spaMM namespace rather than
  using `spaMM::Matern(...)` inside the formula, which fails with "the condition
  has length > 1". This was found by direct test, not by inference.
- Smoke profile runs and writes to `*_SMOKE.csv`.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
- All outputs are new files. The regression gate reports the frozen baselines
  unchanged.
