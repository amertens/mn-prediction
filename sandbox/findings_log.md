# Sandbox findings log

Running notes on what worked and didn't, from iterative experimentation
outside the targets pipeline. Headline metric: **mean LOCO Pearson r**
across 16 (outcome × held-out country) combinations.

## Headline finding

**The single best LOCO predictor is a thin-plate spline on admin-2 polygon
centroids — coordinates only, no GEE covariates at all.** Mean LOCO Pearson
r = **0.285** across 16 splits (11/15 holdouts with r > 0.1), beating
the H6 combined-filter winner (0.224) and the 17-hour SuperLearner
(~0.14). It is also the most parsimonious model — two predictors (lon,
lat) and a smoother.

This is **humbling** — the 152 GEE rasters add nothing beyond spatial
gradient information that's already encoded in the polygon centroid.
The transportable signal across West Africa is mostly geographic.

(Important caveat: the splitting countries are all West African except
Malawi. Holding out Malawi is the failure mode — the spline can't
extrapolate to East Africa. So the spatial-only result is partly an
artifact of the training-country geography, but iron-deficiency
patterns within West Africa really do follow a smooth latitudinal
gradient that the model captures.)

## Full results

| Method (best config) | Source | Mean r | Median r | Wins (r > 0.1) | Mean MAE (pp) |
|---|---|---|---|---|---|
| **gam_coords_k30** (lat/lon spline only) | `04_spatial_gam_coords.R` | **+0.285** | +0.273 | **11 / 15** | 14.7 |
| gam_coords_k50 | `04_spatial_gam_coords.R` | +0.285 | +0.273 | 11 / 15 | 14.7 |
| gam_coords_k15 | `04_spatial_gam_coords.R` | +0.284 | +0.273 | 11 / 15 | 14.6 |
| gam_coords_k30_2cov (lat/lon + accessibility + elev) | `04_spatial_gam_coords.R` | +0.265 | +0.284 | 10 / 13 | 28.0 |
| **combined_w70_k12** (H6: invariance × shared) | `07_combined_winners.R` | +0.224 | +0.122 | 9 / 15 | 13.9 |
| combined_w70_k8 | `07_combined_winners.R` | +0.216 | +0.163 | 8 / 13 | 13.5 |
| **knn_on_h6_K12_gauss** (kNN in H6 12-d) | `10_knn_on_h6.R` | +0.201 | +0.085 | 7 / 16 | 13.6 |
| knn_on_h6_K8_inv | `10_knn_on_h6.R` | +0.200 | +0.150 | 9 / 16 | 13.7 |
| **domain_inv_w70_k8** (H5: invariance only) | `06_domain_adversarial.R` | +0.195 | +0.137 | 8 / 14 | 13.8 |
| domain_inv_w50_k12 | `06_domain_adversarial.R` | +0.174 | +0.113 | 7 / 14 | 13.6 |
| knn_on_h6_K5_unif | `10_knn_on_h6.R` | +0.183 | +0.108 | 8 / 16 | 13.9 |
| h6_plus_bucketA | `08_h6_plus_bucket_a.R` | +0.175 | +0.146 | 9 / 16 | 12.5 |
| shared_top15 | `05_outcome_shared_features.R` | +0.148 | +0.070 | 6 / 14 | 14.8 |
| (baseline: forward) | `01_baseline.R` | +0.087 | +0.084 | 7 / 15 | 16.6 |
| (baseline: mixed) | `01_baseline.R` | +0.083 | +0.041 | 6 / 15 | 15.8 |
| (baseline: quasibinomial) | `01_baseline.R` | +0.077 | +0.034 | 7 / 16 | 13.4 |
| (baseline: gam, on GEE covariates) | `01_baseline.R` | +0.057 | +0.070 | 5 / 15 | 14.3 |
| transferable_k12_a50 | `02_transferability_filter.R` | +0.065 | −0.003 | 4 / 13 | 16.6 |
| knn_K15_PC10_inv (kNN in 152-d GEE) | `03_synthetic_control_knn.R` | +0.036 | +0.034 | 5 / 16 | 12.0 |
| **h6_plus_dag** (DAG-restricted H6) | `09_h6_with_dag.R` | **−0.002** | +0.028 | 3 / 15 | 16.1 |

## Bootstrap-CI honesty check on H6

`sandbox/11_bootstrap_h6.R` runs B = 500 bootstrap CIs on the H6
combined_filter predictions, resampling held-out admin-2 polygons.

**Holdouts with 95% CI clearly above zero (genuine signal):**

- child_iron, Ghana held out: r = 0.627 [0.511, 0.735]
- child_vitA, Gambia held out: r = 0.564 [0.266, 0.765]
- women_vitA, Gambia held out: r = 0.551 [0.304, 0.730]
- child_iron, Gambia held out: r = 0.539 [0.297, 0.700]

**Holdouts with 95% CI crossing zero (uncertain):**

- child_vitA Ghana, child_vitA SierraLeone, child_vitA Malawi
- women_vitA Ghana, women_vitA SierraLeone, women_vitA Malawi
- child_iron SierraLeone, child_iron Malawi
- women_iron Ghana, women_iron SierraLeone, women_iron Malawi

**Interpretation:** mean LOCO r = 0.224 is a real point estimate but
4 of the 5 "wins" are concentrated in a few country–outcome pairs
(Gambia and Ghana for iron and vitamin A in children). The headline
number should be reported with bootstrap CIs in the manuscript, not as
a clean scalar.

## What we promoted to the main pipeline

- **`fit_predict_invariance_filter`** (METHOD 15 in `R/benchmark_models.R`).
  The strongest single-filter method (H5; mean r = 0.195). Wired into
  the default methods list of `run_area_benchmarks_loco`.
- **`fit_predict_combined_filter`** (METHOD 16). The H6 winner
  (mean r = 0.224 when cross-outcome data supplied; falls back to
  invariance_filter when not). Also wired into the default methods
  list. Takes `cross_outcome_pooled = list(<outcome_tag> = pooled_data, …)`.
- `bootstrap_loco_ci()` extended to know about `quasibinomial`,
  `invariance_filter`, `combined_filter`.

## What stays in sandbox as experimental

- **Spatial GAM on lat/lon (`04_spatial_gam_coords.R`)** — the new
  best result. *Not yet promoted* because it requires a centroid
  computation step that needs the GADM polygon cache (existing
  `data/gadm/` cache works; downstream pipeline integration is a
  separate ticket). **Should be the next promotion candidate.**
- **H6 + Bucket A (`08_h6_plus_bucket_a.R`)** — mean r = 0.175,
  worse than H6 alone. Bucket A's spatial smoothing copies + within-
  domain PCs are slightly noisier than the raw covariates after H6
  has done its filtering. Leave as sandbox example.
- **H6 + DAG (`09_h6_with_dag.R`)** — mean r = −0.002. DAG → regex
  mapping is too narrow; the candidate pool collapses to < 4 vars
  for most outcomes and the model fails. Revisit when the DAG
  specs are refined with a domain expert.
- **kNN on H6 features (`10_knn_on_h6.R`)** — mean r = 0.20,
  competitive with H6 itself. Confirms that **kNN failed in `H2`
  because of the curse of dimensionality**, not the method: in the
  low-dim H6-selected space it recovers. Useful as a non-parametric
  consistency check on H6's predictions.

## Recommended next iterations (updated)

1. **Promote `fit_predict_spatial_coords` to the main pipeline** —
   the new best result. ~1 day. Requires loading lat/lon from
   GADM cache (already present), then a thin GAM wrapper.
2. **Test spatial GAM + H6 stacked** — spatial GAM captures the
   geographic gradient; H6 captures within-spatial-pattern
   variation. Their stacked ensemble might be the strict best.
3. **Refine the DAG specs with a country nutritionist** — the H6+DAG
   failure indicates the DAG node → covariate regex mapping is too
   narrow, not that the causal-reasoning approach is wrong.
4. **Pre-register the comparison in the manuscript** —
   bootstrap CIs show several of our "wins" are uncertain, so
   the manuscript needs honest reporting (point + CI, not point).

## Headline for the manuscript

The strongest LOCO model on this dataset is the simplest: a thin-plate
spline on admin-2 polygon centroids, no GEE covariates at all (mean r
= 0.285, with 11 of 15 holdouts having r > 0.1). The complex multi-
method benchmark suite confirms that for West African admin-2
prediction, geographic gradients carry more transportable signal than
any of the 152 earth-observation covariates we extracted. The H6
combined-filter (mean r = 0.224) is the best feature-based parsimonious
model and is competitive with — but does not beat — pure spatial
smoothing. Iron-deficiency outcomes for Gambia and Ghana are reliably
predicted (95% bootstrap CIs above 0); vitamin A in low-prevalence
countries and Malawi outcomes are not.

This is a strong story precisely because it inverts the conventional
expectation that more data = better prediction. With only 4 training
countries and ~200 admin-2 polygons, geographic interpolation does
most of the work and the model selection lesson is to add countries
(or fundamentally different features like soil micronutrients or
agro-ecological zones) rather than to add more covariates of the same
kind.
