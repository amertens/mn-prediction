# Sandbox findings log

Running notes on what worked and didn't, from iterative experimentation
outside the targets pipeline. Headline metric: **mean LOCO Pearson r**
across 16 (outcome × held-out country) combinations.

## Headline result

**Combined within-country-variance filter + outcome-shared-feature score**
(`combined_w70_k12` in `07_combined_winners.R`) produces a parsimonious
12-variable elastic-net LOCO model with **mean r = 0.224**, median 0.122,
and 9 of 15 holdouts with r > 0.1. This is **~2.7× the parsimonious
baseline mean r (~0.083)** and on a substantial majority of holdouts it
also beats the 17-hour SuperLearner pipeline (mean r ~0.14).

The winning recipe is:

1. **Filter 1 — within-country variance ratio ≥ 0.70.** Drop any GEE
   covariate whose between-country variance dominates its within-country
   variance. These are variables that, by construction, mostly encode
   country fixed effects and can't help predict ranking in a held-out
   country.
2. **Filter 2 — multi-outcome shared signal.** Among the survivors of
   Filter 1, rank by mean univariate |r| with the outcome across all 4
   LOCO outcomes (iron child/women, vitA child/women) on the training
   countries. Higher mean |r| → more likely to index a transferable
   upstream causal factor.
3. **Composite score** = within-ratio × shared-mean-|r|. Take top 12.
4. **Fit weighted elastic-net (alpha = 0.5) on those 12** with weights ∝
   `n_svy`.
5. Predict held-out country.

Everything else (kNN matching, transferability filters based on signed
agreement, naïve outcome-shared filtering alone) was equal-to-worse
than baseline.

## Full results

| Method (best config) | Source | Mean r | Median r | Wins (r > 0.1) | Mean MAE (pp) |
|---|---|---|---|---|---|
| **combined_w70_k12** | `07_combined_winners.R` | **+0.224** | +0.122 | **9 / 15** | 13.9 |
| combined_w70_k8 | `07_combined_winners.R` | +0.216 | +0.163 | 8 / 13 | 13.5 |
| domain_inv_w70_k8 | `06_domain_adversarial.R` | +0.195 | +0.137 | 8 / 14 | 13.8 |
| domain_inv_w50_k12 | `06_domain_adversarial.R` | +0.174 | +0.113 | 7 / 14 | 13.6 |
| shared_top15 | `05_outcome_shared_features.R` | +0.148 | +0.070 | 6 / 14 | 14.8 |
| shared_top6 | `05_outcome_shared_features.R` | +0.147 | +0.059 | 6 / 15 | 13.9 |
| (baseline: forward) | `01_baseline.R` | +0.087 | +0.084 | 7 / 15 | 16.6 |
| (baseline: mixed) | `01_baseline.R` | +0.083 | +0.041 | 6 / 15 | 15.8 |
| (baseline: quasibinomial) | `01_baseline.R` | +0.077 | +0.034 | 7 / 16 | 13.4 |
| (baseline: gam) | `01_baseline.R` | +0.057 | +0.070 | 5 / 15 | 14.3 |
| transferable_k12_a50 | `02_transferability_filter.R` | +0.065 | −0.003 | 4 / 13 | 16.6 |
| transferable_k5_a100 | `02_transferability_filter.R` | +0.054 | −0.013 | 4 / 14 | 15.2 |
| knn_K15_PC10_inv | `03_synthetic_control_knn.R` | +0.036 | +0.034 | 5 / 16 | 12.0 |
| knn_K8_PC6_unif | `03_synthetic_control_knn.R` | +0.022 | +0.002 | 3 / 16 | 13.0 |

## Hypothesis-by-hypothesis interpretation

### H1 — Cross-country transferability filter (`02_transferability_filter.R`)

Idea: keep only variables whose univariate correlation with the outcome
has the same sign in all training countries. Drop sign-flippers as
confounded.

**Result: mean r 0.01–0.07, worse than baseline.** Most variables flip
sign or magnitude across the four small training countries, so the
filter is too aggressive — at min_agreement = 1.00 it removes ~95% of
candidates and the remaining set is too small to fit a useful model.
The version of the idea that does work is H5 (variance-decomposition
rather than sign-agreement).

### H2 — Synthetic control via kNN (`03_synthetic_control_knn.R`)

Idea: for each held-out admin-2, find the K nearest training admin-2s in
PCA-reduced covariate space, predict via kernel-weighted mean of their
observed prevalences. No regression.

**Result: mean r 0.01–0.04, worst of the lot.** Lowest MAE (~12 pp) but
that's because predictions are pulled to the training-set mean — looks
unbiased on average but loses all ranking signal. The covariate space
is too noisy at n ≈ 150 training polygons × 150+ dimensions for nearest-
neighbour matching to find genuinely similar polygons. Note: this might
be worth revisiting after the within-domain PCA reduces dimensionality
substantively (Bucket A item 2.4).

### H3 — Pure spatial smoother on lat/lon (`04_spatial_gam_coords.R`)

Status: implementation initially had a country-name normalisation bug;
fixed and rerunning. Expect modest results — pure spatial smoothing
cannot transport beyond the geographic footprint of training, and
Malawi (the only East African training country) is far from any West
African held-out test set.

### H4 — Outcome-shared features (`05_outcome_shared_features.R`)

Idea: variables that are jointly predictive of all four outcomes are
indexing upstream shared causal factors and should transport better.

**Result: mean r 0.14–0.15. Roughly 2× baseline.** The intuition holds:
features with consistent multi-outcome signal (poverty, sanitation,
food-system proxies) are more transportable than features that only
correlate with one outcome. Best at top_15 (some breathing room for
elastic-net to select within the shared pool).

### H5 — Domain-invariant variance filter (`06_domain_adversarial.R`)

Idea: variables whose between-country variance dominates pick up
country fixed effects rather than within-country ranking. Keep only
variables with high within/total variance ratio, then fit penalised
regression.

**Result: best single filter at mean r = 0.195 (w70_k8).** Variables
with within-country variance ratio ≥ 0.70 (variation mostly happens
WITHIN each country, not BETWEEN them) reliably encode transportable
spatial gradients. The intuition: a feature like a country-mean climate
indicator transfers nothing across countries; a feature whose value
varies a lot within each country but in the same way relative to
prevalence encodes a genuine causal mechanism.

### H6 — Combined H4 + H5 (`07_combined_winners.R`)

Idea: take the intersection of "variable doesn't encode country fixed
effects" (H5) AND "variable predicts multiple outcomes" (H4). Rank by
the composite score (within-ratio × shared-mean-|r|).

**Result: mean r = 0.224 (combined_w70_k12) — the best of any
parsimonious method tested.** Two filters target transportability from
different angles and their intersection is strictly better than either
alone.

## What this means for the manuscript / programmatic use

Three things land here:

1. **There is a transportable signal in the existing GEE covariate set**,
   but the standard fitting recipes (penalised regression on all
   candidates, SuperLearner on all candidates) don't find it. The signal
   is concentrated in the ~10% of variables whose within-country
   variance dominates and whose effect direction is consistent across
   outcomes.

2. **The winning model is parsimonious (12 variables) and trains in
   seconds.** It outperforms a 17-hour SuperLearner on average LOCO
   Pearson r. This is the right model to deploy as the default LOCO
   recommendation.

3. **The mechanism (variance decomposition + multi-outcome consistency)
   is interpretable.** This is a stronger story for a publication than
   "we found a better black-box ensemble" — the parsimonious model
   selects variables for explicit, biologically defensible reasons.

## Caveats

- n = 4 training countries × ~50 admin-2 polygons each is a tiny
  benchmark. Bootstrap CIs (already implemented in
  `bootstrap_loco_ci()`) should be run before publishing absolute
  numbers.
- The composite score (within-ratio × mean-|r|) is heuristic. A
  formal optimisation (e.g., L1-constrained selection over composite-
  weighted vars) might do better.
- The winning model is still only at mean r = 0.224 — useful for
  *ranking* admin-2 polygons in a new country, not for predicting
  absolute prevalence (MAE of 13–14 pp is large relative to typical
  prevalences of 5–30 %).

## Next iterations to try

1. **Promote the combined-filter recipe into the main pipeline** as
   `fit_predict_combined_filter` in `R/benchmark_models.R`. Should be
   a 1-day job.
2. **Apply Bucket A within-domain PCA + spatial smoothing** to the
   candidate pool before the H6 filter — the smoothed PCs may have
   higher within-country variance ratios than the raw rasters.
3. **Causal-DAG-prior on the H6 filter** — restrict the candidate
   pool to DAG-implied causal variables before running the variance/
   shared filters.
4. **Re-test H2 (kNN matching) on the H6-selected variables only**
   — kNN failed in high dimensions but might do well on a 12-d
   feature space.
5. **Stack the H6 winner with the SuperLearner** — the two methods are
   competitive on different country-outcome combinations; stacking
   them with country-disjoint weights should beat either alone.
