# WS2b. Covariate-hygiene default flip

## Why this workstream exists

The brief made the `GEE_COVARIATE_HYGIENE` default flip the last step of a
Tanzania rebuild, so that the two changes would not be confounded. Stage 0
established that Tanzania was dropped from the pipeline on 2026-08-26
(`R/config.R:507`, `get_country_config_tanzania_archived_2010()`), so that gate
no longer exists. The flip is therefore evaluated on its own, on the four active
countries.

The flag controls `prune_gee_covariates()` (`R/gee_band_semantics.R:173`), which
removes two classes of column: cross-band `_annual_*` summaries taken over
non-commensurable bands, and static layers with identical per-year copies. The
first class holds quantities such as a soil depth mean averaged with its own
standard deviation, or an FLDAS summary averaging surface pressure with wind
speed. Those numbers do not correspond to a physical measurement whatever their
predictive behaviour, which is the semantic argument for pruning them. This
workstream supplies the empirical half.

## What was run

`scripts/compare_covariate_hygiene.R` refits the area-level LOCO transport
recipe (`AREA_TRANSPORT_RECIPE`, elastic net on the top-30 correlation-screened
predictors, `R/transportability_area.R:34`) twice on the same cached Admin-2
targets, once per variant, and does not rebuild the pipeline.

The script previously wrote in place to
`results/sensitivity/covariate_hygiene_comparison.csv`. That file is now a
frozen baseline, so the script was changed to tag every run with a `scheme`
column and write to a scheme-specific path, and to refuse the old filename
outright. The 2026-08-18 comparison is untouched.

A second script, `scripts/summarize_hygiene_flip.R`, pairs the two variants
within each combination of outcome and held-out country. The original script
reported only means across held-out countries, and a mean over four folds does
not show whether the variants move together.

Outputs:

- `results/sensitivity/covariate_hygiene_comparison_ws2b_2026-08.csv`
- `results/sensitivity/covariate_hygiene_paired_ws2b_2026-08.csv`
- `results/sensitivity/covariate_hygiene_paired_summary_ws2b_2026-08.csv`

Scope: 4 outcomes (child_vitA, women_vitA, child_iron, women_iron) by 4 held-out
countries (Gambia, Ghana, Sierra Leone, Malawi), giving 16 paired folds per
metric.

## Result

The two variants differ in covariate count: 149 predictors under `v1_current`
against 87 under `v2_hygiene`, so pruning removes 62 columns.
(source: `results/sensitivity/covariate_hygiene_comparison_ws2b_2026-08.csv`,
column `n_predictors`, distinct values by `variant`.)

Paired within-fold deltas, computed as `v2_hygiene` minus `v1_current`
(source: `results/sensitivity/covariate_hygiene_paired_summary_ws2b_2026-08.csv`,
one row per `metric`):

| Metric | Folds v2 better | Mean delta | Median delta | Min delta | Max delta |
|---|---|---|---|---|---|
| `rmse_pp` | 12 of 16 | -1.5562 | -0.405 | -13.390 | +0.380 |
| `abs_nat_bias_pp` | 11 of 16 | -1.8219 | -0.270 | -18.910 | +4.320 |
| `mae_pp` | 9 of 16 | -1.3125 | -0.220 | -12.910 | +0.810 |
| `pearson_r` | 9 of 15 | +0.0159 | +0.004 | -0.143 | +0.188 |
| `spearman_r` | 6 of 15 | -0.0125 | -0.011 | -0.153 | +0.110 |
| `n_selected` | 6 of 16 | -1.0000 | 0.000 | -6.000 | +7.000 |

For `rmse_pp`, `abs_nat_bias_pp` and `mae_pp` a lower value is better, so a
negative delta favours the pruned set. For `pearson_r` and `spearman_r` a higher
value is better.

The correlation metrics report 15 folds rather than 16 because `pearson_r` is
absent for women_vitA held out on Gambia
(source: `results/sensitivity/covariate_hygiene_comparison_ws2b_2026-08.csv`,
row `outcome = women_vitA`, `held_out = Gambia`, column `pearson_r`).

The shape of the `rmse_pp` distribution matters more than its mean. The four
folds where pruning is worse are worse by at most 0.380 pp, while the folds
where it helps include child_iron held out on Malawi, where RMSE falls from
30.61 to 17.22 (source:
`results/sensitivity/covariate_hygiene_paired_ws2b_2026-08.csv`, rows with
`metric = rmse_pp`, columns `v1_current`, `v2_hygiene`, `delta`). The downside is
bounded near zero and the upside is large.

## Decision

Flip the default to on.

The grounds are that error metrics improve or are unchanged, ranking metrics are
not degraded in any way the fold count can resolve, and the covariate set is
smaller by 62 columns. The correlation evidence on its own is ambiguous: mean
`pearson_r` moves by +0.0159 and mean `spearman_r` by -0.0125, both far inside
the fold-to-fold spread, and neither would justify a change by itself. The
decision rests on the error metrics, on the bounded downside described above,
and on the semantic argument that the pruned columns are not physical
quantities.

Four held-out countries do not support a p-value, and none is reported.

The flip is committed separately from the evidence above so that the two can be
reverted independently. Because `gee_hygiene_enabled()` is read at call time
rather than captured in a target, area-level and LOCO targets must be
invalidated after the change:

```
targets::tar_invalidate(matches("area_|loco"))
```

The setting is recorded in `pipeline_params` (`R/config.R:674`), so it enters the
target hash for everything downstream and is not invisible to `targets`.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* Both variants run through the same `run_area_transport_loco()`
  (`R/transportability_area.R:203`), which trains on the remaining countries and
  predicts the held-out one. Pruning is applied to the common covariate
  intersection before any fold is formed, identically for every country, so it
  cannot use held-out information.
- *Selection on evaluation folds.* The pruning rule is declared in
  `GEE_BAND_SEMANTICS` (`R/gee_band_semantics.R:74`) from band semantics, not
  fitted from any score, so it does not select on the metric it is judged by. The
  comparison itself is a decision about a default and is reported as such, not as
  an out-of-sample performance claim.
- *Seed coverage.* The only stochastic step is `cv.glmnet` fold assignment,
  seeded from `AREA_TRANSPORT_RECIPE$seed` inside `.tr_fit_predict()`
  (`R/transportability_area.R:176`). Both variants use the same seed, so the
  paired deltas are not contaminated by fold-draw noise.
- *Fold construction.* Leave-one-country-out over the four active countries.
  Countries are the unit of both training and evaluation.
- *Survey weights.* Admin-2 units are weighted by survey n through
  `recipe$weight = TRUE`. This is the area-level design weight, and it is applied
  identically in both variants. The individual-level SuperLearner weighting
  question is not touched by this workstream.
- *Inference.* Sixteen paired folds are reported as counts and ranges. No
  p-value is computed.

**Reproducibility reviewer.**

- No `targets` node was added or changed, so `tar_manifest()` and
  `tar_network()` are unaffected and were not re-run. The one pipeline-visible
  change is the flag default, which flows through the existing
  `pipeline_params` entry.
- Both scripts resolve paths through `here::here()`. There is no absolute path
  and no `setwd()`.
- Both scripts read only from the `targets` store and from files they were given,
  and neither writes to a frozen baseline. `scripts/compare_covariate_hygiene.R`
  stops rather than write the old filename.
- Re-running `scripts/summarize_hygiene_flip.R` on the same input reproduced the
  same summary table, confirmed by running it twice.
- No new file input was introduced, so no stamp target is required.
