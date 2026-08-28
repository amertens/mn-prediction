# WS4. Distributional estimator across every outcome

## What was asked, and what the evidence says

The brief asked to promote the threshold-integrating distributional estimator to
the default, on the strength of its improving vitamin A in the one-country
prototype. The full grid does not support that promotion. The estimator helps
where theory says it should, hurts elsewhere, and on the individual-level
metrics it is worse on average. The recommendation below is to keep it as a
declared alternative for tail-prevalence outcomes rather than to make it the
default, and both arms are written to the same table under a `scheme` column so
either can be selected.

## What was built

`R/corrected/p12_distributional.R` runs two arms under leave-one-Admin-2-out
cross-validation within country:

| Scheme | Model |
|---|---|
| `binary_classifier` | ridge logistic on the individual binary deficiency outcome |
| `distributional_default` | ridge gaussian on the individual log biomarker, with prevalence recovered by integrating the empirical training residual distribution past the cut-point |

Both produce an individual-level predicted probability, which is what makes AUC
and Brier comparable between them, and both aggregate to Admin-2 after the
individual fit rather than before it.

Harmonized inputs are used where WS3 provides them, so this workstream is not
silently re-measuring adjustment heterogeneity: BRINDA CRP+AGP adjusted RBP for
vitamin A, BRINDA CRP+AGP adjusted ferritin for iron. Folate, B12 and zinc have
no inflammation adjustment defined in this pipeline and use their configured
continuous column, which the `harmonized_marker` column records per row.

The covariate set is deliberately the prototype's fixed a-priori eight
(`DIST_COVARIATE_PATTERNS`). Widening it would confound "distributional beats
binary" with "more covariates beat fewer".

Scope: all 24 configured country-outcome cells across the four countries
(source: `results/tables/corrected/distributional_comparison.csv`, one row per
cell per scheme). Every configured outcome turned out to have a usable
continuous biomarker, so no cell was skipped and no gaps file was written; the
script writes one when a cell has to be dropped, and never imputes.

## Result

Twenty-one of the 24 cells carry a noise ceiling distinguishable from
zero. The other three have no predictable between-district variation at all, so
a correlation computed on them is noise; the headline is therefore read on the
21 (source: `results/tables/corrected/distributional_paired_signal.csv`,
grouped by `metric`):

| Metric | Cells distributional better | Mean binary | Mean distributional | Mean delta | Range of delta |
|---|---|---|---|---|---|
| `pearson_r` | 12 of 21 | -0.0964 | -0.0180 | +0.0784 | -0.5714 to +0.8203 |
| `calib_slope` | 14 of 21 | -4.2642 | -2.1626 | +2.1016 | -16.5646 to +43.4202 |
| `auc` | 8 of 21 | 0.6080 | 0.5795 | -0.0285 | -0.3567 to +0.0592 |
| `brier_skill` | 13 of 21 | -0.0176 | -0.0314 | -0.0138 | -0.2568 to +0.0110 |
| `mae_pp` | 11 of 21 | 10.5170 | 10.7235 | +0.2065 | -0.9194 to +2.8259 |
| `r_share` | 4 of 13 | 0.0962 | -0.1254 | -0.2215 | -3.4200 to +0.7600 |

For `mae_pp` a lower value is better; for the rest a higher value is better,
and for `calib_slope` the target is 1 rather than a maximum, so the move from
-4.26 towards -2.16 is an improvement in the sense of being less badly inverted.

The split is consistent: the distributional arm improves the **area-level
correlation and the calibration slope**, and worsens the **individual-level
discrimination, Brier skill and area MAE**. It is not a uniform improvement on
any reading.

The unrestricted version of the same table, covering 22 to 23 cells depending on
which metric is defined (`results/tables/corrected/distributional_paired.csv`),
tells the same story:
`pearson_r` delta +0.0773, `auc` delta -0.0256.

## Where it does help, and why

The prototype predicted the gain concentrates where the cut-point sits off the
median. That holds. The largest gains are in the low-prevalence cells
(source: `results/tables/corrected/distributional_paired_signal.csv`, rows with
`metric = pearson_r`, read against national prevalence in
`distributional_comparison.csv`):

| Cell | National prevalence | Binary r | Distributional r | Delta |
|---|---|---|---|---|
| Ghana women_vitA | 1.2 percent | -0.7594 | 0.0609 | +0.8203 |
| Sierra Leone child_iron | 8.1 percent | -0.6984 | -0.2151 | +0.4833 |
| Sierra Leone women_vitA | 1.2 percent | -0.6045 | -0.2573 | +0.3472 |
| Gambia women_vitA | 1.9 percent | 0.0927 | 0.3041 | +0.2114 |
| Ghana women_b12 | 6.9 percent | -0.1148 | 0.0795 | +0.1943 |

Every one of those is under 10 percent prevalence. The largest losses are
elsewhere: Sierra Leone women_iron at 17.9 percent moves -0.2694 to -0.8408, and
Ghana child_vitA at 13.2 percent moves 0.1392 to 0.0294.

A separate observation the table forces: the mean binary `pearson_r` across
signal cells is **negative** (-0.0964). Within-country leave-one-area-out
prediction from these eight covariates is failing more often than it succeeds,
under either arm. That is consistent with the rest of this project's findings
and is not caused by the choice of arm.

## Recommendation

Do not promote `distributional_default` to the default.

Grounds, all from the tables above: it improves 2 of 6 metric families and
worsens 4; the individual-level metrics that the corrected layer reports for
discrimination and Brier skill both move the wrong way; and `r_share`, which is
the correlation read against what sampling noise permits, is worse in 9 of the
13 cells where it is defined.

Keep it as a declared alternative selected by outcome prevalence. The
prevalence-dependence is not a post-hoc rationalisation: it was predicted by the
prototype before this grid was run and it is what the efficiency argument for
continuous over binary implies. A rule of the form "use the distributional arm
when national prevalence is below roughly 10 percent" is consistent with every
cell in the table above, but four countries and five low-prevalence cells are
not enough to fit a threshold, so the rule is stated as an observation and not
as a tuned parameter.

Both arms are written to one table under a `scheme` column, so adopting either
later is a selection rather than a re-run.

## Reviewer passes

**Statistical reviewer.**

- *Leakage.* Leave-one-Admin-2-out. Every fold's imputation medians,
  standardization centre and scale, glmnet fit, and the residual distribution
  used for the integration all come from training areas only. The held-out
  area contributes to no training quantity, and `train_prev` for the Brier
  null is likewise computed per fold on training rows.
- *Identical treatment of the arms.* Both arms see the same folds, the same
  covariates, the same standardization and the same weights. Only the modelled
  outcome differs.
- *Selection on evaluation folds.* None. The covariate set is fixed a priori
  and is the prototype's, not selected here.
- *Survey weights.* Individual survey weights enter both `cv.glmnet` fits as
  observation weights and the area aggregation as weighted means. Weights that
  are missing or non-positive are set to 1, which is recorded in the code.
- *Evaluation restriction.* Areas with fewer than 8 biomarker reads are dropped
  from the EVALUATION only, never from training, matching the prototype's
  `MIN_NSVY`. The noise ceiling is computed on the same evaluation subset
  the correlation is computed on.
- *Signal restriction.* The 21-cell headline excludes 3 cells whose ceiling is
  indistinguishable from zero. The flag uses the optimistic upper bound, so a
  cell is excluded only when even that bound says there is nothing to predict.
  Both the restricted and unrestricted summaries are written out, and they agree.
- *Seeds.* `cv.glmnet` fold assignment is seeded from `params$seed` and the seed
  is recorded per row.
- *Inference.* Cell counts, means and ranges. No p-value; 21 cells across 4
  countries do not support one, and the cells are not independent within country.
- *Residual distribution.* Empirical rather than Gaussian, because
  `sensitivity/11_distributional_heterosked_transport.R` is the existing reason
  not to assume constant variance.

**Reproducibility reviewer.**

- No `targets` node added, so `tar_manifest()` and `tar_network()` are unchanged
  from WS3 at 845 targets and 2218 edges. The estimator is driven from
  `scripts/run_distributional.R`; wiring it as a target is deferred until an
  arm is actually adopted, so the pipeline does not carry a comparison that
  informs no downstream target.
- The smoke profile runs and writes to `*_SMOKE.csv` filenames.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
- No new file input, so no stamp target is required.
- All outputs are new files under `results/tables/corrected/`. The regression
  gate reports the frozen baselines unchanged.
