# Claims register

Every claim in `docs/SESSION_FINDINGS_FOR_REVIEW.md` Sections 3, 4, 5, 6, 9 and 11,
with the status it holds after the accuracy-and-impact work of September 2026.

**Purpose.** The manuscript is not edited by this branch. This register is the
input to that later edit. A manuscript sentence is revised only where a row
below says `revised` or `withdrawn`.

**Status vocabulary.**

| Status | Meaning |
|---|---|
| `not yet tested` | This branch has not re-measured the claim. |
| `confirmed` | Re-measured and reproduced within rounding. |
| `revised` | Re-measured and different. The new value and its source map are given. |
| `withdrawn` | The claim does not survive. The reason is given. |

**Source map convention** (guardrail 3). Every value carries
`(source: <file>, column <name>, filter <...>)`. A value that cannot be produced
by running code is written `not yet computed`. Values in the "As stated" column
are transcriptions of the session document and carry no independent authority.

**Frozen baseline.** The pre-work state of every table named below is recorded in
`results/tables/frozen_2026-09/MANIFEST.csv` with an md5 per file.

---

## Section 3. Established results carried into this session

The session document labels all of Section 3 `[carried forward]` and states none
were re-verified. Statuses below therefore start from a lower base of evidence
than Sections 4 to 6.

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 3.1 | Individual-level prediction aggregated to Admin-2 | median r 0.516, r_share 0.92 | `frozen_2026-09/area_comparison_all.csv`, columns `pearson_r` and `r_share`, filter `approach == "Individual SL"` | not yet tested | WS1, WS3a |
| 3.2 | Area-level SuperLearner does not beat the national-mean null | MAE 9.31 to 9.84 vs 9.55 pp | `frozen_2026-09/area_comparison_all.csv`, column `mae_pp`, filter on the two Area SL arms and the null | not yet tested | WS3a |
| 3.3 | Covariate-free spatial smoother beats every covariate arm | r 0.304 to 0.393 | `frozen_2026-09/area_comparison_all.csv`, column `pearson_r`, filter `approach == "Spatial only (no covariates)"` | not yet tested | WS4a |
| 3.4 | Covariates beat geography at Admin-2 | 6 of 24 cells | `frozen_2026-09/resolution_comparison.csv`, column `r`, filter `level == "admin-2"`, arms `covariates` against `spatial only` | not yet tested | WS4a |
| 3.5 | Admin-1 covariates against Admin-1 spatial | 0.437 vs 0.209 | `frozen_2026-09/resolution_comparison.csv`, column `r`, filter `level == "admin-1"` | not yet tested | WS4a |
| 3.6 | WHO category accuracy, Admin-2 | 33 percent vs 32 percent null | derivation not located in a committed table | not yet tested | WS8a |
| 3.7 | WHO category accuracy, Admin-1 | 56 percent vs 35 percent null | derivation not located in a committed table | not yet tested | WS8a |
| 3.8 | National prevalence inside the survey 95 percent CI | 24 of 24, mean absolute error 0.96 pp | `frozen_2026-09/national_estimates_all.csv`, columns `obs_lo`, `obs_hi`, `pred_prev`, `abs_diff` | not yet tested | WS2c |
| 3.9 | Predictors surviving FDR control | 0 of 294, all 24 cells | `frozen_2026-09/bivariate_fdr.csv` | not yet tested | WS7a |
| 3.10 | Penalised regression retained predictors | median 0; 13 of 18 cells retain none | `frozen_2026-09/penalized_retained.csv` | not yet tested | not scheduled |
| 3.11 | SuperLearner beats the best of 21 comparators | 0 of 16 LOCO holdouts | `frozen_2026-09/benchmarks_all.csv`, filter `eval_type == "loco"` | not yet tested | not scheduled |
| 3.12 | Positive control, earth observation predicts district education | r 0.48 to 0.71 | no committed table located; Stage 0 searched `results/`, `sensitivity/`, `sandbox_parsimony/` and `docs/` | not yet tested | WS4b |

## Section 4. Where a district estimate's level comes from

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 4.1 | Admin-1 anchor, hard | mean r 0.413, median 0.405, MAE 8.85, mean absolute bias 1.59 | `frozen_2026-09/admin1_arms.csv`, filter `arm == "admin-2 fit + ADMIN-1 benchmark (hard)"` | not yet tested | WS2 |
| 4.2 | Admin-1 anchor, shrunk | mean r 0.318, MAE 9.36, bias 2.76 | same file, filter `arm == "admin-2 fit + ADMIN-1 benchmark (shrunk)"` | not yet tested | WS2 |
| 4.3 | Fit at Admin-1, extrapolate to Admin-2 | mean r 0.230, 22 cells, MAE 8.78, bias 3.18 | same file, filter `arm == "ADMIN-1 fit -> admin-2 (pooled)"` | not yet tested | WS2 |
| 4.4 | National anchor | mean r 0.170, MAE 11.98, bias 5.85 | same file, filter `arm == "admin-2 fit + national benchmark"` | not yet tested | WS2a, WS2c |
| 4.5 | No anchor | mean r 0.164, MAE 10.71, bias 3.24 | same file, filter `arm == "admin-2 fit (LORO), unbenchmarked"` | not yet tested | WS2 |
| 4.6 | The hard Admin-1 anchor beats no anchor | 20 of 24 cells | same file, paired on `country` and `outcome` | not yet tested | WS2b, WS2d |
| 4.7 | The anchor is a one-parameter logit shift, so district ranking is identical in every arm | qualitative | `R/benchmark_area.R`, `benchmark_admin2_to_admin1()` | not yet tested | WS2c |
| 4.8 | National anchoring makes mean absolute bias worse | 5.85 vs 3.24 pp | `frozen_2026-09/admin1_arms.csv`, column `bias_pp` | not yet tested | WS2c, withheld pending the implied-shift audit |
| 4.9 | Hard beats shrunk despite the theoretical risk | 0.413 vs 0.318 | `frozen_2026-09/admin1_arms.csv` | not yet tested | WS2d |
| 4.10 | The anchoring gain is not circular | qualitative argument in Section 4.4 | none | not yet tested | WS2a, WS2b |
| 4.11 | The survey supplies the level; the model supplies the within-region pattern | qualitative | Sections 4 and 6 | not yet tested | WS2d |

## Section 5. The individual-level anchor

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 5.1 | District, proxy arm | 16 cells, mean r 0.154, median 0.126, MAE 8.33 | `frozen_2026-09/individual_anchor.csv`, filter `unit == "district"` and `arm == "proxy"` | not yet tested | WS3a |
| 5.2 | District, questionnaire arm | mean r 0.228, median 0.264, MAE 7.92 | same file, filter `unit == "district"` and `arm == "quest"` | not yet tested | WS3a, WS3b |
| 5.3 | Cluster, proxy and questionnaire arms | 0.146 and 0.229 | same file, filter `unit == "cluster"` | not yet tested | WS3f |
| 5.4 | Questionnaire better in 10 of 16 cells, mean gain +0.075 | 10 of 16, +0.075 | same file, paired on `country`, `outcome` and `unit` | not yet tested | WS3b, WS3c |
| 5.5 | Clears r = 0.4 in 3 of 16 cells | 3 of 16 | same file, column `r` | not yet tested | WS3a |
| 5.6 | In Malawi the questionnaire adds nothing | gains 0.000, 0.000, 0.002, 0.004 | same file, filter `country == "Malawi"` | not yet tested | WS3b |
| 5.7 | Maximum r anywhere across all 64 rows | 0.544 | same file, maximum of column `r` | not yet tested | WS7a |
| 5.8 | The null is not explained by bad linkage, the inflammation adjustment, or overfitting | qualitative | Section 5.3 | not yet tested | WS3a, WS4b |
| 5.9 | Cluster linkage does not help, and helps least in Sierra Leone | Gambia +0.017, Ghana -0.003, Malawi -0.002, Sierra Leone -0.025 | same file, `r` at `unit == "cluster"` minus `r` at `unit == "district"` | not yet tested | WS3f |

## Section 6. A predicted national level cannot substitute for a survey

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 6.1 | The VMNIS LOCO model is competent | Vitamin A preschool random forest, MAE 11.75, Pearson 0.655 | `frozen_2026-09/national_vmnis_loco.csv` | not yet tested | WS6c |
| 6.2 | The null arm's negative Pearson is an artefact, not a finding | qualitative | `frozen_2026-09/national_vmnis_loco.csv`, filter `model == "null"` | not yet tested | WS6b |
| 6.3 | Vitamin A preschool ceiling components | sd country 1.411, method 0.564, residual 0.000, sampling 0.816 | `frozen_2026-09/national_vmnis_ceiling.csv`, row 1 | not yet tested | WS6a |
| 6.4 | Vitamin A preschool ceiling and saturation | r_max 0.818, r_share 0.80 | same file, columns `r_max_report` and `r_share` | not yet tested | WS6a |
| 6.5 | Correction to the record: the model is not 98 percent saturated at a ceiling of 0.66 | 80 percent of 0.818 | same file | not yet tested | WS6a |
| 6.6 | sd_resid hits the boundary at 0.000 in two panels, so r_max is untrustworthy there | 2 of 4 panels | same file, column `sd_resid` | not yet tested | WS6b |
| 6.7 | For Vitamin A NPW, method variance exceeds country variance | 1.996 vs 1.232 | same file, row 3 | not yet tested | WS6b |
| 6.8 | Composition arm errors | MAE 5.81, 8.22, 12.70, 5.58 pp | `frozen_2026-09/national_composition.csv`, column `mae_pp` grouped by `arm` | not yet tested | WS6d |
| 6.9 | Composition arm absolute bias | 3.35, 4.69, 10.03, 0.73 pp | same file, column `bias_pp` | not yet tested | WS6d |
| 6.10 | The predicted national level is on a different scale; the null beats the covariate model | 6 of 8 cells | `frozen_2026-09/national_composition_levels.csv` | not yet tested | WS6c |
| 6.11 | Even a perfect national level buys 0.23 pp, better in only 4 of 8 cells | 5.81 to 5.58 pp | `frozen_2026-09/national_composition.csv` | not yet tested | WS6d |
| 6.12 | Scope limitation: the VMNIS and transport outcomes intersect on vitamin A only | 2 outcomes times 4 countries, 8 cells | Section 6.5 | not yet tested | WS6d |

## Section 9. r_share above one

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 9.1 | Mean r_share by anchoring arm | hard 2.06, shrunk 1.64, Admin-1 fit 1.77, national 1.34, none 1.35 | `frozen_2026-09/admin1_arms.csv`, column `r_share` grouped by `arm` | not yet tested | WS1b, WS1f |
| 9.2 | Measured skill routinely exceeds the estimated ceiling; half of all cells show r at or above r_max | qualitative plus the table above | `frozen_2026-09/area_comparison_all.csv` and `admin1_arms.csv` | not yet tested | WS1 |
| 9.3 | Either the ceiling is biased low or the evaluation is optimistic; unresolved | qualitative | none | not yet tested | WS1a, WS1b, WS1c, WS1e |

## Section 11. The interpretation as it now stands

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 11.1 | Admin-2 is below the resolution these surveys can support | median district 6 to 36 measurements; median r_max 0.098 | `frozen_2026-09/admin1_arms.csv`, column `r_max` over 24 unique cells | not yet tested | WS1a, WS1b, WS4a |
| 11.2 | The constraint is a property of the target, not the predictors | 0.154 to 0.228, and education r 0.48 to 0.71 | `frozen_2026-09/individual_anchor.csv`; the education source is not located | not yet tested | WS3, WS4b |
| 11.3 | Geography carries most of what transportable signal exists | the spatial smoother matches the 294-predictor set; none survive FDR; penalised regression retains a median of zero | `frozen_2026-09/area_comparison_all.csv`, `bivariate_fdr.csv`, `penalized_retained.csv` | not yet tested | WS4a |
| 11.4 | Level and pattern are different problems with different answers | qualitative | Sections 4 and 6 | not yet tested | WS2, WS6 |
| 11.5 | The level must be resolved regionally to be worth anything | 0.164 to 0.413; national anchoring costs bias; an oracle level buys 0.23 pp | `frozen_2026-09/admin1_arms.csv`, `national_composition.csv` | not yet tested | WS2, WS5, WS6d |
| 11.6 | National estimates are the defensible deliverable; district maps are rankings | qualitative | Section 11 | not yet tested | WS8a, WS8b |
