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
| 3.1 | Individual-level prediction aggregated to Admin-2 | median r 0.516, r_share 0.92 | `frozen_2026-09/area_comparison_all.csv`, columns `pearson_r` and `r_share`, filter `approach == "Individual SL"` | **revised** | WS1, WS3a |
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
| 3.12 | Positive control, earth observation predicts district education | r 0.48 to 0.71 | no committed table located; Stage 0 searched `results/`, `sensitivity/`, `sandbox_parsimony/` and `docs/` | **revised** | WS4b |

**3.1 revised.** The stated 0.516 does not reproduce: the committed table gives
**0.524** and `r_share` **1.05** against the stated 0.92 (source:
`frozen_2026-09/area_comparison_all.csv`). More consequentially, Stage 0
established by reading the code what protocol produced it.
`aggregate_admin2_sl()` aggregates `res$yhat_full`, which
`src/analysis/sl_helpers.R` produces via
`origami::make_folds(cluster_ids = id_vec, V = folds)`: a **cluster-blocked
K-fold**, not the region-blocked leave-one-region-out that Section 5 uses. The
preprocessing recipe is also `prep()`ed on all rows before folds are formed.
The comment at `R/area_level_comparison.R:300` describing these rows as
"IN-SAMPLE" is inconsistent with that code and one of the two is wrong.

Measured on the smoke cell, the fold construction alone accounts for most of the
gap: Ghana `child_iron`, proxy arm, district, scores **0.398** under
region-blocked folds and **0.610** under cluster-blocked K-fold (source:
`results/tables/individual_arms_2026-09_SMOKE.csv`). The corrected layer's own
strict-protocol median is **0.031** (`protocol_reconciliation_medians.csv`,
`indiv_region_wt`). The full 2 by 2 across four cells is **not yet computed**;
see Part 3 of `PROJECT_STATUS_2026-09_UPDATE.md`.

**3.12 revised, and the general claim it supports is withdrawn.** Education is
predicted at **0.679 (Ghana), 0.795 (Sierra Leone) and 0.808 (Gambia)** for the
no-education indicator, above the stated 0.48 to 0.71 (source:
`results/tables/positive_control_targets.csv`). The specific claim therefore
holds and is now sourced. The **use** made of it does not: stunting reaches a
median of **-0.122**, improved water **0.014**, improved sanitation **0.007**
and full vaccination **0.164**, all well-measured district quantities. Across 81
DHS indicators the median achieved r is **0.071**
(`results/tables/reliability_skill_curve.csv`). The pipeline predicts
socioeconomic gradients, not health and nutrition outcomes, so a control built on
education does not show that a well-measured target would be predicted well.

**3.3, 3.4 and 3.5: WS4a supplies the resolution evidence.** The sweep does not
re-score the spatial smoother, so 3.3 and 3.4 remain `not yet tested`. Claim 3.5,
that Admin-1 covariates beat Admin-1 spatial 0.437 to 0.209, is not directly
re-measured either, but its premise is confirmed and sharpened: Admin-1 is the
best of four resolutions tested, with mean out-of-fold r **0.315** against
**0.086** at Admin-2 over the 14 cells present at all levels, best in **9 of 14**
cells and better than Admin-2 in **13 of 14**
(source: `results/tables/resolution_sweep.csv`, column `r_oof`).

**Open question 2 answered.** Section 13 asks whether a resolution between
Admin-1 and Admin-2 keeps covariates informative while measuring the target
adequately. Measured across Admin-1, two-way and three-way splits of each region,
and Admin-2, skill declines monotonically and **no intermediate level is a peak**.
An intermediate beats both endpoints in 4 of 14 cells and never decisively.

**A constraint recorded against 11.1.** Sierra Leone (4 Admin-1 regions) and
Gambia (6) cannot support a within-country covariate model at Admin-1 at all
under leave-one-region-out, which would train on 3 and 5 units. The resolution at
which the survey measures the target well is the resolution at which there are
too few units to learn a covariate map.

## Section 4. Where a district estimate's level comes from

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 4.1 | Admin-1 anchor, hard | mean r 0.413, median 0.405, MAE 8.85, mean absolute bias 1.59 | `frozen_2026-09/admin1_arms.csv`, filter `arm == "admin-2 fit + ADMIN-1 benchmark (hard)"` | **confirmed** | WS2 |
| 4.2 | Admin-1 anchor, shrunk | mean r 0.318, MAE 9.36, bias 2.76 | same file, filter `arm == "admin-2 fit + ADMIN-1 benchmark (shrunk)"` | **confirmed** | WS2 |
| 4.3 | Fit at Admin-1, extrapolate to Admin-2 | mean r 0.230, 22 cells, MAE 8.78, bias 3.18 | same file, filter `arm == "ADMIN-1 fit -> admin-2 (pooled)"` | **confirmed** | WS2 |
| 4.4 | National anchor | mean r 0.170, MAE 11.98, bias 5.85 | same file, filter `arm == "admin-2 fit + national benchmark"` | **confirmed** | WS2a, WS2c |
| 4.5 | No anchor | mean r 0.164, MAE 10.71, bias 3.24 | same file, filter `arm == "admin-2 fit (LORO), unbenchmarked"` | **confirmed** | WS2 |
| 4.6 | The hard Admin-1 anchor beats no anchor | 20 of 24 cells | same file, paired on `country` and `outcome` | **withdrawn** | WS2b, WS2d |
| 4.7 | The anchor is a one-parameter logit shift, so district ranking is identical in every arm | qualitative | `R/benchmark_area.R`, `benchmark_admin2_to_admin1()` | **confirmed** | WS2c |
| 4.8 | National anchoring makes mean absolute bias worse | 5.85 vs 3.24 pp | `frozen_2026-09/admin1_arms.csv`, column `bias_pp` | **revised** | WS2c, withheld pending the implied-shift audit |
| 4.9 | Hard beats shrunk despite the theoretical risk | 0.413 vs 0.318 | `frozen_2026-09/admin1_arms.csv` | not yet tested | WS2d |
| 4.10 | The anchoring gain is not circular | qualitative argument in Section 4.4 | none | **withdrawn** | WS2a, WS2b |
| 4.11 | The survey supplies the level; the model supplies the within-region pattern | qualitative | Sections 4 and 6 | **revised** | WS2d |

**4.1 to 4.5 confirmed as reproductions.** All five published arms reproduce from
`anchor_controls.csv` to within 0.004 in mean r (hard anchor 0.409 against 0.413,
MAE 8.91 against 8.85), the difference being rounding before averaging. The
numbers are right; what they mean is not what Section 4 says.

**4.6 withdrawn.** The hard anchor beats no anchor in 22 of 24 cells here, close
to the published 20 of 24. Under a jackknife anchor, where each district's
regional target is computed from the region's other districts so no district
contributes to its own correction, it beats no anchor in **8 of 24** cells, mean
r falls from **0.409 to 0.147** against **0.156** unanchored, and both MAE and
absolute bias become worse (source: `results/tables/anchor_controls.csv`, arm
`5 ADMIN-1 anchor (hard, JACKKNIFE)`). The gain is an artifact of a district
contributing roughly a quarter to a third of the number used to correct it.

**4.10 withdrawn.** Section 4.4's argument that the anchor is not circular rests
on the national anchor as a counterfactual. That control is not valid: a national
anchor built from about 1,000 respondents contains roughly 0.1 percent of any one
district's data, against a quarter to a third for a regional anchor, so its
failure to gain is what a leakage account predicts. The valid control is an arm at
the same resolution using no covariates, and it was not run. It is now: see 4.12.

**4.12 WITHDRAWN as originally stated, 2026-09-21.** The original entry read:
assigning each district its region's design-based survey estimate, with no
covariates at all, reaches mean r 0.516, MAE 7.38 pp and mean absolute bias 0.77
pp, and is better than every covariate arm on every metric.

**That arm had not been jackknifed.** Its regional estimate was computed from all
the region's respondents including the scored district's, which is the same
mechanism that withdrew 4.6. Under a symmetric jackknife (source:
`results/tables/anchor_controls_B1.csv`, arm
`2b flat REGIONAL mean (JACKKNIFE)`):

| Arm | mean r | MAE pp | mean absolute bias pp |
|:---|---:|---:|---:|
| Flat regional mean, as published | 0.516 | 7.38 | 0.77 |
| **Flat regional mean, jackknifed** | **0.076** | **10.87** | **3.41** |
| No anchor (LORO) | 0.156 | 10.77 | 3.21 |
| Admin-1 anchor (hard), jackknifed | 0.147 | 12.33 | 4.35 |

**The corrected claim.** Under a matched out-of-sample control the covariate-free
regional mean does **not** beat the covariate model: it is lower on correlation
(0.076 against 0.156) and indistinguishable on error (10.87 against 10.77 pp).
The Section 4 interpretation that "the covariate pattern subtracts from it" is
withdrawn. What survives from WS2 is narrower and still substantial: the
published anchoring gain is circular, and no arm tested reaches a useful
district-level correlation.

**Which number applies depends on the use case, and both are reported.** For a
district the survey did measure, the un-jackknifed regional mean is what one
would actually compute, but the district's own direct estimate is available and
better. For a district the survey did not measure, which is the deployment case,
the jackknifed number is the honest analogue. WS-B2 settles the choice against
simulated truth rather than against either observed number.

**4.8 released, with its explanation changed.** National anchoring does make
absolute bias worse, by 2.64 pp against no anchor, better on MAE in only 4 of 24
cells. The stated reason, that a single number displaces districts that were
already correct, is not the mechanism. WS2c measured the un-anchored
population-weighted national aggregate against the design-based survey estimate
and found a mean absolute gap of **9.60 pp** over 24 cells, reaching **77.57 pp**
for Sierra Leone child vitamin A, where the model predicts 89.6 percent against a
survey 12.0 percent (source: `results/tables/anchor_implied_shifts.csv`, rows
with `arm == "national"`). The base model's level error varies by region, and a
single national shift cannot correct a region-varying error.

**4.9 not yet tested, and the question is not live.** Hard beats shrunk under the
published anchors. Under the jackknife the hard anchor does not beat no anchor, so
there is no anchoring scheme yet shown to work out of sample for the two variants
to be compared within.

**4.11 revised.** The survey supplies the level, and on this evidence it also
supplies more of the pattern than the covariates do, since the flat regional mean
outperforms the covariate model on correlation as well as on error.

## Section 5. The individual-level anchor

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 5.1 | District, proxy arm | 16 cells, mean r 0.154, median 0.126, MAE 8.33 | `frozen_2026-09/individual_anchor.csv`, filter `unit == "district"` and `arm == "proxy"` | **revised** | WS3a |
| 5.2 | District, questionnaire arm | mean r 0.228, median 0.264, MAE 7.92 | same file, filter `unit == "district"` and `arm == "quest"` | **revised** | WS3a, WS3b |
| 5.3 | Cluster, proxy and questionnaire arms | 0.146 and 0.229 | same file, filter `unit == "cluster"` | **revised** | WS3f |
| 5.4 | Questionnaire better in 10 of 16 cells, mean gain +0.075 | 10 of 16, +0.075 | same file, paired on `country`, `outcome` and `unit` | **revised** | WS3b, WS3c |
| 5.5 | Clears r = 0.4 in 3 of 16 cells | 3 of 16 | same file, column `r` | **revised** | WS3a |
| 5.6 | In Malawi the questionnaire adds nothing | gains 0.000, 0.000, 0.002, 0.004 | same file, filter `country == "Malawi"` | **withdrawn** | WS3b |
| 5.7 | Maximum r anywhere across all 64 rows | 0.544 | same file, maximum of column `r` | **revised** | WS7a |
| 5.8 | The null is not explained by bad linkage, the inflammation adjustment, or overfitting | qualitative | Section 5.3 | **withdrawn** | WS3a, WS4b |
| 5.9 | Cluster linkage does not help, and helps least in Sierra Leone | Gambia +0.017, Ghana -0.003, Malawi -0.002, Sierra Leone -0.025 | same file, `r` at `unit == "cluster"` minus `r` at `unit == "district"` | **confirmed** | WS3f |

**5.1, 5.2 and 5.9: the four-cell re-run has now completed.** The published
numbers were computed on a **contaminated and non-nested** pair of arms and
should not be quoted:

- **Contaminated.** The questionnaire arm saw `gw_wm_whbc` and `gw_gchb`, two
  Ghana haemoglobin measurements, plus thirteen further blood-derived columns
  found from Stata labels (WS7a, WS7b).
- **Non-nested.** The published questionnaire arm applies
  `is_biomarker_column()` to `Xvars_full`, which also strips the MAP sickle-cell
  rasters and the DHS mean-haemoglobin Admin-2 aggregates that the proxy arm
  keeps, because the proxy arm uses `Xvars` unfiltered. The two arms therefore
  differed in more than the questionnaire, so the gap between them is not
  attributable to the questionnaire alone. `allowed_under_arm()` fixes this by
  scoping the filter to the concurrent survey.

**Measured on four cells under the strict protocol** (source:
`results/tables/individual_arms_2026-09_PARITY.csv`, `unit == "district"`,
`protocol == "region_loro"`): proxy mean r **0.274**, questionnaire **0.261**,
questionnaire plus field haemoglobin **0.342**. The questionnaire gain is
**-0.013**, better in 2 of 4 cells, against a published **+0.075**. Section 5's
conclusion that the questionnaire barely beats geospatial proxies is
**strengthened**: on these cells it does not beat them.

**5.9 confirmed.** Cluster minus district gains reproduce the published ordering:
Gambia **+0.018** (published +0.017), Ghana **-0.008** (-0.003), Malawi
**-0.001** (-0.002), Sierra Leone **-0.034** (-0.025). Sierra Leone is again the
worst, which is the falsification Section 5.4 reports.

**5.10 new.** Field haemoglobin, the deployable DHS-standard scenario, adds
**+0.394** for Gambia women's iron (taking that cell to **0.848**) and **+0.027**
for Ghana child iron, while adding nothing for Malawi child vitamin A and costing
**0.150** for Sierra Leone women's iron. Mean **+0.068** over four cells. Its
district MAE is **19.96 pp** against **10.12 pp** for the proxy arm, so it gains
on correlation and loses on level. The pattern is outcome-specific in the
direction physiology predicts.

**A protocol finding that belongs to Section 3 as much as Section 5.** The same
four cells score **+0.319** higher under cluster-blocked K-fold than under
region-blocked folds for the proxy arm, **+0.348** for the questionnaire arm and
**+0.286** for the haemoglobin arm. Fold construction alone accounts for most of
the distance between Section 3's 0.516 and Section 5's 0.154.

**5.4 revised, on the published table's own arithmetic.** The published mean gain
of +0.075 and count of 10 of 16 reproduce exactly. Under the exclusions Section
5.5 asks for and did not perform:

| Subset | Cells | Mean gain | Median gain | Questionnaire better |
|:---|---:|---:|---:|:---|
| As published | 16 | +0.0748 | +0.0080 | 10 of 16 |
| Excluding Malawi | 12 | +0.0993 | +0.0220 | 8 of 12 |
| Excluding Ghana women_iron | 15 | **+0.0375** | **+0.0040** | 9 of 15 |
| Excluding both | 11 | +0.0506 | +0.0200 | 7 of 11 |

(source: `frozen_2026-09/individual_anchor.csv`, `unit == "district"`.) Dropping
the one cell Section 5.5 names halves the mean gain and takes the median to
0.004. The conclusion that the questionnaire adds little is **strengthened**.

**5.6 withdrawn as a finding about Malawi.** Malawi's `GW` domain contains
**zero** columns: `Xvars_full` equals `Xvars` at 1222 columns, against 378 `gw_`
columns for Ghana (measured from the `outcome_data_*` targets). Its questionnaire
arm has **fewer** predictors than its proxy arm in the published table, 1141
against 1147, because the two differ only by the `sd > 0` filter after
imputation. Independently, the WS7a leakage report finds Malawi's proxy and
questionnaire maxima identical to four decimal places in all eight outcomes. The
gains of 0.000, 0.000, 0.002 and 0.004 are the signature of one arm scored twice.
Malawi's questionnaire data exists (242 columns in
`data/IPD/Malawi/Malawi_merged_dataset.rds`, coded `m01` to `m2xx`) but carries
zero labels and no local codebook, so ingestion is held as a to-do.

**5.8 withdrawn.** Section 5.3 excludes three explanations for the null and
concludes the signal is weak "for any predictor set that can be constructed".
WS4b measures the same pipeline on 81 DHS indicators at a median achieved r of
**0.071** against a median reliability of **0.784**, and the micronutrient cells
sit **above** that line. The null is not specific to micronutrient targets, so
it is not evidence about the target. It is evidence about the predictors.

## Section 6. A predicted national level cannot substitute for a survey

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 6.1 | The VMNIS LOCO model is competent | Vitamin A preschool random forest, MAE 11.75, Pearson 0.655 | `frozen_2026-09/national_vmnis_loco.csv` | not yet tested | WS6c |
| 6.2 | The null arm's negative Pearson is an artefact, not a finding | qualitative | `frozen_2026-09/national_vmnis_loco.csv`, filter `model == "null"` | not yet tested | WS6b |
| 6.3 | Vitamin A preschool ceiling components | sd country 1.411, method 0.564, residual 0.000, sampling 0.816 | `frozen_2026-09/national_vmnis_ceiling.csv`, row 1 | **revised** | WS6a |
| 6.4 | Vitamin A preschool ceiling and saturation | r_max 0.818, r_share 0.80 | same file, columns `r_max_report` and `r_share` | **revised** | WS6a |
| 6.5 | Correction to the record: the model is not 98 percent saturated at a ceiling of 0.66 | 80 percent of 0.818 | same file | **revised** | WS6a |
| 6.6 | sd_resid hits the boundary at 0.000 in two panels, so r_max is untrustworthy there | 2 of 4 panels | same file, column `sd_resid` | **withdrawn** | WS6b |
| 6.7 | For Vitamin A NPW, method variance exceeds country variance | 1.996 vs 1.232 | same file, row 3 | **confirmed** | WS6b |
| 6.8 | Composition arm errors | MAE 5.81, 8.22, 12.70, 5.58 pp | `frozen_2026-09/national_composition.csv`, column `mae_pp` grouped by `arm` | **revised** | WS6d |
| 6.9 | Composition arm absolute bias | 3.35, 4.69, 10.03, 0.73 pp | same file, column `bias_pp` | **revised** | WS6d |
| 6.10 | The predicted national level is on a different scale; the null beats the covariate model | 6 of 8 cells | `frozen_2026-09/national_composition_levels.csv` | not yet tested | WS6c |
| 6.11 | Even a perfect national level buys 0.23 pp, better in only 4 of 8 cells | 5.81 to 5.58 pp | `frozen_2026-09/national_composition.csv` | **revised** | WS6d |
| 6.12 | Scope limitation: the VMNIS and transport outcomes intersect on vitamin A only | 2 outcomes times 4 countries, 8 cells | Section 6.5 | **confirmed** | WS6d |

**6.3 revised.** The sampling term is wrong, and by a factor of 4.7. For Vitamin
A / preschool the published `sd_sampling` of 0.816 is the square root of the
**arithmetic mean of a reciprocal**, which the smallest surveys dominate: the
panel holds 34 surveys under n = 50 and a minimum of n = 8, against a median of
373.5, and 37 surveys with prevalence under 2 percent, where the delta method is
clamped at 0.005 and a single survey can contribute up to 301/(n-1). Using the
median survey's variance instead gives **0.172**, and the implied effective
sample size rises from **13** to a plausible value (source:
`results/tables/vmnis_sampling_audit.csv`, columns `sd_samp_from_mean`,
`sd_samp_from_median`, `implied_n_from_mean`). Revised components for that panel:
sd country 1.436, method 0.383, residual **0.703**, sampling **0.171** (source:
`results/tables/national_vmnis_ceiling_revised.csv`, `version == "revised"`).

**6.4 and 6.5 revised.** The corrected ceiling for Vitamin A / preschool is
`r_max_report` **0.869** and `r_max_standardised` **0.893**, against the
published 0.818 and 0.866. Against a best model r of 0.655 the saturation is
therefore about **75 percent**, not 80. The direction of the Section 6.5
correction stands and its size grows: there is more headroom, not less.

**6.6 withdrawn, and its stated cause was wrong.** Two panels report `sd_resid`
exactly 0.000, and Section 6.3 attributes this to a degenerate `lmer` fit. It is
not. The raw residual variance from `lmer` is non-zero in all four panels; the
reported zero is produced by **subtracting an over-large sampling term** from it
and flooring at zero. With the corrected sampling term no panel sits at the
boundary and all four are usable (source:
`results/tables/national_vmnis_ceiling_revised.csv`, columns
`resid_at_boundary`, `sampling_exceeds_resid`, `usable`). No refit with priors
was required, so the `blme` dependency the workstream anticipated is not needed.

**6.7 confirmed.** For Vitamin A / NPW, method variance still exceeds country
variance after the correction: 2.101 against 1.173 (published 1.996 against
1.232). That panel is measuring survey methodology more than country.

**6.8, 6.9 and 6.11 revised in presentation, confirmed in substance.** The arm
MAEs reproduce exactly (5.81, 8.22, 12.70, 5.58). Reported with signed bias
separated from absolute bias, the transported arm's signed bias is **-3.14 pp**
and the oracle's **+0.35 pp** (source:
`results/tables/national_composition_revised.csv`). Four of the eight cells are
near-degenerate: the women's vitamin A cells sit at 1.3 to 2.5 pp true national
prevalence, so there is almost no level to get wrong. Excluding them, the oracle
takes MAE from **9.10 to 8.62 pp**, a gain of 0.48 pp, while taking signed bias
from **-5.31 to +0.98 pp**. The correct statement is that **an oracle national
level removes almost all of the bias and almost none of the error**: the residual
error is pattern, and no national quantity can reach it.

## Section 9. r_share above one

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 9.1 | Mean r_share by anchoring arm | hard 2.06, shrunk 1.64, Admin-1 fit 1.77, national 1.34, none 1.35 | `frozen_2026-09/admin1_arms.csv`, column `r_share` grouped by `arm` | **revised** | WS1b, WS1f |
| 9.2 | Measured skill routinely exceeds the estimated ceiling; half of all cells show r at or above r_max | qualitative plus the table above | `frozen_2026-09/area_comparison_all.csv` and `admin1_arms.csv` | **revised** | WS1 |
| 9.3 | Either the ceiling is biased low or the evaluation is optimistic; unresolved | qualitative | none | **revised** | WS1a, WS1b, WS1c, WS1e |

**9.1 revised.** Under the empirical ceiling and with the production
`r_max > 0.05` guard applied consistently, the arm means become hard **0.75**,
shrunk **0.56**, Admin-1 fit **0.38**, national **0.35**, none **0.33**, and the
medians hard **0.68**, shrunk **0.51**, Admin-1 fit **0.52**, national **0.32**,
none **0.34** (source: `results/tables/r_share_revised_summary.csv`, columns
`mean_share_empirical` and `med_share_empirical`). Two separate corrections
produce the change, and both are needed: applying the `r_max > 0.05` guard that
`add_reliability_columns()` already uses takes the hard arm from 2.06 to 1.17,
and replacing the biased ceiling takes it from 1.17 to 0.75. No arm exceeds 1.

**9.2 revised.** Cells with `r_share > 1` fall from **25 to 7** of 118 arm-cell
rows (source: `results/tables/r_share_revised.csv`, columns
`r_share_analytic` and `r_share_empirical`). The residual seven are concentrated
in the anchored arms, which is consistent with the mechanism WS2 tests and is
not evidence that the empirical ceiling is also biased.

**9.3 revised, and resolved to the first horn.** **The ceiling was biased low.**
WS1c simulates outcomes over the real survey structure with the truth known by
construction: the analytic estimator's mean bias against the attainable
correlation is **-0.161**, the split-half estimator's is **+0.007**
(source: `results/tables/reliability_simulation.csv`, columns `bias_analytic`
and `bias_empirical`). On the real data the analytic median `r_max` is **0.132**
against an empirical **0.613**, with the empirical exceeding the analytic in
**21 of 24** cells (source:
`results/tables/reliability_analytic_vs_empirical.csv`). The design effect that
would reconcile the two has median **0.969** against the assumed **1.5**. This
does not exonerate the anchored arms, whose evaluation may be optimistic for the
separate reason set out in `docs/TARGET_ESTIMAND.md` section 3; WS2 tests that.

## Section 11. The interpretation as it now stands

| ID | Claim | As stated | Source map | Status | Resolved by |
|:---|:---|:---|:---|:---|:---|
| 11.1 | Admin-2 is below the resolution these surveys can support | median district 6 to 36 measurements; median r_max 0.098 | `frozen_2026-09/admin1_arms.csv`, column `r_max` over 24 unique cells | **revised** | WS1a, WS1b, WS4a |
| 11.2 | The constraint is a property of the target, not the predictors | 0.154 to 0.228, and education r 0.48 to 0.71 | `frozen_2026-09/individual_anchor.csv`; the education source is not located | **revised** | WS3, WS4b |
| 11.3 | Geography carries most of what transportable signal exists | the spatial smoother matches the 294-predictor set; none survive FDR; penalised regression retains a median of zero | `frozen_2026-09/area_comparison_all.csv`, `bivariate_fdr.csv`, `penalized_retained.csv` | not yet tested | WS4a |
| 11.4 | Level and pattern are different problems with different answers | qualitative | Sections 4 and 6 | **revised** | WS2, WS6 |
| 11.5 | The level must be resolved regionally to be worth anything | 0.164 to 0.413; national anchoring costs bias; an oracle level buys 0.23 pp | `frozen_2026-09/admin1_arms.csv`, `national_composition.csv` | **revised** | WS2, WS5, WS6d |
| 11.6 | National estimates are the defensible deliverable; district maps are rankings | qualitative | Section 11 | not yet tested | WS8a, WS8b |

**11.1 revised, and the direction of the correction matters.** The stated median
`r_max` of 0.098 does not reproduce: the committed table gives **0.1315**
(source: `frozen_2026-09/admin1_arms.csv`, column `r_max`, 24 unique cells), and
the measured empirical ceiling gives **0.6132** (source:
`results/tables/reliability_analytic_vs_empirical.csv`, column `r_max_emp`).
The accompanying claim that 4 of 24 cells have no detectable signal is also
wrong in both directions: the analytic estimator returns exactly zero in **10**
of 24 cells, and the empirical one in **3**.

The conclusion that Admin-2 is below the resolution these surveys support is
**not** withdrawn, because it does not rest on the ceiling alone: median
district sample sizes of 6 to 36, and 62 of 75 Ghana districts holding a single
cluster, are direct measurements. What changes is the size of the headroom. A
median ceiling near 0.61 rather than 0.10 means the gap between what models
achieve (median r 0.164 unanchored) and what is attainable is **much larger**
than reported, so the reported skill is a smaller fraction of the possible, not
a larger one. The claim that models are near-saturated at Admin-2 does not
survive; the claim that Admin-2 is hard does.

## Session addendum, signal-and-shipping, 2026-09-21

**3.9 revised.** The claim is that 0 of 294 predictors survive FDR control in all
24 cells. That reproduces and is not disputed, but it answers a cell-by-cell
question. Pooled across cells with country as the clustering unit and calibrated
by permutation, **two predictors survive a family-wise correction**:
`soilgrids_ph` (meta z 7.34, family-wise permutation p 0.020) and
`dhs_w_skilled_attendant` (z 7.02, p 0.028), both in the region-partialed family
over 12 cells and 3 countries (source:
`results/tables/predictor_consistency_meta.csv`, columns `z` and `p_perm_fwer`,
filter `family == "region_partialed" & cellset == "shared"`). Nothing survives
marginally. The correct statement is that the cell-by-cell screen was
underpowered, and that pooling recovers two associations rather than none.

**11.3 revised.** The claim is that geography carries most of what transportable
signal exists, supported by the spatial smoother, the FDR screen and penalised
regression retaining a median of zero. The pooled scan sharpens it. What survives
after removing the between-region component is **soil pH** and a cluster of
health-service access measures (skilled attendance z 7.02, antenatal visits 5.76,
improved water 4.19, improved sanitation 4.09). Soil pH is mechanistically
interpretable with the expected sign, since zinc and iron bioavailability falls
as pH rises. Soil zinc itself carries the right sign (z -3.37, more zinc with
less deficiency) and does not survive (family-wise p 0.948). The revised claim is
that the transportable signal is a socioeconomic access gradient plus one soil
chemistry variable, and that it is not nutrition-proximal in the sense the
mechanism hypothesis requires.

**New row A.1.** The analytic random-effects q-values are anti-conservative at
four clusters: they identify up to 22 survivors where the permutation supports at
most two, with the smallest analytic q at 5.9e-11 for a predictor whose
family-wise permutation p is 0.020. Any future pooled analysis in this project
must carry a permutation calibration.

**New row A.2.** Benjamini-Hochberg applied to permutation p-values is
**resolution-limited by design** here. With 500 permutations the smallest
attainable per-predictor p is 0.002 and only one to three predictors per family
reach it, so the smallest attainable BH q is 0.187. The "zero survive under
permutation BH" result is therefore not evidence of absence, and the family-wise
max-statistic calibration is the sound one to quote.

**New row C3.1.** The pre-specified nine-covariate iron-anaemia bridge does not
beat the full 294 within country (mean r 0.027 against 0.364) or across
countries (-0.022 against 0.210), and does not beat a covariate-free jackknifed
regional mean (0.181) (source: `results/tables/iron_bridge.csv`). The two
countries with an exact district crosswalk disagree: Gambia women's iron reaches
0.459 under LOCO with nine covariates against 0.467 with 294, and Sierra Leone is
negative in all four of its cells. Ghana and Malawi discard roughly two thirds of
their prior-round DHS districts in the crosswalk, which weakens the test there
and does not explain Sierra Leone.

**New row D.1.** Six of 24 committed cells fail a calibration gate set at the
larger of 10 pp and twice the survey CI half-width, led by Sierra Leone child
vitamin A at a 77.57 pp gap (source:
`results/tables/calibration_gate_report.csv`). Register row 3.8's claim that
national prevalence lands inside the survey CI in 24 of 24 cells belongs to the
**individual-level aggregate under its own protocol** and not to the area-level
ridge, whose mean absolute national gap is 9.60 pp. Row 3.8 is scoped by
estimator accordingly.

**4.12 further revised.** See the entry above, and the dated correction in
`docs/findings/TWO_READINGS.md`.
