# Project status update, September 2026

Branch `accuracy-and-impact-2026-09`, cut from `main` at `f8b23d9`.

This note separates **changed results** from **new capabilities**. Every number
carries a source map. A number that could not be produced by running code on
this branch is written "not yet computed". The manuscript is not edited by this
branch; `docs/findings/CLAIMS_REGISTER.md` is the input to that later edit.

---

## Part 1. Changed results

### 1.1 A live leak in the arm that carries the central interpretation

`gw_wm_whbc` and `gw_gchb` are Ghana haemoglobin measurements that
`is_biomarker_column()` did not block, because `grepl("Hb", cols)` is
case-sensitive. Measured: `gw_wm_whbc` has n 981, range 7.0 to 17.0 g/dL, median
12.9, and a mean of 12.98 in non-deficient women against 10.71 in iron-deficient
women; it correlated **0.500** with the outcome, the highest of any eligible
column across all 48 cell-and-set combinations
(source: `results/tables/leakage_report.csv`).

A further **thirteen** columns, **fifteen in total on this branch**, were found by deriving the guard from Stata
variable labels: four CRP-and-AGP inflammation flags, two anaemia-malaria-
inflammation statuses, six HbA1c recodes, and a phlebotomist ID
(source: `results/tables/assay_guard_label_comparison.csv`). None names an
analyte.

**Scope.** The primary models use `Xvars`, which excludes the survey domain
entirely, so **no production estimate was affected**. Section 5's questionnaire
arm was.

**Manuscript sentences affected:** Section 5.2 table, Section 5's "10 of 16" and
"+0.075", Section 5's "max r anywhere 0.544", Section 7.3's guard description,
Section 7.4's preserved-column list, Section 7.5's residual-risk paragraph.
Register rows 5.1 to 5.7, 5.9.

### 1.2 The noise ceiling was biased low, and Section 9 is resolved

| Quantity | Published | Measured on this branch |
|:---|---:|---:|
| Median `r_max` across 24 cells | 0.098 | **0.613** (empirical) / 0.132 (analytic, live table) |
| Cells with no detectable signal | 4 of 24 | **3 of 24** empirically, 10 of 24 analytically |
| Mean `r_share`, hard anchor | 2.06 | **0.75** |
| Cells with `r_share > 1` | 25 of 118 | **7 of 118** |

The analytic estimator understates the attainable correlation by a mean of
**0.161**; the split-half estimator by **0.007**
(source: `results/tables/reliability_simulation.csv`). Two causes: a design
effect fixed at 1.5 where the reconciling value has median **0.969**, and a
`max(0, .)` truncation.

**The direction matters.** A higher ceiling means the gap between achieved and
attainable is **larger** than reported. "Models are near-saturated at Admin-2"
does not survive. "Admin-2 is hard" does.

**Manuscript sentences affected:** Section 2.4 in full, Section 9 in full,
Section 11 claim 1, and every reported `r_share`. Register rows 9.1 to 9.3, 11.1.

### 1.3 The anchoring gain does not survive its own control

| Arm | Mean r | MAE pp | Mean absolute bias pp |
|:---|---:|---:|---:|
| No anchor | 0.156 | 10.77 | 3.21 |
| Admin-1 anchor (hard), as published | 0.409 | 8.91 | 1.61 |
| **Admin-1 anchor (hard), jackknife** | **0.147** | **12.33** | **4.35** |
| **Flat regional mean, no covariates** | **0.516** | **7.38** | **0.77** |

(source: `results/tables/anchor_controls.csv`.)

Under a jackknife anchor, where each district's regional target is computed from
the region's other districts, the gain over an unanchored model goes from +0.253
(better in 22 of 24 cells) to **-0.001** (better in 8 of 24). And a predictor
with **no covariates at all** beats every covariate arm on every metric.

The implied-shift audit also found the base model mis-levelled by a mean of
**9.60 pp**, reaching **77.57 pp** for Sierra Leone child vitamin A
(source: `results/tables/anchor_implied_shifts.csv`).

**Manuscript sentences affected:** Section 4.1 table interpretation, Section 4.2,
Section 4.3, Section 4.4 in full, Section 4.5, Section 11 claims 4 and 5.
Register rows 4.6, 4.8 to 4.12, 11.4, 11.5.

### 1.4 The positive control does not license the claim it is used for

Education is predicted at **0.679 to 0.808** across three countries, above the
0.48 to 0.71 the review document states without a locatable source. But stunting
reaches a median of **-0.122**, improved water **0.014** and improved sanitation
**0.007**, all well measured
(source: `results/tables/positive_control_targets.csv`).

Across 81 DHS indicators with a median design-based reliability of 0.784, the
same pipeline achieves a median out-of-fold correlation of **0.071**, while the
24 micronutrient cells reach **0.200** and sit **+0.219 above** the DHS-indicator
line (source: `results/tables/reliability_skill_curve.csv`). The correlation
between reliability and achieved skill across all 105 targets is **0.043**.

**This reverses Section 11 claim 2.** Well-measured district quantities are not
predicted well either, and the micronutrient cells are among the better
performers. The covariates track socioeconomic gradients, not health and
nutrition outcomes.

**Manuscript sentences affected:** Section 3's positive-control row, Section 5.3,
Section 11 claim 2. Register rows 3.12, 5.8, 11.2.

### 1.5 The VMNIS sampling term was wrong by a factor of 4.7

`sd_sampling` 0.816 for Vitamin A preschool is the square root of the
**arithmetic mean of a reciprocal**, which 34 surveys under n = 50 dominate.
Using the median gives **0.172**
(source: `results/tables/vmnis_sampling_audit.csv`). Corrected, `r_max_report`
rises from 0.818 to **0.869**, saturation falls from 80 to about **75 percent**,
and **`sd_resid` at the boundary was not a degenerate `lmer` fit** but the
over-large sampling term being subtracted from a healthy residual. All four
panels become usable.

**Manuscript sentences affected:** Section 6.3 table and both reviewer flags,
Section 6.4's oracle framing, Section 10's VMNIS row, Section 11 claim 5.
Register rows 6.3 to 6.9, 6.11.

### 1.6 The questionnaire gain reverses sign once the arms are clean

Measured on four cells under the strict protocol, with the guard corrected and
the arms nested (source: `results/tables/individual_arms_2026-09_PARITY.csv`,
`unit == "district"`, `protocol == "region_loro"`):

| Arm | Mean r | Gain over proxy |
|:---|---:|---:|
| proxy | 0.274 | |
| questionnaire | 0.261 | **-0.013** |
| questionnaire + field haemoglobin | 0.342 | +0.068 |

The published gain is **+0.075**. Two changes account for the reversal: the
questionnaire arm no longer sees haemoglobin, and the two arms are now nested,
where before the questionnaire arm also lost the MAP sickle-cell rasters and DHS
mean-haemoglobin aggregates that the proxy arm kept.

**Section 5's conclusion is strengthened.** A household questionnaire
administered to the same people does not beat geospatial proxies at all on these
cells.

**Field haemoglobin is the exception, and it is outcome-specific.** It adds
**+0.394** for Gambia women's iron, taking that cell to **0.848**, and nothing
for vitamin A. Its district MAE is 19.96 pp against 10.12 pp for the proxy arm,
so it gains on correlation and loses on level.

**The protocol effect is the larger story.** The same four cells score
**+0.319** higher under cluster-blocked K-fold than under region-blocked folds
for the proxy arm. Fold construction alone accounts for most of the distance
between Section 3's 0.516 and Section 5's 0.154, and Section 3's number should
not be quoted as an out-of-sample individual-level result without saying so.

**Manuscript sentences affected:** Section 3's individual-level row, Section 5.2
table, Section 5.5. Register rows 3.1, 5.1, 5.2, 5.10.

### 1.7 There is no useful resolution between Admin-1 and Admin-2

WS4a swept four levels with the learner, folds and covariates held constant.
Over the 14 cells present at all four, mean out-of-fold r is **0.315** at
Admin-1, **0.166** at a two-way split of each region, **0.093** at a three-way
split and **0.086** at Admin-2 (source: `results/tables/resolution_sweep.csv`).
Admin-1 is the best level in **9 of 14** cells and beats Admin-2 in **13 of 14**;
an intermediate level beats both endpoints in only 4 of 14 and never decisively.

**Open question 2 is answered in the negative.** Section 13 suggests the
crossover lies between Admin-1 and Admin-2. Measured, skill declines
monotonically with resolution and there is no interior peak.

**A structural constraint the sweep surfaced.** Sierra Leone (4 regions) and
Gambia (6) cannot support a within-country covariate model at Admin-1 at all,
because leave-one-region-out would train on 3 and 5 units. Where the survey
supports the estimate there are too few units to learn a map; where there are
enough units the estimate is too noisy. This is the same tension WS2 reports from
the other side, where a covariate-free regional mean beats every covariate arm.

**Manuscript sentences affected:** Section 13 open question 2, Section 3's
Admin-1 comparison row, Section 11 claim 1. Register rows 3.5, 11.1.

### 1.8 Section 5's published questionnaire advantage is substantially one cell

Recomputed from the published table: dropping Ghana women_iron takes the mean
gain from **+0.075 to +0.038** and the median to **+0.004**, with the
questionnaire better in 9 of 15 cells. Malawi's questionnaire arm has **fewer**
predictors than its proxy arm (1141 against 1147) because Malawi has **no `gw_`
block at all**; its four cells cannot show a gain by construction. Excluding
both, 11 cells, mean gain **+0.051**.

Section 5.5 says this check "was not formally checked". It now is, and the
conclusion that the questionnaire adds little is **strengthened**.

**Manuscript sentences affected:** Section 5.2's counts, Section 5.5's fourth
bullet. Register rows 5.4, 5.6.

---

## Part 2. New capabilities

| Capability | Where | What it does |
|:---|:---|:---|
| Admin-2 join lint | `R/lint_admin2_joins.R`, `tests/testthat/test-admin2-join-lint.R` | Ratchet over 69 recorded name-only joins; fails on any new one |
| Unit-count assertions | `R/unit_counts.R`, `tests/testthat/test-unit-counts.R` | One-sided check against counts derived from the store; 0 over-counts in 24 checks |
| Leakage report | `R/leakage_report.R`, target `leakage_report` | Ranks every predictor against the outcome per cell; 24 incoming edges, one per cell |
| Empirical reliability | `R/reliability_empirical.R` | Split-half ceiling with a design-faithful variant |
| Ceiling simulation | `R/reliability_simulation.R` | Estimator bias against a known truth |
| Three-class assay guard | `R/data_prep.R`, `biomarker_column_class()`, `allowed_under_arm()` | `blood_draw` / `hb_field` / `questionnaire`, scoped so arms nest |
| Label-derived guard | `scripts/accuracy_impact/ws7b_label_guard.R` | Ground truth for 865 pipeline columns from Stata labels |
| WHO threshold config | `config/who_thresholds.csv` | Centralised from three hard-coded copies, provenance recorded, all citations marked unverified |
| Fitness-for-use rubric | `docs/FITNESS_FOR_USE.md` | Source-mapped thresholds for what to report and what to suppress |
| Evaluation harness | `harness/` | Folds, locked test cells, scorer that refuses leaky submissions |
| Regression gate | `scripts/accuracy_impact/ws9_regression_gate.R` | 19 of 19 frozen tables unchanged |
| Target estimand note | `docs/TARGET_ESTIMAND.md` | Which ceiling bounds which estimator |

Test suite: **37 tests, 0 failures**, via `Rscript tests/testthat.R`.

---

## Part 3. Not yet computed

| Item | Why |
|:---|:---|
| WS3a's production-stack arm | The 4-cell parity subset ran under the light NNLS stack for both protocols. The learner half of the 2 by 2 is not computed, so the PROTOCOL effect is bounded and the LEARNER effect is not. |
| WS3f cluster-train, district-evaluate | Both variants scored here differ only in the unit of aggregation, matching the published script. Training at the cluster is not computed. |
| WS3d, permutation importance for Ghana women_iron | Not reached. |
| WS6c, method covariates in the LOCO VMNIS model | Not reached. The variance decomposition already shows method variance exceeding country variance for Vitamin A NPW, 2.101 against 1.173. |
| Malawi questionnaire ingestion | Blocked on documentation. 242 columns coded `m01`/`m115a`/`m220h`, zero labels, no codebook locally; the README points to `immpact@cdc.gov`. |
| Côte d'Ivoire as a scored test country | It has rasters and `oos_civ_*` predictions but no biomarker labels. Exported unscored. |

---

## Part 4. Deviations from the plan as approved

1. **WS1, WS2 and WS6 landed in one commit** rather than one per workstream. Not
   corrected, because guardrail 2 forbids history rewriting. The findings files
   remain one per workstream.
2. **One commit message was amended** to correct "nine" to "thirteen" in its
   subject line, where the body already said thirteen. Amending is a form of
   history rewriting; it was done because guardrail 3 forbids a wrong number in a
   commit message and the commit was unpushed. Recorded here rather than left
   silent.
3. **The WS7a lint ships as a ratchet, not a ban.** There are 69 pre-existing
   name-only joins; a test failing on all of them would be switched off.
4. **The WS8c held-out set is 30 percent of regions per country, not one
   region.** One region returned NA for every correlation.
5. **WS5 does not reuse `scripts/run_subsample.R`.** Its stratum vocabulary does
   not express a region-share grid.
