# Fitness for use

What these estimates support, what they do not, and the measured numbers behind
each rule. Written for national partners deciding whether to act on a district
map.

Every threshold below is source-mapped. Where a number could not be produced by
running code on this branch it is written "not yet computed" rather than
estimated.

---

## 1. The short version

| Question | Answer | Basis |
|:---|:---|:---|
| Can I read a national prevalence from this? | **Yes**, where a survey exists. | Section 3 national recovery, 24 of 24 inside the survey 95 percent CI |
| Can I read an Admin-1 (regional) prevalence? | **Yes**, from the survey directly. | WS5: full-survey Admin-1 MAE 0.297 pp |
| Can I read an Admin-2 (district) prevalence? | **No.** | WS5: Admin-2 MAE 5.09 to 13.66 pp at full survey |
| Can I read an Admin-2 **ranking**? | **Only with the caveats in section 4.** | WS4b: median achieved r 0.200 |
| Will a bigger survey fix the district map? | **No.** | WS5: fourfold sample change moves Admin-2 MAE about 1 pp |
| Will better covariates fix it? | **Not established, and the evidence is against it.** | WS4b: 81 DHS indicators reach a median r of 0.071 |

---

## 2. Resolution: where to stop

**Report at Admin-1. Do not report district prevalence.**

The measured gap is large and it is not close:

| Level | Mean absolute error against the full survey | Source |
|:---|---:|:---|
| Admin-1, from the survey | **0.297 pp** | `results/tables/anchoring_design_curve.csv`, `mae_admin1_pp`, all regions and clusters retained |
| Admin-2, best arm | **7.38 pp** | `results/tables/anchor_controls.csv`, arm `2a flat REGIONAL mean (no covariates)`, `mae_pp` |
| Admin-2, anchored covariate model | 8.91 pp | same file, arm `3 ADMIN-1 anchor (hard)` |
| Admin-2, unanchored covariate model | 10.77 pp | same file, arm `1 no anchor (LORO)` |

A district estimate carries roughly **twenty-five times** the error of the
regional estimate it is derived from. No arm tested closes that gap.

**Measurements per unit.** The median district contributes 12 to 36 biomarker
measurements depending on country, and single-cluster districts dominate three
of the four surveys: Gambia 17 of 30, Ghana 62 of 75, Malawi 74 of 87 hold
exactly one survey cluster (measured from the `outcome_data_*` targets;
reproduced by `scripts/accuracy_impact/ws1e_poststrat.R`). A district with one
cluster has no within-district replication at all, so its estimate carries the
full between-cluster variance with nothing to average it against.

**The intermediate resolution question is now answered: there is no useful
intermediate.** WS4a swept Admin-1, Admin-1 cut into two and three spatially
compact parts, and Admin-2, holding the learner, folds and covariates constant.
Over the 14 cells present at all four levels, mean out-of-fold r is **0.315**
(Admin-1), **0.166** (two-way split), **0.093** (three-way split) and **0.086**
(Admin-2); Admin-1 is best in **9 of 14** cells and beats Admin-2 in **13 of 14**
(source: `results/tables/resolution_sweep.csv`, column `r_oof`). Skill falls
monotonically with resolution and no intermediate level is a peak.

Expressed as survey design: skill is highest at a median of **33 measurements per
unit** and has halved by **19.5**. Below roughly **20 measurements per unit**
these covariates carry very little, and Admin-2 sits at a median of **10 to 13**.

**A constraint that follows from the same sweep.** Two of the four countries
cannot support a within-country covariate model at Admin-1 at all: Sierra Leone
has 4 regions and Gambia 6, so leave-one-region-out would train on 3 and 5 units
respectively. Where the survey supports the estimate there are too few units to
learn a map; where there are enough units the estimate is too noisy to trust.

---

## 3. Reliability: which cells to suppress

**Suppress the Admin-2 layer for any cell whose empirical reliability is below
0.30.**

The threshold comes from the measured distribution rather than a convention.
Across 24 cells the empirical `r_max` (split-half, Spearman-Brown corrected,
`results/tables/reliability_empirical.csv`, `scheme == "within"`) is bimodal:
three cells sit at exactly **0.000**, and the next lowest are 0.183 and 0.250.
There is then a gap to 0.384. A cut at 0.30 separates the cells where the
survey's own two halves do not agree about which district is worse from the rest.

Cells failing that cut, and therefore suppressed at Admin-2:

| Country | Outcome | Empirical `r_max` |
|:---|:---|---:|
| Sierra Leone | women_b12 | 0.000 |
| Malawi | women_vitA | 0.000 |
| Malawi | women_b12 | 0.000 |
| Sierra Leone | child_vitA | 0.183 |
| Sierra Leone | child_iron | 0.250 |

**Do not use the analytic ceiling for this decision.** It reports exactly zero in
10 of 24 cells, seven of which the measurement contradicts, and it understates
the attainable correlation by a mean of 0.161 (`docs/findings/WS1_RELIABILITY_CEILING.md`).

---

## 4. How to display an Admin-2 map

**Rankings, with rank intervals, and never a prevalence number.**

1. **Show a rank, not a value.** Median achieved out-of-fold correlation at
   Admin-2 is **0.200** (`results/tables/reliability_skill_curve.csv`,
   `family == "micronutrient"`). A correlation of 0.2 supports a coarse ordering
   and nothing finer.
2. **Show a rank interval.** A point rank from a model at r 0.2 across 14 to 87
   districts is not distinguishable from many neighbouring ranks.
3. **Do not colour by WHO severity band at Admin-2.** The bands are 5, 20 and 40
   percentage points apart for iron (`config/who_thresholds.csv`), and the
   Admin-2 error is 7.4 to 13.7 pp. A single band boundary sits inside the error
   bar.
4. **Say which survey supplied the level.** A district map for a country with a
   survey is that survey's regional estimate applied to districts. A transported
   map for a country without one has no anchor and is a ranking only.

---

## 5. What a survey planner should buy

From WS5 (`docs/findings/WS5_ANCHORING_BUDGET.md`):

- **Buy clusters inside the regions you sample, not more regions.** Admin-1 MAE
  moves from **0.297 pp** at full cluster retention to **5.998 pp** at a quarter,
  a factor of twenty. Region coverage barely matters by comparison: Malawi and
  Sierra Leone are flat within noise across region shares from a quarter to all.
- **Do not buy a bigger survey to improve the district map.** Keeping a quarter
  of clusters instead of all costs a mean of **1.054 pp** in Admin-2 MAE. The
  district error is not budget-limited.

---

## 6. Uncertainty display rules

1. Every national and Admin-1 estimate carries a design-based interval.
2. Every Admin-2 quantity carries a rank interval, not a prevalence interval.
3. Exceedance probabilities against a WHO band are reported at national and
   Admin-1 level only, for the reason in section 4 point 3.
4. Where a threshold row in `config/who_thresholds.csv` is marked
   `source_verified = FALSE`, any band label derived from it is provisional.
   **All eight rows are currently so marked**, and three of them (folate, B12,
   zinc bands) are recorded as project conventions with no known WHO source.
   This should be resolved before any band label is published.

---

## 7. What would change these rules

- WS4a is now computed and found no intermediate optimum, so a finer reporting
  geography is not the route. What would change these rules is a covariate set
  that works at all, not a different aggregation of the current one.
- A covariate set that predicts health and nutrition outcomes rather than
  socioeconomic gradients. WS4b measured the current set at a median r of 0.071
  across 81 DHS indicators, with education and wealth the exceptions at 0.63 to
  0.81 and stunting, water and sanitation near or below zero.
- An anchoring scheme that beats an unanchored model out of sample. The one
  tested does not: under a jackknife it scores 0.147 against 0.156 unanchored
  (`results/tables/anchor_controls.csv`).
