# IL-02. Person-level models: the January figure, IL-01, and the index at person level

27 September 2026. Asked by the PI: why are the person-level models so null when
the January 2026 Ghana deck showed much better accuracy; was January leakage, or
is the current approach too conservative; and what are the person-level AUC and
Brier score of the PCA domain index?

## Answers

1. **The January figure reported in-sample fit, not prediction.** Its AUC and
   Brier skill equal the full-data fit scored on its own training rows, to three
   decimals, in 18 of 18 models. The same fits' own cross-validation gave AUC
   0.54-0.79 and Brier skill 0.00-0.22. The figure cannot be presented.
2. **The rest of the gap to IL-01 is leakage, not conservatism.** On the
   January-era data with one fixed learner, moving the supervised prescreen
   inside the folds costs 0.043 in AUC, and removing identifiers, sampling-design
   columns and blood-derived columns costs 0.061 more. District-blocked folds
   (the "conservative" choice) then cost 0.005.
3. **IL-01 has defects in both directions.** Haemoglobin sits in Malawi's
   "survey" set: women's iron AUC falls from 0.746 to 0.596 without it. Ghana's set
   carries the malaria RDT result and thalassaemia genotype. In the other
   direction, IL-01's scoring makes weak cells look worse than null: a constant
   predictor scores a pooled AUC of 0.40-0.43 under its folds, not 0.50, and 44
   percent of its published rows sit below 0.5.
4. **The index at person level:** mean AUC 0.56 (0.61-0.63 on within-fold
   pairs), Brier skill +0.013 over the national prevalence. No district-level
   predictor can do much better. Knowing the prevalence among each district's
   other respondents gives AUC 0.585 and Brier skill 0.028. Person-level AUC is
   the wrong yardstick for an area-level index.
5. **The January figure, redone on current data (Ghana, section 7).** District
   proxies give person-level Brier skill of 1-9 percent and MSE skill of 4-12
   percent on log concentrations. That is near the between-district ceiling for
   iron and about half of it elsewhere. The questionnaire adds signal only for
   child iron (age) and women's RBP (BMI, pregnancy).
6. **Separate finding: vitamin A.** `targets_v2.csv` was built on 8 September,
   before the 15 September switch to `VITA_RULE=rbp070`, and was not rebuilt.
   Every protocol-v2 vitamin A result reads it, including the dashboard's
   deployment index, so all of them are on the retinol-equivalent rule.

## 1. What the January figure measured

The figure ("Model results by outcome and predictor set", Brier skill 3-67
percent) plots `mn-proxies/results/res_full_bin_GW_Ghana_SL_all.rds`. That file
was built on 22 August 2025 by `src/Ghana/model_performance_bin.R` from eighteen
sl3 fits. The scoring function `auc_pr_brier_from_sl3(pred_source = "auto")` took
`res$sl_fit$predict()` first. Its comment calls these "CV preds", but on an sl3
`Lrnr_sl` a bare `$predict()` is the full-data fit predicting its own training
rows. By February 2026 that line was commented out, but the figure file was never
regenerated.

| Outcome | Set | Figure AUC | Figure Brier skill | In-sample AUC | In-sample Brier skill | CV AUC | CV Brier skill |
|:--|:--|--:|--:|--:|--:|--:|--:|
| Child vitamin A | Survey | 0.966 | 0.383 | 0.966 | 0.383 | 0.654 | 0.064 |
| Women vitamin A | Survey | 0.978 | 0.285 | 0.978 | 0.285 | 0.738 | 0.036 |
| Women B12 | Survey | 0.978 | 0.384 | 0.978 | 0.384 | 0.757 | 0.087 |
| Women folate | Survey | 1.000 | 0.671 | 1.000 | 0.671 | 0.658 | 0.090 |
| Child iron | Survey | 0.885 | 0.379 | 0.885 | 0.379 | 0.787 | 0.217 |
| Women iron | Survey | 0.865 | 0.216 | 0.865 | 0.216 | 0.690 | 0.065 |
| Child iron | All | 0.978 | 0.599 | 0.978 | 0.599 | 0.779 | 0.209 |
| Women B12 | Proxy | 0.886 | 0.209 | 0.886 | 0.209 | 0.701 | 0.100 |

All eighteen rows are in `il02_january_figure_vs_cv.csv`. The fits
themselves used 10-fold cluster-blocked CV: none of the 90 clusters is split
across folds. The fault is in what was scored, not in how the models were fitted.
The B12 slide (AUC 0.689, Brier skill 0.095) came from a different,
cross-validated multi-country set (`compiled_predictions.RDS`, 14 January). In
that set most cells were already near null.

The January fits' own CV was optimistic too. Of the 478 covariates the
six survey fits received, 148 are cluster numbers, household and person ids,
cluster-constant sampling-design columns (PSU weight, listing counts, selection
probabilities), region, team or phlebotomist codes, dates, or blood-derived
columns: the malaria RDT result and referral, sickle and alpha-thalassaemia
genotype, and MRDR selection (`il02_january_covariates.csv`). The supervised
prescreen (`washb_prescreen`, p < 0.2) also ran on all rows before the folds were
drawn.

## 2. From January's cross-validation to IL-01, one factor at a time

The rerun used the January-era Ghana data (the mn-proxies copy of 13 January,
which reproduces n and prevalence for five of the six fits). It kept one learner
stack throughout: mean, lasso, ridge and ranger, with NNLS weights as in
January. Each variant ran over five fold draws. Values are pooled AUC and Brier
skill against the national prevalence, averaged over child vitamin A, child
iron, women's iron, folate and B12.

| Step | AUC | Brier skill |
|:--|--:|--:|
| January figure (in-sample) | 0.939 | 0.407 |
| S1 January survey set, cluster folds, prescreen on all rows | 0.713 | 0.099 |
| S2 prescreen inside the training folds | 0.670 | 0.073 |
| S3 ids, design, region, team, date and blood columns removed | 0.609 | 0.039 |
| S4 same, district-blocked 10-fold | 0.600 | 0.036 |
| S5 same, district-blocked 5-fold (IL-01's folds) | 0.604 | 0.035 |
| S7 S5 without any prescreen | 0.607 | 0.035 |
| IL-01 as published, Ghana survey-only, same five outcomes | 0.607 | 0.031 |
| P1 January proxy set, cluster folds, prescreen on all rows | 0.631 | 0.054 |
| P4 proxies, region and month removed, district 5-fold | 0.615 | 0.040 |

S1 reproduces the saved fits' own CV cell by cell (child iron 0.784 / 0.210
against 0.787 / 0.217). Three steps account for nearly all of the gap:

- **In-sample scoring:** AUC 0.939 to 0.713.
- **Prescreening on all rows:** 0.713 to 0.670. The loss is largest where n is
  smallest (B12 0.771 to 0.684, folate 0.664 to 0.605, n about 475).
- **Identifier, design and blood columns:** 0.670 to 0.609. Kept under district
  folds (S4b), they still lift AUC to 0.666, because cluster numbers and
  sampling-design columns encode geography.

District blocking costs 0.005-0.009 once those columns are gone. Child iron is the
one strong survey cell (0.72 under district folds, carried by the child's age).
Per-outcome values are in `il02_decomposition_cells.csv`.

## 3. IL-01 audit (`scripts/protocol_v2/46_individual_level_models.R`)

**Leaks in the published run** (`il02_il01_survey_columns.csv`):

- **Malawi:** `m432` is women's haemoglobin (range 5.8-17.5; r = -0.73 with the
  `anemia` flag). `m228` is child haemoglobin (`m228 < 11` reproduces the flag
  in 99.9 percent of children). Rerun through IL-01's own code, same folds,
  without them: women's iron falls from AUC 0.746 / Brier skill +0.098 to 0.596 /
  -0.005. Child iron holds (0.784 to 0.781), carried by age (`m07`), height
  (`m232`) and weight (`m231`).
- **Ghana:** `gw_gcmst` (child malaria RDT result; univariate AUC 0.61 for child
  iron), alpha-thalassaemia genotype (`gw_cAlphaThal*`), and a phlebotomist id
  (`gw_pcn`).
- **Sierra Leone:** phlebotomist numbers and urine-specimen flags (weak).
- **Ghana and Gambia:** cluster numbers (`gw_bccn`, `gw_cn`, MICS `HH1`),
  `gw_childid`, region, and cluster-constant sampling-design columns, which are
  the top univariate predictors in most Ghana cells. Under district folds these
  are geography rather than leakage, but they are not questionnaire content.

**The regex is fragile.** It matches `ferr` but not `gw_fer` (ferritin), and
does not list body iron stores (`gw_bis`). `sTfR` is case-sensitive, so it
misses `gw_sTFR`. The anchored patterns (`^sf_`, `^fol`, `^rbp`, `^vit`) can
never match a `gw_`-prefixed name. Applied to the January-era columns, IL-01's
rule lets ferritin through and scores iron AUC 0.97-0.99 (S6). The store data
IL-01 actually read lacks those columns, so the published run escaped this. The
project's tested guard (`is_biomarker_column()` / `allowed_under_arm()` in
`R/data_prep.R`) is the safer rule.

**Also too strict.** `Year` removes women's age (`gw_wAgeYears`), and the
80-column coverage cap then drops `gw_age_days`, BMI, education and sex in
Ghana's women's cells.

**Scoring pushes weak cells below null.** Out-of-fold predictions are pooled
across folds before the AUC is taken, so fold-to-fold differences in the training
mean enter the ranking. A constant predictor scores 0.40-0.43 under these folds
(`il02_decomposition_raw.csv`, learner `mean`). The Brier skill reference is the
full-sample prevalence, which no cross-validated model can match without signal
(-0.005 to -0.007 for a constant). A single fold draw is used. The fix is
within-fold AUC, the training-mean null and replicated draws. The protocol-v2
notes make the same point about the calibrated index.

## 4. The PCA domain index at person level

The index is fitted exactly as in `02_run_benchmarks_v2.R` (`build_cell`, domain
PCs to 80 percent of variance, `make_folds_v2("kfold_district")`, 20 draws). Each
held-out district's out-of-fold prediction is given to its respondents, whose own
deficiency flags are then scored. District-level Spearman reproduces the
protocol's in-fill benchmark cell by cell (Ghana child iron 0.52 against 0.50,
Gambia women's iron 0.75 against 0.74). Sierra Leone is scored
leave-one-district-out, because its 5-fold draws leave 11 training districts,
below the protocol's minimum of 12. Vitamin A is scored under
`VITA_RULE=retinol_equiv` so that respondents and district targets share one
definition (see section 5).

Mean over 19 cells:

| Method | AUC (pooled) | AUC (within-fold) | Brier skill vs national | vs training mean |
|:--|--:|--:|--:|--:|
| Index, calibrated levels (the protocol's product) | 0.529 | 0.608 | +0.009 | +0.014 |
| Index, logistic recalibration on respondents | 0.564 | | +0.013 | +0.018 |
| Index trained on the level target, logistic | 0.564 | 0.626 | +0.015 | +0.020 |
| Ceiling: district's other respondents, shrunk | 0.585 | | +0.028 | |
| Ceiling: cluster's other respondents, shrunk | 0.595 | | +0.027 | |

The within-fold AUC is left blank for the two rows where it was not computed. The
in-sample district prevalence (self included) scores 0.77 / 0.11, but that
ceiling is inflated by each respondent counting towards their own district's
rate. The leave-self-out rows are the honest ceiling.

Best cells, index with logistic recalibration (ceiling in brackets):

| Cell | AUC | Brier skill |
|:--|--:|--:|
| Ghana child iron | 0.675 (0.695) | +0.077 (0.105) |
| Malawi women B12 | 0.676 (0.667) | +0.025 (0.045) |
| Ghana women B12 | 0.654 (0.636) | +0.046 (0.009) |
| Malawi folate | 0.637 (0.729) | +0.037 (0.103) |
| Gambia women's iron | 0.609 (0.573) | +0.037 (0.032) |

The weakest cells are Malawi iron (0.54-0.55), rare vitamin A (women, and
Malawi children at 0.4 percent under the retinol-equivalent rule) and all of
Sierra Leone. Per-cell values are in `il02_index_person_summary.csv`.

What this means: most variation in individual deficiency lies within districts
(and within clusters), so any district-level predictor, including a perfect one,
discriminates individuals only weakly. On person-level AUC the index reaches
about what knowing the district's other respondents would give, and it recovers
about half of that ceiling's Brier skill. Judge it on district-level metrics.
"A coin toss" is the wrong summary: the index is near the attainable ceiling, and
that ceiling is low.

## 5. Separate finding: vitamin A district targets on the retinol-equivalent rule

`results/tables/protocol_v2/targets_v2.csv` was last built on 8 September
(commit 041ec7e). On 15 September the vitamin A default in
`R/brinda_adjustment.R` became `rbp070`, with `retinol_equiv` as the sensitivity
analysis (`docs/survey_report_reconciliation.md` sections 4.2 and 6, commit
672bcb2). The targets were not rebuilt. The benchmarks (`02`, `02b`), external
validation (`scripts/external_validation/04`), the survey planner (`64`),
cross-outcome borrowing (`63`) and the dashboard bundles (`05`, `00_read_targets`)
all read the file.

Evidence: its national vitamin A levels are the reconciliation's
retinol-equivalent figures, not the RBP ones:

| Children | File | Retinol-equivalent (reconciliation) | RBP < 0.70 |
|:--|--:|--:|--:|
| Ghana | 28.0% | 27% | 15% |
| Sierra Leone | 6.7% | 6.5% | 12% |
| Malawi | 0.4% | 0.4% | 10% |

District prevalences recomputed from respondents under `rbp070` correlate with
the file at 0.72-0.83 in Ghana and Sierra Leone and 0.25-0.28 in Malawi, against
1.00 for iron, folate and B12. Recomputed under `retinol_equiv`, every cell
correlates at 1.00. Gambia is unaffected (1.00 under both rules). Iron, folate
and B12 results are not affected. Rebuilding the targets under `rbp070` means
running `01_build_targets_v2.R`, then `02` / `02b`, then the dashboard bundles.
That changes published vitamin A numbers, so it is the PI's decision.
Status, evening of 27 September: RR-14 (`scripts/protocol_v2/rr14_outcome_datasets.R`,
run from another session) rebuilt the four countries' outcome datasets at 16:30
so that the 15 September fixes reach script 01, and `01`/`02` were re-running.
Re-check the vitamin A cells, and section 4, once it lands.

## 7. The January figure, updated to the current data and models

Section 4 scored the index on `targets_v2.csv` as it stood on 27 September
(vitamin A under the retinol-equivalent rule). This section redoes the January
figure itself for Ghana on the current data. It uses the targets store as rebuilt
by RR-14 at 16:30 on 27 September, which moved a handful of women through the
corrected population masks. Outcome definitions are the pipeline's own, with
vitamin A on `rbp070`.

The models are IL-01's person-level SuperLearner (`fit_area_superlearner`: mean,
elastic net and ranger, survey weights, district-blocked inner folds, discrete
pick). Four predictor sets feed it:

- **Survey only:** the respondent's cleaned questionnaire columns, 48 for
  children and 80 for women. The filter is `allowed_under_arm("questionnaire")`,
  plus IL-01's regex with ages exempted from its date patterns, minus the
  identifier, sampling-design, fieldwork-timing and biomarker-module columns
  (RDT, genotypes, phlebotomy).
- **Proxies:** IL-01's 106 district domain components from the 533-column proxy
  table.
- **Survey + proxies:** both sets together.
- **Index:** the PCA domain index, calibrated to respondents (logistic for flags,
  linear for concentrations).

Scoring is outer 5-fold by district, 5 draws averaged, with 95 percent
cluster-bootstrap intervals (`il02_honest_person_level.csv`; figures
`results/figures/il02_person_level_{brier,mse}_{honest,proxy_only}.png`).

| Outcome (prevalence) | Survey | Proxies | Survey + proxies | Index | Ceiling |
|:--|--:|--:|--:|--:|--:|
| *Brier skill, percent* | | | | | |
| Child vitamin A (14.9%) | 0.9 | 1.3 | 1.4 | 1.9 | 0.1 |
| Women vitamin A (1.5%, 15 cases) | 0.1 | -0.5 | -0.1 | -1.8 | 1.5 |
| Child iron (23.9%) | 9.8 | 9.0 | 17.0 | 8.8 | 8.6 |
| Women iron (13.3%) | 0.3 | 2.6 | 1.8 | 1.2 | 1.9 |
| Folate (54.8%) | -0.2 | 3.9 | 5.1 | 1.0 | 13.7 |
| B12 (8.5%) | 0.0 | 6.0 | 4.7 | 5.4 | 13.8 |
| *MSE skill, percent, log concentration* | | | | | |
| Child RBP | 2.9 | 4.7 | 4.8 | 3.6 | 9.1 |
| Women RBP | 13.6 | 6.6 | 15.4 | 5.1 | 13.3 |
| Child ferritin | 15.1 | 11.8 | 25.2 | 11.2 | 16.2 |
| Women ferritin | 2.9 | 6.1 | 7.4 | 3.7 | 8.7 |
| Folate | -0.4 | 8.0 | 9.4 | 2.9 | 14.9 |
| B12 | 0.4 | 10.1 | 7.8 | 8.0 | 19.5 |

**The ceiling.** For a yes/no outcome at prevalence p, the person-level variance
is p(1 - p). Everyone in a district receives the same prediction from a
district-level model, so the best it can do is give each person the district's
true rate. That removes only the variance of the true district rates (tau
squared), so its Brier skill is tau^2 / [p(1 - p)], and for concentrations the
corresponding share of variance. The rest is which people within a district
are deficient, which no district-level information can resolve.

A worked example: true district rates spread with a standard deviation of 7
points around 24 percent (most districts between 10 and 38 percent, a large and
policy-relevant spread) give tau^2 = 0.0049 against p(1 - p) = 0.182, a maximum
Brier skill of 2.7 percent. In Ghana the estimated share (method of moments,
district bootstrap) is about 0 for child vitamin A, 2 percent for women's iron
and vitamin A, 9 percent for child iron and about 14 percent for folate and B12
(wide intervals with about six women per district). For log concentrations it
is 9-20 percent.

**Reading the table.**

- **The district-level models land near their ceiling for iron** (child iron
  8.8-9.0 against 8.6). For folate, B12 and the concentrations they reach half
  or less of it. The SuperLearner on the domain components and the index overlap
  in every cell, with the SuperLearner slightly ahead in most.
- **The questionnaire carries real, physiological signal in two places.** Child
  iron and ferritin (the child's age) and women's RBP (BMI and pregnancy; RBP4
  rises with adiposity). Elsewhere it adds nothing.
- **Where the two carry different information, they add.** Child ferritin reaches
  25 percent with both, because age explains variation within districts and
  proxies explain it between them.

## 8. Places that carry the old account (not edited here)

- `docs/slides/MNF15-talk-2026-09.qmd:580`: the speaker note attributes the January
  figure to folds and biomarker-adjacent columns. It was in-sample scoring.
- `docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd:851` and the full talk:
  the same account, plus "survey-variable AUC up to 0.78 in Ghana and Malawi".
  Malawi women's iron is haemoglobin. Ghana and Malawi child iron stand.
- `docs/findings/SANDBOX_LOG_2026-09.md:1022`: same account.
- `dashboard/R/mod_technical.R:111,220`: "an AUC of about 0.52, close to a coin
  toss" averages all three IL-01 sets on the pooled-AUC scale, where no skill
  reads 0.40-0.43.
- `46_individual_level_models.R`: use the `allowed_under_arm()` guard, drop
  `m228`/`m432`, and score with within-fold AUC, the training-mean null and
  replicated draws.

## Reproduction

| Script | Output (`results/tables/protocol_v2/`) |
|:--|:--|
| `scripts/protocol_v2/67_january_individual_forensics.R` | `il02_january_figure_vs_cv.csv`, `il02_january_covariates.csv` |
| `scripts/protocol_v2/68_individual_decomposition.R` (+ `_lib.R`) | `il02_decomposition_raw.csv`, `_cells.csv`, `_avg.csv` |
| `VITA_RULE=retinol_equiv scripts/protocol_v2/69_index_person_level.R` (+ `_lib.R`) | `il02_index_person_raw.csv`, `il02_index_person_summary.csv` |
| `scripts/protocol_v2/70_il01_column_audit.R` | `il02_il01_survey_columns.csv` |
| `IL_COUNTRY=Malawi IL_OUTS=child_iron,women_iron IL_DROP=m228,m432 scripts/protocol_v2/46_individual_level_models.R` | `individual_level_models_Malawi_child_iron-women_iron_drop-m228-m432.csv` |
| `NW=8 REPS=5 BOOT=1000 IL_CACHE=<local dir> scripts/protocol_v2/71_person_level_honest.R` (+ `_lib.R`) | `il02_honest_person_level.csv`, `il02_honest_survey_columns.csv` |
| `scripts/protocol_v2/72_plot_person_level_honest.R` | `results/figures/il02_person_level_{brier,mse}_{honest,proxy_only}.png` |

Script 71 reads the targets store (`_targets_full`) once at start-up and hands
the workers compact per-outcome matrices. Ten workers each loading the store's
outcome datasets ran the machine out of memory. Its respondent-level
predictions go to `IL_CACHE` and must not be committed. Scripts 67, 68 and the
68 lib read the mn-proxies repository (the saved January fits and the 13 January
data); set `MN_PROXIES` to relocate. The committed 68
tables came from a scratch copy of the same code that differed only in output
paths (61 minutes on 17 workers, 5 draws). Everything else ran from the committed
scripts. No respondent-level rows are written: every table holds aggregate
metrics or column names.
