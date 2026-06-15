# Corrected pipeline — results note

Status of the parallel P1–P8 pipeline (`_targets_corrected.R`). Cross-reference:
[`docs/CRITICAL_REVIEW.md`](CRITICAL_REVIEW.md), methods write-up
[`docs/METHODS_CORRECTED.md`](METHODS_CORRECTED.md), run instructions
[`docs/CORRECTED_PIPELINE_RUN.md`](CORRECTED_PIPELINE_RUN.md).

---

## FULL RUN RESULTS (completed 2026-06-14 23:12)

`CORRECTED_SCOPE=full`, sequential, detached driver; `tar_make` built **289 of
293 targets (4 cached), 0 errored** in ~14 min, then automatically chained HIV
Stage 3 (handoff confirmed: HOLD sentinel deleted, HIV pid 48196 running). The
6 Sierra Leone slices returned no corrected fit (production reference only); the
numbers below are over the **18 slices with corrected metrics**.

| P | Corrected (honest) | Optimistic / production | Headline |
|---|--------------------|-------------------------|----------|
| **P1** CV honesty | spatial-block AUC **0.518**, cluster **0.547** | optimistic **0.571** | leakage inflates AUC by **+0.053** (opt−spatial), 15/18 slices; spatial most conservative |
| **P2** calibration | out-of-fold Platt slope mean **−0.17** | in-sample slope **1.00** (all 18) | in-sample "perfect" calibration is an artefact |
| **P3** decision value | mean lift **1.29×**, reach@20% **0.26**, hit **0.31** | no-targeting = 1.0 | model beats status quo in **13/18** slices |
| **P4** intervals | split-conformal coverage **0.92** (target 0.90); **1,266** district design CIs | Admin-1 broadcast (no guarantee) | coverage restored; no broadcast |
| **P5** error | sampling-adjusted RMSE **15.4 pp** | naive **15.6 pp** | removes 2.5 pp survey sampling noise |
| **P6** trust | **2,080** OK / **1,294** caution / **278** do-not-rely | (none) | unreliable areas flagged |
| **P7** area model | FH (REML) **486** rows + EB fallback **780** | comparator-only | design-aware estimator promoted to primary |
| **P8** integrity | 24 slices, **0 errored**, no cartesian blow-up | (blow-up risk) | by construction |

**Headline:** across the 18 corrected slices the honest spatial-block CV is on
average **0.05 AUC below** the optimistic regime, and in-sample Platt calibration
(slope 1.00) collapses to a near-zero/negative honest out-of-fold slope —
quantitatively confirming the two flagship problems (P1, P2) from the critical
review. Decision-value, design-based interval coverage, sampling-aware error,
trust flags and the partial-pooling area estimator all produced first-class
outputs now rendered in the dashboard "Methods (corrected)" tab.

**Known gap:** the 6 Sierra Leone slices (`sierraleone_*`) returned `NULL` from
`fit_corrected_sl` (no corrected metrics; production reference retained) — to be
debugged and backfilled by resuming the cached run (`tar_make` resume) once the
machine is free; likely a survey-design-column or edge-case in the SL
`outcome_data`. This does not affect the other 18 slices or the HIV handoff.

---

## SMOKE-SLICE VALIDATION (retained for reference)

**Smoke slice:** Gambia × women_iron (n = 1402 individuals, 70 clusters,
6 regions, 30 surveyed districts). Library mean+glmnet+ranger, predictor pool
capped at 400 (unsupervised), V = 5. Full `tar_make` of all 14 targets: **7.3 s**,
sequential. Identity-controlled comparisons (same data + same library) isolate
the *method* effect, not a library difference.

## Headline: honest validation is lower than optimistic (P1)
Out-of-fold ROC-AUC / Brier on identical data and library:

| Scheme | Honest? | AUC | Brier |
|--------|---------|-----|-------|
| random CV + full-data preprocessing (optimistic) | no | **0.592** | 0.209 |
| cluster-block CV (in-fold preprocessing) | yes | 0.579 | 0.212 |
| region/spatial-block CV (in-fold preprocessing) | yes | **0.579** | 0.214 |
| production `cv_perf` (leaky preprocessing, random CV) | no | 0.579 | 0.214 |

The optimistic regime inflates AUC by **+0.013** and Brier-skill ~2–6× relative
to the deployment-honest schemes on this slice (an earlier identical-library run
showed Brier-skill 0.026 optimistic vs 0.004 spatial-block). The gap is the
optimism removed by **P1** (moving imputation, the **supervised prescreen**,
correlation filter and normalization inside each fold) plus spatial-block CV.
Spatial-block is the most conservative — exactly the literature's prescription
(Roberts 2017; Ploton 2020; see [`docs/LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md)).

## Calibration: in-sample is an illusion (P2)
| Calibration | Brier | Intercept | Slope |
|-------------|-------|-----------|-------|
| raw (uncalibrated) | 0.212 | −0.55 | 0.26 |
| in-sample Platt (optimistic) | 0.213 | 0.00 | **1.00** |
| out-of-fold Platt (honest) | 0.216 | −0.83 | **−0.05** |

In-sample Platt reports a *perfect* slope of 1.00 — a fitting artefact. Honest
out-of-fold Platt reveals a slope ≈ **0**: once you stop scoring on the data you
calibrated on, this model has essentially **no calibratable signal** for this
slice. This is the P2 finding made concrete.

## Other corrected targets (smoke values)
- **P5 sampling-aware error:** naive RMSE 22.4 pp vs sampling-adjusted 21.9 pp
  (survey sampling RMSE 4.6 pp removed); Pearson r 0.63. Apparent error is
  partly the survey's own noise.
- **P3 decision value:** reach@worst-20% = 0.23, lift vs no-targeting = 1.15,
  top-quintile hit-rate = 0.17 — modest within-country targeting value, now a
  first-class pipeline target.
- **P4 intervals:** split-conformal (held-out calibration split) half-width
  37.3 pp, empirical held-out coverage 0.87 (target 0.90); 30 district-level
  design-based survey CIs — **no Admin-1 → Admin-2 broadcast**.
- **P6 trust flags:** 32 "OK to use", 3 "Use with caution", 2 "Do not rely"
  (out-of-support districts) across the 37 polygons.
- **P7 area model:** design-aware partial pooling produced per district
  (direct vs partial-pooled vs ML side by side).
- **P8:** enforced by construction in `R_corrected/` — every survey/area join is
  keyed on `(country, Admin2)` (no cartesian blow-up), train/test
  standardization uses training statistics only, and out-of-sample imputation
  uses **training** medians (never the target rows' own values).

## Per-P implementation status
| P | Status | Notes |
|---|--------|-------|
| P1 | ✅ implemented & validated | in-fold preprocessing + cluster & spatial-block CV; honest<optimistic shown. |
| P2 | ✅ implemented & validated | fold-aware Platt; in-sample-vs-OOF contrast shown. |
| P3 | ✅ implemented & validated | targeting/lift as first-class targets. |
| P4 | ✅ implemented & validated | split-conformal (held-out split) + district design CIs; no broadcast. |
| P5 | ✅ implemented & validated | sampling-variance-adjusted RMSE. |
| P6 | ✅ implemented & validated | out-of-support distance + sample-size + width → traffic light. |
| P7 | ⚠️ implemented; FH fell back | Fay-Herriot was **computationally singular** on 30 collinear-covariate areas → empirical-Bayes shrinkage fallback (still design-aware). FH needs better-conditioned covariates (fewer/decorrelated) or BYM2; queued for the full run. |
| P8 | ✅ enforced by construction | joins on (country, Admin2); training-only stats; no OOS-median leakage. |

## Known limitations of the smoke (not bugs)
- **One slice, weak-signal outcome.** Gambia × women_iron has little proxy
  signal, so absolute AUCs hover near 0.5–0.6 and the leakage gap, while
  correctly signed, is small. Larger gaps are expected for higher-signal
  outcomes and bigger predictor pools (the optimistic prescreen overfits more
  when choosing from more candidates) — the full run will quantify this across
  all slices.
- **Discrete SuperLearner + capped pool** are lightweight stand-ins for the full
  stacked SL and complete proxy set; both are configurable for the full run
  (`CORRECTED_SCOPE=full`). The inner learner-selection reuses the fold-fitted
  preprocessing (documented simplification); the **outer** OOF honesty — which
  defines the reported metric — is exact.
- **Decision-value weighting** uses `n_svy` as a population proxy inside the
  pipeline (avoids an external dependency); the dashboard's Decision-value tab
  uses true subgroup population.

## What is staged for the full run
Set `CORRECTED_SCOPE=full` and run unattended once the HIV machine frees (see the
run note). That builds all available country × outcome slices with the fuller
library/pool, refreshes `results/tables/corrected/*.csv` and the dashboard
bundle, and lets the "Methods (corrected)" tab show every slice's honest-vs-
optimistic deltas head-to-head.
