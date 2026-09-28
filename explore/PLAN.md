# `explore/` — hypothesis-generating probes: implementation plan

> **Adaptation note.** The repo has no test suite by design (CLAUDE.md:
> "correctness is checked by running targets and inspecting `results/tables/`").
> So each task below is verified by *running the script and inspecting the
> output table*, not by a failing test. Steps use checkbox syntax for tracking.

**Goal.** Find proxy predictors and modelling approaches that do well for
specific micronutrients, in specific countries, and broadly — as *leads for
future work*, not production claims.

**Architecture.** A thin scorer (`explore/R/harness.R`) reproduces the
protocol-v2 estimands (replicated district in-fill within country;
leave-one-country-out Spearman) but accepts an arbitrary arm function and an
arbitrary predictor matrix. Every probe is a standalone numbered script in
`explore/scripts/` that builds features, screens loosely, and reports the
honest number from the harness next to the loose one. Results land in
`explore/out/`; leads are logged in `explore/FINDINGS.md`.

**Tech stack.** R 4.4.2. Already installed and used here: `glmnet`, `pls`,
`ncvreg`, `kernlab`, `grf`, `hal9001`, `dbarts`, `sva`, `limma`, `corpcor`,
`glasso`, `huge`, `fields`, `terra`, `sf`, `mgcv`, `data.table`. Earth Engine
via the reticulate venv at
`C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe`
(`ee` 0.1.370, initialises against account amertens@berkeley.edu) — Phase 2 only.

**Isolation contract.** `explore/` creates no file outside itself. It *reads*
`data/`, `results/tables/protocol_v2/targets_v2.csv`,
`data/covariates/harmonized/`, `dashboard/data/admin2_boundaries.rds` and
sources `R/protocol_v2.R` read-only. It never writes to `R/`, `scripts/`,
`results/`, `dashboard/`, `docs/`, `_targets.R`, or `metadata/`.

---

## File structure

| File | Responsibility |
|---|---|
| `explore/R/harness.R` | Load cells; `explore_cell()`, `score_infill()`, `score_loco()`; one arm signature for everything |
| `explore/R/features_mechanistic.R` | Nutrient-specific derived features from MapSPAM raw crops, GLW4, ESPEN, soil |
| `explore/R/features_temporal.R` | Cluster phenology phase, lag stacks, Fourier harmonics |
| `explore/R/methods_kernel.R` | Multi-kernel REML-BLUP, kernel builders, multi-trait low-rank BLUP |
| `explore/scripts/NN_*.R` | One probe each; writes `explore/out/NN_*.csv` |
| `explore/FINDINGS.md` | Running log: question, design, result, honest number, verdict |
| `explore/README.md` | What this is, how to run, what it is not |

---

## Task 0: Scaffolding and isolation contract

**Files:** Create `explore/README.md`, `explore/FINDINGS.md`, `explore/.gitignore`

- [ ] **Step 1:** Write `README.md` stating purpose, the isolation contract, how to run a probe, and that nothing here is a production claim.
- [ ] **Step 2:** Write `FINDINGS.md` with the header and an empty entries section.
- [ ] **Step 3:** Commit.

## Task 1: The shared harness

**Files:** Create `explore/R/harness.R`

Mirrors `scripts/protocol_v2/02_run_benchmarks_v2.R` `build_cell()` /
`run_one_draw()` but generalised: an *arm* is any
`function(tr, te, y, X, aux) -> numeric(length(te))`, and `X` is whatever the
probe supplies (raw store, mechanistic block, kernel, embedding).

- [ ] **Step 1:** `exp_load()` — read `targets_v2.csv`, the 578-column shared
      store + metadata, boundaries centroids; apply `drop_near_outcome_v2()`;
      return a list with `TG`, `S`, `MD`, `PREDS`, `domain_of`, `CENT`.
- [ ] **Step 2:** `exp_cell(cn, on, target, cols = NULL, prep = TRUE)` — the
      protocol-v2 cell: pair-key join, `prep_predictors_v2()` (rank-normalise
      within country + median impute), domain representation, `y_nat`/`y_mod`,
      weights `n_eff`, `aux` with lon/lat/Admin1. `cols` selects a feature
      subset; `prep = FALSE` returns raw columns for probes that do their own
      preparation.
- [ ] **Step 3:** `exp_infill(cell, arm, reps = 10, k = 5)` — replicated
      5-fold over districts, out-of-fold prediction, `score_v2()` on the
      natural scale, precision-weighted; returns one row per rep.
- [ ] **Step 4:** `exp_loco(cells, arm, target)` — pool cells across countries
      on the common columns, within-country standardise the outcome, fold by
      country, score Spearman *within each held-out country*.
- [ ] **Step 5:** `exp_baseline_arms()` — `null_train_mean`, `domain_index`
      (the estimator of record), `spatial` so every probe reports its
      comparator on the same folds.
- [ ] **Step 6: Verify** — `Rscript explore/scripts/00_harness_check.R`
      reproduces the domain index's LOCO level Spearman to within Monte-Carlo
      error of the 0.37 on the record, and its in-fill number against
      `results/tables/protocol_v2/benchmarks_v2_cells.csv`. **This is the gate:
      no probe result is trustworthy until the harness reproduces the record.**
- [ ] **Step 7:** Commit.

## Task 2 (Family A): Multi-kernel REML-BLUP

**Files:** Create `explore/R/methods_kernel.R`, `explore/scripts/01_kernel_blup.R`

Why: n=14–87, p=578 is the genomics regime, where the field's answer is not
selection but a relationship matrix plus REML-estimated shrinkage — one
*estimated* hyperparameter, not a tuned one.

- [ ] **Step 1:** `k_linear(X)` = `tcrossprod(X)/ncol(X)`, plus `k_gaussian()`,
      `k_spatial()` (Matérn on centroids). Centre and scale each kernel to
      mean diagonal 1 so variance components are comparable.
- [ ] **Step 2:** `reml_blup(y, Klist, w)` — EM/AI-REML for
      `y = sum_k u_k + e`, `u_k ~ N(0, K_k s2_k)`; returns variance components
      and BLUP predictions for held-out rows via the standard
      `K_te,tr V^-1 (y - mu)` formula.
- [ ] **Step 3:** Arm wrappers: `arm_blup_all` (one kernel, all predictors),
      `arm_blup_multi` (one kernel per domain), `arm_blup_cs` (climate+soil
      only), `arm_blup_spatial_plus_cs`.
- [ ] **Step 4:** Run over all 24 cells × both targets; write
      `explore/out/01_kernel_blup_{infill,loco}.csv` plus a variance-components
      table `01_kernel_varcomp.csv` (the principled replacement for ablation).
- [ ] **Step 5: Verify** — every cell has a finite score, the variance
      components are non-negative and sum sensibly, and the baseline arms in
      the same table match Task 1's gate values.
- [ ] **Step 6:** Commit.

## Task 3 (Family A): AlphaEarth as a kernel, not as columns

**Files:** Create `explore/scripts/02_alphaearth_kernel.R`

Why: `aef_A00..A63` is a learned embedding whose inner product *is* a
similarity. Scoring it as "a domain" collapses 64 dims through one PC
(+0.012/−0.012 on the record — i.e. nothing), which is the one representation
guaranteed to destroy it.

- [ ] **Step 1:** Extract the 64 `aef_*` columns; build a linear kernel and a
      cosine kernel on them.
- [ ] **Step 2:** Score four arms: embedding-kernel BLUP; embedding-PC1 (the
      current representation, as the comparator); embedding + climate/soil
      kernels jointly; climate/soil alone.
- [ ] **Step 3:** Run in-fill and LOCO; write `explore/out/02_alphaearth.csv`.
- [ ] **Step 4: Verify** — the embedding-PC1 arm reproduces the ~0 effect on
      the record, which confirms the comparison is like-for-like.
- [ ] **Step 5:** Commit + `FINDINGS.md` entry.

## Task 4 (Family A): SVA/RUV for the cross-survey level offset

**Files:** Create `explore/scripts/03_ruv_level_offset.R`

Why: `fe_transport_level_offset` records that LOCO transport dies on a
biomarker *level* offset (raw ferritin 6× across countries) that is not assay
(AS-01: same VitMin ELISA in all four surveys) and not adjustment method. That
is the definition of a batch effect, and RUV/SVA is the genomics answer to it.

- [ ] **Step 1:** Assemble the pooled 4-country district matrix per outcome
      with the *unstandardised* level as the target.
- [ ] **Step 2:** Estimate latent factors two ways — `sva::sva()` on the
      predictor matrix with country as the known batch, and RUV using
      predictors with no plausible nutrition pathway as negative controls
      (terrain roughness, distance-to-coast) to estimate the unwanted-variation
      subspace.
- [ ] **Step 3:** Remove the estimated factors, then score LOCO transport on
      the *level* scale — the thing the project currently cannot do at all.
- [ ] **Step 4:** Write `explore/out/03_ruv_level.csv` with, for each outcome,
      transport before and after removal, and the fraction of cross-country
      level variance the factors absorb.
- [ ] **Step 5: Verify** — the correction must not improve *within*-country
      ranking (it should be near-neutral there); if it improves both, it is
      absorbing signal and the probe says so.
- [ ] **Step 6:** Commit + `FINDINGS.md` entry.

## Task 5 (Family A): limma-moderated cross-cell pooling

**Files:** Create `explore/scripts/04_moderated_meta_index.R`

Why: 578 predictors × 24 cells is exactly limma's many-features/many-contrasts
layout. `fe_effective_n` names small effective n as the root cause of unstable
prediction; moderating each cell's per-predictor statistic toward a cross-cell
prior is the direct attack on it.

- [ ] **Step 1:** Per cell, compute the weighted marginal association of every
      predictor with the outcome and its standard error.
- [ ] **Step 2:** Fit an empirical-Bayes hierarchy across cells (limma-style
      moderated t, plus a random-effects meta-analysis per predictor) to get
      shrunk per-predictor weights — pooled over 24 cells, *nested inside the
      fold* so the held-out cell never contributes to its own weights.
- [ ] **Step 3:** Build a "meta-index" from the shrunk weights; score against
      the zero-tuning index on in-fill and LOCO.
- [ ] **Step 4:** Also report which predictors survive moderation at a pooled
      FDR — the hypothesis-generating output the RA annotation work can use.
- [ ] **Step 5:** Write `explore/out/04_moderated_{scores,weights}.csv`.
- [ ] **Step 6: Verify** — the nesting is real: assert the held-out cell's rows
      are absent from the weight estimation (print the row counts).
- [ ] **Step 7:** Commit + `FINDINGS.md` entry.

## Task 6 (Family A): multi-trait low-rank BLUP

**Files:** Create `explore/scripts/05_multitrait.R`

Why: `XO-01` found cross-*nutrient* borrowing null but same-nutrient/other-
population beating the index (0.454 vs 0.403). That is a low-rank structure
across the 24 cells asking to be modelled jointly rather than cell by cell.

- [ ] **Step 1:** Stack all outcomes for a country on shared districts; fit a
      factor-analytic multi-trait BLUP (shared latent "general deficiency"
      factor + nutrient-specific residual) with the kernel from Task 2.
- [ ] **Step 2:** Score each cell's held-out districts, comparing single-trait
      BLUP, multi-trait BLUP, and the index.
- [ ] **Step 3:** Report the estimated genetic-correlation analogue between
      nutrients — which micronutrients share sub-national structure.
- [ ] **Step 4:** Write `explore/out/05_multitrait_{scores,cormat}.csv`.
- [ ] **Step 5: Verify** — the held-out district is held out for *every* trait
      simultaneously, or the multi-trait arm leaks. Assert and print.
- [ ] **Step 6:** Commit + `FINDINGS.md` entry.

## Task 7 (Family B): mechanistic nutrient-specific features

**Files:** Create `explore/R/features_mechanistic.R`,
`explore/scripts/06_mechanistic_build.R`

Why: the index is nutrient-agnostic over 578 generic columns; the literature is
specific. `data/MapSPAM/raw/spam2010V2r0_global_{P,A}_*.csv` carries all 42
crops at 5 arcmin (the project currently uses 5 collapsed group shares).

- [ ] **Step 1:** Build the per-district crop production basket from the raw
      MapSPAM global files, clipped to the four countries' Admin-2 polygons.
- [ ] **Step 2:** Attach a food-composition table (per-crop zinc, iron,
      phytate, provitamin-A carotenoid, folate density) as an explicit CSV in
      `explore/data/food_composition.csv` with a documented source per row.
- [ ] **Step 3:** Derive, per district: **phytate:zinc molar ratio** of the
      basket (the Wessells & Brown mechanism), **provitamin-A carotenoid
      density**, **non-haem iron density and its phytate load**, **folate
      density**, and an **animal-source-food index** from GLW4 livestock plus
      inland-water access.
- [ ] **Step 4:** Write `explore/out/06_mechanistic_features.csv` (one row per
      district, ~10 columns) with a metadata sidecar naming each column's
      source and formula.
- [ ] **Step 5: Verify** — sanity-check against known agronomy: maize/rice
      districts must score high phytate:Zn, cassava districts low zinc and low
      phytate; print the top and bottom 5 districts per feature for eyeballing.
- [ ] **Step 6:** Commit.

## Task 8 (Family B): do 10 mechanistic columns beat 578 generic ones?

**Files:** Create `explore/scripts/07_mechanistic_score.R`

- [ ] **Step 1:** Score, per cell, a nutrient-*matched* mechanistic arm (zinc
      cell gets the phytate:Zn feature, B12 gets ASF, etc.) against the full
      index, the climate+soil index, and a nutrient-*mismatched* mechanistic
      arm as the specificity control.
- [ ] **Step 2:** Run in-fill and LOCO; write `explore/out/07_mechanistic_scores.csv`.
- [ ] **Step 3:** Deep-dive the strongest cell (Malawi B12, 0.70 on the record)
      — does mechanism explain why it is strongest?
- [ ] **Step 4: Verify** — the mismatched control must *not* beat the matched
      arm; if it does, the "mechanism" is a generic agro-ecology proxy and the
      entry says so plainly.
- [ ] **Step 5:** Commit + `FINDINGS.md` entry.

## Task 9 (Family D): n≪p estimators not yet tried

**Files:** Create `explore/scripts/08_np_estimators.R`

Why: the project has tried elastic net, the index, HAL/PCHAL and SuperLearner.
Sparse PLS, supervised PCA and stability selection are the standard
chemometrics/genomics answers at this n and are untried here. This task also
fixes a defect on the record: `XO-01` found the index sums *un-standardised*
axes, so any sd-1 column is under-weighted.

- [ ] **Step 1:** Arms — sparse PLS (`pls`/`spls`), supervised PCA
      (Bair–Tibshirani: marginal screen then PC), `ncvreg` MCP/SCAD, stability
      selection over the elastic net, graphical-lasso decorrelation then ridge,
      and a **standardised-axis index** (the index with each axis scaled to
      sd 1 before summing).
- [ ] **Step 2:** Score all arms on all 24 cells, in-fill and LOCO.
- [ ] **Step 3:** Write `explore/out/08_np_estimators.csv`.
- [ ] **Step 4: Verify** — the standardised-axis index is compared to the
      current index on identical folds so the difference is attributable.
- [ ] **Step 5:** Commit + `FINDINGS.md` entry.

## Task 10 (Family C, Phase 1): fine-temporal features from what is on disk

**Files:** Create `explore/R/features_temporal.R`,
`explore/scripts/09_temporal_phase1.R`

Why: fieldwork windows are 2–3 months per country, so calendar spread is not
the axis — *seasonal phase and lag history* are, and they vary across clusters
because clusters differ in agro-ecology and visit date. The record has only
`_fw`, `_fw_anom`, `_prev3`.

- [ ] **Step 1:** From `data/covariates/cluster/predictors_cluster.csv` and the
      rasters already on disk, build: **months since last harvest at the blood
      draw** (phenology phase), the **annual and semiannual Fourier harmonics**
      (amplitude and phase) of the NDVI/precip series, and a **0–24 month lag
      anomaly stack**.
- [ ] **Step 2:** Fit the lag stack under a *penalised distributed lag* (lag as
      a smooth via `mgcv`, so 24 lags cost 1–2 effective parameters, not 24).
- [ ] **Step 3:** Score at cluster level *and* aggregated to districts, against
      the cluster track's own 0.244/0.368 on the record.
- [ ] **Step 4:** Write `explore/out/09_temporal_phase1.csv`.
- [ ] **Step 5: Verify** — confirm against `metadata/survey_years.csv` that
      every lag window ends before the cluster's fieldwork date (no lookahead).
      Print the min/max lag window per country.
- [ ] **Step 6:** Commit + `FINDINGS.md` entry.

## Task 11 (Family C, Phase 2): targeted GEE extraction — **gated**

**Files:** Create `explore/scripts/10_gee_extract.py`

- [ ] **Step 1:** **Ask the user before running.** Phase 2 runs only for the
      layers Task 10 flags as promising, and only for those.
- [ ] **Step 2:** Extract at cluster buffers (2 km urban / 5 km rural, matching
      the existing convention): monthly CHIRPS/ERA5/MODIS stacks for the 36
      months before each cluster's fieldwork date, and AlphaEarth embeddings at
      the cluster buffer for the survey year.
- [ ] **Step 3:** Write to `explore/out/10_gee_cluster_fine.csv` — never to
      `data/covariates/`.
- [ ] **Step 4:** Re-score Task 10's best arms on the extracted stack.
- [ ] **Step 5:** Commit + `FINDINGS.md` entry.

## Task 12: Synthesis

**Files:** Modify `explore/FINDINGS.md`; create `explore/out/leads.csv`

- [ ] **Step 1:** One table of every probe: loose score, honest in-fill, honest
      LOCO, comparator, verdict (*pre-registrable candidate* / *dead end* /
      *needs data*).
- [ ] **Step 2:** Write the synthesis section — what to pre-register for the
      next country, what is nutrient-specific vs broad, what needs data.
- [ ] **Step 3:** State plainly which leads were chosen post hoc on these four
      countries and therefore cannot be quoted as results.
- [ ] **Step 4:** Commit.

---

## Self-review

- **Spec coverage.** Family A → Tasks 2–6; Family B → Tasks 7–8; Family C →
  Tasks 10–11; Family D → Task 9; shared harness → Task 1; the loose-screen /
  honest-score decision is built into every probe's output columns.
- **Placeholders.** None: every task names exact files, the method, the output
  table, and a concrete verification.
- **Naming consistency.** `exp_load()`, `exp_cell()`, `exp_infill()`,
  `exp_loco()`, `exp_baseline_arms()` are defined in Task 1 and used under
  those names throughout. Arm signature is
  `function(tr, te, y, X, aux) -> numeric(length(te))` everywhere.
- **Risk the plan carries deliberately.** Tasks 2–11 are each capable of
  producing a spurious win at n = 14–87 over 24 cells. That is accepted (this
  is hypothesis generation), and is why every entry carries the honest number,
  a same-folds comparator, and an explicit post-hoc flag in Task 12.
