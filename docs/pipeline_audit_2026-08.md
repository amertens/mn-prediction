# Pipeline audit — 2026-08-18

Scope: correctness, reproducibility and readability of the analysis pipeline
(`_targets.R`, `R/`, `sensitivity/`), carried out alongside the incorporation of
Tanzania as a fifth country. Findings are ordered by consequence.

Everything marked **FIXED** has been changed in this repo; everything marked
**OPEN** is a real finding that was deliberately not changed, with the reason
given.

---

## 1. Silent covariate deletion in every pooled / LOCO analysis — FIXED

**Severity: high.** `build_pooled_dataset()` (`R/transportability.R`) and
`build_area_loco_dataset()` (`R/area_level_comparison.R`) select covariates by a
strict **intersection of column names** across countries. One country that names
a domain differently, or lacks it, deletes that domain **for every country** — no
error, no warning, and the models still fit and still report metrics.

Measured on this repo before the fix (`build_pooled_dataset()` on the stored
`merged_*` targets, outcome `child_vitA`; area-level from the `gee_admin2_*`
targets):

| | 4 countries | + Tanzania |
|---|---|---|
| Individual-level pooled proxy covariates | **585** | **129** |
| — of which `gee_` | 419 | **0** |
| — of which `MAP_` | 30 | 0 * |
| — of which `ihme_` | 116 | 109 |
| — of which `fsec_` | 20 | 20 |
| Area-level (Admin-2) LOCO GEE covariates | 149 | **0** |

\* `MAP_` reads 0 only because the cached `merged_tanzania` predates the Malaria
Atlas merge; the rebuilt dataset carries all 30. The new per-domain log is what
made that distinguishable from a genuine naming mismatch.

Two contributing causes, both now addressed (see `docs/add_new_country.md` §1.1
and `docs/transportability_loco_methods.md` §4.1):

- Tanzania's covariates were extracted through the Earth Engine API and named
  after **EE bands**; the other four are named after **raster filenames**.
  Overlap: 0 of 149. Resolved by `scripts/build_gee_legacy_parity.R`.
- The failure was undetectable from the outputs. `build_pooled_dataset()` now
  logs per-domain, per-country covariate counts and raises a warning when a
  domain contributes zero pooled predictors despite being present in some
  countries.

**Also newly visible:** the pooled model is dominated by whichever country has
the largest survey. For `child_vitA` the pooled sample is Tanzania 6,238 rows vs
486-1,165 for each of the others, and rows are not country-weighted — so in every
LOCO fold that holds out one of the small surveys, Tanzania supplies ~2/3 of the
training rows. Worth a country-weighted sensitivity run.

---

## 2. Untracked file inputs — PARTLY FIXED, mostly OPEN

**Severity: high for reproducibility.** `{targets}` only re-runs a target when a
*tracked* input changes. A `readRDS()` / `read.csv()` inside a function body is
invisible to it, so editing that file on disk leaves the pipeline serving a stale
cached result indefinitely.

The repo already guards the merged datasets this way (`path_merged_*`, with an
explanatory comment at `_targets.R:407`). The same hazard applies elsewhere:

| Reader | File | Status |
|---|---|---|
| `.append_legacy_parity_cols()` (`R/admin2_analysis.R`) | `data/GEE/<ISO3>_legacy_parity_admin2_gee.csv` | **FIXED** — `gee_parity_stamp_<country>` targets (md5, `cue = "always"`) now feed `gee_admin2_*` and `area_covariates_*` |
| `.append_gee_zonal_cols()` | `data/<Country>_GEE_rasters/*.tif` | OPEN |
| `load_dhs_admin1/2()` (`R/data_prep.R:212-285`) | `data/DHS/clean/*.rds` | OPEN |
| `R/cluster_aggregation.R:45,102` | GEE CSV, merged RDS | OPEN |
| `R/conceptual_ablation.R:19-38` | results CSVs | OPEN |
| `R/external_data.R`, `R/oos_prediction.R` | `data/external_cache/*` | OPEN — these are regenerable caches, lower risk |

Not fixed wholesale because each fix invalidates its target and everything
downstream; doing all of them at once would force another full rebuild. The
stamp pattern now in `_targets.R` (md5 + `cue = tar_cue(mode = "always")`) is the
template — it works where `format = "file"` cannot, i.e. when the file is absent
for some countries.

---

## 3. NA-propagating population filter — FIXED

**Severity: medium (latent).** `R/data_prep.R:385` and `R/transportability.R:47`
filtered with

```r
d <- d[d[[pop_col]] == oc$child_flag_val, ]
```

A missing population flag makes the comparison `NA`, and `d[NA, ]` **injects an
all-NA row** rather than dropping it. `R/cluster_aggregation.R:112` already had
the correct form. Both sites now match it:

```r
d <- d[!is.na(d[[pop_col]]) & d[[pop_col]] == oc$child_flag_val, , drop = FALSE]
```

No country currently has NA population flags, so results are unchanged — this is
a trap for the next dataset, not an active error.

---

## 4. Vitamin A inflammation adjustment differed silently by country — FIXED

**Severity: high for the transport claim.** `brinda_rbp_cols()` had no Tanzania
entry, so `apply_brinda_vita_binary()` fell through a warning path and kept the
configured binary. The result was that Tanzania used the DHS-supplied
`rbpadcrp` while the other four used two-marker BRINDA (CRP + AGP) — a
methodological difference presented as a harmonized outcome.

Fixed by making the difference explicit rather than incidental:

- `brinda_adjust_rbp()` gained a **CRP-only** mode (`agp = NULL`) for surveys
  with no AGP.
- `brinda_rbp_cols()` returns a named spec per population with an optional
  `fallback` column, used for rows where the inflammation marker was not
  assayed — so corrected and uncorrected values are not silently mixed in one
  outcome (TDHS 2010 assayed CRP on ~27% of the RBP sample).
- `brinda_country_method()` reports what each country actually received, and the
  method is printed on every run.

Given that the documented dominant LOCO failure mode is a **level offset**, an
undeclared difference in inflammation adjustment is exactly the kind of artifact
that would be misread as biology.

---

## 5. Duplicated outcome-derivation logic — FIXED

The BRINDA VAD derivation existed twice, in `compute_svy_admin2()`
(`R/admin2_analysis.R`) and `apply_brinda_vita_binary()`
(`R/brinda_adjustment.R`), with slightly different logging and error handling.
Two copies of a definition that must agree across the survey-weighted Admin-2
path and the LOCO pooling path is a live drift risk. Both now call one
`brinda_vad_binary()`.

---

## 6. Cross-band summaries over non-commensurable bands — OPEN (documented)

`.append_gee_zonal_cols()` emits `_annual_mean` / `_annual_sd` / `_annual_min` /
`_annual_max` / `_annual_range` **across bands**, whether or not the bands are
commensurable. For the iSDAsoil rasters that averages depth means together with
their standard deviations; for FLDAS it averages distinct physical variables.
The resulting columns are not physically meaningful, but they are in the shared
cross-country vocabulary and are fed to the models for all countries.

Not changed: removing them would alter every existing result. The legacy-parity
extractor reproduces them deliberately, with the reason stated in
`.parity_band_summaries()`, so a new country remains poolable. Worth revisiting
as a covariate-hygiene pass.

---

## 7. Reproducibility checks that came back clean

- **No absolute paths or `setwd()`** in `R/`, `_targets.R` or `sensitivity/`
  (only in comments showing example invocations).
- **No `rm(list = ls())`** in pipeline code (present only in standalone `src/`
  scripts, where it is intentional).
- **Randomness is seeded** at every site that uses it: `set.seed()` precedes each
  `sample()` / `cv.glmnet()` / bootstrap draw, with fixed literal seeds in the
  area-level and benchmark paths.
- **No LOCO label leakage** in `R/transportability_area.R`: median imputation
  (`.tr_prep_X`), correlation screening (`.tr_screen`) and model fitting all use
  training rows only. Within-country centering (`.tr_center_by`) touches held-out
  **covariates** but never held-out labels — permitted domain adaptation under
  the protocol in `docs/transportability_loco_methods.md` §1, and disabled by
  default in `AREA_TRANSPORT_RECIPE` (`center = FALSE`).

---

## 8. Housekeeping

- The working tree is on a **detached HEAD** with ~50 modified/untracked files.
  Nothing here is lost, but the work is not on a branch — worth committing.
- `pipeline_full_rerun.log` (1.3 GB) and `pipeline_full_rerun_initial.log`
  (3.4 GB) sit in the repo root. They are gitignored (`*.log`) so they will not
  be committed, but they are 4.7 GB of local disk.
- A minor readability note: `%||%` is defined twice
  (`R/corrected/00_corrected_utils.R` and `R/transportability_area.R`) and shadows
  `rlang::%||%` with NA-aware semantics. The two definitions are identical and
  the divergence from `rlang` is commented, so this is intentional; flagged only
  because a third, differing copy would be hard to notice.

---

## Verification performed

- `targets::tar_manifest()` builds cleanly with all five countries (890 targets);
  `tar_network()` confirms the new `gee_parity_stamp_*` dependency edges reach
  `gee_admin2_*` and `area_covariates_*`.
- Every edited file parses (`parse()`).
- `brinda_adjust_rbp()` CRP-only mode unit-checked: correction is monotone
  (never lowers RBP), rows with missing CRP are left untouched, and the
  `fallback` path is exercised on real Tanzania data.
- `.append_legacy_parity_cols()` checked for correct join, NA on unmatched areas,
  and graceful no-op when the CSV is absent.
- The legacy-parity extractor is validated against Sierra Leone's raster-derived
  values: **135/135 columns** pass r >= 0.9 and mean within 10%
  (`results/sensitivity/gee_legacy_parity_validation_sierraleone.csv`).
- `scripts/validate_merged_datasets.R` now includes Tanzania and passes
  (16,357 rows x 504 cols; per-prefix NA rates in line with the other countries).
