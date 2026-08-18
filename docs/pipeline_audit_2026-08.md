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

## 2. Untracked file inputs — PARTLY FIXED

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
| `load_dhs_admin2()` via `merge_external_predictors()` | `data/DHS/clean/<country>_<year>_dhs_(custom_)admin2_wide.rds` | **FIXED** — `dhs_stamp_<country>` targets (md5, `cue = "always"`) now feed `merged_ext_*` |
| `R/cluster_aggregation.R:45,102` | GEE CSV, merged RDS | OPEN |
| `R/conceptual_ablation.R:19-38` | results CSVs | OPEN |
| `R/external_data.R`, `R/oos_prediction.R` | `data/external_cache/*` | OPEN — these are regenerable caches, lower risk |

The two highest-risk cases — files that are *regenerated* by other scripts in
this repo, rather than fetched caches — are now fixed. Re-running
`src/DHS/DHS_admin2_aggregation.R` used to leave `merged_ext_*` serving a stale
result indefinitely; it now invalidates the merge.

The remainder are left open deliberately: each fix invalidates its target and
everything downstream, and `data/external_cache/*` is regenerable on demand, so
the payoff is much lower. The stamp pattern in `_targets.R`
(md5 + `cue = tar_cue(mode = "always")`) is the template — it works where
`format = "file"` cannot, i.e. when the file is absent for some countries
(Tanzania has no `_dhs_custom_admin2_wide.rds`, and only Tanzania has a
legacy-parity CSV).

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

## 6. Cross-band summaries over non-commensurable bands — FIXED (opt-in)

`.append_gee_zonal_cols()` emits five cross-band summaries
(`_annual_mean` / `_sd` / `_min` / `_max` / `_range`) for every raster with <= 24
bands. That is correct when the bands are a temporal series of one variable
(`LST_Night_Monthly` -> annual mean night temperature). It is applied
indiscriminately, so it also runs on rasters whose bands are different
quantities. Classified from the actual band names in
`data/Sierra_Leone_GEE_rasters/`:

| Family | Bands | Verdict |
|---|---|---|
| `LST_Night_Monthly` | `2013_01_Mean` … `2013_12_Mean` | temporal — meaningful |
| `WAPOR` | 36 dekadal NPP slices | temporal — meaningful |
| iSDAsoil (x11) | `mean_0_20`, `mean_20_50`, `stdev_0_20`, `stdev_20_50` | location averaged with dispersion |
| `TerraClimate` | `aet, def, pdsi, pet, pr, ro, …` | distinct variables |
| `FLDAS` | `Evap_tavg, Psurf_f_tavg, Qair_f_tavg, …` | distinct variables |
| `LST_Night_Annual_Mean` | `Mean`, `FilledProportion` | Kelvin (~290) averaged with a 0-1 proportion |
| `Productivity` | `Gpp`, `Npp`, `Npp_QC` | averages a **quality-control flag** into the value |
| `LandCoverType` / `LandCoverLayers` | `LC_Type1…5`, class codes | arithmetic on **categorical codes** |
| `GPW_Demographic` | 77 age-sex count bins | a scaled total count |

Consequences that reach reported results:

- The summary is dominated by whichever band has the largest units.
  `gee_fldas_..._annual_mean` = 3640.9 on Sierra Leone, which is surface pressure
  (~98,000 Pa) / 27 — an elevation proxy wearing a climate label.
- Some summaries are **exact duplicates** of a real band:
  `gee_soilzinc_annual_max` is bit-for-bit identical to `gee_soilzinc_mean_0_20`
  (likewise iron, potassium), and `gee_terraclimate_2012_annual_min` is identical
  to `gee_terraclimate_2012_pdsi`.
- Duplicates split variable importance and make the lasso's choice among
  identical columns arbitrary, destabilising the per-fold selected-variable lists
  in `transportability_area_selected_vars.csv` for reasons unrelated to the data.

Scale (Sierra Leone): 543 Admin-2 `gee_` columns, **248 (46%) cross-band
summaries**, **243 columns an exact copy of another**. In the 149-variable
cross-country vocabulary, 62 (42%) are summaries.

**Fix.** `R/gee_band_semantics.R` declares band semantics per family
(`temporal` / `multivariate` / `categorical`) and prunes by NAME — deliberately
not by value, since a value-based filter would drop different columns in
different countries and reintroduce finding #1's country-dependent vocabulary.
Wired into the three predictor-selection choke points: `prune_predictor_cols()`
(individual), `build_area_loco_dataset()` and `assemble_area_transport()` (area),
always applied *after* the intersection so pruning is identical for every
country. An unclassified family keeps its summaries and is reported, so this can
never silently delete a covariate. Diagnostics:
`scripts/check_covariate_hygiene.R`.

**Gated OFF by default** (`GEE_COVARIATE_HYGIENE=true` to enable) and recorded in
`pipeline_params` so the setting is captured in the target hash rather than
living in an invisible env var.

**Measured effect** (`scripts/compare_covariate_hygiene.R`, area-level LOCO
transport recipe, 4 countries, cross-country vocabulary 149 -> 87):

| Outcome | LOCO Pearson r, v1 -> v2 | MAE (pp), v1 -> v2 |
|---|---|---|
| child_iron | 0.207 -> 0.192 (−0.015) | 22.9 -> **20.0** |
| child_vitA | 0.066 -> 0.031 (−0.035) | 10.5 -> **9.7** |
| women_iron | −0.119 -> −0.128 (−0.009) | 21.7 -> **21.3** |
| women_vitA | 0.044 -> 0.050 (+0.006) | 2.29 -> 2.34 |

Removing **42% of the covariates** moves mean LOCO rank correlation by −0.013 and
improves absolute error in 3 of 4 outcomes. With only four held-out countries
neither difference is statistically meaningful — which is the finding: **those 62
columns carry essentially no information**. They can be dropped for
interpretability and parsimony at no measurable cost, and since this project's
documented dominant failure mode is a level/calibration offset (MAE) rather than
ranking, the consistent MAE improvement is the more relevant direction.

**Recommended default: on**, flipped as a deliberate one-line change *after* the
Tanzania rebuild lands, so "Tanzania added" and "vocabulary changed" are not
confounded in the same set of numbers.

---

## 6b. merge_food_security() was not re-runnable — FIXED

**Severity: high.** Found by a Tanzania-scoped smoke run before committing to the
full rebuild; it was **not** Tanzania-specific. `merge_food_security()` failed for
**all five countries** on re-run, i.e. the next `tar_make()` would have died about
an hour in. Two independent causes:

1. **Not idempotent.** The per-country merge scripts
   (`src/<Country>/2_GW_*_data_merge.R`) now call `merge_food_security()`
   themselves, so the `*_merged_dataset.rds` files already carry `fsec_` columns.
   Re-merging made `dplyr::left_join()` suffix the collisions to `.x`/`.y`, after
   which `merge_by_admin()`'s `target[, data_cols]` selected columns that no
   longer existed (`undefined columns selected`). Fixed by dropping stale
   `fsec_` columns before re-merging.
2. **sf geometry column.** Gambia/Ghana/Sierra Leone/Malawi merged datasets carry
   an `sfc_POINT` `geometry` column. `dplyr::left_join()` converts via
   `as_tibble()`, which rejects list/matrix columns outright (`All columns in a
   tibble must be vectors`). Fixed by stashing non-vector columns for the
   duration of the joins and re-attaching by position, guarded on the row count.

Why the cached targets hid it: `merged_fsec_*` was built when the raw RDS files
did *not* yet contain `fsec_` columns, and nothing invalidated it when they
started to. This is finding #2 (untracked inputs) producing a latent break rather
than a wrong number.

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
