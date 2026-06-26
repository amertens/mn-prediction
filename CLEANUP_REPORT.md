# Repo cleanup report — 2026-06-25

Cleanup of the `mn-prediction` repo under two hard rules: **behavior-preserving**
(no change to analysis/pipeline logic, parameters, or outputs — the production
`_targets.R` and the parallel `_targets_corrected.R` remain the analysis of
record and validate byte-identically) and **archive, don't delete** (retired
files moved into `archive/`, reversible). All writes were atomic; no interactive
prompts.

## Outcome at a glance
- **15 throwaway/superseded files archived** into `archive/` (mirrored structure)
  — see [`archive/ARCHIVE_MANIFEST.md`](archive/ARCHIVE_MANIFEST.md).
- **0 production-logic edits.** The active analysis code is already
  header-documented; the apparent "commented-out code" is intentional (formulas,
  citations, decision notes, deliberately-disabled options). Per the conservative
  mandate it was left intact.
- **Both pipelines re-validate with identical DAGs** (498 / 293 targets, matching
  pre-cleanup signatures); the dashboard still constructs headlessly and passes a
  `testServer` smoke.

---

## 1. Inventory & classification

| Area | Classification | Action |
|---|---|---|
| `_targets.R`, `_targets_corrected.R`, `_targets.yaml` | ACTIVE (analysis of record) | Untouched |
| `R/` (31 files) | ACTIVE — all `tar_source("R/")`d by the production pipeline | Untouched (see §3) |
| `src/0-functions.R` | ACTIVE — `source()`d by `_targets.R` | Untouched |
| `R_corrected/` (5), `corrected_driver/` (2) | ACTIVE — corrected pipeline + its driver chain | Untouched (already clean/headed) |
| `dashboard/` (app.R, app_public.R, global.R, `R/` modules) | ACTIVE — deployed app | Untouched |
| `dashboard/data-raw/` active builders (`01_prepare…`, `02_gbd_placeholder`, `_build_area/fh/bym2/sl_admin3_layer`) | ACTIVE — provenance of shipped bundles | Untouched (all already headed) |
| `simplified subset/` | ACTIVE — built this session | Untouched (clean) |
| `_tmp_*.R`, `_run_*.R` (top level) | THROWAWAY probes / ad-hoc runners | **Archived** |
| `app/app.R` | SUPERSEDED by `dashboard/` (original monolithic app) | **Archived** |
| `dashboard/data-raw/_build_new_bundles.R`, `_calc_headline.R` | THROWAWAY (mirror main builder / console probe) | **Archived** |
| `dashboard/data-raw/{smoke_test,test_app_construction,test_deploy_ready,test_endpoints,test_endpoints_v2,test_server}.R` | SMOKE/THROWAWAY dev tests (not deployed) | **Archived** |
| `dashboard/data-raw/_rebuild_admin2_predictions.R` | Borderline (shortcut mirroring canonical builder; possible provenance) | **Flagged, left** |
| `src/` legacy country/DHS/GEE scripts (~60 files) | SUPERSEDED-in-spirit but = upstream data-prep / provenance, regenerates pipeline inputs | **Flagged, left** |
| `national_prediction/` | Separate sub-pipeline (own `_targets.R`); referenced by validation docs | **Flagged, left** |
| `sandbox/`, `sandbox_fe/` | Designated scratch; some are methodology provenance (e.g. `sae_sl_hybrid_prototype.R`) | **Flagged, left** |
| `scripts/` (figure gen, data pulls, benchmark builders) | Utilities / provenance, not throwaway | **Flagged, left** |
| `mn_proxy_tutorial/` | Self-contained teaching artifact | **Left (keep)** |
| `metadata/`, `data/` | Data, not scripts | Untouched |
| Top-level `*.log`, `*.out`, `*.err`, `*.pid`, `.corrected_done` | Build artifacts/logs | **Flagged** (see recommendations) |

When availability/activeness was uncertain, the file was **flagged, not moved**.

---

## 2. Archived files (reasons)

Full per-file reasons in [`archive/ARCHIVE_MANIFEST.md`](archive/ARCHIVE_MANIFEST.md).
Summary: 6 top-level throwaway probes/runners (`_tmp_*`, `_run_*`), the legacy
monolithic `app/app.R`, and 8 `dashboard/data-raw/` throwaway helpers + smoke
tests. None are referenced by either pipeline or the dashboard (verified by a
repo-wide reference grep before moving). Tracked files were moved with `git mv`
(history preserved); untracked files with `mv`. `app/` is now empty.

---

## 3. Simplification & documentation

**No active scripts were edited.** Rationale:

- Every active script already carries a concise header (purpose/inputs/outputs);
  the new `R/` files (`brinda_adjustment.R`, `cluster_aggregation.R`,
  `feature_engineering_constructs.R`, `ingest_new_country.R`) use the full banner
  header style, and the dashboard builders/modules are documented.
- A heuristic scan flagged ~81 "commented-out code" lines, but inspection of the
  two largest concentrations (`R/benchmark_models.R`, `R/sensitivity/mlr3_fitting.R`)
  showed these are overwhelmingly **intentional**: math formulas
  (`v_i = p_i*(1-p_i)/n_svy_i`), literature citations (Stoltzfus 2003, West 2017),
  dated decision notes, and deliberately-disabled learner configs documented for
  provenance. Removing them would degrade documentation, not improve it.
- The only genuinely dead lines found were two commented debug `cat()` calls in
  `R/sensitivity/mlr3_fitting.R` (~614–615). Editing the analysis of record for
  two trivial lines is not worth the risk surface, so they are flagged, not
  removed.

This is the conservative outcome the rules call for: the cleanup value here is in
decluttering the active tree (archiving), not in churning a well-documented,
reproducible analysis codebase.

---

## 4. Verification (behavior preservation)

Captured a baseline before any change, then re-checked after archiving. The DAG
signature is `digest()` of the sorted `(target name, command)` manifest.

| Pipeline (full mode) | Baseline | After cleanup | Result |
|---|---|---|---|
| `_targets.R` | PASS · 498 targets · `5bbcc8ae…` | PASS · 498 · `5bbcc8ae…` | **identical** |
| `_targets_corrected.R` | PASS · 293 targets · `a9a4f042…` | PASS · 293 · `a9a4f042…` | **identical** |

- `tar_validate()` passes for both; `tar_manifest()` DAG signatures match
  baseline exactly → no target added/removed/changed.
- Dashboard: `app.R`, `app_public.R`, `global.R` all parse; `global.R` sources
  all modules and loads data headlessly; `testServer(mod_map_explorer_server)`
  runs (Ghana / women_iron / Admin-2 → 260 rows).
- No heavy `tar_make` or model refits were run (not needed — the DAG is unchanged
  and no sourced file was touched).

---

## 5. Recommended next cleanup (separate, reviewed passes)

1. **`src/` legacy tree (~60 files).** Clear versioned duplicates exist
   (`GW Ghana SL.R` vs `…_V2.R`, `…_bin` vs `…_bin_V2`, `2_GW_Ghana_data_merge.R`
   vs `…_admin2_V2.R`). Worth a dedicated pass that first confirms which merge
   scripts still regenerate current pipeline inputs, then archives the rest.
   Keep `src/0-functions.R` (production dependency).
2. **Build artifacts / logs at top level.** `pipeline_full_rerun.log` (~1.3 GB),
   `pipeline_full_rerun_initial.log` (~3.4 GB), `pipeline_rebuild.{log,out,err}`,
   `_brinda_regen.out` (14 MB), `_bw_regen.out` (26 MB), `_pointcov.{out,err}`,
   `inla_install.log`, `.Rhistory`, `.corrected_driver.pid`. The `*.log` files are
   already git-ignored; the `*.out/*.err` are not. Recommend adding `*.out` /
   `*.err` (or the specific names) to `.gitignore` and deleting the large
   regenerable logs rather than archiving them.
3. **`dashboard/data-raw/_rebuild_admin2_predictions.R`.** Confirm the canonical
   `01_prepare_dashboard_data.R` produces an identical `admin2_predictions.rds`;
   if so, archive the shortcut.
4. **`sandbox/` and `sandbox_fe/`.** Decide whether retained prototypes are
   provenance worth keeping or scratch to archive; they are already segregated.
5. **`national_prediction/` sub-pipeline.** Clarify whether it is active or
   superseded by the main pipeline's OOS/transport targets; document its status
   in its own README.
