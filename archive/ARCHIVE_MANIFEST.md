# Archive manifest

Files retired from the active tree during the 2026-06-25 repo cleanup. Nothing
here was deleted — everything was **moved** (reversible) and the directory
structure is mirrored, so any file can be restored to its original path.

## 2026-08-26 — Tanzania (TDHS 2010) dropped

`src/Tanzania/*` moved to `archive/src/Tanzania/*`, and the `Tanzania` entry
was removed from `get_country_configs()` in `R/config.R` (preserved, unused,
as `get_country_config_tanzania_archived_2010()` in the same file — see the
comment on that function for the full rationale and the restore procedure).

**Reason:** on the 2026-08 UC Davis / BMGF bi-weekly call, project
collaborator Omar flagged the TDHS 2010 dried-blood-spot RBP (vitamin-A)
measurements as unreliable ("completely wrong... totally incorrect... not
even internally consistent"), independent of and consistent with data-quality
anomalies already found in this pipeline (a unit-label error in the survey's
own `TZ61BIOMARKER.DOC`, and a 77-92 percentage-point national-level bias
when Tanzania was held out under leave-one-country-out evaluation). His
recommendation was to use the Tanzania 2023 DHS round instead once its data
is accessible, rather than the 2010 round.

To restore: `git mv archive/src/Tanzania src/Tanzania`, then move
`get_country_config_tanzania_archived_2010()`'s returned `Tanzania = list(...)`
back into `get_country_configs()`'s returned list in `R/config.R`.

**Nothing in this archive is referenced by the production pipeline (`_targets.R`),
the corrected pipeline (`_targets_corrected.R`), or the deployed dashboard.** This
was confirmed by a repo-wide reference grep before moving each file.

To restore a file: `git mv archive/<path> <original-path>` (tracked) or
`mv archive/<path> <original-path>` (untracked).

| Archived path | Original path | Reason |
|---|---|---|
| `_tmp_cluster_cmp.R` | `/_tmp_cluster_cmp.R` | Throwaway probe (cluster-vs-admin2 comparison scratch). `_tmp_` prefix; not sourced by any pipeline. |
| `_tmp_gen_clusterrds.R` | `/_tmp_gen_clusterrds.R` | Throwaway probe that generated a one-off cluster `.rds`. Not in any DAG. |
| `_tmp_loco_cmp.R` | `/_tmp_loco_cmp.R` | Throwaway probe (LOCO comparison scratch). Not in any DAG. |
| `_tmp_verify_oof.R` | `/_tmp_verify_oof.R` | Throwaway probe verifying out-of-fold predictions during dev. Superseded by the corrected pipeline's own checks. |
| `_run_brinda_regen.R` | `/_run_brinda_regen.R` | Ad-hoc one-off runner to regenerate BRINDA outputs. Canonical runner is the targets pipeline. |
| `_run_pipeline.R` | `/_run_pipeline.R` | Ad-hoc `tar_make` wrapper. Canonical entry points are `_targets.R` (via `tar_make()`) and `run_corrected_then_hiv.ps1`. |
| `app/app.R` | `/app/app.R` | Superseded original monolithic Shiny app (single-file, ggplot/patchwork 2×2 map grid, reads `_targets`). Replaced by the modular `dashboard/` (bslib + leaflet + modules). Referenced nowhere; not deployed. |
| `dashboard/data-raw/_build_new_bundles.R` | same | Throwaway dev shortcut; its own header says it "mirrors sections 7d & 7e of `01_prepare_dashboard_data.R`" to avoid a full prep run. Logic lives in the canonical builder. |
| `dashboard/data-raw/_calc_headline.R` | same | One-off console-print probe that computed the Ghana use-case headline numbers; those numbers are now baked into `mod_start_here.R`. |
| `dashboard/data-raw/smoke_test.R` | same | Ad-hoc dashboard smoke test (not a real test suite; never deployed — `deploy.R` does not bundle `data-raw/`). |
| `dashboard/data-raw/test_app_construction.R` | same | Ad-hoc smoke test (headless app construct). Superseded by the `testServer` verification used in dev. |
| `dashboard/data-raw/test_deploy_ready.R` | same | Ad-hoc pre-deploy smoke test. Not part of the deploy bundle. |
| `dashboard/data-raw/test_endpoints.R` | same | Ad-hoc endpoint smoke test, superseded by `test_endpoints_v2.R`. |
| `dashboard/data-raw/test_endpoints_v2.R` | same | Ad-hoc endpoint smoke test (v2). Manual dev tool, not deployed. |
| `dashboard/data-raw/test_server.R` | same | Ad-hoc server smoke test. Manual dev tool, not deployed. |

## Previously removed (before this pass)
For completeness — these were removed in earlier sessions and are **not** in this
archive:

- Four duplicate corrected-pipeline scripts that had been left inside
  `R_corrected/` and were being double-sourced by `tar_source("R_corrected")`;
  the canonical driver scripts were moved to `corrected_driver/`
  (`run_corrected_full.R`, `verify_corrected.R`).
- Several throwaway verify scripts created while building the corrected pipeline.
