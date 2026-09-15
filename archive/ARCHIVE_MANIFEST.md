# Archive manifest

Files retired from the active tree during the 2026-06-25 repo cleanup. Nothing
here was deleted — everything was **moved** (reversible) and the directory
structure is mirrored, so any file can be restored to its original path.

## 2026-08-30 — superseded status, review and exploratory notes

28 Markdown files moved from `docs/` (and `CLEANUP_REPORT.md` from the root)
to `archive/docs/`. Nothing was deleted; every file is recoverable with
`git mv archive/docs/<file> docs/<file>`.

**Reason:** `docs/` held 46 Markdown files, most of them point-in-time session
records rather than reference material, which made the directory hard to
navigate for the two things it now needs to support — running the analysis
pipeline and generating the slides and manuscript. Three categories moved:

- **Superseded status and review snapshots.** `PROJECT_STATUS_2026-06`,
  `PROJECT_STATUS_2026-08_UPDATE`, `PROJECT_STATUS_2026-09_UPDATE`,
  `PROJECT_SUMMARY`, `REVIEW_AND_ROADMAP_2026-06`, `REVIEW_NOTES`,
  `CRITICAL_REVIEW`, `CLEANUP_REPORT`, and three dated `pipeline_review_*`
  files. Each describes a state the repository has since moved past.
- **Completed one-off work plans.** `CORRECTED_PIPELINE_PATCH_PLAN`,
  `CORRECTED_PIPELINE_RUN`, `CORRECTED_PIPELINE_RESULTS`,
  `METHODS_TODO_IMPLEMENTATION_PLAN`, `DIAGNOSTICS_CHANGELOG`,
  `SHAP_KERNEL_RUN`, `audit_analysis_code`, `audit_data_cleaning`.
- **Exploratory findings now folded into the manuscript.**
  `EXPLORATORY_ANALYSES`, the five `exploratory_*` notes,
  `calibration_ladder_and_mimi_divergence`, `ghana_geostat_pilot`,
  `gee_admin3_audit_findings`.

Seven of these are still cited from live code as breadcrumbs
(`R/corrected/p7_area.R`, `p10_nested_loco.R`, `p12_distributional.R`,
`p14_subsample.R`, `src/GEE/extract_gee_cluster.R`,
`docs/manuscript_mcn.qmd`, `simplified subset/methods/aggregate_inference.R`).
Those references were repointed to `archive/docs/` in the same commit, so no
comment now names a path that does not exist.

**Kept in `docs/`** — the 19 files that are reference material rather than
history: data acquisition (`API_DATA_SOURCES`, `manual_data_downloads`,
`biomarker_survey_acquisition`, `micronutrient_survey_candidates`,
`gee_landscape_new_domains`, `ra_new_data_sources_specs`, `RA_data_roadmap`,
`RA_tasks_temporal_alignment`), methods specifications
(`METHODS_CORRECTED`, `AREA_LEVEL_RECIPE_SPEC`, `methods_mrp`,
`transportability_loco_methods`, `feature_engineering_and_sae_notes`),
standing critiques (`survey_weighting_critique`, `dc_h2_brinda_validation`,
`pipeline_audit_2026-08`), `add_new_country`, `LITERATURE_REVIEW`, and the
live `DASHBOARD_ROADMAP`.

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

## Retired 2026-09-02

| File | Archived as | Why |
|---|---|---|
| `dashboard/R/mod_scenarios.R` | same | The "Scenarios" dashboard tab. Retired when the navigation was restructured. Its projections multiplied this project's district estimates by literature-default coverage and effect sizes that were not estimated here, so the output restated the user's own assumptions; its "what-if explorer" mode scaled prevalence by an arbitrary user-chosen shift, which read as a forecast. Replaced by `dashboard/R/mod_survey_planning.R`, which answers the budget question the audience actually has using measured design-curve data from `scripts/accuracy_impact/ws5_anchoring_budget.R`. Restore by moving it back into `dashboard/R/` and re-adding `mod_scenarios_ui`/`mod_scenarios_server` to `dashboard/app.R`. |

## Previously removed (before this pass)
For completeness — these were removed in earlier sessions and are **not** in this
archive:

- Four duplicate corrected-pipeline scripts that had been left inside
  `R_corrected/` and were being double-sourced by `tar_source("R_corrected")`;
  the canonical driver scripts were moved to `corrected_driver/`
  (`run_corrected_full.R`, `verify_corrected.R`).
- Several throwaway verify scripts created while building the corrected pipeline.

## 2026-09-15 — retired modelling code and side experiments

Executed after review of `ARCHIVE_MANIFEST_2026-09-15_PROPOSED.md` (that draft
is folded into this entry and removed). Everything below was **moved** with
`git mv` and is restorable with `git mv archive/<path> <path>`; only regenerable
caches, run logs and stale markers were deleted. Every item was checked for
callers over `_targets.R`, `R/`, `scripts/`, `dashboard/` and `docs/*.qmd`
before moving; none had any.

### Moved

| Original path | Archive path | Why |
|---|---|---|
| `R/bootstrap.R` | `archive/R/bootstrap.R` | cluster bootstrap replaced by conformal intervals; `run_bootstrap_ci()` had no call sites (`B_boot` in `R/config.R` is its only reader and is documented as a no-op) |
| `R/sensitivity/sl_fitting.R` | `archive/R/sensitivity/sl_fitting.R` | the `sl3` fitting path; both DAG fits use `fit_mlr3_models()` |
| `shap_kernel/`, `run_shap_kernel.ps1` | `archive/shap_kernel/`, `archive/run_shap_kernel.ps1` | June 2026 SHAP-kernel experiment; per-district drivers now come from protocol-v2 importance |
| `sandbox_lsff/` | `archive/sandbox_lsff/` | Tang et al. external check; result recorded in `docs/findings/EXTERNAL_CHECK_TANG_LSFF_2026.md` and script 59 |
| `national_prediction/` (9 tracked files; its 55 MB `_targets/` store deleted) | `archive/national_prediction/` | national sub-pipeline superseded by `R/national_vmnis.R` / `R/national_covariates.R` (promoted 2026-08-31) |
| `src/analysis/` (8 files) | `archive/src/analysis/` | pre-`targets` analysis scripts 01–05 + `sl_helpers.R`; every step has an `R/` module |
| `src/combined results/` (15 files) | `archive/src/combined results/` | result compilers and VIM plots for the pre-`targets` outputs |
| `src/test_updated_SL_admin2_code.R`, `src/tutorial on spatial SL.R`, `src/project to_do.R`, `src/spatial_ensemble_example/` | `archive/src/…` | sl3-era scratch, an empty stub, a to-do list, an unused worked example |
| `external docs/aaai12080_am.pdf` | `archive/external docs/` | a reference PDF, not an input |

### Deleted (git-ignored, regenerable)

- `_targets_corrected/` (36 MB): store of the former parallel corrected
  pipeline, whose script is already in `archive/` and whose section now lives
  in `_targets.R`.
- `.corrected_done`, `.corrected_driver.pid`: markers of the June corrected-driver run.
- Root `pipeline_full*`, `pipeline_reported*`, `pipeline_tier1*` logs (36 MB,
  29–30 August full-mode runs): moved, not deleted, to `logs/archive_2026-08/`.

### Held (conditional in the proposal; still in place)

- `data/covariates/harmonized/*.pre_*`, `data/covariates/cluster/*.pre_*`
  (17 files, 19 MB): delete after the LK-02 predictor audit is signed off;
  `pre_LK02/` is the rollback snapshot.
- `.claude/worktrees/silly-swirles-b31d3f/`: removed later the same day
  (`git worktree remove --force`, branch `mert/silly-swirles-b31d3f` deleted;
  it had no commits beyond `signal-audit-and-protocol-v2` and its uncommitted
  edits to `R/area_level_comparison.R` / `R/area_weighted_sl.R` were confirmed
  dead by the owner).
- Ask-first items (`sandbox_parsimony/`, `mn_proxy_tutorial/`,
  `simplified subset/`, `scripts/run_full_pipeline.R`, `scripts/06–08_*.R`,
  `src/dm/`, `src/data_availability_heatmaps.R`, `src/vmnis_heatmaps.R`): untouched.

### Not moved, and why

`data/IPD/<Country>/<Country>_merged_dataset.rds` and
`src/<Country>/2_GW_<Country>_data_merge.R` (outcome source of record for the
DAG and protocol v2); `src/0-functions.R` (sourced by `_targets.R:175`),
`src/0-SL-setup.R`, `src/country_workflow_utils.R`, `src/<Country>/*_workflow.R`,
`src/DHS/`, `src/GEE/`, `src/IHME/`, `src/fuzzy_match_admin2.R`,
`src/food_price_merging.R`; `R/sensitivity/mlr3_fitting.R`,
`R/sensitivity/gp_sensitivity.R`, `sensitivity/` (active DAG sensitivity);
`R/corrected/` (wired into `_targets.R`); `dashboard/`, `harness/`,
`scripts/protocol_v2/`, `scripts/covariates/`, `scripts/accuracy_impact/`,
`scripts/cluster_level/`, `scripts/policy_deck/`, `scripts/signal_probes/`.

README.md and CLAUDE.md were updated to point at the archive paths.
