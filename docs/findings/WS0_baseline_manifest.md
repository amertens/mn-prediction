# WS0. Baseline manifest

Purpose: snapshot every headline result table before any refit in the
`analysis-updates-2026-08` branch, so that each later workstream can be diffed
against a fixed reference rather than against a moving target.

Branch: `analysis-updates-2026-08`, created from the `tanzania-audit` tip at
commit `f3d512d`.

## What was frozen

Script: `scripts/freeze_baselines.R`. The script copies; it does not move, edit
or delete any source file, and it refuses to run a second time unless
`FREEZE_FORCE=true` is set, so a snapshot cannot be silently replaced.

Destination: `results/tables/frozen_2026-08/`, grouped by source location.
Mirroring the full repository-relative path would have nested `results/tables/`
inside itself, so the tree is grouped instead and `MANIFEST.csv` carries the
original path of every file.

| Group | Source location | Files |
|---|---|---|
| `tables` | `results/tables/*.csv` | 23 |
| `tables_corrected` | `results/tables/corrected/*` | 14 |
| `sensitivity` | `results/sensitivity/*.csv` | 6 |
| `sandbox_parsimony_out` | `sandbox_parsimony/out/*.csv` | 9 |
| **Total** | | **52** |

(source: `results/tables/frozen_2026-08/MANIFEST.csv`, column `group`, counted
over all 52 data rows.)

Manifest columns: `frozen_path`, `source_path`, `group`, `why`, `md5`, `bytes`,
`source_mtime`, `frozen_date`. The `md5` is computed on the stored copy rather
than the source, so the manifest verifies what the frozen tree actually holds.
The `why` column records the reason each file was selected, which keeps the
selection rationale attached to the data rather than to this document alone.

`results/tables/corrected/` is frozen as a whole directory rather than as an
enumerated list, so a corrected-methods table added later cannot escape the
freeze by not appearing in a hand-written list.

The previously untracked file `results/tables/area_comparison_all_PREFULL.csv`
is included in the freeze and is committed on this branch.

Every file named in the Stage 0 audit was located. The script reports absent
files explicitly rather than skipping them; it reported none.

## Baselines that could not be frozen because they do not exist

Two items named in the brief have no corresponding artifact anywhere in the
repository. They are recorded here as absent rather than reconstructed.

1. **Survey-subsample cost-of-accuracy outputs.** No file matching
   `*subsample*`, `*cost_of_accuracy*`, `*retention*` or `*sentinel*` exists in
   the repository, and `results/` contains no retention-simulation table. The
   result is described in prose at
   `docs/PROJECT_STATUS_2026-08_UPDATE.md` section 5. Its figures are therefore
   **not yet computed** in the sense of this branch's anti-fabrication rule, and
   no number from that prose is quoted as a result anywhere in this work. WS6
   rebuilds the simulation as a seeded pipeline target before extending it.

2. **Weighted versus unweighted SuperLearner comparison.** Described in prose at
   `docs/PROJECT_STATUS_2026-08_UPDATE.md` section 6. The comparison was produced
   by a monkeypatched refit that was not committed, so no table exists to freeze.
   Its figures are **not yet computed** for the purposes of this branch.

Both absences are consequences of prior work having been reported in prose
without a committed artifact. They are the reason guardrail 3 exists, and they
constrain what later workstreams may cite.

## Reviewer passes

**Statistical reviewer.** WS0 fits no model, constructs no folds, and computes no
estimate, so the leakage, fold-construction, seed-coverage and survey-weight
checks have no surface to act on. The one statistical property that matters here
is that the freeze is a pure copy: the script's only file operations are
`file.copy()` with `overwrite = TRUE` into a directory it has just created, and
`write.csv()` of the manifest. No source file is opened for writing. Checked and
confirmed by reading the script.

**Reproducibility reviewer.** The script adds no `targets` node, so
`tar_manifest()` and `tar_network()` are unaffected and were not re-run for this
workstream. Paths are resolved through `here::here()`; there is no absolute path
and no `setwd()`. The script is deterministic: it performs no sampling and its
only date-dependent output is the `frozen_date` column. Re-running it without
`FREEZE_FORCE=true` stops with an error rather than overwriting, which was
verified by running the script a second time after the freeze: it exited 1 with
the guard message and wrote nothing.

## Consumed by

`scripts/regression_gate.R` (WS9) reads `MANIFEST.csv` and maps each regenerated
table back to its baseline through the `source_path` column, classifying every
change as a new scheme row, a changed baseline number, or unchanged.
