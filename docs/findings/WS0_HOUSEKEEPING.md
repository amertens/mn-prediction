# WS0. Housekeeping and baselines

Branch `accuracy-and-impact-2026-09`, cut from `main` at `f8b23d9`.

## What this workstream did

1. Cut the working branch.
2. Froze the 19 tables behind the claims of Sections 2.3, 3, 4, 5, 6 and 9 of
   `docs/SESSION_FINDINGS_FOR_REVIEW.md` into `results/tables/frozen_2026-09/`,
   with `MANIFEST.csv` recording file, md5, byte count, source modification
   time and freeze date.
3. Created `docs/findings/` and the claims register skeleton.

## Acceptance

| Criterion | Result | Source |
|:---|:---|:---|
| Manifest covers every Stage 0 table | 19 of 19 frozen, 0 missing | `results/tables/frozen_2026-09/MANIFEST.csv`, row count excluding header |
| Register lists every Section 3, 4, 5, 6, 9 and 11 claim | 53 rows | `docs/findings/CLAIMS_REGISTER.md`, count of rows with status `not yet tested` |
| Every register row starts at `not yet tested` | 53 of 53; no row carries `confirmed`, `revised` or `withdrawn` | same file |

## Deviations from the plan as approved

**Branch point.** The prompt named `covariate-harmonisation-and-honest-baselines`.
Stage 0 established that this branch is an ancestor of `main` and that `main` is
one commit ahead, `f8b23d9`, which adds seven lines to `.gitignore` and touches
nothing else. The branch was cut from `main` on instruction, so the gitignore fix
is retained.

**Document location.** `SESSION_FINDINGS_FOR_REVIEW.md` is at `docs/`, not the
repository root as the prompt stated.

**`admin2_reliability()` location.** `R/area_reliability.R:65`, not one of the
files the prompt listed.

## Findings recorded at freeze time

**Two claims have no locatable source table.** Section 3 rows 3.6 and 3.7, the
WHO category accuracy figures of 33 against 32 percent at Admin-2 and 56 against
35 percent at Admin-1, are stated in `docs/manuscript_mcn.qmd` around line 1000
but no committed table under `results/` reproduces them.
`results/tables/resolution_comparison.csv` carries the r and MAE arms at both
levels and is the nearest source, without a category column. WS8a regenerates
the classification from a threshold config and supersedes the need to recover
the original derivation.

**Section 3 row 3.12 has no locatable source at all.** The positive control,
earth observation predicting district education at r 0.48 to 0.71, appears in the
session document and in Section 11 claim 2, which is load bearing. A search of
`results/`, `sensitivity/`, `sandbox_parsimony/` and `docs/` found no table
containing it, and the manuscript does not state it. Under guardrail 3 it is
recorded as unsourced. WS4b re-measures education as one of its pseudo-targets,
which supplies a sourced replacement rather than recovering the original.

**Four discrepancies between the session document and the tables it cites**, all
carried into the register and left at `not yet tested`:

| Document | Live value | Source |
|:---|:---|:---|
| Section 2.4 median r_max 0.098 | 0.1315 | `frozen_2026-09/admin1_arms.csv`, column `r_max`, over 24 unique country-outcome cells |
| Section 2.4 "4 of 24 cells have no detectable signal" | 10 of 24 cells have `r_max` exactly 0.000 | same |
| Section 3 Individual SL median r 0.516 | 0.524 | `frozen_2026-09/area_comparison_all.csv`, column `pearson_r`, filter `approach == "Individual SL"` |
| Section 3 Individual SL r_share 0.92 | 1.05 | same file, column `r_share`, same filter |

The Section 2.4 pair matters most. The r_max distribution over the 24 cells is
bimodal rather than uniformly low: 10 cells sit at exactly 0.000, 11 sit below
0.05, and 9 sit above 0.55. A median of 0.13 summarises a distribution that has
no mass near its own median. WS1 treats this as the primary evidence.

## Reviewer pass, statistical

WS0 fits no model and computes no estimate, so leakage, fold construction,
prescreen placement, seed coverage and survey-weight handling are not engaged.
Three points still apply.

**Freeze completeness against the comparisons the later workstreams need.** The
guardrail requires that every later comparison be between arms scored on
identical cells and identical folds. The freeze was chosen to make that
checkable rather than to cover every table under `results/tables/`. It includes
both sides of each planned comparison: `admin1_arms.csv` for WS2,
`individual_anchor.csv` and `protocol_reconciliation.csv` for WS3,
`national_vmnis_ceiling.csv` and `national_composition.csv` for WS6, and
`cv_honesty_compare.csv` for the Section 2.3 optimism figure that WS8c must cite.
Tables not frozen are those no workstream regenerates. WS9 will fail its
regression gate if a later workstream writes to a table absent from the manifest,
which is the intended detection path.

**The frozen values are not independently verified.** Freezing records the
starting state; it does not endorse it. Every value in the register's "As stated"
column is a transcription and is marked as carrying no independent authority.

**No stochastic call was made.** No seed is required.

## Reviewer pass, reproducibility

**Targets graph.** WS0 adds no target. `scripts/accuracy_impact/ws0_freeze.R` is
a standalone script, consistent with the existing `scripts/covariates/` scripts
that write to `results/tables/` outside the graph. `tar_manifest()` and
`tar_network()` are therefore unchanged, and were not re-run for this workstream.
WS7a adds the first new target and will exercise both.

**Stamp targets.** The freeze reads only files already tracked in the repository.
No untracked input is consumed, so no stamp target is required.

**Paths.** The script uses `here()` throughout, contains no absolute path and no
`setwd()`. It creates `results/tables/frozen_2026-09/` and mirrors the source
directory structure beneath it, so `corrected/protocol_reconciliation.csv` does
not collide with a same-named table at the top level.

**Determinism.** Re-running the script overwrites the frozen copies with
byte-identical content and rewrites `MANIFEST.csv` with the same md5 values and a
refreshed `frozen_date`. It is safe to re-run and does not touch any file under
`results/tables/` outside `frozen_2026-09/`.

**Environment.** Run with `R_USER` and `HOME` set to the OneDrive Documents
directory so `~/.rdhs.json` resolves, as the session document requires. R 4.4.2.
The only package dependencies are `here` and `tools`.

**Smoke profile.** Not applicable. WS0 has no compute-bearing step, and the
`PROFILE=smoke` cell of Ghana `child_iron` is introduced in WS7a.
