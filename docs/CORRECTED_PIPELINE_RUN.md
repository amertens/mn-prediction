# Corrected pipeline — launch / driver note

The parallel, corrected-methods pipeline (`_targets_corrected.R`) implements the
P1–P8 fixes from [`docs/CRITICAL_REVIEW.md`](CRITICAL_REVIEW.md) **alongside**
the production pipeline. It writes to its **own store** (`_targets_corrected/`)
and never modifies `_targets.R` or the production store.

> **Compute guardrail.** While the HIV simulation is running (~6 workers), do
> **not** launch the `full` scope. The `smoke` scope (1 slice, ~7 s, sequential)
> is safe to run anytime. Launch `full` only once the machine frees.

## What it reuses (no heavy recompute)
It reads already-built **production** objects from `_targets_full/objects/`
(`outcome_data_*`, `svy_admin2_*`, `gee_admin2_*`) via `read_prod()`. It does
**not** rerun GEE extraction or data merging. Override the source store with env
`PROD_STORE` if needed (default `_targets_full`).

## Files
- `_targets_corrected.R` — pipeline definition (separate store).
- `R_corrected/00_corrected_utils.R` — fold builders, in-fold preprocessing, learners, metrics, `read_prod()`.
- `R_corrected/p1_fitting.R` — P1 leakage-free SL (honest cluster + spatial block + optimistic).
- `R_corrected/p2_p6_methods.R` — P2 OOF calibration, P3 decision value, P4 intervals, P5 sampling-aware error, P6 trust flags.
- `R_corrected/p7_area.R` — P7 design-aware partial-pooling area estimator.
- `R_corrected/comparison.R` — assembles corrected-vs-production bundle + CSVs.
- Outputs: `results/tables/corrected/*.csv` and `dashboard/data/methods_comparison.rds` (feeds the dashboard "Methods (corrected)" tab).

## Scope (env `CORRECTED_SCOPE`)
| Scope | Slices | Library | Pool cap | Use |
|-------|--------|---------|----------|-----|
| `smoke` (default) | `gambia_women_iron` | mean+glmnet+ranger | 400 | validation, ~7 s |
| `full` | all available country × outcome | +rpart | 800 | unattended publication run |

## Validate (no build — safe anytime)
```r
Rscript --vanilla -e 'targets::tar_manifest(script="_targets_corrected.R")'
Rscript --vanilla -e 'targets::tar_validate(script="_targets_corrected.R")'
```

## Smoke run (one slice, ~7 s, sequential — safe anytime)
```bash
CORRECTED_SCOPE=smoke Rscript --vanilla -e \
  'targets::tar_make(script="_targets_corrected.R", store="_targets_corrected")'
```

## Full run — UNATTENDED, only once the HIV machine is free
Keep parallelism minimal (the per-slice SL fit is already cheap; the slice count
is the multiplier). Recommended: sequential or ≤2 workers.

```bash
# sequential (lowest footprint; recommended while unsure of free capacity)
CORRECTED_SCOPE=full Rscript --vanilla -e \
  'targets::tar_make(script="_targets_corrected.R", store="_targets_corrected")'

# optional light parallelism (≤2 workers) once the machine is fully free
CORRECTED_SCOPE=full Rscript --vanilla -e \
  'library(future); plan(multisession, workers=2);
   targets::tar_make_future(script="_targets_corrected.R",
                            store="_targets_corrected", workers=2)'
```
Estimated cost: ~22 slices × a few seconds each for the SL fit at pool cap 800
(ranger dominates) → on the order of minutes, not hours. Memory is modest
(one slice in memory at a time; `memory="transient"`, `garbage_collection=TRUE`).

### Run in the background and detach (Windows PowerShell)
```powershell
$env:CORRECTED_SCOPE="full"
Start-Process -NoNewWindow -FilePath "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" `
  -ArgumentList '--vanilla','-e','targets::tar_make(script="_targets_corrected.R", store="_targets_corrected")' `
  -RedirectStandardOutput "corrected_full.log" -RedirectStandardError "corrected_full.err"
```

## Checkpoint / resume
`{targets}` is resumable by design: every target is cached in
`_targets_corrected/`. If the run is interrupted (or you change one helper),
**just rerun the same `tar_make`** — completed targets are skipped and only
stale/new ones rebuild. Useful commands:
```r
targets::tar_progress(store="_targets_corrected")     # what built / errored
targets::tar_errored(store="_targets_corrected")      # which targets failed
targets::tar_outdated(script="_targets_corrected.R",  # what would rebuild
                      store="_targets_corrected")
```
Per-slice failures do not abort the others’ cached results; fix the cause and
rerun to fill gaps.

## After the full run
1. The comparison bundle + CSVs refresh automatically
   (`results/tables/corrected/`, `dashboard/data/methods_comparison.rds`).
2. Launch the dashboard and open **Methods (corrected)**:
   `R -e "shiny::runApp('dashboard')"`.
3. See [`docs/CORRECTED_PIPELINE_RESULTS.md`](CORRECTED_PIPELINE_RESULTS.md) for
   how to read the head-to-head deltas.
