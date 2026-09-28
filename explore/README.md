# `explore/` — exploratory, hypothesis-generating probes

Started 2026-09-28, at the user's direction. **Nothing in this folder is a
production claim.** It exists to generate leads for future work: proxy
predictors and modelling approaches that might do well for specific
micronutrients, in specific countries, or broadly.

It is explicitly licensed to overfit. Probes try many models, many feature
sets, and post-hoc choices on the same four countries. That is the point — but
it means **no number in `explore/out/` may be quoted as a result**, only as a
candidate to pre-register and test on a country these models have never seen.

## Isolation contract

`explore/` creates no file outside itself. It *reads*:

- `results/tables/protocol_v2/targets_v2.csv` (the 24 country × outcome cells)
- `data/covariates/harmonized/predictors_admin2_shared.csv` (+ metadata)
- `data/`, `dashboard/data/admin2_boundaries.rds`
- `R/protocol_v2.R`, sourced read-only

It never writes to `R/`, `scripts/`, `results/`, `dashboard/`, `docs/`,
`metadata/` or `_targets.R`, and is not part of the `{targets}` DAG.

## The gate

`explore/R/harness.R` reuses the project's own cell construction,
rank-normalisation (fix 3), domain representation (fix 4), fold seeds (fix 1)
and scoring. What is new is only that the *arm* and the *feature matrix* are
arguments.

`scripts/00_harness_check.R` is the gate: it runs the protocol's own arms
through the harness and compares cell by cell against
`results/tables/protocol_v2/benchmarks_v2_cells.csv`.

**Status: CLEAN** (2026-09-28). 98.8% of 376 cell × arm comparisons within
0.02; mean |diff| 0.0009; the LOCO transport rows reproduce exactly
(median Spearman 0.287 level / 0.192 prevalence, max |diff| 0.0000). The
residual differences are all in the constant-prediction `null_train_mean` arm
under the region estimand on the prevalence target, where Spearman is
tie-dominated.

Two configuration facts the gate exposed, worth knowing before writing any
probe:

1. The headline record uses `V2_PREDICTOR_TIERS = "open,survey_public"`
   (no DHS) throughout, so `exp_load()` defaults to the same set. Without it
   a probe is quietly computed on a larger predictor set than the record it is
   compared against.
2. Under LOCO the domain axes must be rebuilt **inside the fold** on the
   pooled matrix with `sign_rows = training countries`. Building them per
   country lets each country learn its own PC1 sign, so a domain score means
   the opposite thing in two countries and transport is destroyed. (First
   version of the harness got this wrong; it cost up to 0.76 Spearman.)

## Running a probe

```bash
Rscript explore/scripts/00_harness_check.R        # the gate, run it first
Rscript explore/scripts/01_kernel_blup.R          # etc.
```

`EXP_REPS` overrides the in-fill replication count (default 10, the protocol's).

## Layout

| Path | What |
|---|---|
| `PLAN.md` | The implementation plan, task by task |
| `FINDINGS.md` | Running log: question, design, result, honest number, verdict |
| `R/harness.R` | `exp_load()`, `exp_cell()`, `exp_infill()`, `exp_region()`, `exp_loco()` |
| `R/features_*.R` | Feature builders (mechanistic, temporal) |
| `R/methods_*.R` | Estimators (kernel BLUP, etc.) |
| `scripts/NN_*.R` | One probe each |
| `out/` | Result tables |

## Reference numbers to beat

From `benchmarks_v2_cells.csv`, the domain index (the estimator of record):

| estimand | target | domain index | spatial | n cells |
|---|---|---|---|---|
| in-fill | level | 0.394 | 0.457 | 18 |
| in-fill | prev | 0.298 | 0.307 | 18 |
| region | level | 0.339 | 0.436 | 18 |
| region | prev | 0.222 | 0.240 | 18 |
| transport (LOCO) | level | 0.287 | — | 22 |
| transport (LOCO) | prev | 0.192 | — | 22 |

The pre-registered climate+soil index transports at 0.368 on the level
(`transport_domains_climate_soil`), better than the full index — chosen on
these same four countries, so it is a prediction for the next one, not a
result. Sierra Leone's 14 districts cannot be folded for in-fill, which is why
the in-fill estimands have 18 cells and transport has 22.
