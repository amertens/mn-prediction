# Manuscript review notes — `manuscript_mcn.qmd`

Critical review and corrections, 2026-06-15. The manuscript prose had been written
against **earlier** versions of the result tables and had drifted out of sync; all
headline numbers were recomputed from the current `results/tables/*.csv` and the
`.qmd` source corrected. **Re-render to refresh the `.docx`:** `quarto render docs/manuscript_mcn.qmd`.

---

## 1. Corrections already applied (factual sync to current tables)

| Claim (old → new) | Location | Source table |
|---|---|---|
| Area-SL within-country "r > 0.9 in 27 of 30" → **16 of 21 (MSE), 18 of 23 (NLL)** | abstract, Results §within, Discussion | `area_comparison_all.csv` (24 combos total, not 30) |
| Individual-SL median r 0.55 → **0.50** | abstract, Results | `area_comparison_all.csv` |
| LOCO median 0.18 → **0.21**; range −0.22→0.53 → **−0.14→0.53** | abstract, Results, Discussion | `transportability_area_loco_metrics.csv` |
| LOCO MAE "5–46 pp" → **1.5–44 pp** | abstract, Discussion | same |
| Iron median 0.31 → **0.28**; women-vitA −0.05 → **−0.12** | Results | same |
| Benchmark "9 of 16, SL won 7" → **best benchmark beat SL on 10 of 13 estimable splits (SL won 3)** | abstract, Results | `benchmarks_all.csv` |
| Negative Brier skill "6 of 24" → **13 of 23** (+ noted recalibration fixes all) | Results, Limitations | `diagnostics_binary.csv` / `_calibrated.csv` |
| National est. "within ±5pp for 22 of 24" → **all 24** (median abs diff 0.4 pp) | Results | `national_estimates_all.csv` |
| SierraLeone child-iron SL-vs-OLS anecdote 0.37 → **0.51**; fixed women-iron / child-vitA win examples | Results, Discussion | `benchmarks_all.csv` |
| Removed stale "benchmark targets still running" placeholder note | `@tbl-benchmarks` caption | table is now fully populated |

**Net effect on the thesis:** the benchmark result is now *stronger* for the paper's
argument, not weaker — the prescreened ensemble loses to the best SAE/regression
benchmark on **10 of 13** estimable LOCO splits (old text understated this as a close 9–7).

## 2. Framing / scientific edits applied

- **Effective-n caveat (new paragraph).** Results §Within-country now states the
  area-level effective n is the number of admin-2 polygons (14 / 30 / 75 / 87), not
  the underlying individuals, so the within-country CV r is an optimistic upper bound
  (near-perfect smallest-country values reflect low-DoF overfitting). Reinforced in
  §Limitations. Rationale: aligns framing with the over-parameterization mechanism the
  paper already uses to explain LOCO failure.
- **"West Africa" → "sub-Saharan Africa."** Title, short title, abstract, Introduction,
  Discussion summary, Limitations, and external-validation paragraph no longer mislabel
  all four countries as West African (Malawi is Southern/Eastern). The five legitimate
  uses of "West African geography" — where Malawi-as-out-of-region is the point — were
  intentionally kept.

## 3. New figure

- `scripts/fig_loco_scatter.R` → `results/figures/loco_scatter.png`. Observed vs
  predicted admin-2 prevalence under LOCO, 4 outcomes × 4 held-out countries, y = x
  line, point size = survey n. Panel Pearson r uses the **unweighted** definition so it
  matches `@tbl-loco` exactly (verified). The `fig-loco-scatter` chunk now points at the
  real PNG (previously a non-existent placeholder rendered nothing).

---

## OPEN ITEMS — need author decision before submission

### A. Two different `sl_prescreened` numbers exist (data hygiene) — IMPORTANT
`benchmarks_all.csv` and `sl_prescreened_main.csv` report **different** LOCO Pearson r
for the same outcome × held-out country (different `n_train`; e.g. Ghana child-iron
0.336 vs 0.435), and some combos are **missing** from `benchmarks_all` (e.g. child-vitA
Ghana, child-iron Sierra Leone has it but others don't).
- Decide which file is canonical for "the SuperLearner," and regenerate so all 16
  main combos are present from one source.
- Consequence if unfixed: the `@tbl-benchmarks` pivot renders with confusing gaps —
  `fh` covers only 10/16 splits and `bym2` only 5/16, so the printed table will be
  sparse and the head-to-head tally depends on which cells happen to be populated.

### B. Unverifiable Methods calibration "before/after" slopes
Methods reports e.g. "Malawi child vitamin A slope moved from −4.30 to −0.15." The raw
−4.30 checks out (`diagnostics_binary.csv`), but `diagnostics_binary_calibrated.csv`
reports post-calibration slope = 1.0 **by construction**, so the "−0.15" cannot be
reproduced from the available tables. Confirm where the post-calibration slope values
come from (an independent held-out recalibration assessment?) or update them.

### C. `@tbl-loco` method label — confirm it is the ensemble
The caption says "area-level SuperLearner," but the underlying
`transportability_area_loco_metrics.csv` has an `n_selected` column (5–16 vars/split)
that looks like a forward/L1-selected linear model rather than the full ensemble.
Confirm the label matches what produced the file. (The new scatter figure and the table
are mutually consistent — both come from the same predictions/metrics pair — so this is
a labeling question, not a new inconsistency.)

### D. Possible overstatement to revisit (optional)
The within-country "excellent / near-perfect" framing still leads the Results even with
the new effective-n caveat. Consider whether the headline should be tempered further
given both within- and cross-country results rest on the same 14–87 effective points.

### E. Housekeeping placeholders
- Repository URL: `<repository-URL placeholder>` (two spots: Methods §Evaluation and
  Data availability).
- Grant number: `INV-XXXXXX` (Acknowledgements).
- "Co-authors to be added on revision" (Author contributions).

### F. `@tbl-cohort` chunk fragility (low priority)
The "Observed prevalence" column is computed by string-matching `country`+`outcome`
against `diag_binary` inside an `ifelse`; it is fragile if either table's row set
changes. Consider a clean `left_join` instead.

---

## How to regenerate
```sh
Rscript scripts/fig_loco_scatter.R          # rebuild results/figures/loco_scatter.png
quarto render docs/manuscript_mcn.qmd       # rebuild docs/manuscript_mcn.docx
```
