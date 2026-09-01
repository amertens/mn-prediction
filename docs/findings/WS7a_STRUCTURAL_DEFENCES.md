# WS7a. Structural defences

Three guards, built before the analytic workstreams because Section 8 of
`docs/SESSION_FINDINGS_FOR_REVIEW.md` records ten instances of one defect class
and Section 7 records three iterations of another.

The workstream found an eleventh instance of the leak class. It is described in
full below and is the reason this file is longer than a housekeeping note.

---

## 1. The headline finding: two haemoglobin columns inside the no-blood-draw arm

**`gw_wm_whbc` and `gw_gchb` are haemoglobin measurements that
`is_biomarker_column()` did not block.** Both are Ghana columns. Both were
available to the questionnaire arm of the individual anchor, the arm whose
stated purpose is to answer what a household survey *without* a blood draw buys.

### How they were found

Not by reading the guard. By the ranking the guard exists to survive. The
leakage report ranked every surviving predictor against the outcome in all 24
cells and both predictor sets, and `gw_wm_whbc` came out at the top of the whole
run: |r| 0.500 against Ghana women's iron deficiency, on 981 complete pairs,
with the second-ranked column in that cell at 0.183
(source: `results/tables/frozen_2026-09/../leakage_report.csv` as regenerated,
columns `abs_r` and `n_pairs`, filter `country == "Ghana"`,
`outcome == "women_iron"`, `predictor_set == "questionnaire"`). A column
standing 0.317 clear of everything else in its cell is the signature Section 7.2
describes.

### What they are

Measured, not inferred from the names:

| Column | n | min | median | max | distinct values | mean by outcome |
|:---|---:|---:|---:|---:|---:|:---|
| `gw_wm_whbc` | 981 | 7.0 | 12.9 | 17.0 | 74 | 12.98 not deficient, 10.71 deficient |
| `gw_gchb` | 1159 | 6.2 | 11.4 | 15.0 | 82 | not computed |

Both distributions are haemoglobin in g/dL. The women's column centres at 12.9
and the children's at 11.4, which is the expected ordering and the expected
separation. The name of the first decodes as woman-haemoglobin-concentration and
the second as child-haemoglobin.

### Why the guard missed them

`is_biomarker_column()` matched haemoglobin with `grepl("Hb", cols)`,
**case-sensitively**. The comment in `R/data_prep.R` states the reason: a
case-insensitive match would also catch the lower-case household block,
`gw_hBuy*`, `gw_hBreadType`, `gw_hBirds*`, where `h` is a household prefix
followed by an upper-case word. That reasoning is sound and the consequence was
that `whbc` and `gchb`, which spell the token in lower case, passed.

### The fix

A second pattern, `^gw_.*hb`, case-sensitive and scoped to the survey domain.
Verified against all 4,906 columns across the four surveys: it matches exactly
four, `gw_bs2_hb`, `gw_bs3_hba1c`, `gw_gchb` and `gw_wm_whbc`, every one of them
draw-derived. It matches no column of the household block, and every exposure
and identifier listed in Section 7.4 as deliberately preserved remains
preserved. The `^gw_` scope keeps it away from the MAP and map2 sickle-cell and
HbC allele-frequency rasters, which are geospatial covariates rather than
measurements taken from these respondents.

Guarded column counts across the 4,906: 65 in the survey domain and 12 outside
it, 77 in total, against 75 before the fix.

### What it changes, and what it does not

Removing the column drops that cell's maximum eligible correlation from
**0.5000 to 0.1830**, and the maximum across all 48 cell-and-set combinations
from **0.5000 to 0.3881** (source: `results/tables/leakage_report_summary.csv`,
column `max_abs_r`, before and after the fix).

The scope is narrower than it first appears, and stating the limit matters:

- **The primary models are unaffected.** They use `Xvars`, which is
  `Xvars_all` minus the entire `GW` domain (`R/data_prep.R`, `Xvars_no_gw`). No
  `gw_` column can ever enter it, so no production estimate in this project used
  either haemoglobin column.
- **The individual anchor's questionnaire arm is affected.** That arm uses
  `Xvars_full` with the guard applied, so it saw both columns. Section 5's
  questionnaire results for Ghana are therefore built on a predictor set that
  included a blood measurement. WS3 rescores that arm.
- **Ghana only.** The name-independent distribution scan found no comparable
  column in Gambia, Sierra Leone or Malawi, and no other cell's top-ranked
  column resembles haemoglobin.

Claims register rows 5.2, 5.4 and 5.7 are affected. They remain at `not yet
tested` because WS3 supplies the replacement numbers; WS7a establishes only that
the existing ones were computed on a contaminated arm.

### Residual risk this does not close

The name-independent scan looked for a haemoglobin-shaped distribution. It
surfaced six Gambia columns that are anthropometry, `gw_cWeight` and `gw_weight`
(child weight, kg), `gw_cMUAC` and `gw_muac` (mid-upper-arm circumference, cm),
and `gw_an2a` and `gw_an2b` (the raw anthropometry pair). Those are correctly
retained. A blood-derived column whose distribution resembles neither
haemoglobin nor its outcome would still pass. One column deserves a label check
in WS7b: `gw_wAnemMalarInflamm`, the top-ranked column for Sierra Leone women's
iron at |r| 0.358, whose name suggests a derived anaemia, malaria or
inflammation status rather than an exposure.

---

## 2. A defect in the detector itself, found before it could mislead

The first run of the leakage report ranked `gw_cFormMilkFreq` at |r| **0.8165**
on the smoke cell, above the 0.80 of `gw_cAnemiaYN`, the worst leak Section 7
found. It is not a leak. The column has **four** complete pairs with the
outcome. In that one cell, 134 of 2,253 guarded predictors carry fewer than 50
complete pairs.

A detector whose top rank is occupied by four points of noise is worse than no
detector, because the real leak sits below the noise and the reader learns to
skip the report. The report now carries an effective-n floor: a column below
`min_n = 100` complete pairs is ranked and shown but never flagged, and
`n_pairs` appears on every row. Across the full run, 370 column-cell
combinations fall below the floor
(source: `results/tables/leakage_report_summary.csv`, column
`n_under_measured`, summed).

Without that floor the eleventh instance would have been the second entry on a
list headed by an artefact.

---

## 3. The join lint

**Design change from the plan, forced by measurement.** The prompt asked for a
test that fails on any `by = "Admin2"`. There are 57 bare name-only joins and 12
composite joins lacking Admin1 across `R/`, `scripts/` and `dashboard/`
(source: `tests/testthat/admin2_join_baseline.csv`, column `kind`). A test that
fails on all 69 fails on first contact and is switched off within a day, which
leaves the codebase less protected than before. The lint is therefore a
**ratchet**: it fails on any site absent from a recorded baseline. Existing sites
are grandfathered and visible; a new one cannot appear without either a fix or a
deliberate, reviewable baseline amendment.

The baseline is keyed on file and trimmed source line rather than line number,
so edits elsewhere in a file do not invalidate it.

### Joins the lint caught

**In its own first version, two classification errors**, both found by the
self-test that runs the scanner over a fixture of known cases:

1. `by = c("Admin1", "Admin2")`, the *correct* pair-key form, matched the
   positional-join pattern, because the character class before the `Admin2`
   token spans the `Admin1` argument. Five sites were being reported as defects
   for using the very key the project migrated to. Fixed by exempting any line
   naming Admin1.
2. `by = c("country", "Admin2")` and `by = c("country_label", "Admin2")` were
   labelled "positional" when they are composite. They now carry their own kind,
   `composite_without_admin1`, and are **not** exempted. `country` separates
   Malawi from Ghana, and the Malawi fan is *within* Malawi, where six district
   names each occur in more than one region. Twelve sites are in this class,
   including `R/corrected/p2_p6_methods.R` and three dashboard modules.

A lint that reports the correct pattern as a defect is a lint that gets
disabled, so both corrections matter more than their size suggests.

**No new defect was found in the codebase.** All 69 sites predate this branch.
The highest-risk one for this work, `R/corrected/p9_protocol_reconciliation.R`
with six sites, produces the protocol table WS3 depends on; the unit-count check
below shows it did not fan.

---

## 4. Unit counts

The reference is derived from the built store rather than declared, so the test
tracks the data instead of agreeing with whatever the author believed:

| Country | Analytic Admin-2 units | Outcomes |
|:---|---:|---:|
| Gambia | 30 | 4 |
| Ghana | 75 | 6 |
| Sierra Leone | 14 | 6 |
| Malawi | 87 | 8 |

(source: `results/tables/admin2_unit_reference.csv`, columns `analytic_units`
and `n_outcomes`.) These match the counts the prompt gave and the counts in
`results/tables/frozen_2026-09/corrected/admin2_unit_counts.csv`, which records
Malawi's GADM level 2 at 256 polygons under 243 distinct names.

The assertion is one-sided. A table may hold fewer units than the survey
supports, because cells drop districts with no outcome data and aggregation
drops units below a minimum sample size. It may never hold more, and every
mechanism that raises the count is a defect.

**Result: 0 over-counts in 24 checks** across six tables
(source: `results/tables/admin2_unit_check.csv`, column `status`). Two rows sit
legitimately below reference: `individual_anchor.csv` reports 72 Ghana and 83
Malawi districts against 75 and 87, which is the `n >= 5` aggregation filter in
`scripts/covariates/18_individual_anchor.R`.

---

## Acceptance

| Criterion | Result |
|:---|:---|
| All three tests present and passing on the current data | 28 tests pass, 0 fail, 0 skip (`Rscript tests/testthat.R`) |
| Leakage report runs as a target and its output is a file | `leakage_report` is in the manifest with 24 incoming edges, one per `outcome_data_*` target |
| FINDINGS documents any join the lint caught | Section 3 above |

---

## Reviewer pass, statistical

**Leakage.** This is the workstream whose subject is leakage, and it found a
live instance. The measurement is a marginal Pearson correlation between each
predictor and the binary outcome at the individual level. That is the same
statistic Section 7.2 used, chosen for comparability with the recorded history
rather than for power. It will not detect a column that leaks only in
combination with another, nor a non-monotone relationship. The report is a
detector, not a proof of cleanliness, and the file says so.

**Selection on evaluation folds.** The report is computed on all rows of a cell,
with no folds involved. It is a data-screening artefact and is never used to
select predictors: nothing downstream reads it to decide what enters a model. If
it were used that way it would become an all-data prescreen of exactly the kind
Section 2.3 measures at +0.182 optimism. Its only sanctioned uses are the
testthat assertion and human inspection.

**Survey weights.** The correlations are unweighted. For a leak detector this is
the right choice: leakage is a property of the recorded values, not of the
population they represent, and weighting would add variance without adding
detection.

**Fold construction and seed coverage.** No model is fitted and no stochastic
call is made in any of the three guards, so no seed is required. The one
randomised element in the suite is the lint self-test fixture, which is a fixed
literal.

**Identical cells and identical folds.** The before-and-after comparison of the
guard fix is between two runs of the same code over the same 24 cells and the
same rows, differing only in whether two columns are present. The 0.5000 to
0.3881 change is therefore attributable to the columns and nothing else.

**A limit on the effective-n floor.** Setting `min_n = 100` means a genuinely
leaked column measured on 80 respondents would be reported but not flagged, and
the test would pass. The floor was set from the measured failure at n = 4 rather
than from a power calculation. The 370 under-measured column-cells are carried
in the report so an inspection can reach them.

## Reviewer pass, reproducibility

**Targets graph.** `tar_manifest()` lists 845 targets and includes
`leakage_report`. `tar_network(targets_only = TRUE)` shows **24 edges into it,
all from `outcome_data_*`**, one per country-outcome cell across all four
countries. That is the property the acceptance criterion asks for: a newly
ingested country adds an `outcome_data_*` target, which becomes an edge, which
rebuilds the report over its columns.

**Store staleness, measured rather than assumed.** Changing
`is_biomarker_column()` alters a function that every `outcome_data_*` target
calls, so the invalidation question had to be answered rather than reasoned
about. Under the store's own build settings, `PIPELINE_MODE=full` and
`GEE_COVARIATE_HYGIENE=true`:

| State | Outdated targets |
|:---|---:|
| Before the guard change | 806 of 845 |
| After the guard change | 806 of 845 |

**The guard change invalidates nothing that was not already invalid.** All 24
`outcome_data_*` targets were already outdated on this branch before this
workstream began. This is a known condition rather than a new one: the store
predates the assay-guard fix in commit `fe30fa8`, which is precisely why
`scripts/covariates/18_individual_anchor.R` applies the guard at read time and
says so in its comments. The same mechanism carries the WS7a fix, so WS3 can
rescore Section 5 without a store rebuild.

A first measurement of 809 was taken under the default `PIPELINE_MODE=fast` and
is not comparable, because the mode enters `pipeline_params` and invalidates
almost everything by itself. The table above uses the store's own settings for
both arms.

**Stamp targets.** No untracked file input is consumed. The report reads
`outcome_data_*` targets; the lint reads tracked source; the unit-count
reference reads `svy_admin2_*` targets.

**Paths and determinism.** All three guards use `here()`, contain no absolute
path and no `setwd()`. Re-running any of them reproduces its artefact.

**Smoke profile.** `PROFILE=smoke` restricts the leakage report to Ghana
`child_iron` and writes `_SMOKE`-suffixed artefacts, following the convention in
`scripts/run_subsample.R`. Verified: the smoke run completes and was in fact the
run that first surfaced the `gw_cFormMilkFreq` artefact that forced the
effective-n floor. Full run over 24 cells and both predictor sets completes in
well under a minute.

**Test suite entry point.** `Rscript tests/testthat.R`. The helper sources only
`R/lint_admin2_joins.R` and `R/unit_counts.R` and extracts
`is_biomarker_column()` from `R/data_prep.R` by parsing, rather than calling
`targets::tar_source("R")`, which pulls in the modelling stack and takes long
enough to stop anyone running the tests. Tests that need the guard artefacts
skip with an instruction when those are absent, so the suite runs on a clean
checkout without a built store.
