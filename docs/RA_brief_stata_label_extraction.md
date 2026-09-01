# RA brief — recover Stata variable labels and re-derive the assay guard

**Status:** ready to start. **Estimated effort:** 2–4 days.
**Deliverable:** a label dictionary, plus a labels-versus-names disagreement
report for `is_biomarker_column()`.

## Why this matters (read this first)

`is_biomarker_column()` in `R/data_prep.R` decides which survey columns are
"blood draw derived" and must be excluded from any predictor set. If it lets one
outcome-derived column through, a model predicts the outcome from itself and the
project reports a spurious success.

**It has already failed twice.** The individual-level anchor analysis returned
r = 0.973, then r = 0.986, before the leak was found; the honest figure is
0.228. The columns that caused it were not named after any analyte —
`gw_cAnemiaYN` (|r| 0.80), `gw_cAnemiaCat` (0.78), `gw_bs2` (0.69), `gw_bis`
(0.67), `gw_cID_NoAdj` (0.57). No pattern over analyte names could have caught
them.

The guard is now a five-class regular expression over **variable names**, and
the project's own review states the decisive limitation plainly: the merged data
carries **no variable labels at all** — `attr(x, "label")` is `NULL` for every
column tested — so columns such as `gw_bis`, `gw_rpb1` and `gw_bs2` were
classified as blood-sample fields *by name pattern and correlation*, not by
documentation. The review's own recommended fix is exactly this task: obtain the
original Stata files with labels intact and re-derive the guard from labels
rather than names.

Two failure modes are currently undetectable:

- **Under-blocking.** A leaked column with a merely moderate correlation would
  not stand out in the correlation ranking that found the others.
- **Over-blocking.** `gw_wFFAnemia` is *assumed* to be a "beliefs about
  fortified food" item because of its `gw_wFF*` siblings. If it is actually a
  measured status it should be blocked. Conversely, blocking every column
  containing `NoAdj` may be discarding legitimate unadjusted exposures.

A fifth country (Tanzania) is planned and will arrive with its own naming
conventions, so a name-based guard has to be re-verified from scratch each time.
A label-based guard does not.

## Tasks

### 1. Locate the original files

Find the original Stata `.dta` (or SPSS `.sav`) files for all four surveys:
Gambia 2021, Ghana 2017, Malawi 2015, Sierra Leone 2013. Start from
`src/<Country>/` — those scripts are the provenance for `data/IPD/*` and will
name the files they read. Record for each survey the exact file, its version or
date, and where it came from, in a provenance table. If a file cannot be found,
say so explicitly rather than substituting a different one.

### 2. Extract the label dictionary

For every column in every survey file, extract:

| Field | Notes |
|:---|:---|
| `survey` | country and year |
| `variable` | column name as it appears in the raw file |
| `label` | the Stata variable label |
| `value_labels` | the value-label set, if any, as `code=label` pairs |
| `type` | numeric / string / labelled |
| `n_nonmissing` | count |

Use `haven::read_dta()` and read the `label` attribute of each column; do not
retype labels by hand. Write to `data/labels/<country>_variable_labels.csv` and
a combined `data/labels/all_variable_labels.csv`.

Note that many columns are renamed between the raw file and the merged dataset.
Where the merge scripts rename, carry **both** names so the dictionary can be
joined to the analysis data.

### 3. Classify from labels, independently

Working **from the labels only, without looking at the current guard**,
classify every column in each survey's `Xvars_full` into:

- `blood_draw_derived` — measured from the blood sample, or derived from it
  (this includes derived status flags such as anaemia category, and unadjusted
  twins of adjusted analytes)
- `blood_sample_admin` — fieldwork bookkeeping about the draw itself (sample
  IDs, tube numbers, whether a sample was taken)
- `not_blood` — everything else
- `unclear` — the label is absent or too terse to decide

Doing this blind matters. If you read the regex first you will tend to
reproduce it, and the whole value of the task is an independent opinion.

### 4. Produce the disagreement report

Join your classification to `is_biomarker_column()`'s verdict on the same
columns and report every disagreement, in two tables:

- **Guard says keep, labels say blood-derived.** These are candidate live
  leaks. Rank by absolute correlation with each outcome so the highest-risk
  ones are obvious.
- **Guard says block, labels say not blood.** These are over-blocks costing
  legitimate predictors. Check `gw_wFFAnemia`, `gw_wHeardAnemia` and the
  `NoAdj` family specifically — the review names all three as assumed rather
  than verified.

For every disagreement give the variable, its label, both verdicts, and a
one-line recommendation.

### 5. Propose a label-based guard

Write the classification rules you would use if labels were available for every
country, and note which of the current regex's five classes become unnecessary.
Do **not** edit `R/data_prep.R` — propose, and let the analyst make the change,
because altering the guard invalidates cached results and changes published
numbers.

## Acceptance criteria

- Label dictionaries for all four surveys, with a provenance table naming each
  source file.
- An independent, blind classification of every `Xvars_full` column.
- A disagreement report with every case listed, not a summary count.
- An explicit statement of coverage: what fraction of columns have a usable
  label, per survey. If a survey's file has no labels either, that is a
  legitimate and important negative result — report it rather than inferring
  labels from names.

## Why blind classification, restated

The point of this task is to produce evidence the current guard cannot produce
about itself. A regex over opaque names cannot be validated by more careful
reading of the same names. Only an independent source of truth — the labels the
survey team wrote — can confirm or refute it, and only if it is consulted
independently.
