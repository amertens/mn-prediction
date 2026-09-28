# `explore/` findings log

Running record of the hypothesis-generating probes in this folder. One entry
per probe: the question, the design, the result, the honest number next to the
loose one, and a verdict.

**Verdicts.** *candidate* = worth pre-registering for a country these models
have never seen · *dead end* = tried, does not work, recorded so nobody
re-tries it · *needs data* = the idea is untested because the data is not here.

**Everything here is post hoc on the same four countries** unless an entry says
otherwise. See `README.md` for the isolation contract and the reference numbers.

---

## GATE (2026-09-28) — the harness reproduces the record

**Question.** Can a generalised scorer, which accepts an arbitrary arm and an
arbitrary feature matrix, reproduce the numbers in
`results/tables/protocol_v2/benchmarks_v2_cells.csv`?

**Design.** `explore/R/harness.R` reuses the project's own `exp_cell`
construction, `prep_predictors_v2` (fix 3), `domain_representation_v2` (fix 4),
`make_folds_v2` seeds (fix 1) and `score_v2`. Run the protocol's own arms
(`null_train_mean`, `spatial`, `domain_index`) over all cells, both targets,
all three estimands, and compare cell by cell against the record.

**Result. CLEAN.** 98.8% of 376 cell × arm comparisons within 0.02; mean
|diff| 0.0009. In-fill and region reproduce to three decimals
(level 0.394 / 0.457 / −0.199; prev 0.298 / 0.307 / −0.195). LOCO transport
reproduces **exactly** (level 0.287, prev 0.192, max |diff| 0.0000). Residual
differences are confined to the constant-prediction null arm under the region
estimand on prevalence, where Spearman is tie-dominated.

**Two traps it caught**, both of which would have silently corrupted every
probe downstream:

1. **Predictor tiers.** The headline record runs
   `V2_PREDICTOR_TIERS = "open,survey_public"` (no DHS) for *all three*
   estimands, not just transport. A probe on the default (all tiers) is
   computed on a bigger predictor set than the record it is compared to.
   `exp_load()` now defaults to the headline set.
2. **Where the domain axes are built under LOCO.** Building them per country
   and intersecting the columns lets each country learn its own PC1 sign, so a
   domain score means the opposite thing in two countries. The first harness
   did this and lost up to **0.76 Spearman** on a single cell (Gambia
   women_iron prev: −0.236 against the record's +0.524). The axes must be
   rebuilt inside each fold on the pooled matrix with
   `sign_rows = training countries`, as `02b_merge_and_loco.R` does.

**Verdict.** Gate passed; probes may proceed. Trap 2 is worth carrying into
any future cross-country work in the main pipeline too — it is a silent,
sign-flipping failure that looks like a null result rather than a bug.

---
