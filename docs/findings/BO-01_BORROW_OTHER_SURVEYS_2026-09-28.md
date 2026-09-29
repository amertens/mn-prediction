# BO-01: borrowing the other surveys when filling in a surveyed country (28 Sep 2026)

**Question.** The in-fill index learns its weights from the 24 to 70 training
districts of one survey. The same outcome was measured in up to three other
surveys. Does adding their evidence to the weights improve the in-country
ranking?

**Answer. No. The pre-registered test FAILS.** Borrowing lowers the mean in-fill
Spearman from 0.399 to 0.377 (gain -0.023). It is better in 8 of 18 cells,
worse in 8, and tied in 2. The bar was a gain of at least +0.03 and wins in at
least 12 of 18. The secondary settings point the same way.

## Design (fixed before running; see the script header)

- `index_borrow`: z_total = z_own + sum over donors of z_C. Here z_own is
  `.index_weights_v2` on the target country's training districts, and z_C is
  the same function on all of a donor country's districts. Donors are every
  other country that measured the outcome and has at least 12 districts. The
  prediction is rescaled exactly as `arm_domain_index_v2` does it (rho = 1).
- Components have to mean the same thing in every country, and 02b is how the
  project arranges that. It rank-normalises within each country, keeps the
  columns common to all the countries, stacks them, and fits one PC basis per
  domain on the training rows (`sign_rows` = the target's training districts
  plus every donor district). BO-01 copies that construction.
  Per-country bases cannot be summed. Across countries, same-named PC1s
  have a median loading cosine of 0.70 to 0.74, and 22% point the opposite
  way. PC2 and later are unrelated (median cosine 0.03, 45% opposite)
  (`bo01_borrow_alignment.csv`).
- `domain_index_common`: the country's own weights on that same pooled basis.
  This arm separates the change of representation from the borrowing itself.
- Malawi zinc (both populations) has no donor, so it is a tie by construction.
- Primary: level target, in-fill 5-fold x 10 draws (`make_folds_v2` rep_id
  1..10), and the 18 cells with a finite benchmark Spearman.

**Checks.**
- `domain_index` reproduces `benchmarks_v2_cells.csv` in 96 of 96
  cell-settings, with a largest difference of 0.0000. That covers in-fill and
  region, level and prevalence.
- 748 leakage checks passed. Replacing the target's held-out outcomes with
  arbitrary values never moved a prediction.

## Primary result: level, in-fill (mean Spearman over 10 draws)

| Country | Outcome | Donors | Own index | Own, pooled basis | Borrow | Gain |
|---|---|---|---|---|---|---|
| Gambia | child vit A | GH, MW, SL | 0.706 | 0.713 | 0.686 | -0.021 |
| Gambia | women vit A | GH, MW, SL | 0.697 | 0.687 | 0.726 | +0.028 |
| Gambia | child iron | GH, MW, SL | 0.369 | 0.386 | 0.379 | +0.010 |
| Gambia | women iron | GH, MW, SL | 0.668 | 0.675 | 0.693 | +0.025 |
| Ghana | child vit A | GM, MW, SL | 0.277 | 0.287 | 0.311 | +0.034 |
| Ghana | women vit A | GM, MW, SL | 0.464 | 0.460 | 0.367 | -0.098 |
| Ghana | child iron | GM, MW, SL | 0.529 | 0.523 | 0.539 | +0.010 |
| Ghana | women iron | GM, MW, SL | 0.403 | 0.395 | 0.334 | -0.068 |
| Ghana | women folate | MW, SL | 0.385 | 0.383 | 0.442 | +0.058 |
| Ghana | women B12 | MW, SL | 0.567 | 0.574 | 0.582 | +0.015 |
| Malawi | child vit A | GM, GH, SL | 0.230 | 0.253 | 0.235 | +0.004 |
| Malawi | women vit A | GM, GH, SL | 0.302 | 0.313 | 0.254 | -0.048 |
| Malawi | child iron | GM, GH, SL | 0.352 | 0.345 | 0.115 | -0.237 |
| Malawi | women iron | GM, GH, SL | 0.244 | 0.226 | 0.214 | -0.029 |
| Malawi | women folate | GH, SL | 0.412 | 0.400 | 0.359 | -0.053 |
| Malawi | women B12 | GH, SL | 0.694 | 0.692 | 0.658 | -0.036 |
| Malawi | child zinc | none | -0.105 | -0.105 | -0.105 | 0 |
| Malawi | women zinc | none | -0.007 | -0.007 | -0.007 | 0 |
| **Mean** | | | **0.399** | **0.400** | **0.377** | **-0.023** |

**Verdict: FAIL.** The gain was -0.023 against a bar of +0.03, with 8 wins
against a bar of 12.

The loss comes from borrowing, not from the representation:

- Columns: the own index uses 351 to 367 predictor columns (66 to 87
  components). The common set keeps 329 to 331 (86 to 92 pooled components).
- Representation: the own weights on the common, pooled basis score 0.400
  against 0.399.
- Borrowing: on that same basis, borrowing costs 0.023 and is better in only
  6 of the 16 cells that had donors.
- Without Malawi child iron (-0.237), the mean gain is still -0.010.

## Secondary (not used for the verdict)

| Setting | Own index | Borrow | Gain | Better in |
|---|---|---|---|---|
| Prevalence, in-fill | 0.301 | 0.292 | -0.009 | 7 of 18 |
| Level, region hold-out | 0.385 | 0.366 | -0.019 | 7 of 18 |
| Prevalence, region hold-out | 0.283 | 0.279 | -0.004 | 8 of 18 |
| Pairs in the survey's order, level in-fill | 63.8% | 63.1% | -0.7 pt | 7 of 18 |

By country (level, in-fill):

| Country | Gain | Better in |
|---|---|---|
| The Gambia | +0.011 | 3 of 4 |
| Ghana | -0.008 | 4 of 6 |
| Malawi | -0.050 | 1 of 6 with donors |

The prior expectation was half right. The Gambia is the only country that
gains in-fill (+0.011 level, +0.031 prevalence). Under region hold-out,
however, The Gambia loses (-0.030 level, 0 of 4), and the overall region gain
is negative.

## Why it does not help

- **The countries' weights disagree.** Across components, the correlation
  between a country's own weights and the donors' summed weights averages 0.19
  (range -0.07 for Malawi child iron to 0.32).
- **The donors outweigh the home survey.** The fixed-effect sum gives the
  donors 1.0 to 2.8 times the weight of the home survey. Where the two
  disagree, the home signal is diluted.
- **Agreement tracks the gain.** Across the 16 cells with donors, the gain
  rises with that agreement (Pearson 0.56, Spearman 0.37). This is
  descriptive, not tested.

## Caveats

- The fixed-effect sum assumes that every country has the same association.
  A down-weighted or random-effects version was not pre-registered and was not
  run.
- Only three countries can be scored as targets. Sierra Leone's 14 districts
  leave fewer than 12 training rows per fold, so it serves only as a donor.
- Two of the 18 cells (zinc) cannot borrow and are ties by construction.
- All of this is under the headline tiers (open + survey_public), on
  `targets_v2.csv` as rebuilt on 27 September.

**For a policy audience.** Adding other countries' survey results to a
country's own survey did not improve its district rankings, so maps for a
surveyed country should keep relying on that country's own survey.

## Files

- Script: `scripts/protocol_v2/75_borrow_other_surveys_infill.R` (the
  pre-registration is in the header).
- Tables, all in `results/tables/protocol_v2/`:
  - `bo01_borrow_cells.csv`: per cell, target and estimand; column counts,
    donors, weight diagnostics and pairs.
  - `bo01_borrow_summary.csv`: the verdict, secondary settings and per-country
    rows.
  - `bo01_borrow_raw.csv`: one row per draw and arm.
  - `bo01_borrow_reproduction.csv`
  - `bo01_borrow_alignment.csv`
- Runtime: 9 minutes, one R process.
