# NX-01. Is each transported map specific to its nutrient?

28 September 2026. Script `scripts/protocol_v2/76_nutrient_specificity_swap.R`;
tables `results/tables/protocol_v2/nx01_per_source.csv`, `nx01_per_cell.csv`,
`nx01_summary.csv`.

## Question

The transported maps for different nutrients agree closely (about 0.93 between
outcomes in Côte d'Ivoire). The talk's "which deficiencies can it map" table
could therefore be describing one general deprivation map instead of separate
nutrient maps. This is a weight-swap negative control.

## Design (fixed before the run)

Leave-one-country-out, as in 02b estimand C. The level target is primary and
the tiers are open plus survey_public. For each of the 22 cells (held-out
country H, outcome B), B's districts in H are scored with:

- **own**: the index trained on B in the other countries. This is the committed benchmark.
- **swap**: for each outcome of a different nutrient measured in at least 2
  countries other than H, the index trained on that outcome outside H and
  applied to H.
- **same nutrient, other population** (for example child iron for women's
  iron): reported separately, not counted as a swap.
- **generic**: the sum of every eligible outcome's weight vector on one shared
  basis, a single "deficiency in general" index.

Zinc is measured only in Malawi, and selenium and iodine are not in
`targets_v2.csv`, so no source for those nutrients was ever eligible.

**Rule.** A cell is nutrient-specific if own minus the best swap is at least
0.03 (level Spearman). The nutrient table is defended for a nutrient if a
strict majority of its cells are specific. Own minus the mean swap is reported
as the second reading, because a maximum over noisy swaps is biased upward.

**Reproduction.** Own matches `benchmarks_v2_cells.csv` (country,
domain_index) in all 22 cells to 4e-16, on both the level and prevalence
targets.

## Result (level target)

| Nutrient | Cell | Own | Best swap (source) | Mean swap | Generic | Same nutr., other pop. | Specific |
|---|---|---|---|---|---|---|---|
| B12 | Ghana | 0.60 | 0.54 (women vitA) | 0.27 | 0.51 | | yes |
| B12 | Malawi | 0.44 | 0.31 (child vitA) | 0.17 | 0.27 | | yes |
| B12 | Sierra Leone | 0.56 | 0.60 (child vitA) | 0.34 | 0.53 | | no |
| Folate | Ghana | -0.07 | 0.32 (women iron) | 0.27 | 0.30 | | no |
| Folate | Malawi | 0.08 | 0.12 (women iron) | -0.06 | -0.06 | | no |
| Folate | Sierra Leone | 0.19 | -0.07 (women iron) | -0.18 | -0.10 | | yes |
| Iron, child | Gambia | 0.37 | 0.50 (women vitA) | 0.43 | 0.45 | 0.26 | no |
| Iron, child | Ghana | 0.54 | 0.58 (child vitA) | 0.32 | 0.56 | 0.57 | no |
| Iron, child | Malawi | -0.08 | 0.02 (women vitA) | -0.09 | -0.06 | 0.11 | no |
| Iron, child | Sierra Leone | 0.24 | 0.28 (women B12) | 0.16 | 0.20 | 0.11 | no |
| Iron, women | Gambia | 0.60 | 0.71 (women vitA) | 0.64 | 0.64 | 0.59 | no |
| Iron, women | Ghana | 0.36 | 0.31 (child vitA) | 0.17 | 0.27 | 0.37 | yes |
| Iron, women | Malawi | 0.22 | 0.12 (women vitA) | 0.04 | 0.08 | 0.07 | yes |
| Iron, women | Sierra Leone | -0.12 | 0.24 (women B12) | 0.06 | 0.09 | 0.01 | no |
| VitA, child | Gambia | 0.62 | 0.65 (women B12) | 0.57 | 0.64 | 0.67 | no |
| VitA, child | Ghana | 0.36 | 0.37 (child iron) | 0.21 | 0.35 | 0.32 | no |
| VitA, child | Malawi | 0.26 | 0.27 (child iron) | 0.20 | 0.25 | 0.27 | no |
| VitA, child | Sierra Leone | 0.07 | 0.38 (women iron) | 0.21 | 0.14 | 0.30 | no |
| VitA, women | Gambia | 0.75 | 0.74 (women B12) | 0.67 | 0.73 | 0.73 | no |
| VitA, women | Ghana | 0.31 | 0.48 (child iron) | 0.25 | 0.40 | 0.43 | no |
| VitA, women | Malawi | 0.23 | 0.26 (women iron) | 0.18 | 0.23 | 0.24 | no |
| VitA, women | Sierra Leone | 0.02 | 0.20 (women B12) | 0.06 | -0.05 | -0.02 | no |
| **All 22** | | **0.30** | **0.36** | **0.22** | **0.29** | **0.31** (16 cells) | **5** |

## Verdict against the pre-registered rule

| Nutrient | Specific (best-swap test) | Verdict | Second reading (mean swap) |
|---|---|---|---|
| B12 | 2 of 3 | **defended** | 3 of 3 |
| Iron | 2 of 8 | not defended | 4 of 8 (not a majority) |
| Vitamin A | 0 of 8 | not defended | 6 of 8 |
| Folate | 1 of 3 | not defended | 2 of 3 |
| All | 5 of 22 | | 15 of 22 |

**B12 is a nutrient map.** Its own index beats every other nutrient's index
in Ghana and Malawi. It beats the generic index in all three cells, and its
maps agree least with the swaps (mean Spearman 0.52 between own and swap
maps).

**Iron and vitamin A are one map.** The vitamin A indices rank iron districts
as well as iron's own index does (0.27 and 0.28 against 0.27). The iron
indices rank vitamin A districts slightly better than vitamin A's own index
(0.39 and 0.34 against 0.33). The own and swap maps agree at 0.68 (iron) and
0.74 (vitamin A). For these two nutrients the honest wording is: one general
deprivation map, which ranks this nutrient about as well as any other.

**Folate** is not defended. Its own index barely transports (0.07), so the
question does not really arise.

**Generic.** Across all 22 cells the generic index equals own (0.29 against
0.30; own is ahead in 12 of 22). The borrowed same-nutrient index from the
other population is slightly ahead of own (0.31 against 0.30, ahead in 9 of
16). This agrees with XO-01.

**The two readings disagree** for vitamin A and folate, where the verdict flips
under the mean reading. The registered verdict is the best-swap reading. The
mean reading is not clean either. Every mean includes the folate index, which
is a poor map for everything (mean Spearman at most 0.14 as a swap), so part
of what the mean reading measures is "not the folate map" rather than
specificity. For iron and vitamin A cells, the best swap was the other of the
two nutrients in 11 of 16 cells and B12 in the other 5. The iron and vitamin A
swaps were trained on the same three countries as own.

## Secondary: prevalence target (not used for the verdict)

5 of 22 cells are specific (iron 2 of 8, vitamin A 0 of 8, folate 2 of 3, B12
1 of 3). B12 loses its majority on prevalence: in Malawi the child iron
index ranks B12 prevalence better (0.36 against 0.22). The B12 result
therefore rests on the level target and on three cells.

## Caveats

- Three cells per nutrient for B12 and folate. Sierra Leone has 14 districts,
  where a Spearman correlation has a standard error of about 0.27. Dropping
  Sierra Leone gives 4 of 16 specific, with the same verdict for every
  nutrient.
- The best-swap test is strict by design: the maximum of 4-5 noisy
  correlations is biased upward.
- The number of training countries differs. Folate and B12 exist in three
  countries, so for B12 cells the iron and vitamin A swaps trained on three
  countries against own's two. That favours the swaps, so the B12 result is
  conservative.
- The generic index needs one shared basis, learned on the stacked training
  rows of all six outcomes. Its per-outcome weights therefore differ slightly
  from each outcome's own fit.

## For a policy audience

Except for vitamin B12, the maps find districts that are deprived in general
rather than short of one particular nutrient, so the iron map and the vitamin A
map point to largely the same places.
