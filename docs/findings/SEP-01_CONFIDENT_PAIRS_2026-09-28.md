# SEP-01: accuracy on the pairs the survey can tell apart (28 Sep 2026)

**Question.** The talk's "pairs" score counts every pair of districts, including
pairs whose survey difference sits inside the survey's own sampling noise. How
well do the methods do on the pairs the survey clearly separates?

**Answer.** When the survey clearly separates two districts, the model puts them
in the survey's order 83 times in 100 (regional figures 80, neighbour map 83).
Only about one pair in six qualifies. Every method gains 16 to 18 points on
these pairs, so the rise mostly shows that they are easier pairs. It is not
evidence that the model does especially well on them.

## Design (fixed before running; see the script header)

- District standard error. For the level target: the respondent-level SD
  (`sd_level`, which is an SD and not an SE) divided by sqrt(`n_eff_cont`).
  For prevalence: sqrt(p(1-p)/`n_eff`). The 8 one-respondent districts get an
  infinite SE.
- A pair is confident when |difference| > 1.96 x sqrt(SE_i^2 + SE_j^2).
- Arms are the same held-out predictions as the deck, averaged over the ten
  in-fill draws. The model and the neighbour map skip tied predictions. The
  fair regional figure (script 39) counts same-region pairs as half.
- The probability-weighted score weights each pair by 2 Phi(|d|/SE) - 1.
- The level target is primary. The main cell set is the deck's 14 measurable
  in-country combinations.

**Reproduction.** All 104 checks pass, with a largest difference of 7e-16:
- all-pairs model scores match `v3_percell_pairs.csv` (18 of 18);
- fair regional scores match `v6_fair_regional_pairs.csv` (14 of 14);
- Spearman for the model and the neighbour map matches `benchmarks_v2_cells.csv`
  for both targets.

## Results: level, 14 combinations (mean over combinations)

| | model | regional (fair) | neighbour map | coin |
|---|---|---|---|---|
| All pairs (the deck) | 66.1 | 64.6 | 65.8 | 50 |
| Confident pairs only | 83.4 | 80.1 | 82.8 | 50 |
| Probability-weighted | 71.2 | 69.3 | 70.8 | 50 |

- 16.6% of pairs are confident. By combination this runs from 3% (Malawi
  children's iron, 105 pairs) to 35% (Gambia women's iron).
- On confident pairs the model beats the regional figure in 10 of 14
  combinations and the neighbour map in 7 of 14.
- The model's biggest loss is Ghana children's vitamin A: 74 against 93 for the
  regional figure.
- Pooled over all confident pairs, rather than averaged over combinations: model
  82, regional 80, neighbour map 84.
- Prevalence (secondary) on confident pairs: model 80, regional 74, neighbour
  map 81. 11.8% of pairs are confident.
- All 18 combinations, level, on confident pairs: model 79, regional 77,
  neighbour map 78.

## Caveats

- Confident pairs are a selected, easier subset for every method. Their median
  survey gap is 2.5 times that of all pairs. They are also half as likely to
  fall within one region (3.9% vs 8.4%), so the regional figure loses fewer
  half-pairs on them.
- The model's lead over the regional figure grows from 1.5 to 3.3 points. On a
  10-of-14 split that is weak evidence, and the neighbour map is level with the
  model.
- The prevalence SE is zero for districts at 0% or 100%, and many districts sit
  there. That overstates confidence for their pairs.
- No uncertainty interval was pre-specified or computed.

Files: `scripts/protocol_v2/77_confident_pairs.R`,
`results/tables/protocol_v2/sep01_{cells,summary,districts,reproduction}.csv`.
Only the summary's NA handling was changed after the first run: two prevalence
cells have no confident pair, and neither is among the 14. The per-cell results
were identical across both runs.
