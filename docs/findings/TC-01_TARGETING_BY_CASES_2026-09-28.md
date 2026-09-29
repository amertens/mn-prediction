# TC-01: targeting by predicted cases, under two budget rules (28 September 2026)

Script: `scripts/protocol_v2/74_targeting_by_cases.R` (design pre-registered in its header before any result).
Tables: `results/tables/protocol_v2/tc01_targeting_by_cases{,_cells,_summary}.csv`.
Figure: `results/figures/mnf15_v6/v6_targeting_rules.png` (script `scripts/policy_deck/43_mnf15_v6_targeting_and_transport.R`).

## Why

Script 12 ranks districts by predicted rate and counts the deficient people in the worst-ranked fifth. In a country held out, that fifth is more deficient than the country but holds only 17% of its deficient people, less than a random fifth, because the high-rate districts are small. The question is which ranking a programme should use, given how its budget is set.

## Design

Script 12's units, cells, folds and rate ranking, reproduced exactly: in-fill capture matches in 180 of 180 draws, transport in 22 of 22 (largest difference 5e-16). Measurable combinations only: 14 in-country (Sierra Leone has no in-country test) and 16 with the country held out.

- **Budget of districts** (choose a fifth of them): share of deficient people inside the chosen districts. Compared: random, the model's rate, population alone (no model), the model's expected cases (rate x population), perfect knowledge.
  - The in-country rate comes from the calibrated index.
  - For a held-out country, the rate is anchored at the country's national prevalence, with its spread shrunk by the training countries as in AR-01.
- **Budget of people** (reach a fifth of the target population, the last district in part): share of deficient people reached. Compared: random (20%), the survey's regional figures (in-country only), the model's rate, perfect knowledge.
- **Pre-registered readings:**
  - Budget of districts: expected cases must beat population alone on the mean and in a majority of cells.
  - Budget of people: the model's rate must beat 20% on the mean and in a majority of cells.

## Results (mean over combinations, share of deficient people reached)

| Rule | Inside a surveyed country (14) | Country with no survey (16) |
|---|---|---|
| **Budget: a fifth of the districts** | | |
| Random | 20% | 20% |
| Model, highest rates | 25% | 18% |
| Most populous districts (no model) | 48% | 46% |
| Model, most expected cases | 49% (better than population in 8 of 14: PASS, a tie in practice) | 44% (2 of 16: FAIL) |
| Perfect knowledge | 64% | 61% |
| **Budget: a fifth of the people** | | |
| Random | 20% | 20% |
| Survey's regional figures, highest rates | 27% | n/a |
| Model, highest rates | 31% (above 20% in 13 of 14: PASS) | 25% (9 of 16: PASS, weakly) |
| Perfect knowledge | 49% | 47% |

## Reading

- **When the budget is a number of districts,** choose the most populous ones. Population alone reaches 46 to 48% of deficient people, and the model adds nothing to it.
  - In a new country, the shrunk rate is often flat: the training countries' own cross-border ranking was at or below zero for B12 and several iron combinations.
  - So expected cases there are close to population alone.
  - Ranking districts by rate for this budget is the wrong rule (18 to 25%).
- **When the budget is per person** (supplements, a fixed number of people to screen), rank by the model's rate.
  - Inside a surveyed country it reaches 31% of deficient people with a fifth of the population, against 27% for the survey's regional figures and 20% at random.
  - In a country with no survey it reaches 25%.
- **Talk wording:** "The worst fifth of districts holds 25% of deficient people" answers a question no programme asks. Use the per-person budget instead: with a fifth of the people covered, the model reaches 31 in 100 deficient people inside a surveyed country, 25 with no survey, against 20 at random.

## Caveats

- Capture is computed among surveyed districts only, where the true prevalence is known.
- The regional figure is fold-based, with the handicap described in script 39. Multiplying by population dilutes it for the district budget, but not for the people budget.
- The rho_train shrinkage follows AR-01 as pre-registered. A less-shrunk spread was not tested.
