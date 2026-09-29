# The dietary block was never given a fair test: LOCO's 4-country intersect discards it

*2026-09-29. Follow-up to the question "are dietary intake or food security
variables important predictors in any models?" — the first answer was "no, with
one exception", and that answer was measuring the wrong thing.*

---

## The claim, corrected twice

**First statement (wrong):** three dietary domains — Market prices (RTFP),
Adult nutrition, Dietary inadequacy (MIMI) — "were never ablated", implying the
ablation table was stale.

**Second statement (also wrong):** they are excluded by the leakage filter.
They are not. All three are in `predictors_admin2_shared.csv` and all three
survive `drop_near_outcome_v2()`: Adult nutrition 3 columns, Dietary inadequacy
6, Market prices 6. Only `Nutrition status (MODELLED SURFACE)` (1 column) is
dropped as near-outcome, correctly.

**What is actually true.** `scripts/protocol_v2/23_domain_ablation_loco.R:81`
pools on

```r
common <- Reduce(intersect, lapply(cl, function(z) colnames(z$X)))
```

columns present in **all four countries**. Any predictor missing in one country
vanishes from every LOCO analysis. Re-running the ablation raised it from 21 to
25 domains — and MIMI and RTFP are *still* absent, because they are structurally
ineligible, not stale.

## What the intersect costs, and who pays

| | |
|---|---|
| predictors after the leakage filter | 533 |
| present in all 4 countries (LOCO-usable) | 482 (90%) |
| present in 3 / 2 / 1 | 33 / 10 / 8 |
| **lost to the intersect** | **51 (10%)** |

The 10% is not spread evenly. It is almost entirely the dietary and food block:

| domain | lost | of | % |
|---|---|---|---|
| Dietary inadequacy (MIMI) | 6 | 6 | **100%** |
| Market prices (RTFP) | 6 | 6 | **100%** |
| Adult and maternal mortality | 9 | 9 | 100% |
| Household diet and consumption (HCES) | 13 | 15 | **87%** |
| Food prices and supply | 8 | 12 | **67%** |
| Agricultural production | 2 | 29 | 7% |
| *Climate and weather* | **0** | 60 | **0%** |
| *Soil characteristics* | **0** | 45 | **0%** |
| *Satellite embedding* | **0** | 64 | **0%** |

The remotely-sensed domains are global rasters, so they are complete by
construction. Dietary data is survey-derived and country-specific by nature, so
it is exactly the block the intersect deletes.

## So the dietary ablation result does not mean what it looked like

The earlier table showed `Household diet and consumption (HCES)` scoring
**−0.134 alone** under LOCO and −0.008 to drop — apparently worse than useless.
That number was computed on **2 of its 15 columns**. MIMI and RTFP were
computed on **zero of theirs** and do not appear at all. Climate, soil and
satellite were computed on 100% of theirs.

This is not a fair comparison, and "diet does not predict micronutrient status"
was never a conclusion the design could support. The correct statement is:
**the project's headline estimand cannot see the dietary block.**

Per-country coverage, which is the underlying cause:

| domain | Gambia | Ghana | Malawi | Sierra Leone |
|---|---|---|---|---|
| Household diet (HCES) | 93% | **20%** | 100% | 100% |
| Dietary inadequacy (MIMI) | 0% | **100%** | 0% | 0% |
| Market prices (RTFP) | 100% | **0%** | 100% | **0%** |
| Food prices and supply | 96% | 81% | 85% | 93% |
| Infant and young child feeding | 100% | 99% | 99% | 100% |
| Adult nutrition | 100% | 100% | 100% | 100% |

MIMI exists for one country, so under leave-one-country-out it is useless by
construction whichever fold it is in.

## The highest-leverage fix is Ghana, and it is already half-built

Which country's gap, if filled, returns the most columns to LOCO:

| fill | columns restored |
|---|---|
| **Ghana** | **22** |
| Malawi | 5 |
| Sierra Leone | 4 |
| Gambia | 2 |

**12 of the 13 lost HCES columns are lost to Ghana alone.** The cause is known
and specific: the production builder reads
`data/LSMS/g7aggregates_hhlevel.dta`, which carries no item detail, so Ghana
yields 3 of 15 HCES columns. `explore/scripts/25` and `29` already produce the
missing indicators from GLSS7 section 9b (528,678 household-item rows, 13,924
households, 484 items, brand labels resolved by item-code block). They are
blocked on one artefact — the GSS 216-district list in official per-region
order — and on nothing else.

So the sequence is: **GSS district list → Ghana HCES at Admin-2 → 22 columns
restored → the dietary block becomes testable under LOCO for the first time.**

## What was learned from the domains that could be re-tested

Re-running the ablation did test four previously-absent domains (they are
4-country complete). None is load-bearing, but one is interesting:

| domain | Δ when dropped | alone (level) |
|---|---|---|
| **Adult nutrition** | 0.000 | **0.197** |
| Fertility, reproductive health | +0.002 | 0.055 |
| Healthcare access | −0.002 | 0.058 |
| Child mortality | 0.000 | −0.052 |

Adult nutrition transports at 0.197 on its own — 5th best of 25 solo scores,
above Anaemia (0.156) and Infection (0.090) — yet removing it costs nothing.
That is the signature of a domain fully redundant with the others, not a weak
one. Three columns of adult BMI/height carry most of what the larger
socio-economic blocks carry.

Also worth noting from the re-run: `Agricultural production, land use` improved
to **0.242 alone** (from 0.171) and remains the only food-related domain that is
load-bearing (Δ +0.012). But it is food *production* from remote sensing, not
intake.

## Recommendations

1. **Report LOCO results with a coverage column.** A domain's LOCO score is
   uninterpretable without the fraction of its columns that survived the
   intersect. Three of the dietary domains currently show scores computed on a
   minority or none of their columns, with nothing marking that.
2. **Add a guard** to `23_domain_ablation_loco.R` that refuses to report a
   domain whose surviving-column share is below, say, 50%, or reports it
   explicitly as "not evaluable".
3. **Prioritise the GSS district list.** It is the gate on 22 columns, and on
   the only honest test of whether diet predicts micronutrient status here.
4. **Evaluate the dietary block on the within-country estimands** (in-fill and
   region), which do not require the intersect and can use all 15 HCES columns
   in Malawi, Sierra Leone and Gambia today. That is a test that can be run
   without waiting for Ghana.
