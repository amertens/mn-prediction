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

---

## DC-01 — the within-country test, run: diet predicts, but adds nothing

*2026-09-29 · `explore/scripts/30_within_country_diet.R` · answers recommendation 4 above*

Three predictor sets per country × outcome × target, identical rows and
identical folds, scored with the estimator of record (zero-tuning domain
index): **full** (533 columns), **nodiet** (467, the six dietary domains
removed), **diet** (66, those domains alone). Agricultural production is
remotely sensed and already known load-bearing, so it stays in `nodiet` —
this measures the survey-derived dietary block specifically.

### The dietary block is a strong predictor on its own

| estimand | full (533 cols) | nodiet (467) | **diet alone (66)** |
|---|---|---|---|
| in-fill | 0.372 | 0.374 | **0.288** |
| region | 0.316 | 0.322 | **0.251** |

Sixty-six dietary columns reach **77% of what all 533 reach** in-fill, and 79%
under region extrapolation. Whatever else is true, "diet does not predict
micronutrient status" is not.

### But it adds nothing on top of everything else

| estimand | cells where full > nodiet | median Δ |
|---|---|---|
| in-fill | 14 of 36 | **−0.0030** |
| region | 15 of 36 | **−0.0027** |

A clean null. The dietary block is **redundant with the rest, not weak** — the
same signature Adult nutrition shows under LOCO (0.197 alone, 0.000 to remove).
The remotely sensed climate/soil/satellite block already carries the
information diet carries.

### The LOCO verdict was an artefact, and here is the size of it

Identical domains, scored alone:

| domain | LOCO (alone) | **within-country (alone, level)** | columns it got under LOCO |
|---|---|---|---|
| Food prices and supply | 0.011 | **0.316** | 4 of 12 |
| Infant and young child feeding | 0.042 | **0.306** | 24 of 24 |
| Dietary inadequacy (MIMI) | *not evaluable* | **0.302** | 0 of 6 |
| Household diet (HCES) | **−0.134** | **0.263** | 2 of 15 |
| Adult nutrition | 0.197 | 0.244 | 3 of 3 |
| Market prices (RTFP) | *not evaluable* | 0.155 | 0 of 6 |

Household diet goes from −0.134 to +0.263; MIMI, which LOCO literally cannot
see, scores 0.302 within Ghana on six columns. The two domains LOCO evaluated
in full (IYCF, Adult nutrition) are the two whose LOCO and within-country
scores are closest — exactly what the coverage explanation predicts.

### Gambia is the standout, and Ghana is where diet hurts

| country × target | diet adds | median Δ | diet alone |
|---|---|---|---|
| **Gambia** level (region) | **3 of 4** | **+0.0082** | **+0.632** |
| **Gambia** prev (region) | **3 of 4** | **+0.0106** | **+0.588** |
| Gambia level (in-fill) | 1 of 4 | −0.0024 | +0.602 |
| Malawi level (in-fill) | 4 of 8 | +0.0013 | +0.210 |
| Malawi prev (in-fill) | 5 of 8 | +0.0050 | +0.141 |
| **Ghana** level (region) | **0 of 6** | **−0.0212** | +0.310 |
| Ghana level (in-fill) | 0 of 6 | −0.0075 | +0.339 |

Gambia is the only country where the dietary block **reliably adds** — 3 of 4
cells on both targets under region extrapolation — and its solo score, 0.602
in-fill and 0.632 region, is far above the pooled full-model median. Gambia
also has the most dietary columns of any country (58, against Ghana's 43).

Ghana is the mirror image: diet never adds and actively hurts under region
extrapolation (−0.021). Ghana is the country whose HCES block is crippled — 12
of the 13 lost HCES columns are lost to Ghana alone. A block that is 20%
populated is worse than no block, because the domain axis is built from
whatever happens to be present.

### Sierra Leone cannot be evaluated at all

All 72 Sierra Leone rows are NA. With 14 districts, 5-fold CV leaves ~11 in
training, below the harness's `length(tr) < 12` floor, so every fold is
skipped; leave-one-region-out over 4 provinces fails the same way. This is the
same constraint the chiefdom analysis hit from the other side — Sierra Leone
has the best-measured areas in the project and too few of them to fit anything
within country.

### What this changes

1. **Retract the dietary verdict.** "Diet is not an important predictor" was
   measuring the 4-country intersect, not diet. The honest statement is that
   the dietary block predicts about as well as everything else and is
   redundant with it.
2. **A dietary-only model is a viable cheap instrument.** 66 survey-derived
   columns reach 77% of a 533-column model that needs the full GEE stack. For a
   programme that already runs an HCES or DHS but has no remote-sensing
   pipeline, that is the more useful result than the marginal-contribution null.
3. **Gambia deserves a closer look.** Diet alone at 0.602/0.632 is the highest
   single-block score seen anywhere in this folder. Whether that is real or a
   30-district artefact is worth one probe.
4. **Fixing Ghana's HCES block is now doubly motivated** — it restores 22
   columns to LOCO *and* removes the one country where the dietary block is
   actively harmful.
