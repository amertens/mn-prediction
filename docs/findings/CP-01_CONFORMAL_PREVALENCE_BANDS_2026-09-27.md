# CP-01. Calibrated 90% bands for the planning prevalence, by split conformal on out-of-fold residuals

27 September 2026. Design written before the script ran (the post-audit rule).
Script: `scripts/protocol_v2/66_conformal_prevalence_bands.R`. Results:
`results/tables/protocol_v2/conformal_prev_cells.csv` (+ `_residuals`,
`_districts`). Dashboard integration: builder 05 attaches the calibrated
columns to `admin2_index.rds`; the map's exceedance layer and the district
profiles switch from the stability quantities to these.

## Why

The dashboard's prevalence bands and WHO-threshold exceedance chances (UE-01)
come from refitting on resampled training districts. They measure stability,
and the decks say so — but VZ-01 showed what stability means for coverage
(38% for rank intervals under a held-out country), and the survey-planning
tools this project is heading toward (the LQAS-style calculator) would consume
these bands arithmetically. Before anything multiplies them, they need to
cover.

## Target of coverage — stated up front

The bands are calibrated to cover **the survey's own measured district value**
(the weighted district prevalence the protocol scores everything against),
not an unobservable noiseless truth. That is the checkable claim, it includes
the survey's own sampling noise, and it is the quantity a confirmatory survey
visit would produce. Nothing here claims coverage for a country outside the
four (the transport analogue is future work; CIV shows no prevalence product,
so nothing there needs this).

## Design

Per cell (country × outcome, all cells with ≥ 12 surveyed districts,
including the Malawi selenium/iodine cells via the shared targets reader):

1. **Out-of-fold predictions of the deployed quantity.** 5-fold by district,
   10 draws (the protocol's in-fill folds, `make_folds_v2`). Per fold: fit the
   calibrated index (`domain_index_cal`, IS-01) on the training districts,
   apply to every district of the country, re-solve the national anchor on
   population exactly as deployment does, and keep the held-out districts'
   anchored prevalences. Average over draws: one honest p-hat per surveyed
   district. (The deployed map uses the full-data fit; the draw-averaged
   out-of-fold predictor is its closest honest stand-in, and the direction of
   the approximation — fold fits see 80% of the data — makes the bands
   conservative, not optimistic.)
2. **Residuals.** e_i = p_survey,i − p-hat_i, one per surveyed district.
3. **Band: split conformal, constant width per cell.** Half-width = the
   ceil(0.9 · (n+1))-th smallest |e|; with small n this hits max|e| and the
   band is honestly wide. Bands are point ± half-width on the probability
   scale, clamped to [0, 1]. Width in percentage points is the legible unit
   this project already reports (AR-01, MAE).
4. **Exceedance: conformal predictive distribution.** For any district and
   threshold t, P(≥ t) = (#{e_j : p-hat_i + e_j ≥ t} + 0.5) / (n + 1) — the
   split-CPS estimate from the same residual set, calibrated by construction
   against the same target. Computed for the WHO "moderate-plus" and "severe"
   cuts where they exist.
5. **Checks, reported per cell:**
   - leave-one-out coverage of the 90% band (each district's residual against
     the quantile of the others') — should sit near 0.90;
   - the UE-01 **stability** band's empirical coverage of the survey value on
     the same districts — the number that motivated this work;
   - median |e| and the conformal half-width, in points.

## Pre-stated readings

- The stability bands under-cover (that is the premise); the conformal bands'
  LOO coverage lands in 0.85–0.95 in most cells.
- Widths will be wide where the model is weak (Malawi vitamin A, zinc) and
  the dashboard will show them anyway: a wide honest band is the product
  working, not failing.
- The calibrated exceedance chances are less extreme than the stability ones
  (they inherit model error, not just refit spread). Where the two disagree,
  the calibrated one is displayed and the stability one retires from the
  prevalence side (it remains the right object for rank firmness).

## Results (same day)

27 cells, all with LOO coverage 0.90-0.93 (the conformal machinery does what
it says). The pre-stated readings all held, the first one more dramatically
than expected:

- **The stability bands covered a mean of 8% (0-23% per cell)** — the
  prevalence-side stability bands were near-meaningless as uncertainty, worse
  than the rank-side 38% of VZ-01. They are retired from the prevalence
  displays as of this note.
- **Honest widths are wide: median half-width 26 pp.** Iron cells sit at
  ±30-37 pp; Ghana folate ±49; Malawi selenium ±56-59 (a cell that ranks
  better than the survey's regional averages, 0.43-0.45 vs 0.24-0.32, yet
  carries the widest level band; corrected 27 Sep: it is not the strongest
  ranking cell, which is Malawi women's B12 at 0.70);
  the rare-prevalence vitamin A cells are ±2-13 pp because the outcome is
  near zero. Much of the width is the coverage target's own noise: the bands
  cover the survey's measured district value, which mostly rests on one or
  two clusters (the ceiling analysis quantifies that share).
- **Calibrated exceedance is less extreme but not useless.** Ghana child
  vitamin A: districts at >=80% chance of sitting at or above the WHO
  moderate line fall from 204 of 260 (stability) to 67 (calibrated); the
  honest uncertain middle grows from 51 to 193. Confident calls survive where
  the point sits far from the line (Tamale: band 1-59%, chance >= the 10%
  line still 98%).

Consequence for the LQAS-style calculator: its main blocker (an uncalibrated
prior) is cleared — the CPS residual sets in `conformal_prev_residuals.csv`
are exactly the calibrated predictive distribution it needs. The remaining
blockers are the per-district design-effect problem and the missing WHO
conventions for selenium/iodine.

## What changes on the dashboard

`admin2_index.rds` districts gain `prev_cal_lo`, `prev_cal_hi`,
`p_modplus_cal`, `p_severe_cal`, `conformal_n`. The map's exceedance layer,
the district click-through and the district profiles switch to the calibrated
quantities with the coverage check quoted; the rank stability range keeps its
stability label (ranks stay class-probability territory, as the deck's
calibration slide argues). Unsurveyed districts get the same cell width —
stated as an extrapolation of the surveyed districts' error, the standard
split-conformal reading.
