# WS-I. The cross-country estimator: anchored level, covariate tilt

`results/tables/loco_decomposition.csv`,
`results/tables/loco_decomposition_summary.csv`. 22 cells, tilt swept rather
than tuned, seed 20260926.

## Leave-one-country-out does not fail the way leave-one-region-out fails

Within a country, covariate models lose to a covariate-free shrinkage rule but
the loss is measured in fractions of a correlation point. Across countries the
failure is of a different kind, and the earlier `relative_target_loco`
experiment had already isolated it:

| Arm | mean r | MAE pp | \|bias\| pp |
|:---|---:|---:|---:|
| absolute — covariates, no anchor | 0.201 | 13.50 | 10.67 |
| relative + **true** national anchor | 0.140 | 10.81 | 7.00 |
| relative + **blind** predicted anchor | 0.141 | 14.23 | 11.01 |
| **anchor only** — flat national, no covariates | — | **9.85** | **5.68** |

Two readings follow, neither of them about the learner.

1. Applying the held-out country's national prevalence **flat to every
   district** beats every covariate model on error.
2. The two relative arms differ *only* in whether the anchor is known or
   predicted, and they differ by **3.4 pp**. That gap is the LOCO problem.

## The decomposition, scored

`loco_anchor_tilt()` estimates the two failing quantities separately:

```
prediction = anchor + tilt * spread * (within-country covariate rank signal)
```

Scored over 22 cells with the harmonised predictor set, sweeping the tilt:

| Anchor | Tilt | r | MAE pp | \|bias\| pp |
|:---|---:|---:|---:|---:|
| known | **0.00** | n/a | **8.94** | **1.02** |
| known | 0.20 | 0.080 | 9.02 | 1.04 |
| known | 0.35 | 0.082 | 9.26 | 1.05 |
| known | 0.50 | 0.080 | 9.61 | 1.10 |
| known | 0.75 | 0.078 | 10.49 | 1.33 |
| known | 1.00 | 0.080 | 11.53 | 1.73 |
| predicted | 0.00 | n/a | 17.55 | 14.18 |
| predicted | 1.00 | 0.074 | 18.69 | 14.31 |
| *un-anchored transported model* | — | 0.082 | 17.49 | — |

**The `r` column is constant by construction and must not be read as a
trade-off.** `anchor + tilt * spread * z` is an AFFINE transform of the
covariate signal `z` whenever the anchor is a scalar, and Pearson correlation is
invariant under affine transforms with positive scale. So the tilt cannot change
the ranking at any setting. Verified on these results: `r` is *exactly* constant
across the sweep in **9 of 22 cells** and near-constant in the rest, with a
median spread of **0.0034**. The residual movement is entirely the
`pmin(pmax(., 0), 1)` clipping, which is only weakly monotone at the boundary --
which is why Sierra Leone child_iron, whose raw predictions run off the ends,
varies most at 0.076.

This was flagged by an independent implementation on another branch
(`solve_anchor_shift()`), whose author states it exactly: a monotone shift is
"MATHEMATICALLY INCAPABLE of improving RANKING". The same caution applies to any
top-k agreement statistic computed after such a shift -- a 100 percent overlap
would be guaranteed by construction rather than measured.

**What the sweep therefore shows** is narrower than first written, and the
conclusion is unchanged but better founded: the tilt is a **pure error knob**.
It cannot buy ranking, so the only reason to raise it is if it lowers error --
and pooled it does not.

**The anchor is worth 8.6 pp. Bias collapses from 14.18 to 1.02 on the anchor
alone.** The covariate tilt is worth nothing pooled, and every increase in it
makes the error worse.

## Where the tilt does earn its place

| Held out | best tilt | MAE pp | r (constant in tilt) |
|:---|---:|---:|---:|
| Gambia | **0.50** | 9.55 | 0.309 |
| Ghana | **0.35** | 11.75 | 0.267 |
| Malawi | 0.00 | 9.64 | — |
| Sierra Leone | 0.00 | 4.05 | — |

"Best tilt" here is an **error minimum only**; the r column is shown for context
and is the covariate signal's own correlation, identical at every tilt.

The same two countries that beat the null within-country, and the same two the
simulated tournament ranks highest. **Three independent analyses, one ordering.**

## Why the weights are fixed and not learned

Within a country there are 14 to 87 districts to fit a stack on. Across
countries there are **four**, so a leave-one-country-out stack fits its weights
on three points. Transportability is not estimable at that n, and a stack that
pretended otherwise would be fitting noise and presenting it as adaptivity. The
tilt is a declared parameter and the sweep is reported rather than optimised
away.

## What this does not fix

The level offset is a cross-survey calibration problem — raw ferritin differs
sixfold across these countries — not a shrinkage problem inside one. The
decomposition **routes around** it by taking the level from an anchor; it does
not repair it. A country with no national estimate at all is served by the
`predicted` rows above, at 17.55 pp, which is no better than the un-anchored
model.

That makes the highest-leverage open question in the project not "which
learner" but **how well a country's national level can be estimated without a
survey** — which is the VMNIS national model.

## Reviewer pass

**Covariates usable in every country.** The transported model is restricted to
predictors that are present and non-constant in all four, so the fold does not
silently change the predictor set.

**Median imputation from training rows only.**

**The `spread` is borrowed.** A transported ranking has no scale of its own, so
the within-country SD is taken from the training countries' median. That is an
assumption, it is most likely to be wrong for a country whose districts are
unusually equal or unequal, and it is recorded per row as `spread_assumed`.

**The blind anchor here is the weakest honest one** — the training countries'
precision-weighted mean. VMNIS supplies a better anchor where it covers the
country, and the gap between the two columns is the value of having a national
estimate at all.

**Additive.** Both output tables are new.

## Addendum: an extrapolation diagnostic this analysis lacks

An independent implementation on another branch reports that its Sierra Leone
LOCO predictions were **extrapolations rather than shrinkage** -- Sierra Leone's
predictors falling outside the training countries' joint range, with four of
fourteen districts receiving negative predicted prevalence. That work measures
this with Stahel-Donoho outlyingness, having first verified that halfspace depth
saturates at zero for essentially every point at p=16, n=180.

**That specific finding does not reproduce here, and the disagreement is
informative rather than contradictory.** Under the production harmonised
predictor set, checking how many held-out districts fall outside the training
per-covariate range:

| Held out | districts with >10% of covariates out of range | negative predictions |
|:---|---:|---:|
| **Malawi** | **100.0%** | 0 |
| Ghana | 30.7% | 0 |
| Gambia | 0.0% | 0 |
| **Sierra Leone** | **0.0%** | 0 |

No negative predictions occur in any of 16 cells, and Sierra Leone is the least
extrapolated country by this measure while Malawi is fully extrapolated.

**This cannot refute the other result**, for a reason that is the whole point of
their method choice: a marginal range check asks whether each covariate is
individually in range, while a depth measure asks whether the *joint* point is
plausible. A district can sit inside every marginal range and far outside the
joint distribution. The configurations also differ -- 294 harmonised predictors
and a screened ridge here, against roughly 16 predictors on the teaching subset
there.

What this check does establish independently is that **Malawi is 100 percent
extrapolated** under the production covariate set, which nothing in this project
had noticed and which is a candidate explanation for Malawi's weak transport
(r 0.095, and the lowest tilt benefit of the four countries).

Porting a proper joint-outlyingness diagnostic to the production LOCO arm is the
right follow-up. It is not done here.
