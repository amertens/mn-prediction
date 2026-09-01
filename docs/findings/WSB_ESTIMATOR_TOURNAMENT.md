# WS-B. The shipping estimator, decided against simulated truth

`results/tables/estimator_tournament_truth.csv`, 8 cells spanning all four
countries and both nutrient families, 100 replicates per setting, four levels of
covariate signal, seed 20260922.

## Why the yardstick had to change

Every estimator comparison in this project scored predictions against the
survey's district estimates. WS1 measured the reliability of that yardstick at a
median of 0.613, so roughly forty percent of its variance is sampling noise. Two
consequences followed and neither was visible from inside that frame:

1. The **direct district estimate** has a correlation of exactly 1.000 with the
   yardstick by construction, so it was never treated as a candidate estimator.
2. Any estimator that **shrinks** is penalised for failing to reproduce the
   yardstick's own noise, which is precisely backwards.

Simulating the truth removes both problems.

## The result

Pooled over 8 cells, scored against truth:

| ρ | Metric | direct | **eb_blend** | ridge | region_tilt | flat_region |
|---:|:---|---:|---:|---:|---:|---:|
| 0.20 | r | **0.724** | 0.583 | 0.188 | 0.066 | 0.063 |
| 0.20 | MAE pp | 6.71 | **6.25** | 9.52 | 9.85 | 9.34 |
| 0.35 | r | **0.718** | 0.636 | 0.324 | 0.149 | 0.143 |
| 0.35 | MAE pp | 6.84 | **6.05** | 9.00 | 9.28 | 8.75 |
| 0.60 | r | **0.726** | 0.718 | 0.522 | 0.320 | 0.315 |
| 0.60 | MAE pp | 6.79 | **5.54** | 7.70 | 7.71 | 7.28 |

**The EB blend wins on error at every level of covariate signal.** The direct
estimate wins on correlation at every level. Head to head, the blend beats the
flat regional mean on correlation in 8, 8, 8 and 7 of 8 cells and on error in 7,
6, 6 and 6 of 8.

**The covariate arms lose everywhere in the measured range.** The ridge reaches
0.188 to 0.324 at the ρ this project has actually measured, against 0.583 to
0.636 for the blend, and carries a systematic **-3.1 to -3.3 pp bias** that is
the same defect WS-D's gate catches in six real cells.

## What was confirmed, and what was not

**Confirmed.** The EB blend is the best available estimator of district
prevalence, it dominates the covariate-free regional mean that the previous
session recommended, and it dominates every covariate model.

**Not confirmed.** The smoke run on one cell had the blend overtaking the direct
estimate on correlation at ρ 0.60. The eight-cell run does not reproduce that:
the direct estimate wins on correlation at all four ρ. The blend's advantage is
on error, not on ranking.

**A limitation of the simulator that biases against the blend.** The simulated
truth is `sqrt(rho) * z + sqrt(1-rho) * e`, where `z` is the covariate axis and
`e` is independent district noise. It contains **no explicit regional block
structure**. WS-E2 measured that about **40 percent** of real district-level
variance sits between regions, so in the simulated world the region mean carries
less information than it does in reality, and an estimator that shrinks toward
it is handicapped. The blend's correlation is therefore an underestimate of what
it would achieve on real data, and the direct estimate's advantage on ranking is
an overestimate. Adding a region random effect to the generator is the fix and
is not done here.

## The shipping decision

**Ship the empirical Bayes blend for prevalence.** Report district rankings from
it with rank intervals, and read the ranking with the caveat above.

The two claims are separable and both are honoured in
`docs/FITNESS_FOR_USE.md`: the district layer is a ranking rather than a
prevalence, and where a prevalence is needed the blend is the estimator with the
lowest error against truth.

## The defect that nearly shipped

The first implementation estimated `tau2` as `var(p) - mean(v_d)`, floored at
zero. That is the same estimator WS1 showed is biased low, and it returned
**tau2 = 0 in fourteen of twenty-four cells**, forcing `lambda = 0` everywhere
and collapsing the blend into the flat regional mean, the arm WS-B1 had just
withdrawn at a jackknifed mean r of 0.076.

Sourcing `tau2` from this session's own split-half reliability instead,
`tau2 = lambda_emp / (1 - lambda_emp) * vbar`, leaves `lambda` at zero in only
**three** cells, and those three are exactly the cells WS1a found with zero
empirical reliability. Measured on the shipped output: `tau2` comes from the
reliability route in **21 of 24** cells (source:
`results/tables/shipped_estimates_summary.csv`, column `tau2_source`).

## What shipped

`results/deliverables/district_estimates.csv` and `.rds`, 24 cells:

- **19 Admin-2 layers released, 5 suppressed** at the reliability cut of 0.30.
  The five are Malawi women_vitA and women_b12, Sierra Leone women_b12,
  child_vitA and child_iron, which is exactly the list
  `docs/FITNESS_FOR_USE.md` section 3 specifies.
- **6 cells carry `calibration_failed`** from the WS-D gate and are marked
  rather than deleted, so a consumer can see why a cell is absent from a map.
- Median `lambda` ranges 0.000 to 0.772 across cells; the median shift the blend
  applies ranges 0.00 to 9.69 pp.
- Rank intervals are wide, with a median width of 42 to 66 ranks in the
  87-district and 75-district countries. That is the honest reading: the
  district ranking separates the ends of the distribution and not the middle.

`results/deliverables/national_regional.csv`, 384 rows of national and Admin-1
prevalence with standard errors and WHO band labels, each carrying
`band_source_verified` from `config/who_thresholds.csv`, which is **FALSE for
every row** because no threshold citation in this project has been verified
against a source document.

## Reviewer pass, statistical

**Shared-noise check.** The blend's regional target is jackknifed, so a
district contributes zero of its own sampling noise to the quantity it is shrunk
toward. Without that the estimator would inherit the defect that withdrew
registers 4.6 and 4.12; a test pins it.

**The direct estimate is not shared-noise free in the ordinary sense**: it *is*
the observed estimate, so scoring it against the observed estimate is
meaningless, which is exactly why the tournament scores against truth.

**Identical replicates across estimators.** All five estimators are scored on
the same simulated draw within a replicate, so the comparison is paired.

**Eight cells, not twenty-four.** The tournament costs about half an hour per
cell because every replicate refits the ridge and the tilt once per region. The
subset is named in the script and spans all four countries and both nutrient
families. The estimator ranking is stable across all eight and across all four
ρ, which is the property the shipping decision needs.

**Seeds.** 20260922 for the tournament, 20260924 for the rank-interval
bootstrap.

## Reviewer pass, reproducibility

**Additive.** `estimator_tournament_truth.csv`, `shipped_estimates_summary.csv`
and everything under `results/deliverables/` are new.

**Dashboard untouched**, per the earlier session's constraint.

**Tests.** Seven cover the estimator, including one that fails if the moment
route for `tau2` is reinstated. Suite is 72 passing.

**Runtime.** About four hours for eight cells at 100 replicates and four ρ.

---

# Addendum, late September 2026: the tau2 route, and a target that is noisier than no target

## A confound in the arm comparison, found and fixed

The three empirical Bayes arms were not using the same definition of `tau2`:

| arm | tau2 |
|:---|:---|
| `eb_blend` | `var(p) - vbar`, the **total** between-district variance |
| `eb_covariate` | `var(p - ridge) - vbar`, residual after its own target |
| `eb_stack` | `var(p - stack) - vbar`, residual after its own target |

For a shrinkage toward target `m` the model is `theta_d = m_d + u_d`, so `tau2`
is `Var(p - m) - vbar`. Using the total variance overstates it by the
between-region component, which WS-E2 puts at about 40 percent of district
variance.

The consequence was not cosmetic. It meant that any comparison of `eb_stack`
against `eb_blend` conflated **a better target** with **a better tau2 route** --
which is the single comparison this tournament exists to make. All three arms
now use residual-after-own-target, and the total-variance route is retained as a
scored arm, `eb_blend_totvar`, so its cost is measured rather than assumed.

## The correction's own justification was half wrong

The obvious reading of the above is that `eb_blend` had been shrinking too
little everywhere. Measured across the 24 shipped cells, that is true in about
half of them and false in the rest:

| `tau2_residual / tau2_total` | cells |
|:---|---:|
| both collapse to the floor | **14** |
| residual smaller, 0.00 to 0.76 | 5 |
| residual larger, 1.19 to 1.60 | 4 |
| residual very much larger | **2** |

Sierra Leone women_iron has a ratio of **61**, and Sierra Leone child_vitA of
**107,004**.

## What that means, and it is a real limitation of the shipped estimator

A residual variance larger than the total says the jackknifed region mean is a
**worse** predictor of a district than the grand mean is. That is not a coding
error. Sierra Leone has 14 districts across 4 regions, so a district's
jackknifed regional target is computed from **two or three other districts** and
carries substantial sampling error of its own.

The residual-variance formula is correct only when the target is estimated
without error. Where districts are scarce that assumption fails, and the shipped
estimator is then shrinking toward a quantity noisier than the national mean it
could have used instead.

This is a limitation to state rather than a defect to patch here. The fix would
be to propagate the target's own variance into `tau2` -- shrinking toward the
region mean only as far as the region mean is itself known -- and to let the
national mean take weight where it is not. That is exactly what the `eb_stack`
arm's candidate set makes available, and it is one more reason to score that arm
rather than reason about it.

## The shipped estimator does not need re-issuing on this

`tau2` comes from split-half reliability in **21 of 24** shipped cells, with a
median lambda of 0.478. Both moment routes -- total and residual -- collapse to
the floor in **14 of 24**, which is the degeneracy WS-B3 already documented and
the reason the reliability route was adopted in the first place. The route this
addendum corrects is a fallback that the shipped output mostly does not take.
