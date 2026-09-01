# The target estimand, and which ceiling bounds which estimator

Written for WS1f. The purpose is to stop one number being used to bound a
quantity it does not bound, which is how Section 9 of
`docs/SESSION_FINDINGS_FOR_REVIEW.md` came to report measured skill at twice its
own estimated ceiling.

## 1. The estimand

**The district population prevalence.** For district *d* and a given
country-outcome cell, the quantity every model in this project is trying to
estimate is

> the proportion of the district's target population (preschool children, or
> non-pregnant women of reproductive age, as the outcome specifies) that would be
> classified deficient if the whole population were measured under this
> project's assay and inflammation-adjustment protocol.

Three things follow immediately, and each has been a source of confusion.

**It is a population quantity, not a sample quantity.** The survey's
design-based estimate `svy_prev` is an estimate of it, and carries sampling
error. A model scored against `svy_prev` is scored against a noisy yardstick,
which is what makes a reliability ceiling necessary at all.

**It is protocol-specific.** VMNIS records what national surveys report under
*their own* assays and adjustments, which is a different quantity. Section 6.4
attributes the composition failure to a weak model; the more accurate statement
is that it is a model of a different estimand. This is the same cross-survey
level offset already documented as defeating leave-one-country-out transport.

**It is not the achieved sample's mean.** The unweighted mean over the
respondents a survey happened to interview in a district is a third quantity
again. Where survey weights vary within a district the three differ; where they
do not, they coincide. See section 4.

## 2. The two ceilings

A reliability ceiling answers: given that the yardstick is noisy, how large can
a correlation against it be? Write `pi_d` for the estimand and `p_d` for the
survey estimate. The ceiling is

> `r_max = sqrt(lambda)`, `lambda = Var(pi_d) / Var(p_d)`

and it bounds `cor(pi_d, p_d)`, the correlation an oracle predictor of the truth
would score against the observed yardstick.

**The analytic ceiling** (`admin2_reliability()`, `R/area_reliability.R`)
estimates `lambda` by subtracting an assumed binomial sampling variance,
inflated by a design effect fixed at 1.5, from the observed between-district
variance, then truncates at zero.

**The empirical ceiling** (`split_half_reliability()`,
`R/reliability_empirical.R`) estimates the same `lambda` by splitting each
district's respondents into halves, estimating the prevalence twice, correlating
the halves across districts, and applying the Spearman-Brown correction. It
assumes no design effect and no binomial model.

### Which is trustworthy

Measured, not argued. WS1c simulates outcomes over the real survey structure
with the truth known by construction, so the attainable correlation `r_oracle`
is observable:

| Estimator | Mean bias against `r_oracle` |
|:---|---:|
| Analytic, deff 1.5 | **-0.161** |
| Empirical split-half | **+0.007** |

(source: `results/tables/reliability_simulation.csv`, columns `bias_analytic`
and `bias_empirical`, over four cells and four reliability settings.)

The analytic ceiling is biased low, and severely so where reliability is
genuinely low: at the lowest setting it returned exactly zero in 56 to 97
percent of replicates while the oracle correlation was 0.16 to 0.37. On the real
data the analytic estimator reports `r_max` exactly 0.000 in **10 of 24** cells
against **3 of 24** for the empirical one
(source: `results/tables/reliability_analytic_vs_empirical.csv`).

**Use the empirical ceiling.** The analytic one may be quoted for continuity
with earlier documents, and must then be labelled as a lower bound.

## 3. What each ceiling does and does not bound

This is the section that matters for Section 9.

**Both ceilings bound an estimator whose predictions are independent of the
sampling realisation.** An area-level model that sees only Admin-2 covariates
qualifies: its prediction for a district would be the same whichever respondents
the survey had drawn. For such an estimator, a measured correlation above
`r_max` is evidence that `r_max` is wrong, or that the evaluation is optimistic.

**Neither ceiling bounds an estimator that sees respondent-level or
cluster-level information and is then averaged over the realised sample**,
unless its predictions are post-stratified to a fixed population composition.
Such an estimator's district aggregate depends on who was interviewed, and so
does the observed district prevalence. The two share a sampling realisation, and
shared sampling noise contributes correlation that no ceiling on transportable
skill was ever meant to cover.

**Consequence for the anchoring arms.** An arm that rescales predictions so their
regional mean equals the survey's own regional estimate has, by construction,
used the survey. Its correlation against the survey's district estimates is not
bounded by a ceiling derived on the assumption of independence. WS2 tests how
much of the anchoring gain this accounts for; until then, `r_share` for the
anchored arms should be read as a descriptive ratio and not as a fraction of an
attainable maximum.

## 4. Post-stratification, and why it changes nothing here

The general concern in section 3 is real. Its empirical force in this project is
small, for a reason specific to these surveys.

The primary predictor set `Xvars` is `Xvars_all` minus the entire `GW` survey
domain (`R/data_prep.R`), so it holds **no respondent-level covariate**.
Predictions are constant within a cluster. The realised sample can therefore
reach the district aggregate only through the weights clusters receive. Measured
across the four surveys, the median within-district coefficient of variation of
the survey weight is:

| Country | Districts | Districts with one cluster | Median within-district weight CV |
|:---|---:|---:|---:|
| Gambia | 30 | 17 | 0.0000 |
| Ghana | 75 | 62 | 0.0000 |
| Malawi | 87 | 74 | 0.0000 |
| Sierra Leone | 14 | 0 | 0.5164 |

(source: measured from the `outcome_data_*` targets, first two outcomes per
country; reproduced by `scripts/accuracy_impact/ws1e_poststrat.R`.)

In three of the four countries the survey weight is **constant within a
district**, so the weighted and unweighted district means of any quantity are
identical and post-stratification is an exact no-op. Only Sierra Leone, which
has no single-cluster district and a within-district weight CV above 0.5, can
show a difference at all.

This does not weaken section 3. It says that for the **proxy** arm the channel is
closed by the data. The **questionnaire** arm of the individual anchor does carry
respondent-level covariates, and there the concern applies with full force.

## 5. Practical rules

1. Score against `svy_prev`, the design-based estimate of the estimand, not
   against the achieved sample mean.
2. Quote `r_max_emp` from `results/tables/reliability_empirical.csv`. If the
   analytic ceiling is quoted, label it a lower bound and give the empirical one
   beside it.
3. Report `r_share` as a **median** across cells, never a mean. It is a ratio
   whose denominator can approach zero, and a mean over such ratios is dominated
   by whichever cell has the smallest denominator. Section 9's mean of 2.06 for
   the hard anchor includes a cell whose `r_max` is 0.030 and whose `r_share` is
   consequently 13.61.
4. Suppress `r_share` entirely where the denominator is below a stated cut. The
   production code uses `r_max > 0.05` in `add_reliability_columns()`; the
   anchoring script did not, which is why the two disagree.
5. Do not use `r_share` to compare an anchored arm against an unanchored one.
   The anchored arm has seen the survey; the ceiling assumes it has not.
