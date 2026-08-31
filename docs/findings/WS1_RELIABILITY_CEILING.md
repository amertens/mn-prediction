# WS1. Validating the reliability ceiling

Section 9 of `docs/SESSION_FINDINGS_FOR_REVIEW.md` states the question and leaves
it open: measured skill exceeds the estimated ceiling in about half of all cells
and by a factor of two for the best anchoring arm, so either the ceiling is
biased low or the evaluation is optimistic.

## The one-paragraph answer

**The ceiling was biased low, and that is the larger part of the explanation.**
The analytic estimator understates the attainable correlation by a mean of
**0.161** when measured against a known truth, while the split-half estimator
understates it by **0.007** (source: `results/tables/reliability_simulation.csv`,
columns `bias_analytic` and `bias_empirical`, four cells by four reliability
settings). On the real data the analytic median `r_max` is **0.1315** against an
empirical **0.6132**, with the empirical larger in **21 of 24** cells (source:
`results/tables/reliability_analytic_vs_empirical.csv`). Two mechanisms account
for it: a design effect fixed at 1.5 where the value that would reconcile the
formula with the measurement has median **0.969**, and a `max(0, .)` truncation
that reports `r_max` as exactly zero in **10 of 24** cells against **3 of 24**
empirically. The evaluation is **also** optimistic, separately and for the
anchored arms only, and WS2 measures that; it is not needed to explain Section 9.
Under the empirical ceiling and with the `r_max > 0.05` guard the production code
already applies, no arm has a mean `r_share` above 1 and the count of cells above
1 falls from **25 to 7** of 118 (source: `results/tables/r_share_revised.csv`).

## WS1a. Empirical split-half reliability

`results/tables/reliability_empirical.csv`, 200 random splits per cell, seed
20260901, Spearman-Brown corrected.

Two schemes are computed and the difference matters. The **within** scheme splits
respondents inside each cluster, so both halves draw on every cluster and
between-cluster variation is not reproduced; it is an **upper bound**. The
**cluster** scheme assigns whole clusters to halves and is design-faithful, but
needs a district to hold at least two clusters, which the median district does
not in three of the four countries. Medians across cells:

| Scheme | Median `r_max` | Cells | Note |
|:---|---:|---:|:---|
| Analytic, deff 1.5 | 0.1315 | 24 | the published ceiling |
| Empirical, within | 0.6132 | 24 | upper bound |
| Empirical, cluster | 0.4689 | 23 | design-faithful, median 13 districts per cell |

The design-faithful scheme sits below the within scheme, as it must, and still
far above the analytic one. The honest statement of the ceiling is therefore a
range of roughly **0.47 to 0.61**, against a published **0.13**.

Three cells return an empirical `r_max` of exactly zero: Sierra Leone
`women_b12`, Malawi `women_vitA` and Malawi `women_b12`, each with the
Spearman-Brown correction negative in 52 to 100 percent of splits. Those cells
genuinely carry no recoverable between-district signal. The analytic estimator
says the same of ten cells, seven of which the measurement contradicts.

## WS1b. Where the analytic estimator goes wrong

The formula is
`lambda = max(0, Var(p) - mean(deff x p(1-p)/(n-1))) / Var(p)`, with `r_max` its
square root (`R/area_reliability.R`). Three defects, in order of size.

**1. The design effect is assumed, and the assumption is too large.** Solving the
formula for the design effect that would reproduce the measured reliability gives
a median of **0.969** across cells, against the assumed **1.5** (source:
`results/tables/reliability_analytic_vs_empirical.csv`, column `implied_deff`).
Assuming 1.5 subtracts roughly half again as much sampling variance as the data
support, which pushes `lambda` down and, in eleven cells, through zero. The file
comment already concedes the design effect cannot be estimated at this resolution
because most districts hold a single PSU; the measurement now shows the chosen
substitute is wrong in a known direction.

**2. The estimator and the estimate do not match.** `svy_prev` is a design-based
`srvyr::survey_mean()` over PSU and strata, while `n_svy` is `dplyr::n()`, the
raw unweighted respondent count (`R/admin2_analysis.R`). The variance formula is
a simple-random-sample binomial one applied to a weighted estimator.

**3. The truncation is not symmetric, and it feeds a denominator.** `lambda` is
floored at zero before the square root, so `r_max` cannot fluctuate below zero
while a correlation can. Ten of 24 cells sit exactly on that floor. Every
`r_share` divides by this quantity, and the anchoring script divides whenever
`r_max > 0` with no lower guard, which is how one cell reached `r_share` **13.61**
on an `r_max` of **0.030**.

## WS1c. Simulation with a known truth

`results/tables/reliability_simulation.csv`. The real survey structure is kept
(actual districts, clusters and respondents per cluster) and only the outcome is
simulated, so the quantity the ceiling should equal, `cor(true, observed)`,
becomes observable.

| Between-district SD (logit) | Attainable `r` | Analytic | Empirical | Analytic bias | Empirical bias |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.291 | 0.084 | 0.306 | **-0.207** | +0.015 |
| 0.5 | 0.603 | 0.331 | 0.599 | **-0.272** | -0.004 |
| 1.0 | 0.831 | 0.717 | 0.841 | **-0.114** | +0.010 |
| 1.5 | 0.913 | 0.862 | 0.920 | **-0.050** | +0.007 |

(means over the four cells, from the per-cell rows in the same file.)

Two things follow. The bias is **largest exactly where it matters**: at the low
end, where these cells actually sit. And the analytic estimator returns exactly
zero in **56 to 97 percent** of replicates at the lowest setting, while the
attainable correlation there is 0.16 to 0.37. The statement that four cells
"have no detectable signal above sampling noise" describes the estimator, not the
data.

## WS1d. Noise injection, as a test

`tests/testthat/test-reliability-ceiling.R`, five tests, running on the Ghana
2014 DHS district education indicator. They pin: that injected binomial noise
reduces the oracle correlation monotonically; that the analytic ceiling is
biased downward; that the split-half estimator recovers the attainable
correlation within 0.15; and that the analytic estimator returns exactly zero in
a regime where real signal remains. The suite runs in seconds and skips cleanly
when the DHS file is absent.

## WS1e. Post-stratification

Run branch, with a reframe the data forced.

The concern is real in general: an estimator that sees respondent-level
information and is averaged over the realised sample shares a sampling
realisation with the yardstick, and no ceiling on transportable skill covers
that. For the **proxy** arm in these surveys the channel is closed by the data.
`Xvars` excludes the whole `GW` domain, so predictions are constant within a
cluster, and the realised sample can enter only through cluster weights.
Measured within-district coefficient of variation of the survey weight:

| Country | Districts | Single-cluster districts | Median within-district weight CV |
|:---|---:|---:|---:|
| Gambia | 30 | 17 | 0.0000 |
| Ghana | 75 | 62 | 0.0000 |
| Malawi | 87 | 74 | 0.0000 |
| Sierra Leone | 14 | 0 | 0.5164 |

In three countries the weight is **constant within a district**, so weighted and
unweighted district means are identical and post-stratification is an exact
no-op. On the Ghana smoke cell all six aggregation-by-target combinations return
the same `r` to four decimal places, which is that identity showing up as a
measurement (source: `results/tables/poststratified_aggregation_SMOKE.csv`).

Only Sierra Leone can show a difference, and the full sweep was redirected to it
rather than spending hours confirming an identity in the other three. WorldPop
holds total population only for these countries, with no age-sex structure, so
the composition used is the survey's own weighted composition, labelled as such.

**In Sierra Leone the difference is negligible, and it has the wrong sign to
support the concern.** Over its six cells, post-stratified minus realised
aggregation gives a mean change in `r` of **+0.0027** against the district sample
mean, **+0.0025** against `svy_prev` and **+0.0040** against the survey-weighted
district mean, with the post-stratified version *lower* in only one or two of six
cells (source: `results/tables/poststratified_aggregation_SLE.csv`). If the
realised-sample channel were inflating the correlation, removing it would lower
`r`. It raises it, by an amount indistinguishable from zero.

**Conclusion for WS1e.** Aggregating over the realised sample rather than a fixed
population composition contributes nothing measurable to the proxy arm's
correlation, in the one country where it is capable of contributing anything.
This eliminates the second of the two candidate explanations for Section 9 as
far as the proxy arm is concerned, and leaves the biased ceiling as the
explanation.

**This does not dispose of the concern for the questionnaire arm**, which does
carry respondent-level covariates. WS3 handles that arm.

## WS1 acceptance: `r_share` recomputed for every arm

`results/tables/r_share_revised.csv` and `..._summary.csv`.

| Arm | Published mean | With the 0.05 guard | Under the empirical ceiling | Median, empirical |
|:---|---:|---:|---:|---:|
| Admin-1 anchor (hard) | 2.06 | 1.17 | **0.75** | 0.68 |
| Admin-1 anchor (shrunk) | 1.64 | 0.96 | **0.56** | 0.51 |
| Admin-1 fit to Admin-2 | 1.77 | 0.77 | **0.38** | 0.52 |
| National anchor | 1.34 | 0.68 | **0.35** | 0.32 |
| No anchor | 1.35 | 0.69 | **0.33** | 0.34 |

Both corrections are needed. Applying the `r_max > 0.05` guard that
`add_reliability_columns()` already uses removes the cells with a near-zero
denominator and takes the hard arm from 2.06 to 1.17. Replacing the biased
ceiling takes it from 1.17 to 0.75. Cells with `r_share > 1` fall from **25 to
7** of 118.

The published mean of 2.06 was never a meaningful summary: it is a mean of
ratios whose denominator can be 0.030.

## Reviewer pass, statistical

**Fold construction.** WS1a fits no model, so folds do not arise. The
correlations are between two halves of the same survey, across districts.

**Survey weights.** Split-half prevalences are survey-weighted, matching
`svy_prev`, which is the estimand. The analytic estimator's mismatch between a
weighted estimate and an unweighted `n` is reported above as defect 2 and is not
reproduced here.

**Seed coverage.** WS1a seed 20260901, 200 splits. WS1c seed 20260902, 100
replicates per setting, with the inner split-half seeded at `seed + r` so
replicates do not share a split. WS1e seed 20260903. WS1d seeds are fixed
literals per test.

**Identical cells.** The analytic and empirical ceilings in
`reliability_analytic_vs_empirical.csv` are computed on the same 24 cells, and
the analytic one is recomputed inside the same run rather than read from the
published table, so a difference cannot come from a stale input. The `r_share`
recomputation joins on country and outcome and preserves all 118 arm-cell rows.

**A limit worth stating.** The **within** scheme is an upper bound by
construction. If the true reliability were near the design-faithful 0.47 rather
than the 0.61 the primary scheme reports, every conclusion above holds with
smaller margins: the analytic ceiling would still be biased low by roughly 0.34
in the median cell, and the hard anchor's `r_share` would be about 0.87 rather
than 0.75, still below 1. No conclusion in this workstream turns on which of the
two empirical schemes is used.

**What this does not establish.** That the anchored arms' evaluation is honest.
The ceiling being biased low removes the paradox; it does not show the anchoring
gain is real. WS2 tests that directly and finds it is not.

## Reviewer pass, reproducibility

**Targets graph.** WS1 adds no target. Its three scripts sit in
`scripts/accuracy_impact/` alongside the existing `scripts/covariates/` scripts
that write to `results/tables/` outside the graph. `tar_manifest()` is unchanged
from the 845 targets WS7a verified.

**Stamp targets.** WS1d reads
`data/DHS/clean/Ghana_2014_dhs_custom_admin2_wide.rds`, which is a tracked file,
and skips with an explanatory message when it is absent. No untracked input is
consumed.

**Paths and determinism.** All scripts use `here()`, contain no absolute path
and no `setwd()`. Every stochastic call is seeded and re-running reproduces the
tables.

**Smoke profile.** `PROFILE=smoke` runs Ghana `child_iron` with `B = 50` for
WS1a and `R = 20` for WS1c; `PROFILE=sierraleone` runs WS1e on the only country
where its answer can be non-trivial. All three were exercised before the full
runs. Full runtimes measured: WS1a and WS1b together about four minutes over 24
cells and both schemes; WS1c about twelve minutes over four cells and four
settings; WS1e about five minutes per cell, which is why it was redirected.

**Test suite.** `Rscript tests/testthat.R` reports 37 passing, 0 failing, 0
skipped, up from 28 after WS7a.
