# Andy Kim: starter project on the simplified subset

A sequenced plan to apply the ideas from the Mertens, Hubbard, and Kim meeting to
the demo dataset in this folder. It is written so you can start without learning
the full project pipeline. Each task connects to your existing work on
sensitivity-constrained estimation, quintile (ranking) targeting, and
parsimonious DHS proxies.

## The dataset

This folder holds every micronutrient-deficiency outcome for the four study
countries (Gambia, Ghana, Sierra Leone, Malawi), aggregated to three geographic
levels and stored in wide format (one row per area, with a small block of
columns for each outcome plus a shared set of 16 predictors). See
[`README.md`](README.md) and the
column-by-column [`data_dictionary.csv`](data_dictionary.csv).

| File | Geographic level | Rows |
|---|---|---|
| `data/mn_cluster.csv` | survey cluster (the survey's primary sampling unit, a small group of neighbouring households) | 323 |
| `data/mn_admin2.csv` | district (Admin-2) | 206 |
| `data/mn_admin1.csv` | region or province (Admin-1) | 53 |

For each outcome there are four columns, following a fixed naming pattern:

- `prev_<outcome>`: the survey-weighted prevalence of the deficiency in that area
  (a proportion between 0 and 1). This is the quantity you model.
- `n_<outcome>`: the number of people surveyed for that outcome in the area (the
  denominator). Larger n means a more reliable prevalence estimate, so n is a
  reasonable weight.
- `ndef_<outcome>`: the number found deficient (the numerator). The unweighted
  prevalence is `ndef / n`.
- `var_<outcome>`: the estimated sampling variance of the prevalence, that is, how
  much the area's prevalence estimate would bounce around if the survey were
  repeated. Districts with bigger samples have smaller variance. Using `1 / var`
  as a weight lets the better-measured districts count more. This style of weight
  is what an area-level small-area model (for example Fay-Herriot) expects.

Outcomes present, and where they were measured: `women_iron`, `women_vitA`,
`child_iron`, `child_vitA` in all four countries; `women_b12` and `women_folate`
in Ghana, Malawi, and Sierra Leone; `women_zinc` and `child_zinc` in Malawi only.
A column is left blank (NA) for a country that did not measure that outcome.

The 16 predictors are geospatial proxies (rainfall, temperature, vegetation
greenness, land use and human footprint, and malaria indicators). They are the
same for every outcome. They were chosen by ranking the full proxy set on how
much each helped predict iron deficiency, keeping only proxies available in all
four countries and dropping near-duplicates. The data dictionary describes each.

`run_superlearner_example.R` fits a SuperLearner for one outcome and then runs a
short two-outcome demonstration. Run it first to confirm your setup works.

## Two ways a model is tested here

Both evaluation setups appear throughout the tasks, so it helps to fix the terms:

- Within-country cross-validation: split one country's areas into folds, train on
  most of them, predict the held-out fold, and rotate. This measures how well the
  model works inside a country that already has survey data.
- Leave-one-country-out (LOCO): train on three countries and predict the fourth,
  which contributed no data to the fit. This imitates the real goal of predicting
  a country that has never been surveyed, and it is much harder. The gap between
  the two is the central story of the project.

## Task 1: decision-focused targeting of the highest-burden districts

Why this matters. For program planning the question is usually not the exact
prevalence in a district but whether a district is among the worst affected, so
it can be prioritised for fortification, supplementation, or a follow-up survey.
A standard model is trained to make the average prediction error small across all
districts, which spreads attention evenly. We can instead tune the model so it is
good at the specific job of flagging the highest-burden districts. This is the
same structure as your sensitivity-constrained and quintile work, now applied to
geographic areas.

How a SuperLearner makes this possible. A SuperLearner first fits several base
models (for example a linear model, a penalised regression, and a random forest),
then combines them. The base models are trained on ordinary error criteria, but
the combination step is free to optimise a different objective. So you can leave
the base models as they are and only change how they are blended.

What to do. Work at the district level (`mn_admin2.csv`) and start with iron in
women (`prev_women_iron`). Define the planning target two ways: districts above a
WHO public-health threshold, and separately the top 20 percent most affected
districts within each country. Fit the base models, then choose the combination
weights to maximise how many truly high-burden districts the model catches
(recall, also called sensitivity) while holding the false-alarm rate below a set
limit. Evaluate both within-country and with LOCO.

What success looks like. A table or curve showing the share of high-burden
districts correctly flagged, within a country and when transported to a new one,
compared with an ordinary error-minimising model. The reason to start here: the
data supports ranking districts better than it supports pinning down exact levels,
this matches the real planning need, and it reuses the base-model fits.

## Task 2: anchoring transported predictions to a known national prevalence

Why this matters. When a model trained in some countries is applied to a new one,
the ranking of districts often survives but the overall level drifts, because the
model pulls every prediction toward the average of the countries it learned from.
Many countries that lack district-level survey data still have a national
prevalence number from a past survey or from a compilation such as VMNIS, BRINDA,
or DHS anaemia. We can use that single national number to correct the level.

What to do. In the LOCO loop, treat the held-out country's national prevalence as
known (for the demo, compute it as the sample-size-weighted average of
`prev_<outcome>` across that country's districts, standing in for the external
number you would have in practice). Then apply a one-number correction that
shifts all predictions so their population-weighted average equals the national
value. Doing the shift on the log-odds scale, or with isotonic regression, keeps
the district ordering exactly as it was. Report the ranking metrics (which should
not change) and the level metrics (mean absolute error and bias) before and
after.

What success looks like. A clear measure of how much of the transport level error
a single national anchor removes. This addresses Alan's calibration suggestion
directly and tests the project finding that transport mostly fails on level, not
on ranking.

## Task 3: a proper model for proportions, with sampling weights

Why this matters. Prevalence is a proportion bounded between 0 and 1. A plain
linear model can predict impossible values below 0 or above 1 and assumes every
district's estimate is equally noisy. A binomial or quasi-binomial model respects
the 0 to 1 range and the fact that proportions near 0 or 1 vary less than
proportions near one half. On top of that, districts measured with larger samples
deserve more influence.

What to do. Refit the area-level models for a chosen outcome using a binomial or
quasi-binomial family on `prev_<outcome>`, weighting each district by its sample
size `n_<outcome>` or by the precision weight `1 / var_<outcome>`. Compare
calibration and stability against the plain linear version from Task 1. The
`var_<outcome>` column is included precisely so this weighting is easy.

What success looks like. Evidence on whether the proper proportion model gives
better-calibrated and more stable predictions, especially for low-prevalence
outcomes such as vitamin A or B12 in women.

## Task 4: uncertainty on each estimate, and where to survey next

Why this matters. Planners need to know how confident an estimate is, and an
honest uncertainty estimate is itself a useful product. There are two accessible
ways to attach an interval to each district's prediction. A bootstrap resamples
the clusters or districts, refits, and watches how much the predictions move. A
conformal interval uses the model's past prediction errors to build an interval
with a guaranteed long-run coverage rate, without assuming a particular error
distribution.

Once each district has both a predicted prevalence and a measure of uncertainty,
you can rank districts by where a new survey would be most informative. A district
that is both likely high-burden and highly uncertain is a strong candidate,
because measuring it resolves the most doubt about a place that may need action.
This is a value-of-information ranking, and it answers a concrete funder question:
given a fixed budget, where should the next biomarker survey go.

What to do. Produce a per-district interval (start with a cluster or district
bootstrap, or a conformal interval; cross-validated loss-based inference is a
later refinement). Then rank districts by predicted risk multiplied by prediction
uncertainty.

What success looks like. A ranked list of districts that would most improve the
map if surveyed next, demonstrated on at least one outcome.

## Task 5: borrowing strength across deficiencies (sequential modelling)

Why this matters. The deficiencies tend to occur together: a district high in one
is often high in others. If you have already estimated one deficiency well, that
estimate carries information about the others. Modelling the outcomes one at a
time, each conditional on the ones before it, is a standard way to use that shared
information. Statistically this is the forward factorisation of a joint
distribution, the same logic used in g-computation: you build the joint behaviour
of several variables by chaining one-step-ahead models. The wide layout of this
dataset, where each outcome is its own column in the same row, is exactly what
makes the chaining straightforward.

What to do. Order the outcomes from best measured to least measured (iron, then
vitamin A, then B12 and folate, then zinc). Predict the first outcome from the 16
proxies. Then predict the next outcome from the proxies plus the first outcome,
the next from the proxies plus the two before it, and so on. When a country did
not measure an upstream outcome, add a 0/1 indicator marking it missing and
impute a value, so the model can still run. The sequential demonstration inside
`run_superlearner_example.R` shows the mechanics on a single pair.

A caution that matters for transport. Within a country you can condition on the
observed upstream outcome. When you transport to a new country you do not have the
upstream truth, so you must feed in your model's prediction of it. Because that
prediction is itself imperfect, the benefit across countries will be smaller than
within a country, and it can even hurt if the upstream model transports poorly.
For this reason, test the within-country version first and treat the cross-country
version as exploratory.

What success looks like. A comparison, per outcome, of prediction accuracy with
proxies only versus proxies plus the earlier outcomes, reported separately for
within-country and for transport. The clearest expected gains are for the
data-poorest outcomes (zinc, B12), which have the most to borrow.

## Optional connection to your parsimonious DHS work

Your existing project distils the DHS down to roughly four highly informative
proxy variables that recover much of the survey's signal. That is the same goal as
this project's aim of finding cheap, widely available proxies for deficiency. A
natural crossover is to aggregate your four-variable distillation to the district
level and test whether it adds ranking signal on top of, or can stand in for, the
16 geospatial proxies used here.
