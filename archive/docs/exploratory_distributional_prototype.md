# Exploratory prototype: distributional (continuous-biomarker) vs binary modelling

**Script:** `sensitivity/10_distributional_vs_binary_prototype.R`
**Output:** `results/sensitivity/distributional_prototype_metrics.csv`
**Date:** 2026-06-30. Prototype only — small fixed covariate set, single country.

## Question

Does modelling the **continuous** biomarker distribution and integrating to the
deficiency cut-point recover statistical power versus **dichotomising first**?
And how much of any gain is simply about modelling **individuals** instead of
noisy **area proportions**?

## Design

One country (Ghana), three outcome cells spanning a range of prevalence, an
identical a-priori 8-covariate set (soil Zn/Fe/N/P, accessibility, elevation,
population density, built surface), and identical **leave-one-Admin2-out** CV.
All preprocessing (median impute + standardize) is fit on training areas only.

- **S1** — ridge on the survey-weighted district *proportion* (logit). The
  status-quo "predict area-level prevalence from proxies" analog.
- **S1b** — ridge *logistic* on the individual *binary* deficiency outcome.
  Fair binary baseline that uses individual data.
- **S2** — ridge on the individual *continuous* log-biomarker; prevalence
  recovered by integrating the empirical-residual distribution past the cut-point.

Evaluation: Pearson/Spearman/MAE of predicted vs survey-weighted district
prevalence, restricted to areas with survey n ≥ 8. The S2−S1b Pearson contrast
carries a 2000-rep area bootstrap CI.

## Results (leave-one-area-out, Ghana)

| Cell | National prev | S1 (area prop) r | S1b (indiv binary) r | S2 (distributional) r | S2−S1b gain [95% CI] |
|---|---|---|---|---|---|
| child vitamin A | 23.8% | −0.16 | 0.21 | 0.22 | +0.01 [−0.06, 0.09] |
| child iron | 22.0% | 0.48 | 0.48 | 0.52 | +0.04 [−0.04, 0.10] |
| **women iron** | **12.3%** | −0.18 | 0.06 | **0.18** | **+0.12 [0.003, 0.24]** |

(glmnet lambda selection is stochastic; point estimates move ~±0.03 run to run,
but the qualitative pattern is stable.)

## Findings

1. **The dominant lever is individual- vs area-level modelling, not
   continuous-vs-binary.** Regressing the raw weighted *district proportion* on
   covariates (S1) is *actively harmful* at this n — strongly negative for both
   lower-prevalence cells (−0.16, −0.18), because each area's proportion is a
   noisy binomial (some areas n≈6). Any model fit on individuals (S1b or S2)
   recovers a positive signal. **This directly questions the current framing of
   the area-level regression step as primary** for the proxy→prevalence map: the
   area aggregation should happen *after* an individual-level fit, not before it.

2. **The distributional gain is real but conditional — it concentrates in the
   tail.** For moderate-prevalence outcomes (~22–24%, cut-point near the middle
   of the distribution) S2 ≈ S1b. For the low-prevalence outcome (women iron,
   12.3%, cut-point in the lower tail) S2 roughly triples the binary correlation
   (0.06 → 0.18) with a bootstrap CI excluding zero, and modestly lowers MAE.
   This matches theory: the efficiency of continuous-over-binary grows as the
   threshold moves off the median — i.e. exactly the rare/common outcomes where
   the project currently has little to no signal.

## Implication for the analysis

- Refit the proxy→prevalence step at the **individual level** (logistic or
  distributional), aggregating to Admin-2 *after* fitting — do not regress raw
  area proportions. This alone flips several cells from negative to positive.
- Reserve the full **distributional (GAMLSS / location–scale) treatment** for the
  **low- and high-prevalence outcomes**, where it adds the most; near-median
  outcomes get little from it.
- Caveats before promoting beyond prototype: small fixed covariate set (absolute
  r's are below the pipeline's full-feature numbers by design), single country,
  and stochastic lambda (fix or average for a production result).

## Extensions (1a heteroskedastic, 1b transport)

**Script:** `sensitivity/11_distributional_heterosked_transport.R`
**Outputs:** `results/sensitivity/distributional_within_heterosked.csv`,
`results/sensitivity/distributional_transport_childiron.csv`

Added a third strategy **S3** = heteroskedastic distributional model (log-variance
also regressed on covariates, standardized-residual integration), and a
cross-country **leave-one-country-out (LOCO) transport** test on child iron
(Gambia + Ghana + Sierra Leone, biomarker harmonised to the log-ferritin scale,
common `gee_` covariates, cut-point `log 12`).

### 1a — Heteroskedastic σ adds essentially nothing

Within-country Ghana, S3 vs S2 (homoskedastic): child vitA 0.194 vs 0.197,
child iron 0.512 vs 0.512, **women iron 0.174 vs 0.169** (Pearson). The
log-variance layer does not widen the tail advantage. **Conclusion: the tail
gain comes from modelling the continuous response at all, not from area-varying
spread — a full GAMLSS scale model is not worth the complexity; the simple
homoskedastic empirical-residual integration captures the gain.**

### 1b — The distributional model is a within-country tool, not a transport fix

LOCO transport, child iron (Spearman ranking / MAE level on the held-out country):

| Held-out | S1b (binary) Spearman | S2/S3 (distrib) Spearman | S1b MAE | S2 MAE |
|---|---|---|---|---|
| Gambia | 0.519 | 0.53–0.54 | 0.454 | 0.51 |
| Ghana | 0.539 | 0.522 | 0.163 | 0.139 |
| Sierra Leone (n=14) | −0.222 | −0.34 to −0.35 | 0.146 | 0.17 |

- Distributional transports **ranking about as well as binary, within noise** —
  neither is clearly better on any held-out country.
- Its **level (MAE) is sometimes worse** (Gambia 0.51 vs 0.45) because the
  predicted prevalence inherits the training-country ferritin level directly —
  an empirical reproduction of the documented cross-survey ferritin level-offset
  transport failure.
- **Sierra Leone transport fails for every strategy** (negative ranking), and the
  distributional form does not rescue it.

**Conclusion: sketch A (distributional modelling) helps the *within-country*
ranking / value-of-information deliverable for low/high-prevalence outcomes. It
does not improve cross-country transport and can worsen transported level, so it
is not a remedy for the transport problem — that still needs the level-offset
modelling of sketch B and/or national-anchor calibration.**

(Malawi omitted from the transport test — no `gw_child_flag`; child/women share
one biomarker column and must be split on the text `population` field. A clean
add for a follow-up.)
