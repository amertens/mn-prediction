# Literature benchmarking — methods this project should learn from

**Date:** 2026-06-14
**Purpose:** Position the micronutrient prediction pipeline against external work
with the same methodological goals (spatial/survey prediction, small-area
estimation, transportability, decision-focused evaluation) and distil concrete
lessons, cross-referenced to the prioritized recommendations **P1–P8** in
[`docs/CRITICAL_REVIEW.md`](CRITICAL_REVIEW.md).
**Compute note:** desk research only (literature + web verification); no models
were run.

> **Citation policy.** Every reference below was checked against a primary or
> authoritative source (publisher page, arXiv, PubMed/PMC) before inclusion.
> Identifiers (DOI / arXiv id / PMID / PMCID) are given so claims are traceable.
> Where a metadata field could **not** be independently confirmed it is marked
> **⚠** in the verification table at the end — nothing is asserted that was not
> checked.

---

## How this maps to our findings (one-paragraph orientation)

Our review found that (a) within-country cross-validated performance is
optimistic because preprocessing/prescreening is fit outside the resampling
(**P1**); (b) calibration is assessed in-sample (**P2**); (c) intervals likely
under-cover and ignore survey design (**P4**); (d) accuracy is judged against
noisy survey "truth" (**P5**); (e) the SAE benchmarks are comparators rather
than a principled primary estimator (**P7**); and (f) AUC/MSE hide the model's
actual, highly heterogeneous decision value (**P3/P6**). The external literature
below independently validates each of these concerns and supplies established
solutions.

---

## Theme 1 — Honest spatial validation & data leakage

**Supports:** P1 (preprocessing inside resampling; re-report within-country
metrics), and validates our **LOCO design** as the deployment-honest instinct.

- **Roberts et al. (2017)**, *Cross-validation strategies for data with
  temporal, spatial, hierarchical, or phylogenetic structure*, **Ecography**
  40(8):913–929, DOI [10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881).
  The foundational treatment: when data are structured, random CV "regularly
  ignores" that structure and **seriously underestimates predictive error**;
  the fix is *block* CV split along the structure. → Our cluster-blocked CV is
  the right idea but insufficient (spatially adjacent clusters still leak); add
  **spatial-block CV** as a second, more conservative estimate.

- **Ploton et al. (2020)**, *Spatial validation reveals poor predictive
  performance of large-scale ecological mapping models*, **Nature
  Communications** 11:4540, DOI
  [10.1038/s41467-020-18321-y](https://doi.org/10.1038/s41467-020-18321-y).
  The cautionary tale: a non-spatial CV implied the model explained >50% of
  biomass variance; spatial CV revealed **quasi-null** predictive power. → This
  is *exactly* the gap we expect between our (leaky) within-country numbers and
  the honest LOCO numbers, and the reason to lead with LOCO.

- **John, Saurette & Heung (2025)**, *The problematic case of data leakage: a
  case for leave-profile-out cross-validation in 3-dimensional digital soil
  mapping*, **Geoderma** 455,
  [ScienceDirect S0016706125000618](https://www.sciencedirect.com/science/article/pii/S0016706125000618).
  Defines leakage precisely as **overlap between data used for fitting/tuning
  and for testing**, biasing metrics so they are "uninformative regarding the
  model's ability to generalize." → Direct analogue of our **P1** finding: in
  `R/mlr3_fitting.R` imputation, NZV, the supervised `washb_prescreen` (line 307)
  and the correlation/normalize recipe (line 341) are fit on the full data
  **before** folds (line 409). The remedy is structural: fit the entire
  preprocessing graph **inside** each training fold (e.g. an `mlr3pipelines`
  graph or per-fold `recipes::prep`).

- **Brenning & Suesse (2026)**, *Aligning Validation with Deployment in Spatial
  Prediction: Target-Weighted Cross-Validation*, **arXiv:2603.29981**
  ([abs](https://arxiv.org/abs/2603.29981)). Validation is only honest if the
  validation tasks match the *deployment* distribution; TWCV reweights
  validation losses (with spatially-buffered resampling) to the target domain.
  → For us, "deployment" = unsurveyed districts/countries; this is the formal
  justification for weighting evaluation toward the prediction domain rather
  than the surveyed sample.

- **Linnenbrink, Nowosad & Meyer (2026)**, *Moving beyond spatial and random
  cross-validation in environmental modelling: a call for prediction-domain
  adaptive evaluation*, **arXiv:2605.13689**
  ([abs](https://arxiv.org/abs/2605.13689)). Neither random nor spatial CV is
  universally correct; the appropriate scheme depends on whether prediction is
  interpolation or extrapolation relative to the sample. → Argues we should
  report **both** an interpolation-style (cluster/random) and an
  extrapolation-style (spatial-block / LOCO) estimate and label which regime a
  given decision falls in — feeding our out-of-support flags (**P6**).

**Lesson (Theme 1):** within-country CV must be spatial-block (or
target-weighted to the deployment domain), and **all** preprocessing and
supervised screening must live inside the resampling — our **P1** fix exactly.
LOCO is the correct deployment-honest instinct and should be reported as the
headline, with explicit interpolation-vs-extrapolation framing.

---

## Theme 2 — Design-based small-area estimation for survey prevalence mapping

**Supports:** P4 (design-based uncertainty), P5 (sampling-error-aware truth),
P7 (partial-pooling area/cluster model as a *primary* estimator), and the
benchmark suite (`benchmarks_all`).

- **Dong, Wu, Li & Wakefield (2025)**, *Toward a Principled Workflow for
  Prevalence Mapping Using Household Survey Data*, **arXiv:2504.16435**
  ([abs](https://arxiv.org/abs/2504.16435)); also **J. Survey Statistics and
  Methodology** (2025), DOI
  [10.1093/jssam/smaf048](https://doi.org/10.1093/jssam/smaf048). A concrete,
  reproducible workflow that **acknowledges the sampling design explicitly** and
  models spatial dependence, with model-choice and interpretation guidance for
  data-sparse LMIC settings. → A ready-made template for restructuring our SAE
  path around design-based estimation (**P4/P7**).

- **[Authors per arXiv] (2026)**, *Mapping Subnational Vulnerability to
  Inadequate Micronutrient Intake using a Bayesian Small Area Estimation
  Framework*, **arXiv:2604.14971** ([abs](https://arxiv.org/abs/2604.14971)).
  Most on-point comparator: applies Bayesian SAE to **micronutrient** inadequacy
  at admin-2 from household surveys; in simulation a **cluster-level
  Beta-binomial** model performed best overall, and an **area-level joint
  smoothing** model was the best among design-aware options. → Strong evidence
  that a partial-pooling cluster/area model should be our **primary** estimator
  (**P7**), not merely a benchmark, and points to specific model forms to adopt.

- **Wakefield (2020)**, *Small Area Estimation for Disease Prevalence Mapping*,
  **International Statistical Review**, DOI
  [10.1111/insr.12400](https://doi.org/10.1111/insr.12400) (PMID 36081593;
  PMC9451141). Review of design-based vs model-based, area-level vs unit-level
  SAE, illustrated on **Malawi 2015–16 DHS** (the same survey family we use). →
  Clarifies the area-vs-unit choice that our pipeline currently leaves implicit
  (individual-SL-aggregated vs area-level) and that we flagged cannot be cleanly
  compared today.

- **Mercer, Wakefield, Pantazis, Lutambi, Masanja & Clark (2015)**,
  *Space–Time Smoothing of Complex Survey Data: Small Area Estimation for Child
  Mortality*, **Annals of Applied Statistics** 9(4):1889–1905, DOI
  [10.1214/15-AOAS872](https://doi.org/10.1214/15-AOAS872). The canonical
  design-based SAE pattern: compute **design-based direct (weighted) estimates +
  design variances** per area, then smooth them with a spatial model. → This is
  the principled replacement for (i) our naïve weighted means and (ii) the
  conformal delta-method SE that currently ignores strata/PSU (**P4**); the
  `SUMMER` R package implements it.

- **[LSHTM group] (2025)**, *Predicting risk of inadequate micronutrient intake
  with transferable machine learning models*, **Scientific Reports**, DOI
  [10.1038/s41598-025-26179-7](https://doi.org/10.1038/s41598-025-26179-7).
  Demonstrates **cross-country model transfer** (Ethiopia↔Nigeria) for
  micronutrient risk using routinely available predictors — the same goal as our
  LOCO/OOS work — and reports predictor consistency across countries. → External
  precedent that our transport ambition is sound; useful comparator for honest
  framing of transfer performance.

**Lesson (Theme 2):** make a **design-aware partial-pooling** area/cluster-level
model a *primary estimator* (Beta-binomial cluster model or area-level smoothing
of design-based direct estimates), and let the survey design enter **both** the
point estimate and its uncertainty — our **P4/P7**, with **P5** following
naturally because design variances quantify the noise in the "truth."

---

## Theme 3 — Calibrated uncertainty & the gold-standard geostatistical template

**Supports:** P2 (out-of-sample calibration), P4 (calibrated, design-based
intervals), and the dashboard's currently-stubbed GBD anchoring.

- **Ferreira LZ et al. (2020)**, *Geospatial estimation of reproductive,
  maternal, newborn and child health indicators: a systematic review of
  methodological aspects of studies based on household surveys*, **International
  Journal of Health Geographics** 19:41, DOI
  [10.1186/s12942-020-00239-9](https://doi.org/10.1186/s12942-020-00239-9)
  (PMID 33050935; PMC7552506). Across 82 studies, **out-of-sample validation was
  reported in only ~56%**, and **>25% reported no uncertainty measure** at all.
  → Doing honest OOS validation (P1) *and* calibrated, design-based uncertainty
  (P4) would put this project **ahead of the published field** — a framing worth
  stating explicitly in the manuscript.

- **Burstein et al. (2019)**, *Mapping 123 million neonatal, infant and child
  deaths between 2000 and 2017*, **Nature** 574(7778):353–358, DOI
  [10.1038/s41586-019-1545-0](https://doi.org/10.1038/s41586-019-1545-0). The
  IHME Local Burden of Disease template: Bayesian model-based geostatistics
  producing **posterior surfaces** that are **aggregated to admin units by
  summarizing posterior draws**, and **calibrated/anchored to GBD** national
  totals. → The reference design for (i) propagating predictive uncertainty to
  admin-2 by draw-aggregation rather than broadcasting an Admin-1 CI downward
  (our **P4** anti-pattern), and (ii) the **GBD anchoring** the dashboard's GBD
  comparison tab currently stubs.

**Lesson (Theme 3):** move from a single point + broadcast interval to
**posterior/predictive draws aggregated to admin units**, report calibrated
design-based intervals, validate calibration **out-of-sample** (P2), and
consider **GBD-anchoring** national totals (replacing the placeholder).

---

## Theme 4 — Decision-focused evaluation for stakeholders

**Supports:** Part D and our new **Decision value** dashboard tab (P3), and the
trust/fitness flags (P6).

- **Vickers & Elkin (2006)**, *Decision Curve Analysis: A Novel Method for
  Evaluating Prediction Models*, **Medical Decision Making** 26(6):565–574, DOI
  [10.1177/0272989X06295361](https://doi.org/10.1177/0272989X06295361)
  (PMID 17099194). Net benefit across action thresholds, compared to
  **treat-all / treat-none** — accuracy metrics "do not address clinical
  consequences." → The formal backbone for upgrading our targeting/enrichment
  panel into a **net-benefit / decision-curve** view at an explicit action
  threshold.

- **Singh, Shah & Vickers (2023)**, *Assessing the net benefit of machine
  learning models in the presence of resource constraints*, **JAMIA**
  30(4):668–673 (PMID 36810659; PMC10018264). Extends decision curves to the
  **resource-constrained** case (you can only act in *k* units) — precisely the
  ministry's "fund the top-N districts" problem our targeting panel models. →
  Gives a principled threshold/capacity formulation for the Decision-value tab.

- ***Evaluation of performance measures in predictive AI models to support
  medical decisions: overview and guidance*** (2025), **The Lancet Digital
  Health**,
  [PIIS2589-7500(25)00098-6](https://www.thelancet.com/journals/landig/article/PIIS2589-7500(25)00098-6/fulltext).
  Of five performance domains (discrimination, calibration, overall,
  classification, **clinical utility**), **clinical utility — does the model
  lead to better decisions than not using it — is the decision-relevant one**,
  with net benefit the standard measure. → Authoritative support for our core
  thesis that AUC/MSE are insufficient and decision utility must lead (P3).

- **Griffin, Claxton & Sculpher (2008)**, *Decision analysis for resource
  allocation in health care*, **Journal of Health Services Research & Policy**
  13(Suppl 3):23–30, DOI
  [10.1258/jhsrp.2008.008017](https://doi.org/10.1258/jhsrp.2008.008017).
  Resource-allocation decision analysis and **value of information**: the value
  of acting on imperfect evidence vs the status quo, and the expected cost of
  the wrong allocation. → Framework for a **value-of-information** framing in the
  dashboard: expected gain from acting on the model vs the national-average
  status quo, and expected cost of misallocation where the model is unreliable
  (ties to our fitness-for-purpose scorecard, P6).

- ***Qualitative evaluation of the use of modelling in resource allocation
  decisions for HIV and TB*** (2024), PMC10806894 (medRxiv preprint
  [2023.04.11.23288405](https://www.medrxiv.org/content/10.1101/2023.04.11.23288405));
  evaluates the Optima allocative-efficiency tool. Headline: uptake hinges on
  **stakeholder co-design** and timing analyses to real budget/strategy
  decisions. → Methodological correctness is necessary but not sufficient; the
  Decision-value tab should be **co-designed with ministry users** and tied to
  actual planning cycles.

**Lesson (Theme 4):** formalize the decision tab with **net-benefit /
decision-curve analysis at the action threshold** (resource-constrained
variant), add a **value-of-information** framing (expected value of acting on
the model vs status quo; expected cost of misallocation), and **co-design** with
stakeholders for uptake.

---

## What to adopt — mapped to P1–P8

| Our rec | Adopt from the literature | Concrete step |
|---|---|---|
| **P1** preprocessing inside CV; re-report within-country metrics | Roberts 2017; Ploton 2020; John–Saurette–Heung 2025; Brenning–Suesse 2026; Linnenbrink 2026 | Wrap impute/NZV/**prescreen**/corr/normalize in an `mlr3pipelines` graph cross-validated per fold; add **spatial-block CV**; report interpolation *and* extrapolation estimates. |
| **P2** out-of-fold calibration | Lancet Digit Health 2025 (calibration domain); Ferreira 2020 | Fit/evaluate Platt on disjoint folds (`calibrate_binary_oof`); re-issue `diagnostics_binary_calibrated.csv`; assess calibration only OOS. |
| **P3** lead with decision value | Vickers–Elkin 2006; Singh–Shah–Vickers 2023; Lancet Digit Health 2025 | Promote net-benefit/decision-curve as the headline metric in dashboard + manuscript; keep AUC/MSE secondary. |
| **P4** design-based, calibrated intervals; stop Admin-1→Admin-2 broadcast | Mercer 2015; Wakefield 2020; Dong–Wakefield 2025; Burstein 2019 | Compute design-based direct estimates + design variances (`srvyr`/`SUMMER`), smooth spatially, and **aggregate posterior draws** to admin units. |
| **P5** sampling-error-aware "truth" | Mercer 2015; Wakefield 2020; arXiv 2604.14971 | Use design variances to denoise error metrics (or inverse-variance weight / stratify by `n_svy`). |
| **P6** out-of-support / trust flags | Brenning–Suesse 2026; Linnenbrink 2026; Griffin 2008 (VoI) | Flag extrapolation/low-support districts; show expected misallocation cost where unreliable. |
| **P7** partial-pooling SAE as **primary** estimator | arXiv 2604.14971 (Beta-binomial / area smoothing); Dong–Wakefield 2025; Wakefield 2020 | Make a design-aware cluster/area model the primary estimator; demote the ML stack to one comparator in one harness. |
| **P8** fix latent leakage; key joins on (country, Admin2) | John–Saurette–Heung 2025 (leakage definition) | Center transport features by training-pooled means only; persist `X_train`; join on `(country, Admin2)` throughout. |

---

## Headline lessons
1. **Honest validation first.** The field's own cautionary tales (Ploton 2020)
   show non-spatial CV can invert conclusions; our P1 leakage must be fixed and
   LOCO led, with spatial-block CV added.
2. **Design-aware partial-pooling SAE should be the primary estimator**
   (arXiv 2604.14971; Wakefield 2020; Mercer 2015) — not a benchmark — with the
   survey design in both estimate and uncertainty.
3. **Propagate draws, calibrate out-of-sample, consider GBD anchoring**
   (Burstein 2019); doing OOS validation + uncertainty already beats ~half the
   published RMNCH literature (Ferreira 2020).
4. **Decision utility leads, not AUC/MSE** (Lancet Digit Health 2025;
   Vickers–Elkin 2006; Singh 2023): formalize the Decision-value tab with
   net-benefit at the action threshold and a value-of-information framing, and
   co-design it with stakeholders (PMC10806894).

---

## Citation verification log

| Ref | Identifier | Verified | Notes / flags |
|---|---|---|---|
| Roberts et al. 2017, Ecography 40(8):913–929 | DOI 10.1111/ecog.02881 | ✅ publisher + Semantic Scholar | — |
| Ploton et al. 2020, Nat Commun 11:4540 | DOI 10.1038/s41467-020-18321-y | ✅ nature.com | — |
| John, Saurette & Heung 2025, Geoderma 455 | PII S0016706125000618 | ✅ ScienceDirect + search | Numeric DOI not separately captured (paywall 403); cite via PII. |
| Brenning & Suesse 2026 | arXiv:2603.29981 | ✅ arXiv | User's short title corrected to full title. |
| Linnenbrink, Nowosad & Meyer 2026 | arXiv:2605.13689 | ✅ arXiv | — |
| Dong, Wu, Li & Wakefield 2025 | arXiv:2504.16435; DOI 10.1093/jssam/smaf048 | ✅ arXiv + OUP | Also published in JSSAM. |
| Micronutrient Bayesian SAE 2026 | arXiv:2604.14971 | ✅ arXiv | **⚠ author names not captured** (cite by arXiv id). |
| Wakefield 2020, Int Stat Rev | DOI 10.1111/insr.12400; PMC9451141; PMID 36081593 | ✅ Wiley + PMC | — |
| Mercer et al. 2015, Ann Appl Stat 9(4):1889–1905 | DOI 10.1214/15-AOAS872 | ✅ Project Euclid + PubMed | — |
| Transferable ML micronutrient 2025, Sci Rep | DOI 10.1038/s41598-025-26179-7 | ✅ nature.com | **⚠ exact author list not captured** (LSHTM group; biorxiv 2025.04.08.647715). |
| Ferreira LZ et al. 2020, IJHG 19:41 | DOI 10.1186/s12942-020-00239-9; PMC7552506; PMID 33050935 | ✅ PMC | First author Ferreira, L.Z.; many co-authors. |
| Burstein et al. 2019, Nature 574:353–358 | DOI 10.1038/s41586-019-1545-0 | ✅ nature.com | — |
| Vickers & Elkin 2006, Med Decis Making 26(6):565–574 | DOI 10.1177/0272989X06295361; PMID 17099194 | ✅ SAGE + Scholar | — |
| Singh, Shah & Vickers 2023, JAMIA 30(4):668–673 | PMID 36810659; PMC10018264 | ✅ OUP + PMC | **⚠ numeric DOI not separately captured**; PMID/PMC confirmed. |
| LDH performance-measures guidance 2025 | PII S2589-7500(25)00098-6 | ✅ thelancet.com | **⚠ author list not captured**; title/venue confirmed. |
| Griffin, Claxton & Sculpher 2008, J Health Serv Res Policy 13(S3):23–30 | DOI 10.1258/jhsrp.2008.008017 | ✅ SAGE | Distinct from their "Dangerous omissions" paper. |
| Optima HIV/TB modelling evaluation 2024 | PMC10806894; medRxiv 2023.04.11.23288405 | ✅ PMC + medRxiv | **⚠ journal venue not independently confirmed** beyond PMC. |

*All 17 references had their identity (title + venue/repository + identifier)
confirmed against a primary source. Items marked ⚠ have a single secondary
metadata field (author list, numeric DOI, or venue) that could not be confirmed
without paywalled access; their identifiers are nonetheless verified.*
