# Subnational micronutrient deficiency prediction from proxy data: project overview for external review (June 2026)

**Purpose.** This document is a self-contained overview of the project, written to
solicit feedback from external experts. It assumes no prior familiarity with our
goals or approach, but it is written for readers who are knowledgeable about
micronutrient deficiency assessment, geospatial / small-area modelling, and the
workings of the IHME Global Burden of Disease (GBD) enterprise. It describes the
data and data sources, the data-cleaning and harmonisation decisions, the methods
we have tried and what we have found scientifically, our relationship to GBD, and
the ideas we are weighing for extending the work. A list of specific questions for
reviewers is at the end (Section 9).

Deeper technical detail lives in the manuscript and its supplementary appendix
(`manuscript_mcn.qmd`, `manuscript_mcn_appendix.qmd`), the full data lineage
(`data_sources_appendix.qmd`), and the exploratory findings logs
(`exploratory_*.md`); this overview summarises and connects them.

---

## 1. Executive summary

Micronutrient deficiencies (iron, vitamin A, zinc, folate, vitamin B12) are
measured by blood biomarkers, which are expensive and rarely collected. Many
countries have at most one national micronutrient survey in the last decade, and
those surveys do not measure every nutrient. The goal of this project is to test
whether **routinely available proxy data** (household-survey indicators, satellite
and climate layers, modelled disease burden, food security and prices, soil
chemistry) can predict **subnational (district-level) micronutrient deficiency**,
both in countries with outdated data and, by transport, in countries with none.

We assembled individual-level biomarker data from four sub-Saharan African surveys
(The Gambia, Ghana, Sierra Leone, Malawi) and linked them to roughly 200
harmonised household-survey proxies and more than 2,500 earth-observation
covariates. We compared machine-learning ensembles (SuperLearner) against
canonical small-area-estimation (SAE) methods, evaluated honestly with
leakage-free, spatially blocked cross-validation within countries and
leave-one-country-out (LOCO) holdout across them.

Headline findings, stated honestly:
- The **binding constraint is the data, not the model**. The proxy predictors are
  almost all resolved at the district level, so the effective sample size is the
  number of districts (14 to 87 per country), not the thousands of individuals.
- Under honest evaluation, **within-country district ranking is modest** and
  **cross-country transport is weak**, failing mostly on the absolute *level*
  (predictions are pulled toward the training-country mean) rather than on
  ranking. **National prevalence estimates, by contrast, are accurate** (within
  about 5 percentage points).
- A multi-learner **SuperLearner does not reliably beat simple SAE or a penalised
  regression**; the best transportable model we found is a thin-plate spatial
  spline on district centroids plus a few soil-micronutrient features.
- We have therefore **reframed the deliverable** toward what the data can support:
  correctly flagging the highest-burden districts, honest uncertainty, and
  value-of-information (where to survey next), rather than precise absolute
  prevalence everywhere.

We want external feedback on whether this framing is right, on the
harmonisation and modelling choices, and on how best to position the work
relative to GBD.

## 2. Background and goal

Two methodological traditions try to fill the micronutrient data gap. The
**small-area-estimation** tradition (Fay-Herriot, Bayesian disease mapping such as
BYM2) builds structured spatial smoothers from survey-derived area prevalence and
known sampling variances. The **machine-learning** tradition trains flexible
models on satellite-derived covariates to map indirect indicators of nutritional
status. These two have rarely been compared head to head on the same biomarker
data at the same geographic scale. We do that, and we ask two questions:

1. **Within-country**: trained on a country's own data, how well do proxy models
   reproduce district-level deficiency?
2. **Across countries (transport)**: trained on some countries, do predictions
   transport to a held-out country using only proxies that exist everywhere?

The project has two aims: a **variable-importance** aim (which proxy domains carry
signal) and a **parsimonious, transportable model** aim (a deployable model for
countries without biomarker data).

## 3. Data and data sources

### 3.1 Biomarker (outcome) surveys

Four nationally representative micronutrient or micronutrient-containing biomarker
surveys, one per country, each a different year with a partly different panel:

- The Gambia, 2018; Ghana, 2017; Sierra Leone, 2013; Malawi, 2015 to 2016.
- Biomarkers were extracted for women of reproductive age (15 to 49 y) and
  children 6 to 59 months where measured; one to two surveys also include adult
  men. Sample sizes are roughly 750 to 1,400 individuals per outcome per country.
- Each survey is georeferenced at the survey-cluster (primary sampling unit) level
  and carries survey design information (weights, PSUs, strata).
- Iron is measured by serum ferritin in all four countries; vitamin A by
  retinol-binding protein (a proxy for serum retinol); vitamin B12, folate and
  zinc by serum concentration in a subset of countries. The cross-country
  comparability of these biomarkers is itself debated, which is one reason we
  treat absolute levels cautiously.

Outcomes available in all four countries: iron and vitamin A (women and
children). B12 and folate: Ghana, Malawi, Sierra Leone. Zinc: Malawi only.

### 3.2 Predictor domains

Predictors come from eleven domains. With one exception (the survey biomarker
matrix, deliberately excluded) all are public or routinely available, so a fitted
model can be applied where no biomarker survey exists.

1. **Survey biomarker matrix** (`gw_`): the original survey items, **excluded**
   from the deployable proxy-only models; retained only for benchmark models that
   deliberately use full survey information.
2. **DHS** (`dhs<year>_`): household and child indicators aggregated to the
   administrative unit (education, assets, water/sanitation, child health,
   malaria interventions). Stamped with each survey's year, so DHS columns do not
   share names across countries.
3. **MICS** (`mics_`): additional household/child aggregates where available.
4. **IHME GBD geospatial estimates** (`ihme_`): modelled related conditions
   (stunting, wasting, diarrhoea, anaemia, under-five mortality) at admin-1/2.
5. **LSMS** (`lsms_`): World Bank living-standards aggregates; wired for Ghana.
6. **Malaria Atlas Project** (`MAP_`): falciparum/vivax incidence, parasite rate,
   mortality, reproductive number; intervention coverage (ITN, IRS, effective
   antimalarial treatment); and inherited blood disorders (HbS, HbC, G6PD, Duffy)
   that confound ferritin-based iron assessment.
7. **Google Earth Engine rasters** (`gee_`): (a) cluster-buffer summaries at
   10/25/50 km around each cluster GPS (monthly precipitation, NDVI, land-surface
   temperature; plus population density, built-up surface, settlement model,
   global human modification, grassland fraction); and (b) admin-2 polygon zonal
   means for the full raster set including soil chemistry, climate, vegetation and
   agriculture, settlement, accessibility, demography, land cover, livestock.
8. **WFP food prices** (`wfp_`).
9. **FluNet** influenza surveillance (`flunet_`).
10. **Food security** (`fsec_`): HungerMap / HFID household food-insecurity and
    IPC / Cadre Harmonise phase classifications.
11. **SoilGrids / ISRIC** (`soil_`, `gee_a2_Soil*`): remotely sensed soil
    micronutrient and chemistry layers (zinc, iron, calcium, magnesium,
    phosphorus, potassium, sulphur, nitrogen, carbon, CEC, aluminium).

Two predictor sets are built: an **individual-level set** of roughly 200
harmonised household proxies shared across the DHS surveys, and an **area-level
set** of more than 2,500 earth-observation aggregates at admin-2; for
cross-country work the latter is restricted to the 152 variables present in at
least three countries (the common-vocabulary set).

### 3.3 Geographic and temporal scope

Administrative units are GADM admin-1 (region) and admin-2 (district); admin-3
(chiefdom) for Sierra Leone in a resolution sensitivity analysis. The pooled
admin-2 dataset is about 206 districts (Gambia 30, Ghana 75, Sierra Leone 14,
Malawi 87). The proxy predictors exist for all of sub-Saharan Africa, so the
fitted models can in principle be applied beyond the four survey countries.

## 4. Data cleaning, linkage, and harmonisation

### 4.1 Outcome definitions and inflammation adjustment

Each continuous biomarker is converted to a binary deficiency indicator using WHO
or BRINDA-recommended thresholds; harmonising these across surveys is a central
step.

- **Vitamin A**: retinol-binding protein below 0.70 µmol/L (women and children).
- **Iron**: serum ferritin below 15 µg/L in women, below 12 µg/L in children. The
  Gambia uses the equivalent cut-points on the log-ferritin scale; Ghana uses
  Thurnham-adjusted ferritin; Sierra Leone and Malawi use the standard cut-points.
- **Folate**: serum folate below 10 nmol/L (Malawi's native rule, below 3 ng/mL,
  is harmonised to this for cross-country comparability).
- **Vitamin B12**: serum B12 below 148 pmol/L.
- **Zinc**: serum zinc below WHO age/sex/fasting-specific cut-points (Malawi).

Ferritin (and to a lesser degree retinol-binding protein) rises with
inflammation, so where C-reactive protein and alpha-1-acid glycoprotein were
measured, ferritin is corrected (BRINDA regression or Thurnham correction factors,
depending on the country's original processing) before thresholding. A practical
consequence is that the iron-deficiency definition is **not perfectly identical
across countries**, which contributes to between-country prevalence differences
and to the unreliability of transporting absolute levels.

### 4.2 Predictor linkage

- **Cluster-level linkage**: earth-observation cluster-buffer summaries and
  Malaria Atlas rasters are extracted at each cluster's GPS point and joined on
  cluster ID; these predictors are therefore constant within a cluster.
- **Area-level linkage**: sources available only as area aggregates (DHS, MICS,
  IHME, LSMS, food security, GEE admin-2 zonal means) are joined by administrative
  unit. Because survey and external admin names are spelled inconsistently, names
  are matched with fuzzy string matching (Jaro-Winkler at admin-1, string-distance
  methods at admin-2) with conservative thresholds; unmatched units are left
  missing rather than guessed. **Every survey-to-area join is keyed on the pair
  (country, admin-2), never on admin-2 name alone**, to prevent cross-country
  mismatches.

### 4.3 Temporal alignment

For multi-version sources, the version closest to each country's survey year is
selected; climate uses the survey year directly, slowly changing layers (soil)
use the most recent global version, and multi-band annual rasters are collapsed to
the survey-year annual mean.

### 4.4 Cleaning conventions

All `gw_` survey items are removed from the proxy pool (the proxy-only rule);
high-cardinality categoricals are dropped or collapsed; geometry and raw temporal
columns are removed from the feature matrix after they have been used for linkage;
missing values are retained as NA in the merged table and imputed downstream
inside model fitting (strictly within each training fold in the corrected
methods). Survey design is handled with weights and a design-effect assumption of
1.5 where the cluster structure does not support a direct estimate.

### 4.5 The effective-sample-size problem (important context)

Almost every predictor is resolved at the admin-2 level (DHS, MICS, IHME, LSMS,
food security, soil, and the admin-2 GEE means), and even the cluster-buffer GEE
layers are nearly constant within a district. As a result **close to 100% of
predictor variance is between-area**, and the effective sample size for modelling
is the number of districts (14 to 87 per country), not the underlying individuals.
This single fact explains much of what follows: fitting hundreds of covariates
against tens of independent points is heavily over-parameterised, the apparent
precision from a thousand individuals is pseudo-replication, and transfer is
fragile. It is also why we treat the area as the unit of analysis.

## 5. Methods

### 5.1 Estimation strategies compared

We compared strategies from both traditions, in two families. The **area-level
SAE family** (the primary analysis) fits one row per admin-2 polygon with
survey-weighted prevalence as the target: an area-level SuperLearner on the
polygon-aggregated covariates; frequentist baselines (a country-mean floor, OLS,
cross-validated elastic-net); Bayesian areal models (Fay-Herriot, and BYM2 fit by
integrated nested Laplace approximation); and transport-aware comparators (a
country random-intercept model, a two-stage country-mean-plus-deviation model, a
forward selection that optimises a leave-one-training-country-out criterion, a
thin-plate spatial-spline GAM, a data-source group-lasso, and a stacked ensemble).
The **individual-level SuperLearner aggregated to admin-2** (a sensitivity
analysis) fits a person-level ensemble (mean, penalised GLMs, random forest,
gradient-boosted trees, dropout BART, with class-weighted learners for rare
outcomes) with cluster-blocked folds, then averages predictions to the district.

### 5.2 Honest-evaluation ("corrected") methods, the analysis of record

Initial benchmarking showed that an analysis at this scale can flatter itself, so
our headline results come from a set of leakage-free corrected methods that are
the analysis of record (the original pipeline is kept as a sensitivity
comparison): all preprocessing inside each training fold; folds blocked by survey
cluster and by spatial region; out-of-fold (not in-sample) recalibration;
separation of survey-sampling error from model error in admin-2 error;
split-conformal and design-based prediction intervals; decision-value scoring of
district ranking; out-of-support trust flags; and a partial-pooling
(empirical-Bayes) area estimator. All joins and standardisation are
training-only by construction.

### 5.3 Evaluation design

Within-country accuracy uses 5-fold (and cluster/spatial-block) cross-validation
with metrics on out-of-fold predictions only. Cross-country transport uses
leave-one-country-out holdout. Metrics: Pearson and Spearman correlation (spatial
ranking), MAE and RMSE in percentage points and bias (absolute accuracy), and
calibration slope; for binary classifiers, ROC-AUC, Brier score and skill,
expected calibration error, and calibration intercept/slope. We also report the
sampling-adjusted RMSE, split-conformal coverage, and decision-value measures.

## 6. Scientific findings so far

1. **The binding constraint is the data, not the model** (four surveys, ~206
   districts, individually weak proxies, disputed biomarker comparability).
2. **Effective sample size is the number of areas (14 to 87)**, not individuals;
   ~100% of predictor variance is between-area. This drives weak, unstable
   transfer.
3. **Honest evaluation lowers performance.** Replacing random folds and full-data
   preprocessing with cluster/spatial-block folds and strict in-fold preprocessing
   moves within-country discrimination toward chance (for example Gambia child
   iron AUC about 0.57 down to 0.48), and the directly-fit area ensemble's
   apparent excellence was an in-sample artifact.
4. **Within-country ranking is modest; national estimates are good.** The
   individual-level SuperLearner aggregated to admin-2 ranks districts at a median
   Pearson r near 0.53 (the area ensemble near 0.12 under honest evaluation), but
   national prevalence is accurate to within about 5 percentage points.
5. **Transport fails mainly on level, not ranking.** Predictions are pulled toward
   the training-country mean; iron transports best, vitamin A essentially does
   not, and Malawi (the only non-West-African country) anti-correlates when held
   out. A one-parameter calibration to a known national prevalence removes most of
   the level error (mean MAE about 20.9 down to 9.6 pp in our prototype) without
   changing ranking.
6. **The ensemble does not beat simple methods.** A plain elastic-net, OLS, or
   Fay-Herriot matches or beats the SuperLearner on most LOCO splits. The best
   transportable model is the simplest: a thin-plate spatial spline on district
   centroids plus a few soil-micronutrient features (mean LOCO r about 0.33),
   beating all 152 earth-observation covariates. Soil micronutrient content is the
   closest covariate we have to a causal pathway (soil to crop to intake to
   biomarker).
7. **Calibration and uncertainty are now honest.** Several binary models had
   negative Brier skill when raw, repaired by out-of-fold recalibration;
   split-conformal intervals cover roughly 87 to 100 percent against a 90 percent
   target; honest area intervals remain wide (Fay-Herriot tightens them 19 to 33
   percent).
8. **Preprocessing lessons** (from a dedicated feature-engineering investigation):
   rank/quantile normalisation beats z-scoring for the linear learners; aggressive
   dimensionality reduction matches the full feature set within country; median
   imputation suffices; unsupervised PCA loadings are not transportable (components
   flip sign across countries); and predictor domains should be selected by outcome
   biology (vitamin A transfers on environment, iron on malaria and modelled
   health burden) rather than pooled.

The overall reframe: because **ranking survives where absolute level does not**,
the defensible deliverables are (a) flagging the highest-burden districts and (b)
value-of-information sampling, plus honest national estimates anchored to known
totals, rather than precise absolute prevalence everywhere.

## 7. Relationship to GBD / IHME and the complementary data regime

### 7.1 Per-nutrient covariate and method mapping (GBD vs this project)

| Nutrient | GBD quantity and case definition | GBD data and model | GBD key covariates / drivers | This project | Where we depart / the opportunity |
|---|---|---|---|---|---|
| **Vitamin A** | Deficiency prevalence (serum retinol < 0.70 µmol/L) and vision-loss sequelae | WHO VMNIS + DHS + surveys, plus modelled vitamin A supplementation coverage; GBD 2017 DisMod-MR 2.1, GBD 2019+ spatio-temporal Gaussian process regression (ST-GPR) | Vitamin A supplementation coverage (a strong driver), location-level stunting, socio-demographic index / income | Measured RBP-based deficiency, women and children, 4 countries, Admin-2, proxy ML + SAE | We use a measured biomarker at district resolution; GBD leans on supplementation-coverage covariates, the suspected driver of its modelled decline. We can test whether biomarker VAD tracks that decline. |
| **Iron** | Iron-deficiency anaemia, via the anaemia envelope (not a direct ferritin model) | Haemoglobin distributions (DHS and other Hb surveys) give total anaemia by severity; anaemia split to causes via cause-specific haemoglobin shifts, with a residual (~10%) assigned mostly to dietary iron deficiency | Haemoglobin survey data, cause-specific Hb shifts, malaria and haemoglobinopathy prevalence, SDI | Measured serum ferritin (inflammation-adjusted) deficiency, Admin-2 | We estimate iron deficiency directly from the exposure (ferritin), not by attributing anaemia. Our biomarker fractions could inform their iron-attribution step. |
| **Zinc** | Risk-factor exposure = prevalence of inadequate dietary zinc intake (NOT a biomarker) | FAO food-balance-sheet zinc availability vs IZiNCG physiological requirements (Miller absorption equation); diet surveys; ST-GPR | Lag-distributed income, dietary energy availability; historically correlated with stunting | Measured serum zinc (Malawi only), Admin-2 | A fundamentally different signal: measured biomarker vs diet-supply proxy. Where they diverge is scientifically important; we can flag where diet-supply mis-states biomarker status. |
| **Folate** | Not a standard deficiency-prevalence cause; enters mainly via neural-tube-defect burden | No subnational biomarker prevalence product | n/a | Measured serum folate, women, 3 countries, Admin-2 | We produce a biomarker deficiency map GBD does not. |
| **Vitamin B12** | Not a standard GBD cause or risk | n/a | n/a | Measured serum B12, women, 3 countries | Fills a gap GBD does not cover. |
| **Iodine** | Iodine deficiency modelled (urinary iodine, goiter) | Survey-based, DisMod / ST-GPR | Salt iodisation coverage | One country only; not a focus | GBD is stronger here; not a comparative priority for us. |

### 7.2 How IHME GBD estimates micronutrient deficiencies

GBD treats micronutrient problems in two different machines. **Vitamin A and
iodine deficiency** are modelled as *causes* (deficiency prevalence plus
sequelae). **Zinc deficiency** is modelled as a *dietary risk factor* (inadequate
intake), and **iron** enters through the **anaemia** cause-attribution machinery,
not as a direct serum-ferritin model. **Folate and B12** are not standard
stand-alone GBD deficiency-prevalence products.

The estimation engines are **DisMod-MR 2.1** (a Bayesian compartmental
meta-regression enforcing consistency among incidence, prevalence, remission and
mortality) and increasingly **ST-GPR** (spatio-temporal Gaussian process
regression) for exposure and prevalence. Both ingest sparse, heterogeneous data
and produce a complete cube of estimates for every location, year, age and sex,
by **borrowing strength across a geographic hierarchy** (global to super-region to
region to country to, where modelled, subnational) and **across time**, with
study-level covariates (adjusting for measurement and case-definition differences)
and country-level covariates (socio-demographic index, lag-distributed income,
maternal education, healthcare access and quality, and cause-specific drivers).

Per nutrient:
- **Vitamin A**: serum retinol below 0.70 µmol/L, from WHO VMNIS, DHS and other
  surveys, with modelled vitamin A supplementation coverage as a covariate; GBD
  2017 used DisMod-MR 2.1, GBD 2019 moved to ST-GPR and added stunting. Reported
  burden fell sharply between rounds, driven largely by methods changes (a switch
  to MR-BRT with 10% trimming, reduced relative risks, removed outcomes), and the
  authors caution the change may be partly artifact
  ([basis for GBD 2017/2019 vitamin A and zinc changes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9991746/)).
- **Iron**: total anaemia prevalence by severity from haemoglobin distributions
  (DHS and other Hb surveys) is split to causes using cause-specific haemoglobin
  shifts; a residual (~10%) is assigned mostly to dietary iron deficiency
  ([GBD anaemia methods](https://ashpublications.org/blood/article/123/5/615/32839/A-systematic-analysis-of-global-anemia-burden-from);
  [GBD 2021 anaemia](https://www.thelancet.com/journals/lanhae/article/PIIS2352-3026(23)00160-6/fulltext)).
- **Zinc**: not a biomarker; the prevalence of inadequate dietary zinc intake from
  FAO food supply versus IZiNCG requirements
  ([Wessells and Brown 2012](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0050568)),
  modelled with ST-GPR.
- **Iodine** from urinary iodine and goiter with salt iodisation coverage; **folate
  and B12** have no standard subnational biomarker prevalence product.

### 7.3 How our data and methods could augment, improve, or provide alternatives

We are a different data regime: cross-sectional measured biomarkers in four
countries plus universal open proxies, borrowing strength **across space**, versus
GBD's longitudinal multi-source modelling that borrows **across time and
neighbouring locations**. That difference is the source of both our limitations
and our distinct value:
1. **Independent measured-biomarker check** on exposure-driven estimates,
   especially vitamin A (whose modelled decline leans on supplementation
   covariates) and the diet-supply-based zinc estimate.
2. **An alternative serum-zinc signal** versus GBD's food-supply proxy.
3. **Direct ferritin-based iron deficiency** that could constrain GBD's
   iron-attributable fraction of anaemia.
4. **Filling gaps** GBD does not model (subnational folate and B12).
5. **Admin-2 resolution** for within-country targeting.
6. **A different, contributable covariate space** (soil micronutrients, climate,
   vegetation, malaria) that could be offered as candidate covariates for GBD's
   framework.
7. **Two-way benchmarking**: our national-anchor calibration makes it easy to
   reconcile our subnational maps with GBD or VMNIS national totals, and to
   disaggregate a GBD national number to districts.
8. **Value-of-information targeting** of where new surveys would most improve the
   evidence base for everyone.

We position this as triangulation and augmentation, not replacement. The
scientifically interesting frontier is where the two disagree (the vitamin A
decline, the direction of iron, diet-supply versus biomarker zinc).

## 8. Ideas for extending the analysis

Methodological (prototyped or planned; full plan in
`METHODS_TODO_IMPLEMENTATION_PLAN.md`, with honesty flags on novelty):
- **Decision-focused targeting**: tune the ensemble to correctly rank the
  highest-burden districts (a top-N / sensitivity-constrained objective in the
  ensembling step) rather than minimise average error.
- **National-anchor calibration** (prototyped): benchmark transported predictions
  to a known national prevalence to remove the level bias.
- **Proper binomial loss** (prototyped) for the bounded prevalence response.
- **Aggregate uncertainty and value-of-information** (prototyped): split-conformal
  and bootstrap intervals, partial-identification bands for transported aggregates,
  and a "where to sample next" map. The rigorous loss-based / transport-EIF
  inference is documented but flagged for statistician design.
- **Sequential / joint multi-outcome modelling** to borrow strength across
  co-occurring deficiencies (regressor chains; multivariate Fay-Herriot).

Data and scope extensions:
- **More countries and biomarker panels** (e.g. Kenya, Tanzania, Ethiopia national
  micronutrient surveys; DHS rounds with haemoglobin), to widen the LOCO panel
  beyond West Africa and address the "cross-sectional, four-country" limitation.
- **An anaemia validation track** on DHS haemoglobin (large samples, nearly every
  sub-Saharan country) to validate the methodology on a trusted, abundant outcome
  (noting anaemia is multifactorial, so it validates the pipeline, not
  micronutrient biology).
- **A formal GBD comparison** once GBD Results Tool data is sourced, focused on the
  vitamin A and iron divergences.
- **Finer resolution** where survey design allows (the Sierra Leone admin-3
  prototype) and a **public, policymaker-facing dashboard**.

## 9. Questions for external reviewers

**On the data and harmonisation:**
1. Is harmonising iron deficiency across BRINDA, Thurnham, and unadjusted
   ferritin defensible, or does residual non-comparability undermine cross-country
   pooling? Is there a better common definition?
2. How seriously should we take the disputed comparability of these biomarkers,
   and does it change the conclusions we can draw?
3. Are there proxy data domains we are missing that carry real micronutrient
   signal (beyond soil, climate, malaria, DHS)?

**On the methods:**
4. Given effective n equals the number of districts, is the area-level SAE the
   right primary analysis, and are our honest-evaluation choices (cluster/spatial
   block CV, in-fold preprocessing) adequate to avoid optimism?
5. Is the spatial-spline-plus-soil parsimonious model a credible transportable
   result, or are we over-reading a geographic gradient?
6. Is the reframe toward ranking / top-N targeting and value-of-information the
   right response to weak transport, and is the national-anchor calibration a
   sound way to handle the level bias?

**On the GBD relationship:**
7. Is the complementary-data-regime framing accurate, and which comparison
   (vitamin A decline, iron direction, biomarker vs diet-supply zinc) is most
   worth pursuing?
8. Could our covariates or biomarker estimates realistically feed into or
   constrain GBD's framework, and if so, how?

## 10. Pointers to deeper documentation

- Methods and results: `manuscript_mcn.qmd` and the supplementary appendix
  `manuscript_mcn_appendix.qmd`.
- Full data lineage (exact sources, fuzzy-match thresholds, per-domain
  missingness): `data_sources_appendix.qmd`.
- Exploratory findings: `exploratory_loco_findings.md`,
  `exploratory_fe_findings.md`, `exploratory_fe_sae_findings.md`.
- Methodological roadmap and implementation plan: `REVIEW_AND_ROADMAP_2026-06.md`,
  `METHODS_TODO_IMPLEMENTATION_PLAN.md`.
- A small shareable example dataset and analysis recipes: the `simplified subset/`
  folder.
