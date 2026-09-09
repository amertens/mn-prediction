# Brief for a 15-minute policy deck (written 2026-09-09 for a separate Claude Code session)

Build a concise slide deck (12 to 15 slides, 15 minutes) for a non-technical nutrition policy audience from the
results already on disk in this repository. Do not run the pipeline, `targets::tar_make()`, or any script under
`scripts/` that fits models; every number below is in a CSV under `results/tables/protocol_v2/` or
`results/tables/cluster_level/`, and you only read those. Make new figures with R (ggplot2; Rscript is at
`C:/Program Files/R/R-4.4.2/bin/Rscript.exe`) saved as PNG under `results/figures/policy_deck/`, and assemble the deck
as a Quarto file `docs/slides/MN-proxy-policy-deck-2026-09.qmd` rendered to PowerPoint (match the header of
`docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd`, which also shows how tables and the Ghana map are built from the
result tables). If Quarto is unavailable, build the .pptx with pptxgenjs from the same PNGs.

## Audience and style

- Nutrition policy makers and programme managers. No equations, no acronyms without a gloss, no variable codes on
  any slide (use the plain-language names below), no tables where a figure will do, at most one figure per slide,
  one message per slide in the title.
- Figures must read at slide size: base font 18 pt or larger, at most six categories per chart, direct labels rather
  than legends where possible, a grey band or line for "chance" wherever a correlation is shown, consistent colours
  across the deck (one colour for the proxy model, one for the survey baseline, one for the DHS-style geostatistical
  model, one for everything else).
- The score used throughout is the rank correlation between predicted and survey-measured district values
  (Spearman), where 1 is a perfect ranking and 0 is chance. Chance is not exactly zero here: the 95th percentile of a
  no-signal null is 0.08 for districts and 0.16 for regions (`transport_null_calibration.csv`). Say "ranking accuracy"
  on slides, not "Spearman".

## Very brief methods (two slides at most)

1. Four national micronutrient surveys (The Gambia 2018, Ghana 2017, Malawi 2015-16, Sierra Leone 2013;
   `metadata/survey_years.csv`) give district-level values for vitamin A, iron, folate and B12 in children and women
   (six outcomes; 206 surveyed districts: Gambia 30, Ghana 75, Malawi 87, Sierra Leone 14). Districts are below the
   surveys' design resolution, so the survey's own district values are noisy; the reliability ceiling per cell is in
   `variance_components_ceiling.csv`.
2. 454 public predictors in 24 domains (satellite imagery, climate, soil, crops, livestock, malaria, market prices,
   modelled child-nutrition surfaces, and aggregates of the DHS household survey) are attached to every district,
   including districts with no survey (`data/covariates/harmonized/predictors_admin2_shared_metadata.csv` for the
   list; `results/tables/protocol_v2/variable_sheet.csv` for plain descriptions).
3. The model of record is deliberately simple: each domain is summarised into a few components, each component is
   weighted by how strongly it tracks the outcome in the training districts, and the weighted sum ranks districts.
   There is nothing to tune, which is why it beats every more flexible method at 14 to 87 districts per country.
4. Three questions, each tested honestly by holding out what would be unknown: a district inside a surveyed country
   ("in-fill"), a whole region of a surveyed country ("region"), and a country with no survey at all ("new country",
   leave-one-country-out). The baselines are what a programme could do without the model: the survey's own regional
   average (computed without the district being scored) and a spatial smoother over neighbouring districts.

## Model comparison (numbers to plot; source files in brackets)

Mean ranking accuracy over the country-outcome cells (`benchmarks_v2_summary.csv`; level target = biomarker
concentration, prev = prevalence below the cut-off). In-country figures are over 18 scored cells (Sierra Leone's six
score NA with 14 districts); new-country figures over 22.

| Method | In-fill level / prev | Region level / prev | New country level / prev |
|---|---|---|---|
| Proxy index (model of record) | 0.40 / 0.29 | 0.38 / 0.27 | 0.28 / 0.21 |
| Spatial smoother, no covariates | 0.39 / 0.26 | 0.38 / 0.23 | cannot run |
| Smoother + covariates | 0.40 / 0.25 | 0.40 / 0.23 | cannot run |
| Penalised regression on all columns | 0.33 / 0.15 | 0.25 / -0.05 | not run |
| Penalised regression on domain components | 0.33 / 0.15 | 0.27 / 0.06 | 0.29 / 0.18 |
| Survey's own regional average (jackknifed) | 0.31 / 0.22 | cannot run | cannot run |
| Climate + soil only index (chosen after the fact) | | | 0.37 / 0.27 district, 0.45 regional (`nested_domain_selection.csv`, `climate_soil_admin1.csv`) |
| Twenty-predictor equal-weight composite | 0.36 / 0.26 | 0.31 / 0.23 | 0.30 / 0.20 (`weight_sources_summary.csv`, arm sparse20) |

DHS-style geostatistical model (INLA-SPDE spatial field + covariates, fitted at survey clusters, scored on the same
district folds; `results/tables/cluster_level/mbg_comparison_cells.csv`, 24 cells): in-fill level 0.27 against 0.38
for the index on the same cells (index better in 20 of 24), prevalence 0.21 against 0.22 (tie), region 0.28 against
0.36; the spatial field alone 0.20 / 0.10. It has lower absolute error on prevalence (9.2 against 10.7 percentage
points) because it shrinks toward the mean, and it cannot run for a new country at all.

Other methods tried and beaten by the index: SuperLearner ensembles under four losses (0.19 to 0.27 against 0.28 to
0.29 on the same cells, `sl_*_scores.csv`), random forest (0.25), highly adaptive ridge (0.21), ridge / lasso /
elastic-net weightings of the same components (`weight_sources_summary.csv` when the WS-01 rerun has landed; if the
file lacks arms other than sparse*, say "pending"), individual-level prediction (AUC 0.51 to 0.53,
`individual_level_models.csv`). Message: at these sample sizes, simpler wins.

## Variable importance (one or two slides)

- Existing figure: `results/figures/protocol_v2/index_importance_top10_level.png` (six panels, ten bars each, coloured
  by data source). Too dense for a policy slide; redraw with the top five per outcome, plain-language names, and
  three colour groups: remotely sensed environment (satellite, climate, soil, vegetation, water), household survey
  aggregates (DHS), modelled or administrative surfaces (IHME, livestock, malaria, prices). Data:
  `index_importance_top.csv` (columns: outcome, target, rank, column, beta_std, incountry_sign_agree, incountry_fits;
  use target == "level"; beta_std positive means more deficiency).
- Plain-language names for the codes that appear in the top lists: lcover_grass_frac_t0 = grassland cover;
  npp_gpp_t0 and npp_npp_t0 = vegetation productivity; ihme_stuntingprevalence / ihme_underweightprevalence =
  modelled child stunting / underweight; wpop_share_under5 = share of population under five; wpop_dependency_ratio =
  dependency ratio; dhs_w_primary_edu = women with primary schooling; dhs_w_no_education = women with no schooling;
  fprice_staple_rel = relative staple food price; glw_tlu_per_capita = livestock per person; glw_sheep_km2 = sheep
  density; glw_cattle_km2 = cattle density; glw_pigs_km2 = pig density; glw_ruminant_share = ruminant share of
  livestock; dhs_w_owns_house = women owning their home; dhs_w_working = women in paid work; dhs_hh_cows = household
  cattle ownership; spam_share_cereals = cereal share of cropland; map_blooddisorders...hbc... = haemoglobin C gene
  frequency; soil_phosphorus_stdev_0_20 = soil phosphorus variability; ihme_*anemia = modelled anaemia;
  dhs_FP_CUSA_W_MOD = modern contraceptive use; espen_sth_cov_mean = deworming coverage; dhs_w_health_insurance =
  health insurance; dhs_AN_NUTS_W_THN = thin women; dhs_CN_NUTS_C_HA2 / _WH2 = stunted / wasted children (survey);
  tclim_pdsi_t0 = drought index; wdist_perm_km_mean = distance to permanent water; wdist_coast_km_mean = distance to
  coast; lcover_crops_frac_t0 = cropland cover; wapor_sd_t0 = vegetation seasonality; grassland_frac = grassland
  share; elevation = elevation; map_sy_pf_* = malaria transmission; spam_share_oilcrops = oil-crop share.
- The pattern to state: one gradient recurs in every outcome. Districts that are drier and grassier, less productive,
  more pastoral, poorer (less home ownership, fewer women in work, less schooling) and more stunted rank worst on all
  six deficiencies. Livestock density ranking positive for iron and B12 is ecological, not dietary (pastoral zones are
  the poor, dry zones); say so, and do not present any weight as a cause.
- Two importances that disagree usefully (`index_importance_domains.csv`, scope == "pooled", target == "level", for
  shares; `domain_ablation_loco_summary.csv` for the cost of dropping a domain from new-country transport): satellite
  imagery, climate and soil together carry 40 to 50 percent of the model's variance and soil and climate are also
  the domains whose removal hurts a new country most (-0.024, -0.015); the household-survey domains hold several of
  the largest single weights yet removing all 130 DHS columns IMPROVES new-country accuracy by 0.056 (0.275 to
  0.331, `source_ablation_loco_summary.csv`). A scatter of "share of the model" (x) against "cost of dropping for a
  new country" (y), one dot per domain, with the four quadrants labelled, is the figure.

## Transportability (two or three slides)

- Rankings cross borders; absolute levels do not (raw ferritin differs six-fold between surveys). Say "the model
  tells a new country which districts are worst, not how bad".
- New-country accuracy: full index 0.28 (positive in 17 of 22 cells), climate + soil 0.37 (22 of 22) at district
  level and 0.45 at regional level, twenty-predictor composite 0.30 (20 of 22); chance line 0.08 (district) / 0.16
  (region). The climate + soil pair was chosen on these same four countries, so label it "to be confirmed on the next
  country" (it is pre-registered).
- Each added training country buys about 0.05 (`training_country_curve.csv`, mean spearman by n_train_countries for
  arm domain_index, level): a simple line chart.
- What the ranking is worth: directing effort to the worst fifth of districts by the model reaches 22 percent of the
  deficient people against 20 percent for the survey's regional averages and 14 percent for no information
  (`nce_targeting_summary.csv`, in-fill, mean_capture; oracle 47 percent). Risk bands: transported rankings put a
  district in the right WHO band 52 to 54 percent of the time and within one band 91 to 92 percent
  (`risk_category_accuracy_summary.csv`, scheme bands_5_20_40). A survey designed to sample regions and let the model
  rank districts gives a district error of 10.3 points against 12.3 for regional averages alone and 19.5 for a thin
  district survey of the same size (`anchor_and_rank_summary.csv`; see the sandbox log entry AR-01).
- Map slide: reuse the Ghana held-out map from the existing deck (survey rank vs predicted rank, side by side), it is
  the single most persuasive image.

## Figure ideas (simple, one per slide)

1. Lollipop chart: methods on the y axis, ranking accuracy on x, three small panels (in-fill, region, new country),
   grey band = chance, the proxy index highlighted. Replaces the model table.
2. "Simpler wins" strip: a single row of dots for index, ridge, SuperLearner, random forest, elastic net, lasso,
   individual-level, ordered by accuracy, with the tuned methods in grey.
3. Top-five-predictors chart per outcome, plain names, three colour groups, sign shown by bar direction.
4. Domain scatter: share of the model vs cost of dropping for a new country.
5. Learning curve: accuracy in a new country vs number of training countries.
6. Targeting bars: share of deficient people reached by the worst fifth of districts under model / survey average /
   no information / perfect knowledge.
7. Ghana map pair (survey vs predicted rank).
8. "Twenty public layers" slide: the child iron composite as an icon list grouped by source, with signs
   (`index_importance_top.csv`, outcome child_iron, target level, rank <= 20; all twenty replicate in every country).

## Caveats that must appear (one closing slide)

- Scores are rank correlations over few districts; differences under 0.03 between methods are ties.
- In-country, a spatial smoother with no covariates does almost as well: the model's value is in unsurveyed regions
  and countries, and in ranking, not in absolute prevalence.
- Prevalence numbers from the model are not calibrated; the geostatistical model is better for a number to publish.
- The climate + soil result is post hoc until the next country confirms it.
- Person-level prediction does not work with these data and is not proposed.

## Where the narrative lives

`docs/findings/SANDBOX_LOG_2026-09.md` (entries MB-01, WS-01, WS-02), `docs/manuscript_mcn_v2.qmd` (Results and
Supplement S4), `docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd` (the technical deck this one simplifies).
