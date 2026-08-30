# sandbox_parsimony — improving the area-level (GEE-only) and LOCO models

Exploratory sandbox. Nothing here is wired into `_targets.R`; it reads the
cached `_targets_full/objects/` store directly so experiments run in seconds
instead of hours.

Start with [`FINDINGS.md`](FINDINGS.md). `out/SUMMARY.txt` holds every headline
table; regenerate it with `Rscript sandbox_parsimony/R/99_summary.R`.

```r
# from the project root, in order
Rscript sandbox_parsimony/R/00_assemble.R          # pooled Admin-2 dataset (once, ~1 min)
Rscript sandbox_parsimony/R/01_noise_audit.R       # how much signal is even there
Rscript sandbox_parsimony/R/07_join_diagnostics.R  # rule out a broken join first
Rscript sandbox_parsimony/R/04_within_country.R    # within-country bake-off (~45 min)
Rscript sandbox_parsimony/R/13_combine_space.R     # + geography (produces the winner)
Rscript sandbox_parsimony/R/08_spatial_models.R    # binomial GAM + spatial field (slow)
Rscript sandbox_parsimony/R/05_loco.R              # LOCO round 1
Rscript sandbox_parsimony/R/09_loco_spatial.R      # LOCO round 2
Rscript sandbox_parsimony/R/11_loco_headline.R     # LOCO, comparable to benchmarks_all.csv
Rscript sandbox_parsimony/R/12_scale_sweep.R       # logit vs raw response scale
Rscript sandbox_parsimony/R/06_pool_composition.R  # cost of an incomplete GEE country
Rscript sandbox_parsimony/R/10_hygiene_effect.R    # GEE_COVARIATE_HYGIENE on vs off
Rscript sandbox_parsimony/R/14_decision_accuracy.R  # is it accurate enough to use?
Rscript sandbox_parsimony/R/15_level_and_target.R   # spatial level; biomarker vs indicator
Rscript sandbox_parsimony/R/16_var_importance.R     # what the 16 variables contribute
Rscript sandbox_parsimony/R/20_recommended_recipe.R # the two recipes, runnable
Rscript sandbox_parsimony/R/23_travel_time.R        # downloads 406 MB once
Rscript sandbox_parsimony/R/21_extend_covariates.R  # extended covariate table
Rscript sandbox_parsimony/R/17_admin1_layer.R       # Admin-1 + adaptive selection
Rscript sandbox_parsimony/R/18_mixed_level_check.R  # can levels be mixed?
Rscript sandbox_parsimony/R/19_biomarker_distribution.R
Rscript sandbox_parsimony/R/24_national_recovery.R  # national recovery + representativeness
Rscript sandbox_parsimony/R/30_vmnis_national.R     # VMNIS panel (reads a 522 MB .dta once)
Rscript sandbox_parsimony/R/31_national_model.R     # honest national model
Rscript sandbox_parsimony/R/99_summary.R            # collate everything
```

## Files

| Script | What it answers |
|---|---|
| `00_assemble.R` | Builds one row per country × Admin-2 per outcome: `svy_prev`, `svy_prev_se`, `n_svy` + every `gee_`/`ihme_`/`MAP_`/`fsec_` covariate + centroid lon/lat. Caches to `out/pooled_all.rds`. |
| `01_noise_audit.R` | Variance decomposition of the Admin-2 outcome into signal vs survey sampling noise. Produces `r_max`, the largest Pearson r any model could achieve. Also checks what the degenerate design SEs do to Fay-Herriot shrinkage. |
| `02_features.R` | Four parsimonious predictor sets: a priori curated (~16 constructs), correlation-cluster representatives, per-domain PCs, and the production top-K screen for comparison. |
| `03_core.R` | Shared CV harness, model fitters, and metrics (including precision-weighted r and `r_share` = r / r_max). |
| `04_within_country.R` | Within-country bake-off, repeated 5-fold CV, 15 model/feature combinations, including the null and spatial-only benchmarks the production pipeline lacks. |
| `05_loco.R` + `05a_loco_fns.R` | LOCO round 1: production recipe, the two centering variants, parsimonious sets; level and pattern reported separately. |
| `06_pool_composition.R` | What the pooled covariate intersection costs when a country with an incomplete GEE extraction joins. |
| `07_join_diagnostics.R` | Moran's I on covariates and outcome + Admin-2 key match rates, to rule out a broken join before blaming the model. |
| `08_spatial_models.R` | Binomial GAM on survey counts with a `s(lon, lat)` field on a parsimonious covariate set. **Negative result** — it loses to the simpler models. |
| `09_loco_spatial.R` + `09_loco_fns.R` | LOCO round 2: spatial + centering + within-country z-scored target. |
| `10_hygiene_effect.R` | Scores the existing `GEE_COVARIATE_HYGIENE` flag, which is written but off by default. |
| `11_loco_headline.R` | One LOCO table on exactly the cells `results/tables/benchmarks_all.csv` uses, so the numbers are directly comparable. |
| `12_scale_sweep.R` | `logit` vs raw prevalence as the response scale. |
| `13_combine_space.R` | Gives the winning covariate models geography, two ways (lon/lat as predictors; spatial smoothing of residuals). Produces the best within-country model. |
| `14_decision_accuracy.R` | Scores the winning model the way a programme consumes it: WHO severity band, worst-20% recall, and MAE, all against the national-mean map as the no-model baseline. |
| `15_level_and_target.R` | Reliability at cluster / Admin-2 / Admin-1, and a head-to-head between predicting the deficiency indicator and predicting the continuous biomarker. |
| `16_var_importance.R` | Permutation importance of the 16 a-priori constructs, plus the note on how they were chosen. |
| `20_recommended_recipe.R` | `area_model_v2()` and `area_transport_v2()` — the two recipes as functions, plus a demo run. |
| `17_admin1_layer.R` + `17_admin1_layer_fns.R` | Admin-1 layer: aggregates outcomes and covariates up a level, and runs the pre-specified vs data-adaptive predictor bake-off (selection strictly inside each training fold, plus a deliberately leaky variant to size the optimism). |
| `18_mixed_level_check.R` | Is it appropriate to mix Admin-1 and Admin-2 across countries? MAUP check, level/country confounding, and pooled LOCO under three unit schemes. |
| `19_biomarker_distribution.R` | Location-scale model of the district biomarker distribution vs the deficiency indicator. **Negative result** — modelling the scale does not help. |
| `21_extend_covariates.R` | Rebuilds the Admin-2 covariate table with the diet-quality pathways: travel time, MapSPAM crop mix, Köppen/AEZ agro-ecology, and the rich GEE families. Also drops GADM water bodies and collapses duplicate Admin-2 names. |
| `22_fetch_tanzania_access.R` | Attempt to fetch Tanzania accessibility from MAP. Kept for reference; the service returned an empty raster. |
| `23_travel_time.R` | Travel time to cities for all five countries from one consistent global product (406 MB download, cached). |
| `24_national_recovery.R` | How well the LOCO models recover national prevalence, and whether the surveyed districts are representative of the country. |
| `30_vmnis_national.R` | Builds a national panel from the full WHO VMNIS dump (187 countries) — the 8-country CSV in `results/external/` is a small subset. |
| `31_national_model.R` | Honest leave-one-COUNTRY-out national prevalence model on VMNIS + World Bank indicators. No benchmarking, no anchor. |
| `32_vmnis_heterogeneity.R` | Models VMNIS survey-method heterogeneity within country; recomputes the ceiling with method effects removed; permutation importance of the national predictors. |
| `33_vmnis_standardised.R` | Residualises the VMNIS outcome out-of-fold and predicts the method-standardised quantity. |
| `34_rescore_spatial_plus_soil.R` | Re-scores the leading benchmark method with its features selected inside the LOCO loop. |
| `35_hierarchical_biomarker.R` | Hierarchical individual-biomarker model with district random effects. **Negative result.** |
| `merge_tanzania_chunks.R` | Reassembles the chunked Earth Engine downloads and restores band names. |
| `99_summary.R` | Collates every headline table into `out/SUMMARY.txt`. |

`python/` holds the Earth Engine fetchers (`fetch_tanzania_rasters.py` is the
one on the working path; it needs the venv described in FINDINGS §21).

Outputs land in `sandbox_parsimony/out/`. The `.rds` caches and run logs are
gitignored; the result CSVs are not.

## Evaluation conventions used here

* **Repeated** 5-fold CV (fixed seeds), not a single split — several published
  differences are smaller than the fold-to-fold SD.
* Every table carries a **null** (country mean) and a **spatial-only** row.
* `r_max` and `r_share` put the correlations on a scale where "0.30" can be
  read as good or bad depending on how noisy the outcome is.
* For LOCO, level (bias, RMSE) and pattern (Spearman, level-removed RMSE) are
  never collapsed into one number.
