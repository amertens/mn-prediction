# Critical review: area-level (GEE-only) and LOCO models

Everything below is measured on the cached `_targets_full` store — same
outcomes, same covariates, same districts the production pipeline uses. Scripts
and output CSVs are in this folder; see `README.md` for how to re-run.

---

## Summary

The area-level and LOCO models are built sensibly, but three things are working
against them, in descending order of size:

1. **The outcome is noise-dominated and nothing reports that.** A median
   Admin-2 district contributes 6–36 biomarker measurements. Once sampling
   variance is subtracted from the between-district spread, several
   country × outcome cells have *no detectable signal at all* — the model is
   being scored against a target that is essentially white noise.
2. **Geography is left out of both the model and the benchmark table.** A
   covariate-free lon/lat smoother scores *higher* than the 237-predictor
   production recipe. The pipeline's own `benchmarks_all.csv` already shows the
   spatial methods leading the LOCO table; the production recipe is not one of
   them.
3. **The `logit` response scale costs accuracy.** 11–89% of districts are
   observed at exactly 0% and get clamped to 0.005, i.e. logit = −5.3. Fitting
   on the raw prevalence scale instead is worth up to +0.17 LOCO Spearman.

On the specific question of **parsimony**, the answer splits:

* **Within country: yes, decisively.** 16 a-priori mechanistic variables plus
  longitude and latitude, in a random forest, beat the 237–487-predictor
  production recipe in **14 of 16** country × outcome cells (paired mean
  ΔSpearman +0.063, t = 4.9). 18 predictors instead of 237.
* **Across countries (LOCO): no.** Every parsimonious set underperformed ridge
  on the full pool. What helps LOCO is not fewer predictors but a different
  *target* — see §6.

Two things I expected to help and that **did not**, reported because the
negative result is useful: a binomial likelihood on the survey counts, and a
thin-plate spline field instead of raw coordinates. Both are more principled
than what they replace; at 30–170 districts both lose to the simpler model
(§5b).

---

## 1. The Admin-2 outcome is mostly sampling noise

`R/01_noise_audit.R`, `out/noise_audit.csv`

Decompose the between-district variance of `svy_prev` into signal and sampling
noise: λ = (Var_obs − Var_sampling) / Var_obs, and r_max = √λ is the highest
Pearson r *any* model could reach, even one that predicted the truth exactly.

Because most districts contain a single PSU, `svy_prev_se` collapses to ~1e-17
(57–100% of districts) and cannot be used to estimate a design effect. λ is
therefore reported as a band over deff ∈ {1.0, 1.5, 2.0} with a district
bootstrap.

Reading the published within-country CV correlations against that ceiling:

| outcome | country | published CV r | r_max (deff 1.5) | r_max upper | % of ceiling |
|---|---|---|---|---|---|
| child_iron | ghana | 0.609 | 0.672 | 0.790 | **77–91%** |
| women_iron | gambia | 0.538 | 0.706 | 0.796 | **68–76%** |
| child_vitA | gambia | 0.426 | 0.671 | 0.823 | **52–63%** |
| child_iron | malawi | 0.038 | 0.199 | 0.540 | 7–19% |
| child_iron | gambia | 0.087 | 0.565 | 0.753 | 12–15% |
| women_vitA | ghana | −0.039 | 0.000 | 0.055 | — |
| women_vitA | malawi | 0.176 | 0.000 | 0.000 | — |

Two different conclusions live in that table, and the pipeline currently
presents them identically:

* Ghana child-iron and Gambia women-iron are **close to saturated**. More
  covariates or fancier learners will not move them much.
* Malawi has **genuine headroom** (7–19% of ceiling).
* Ghana and Malawi **women_vitA have nothing to predict** — even the optimistic
  bound on r_max is 0.06 and 0.00. Reporting r = −0.039 there reads as a model
  failure when it is a data limit.

**Recommendation.** Publish `r_max` beside every reported correlation, report
r / r_max, and suppress (or explicitly flag) cells whose ceiling is
indistinguishable from zero.

## 2. Fay–Herriot is shrinking almost nothing

`out/fh_shrinkage_check.csv`

`dashboard/data-raw/_build_fh_layer.R:84` takes the sampling variance from
`svy_prev_se^2`, treats only `D <= 0` as bad, then floors at `1e-5`. A
degenerate SE of 1e-17 is finite and positive, so it survives the guard and
becomes D = 1e-5 — an implied SE of 0.32 pp for a district measured on six
children.

The consequence is the EBLUP weight γ = A/(A+D):

| | as coded | with D = deff·p(1−p)/n |
|---|---|---|
| mean γ (weight on the raw direct estimate) | 0.71 | 0.26 |
| % districts with γ > 0.95 (i.e. no shrinkage at all) | 29% | 0% |

Fay–Herriot exists to shrink unreliable direct estimates toward a model. As
coded, the EBLUP keeps about 71% of the raw direct estimate where it should keep
about 26% — and for 29% of districts it applies effectively no shrinkage at all,
precisely the districts that most need it.

**Fix.** Compute D from `n_svy` unconditionally (design-effect inflated), or
treat any SE below, say, `0.5 * sqrt(p(1-p)/n)` as degenerate and replace it.
The same `svy_prev_se` is read by `_build_bym2_layer.R:161` and
`R/corrected/p2_p6_methods.R:92`.

## 3. Evaluation gaps

* **No null and no spatial-only benchmark.** Adding them changes the ranking
  (§5). Without them there is no way to tell whether the GEE block is doing
  anything.
* **Single-split CV.** Fold-to-fold SD of Pearson r across repeats has median
  0.071 and reaches 0.27. Several published between-method differences are
  inside that.
* **Degenerate predictions score as `NA`, not as zero skill.** `spatial_coords`
  returns a constant for child_iron/Malawi, so `pearson_r` is NA and
  `mean(..., na.rm = TRUE)` silently averages that method over 15 cells while
  others get 16. A method that occasionally collapses is rewarded.
* **`spatial_plus_soil` is scored on the metric its features were chosen by.**
  Its eight locked SoilGrids variables come from
  `archive/sandbox/13_univariate_loco_ranking.R`, ranked by LOCO r on the *same*
  four outcomes and four held-out countries it is then evaluated on. Its
  ρ = 0.305 — the current leader of `benchmarks_all.csv` — is optimistic. The
  honest spatial comparator is `spatial_coords` (ρ = 0.267; reproduced here at
  0.274).
* **Joins are fine.** Worth stating because it rules out a common explanation:
  Admin-2 keys match at 100% for every country × outcome, and Moran's I on the
  covariates is 0.38–0.73, as it must be for raster aggregates.

## 4. The covariate pool

For the five-country `child_vitA` intersection, the 237 "common predictors" are:

* ~99 SoilGrids columns — 10 nutrients × ~9 near-identical summaries of one
  static raster,
* ~108 IHME columns — mostly HIV and male-circumcision *counts*, which are
  population-size proxies,
* ~27 MAP columns.

Screening "top 30 by |r|" on that pool spends most of its budget on duplicates
of one or two underlying constructs.

**Tanzania shrinks the pool for everyone.** *(Since fixed — see §21, which
fetched the missing rasters through the Earth Engine API and recovered six of
these families for every country. The diagnosis below is what was true before
that.)* `data/Tanzania_GEE_rasters/` did not exist; Tanzania had 20 of the 35
GEE families the other four share. Because
the pooled model intersects on covariate *names*, adding it deletes — for every
country — `accessibility` (travel time to cities, the standard market-access
proxy), `lst`, `terraclimate` (precipitation, PET, PDSI, soil moisture, VPD),
`productivity` (GPP/NPP), `dailyevi`, `lai8days`, `landcovertype` and
`landcoverlayers` (crop/tree/urban cover fractions), `esa` WorldCereal, `wapor`,
`ghsbuilts`, `ghspop`, `fldas`, `aerosoloptical`, `atmosphere`. What survives is
soil chemistry plus a few static urbanicity layers.

Measured cost (`out/pool_composition.csv`, same four countries, same held-out
sets, only the covariate pool differs): child_vitA LOCO mean ρ 0.095 → 0.134.
For the production top-30 recipe specifically, ρ −0.106 → +0.089.

**Recommendation.** Complete the Tanzania extraction, or keep Tanzania out of
the pooled/LOCO covariate intersection until it is complete. Also worth logging
the intersection size per run — this failure is silent.

`R/gee_band_semantics.R` already diagnoses the redundancy correctly and its
pruning is written; it is gated off by `GEE_COVARIATE_HYGIENE`. §8 measures it.

## 5. Within-country bake-off

`R/04_within_country.R`, `out/within_country_bakeoff.csv` — repeated 5-fold CV,
5 repeats, fixed seeds, 16 country × outcome cells.

Restricted to the 8 cells with a detectable ceiling (r_max ≥ 0.35):

| model | predictors used | ρ | r | r/r_max |
|---|---|---|---|---|
| `curated16_rf_xy` | **16 a-priori variables + lon/lat**, random forest | **0.374** | 0.376 | 0.58 |
| `decorr20_rf_xy` | 20 cluster representatives + lon/lat, RF | 0.362 | 0.394 | 0.60 |
| `decorr20_rf` | 20 cluster representatives, RF | 0.356 | 0.385 | 0.59 |
| `spatial_only` | **none** — lon/lat spline | 0.347 | 0.341 | 0.53 |
| `curated16_ridge` | 16 a-priori variables, ridge | 0.318 | 0.323 | 0.50 |
| `domainPC2_ridge` | 2 PCs per conceptual domain (~22 features) | 0.315 | 0.316 | 0.49 |
| `PROD enet_screen30` | top-30 screen, elastic net (production) | 0.310 | 0.320 | 0.49 |
| `binom_screen30` | same, binomial likelihood on counts | 0.296 | 0.310 | 0.48 |
| `PROD ridge_all` | all 237–487, ridge | 0.230 | 0.209 | 0.32 |

`curated16_rf_xy` versus the production recipe, paired cell by cell:
mean ΔSpearman **+0.063**, SE 0.016, **t = 4.0**, better in **8 of 8**
signal-bearing cells (and 14 of 16 across all cells, t = 4.9). This is the one
comparison in this review that is unambiguous.

The 16 variables are NDVI, elevation, population density, nightlights, global
human modification, settlement degree, built surface, precipitation, four
SoilGrids topsoil nutrients (carbon, zinc, iron, nitrogen), *Pf* parasite rate,
ITN use, and HbS and G6PD allele frequency — chosen from the causal story
before looking at any performance number.

Also worth taking from this table:

* **A covariate-free spatial smoother matches the production recipe** (0.347 vs
  0.310). That does not make the covariates useless — it means geography must be
  *in* the model rather than treated as an alternative to it. Adding lon/lat
  lifts the forest from 0.356 to 0.362–0.374.
* **Ridge on the whole pool is the worst covariate model tested** (0.230). The
  redundancy in the pool is actively harmful when nothing selects against it.
* Models capture roughly **half to three-fifths** of the achievable signal
  (r/r_max ≈ 0.5–0.6).

### 5b. Two principled changes that did not pay off

* **Binomial likelihood on the counts.** `binom_screen30` (0.296) vs its
  Gaussian-on-logit twin `PROD enet_screen30` (0.310); `curated16_binom` (0.305)
  vs `curated16_ridge` (0.318). Statistically it is the right model for a
  proportion built from n Bernoulli draws and it needs no clamping — but at
  30–170 districts it does not win. (It is still the right fix for the *LOCO*
  scale problem in §7, where the clamp does real damage.)
* **Binomial GAM with an `s(lon, lat)` field** (`R/08_spatial_models.R`). On the
  10 cells it was run on, the best GAM variant scored ρ = 0.298 against
  `spatial_only` 0.354 and `decorr20_rf` 0.350 on the same cells, with worse
  RMSE (17.7 vs 15.0–15.2 pp). Two structural problems: the spatial basis is
  silently dropped when a training fold has fewer than 25 districts (every fold
  in The Gambia), and deflating the trials by the design effect shrinks an
  already-small effective sample further.

The general lesson is that at this sample size, model *machinery* is not the
binding constraint — predictor selection and getting geography into the model
are.

## 6. LOCO

`R/11_loco_headline.R`, `out/loco_headline.csv` — 4 outcomes × 4 held-out
countries, the same cells as `benchmarks_all.csv`. Tanzania excluded (§4, and it
is 4,000 km away, so a spline can only extrapolate into it).

| variant | ρ | ρ SD | RMSE pp | level-removed RMSE pp | \|level bias\| pp |
|---|---|---|---|---|---|
| `zscore ridge_all` | **0.261** | 0.24 | 17.8 | 13.0 | 8.7 |
| `zscore screen30` | 0.190 | 0.29 | 18.6 | 14.0 | 8.7 |
| `centered_own ridge_all` | 0.187 | 0.26 | 19.7 | 11.8 | 13.1 |
| `PROD enet30` (production) | 0.154 | 0.30 | 14.5 | 10.0 | 8.9 |
| `PROD ridge_all` | 0.140 | 0.26 | 16.1 | 10.8 | 9.2 |
| `spatial_tps` (logit scale) | 0.106 | 0.33 | 16.5 | 10.5 | 11.5 |
| `spatial_tps` (raw scale) | 0.274 | — | 16.5 | — | 11.8 |

`zscore` = train on the **within-country z-score of the transformed
prevalence**, so between-country level and spread are structurally unlearnable
and the model can only carry spatial *pattern*; the level and spread come from
an anchor.

Paired cell-by-cell against the production recipe (`zscore ridge_all` −
`PROD enet30`, 16 cells): mean ΔSpearman **+0.107**, SD 0.256, SE 0.064,
t = 1.7, better in **11 of 16** cells. Suggestive and consistently positive, not
statistically decisive at this sample size — 16 country × outcome cells is
simply not much to compare methods on, which is itself worth saying in any
write-up of the LOCO results.

It is, though, the only variant that is positive for **all four** outcomes:

| variant | child_iron | child_vitA | women_iron | women_vitA |
|---|---|---|---|---|
| `zscore ridge_all` | 0.293 | **0.378** | 0.284 | **0.088** |
| `centered_own ridge_all` | 0.298 | 0.164 | 0.238 | −0.002 |
| `PROD enet30` | 0.246 | 0.089 | 0.224 | 0.055 |
| `PROD ridge_all` | 0.224 | 0.206 | 0.136 | −0.056 |
| `spatial_tps` (logit) | **0.363** | 0.243 | −0.014 | −0.169 |

The spatial spline is the strongest single method on child iron and vitamin A
and the *worst* on the women's outcomes — a reminder that "extrapolate a smooth
surface from West Africa into Malawi" is a fragile mechanism, not a transported
relationship.

**Level and pattern must be reported separately.** Supplying the held-out
country's *national* prevalence as the level — realistic, since a national
figure is nearly always published even when district estimates are not —
removes most of the catastrophic LOCO errors:

| outcome | held out | RMSE with train-mean level | with national anchor | saved |
|---|---|---|---|---|
| women_folate | sierraleone | 49.8 | 9.8 | **40.0 pp** |
| child_iron | gambia | 36.6 | 12.7 | **23.9 pp** |
| women_iron | gambia | 18.0 | 11.4 | 6.6 pp |
| everything else | | | | ≈ 0 |

That is the same conclusion as the existing note that LOCO fails on cross-survey
biomarker *level* offsets — but it makes it actionable: transport the map, take
the level from a national anchor, and report the two error components apart.

**The centering guard is over-conservative but nearly harmless.**
`R/transportability_area.R:~215` centers the held-out country's covariates on
the *training* means, with a comment that using its own means "would leak
held-out information". Only the held-out *covariates* are involved, and those
are observable everywhere without a survey — that is the premise of the whole
method, so it is not leakage. In practice, for a linear model, the two differ by
a constant, so correlations are identical and only the level/RMSE moves. Worth
correcting for clarity; not a source of lost accuracy.

## 7. The response scale is costing accuracy

`R/12_scale_sweep.R`, `out/scale_sweep.csv`

`AREA_TRANSPORT_RECIPE$scale = "logit"` with a clamp at 0.005/0.995. But
11–89% of districts are observed at exactly 0% — 84% for women_vitA in Ghana,
89% in Malawi — and every one of them becomes an extreme point at logit = −5.3
that the fit then chases.

Same models, same cells, only the response scale changes (LOCO ρ):

| model | raw prevalence | logit | gain |
|---|---|---|---|
| `spatial_tps` | 0.274 | 0.106 | **+0.168** |
| `centered_own ridge_all` | 0.248 | 0.187 | +0.061 |
| `PROD enet30` | 0.202 | 0.154 | +0.048 |
| `ridge_all` | 0.191 | 0.140 | +0.051 |
| `curated16_ridge` | −0.020 | −0.089 | +0.069 |
| `centered_own decorr20` | 0.150 | 0.163 | −0.013 |

**Recommendation.** Drop the logit default. Better still, fit the binomial
likelihood on the counts (`cbind(deficient, n - deficient) ~ X`, `glmnet`
`family = "binomial"` or `mgcv::gam`): it needs no clamp, handles the
zero-count districts natively, and weights each district by its own information
instead of treating logit(p) as homoskedastic. `fit_predict_quasibinomial`
already exists in `R/benchmark_models.R:943` — it just is not the default.

## 8. Turning on `GEE_COVARIATE_HYGIENE`

`R/10_hygiene_effect.R`, `out/hygiene_effect.csv`

The flag drops 52 of 237 (child_vitA) and 98 of 366 (child_iron) columns: the
cross-band `_annual_*` summaries over non-commensurable bands that
`R/gee_band_semantics.R` documents (TerraClimate averaging `aet`, `pdsi` and
surface pressure into one number; SoilGrids summaries that are bit-identical
copies of a real band). Paired on the same cells:

| eval | model | ρ off | ρ on | RMSE off | RMSE on |
|---|---|---|---|---|---|
| LOCO | `enet_screen30` | 0.212 | **0.228** | 30.9 | **26.3** |
| LOCO | `ridge_all` | 0.222 | **0.238** | 29.2 | 29.7 |
| LOCO | `decorr20_enet` | 0.133 | 0.052 | 24.5 | 25.0 |
| within | `enet_screen30` | 0.264 | **0.276** | 16.6 | 16.5 |
| within | `ridge_all` | 0.162 | **0.177** | 16.5 | 16.4 |
| within | `decorr20_enet` | 0.160 | 0.155 | 16.8 | 16.6 |

Small and positive for the screened and ridge models, negative for
`decorr20_enet` — which already de-duplicates by construction, so pruning first
just narrows what its clustering has to work with. The largest single effect is
LOCO RMSE for the production recipe: 30.9 → 26.3 pp.

The accuracy case is modest; the better argument for flipping the default is
that the pruned columns are provably meaningless (an average of surface
pressure and a drought index) and that removing exact duplicates stops
variable importance from splitting across nine copies of the same raster,
which is what currently destabilises the per-fold selected-variable lists.

## 9. The two recipes, written out

`R/20_recommended_recipe.R` implements both as functions that could be lifted
into `R/` with minimal editing.

**`area_model_v2()`** — within-country map. Random forest on the 16 a-priori
constructs plus longitude and latitude, on the logit scale with survey-n case
weights, falling back to correlation-cluster representatives for any construct a
country is missing.

**`area_transport_v2()`** — LOCO. Trains on the within-country z-score of
prevalence so level and spread are structurally unlearnable, then places the
held-out country using a national anchor.

Demo output on child iron (`Rscript sandbox_parsimony/R/20_recommended_recipe.R`):

```
--- area_model_v2, within Ghana, 5-fold CV ---
  r = 0.494, rho = 0.543, RMSE = 17.0 pp | ceiling r_max = 0.67 (74% of it)

--- area_transport_v2, LOCO, with and without a national anchor ---
  gambia       rho = +0.532 | RMSE anchored  13.2 pp vs unanchored  46.2 pp
  ghana        rho = +0.463 | RMSE anchored  16.7 pp vs unanchored  17.8 pp
  sierraleone  rho = +0.223 | RMSE anchored  12.1 pp vs unanchored  26.3 pp
  malawi       rho = +0.030 | RMSE anchored  18.9 pp vs unanchored  22.3 pp
```

Compare the production LOCO recipe on child iron: mean ρ 0.246, and
`transportability_area_loco_metrics.csv` reports Gambia at RMSE 45.8 pp with a
−43.9 pp national bias. Anchoring turns that into 13.2 pp.

Caveats worth carrying forward:

* the anchor in the demo is the held-out country's own weighted mean, i.e. the
  best case. In production it would be a published national figure, which
  carries its own error — small relative to a 44 pp offset, but not zero;
* **lon/lat flatter the model when whole regions are unsurveyed.** With random
  CV folds a district's neighbours are in the training set, so coordinates let
  the forest interpolate. That matches the real use case (an unsurveyed district
  surrounded by surveyed ones) but not the case where an entire Admin-1 region
  is unsurveyed. Stress-test with leave-one-Admin1-out CV before promising
  accuracy there;
* the `se` returned by `area_model_v2()` is forest uncertainty only. It does not
  include district sampling error and must not be presented as a credible
  interval for the underlying prevalence.

## 10. Is the map actually accurate enough to use?

`R/14_decision_accuracy.R`, `out/decision_accuracy.csv`

"Beats the other models" and "good enough to target with" are different claims.
Scored the way a programme would consume it — WHO severity band, and whether it
finds the worst-off districts — against the comparator a programme already has
without any model at all: **paint the national prevalence on every district**.

Averaged over the 11 cells with detectable signal:

| | model | national-mean map |
|---|---|---|
| WHO severity band correct | 60.2% | 53.6% |
| MAE | 13.9 pp | 15.4 pp |
| Worst-20% districts recovered | 27.1% | 20% (random) |
| Spearman | 0.287 | — |

**In absolute terms this is a weak map.** It beats the no-model baseline, but by
1.5 pp of MAE and 6.6 points of WHO accuracy, and it finds 27% of the worst
districts where random targeting finds 20% — a lift of **1.36×**. Some cells are
genuinely useful (Ghana child iron: WHO 40% vs 21%, lift 2.3×, MAE −2.6 pp,
ρ = 0.55; Gambia women's iron: WHO 80% vs 67%, MAE −4.6 pp, ρ = 0.67; Ghana
women's B12: WHO 54% vs 17%, lift 2.1×). Several are *worse* than the national
mean (Tanzania and Ghana child vitamin A).

Caveat on the WHO metric: where a whole country sits in one band it is trivially
high for everything — Gambia child iron scores 93% for both model and null
because every district is "Severe". The targeting lift is the more honest
number, and it is 1.36×.

## 11. The 16 variables, and which of them matter

`R/16_var_importance.R`, `out/var_importance.csv`

**How they were chosen.** Written down from the causal story before any
performance number was looked at — one variable per distinct mechanism, no two
measuring the same thing:

| pathway | variables |
|---|---|
| diet quality / agro-ecology | NDVI, precipitation (TRMM), soil organic carbon, soil zinc, soil iron, soil nitrogen |
| market access / poverty | nightlights, population density, built surface, settlement degree, global human modification, elevation |
| inflammation & blood loss | *Pf* parasite rate, ITN use |
| inherited anaemia | HbS allele frequency, G6PD allele frequency |

The binding constraint was availability in all five countries, and that is **not
a neutral constraint**: the five-country intersection has already lost travel
time to cities, land surface temperature, TerraClimate, cropland fraction and
vegetation dynamics (§4) — all with a better claim on the diet-quality pathway
than soil chemistry has. This is the best available set, not the best
conceivable one.

Permutation importance over the 11 signal-bearing cells (relative to the
strongest variable in each cell; `share_top5` = fraction of cells where it lands
in that cell's top five):

| construct | mean rel. importance | share top-5 |
|---|---|---|
| soil zinc | 0.64 | 0.73 |
| settlement degree | 0.53 | 0.64 |
| longitude | 0.48 | 0.45 |
| latitude | 0.43 | 0.36 |
| human modification | 0.37 | 0.36 |
| built surface | 0.35 | 0.27 |
| HbS allele frequency | 0.34 | 0.36 |
| nightlights | 0.32 | 0.27 |
| soil iron | 0.32 | 0.27 |
| precipitation | 0.30 | 0.18 |
| soil organic carbon | 0.30 | 0.18 |
| population density | 0.26 | 0.09 |
| *Pf* parasite rate | 0.25 | 0.18 |
| ITN use | 0.22 | 0.27 |
| G6PD allele frequency | 0.21 | 0.27 |
| elevation | 0.20 | 0.09 |
| soil nitrogen | 0.16 | 0.00 |
| **NDVI** | **0.08** | **0.00** |

Two things to be honest about:

* **NDVI is dead last and never makes a top five.** The single most-used
  remote-sensing covariate in this literature is contributing nothing here.
* **Soil zinc ranks first for outcomes it has no mechanism for** — vitamin A,
  B12, folate. Together with HbS allele frequency (a smooth continental
  surface) and longitude/latitude in the top four, the honest reading is that
  much of the "covariate" signal is **geography wearing a covariate label**.
  That is consistent with §5, where a covariate-free spatial smoother nearly
  matches the full model. Do not present soil zinc as a mechanism.

The mechanistically clean hits do show up where they should: *Pf* parasite rate
is the top variable for Malawi child zinc and Malawi women's folate, and HbS for
Malawi women's zinc.

## 12. Would a different spatial level help? Admin-1 yes, cluster no

`R/15_level_and_target.R`, `out/level_reliability.csv`

The number that limits a proxy-only model is not the number of individuals — it
is the number of **design points at which the covariates take distinct values**.
Checked directly:

* every proxy covariate in `merged_*`, including the ones named `_10km`, is
  **constant within Admin-2** (Ghana: 75 distinct values across 2,450 rows,
  between-district variance share 1.00);
* the survey coordinates are effectively Admin-2 too (Ghana: 90 distinct
  lon/lat for 75 districts, median 1 per district);
* so fitting at individual or cluster level cannot add information. The
  effective n is the district count either way. This is the mechanism behind
  the existing note that the effective n is the number of areas.

Reliability by level (mean over cells, deff 1.5):

| level | design points (median) | n per unit (median) | ceiling r_max, prevalence | ceiling r_max, biomarker mean |
|---|---|---|---|---|
| Admin-1 | 16 | 37 | **0.59** | **0.76** |
| Admin-2 | 75 | 13 | 0.31 | 0.63 |
| cluster | 90 | 10 | 0.22 | 0.57 |

**Going up to Admin-1 nearly doubles the achievable correlation; going down to
cluster level makes it worse.** The pipeline is currently mapping at a finer
resolution than these surveys can support. Going finer than Admin-2 would only
pay off after re-extracting covariates at true cluster coordinates — and even
then the surveys have just 60–105 GPS clusters each, so the ceiling would still
fall, not rise.

The tension is real, not a slam dunk: Admin-1 gives a far more reliable estimate
of a far less useful unit, and only ~16 units to fit on. The defensible answer
is to **publish both**, with the Admin-2 map carrying its reliability explicitly.

## 13. Predict the biomarker or the deficiency indicator?

`R/15_level_and_target.R`, `out/continuous_vs_binary.csv`,
`out/head_to_head_target.csv`

Thresholding a biomarker at a WHO cutoff throws away how far from the cutoff
each person sits, and it costs a lot: at Admin-2 the district mean log-biomarker
has a ceiling of **r_max = 0.63** against **0.31** for the district prevalence.
Twice the achievable signal. The district mean is also tightly coupled to the
prevalence it implies (|r| = 0.43–0.90 across cells, median ≈ 0.7).

Head-to-head — same learner, same folds, both scored against the observed
Admin-2 prevalence, 18 cells:

| route | mean ρ gain vs indicator | cells better | MAE change |
|---|---|---|---|
| biomarker → normal conversion | **+0.043** (SE 0.029, t = 1.5) | 10/18 | −2.0 pp (worse) |
| … + per-fold intercept shift | +0.011 | 9/18 | −0.1 pp (parity) |
| … + per-fold linear calibration | +0.004 | 7/18 | −0.1 pp (parity) |

So it is **promising but not yet a win**. The raw conversion improves the
ranking and damages the level; both recalibrations repair the level and give
most of the ranking gain back. (I expected the intercept-only shift to leave
Spearman untouched — it does not, because the shift is re-estimated per fold, so
the assembled out-of-fold vector is only piecewise monotone.)

Where it helps is interpretable: the gain is largest where the indicator route
has almost no signal — Ghana women's vitamin A −0.08 → +0.24, Malawi child iron
+0.03 → +0.27, Malawi women's vitamin A +0.00 → +0.20, Ghana child vitamin A
+0.18 → +0.35 — and it is slightly negative where the indicator already works
(ρ ≥ 0.4: mean −0.013). Correlation between baseline ρ and the gain is −0.49.
Part of that pattern is regression to the mean and should not be over-read.

The conclusion is not "patch the conversion harder". It is that the two-step
mean → normal-CDF → prevalence route is the wrong shape: it assumes one
within-district SD and a normal distribution. The right version models the
**district biomarker distribution** (location *and* scale, lognormal or
skew-normal), or fits one hierarchical model to individual biomarkers with
district random effects and derives prevalence from the posterior. That is a
real piece of work, and §12 says it should be attempted at Admin-1 as well as
Admin-2.

---

# Follow-up round: acting on the recommendations

## 14. Closing the diet-quality covariate gap

`R/21_extend_covariates.R`, `R/23_travel_time.R`

Section 4 named the pathways missing from the five-country intersection.
Most of the gap was closeable without leaving the machine:

| added | source | coverage |
|---|---|---|
| **travel time to cities** (mean, median, log-mean, p90) | Nelson et al. 2019 via `geodata::travel_time()`, 406 MB global 30-arcsec GeoTIFF, CC BY 4.0 | all 5 |
| **crop mix** — production and area share of cereals, roots, pulses, oilcrops, vegetables, plus a yield proxy | MapSPAM, already in `data/MapSPAM/` and previously unused by the area model | all 5 |
| **agro-ecology** — Köppen-Geiger and AEZ class shares + a Shannon diversity index | global rasters already in `data/` | all 5 |
| accessibility, LST, TerraClimate, productivity, EVI, LAI, land-cover fractions, ESA WorldCereal | `area_covariates_*` on disk | 4 (not Tanzania) |

Five-country common covariates: **243 → 263**. For the four countries with a
complete GEE extraction the gain is much larger, because the rich families come
back in full.

Two notes on how this was done:

* **One product, not four.** Gambia, Ghana, Sierra Leone and Malawi already have
  `Accessibility_<Country>_2019.tif` (Weiss et al. 2018) on disk; Tanzania has
  none. Rather than mix a MAP-derived column for four countries with something
  else for the fifth — the exact failure mode §4 documents — travel time was
  taken from a single global product for all five. Its correlation with the
  existing per-country MAP clips is only 0.41–0.73, which confirms they are not
  interchangeable.
* Attempting to fetch Tanzania's accessibility from MAP's own service
  (`malariaAtlas::getRaster`) returned an empty raster; the code is kept in
  `R/22_fetch_tanzania_access.R` but is not on the working path.

### 14b. Two data-integrity bugs found while doing this

**GADM water bodies are being used as prediction units, and two are in the
survey data.** Malawi has 4 inland-water Admin-2 polygons and Tanzania 5.
Covariates are extracted over open water — NDVI, soil chemistry and LST of a
lake. Worse:

* **`Lake Malawi`** appears in the Malawi survey tables with `n_svy` 21–33;
* **`Lake Victoria`** appears in the Tanzania survey tables with `n_svy` 60–110,
  making it one of Tanzania's larger "districts".

These rows are in both the fitting data and the map.

**CORRECTED in §32.** Dropping them was the wrong response: checking the
coordinates shows all four cluster locations sit within the DHS/MICS GPS
displacement bounds of the nearest land district, so they are real lakeside
communities displaced into the water, not errors. They are now snapped back to
land rather than discarded. Water polygons are still dropped from the covariate
and prediction side.

**Duplicate Admin-2 names cause silent cartesian joins.** Beyond the lakes,
Malawi has genuinely repeated Traditional Authority names across regions
(`TA Lundu`, `TA Malemia`, `TA Ngabu`, `TA Pemba`). Every join in this project
is keyed on the Admin-2 *name*. Unguarded, that inflated the Malawi covariate
table from 256 rows to **29,183** in my first run. In production the same
mechanism turns Tanzania's `inner_join(svy_admin2, gee_admin2)` from 163 rows
into 167. Small there, but it is silent, and it is the same defect.

## 15. Is it appropriate to mix Admin-1 and Admin-2 across countries?

`R/18_mixed_level_check.R`, `out/maup_check.csv`, `out/mixed_level_loco.csv`

Short answer: **no, not to equalise unit counts** — though the reasoning matters
more than the verdict.

**(a) The covariate–outcome relationship is mostly but not entirely
level-invariant.** Over 320 variable × country × outcome comparisons, the
correlation between r(Admin-1) and r(Admin-2) is **0.87** — reassuring. But
**14% flip sign** between levels, and the standardised slope shrinks about 20%
at Admin-1 (1.42 vs 1.74 pp per SD). Stability varies a lot by variable: *Pf*
parasite rate and NDVI are the least stable, elevation and the blood-disorder
frequencies the most. And by outcome: child zinc is unstable (r = 0.50 between
levels), women's folate flips sign in 25% of comparisons. A model fitted on
mixed levels is therefore fitting a blend of two related-but-different
relationships. This is the modifiable areal unit problem, and it is the
substantive objection.

**(b) It confounds level with country and with reliability.** Under the obvious
rule (Admin-1 where a country has ≥15 regions), Ghana, Malawi and Tanzania go
Admin-1 while Gambia and Sierra Leone stay Admin-2. Level, country and noise
then vary together, so a LOCO model that "fails to transport" to Sierra Leone
could just be meeting a different measurement-error regime.

**(c) It does not equalise the thing that actually matters.** Unit *counts* do
get more even (14–31 instead of 14–162). But people per unit gets *less* even:

| country | Admin-1 units / n per unit | Admin-2 units / n per unit |
|---|---|---|
| Gambia | 6 / 185 | 30 / 24 |
| Ghana | 16 / 45 | 75 / 12 |
| Malawi | 27 / 27 | 86 / 8 |
| Sierra Leone | — (only 4 regions) | 14 / 47 |
| Tanzania | 31 / 248 | 162 / 36 |

All-Admin-2 spans 8–47 people per unit; the mixed scheme spans 24–248. Since
reliability is driven by people per unit, mixing makes the heteroskedasticity
worse, not better.

**(d) Empirically it buys nothing.** Pooled LOCO under the three schemes:
all-Admin-2 ρ = −0.024, mixed ρ = −0.020, all-Admin-1 ρ = 0.050 (4 countries
only). No scheme transports; mixing is indistinguishable from Admin-2.

**What to do instead.** If the goal is that each country contributes comparable
information, address that directly rather than by changing the geography:

* fit at the finest **common** level and carry an explicit measurement-error
  term, so each unit is automatically weighted by its own reliability — a
  Fay–Herriot / precision-weighted model does this by construction, and it is
  already in the codebase;
* or weight rows by country so no country dominates;
* and **decouple the modelling level from the reporting level**. Fitting at
  Admin-2 everywhere and *reporting* Admin-1 where the survey cannot support
  Admin-2 is coherent; fitting on a mixture is not.

The one place mixing is defensible is a deliberately per-country product, where
each country's map is fitted and reported at whatever level its own survey
supports and no pooled model is fitted across them.

## 15b. Do the new covariates actually carry signal?

Univariate |r| with Admin-2 prevalence, averaged over the 11 signal-bearing
cells, with the share of the achievable ceiling each variable reaches on its own:

| new variable | mean \|r\| | max \|r\| | share of ceiling |
|---|---|---|---|
| travel time to cities (log-mean) | 0.228 | 0.627 | 0.36 |
| travel time to cities (median) | 0.223 | **0.696** | 0.34 |
| MapSPAM pulses share | 0.220 | 0.498 | 0.36 |
| MapSPAM cereals share | 0.209 | 0.543 | 0.32 |
| MapSPAM oilcrops share | 0.205 | 0.469 | 0.32 |
| MapSPAM roots share | 0.182 | 0.451 | 0.29 |
| MapSPAM vegetables share | 0.177 | 0.446 | 0.29 |
| Köppen diversity | 0.143 | 0.363 | 0.22 |
| AEZ diversity | 0.121 | 0.226 | 0.19 |
| MapSPAM yield proxy | 0.107 | 0.270 | 0.18 |

Yes. Travel time alone reaches a third of the achievable correlation on average
and 0.70 in its best cell — comparable to the strongest variables in the
original curated set (§11), and mechanistically far easier to defend than soil
zinc. Crop-mix shares sit just behind it. These belong in the model.

## 16. Modelling the biomarker distribution (next step 2)

`R/19_biomarker_distribution.R`, `out/biomarker_distribution.csv`

§13 concluded the two-step mean → normal-CDF → prevalence route was the wrong
shape because it assumes a single pooled within-district SD. So I modelled the
scale too: predict μ and σ from the proxies, then p = Φ((cut − μ)/σ). Also a
variant replacing the normal with the pooled empirical distribution of
within-district standardised biomarker values.

18 country × outcome cells, same folds and learner, all scored on the observed
Admin-2 prevalence:

| route | ρ | MAE pp | mean abs bias pp |
|---|---|---|---|
| mean + **pooled** SD | **0.276** | 12.13 | 5.47 |
| mean + **modelled** SD (location–scale) | 0.252 | 12.09 | 5.23 |
| location–scale, empirical shape | 0.252 | 11.80 | 5.43 |
| deficiency indicator (production) | 0.228 | **10.39** | **3.53** |

Paired against the indicator route: pooled-SD **+0.048 ρ** (t = 1.94, 11/18
cells); location–scale **+0.023** (t = 0.83).

**My proposed fix did not work.** Modelling σ makes things worse, not better —
in hindsight for a clear reason: the district SD is itself estimated from ~13
people, so it is mostly noise, and the proxies have little purchase on it.
Predicting a noisy quantity and dividing by it adds variance.

So the honest state of §13 stands: the biomarker route buys a modest,
borderline-significant ranking gain and costs ~1.7 pp of MAE. The remaining
idea — a hierarchical model on individual biomarkers with district random
effects, deriving prevalence from the posterior rather than from a plug-in
conversion — is untested here.

**Where it is worth using is now clear**, and it is a usable rule. Splitting the
18 cells by whether the indicator route has any signal at all:

| | cells | indicator ρ | biomarker ρ | gain |
|---|---|---|---|---|
| indicator route has no signal (ρ < 0.2) | 9 | 0.040 | **0.129** | **+0.090** |
| indicator route works (ρ ≥ 0.2) | 9 | 0.417 | 0.422 | +0.006 |

So: **fit both, and use the biomarker route only for the outcomes where the
prevalence model has failed.** That is most of the low-prevalence outcomes —
vitamin A in women, child iron in Malawi — where thresholding leaves too few
deficient people per district for the proportion to carry information.

## 16b. Production fixes applied

Four of the quick fixes are now in the repo, not just recommended:

| change | file |
|---|---|
| Fay–Herriot no longer accepts degenerate design SEs: any SE below 0.25·p(1−p)/n is replaced, not just SE ≤ 0 | `dashboard/data-raw/_build_fh_layer.R` |
| same guard for the BYM2 precision multiplier | `dashboard/data-raw/_build_bym2_layer.R` |
| new `admin2_reliability()` / `add_reliability_columns()`; every area-level comparison table now carries `r_max`, `r_share` and a `signal` flag | `R/area_reliability.R` |
| the area-level comparison gains a **National mean (null)** row and a leave-one-out **Spatial only (no covariates)** row | `R/area_level_comparison.R` |

Effect on the Ghana child-iron table, which previously showed only the three
model rows:

| approach | MAE pp | Pearson r | r_share |
|---|---|---|---|
| Individual SL | 9.52 | 0.748 | 1.11 |
| Area SL (MSE) | 11.35 | 0.590 | 0.88 |
| Area SL (NLL) | 13.97 | 0.544 | 0.81 |
| **National mean (null)** | **15.23** | — | — |
| **Spatial only (no covariates)** | **12.59** | 0.558 | 0.83 |

Two things this now makes visible that were invisible before: the covariate
models beat the no-model baseline by ~4 pp of MAE here (they earn their keep in
this cell), and the spatial-only row is within 0.03 of the Area SL, so most of
what the 500-covariate model is doing is reproducible from geography alone.
Note also `r_share = 1.11` for the individual SL — exceeding the estimated
ceiling, which is a useful alarm: either deff = 1.5 is too aggressive for this
cell or that model's CV is optimistic.

## 17. The Admin-1 layer — my own recommendation was wrong

`R/17_admin1_layer.R`, `out/admin1_vs_admin2.csv`

Recommendation 10 said to publish an Admin-1 layer because the ceiling roughly
doubles. The ceiling does double. **The predictions get worse anyway.**

| level | median units | r_max | ρ (curated + lon/lat RF) | share of ceiling |
|---|---|---|---|---|
| Admin-1 | 27 | 0.571 | 0.064 | **0.11** |
| Admin-2 | 75 | 0.349 | 0.137 | **0.35** |

Admin-1 wins in only **5 of 16** paired country × outcome cells. The reason is
the trade-off §12 flagged and I then under-weighted: Admin-1 has more signal per
unit but only 16–32 units to learn the covariate relationship from, and Gambia
(6 regions) and Sierra Leone (4) cannot be fitted at all. Model error swamps the
reliability gain.

**Revised position.** Admin-1 is the better level at which to *report a direct
survey estimate* — the survey supports it and Admin-2 it does not. It is the
worse level at which to *fit a covariate model*. Those are different uses, and
conflating them is what made recommendation 10 wrong. Fit at Admin-2; report
Admin-1 direct estimates alongside, and label them as design-based rather than
modelled.

## 18. Data-adaptive predictor selection vs the pre-specified 16

Same script. Selection is done strictly **inside** each training fold, with a
deliberately leaky variant that selects once on all the data for contrast.

Paired against `curated` (the 16 a-priori constructs):

| | Admin-1 (16 cells) | Admin-2 (26 cells) |
|---|---|---|
| adaptive lasso | +0.087 (t = 1.38, 9/16) | −0.007 (t = −0.33, 13/26) |
| adaptive top-20 by \|r\| | +0.022 (t = 0.30, 8/16) | −0.012 (t = −1.03, 11/26) |
| **LEAKY top-20 (selected once on all data)** | **+0.320 (t = 4.34, 14/16)** | **+0.147 (t = 4.21, 21/26)** |

**Two answers, and the second matters more than the first.**

*(To be precise about what the leaky row measures, since the wording invites a
misreading: there is NO leakage when selection happens inside the folds. The
honest adaptive rows above are unbiased. `LEAKY_top20_xy` selects once on all
the data and then "cross-validates" the fixed set -- it is a deliberate control
showing the cost of the common mistake, not a property of adaptive selection.)*

Letting the data choose does **not** reliably beat the pre-specified set. At
Admin-2 it is slightly worse; at Admin-1 the lasso is nominally ahead but at
t = 1.4 over 16 cells that is not a result. Given the pre-specified set is also
interpretable and stable across folds, there is no case for switching.

The selection-leakage effect is **larger than any real modelling difference in
this entire review**: selecting features once on the full data inflates Spearman
by **+0.30 at Admin-1** and **+0.16 at Admin-2** (t = 6.9 and 4.9). At Admin-1
the leaky number is ρ = 0.384 against an honest 0.15 — it would read as "letting
the data pick variables more than doubled performance", and all of it is
illusory. Any variable-selection step in this project must sit inside the
resampling loop, and any past result where screening happened first should be
re-checked (§3 already flags `spatial_plus_soil` as one such case).

Two cautions that survive correct in-fold selection:

* **Winner's curse on model choice.** Picking `curated16_rf_xy` as "best" after
  comparing ~20 variants on the same CV makes that headline optimistic. The
  PAIRED comparison against the production recipe (+0.063, t = 4.9) is sound
  because both pipelines were fixed in advance; "best of 20" is not.
* **Outcome-free screens are exempt.** `decorr_reps` clusters on X only, so it
  is safe anywhere. `screen_topK` and the lasso use y and must stay in-fold.

## 19. How accurately do the LOCO models recover NATIONAL prevalence?

`R/24_national_recovery.R`, `out/national_recovery.csv`

Aggregating each LOCO model's district predictions to a national figure and
comparing with the survey's own design-based national estimate — the quantity
these surveys are actually powered for:

| model | mean abs error | median | worst |
|---|---|---|---|
| PROD enet30 | **5.4 pp** | 4.2 | 15.6 |
| zscore ridge (best at spatial ranking) | 7.9 pp | 8.2 | 17.8 |
| centered_own ridge | 12.1 pp | 8.8 | 39.6 |
| train-country mean (null) | 12.3 pp | 10.0 | 31.2 |
| spatial_tps | 18.1 pp | 11.1 | 47.2 |

**Not accurately enough to use.** The best model misses by ~5 pp on average
against national prevalences of 1–38%. Per cell, iron is off by factors of 2–3:

| outcome | country | design-based | LOCO predicted |
|---|---|---|---|
| child_iron | Malawi | 10.6% | 28.4% |
| women_iron | Ghana | 7.6% | 19.0% |
| women_iron | Gambia | 28.2% | 10.6% |
| child_iron | Gambia | 37.5% | 24.8% |

Only women's vitamin A is close, and only because every country sits at 1–2%.

Note the ordering flips against §6: `zscore ridge` transports the spatial
*pattern* best but recovers the national *level* worse, because it discards the
level by construction and rebuilds it from the training countries. That is the
level/pattern split from §6 showing up again — a model can be good at one and
bad at the other, and a single metric hides it.

## 20. What follows from the surveys not being designed for Admin-2

Admin-2 is an **unplanned domain**: the weights are built to make national (and
usually Admin-1) estimates unbiased, and whatever sample lands in a district is
incidental. Three consequences, in order of how much they should change what you
do:

**(a) National prevalence should not come from the model at all.** The survey's
design-based national estimate is unbiased and precise (n in the thousands). Any
model can only add error to it — §19 shows the best LOCO model adding ~5 pp. So
the model's job is *disaggregation*: distributing a known national total across
districts, not estimating the total. This is the same conclusion the LOCO
anchoring result reached from the other direction (§6), and it argues for
**benchmarking**: constrain the Admin-2 predictions so their population-weighted
aggregate reproduces the design-based national (or Admin-1) figure. That is
standard SAE practice and it is not currently done.

**(b) Weights do not guarantee representativeness within a district.**
Calibration operates at the design strata (typically region × urban/rural), not
at Admin-2, so a district whose sampled clusters happen to be urban-skewed is
biased, not merely noisy — and the `r_max` ceilings in §1 account only for
*variance*, so they are if anything optimistic. Model-based SAE fixes the
variance part by borrowing strength; it does **not** fix a domain-level bias,
and Fay–Herriot will happily shrink toward a model fitted on biased direct
estimates.

**(c) Which districts got sampled is, encouragingly, not badly selective.**
Standardised mean differences between surveyed and unsurveyed districts on the
curated covariates (`out/surveyed_representativeness.csv`):

| country | surveyed / unsurveyed | max \|SMD\| | vars with \|SMD\| > 0.25 | worst |
|---|---|---|---|---|
| Ghana | 75 / 185 | 0.33 | 1 of 18 | soil zinc |
| Malawi | 86 / 153 | 0.24 | 0 of 18 | soil carbon |
| Gambia | 30 / 7 | **0.96** | 8 of 18 | NDVI |
| Sierra Leone | 14 / 0 | — | — | all surveyed |

For Ghana and Malawi — where most districts are unsurveyed and the extrapolation
actually matters — the surveyed districts are close to exchangeable with the
unsurveyed ones on observables. That is genuine support for the
predict-unsurveyed-districts use case. Gambia is badly non-exchangeable but has
only 7 unsurveyed districts, so little rides on it.

## 21. The missing GEE layers — fetched, not delegated

`sandbox_parsimony/python/fetch_tanzania_rasters.py`

The Earth Engine API works from this machine: the refresh token already in
`~/.config/earthengine/credentials` authenticates non-interactively, so no
student time is needed. `earthengine-api` was installed into a scratch venv
rather than any existing conda environment.

Twelve GeoTIFFs (~38 MB total) were written to `data/Tanzania_GEE_rasters/`,
which is the directory the pipeline already looks in — so this fixes the cause
rather than working around it. Tanzania now takes the raster path like every
other country instead of the legacy-CSV fallback.

Effect on Tanzania's own vocabulary: **102 → 176 GEE variables, 20 → 27
families**. Effect on the five-country intersection that the pooled and LOCO
models actually use:

| | stems | families |
|---|---|---|
| 4-country common (the achievable target) | 224 | 35 |
| 5-country common, before | 102 | 20 |
| **5-country common, after** | **139** | **26** |

Six families recovered for *every* country: **accessibility, dailyevi,
lai8days, landcoverlayers, productivity, terraclimate**. That closes about 30%
of the gap to the 4-country target.

Three production bugs had to be fixed to make this land, and each would have
silently wasted the download:

* `.append_gee_zonal_cols()` strips the country name out of derived variable
  names using a **hardcoded list that did not contain Tanzania** — every new
  column would have carried `_tanzania_` and failed the exact-name intersection.
  Now driven off `get_country_configs()` so adding a country cannot leave it
  behind.
* `extract_gee_admin2()` treated rasters and the legacy-parity CSV as
  alternatives. Giving Tanzania rasters would have **discarded its 134
  CSV-derived columns**, shrinking its vocabulary from 102 to 42. Now takes the
  union.
* the water-body and duplicate-key fixes from §14b, without which the new
  extraction reproduced the lake polygons.

What is still missing for Tanzania: `esa` WorldCereal, `fldas`, `ghsbuilts`,
`ghspop`, `wapor`, `landcovertype`, `aerosoloptical`, `atmosphere`, and an
`lst` naming that matches the other countries' year-position convention. Those
need the original export script (mn-proxies repo) to reproduce exactly; the
same API route works for them.

**Caveat on years.** Tanzania's survey is 2010, so annual composites use 2010 —
but Copernicus land cover only starts in 2015 and the MAP accessibility surface
is 2015, while the other countries' equivalents are 2016–2019. The filenames
carry the year so the mismatch stays visible rather than buried.

## 22. National benchmarking, implemented

`R/benchmark_area.R`, wired into `fit_area_level_model()` behind
`AREA_BENCHMARK_NATIONAL=true`

A single shift on the **logit** scale, solved so the population-weighted
aggregate of the district predictions equals the design-based national
prevalence. Chosen over the two textbook alternatives deliberately: ratio
benchmarking can push predictions above 1, difference benchmarking below 0 and
distorts low-prevalence districts hardest, while the logit shift is strictly
monotone — so it moves **only the level**, which is the component the model is
bad at, and leaves the ranking, which it is comparatively good at, untouched.

Measured over 24 country × outcome cells with population weights from the GEE
population layer:

| | before | after |
|---|---|---|
| national error vs design-based estimate | 4.33 pp | **0.00 pp** |
| district MAE vs direct estimates | 9.99 pp | 10.55 pp |
| mean \|change in district Spearman\| | — | **0.0000** |

End-to-end check on Ghana child iron: aggregate 19.1% → 20.2% (the design-based
figure), Spearman(before, after) = 1.000000.

**The district MAE gets slightly worse, and that is the point.** The survey's own
district estimates aggregate **4.50 pp** away from the national figure the *same
survey* produces — almost exactly the model's unbenchmarked drift of 4.33 pp. So
most of that drift is inherited from the direct estimates, not manufactured by
the covariate model, and forcing agreement with the trustworthy national number
necessarily increases disagreement with the untrustworthy district ones. That
4.50 pp gap is the cleanest single measurement of the unplanned-domain problem
in this whole review.

Left off by default so existing outputs do not change silently. `svy_admin2`
now carries the design-based national prevalence as an attribute, computed in
`compute_svy_admin2()` where the individual data is already in hand.

## 23. An honest national model, with no benchmarking — the VMNIS panel

> **Numbers superseded by §33.** Computed under complete-case with no
> imputation, which discarded a third of the panel. Conclusions unchanged.

`R/30_vmnis_national.R`, `R/31_national_model.R`

If the goal is a national prevalence for a country with no survey, aggregating a
district model is the wrong shape: it is trained on 4 countries to resolve
*within*-country pattern. Better to make the country-survey the unit of
observation. WHO's VMNIS supports that.

`results/external/vmnis_national_full.csv` holds only the 8 countries the puller
was pointed at. The full long-format dump (522 MB, already on disk) has **187
countries**, of which **102** have nationally representative surveys with a
deficiency prevalence. Usable panels:

| micronutrient | population | surveys | countries |
|---|---|---|---|
| Vitamin A | preschool children | 528 | **70** |
| Zinc | preschool children | 219 | 30 |
| Folate | non-pregnant women | 184 | 28 |
| Vitamin B12 | non-pregnant women | 128 | 21 |

*(One trap: `Representativenessname` names the SUBNATIONAL unit and is blank
exactly when a survey is national — filtering on it is backwards and returns 23
rows from 1 country. The `Representativeness` code is the right column.)*

Joined to 17 a-priori World Bank indicators (food supply, poverty, market
access, health system, demography), evaluated **leave-one-country-out** so the
model has never seen the held-out country's prevalence at any date:

| model | MAE | RMSE | Pearson | Spearman |
|---|---|---|---|---|
| random forest | **11.60 pp** | 16.59 | 0.187 | 0.231 |
| ridge | 11.99 pp | 16.76 | 0.147 | 0.136 |
| global mean (null) | 13.06 pp | 18.29 | — | — |

Best panel, vitamin A in preschool children (75 surveys, 53 countries):
**MAE 12.13 pp, r = 0.648**.

**Read that against the ceiling, not against zero.** Applying the §1 reliability
logic at national level — using countries with repeat national surveys of the
same panel to measure how repeatable a national prevalence is:

| panel | within-country SD (logit) | λ | r_max |
|---|---|---|---|
| Vitamin A / preschool | 1.12 | 0.44 | **0.66** |
| Zinc / preschool | 0.69 | 0.59 | 0.77 |
| Folate / non-pregnant women | 0.74 | 0.81 | 0.90 |

Vitamin A preschool achieves **r = 0.648 against a ceiling of 0.66 — 98% of what
is achievable**. The model is close to saturated; the 12 pp MAE is the target
being noisy, not the covariates being weak. Repeat national surveys of the same
country disagree by ~1.1 on the logit scale (roughly a factor of 3 in odds),
mixing genuine temporal change with assay, cutoff and inflammation-adjustment
heterogeneity — VMNIS records "Not Specified" adjustment for 290 of 528 vitamin
A surveys. This is the cross-survey level-offset problem from §6, at global
scale.

**Honest bottom line for unbenchmarked national prediction.** Two routes, and
they are not equivalent:

* aggregate the Admin-2 LOCO model → ~5.4 pp mean absolute error, but only
  demonstrated across 4 similar countries with harmonised biomarker definitions;
* the VMNIS national model → ~11.6 pp, but across 53 countries with no
  harmonisation and near the achievable ceiling.

The first number is smaller mostly because its test set is easier. Neither is
accurate enough to replace a survey; the VMNIS route is the more honest estimate
of what "predict a country we have never surveyed" actually costs, and its
correlation (r ≈ 0.65 for vitamin A) is good enough to rank countries even
though the level is not good enough to quote.

## 24. Modelling the VMNIS survey heterogeneity, and what predicts nationally

`R/32_vmnis_heterogeneity.R`

**The raw method comparison is confounded.** Marginally, unadjusted surveys
report 12.9% vitamin A deficiency and inflammation-adjusted ones 21.4% — but
countries that ran adjusted surveys are not a random sample, so the contrast
mixes method with setting and era. The method effect is only identified *within*
a country, so it is estimated with country fixed effects.

**Removing method effects raises the ceiling.** `r_max` is defined by how much
repeat national surveys of the same country disagree; taking the estimated
method and biomarker effects out of that disagreement first:

| panel | surveys | countries (repeated) | within-SD raw → adj | **r_max raw → adj** |
|---|---|---|---|---|
| Vitamin A / preschool | 115 | 69 (30) | 1.39 → 1.23 | **0.505 → 0.633** |
| Zinc / preschool | 41 | 30 (7) | 0.69 → 0.69 | 0.767 → 0.767 |
| Folate / non-pregnant women | 38 | 28 (9) | 0.83 → 0.83 | 0.879 → 0.879 |

For vitamin A — the only panel with enough repeated countries to identify the
contrast — roughly **a quarter of the apparently-unexplainable within-country
variance is survey method, not real change or sampling**. Zinc and folate have
7 and 9 repeated countries, too few for the fixed effects to move anything;
their unchanged rows are "no information", not "no effect".

*(Note `r_max` for vitamin A is 0.505 here versus 0.66 in §23. The difference is
the grouping: §23 averaged to one row per country-year, this groups by
country-year-method-biomarker and so retains more within-country variance. The
ceiling estimate is sensitive to how repeat surveys are collapsed, which is
worth stating whenever it is quoted.)*

**But modelling heterogeneity barely improves prediction:**

| panel | wdi only | + heterogeneity terms | standardised prediction |
|---|---|---|---|
| Vitamin A / preschool | ρ 0.549, MAE 12.22 | **ρ 0.590**, MAE 12.14 | ρ 0.528, MAE 12.36 |
| Zinc / preschool | ρ 0.304, MAE 14.58 | ρ 0.296, MAE 14.49 | ρ 0.297, MAE 14.66 |
| Folate / NPW | ρ 0.652, MAE 13.32 | ρ 0.630, MAE 13.58 | ρ 0.633, MAE 13.45 |

That is the expected result once stated plainly: you cannot know what method a
*future* survey will use, so knowing the method effect does not help predict the
number that survey will report. **The payoff is a better-defined target, not a
better fit.** If both the training outcomes and the prediction are standardised
to one method, the quantity being predicted is comparable across countries and
its ceiling is higher (0.633 vs 0.505). The coherent next step is to residualise
the *outcome* — out-of-fold — and evaluate against the standardised outcome,
which is not done here.

### Which national predictors matter

Permutation importance, averaged over the three panels (relative to the top
variable in each):

| rank | variable | mean rel. | rank | variable | mean rel. |
|---|---|---|---|---|---|
| 1 | under-5 mortality | 0.52 | 7 | agricultural land % | 0.30 |
| 2 | electricity access | 0.50 | 8 | year | 0.30 |
| 3 | food production index | 0.50 | 9 | sanitation | 0.30 |
| 4 | life expectancy | 0.43 | 10 | basic water | 0.24 |
| 5 | health expenditure pc | 0.38 | … | | |
| 6 | GDP per capita | 0.31 | **20** | **undernourishment prevalence** | **0.07** |

Two things stand out:

* **The panels split along interpretable lines.** Vitamin A is predicted by a
  child-health composite (under-5 mortality, life expectancy, water,
  electricity); folate by the food system (food production index, agricultural
  land, GDP); zinc by infrastructure (electricity, sanitation, food production).
  That is a coherent story, and a more defensible one than the Admin-2
  importance table where soil zinc topped outcomes it has no mechanism for
  (§11).
* **"Prevalence of undernourishment" — the headline food-security indicator —
  ranks last of 20.** Calorie sufficiency is not micronutrient sufficiency, and
  the data say so plainly. Broad development indicators beat the nutrition-
  specific one.

## 25. Predicting the method-standardised national prevalence (step 1)

> **Numbers superseded by §33.** Same complete-case caveat; the standardisation
> gain is larger once no surveys are dropped.

`R/33_vmnis_standardised.R`, `out/vmnis_standardised.csv`

§24 established that method terms as *predictors* barely help, because you cannot
know what method a future survey will use. So the correction is applied to the
**outcome** instead: for each held-out country, method and biomarker offsets are
estimated on the training countries only (with country fixed effects), then
subtracted from both the training outcomes and the held-out country's outcome.
That uses the held-out survey's method *label* — observable metadata — never its
prevalence.

| panel | target | MAE | Pearson | Spearman | folds identified |
|---|---|---|---|---|---|
| Vitamin A / preschool | raw (control) | 12.20 pp | 0.622 | 0.548 | 53 |
| Vitamin A / preschool | **standardised** | **11.04 pp** | 0.523 | 0.542 | 53 |
| Zinc / preschool | either | 14.58 pp | 0.222 | 0.304 | **0** |
| Folate / NPW | either | 14.36 pp | 0.576 | 0.540 | **0** |

**Vitamin A: MAE falls 1.16 pp with rank agreement unchanged.** Pearson drops,
which is the expected consequence of standardising rather than a loss —
removing method-induced spread leaves a less variable target, so the same
absolute accuracy yields a lower correlation. The two rows are scored against
**different targets**, so MAE is the meaningful comparison and r is not
comparable across them.

**Zinc and folate could not be corrected at all**: `folds_identified = 0` means
no training country had surveys under more than one method, so the contrast is
unidentified and the offsets are zero by construction. Their identical rows are
"no information", not "no effect" — the same caveat as §24.

## 26. National benchmarking turned on (step 2)

`AREA_BENCHMARK_NATIONAL` now defaults to **true** in `fit_area_level_model()`,
and `compare_admin2_approaches()` reports the level and the pattern as separate
columns rather than collapsing them into one RMSE.

The comparison table for Ghana child iron now reads:

| approach | MAE pp | Pearson r | r_share | national err pp |
|---|---|---|---|---|
| Individual SL | 9.52 | 0.748 | 1.11 | +3.71 |
| Area SL (MSE) | 11.35 | 0.590 | 0.88 | +2.04 |
| Area SL (NLL) | 13.97 | 0.544 | 0.81 | −4.80 |
| National mean (null) | 15.23 | — | — | — |
| Spatial only (no covariates) | 12.59 | 0.558 | 0.83 | — |

with `national_design_pp = 20.2` carried alongside. The unbenchmarked vector is
always retained as `area_pred_prev_unbenchmarked`, so nothing is lost and the
older behaviour is one environment variable away.

## 27. Finishing the Tanzania extraction (step 3)

`python/fetch_tanzania_rasters2.py`, `python/fetch_tanzania_rasters3.py`,
`R/merge_tanzania_chunks.R`

Tanzania's contribution to the shared vocabulary, measured against the
4-country target of 224 stems / 35 families:

| | stems | families |
|---|---|---|
| original | 102 | 20 |
| after pass 1 (§21) | 139 | 26 |
| **after pass 2** | **176** | **31** |

Five more families recovered for every country: `esa`, `fldas`, `ghspop`,
`landcovertype`, `lst`. **79% of the gap closed.**

Four things had to be got right, each of which would otherwise have silently
wasted the download:

* **ESA WorldCereal seasons live in a flat collection** filtered on a `season`
  property, not as per-season collections — the obvious path returns "asset not
  found" for every season.
* **Earth Engine caps a direct download at 50 MB.** FLDAS (143 MB) and the
  12-month LST stack (90 MB) had to be fetched in band chunks and merged.
* **Band names do not survive `writeRaster` to a plain GeoTIFF**, and
  `.append_gee_zonal_cols()` derives one column per band *from the band name*.
  Left unset, FLDAS arrived as `gee_fldas_..._chunk_0_1` and matched nothing.
* **GHS population carried a nodata fill** that averaged to −10.96 people per
  cell until masked.

Still outside the intersection: `ghsbuilts` and `wapor`. WAPOR is not in the
public Earth Engine catalogue at all — MODIS NPP stands in under the WAPOR
filename so the column name matches, which means **that column is a different
product for Tanzania than for the other four countries** and should be treated
as such. ESA WorldCereal band names carry the AEZ tile id (Ghana's read
`32121_TC-MAIZE-MAIN_...`), so they are country-specific by construction and
will never match across countries however they are fetched.

## 28. Re-scoring `spatial_plus_soil` (step 4)

`R/34_rescore_spatial_plus_soil.R`, `out/spatial_plus_soil_rescored.csv`

The eight locked SoilGrids features in `benchmarks_all.csv`'s leading method were
ranked by LOCO r on the same four outcomes and four held-out countries the method
is scored on. Moving that selection inside the loop — ranking candidates on the
training countries only, then fitting the identical model:

| variant | ρ | Pearson | RMSE pp |
|---|---|---|---|
| locked features, as published | 0.316 | 0.332 | 18.3 |
| spatial only, no soil | 0.274 | 0.286 | 16.5 |
| **features selected in-fold** | **0.117** | 0.143 | 18.2 |

*(The locked reproduction lands at 0.316 against the published 0.305, which
validates the setup.)*

* **Optimism from locking on the test metric: +0.197** (SE 0.104, t = 1.89,
  better in 10 of 14 cells).
* **Honest value of the soil block over pure spatial: −0.172** (t = −1.62,
  better in only 4 of 14 cells) — the soil features make things *worse*.

The published lead of `spatial_plus_soil` over `spatial_coords` is entirely a
selection artefact. It should be withdrawn from the headline comparison, and
`spatial_coords` treated as the spatial benchmark.

## 29. The hierarchical individual-biomarker model (step 5)

`R/35_hierarchical_biomarker.R`, `out/hierarchical_biomarker.csv`

The last untested route to the continuous target's extra signal: fit
`log(biomarker)_ij = X_j β + u_j + e_ij` to individuals, then derive district
prevalence. Its promised advantage was σ estimated from every individual in the
country rather than the ~13 in one district.

| route | MAE pp | Pearson | Spearman |
|---|---|---|---|
| indicator (production) | **10.83** | 0.207 | 0.225 |
| two-step, mean + pooled SD | 12.58 | 0.280 | **0.279** |
| hierarchical individual | 15.51 | 0.207 | 0.187 |

**It loses to both simpler routes.** Two bugs had to be fixed before that
verdict was trustworthy, and both are worth recording:

* For an **unsurveyed** district the predictive spread is `√(σ² + τ²)`, not σ —
  integrating the random effect out of the mean but not the variance makes every
  district look alike. Correcting it was worth ~10 pp of MAE.
* **`lmer` treats `weights` as precision, and the survey weights are not on a
  common scale** — Malawi's median is 790,000 (raw DHS 1e6 convention) against
  Ghana's 0.98. Each Malawian was being counted as ~790,000 observations. This
  affects any use of these weights other than a weighted average, and Malawi is
  8 of the 18 cells here.

One caveat on the verdict: the hierarchical model's mean structure is linear
(`lmer`) while both competitors use a random forest, so the comparison is not
clean. Its better σ did not compensate for a weaker mean model, and the mean
structure is what drives ranking. A fair rematch would put a forest in the mean.

**Standing recommendation is unchanged:** use the two-step route only for
outcomes where the indicator model has no signal, where §16 measured it at
+0.09 ρ.

## 30. Suppressing zero-ceiling cells (step 6)

`area_comparison_all` now blanks `pearson_r` and `r_share` for country × outcome
cells whose `signal` flag is FALSE — that is, where even the *optimistic* bound
on r_max is below 0.15, so the entire between-district spread in the survey
estimate is explainable by sampling noise. Raw values are preserved in
`pearson_r_raw` / `r_share_raw`, and the count is logged.

The point is the one from §1: reporting r = −0.04 for Ghana women's vitamin A
reads as a model failure when it is a data limit, and the table was presenting
that identically to Ghana child iron sitting at 91% of its achievable ceiling.

## 31. Reconciling the two VMNIS ceilings — and a correction to §23

`R/36_vmnis_ceiling_resolution.R`, `out/vmnis_ceiling_resolution.csv`

§23 reported r_max = 0.66 for vitamin A and §24 reported 0.505 for the same
panel. Nothing about the data differed — only how repeat surveys were grouped
before the within-country variance was taken:

* **one row per country-year** averages method differences away, counting them
  as *signal* → 0.663
* **one row per country-year-method** retains them, counting them as *error*
  → 0.572

**This is a VMNIS-only issue.** The Admin-2 ceiling uses the binomial sampling
variance and never involves repeat surveys, so nothing there is affected.

### The fix: estimate the components, do not pick a grouping

VMNIS records `Samplesize` for ~85% of surveys, so the sampling term can be
computed directly — exactly as at Admin-2 — instead of being bundled into a
residual. Everything else goes into a variance-components model:

`logit(prev) = country + method + year trend + residual + sampling`

Vitamin A, 528 surveys, 70 countries (logit-scale SDs):

| component | SD | what it is |
|---|---|---|
| country | **1.411** | between-country signal — what a national model predicts |
| method | 0.564 | inflammation adjustment + biomarker |
| residual | 0.550 | real temporal change + unexplained |
| sampling | 0.479 | binomial, from `Samplesize` |

Two ceilings follow, for two genuinely different estimands:

| panel | legacy pair | **predict what a survey reports** | **predict a standardised prevalence** |
|---|---|---|---|
| Vitamin A / preschool | 0.572 – 0.663 | **0.837** | **0.888** |
| Zinc / preschool | 0.767 | 0.892 | 0.960 |

The gap between the two columns is exactly what §25's outcome standardisation
buys, now stated as a number rather than appearing as a side effect of a
grouping choice. Folate drops out — too few countries for the components to be
identified.

### The correction

§23 said the national model was at "98% of achievable" (r = 0.648 against
r_max = 0.66). That rested on one particular estimator. Across defensible
estimators the vitamin A ceiling spans **0.57 to 0.89**, so the same model is
somewhere between *at* the ceiling and at *73%* of it.

Two reasons the variance-components ceilings are higher: a secular **year trend
is removed** (a predictable trend is not noise), and **method variance is
separated** rather than lumped into the error. Both are legitimate for a model
that includes year and can standardise method — but they are choices, and they
move the answer a lot.

**So: the ceiling should always be quoted with its estimator and estimand, never
as a bare number.** The honest summary for the national model is that it
captures most of the between-country signal and that the remaining headroom is
uncertain, not that it is 98% saturated.

## 32. The lake clusters are displaced GPS, not bad data

`snap_water_to_land()` in `R/admin2_key_hygiene.R`

§14b treated the survey rows assigned to `Lake Malawi` and `Lake Victoria` as
errors and dropped them. That was the wrong call. Checking the actual
coordinates against the polygons:

| country | cluster coordinate | inside lake polygon | nearest land district | distance |
|---|---|---|---|---|
| Malawi | 34.2052, −11.2322 | yes | TA Musisya | **0.00 km** |
| Malawi | 34.2595, −11.6653 | yes | TA Timbiri | **3.84 km** |
| Malawi | 34.3092, −11.6263 | yes | TA Timbiri | **7.95 km** |
| Tanzania | 31.7598, −2.6308 | yes | Chato | **0.36 km** |

Four distinct cluster locations, every one within the DHS/MICS confidentiality
displacement bounds (2 km urban, 5 km rural, 10 km for a random 1% of rural
clusters). These are **real lakeside communities whose GPS was displaced into
the water**, not people living on a lake and not a coding error.

Dropping them discards 87 individuals in Malawi and 175 in Tanzania — the
Tanzanian one is among that survey's larger clusters. They are now **snapped to
the nearest land Admin-2 polygon** instead, in the `merged_*` target so every
downstream consumer sees the correction:

* Malawi 87 rows → TA Musisya (27), TA Timbiri (60)
* Tanzania 175 rows → Chato (175)

A `max_km = 12` guard leaves anything beyond the displacement limit alone rather
than moving it to a district it may have nothing to do with; nothing hit that
guard here.

Water polygons are **still dropped from the covariate and prediction side** —
a lake is not a populated prediction unit, and its zonal statistics are NDVI and
soil chemistry of open water. Only the survey side snaps.

---

## 33. Rerunning the national model on the project's own covariate panel

`R/37_national_covariates.R`, `R/31_national_model.R`, `R/33_vmnis_standardised.R`

§23 used 17 a-priori WDI indicators pulled from the World Bank API, because 31
was written as a self-contained viability check and reached for its own
covariates. That was never a considered parsimony decision: the project already
holds a national country x year panel
(`mn-proxies/data/national/combined_dataset.dta`) and nothing was reading it.
Two questions, tested together — does the wider pool help, and was the
complete-case rule costing more than the covariates were worth.

### What the panel actually holds

5,433 columns is a misleading headline. It is a full country x year grid, but
most columns are DHS STATcompiler indicators that exist only in DHS survey
years, so the median column is **5.5% populated** over 1990-2022. Carrying the
nearest value within +/- 5 years lifts that to 104 columns at 90% coverage, 280
at 70%, 1,375 at 50%. 61 countries, 57 of which map to ISO3.

Thirteen measured anaemia/haemoglobin columns are dropped by default —
predicting a micronutrient deficiency from a measured deficiency is a different
and much weaker claim. Supplementation and dietary columns ("received iron
tablets", "vitamin A-rich foods", "micronutrient powder") are kept: those are
interventions and exposures, and mechanistically among the most relevant proxies
in the file. The distinction is only visible in the Stata variable labels — the
column names are opaque codes (`AN_ANEM_C_ANY`, `ML_HEMO_C_HL8`), so a
name-based filter silently keeps all of them.

### Three changes to the fitting, all fold-internal

* **Surveys are no longer dropped for missing covariates.** Complete-case now
  applies to the outcome only. Against ~1,600 candidate columns the old rule
  would have discarded essentially every survey.
* **KNN imputation** (k = 5; distance is mean squared difference over the
  dimensions two rows both observe, which is the only usable definition at this
  sparsity). Donors and the standardisation defining the distance come from
  training rows only.
* **Missingness indicators** for columns whose training missingness falls
  between 2% and 98%. "This country never reported this indicator" is patterned
  by survey programme, era and capacity, so it is at least a candidate signal
  rather than something to paper over with a median.
* **Per-fold prescreen** by absolute Spearman correlation with the training
  outcome. With ~1,600 candidates and ~50 folds this is a screening problem, and
  the screen has to sit inside the fold or it sees the held-out country.

### Results: vitamin A, preschool children, leave-one-country-out

| covariates | countries | surveys | MAE pp | Pearson |
|---|---|---|---|---|
| 17 WDI, complete-case (§23 as published) | 53 | 75 | 12.13 | 0.648 |
| **17 WDI + KNN + indicators** | **69** | **115** | **11.75** | **0.655** |
| 17 WDI, common support | 26 | 35 | 16.64 | 0.398 |
| 1,587 panel columns, common support | 27 | 36 | 18.64 | 0.069 |
| union of both, common support | 26 | 35 | 15.47 | 0.441 |

"Common support" restricts to the 54 countries both sources cover, so the
covariate contrast is not confounded with the country sample. Without it the
panel is scored on a smaller, more African, higher-burden subset (mean
prevalence 39.9% against 24.9%) and the comparison means nothing.

And the §25 standardisation, on the same rerun:

| target | surveys | countries | MAE pp | folds identified |
|---|---|---|---|---|
| raw outcome (control) | 115 | 69 | 12.46 | 69 |
| **method-standardised** | 115 | 69 | **10.81** | 69 |

### Four things follow

**Not dropping surveys was the whole win.** 75 surveys became 115 and 53
countries became 69; every metric improved. That single change matters more
than the covariate question it was run to answer.

**The panel does not beat 17 hand-picked indicators.** On identical countries,
panel-only is clearly worse (18.64 vs 16.64 pp). The screen was swept from 10 to
150 predictors and the ordering never changed, so it is not a screening
artifact. The union is the best of the three, so the panel adds something on top
of WDI — it just cannot carry the model alone.

**Country coverage dominates covariate count.** Holding the covariates fixed at
the same 17 WDI variables, restricting to the panel's countries costs **4.9 pp
of MAE**. Nothing on the covariate side moved the answer by half that. Seventeen
variables covering ~200 countries beat 1,587 covering 57, and that — not
parsimony — is the reason the original choice was right.

**The missingness indicators never fired.** Not one of the ~1,400 panel
indicators survived the top-60 screen in any fold, while all 5-6 WDI ones were
used. DHS reporting patterns carry no rank information about deficiency here,
which was the main argument for keeping the sparse columns.

**Method standardisation needs the full country set.** Every common-support run,
and both the zinc and folate panels, report `folds_identified = 0`: no training
country had surveys under more than one method, so the contrast is unidentified
and the offsets are zero by construction. Their identical raw and standardised
rows are "no information", not "no effect" — the same caveat as §24 and §25.

### Correction to §23 and §25

Both sections' numbers were computed under complete-case with no imputation.
They remain accurate for that configuration; the current numbers are the table
above. The direction of every §23/§25 conclusion is unchanged — the national
model ranks countries usefully and cannot be quoted on level, and standardising
the outcome buys MAE while leaving rank agreement alone. The standardisation
gain is now larger (1.65 pp, from 1.16 pp) because it is measured on 69
identified folds rather than 53.

---

## What I would change, in order

1. **Report the noise ceiling.** `r_max` next to every correlation; suppress
   cells where it is indistinguishable from zero. Cheapest change, biggest
   effect on how the results read.
2. **Fix the Fay–Herriot sampling variance.** One-line guard; the EBLUP
   currently keeps ~71% of the raw direct estimate where it should keep ~26%.
3. **Add null and spatial-only rows to every results table**, and switch to
   repeated CV.
4. **Put geography in the model.** Adding lon/lat to the area-level predictors
   is a two-line change and is where most of the measured gain comes from.
5. **Adopt a parsimonious within-country predictor set** — 16 curated variables
   plus coordinates, in a random forest, beat the production recipe in 14 of 16
   cells. Better accuracy, an interpretable model, and variable importance that
   does not split across nine copies of the same raster.
6. **Move the LOCO fits off the logit scale** (§7). Within country the scale
   barely matters; for transport it is worth up to +0.17 Spearman.
7. **For LOCO, change the target, not the predictor count**: train on the
   within-country z-score and take the level from a national anchor. Report
   pattern and level error separately.
8. **Fix the Tanzania GEE extraction** or exclude it from the pooled covariate
   intersection; log the intersection size each run.
9. **Re-score `spatial_plus_soil` with features selected inside the LOCO loop**,
   or drop it from the headline comparison.

### Larger pieces of work, in expected-value order

10. ~~Publish an Admin-1 layer alongside Admin-2~~ — **tested, and I was
    wrong** (§17). The ceiling doubles but achieved accuracy halves, because
    16–32 regions is not enough to fit a covariate model on. Report Admin-1
    *direct survey* estimates alongside the Admin-2 map, labelled design-based;
    do not fit the model there.
11. **Benchmark the Admin-2 predictions to the design-based national estimate**
    (§19, §20). These surveys are powered for national, not Admin-2. The model
    should disaggregate a known national total, not re-estimate it — the best
    LOCO model misses national prevalence by ~5 pp on average and by a factor
    of 2–3 on iron. This is the single most consequential change left.
12. ✅ **Extend the covariate set on the diet-quality pathway** (§14, §21) —
    done. Travel time, MapSPAM crop mix and agro-ecology added for all five
    countries; six GEE families recovered for everyone by fetching Tanzania's
    missing rasters. Remaining: `esa`, `fldas`, `ghsbuilts`, `ghspop`, `wapor`,
    `landcovertype`, and an `lst` naming that matches the other countries.
13. **Model the district biomarker distribution** (§13, §16) — but not the way
    I proposed. Modelling the scale failed; what remains untested is a
    hierarchical model on individual biomarkers with district random effects.
    Use the biomarker route today only for outcomes where the indicator model
    has no signal, where it is worth +0.09 ρ.
14. **Re-extract covariates at true cluster coordinates** (§12) — worth much
    less than it sounds: the surveys carry only 60–105 GPS clusters, so the
    design-point count barely rises and per-unit reliability falls. Last.

### What NOT to do

* **Do not fit at individual or cluster level hoping for more power** (§12).
  Every proxy covariate is constant within Admin-2, so the effective sample size
  is the district count no matter what level the model is fitted at. "Fit at the
  finest common level" means **Admin-2** -- the finest AREAL level every country
  supports -- not cluster or individual.
* **Do not mix Admin-1 and Admin-2 across countries to equalise unit counts**
  (§15). It blends two different covariate-outcome relationships (14% of
  variables flip sign between levels), confounds level with country and with
  reliability, makes people-per-unit LESS even rather than more, and buys
  nothing empirically.
* **Do not select variables outside the resampling loop** (§18). Doing it once
  on the full data inflates Spearman by +0.30 at Admin-1 and +0.16 at Admin-2 --
  larger than any genuine modelling difference found anywhere in this review.
* **Do not reach for a more elaborate learner.** A binomial GAM with a spatial
  field, the correct likelihood on survey counts, and a thin-plate spline all
  lost to a random forest on 18 predictors (§5b). Model machinery is not the
  binding constraint.
* **Do not present soil zinc as a mechanism** (§11). It ranks first for vitamin
  A, B12 and folate, where it has none; it is a geographic gradient in disguise.
