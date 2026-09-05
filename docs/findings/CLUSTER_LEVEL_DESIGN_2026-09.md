# Modelling at the cluster level: linkage, feature engineering, and what the field does

*Design note, 4 September 2026. A brainstorm in response to the question "how
would GEE extractions get linked to clusters, and what should the features
be, without overwhelming the model?" Nothing here has been run; the facts
about our own data were checked on disk today.*

## 1. What we have, and a correction

**Cluster GPS exist for all four surveys**, undisplaced as far as the files
show (these are research surveys, not DHS releases):

| Survey | GPS file | Clusters | Keys carried |
|---|---|---|---|
| Gambia 2018 | `data/IPD/Gambia/Gambia_GMS_GPS_cleaned.csv` | 71 | LGA, district code, ward code |
| Ghana 2017 | `data/IPD/Ghana/Ghana_GMS_GPS_cleaned.csv` | 92 | region, EA code, district, urban/rural, EA population |
| Malawi 2015-16 | `data/IPD/Malawi/Malawi_GMS_GPS_cleaned.csv` | 105 | cluster number only |
| Sierra Leone 2013 | `data/IPD/Sierra Leone/Sierra Leone_GMS_GPS_cleaned.csv` | 60 | region, village |

**A cluster-buffer extraction already exists.** `data/GEE/DHS GEE merge.Rmd`
builds `sf::st_buffer` polygons of 5, 10, 25 and 50 km around cluster points
in a projected CRS and takes raster means inside them; its output for Ghana
(`data/GEE/GMS_gee_merge_GH_2017.csv`, 92 clusters) carries TRMM rainfall,
population, NDVI and temperature at 10/25/50 km. The current Admin-2
pipeline uses the same operation with district polygons instead of buffers
(`reduceRegions` in GEE, cached per layer under `data/GEE/.cache_*`). So the
linkage is not a new mechanism: **the buffer polygons replace the district
polygons as the feature collection sent to the same reducers.** Fieldwork
dates are now on disk for every cluster (`fieldwork_windows_cluster.csv`,
FW-01), so time-matching is possible at the cluster too.

**Correction to what I said earlier.** I described cluster-level modelling as
raising the effective n "from 146 units to several hundred per country".
The cluster counts above sum to 328 against 146 units on the consistent
district rung (208 at the Admin-2 rung). That is a 1.6-2.2x gain, not an
order of magnitude. The larger gains are elsewhere:

- **resolution of the covariate, not of the outcome.** In Ghana, Malawi and
  Gambia most districts hold one cluster, so the district outcome *is* a
  cluster outcome already, but the covariate is averaged over a district the
  cluster may not resemble. Extracting at the cluster removes that
  ecological averaging.
- **honest variance structure.** VC-01 showed 21% of the published
  reliability ceiling is cluster effect and that child zinc has no district
  variance at all. A cluster-level model carries the cluster as a level of
  its own instead of folding it into "district geography".
- **Sierra Leone stops being 14 rows.** Its 60 clusters sit in 14 districts;
  at the cluster level it contributes as much as Gambia.

## 2. How the field does the linkage

Three families of practice, all using the cluster point rather than the
district polygon:

1. **DHS Program geospatial covariate datasets** (Mayala et al. 2018). One
   value per cluster per layer, the mean inside a buffer of **2 km for urban
   and 10 km for rural clusters** (matched to the DHS displacement radii),
   about 30 layers, time-varying layers taken for the survey year, static
   layers once. No within-buffer SD, no seasonality columns. This is the
   most-used format in the nutrition-epidemiology literature and the one
   collaborators will recognise.
2. **Model-based geostatistics (IHME Local Burden of Disease, Malaria Atlas
   Project).** Covariates on a 5 x 5 km grid; the cluster takes the pixel
   value (or a small buffer). Dynamic layers are annual for the survey year
   (LBD: EVI, LST, rainfall, night lights, accessibility, urbanicity,
   population; Osgood-Zimmerman et al. 2018, Weiss et al. 2019). MAP goes
   further for infection outcomes: **synoptic (long-run) monthly means plus
   survey-month-matched anomalies with 0-3 month lags** (Bhatt et al. 2015),
   and seasonality as **temporal Fourier components** (annual mean,
   amplitude and phase; Scharlemann et al. 2008).
3. **Household-survey economics.** Cluster linked to the **nearest market**
   or an inverse-distance mean of the nearest few, with distance kept as a
   covariate; price features are relative and seasonal (below).

Our undisplaced GPS means the 2/10 km buffers are a choice, not a necessity.
The defensible default is the DHS convention (2 km urban, 5 km rural, with a
10 km rural sensitivity) because a cluster's respondents live within that
radius and it keeps us comparable with the literature. 25 km and 50 km
buffers, which the old Ghana merge also produced, average away exactly the
local variation the cluster level is meant to recover; keep at most one
"context" radius (25 km) for layers where the neighbourhood matters
(markets, health facilities, night lights), not for soil or rainfall.

## 3. Feature engineering: a recipe that stays small

The question was whether each layer needs a mean, an SD, a seasonality and a
deviation from an annual or national mean. Most of those are unnecessary for
this project, and the domain-PC step already absorbs the count within a
domain. The recipe:

| Layer type | Spatial summary | Temporal summaries | Columns per layer |
|---|---|---|---|
| Static (soil, elevation, AEZ class, land cover) | buffer mean; SD only for elevation (ruggedness) | none | 1 (2 for elevation) |
| Slow (population, built surface, night lights, travel time) | buffer mean (+ 25 km context for night lights and travel time) | survey year | 1-2 |
| Dynamic climate (LST day/night, rainfall, EVI, soil moisture) | buffer mean | climatology mean (2001-2020), seasonal amplitude and peak month, survey-year annual mean, **fieldwork-window anomaly with 0-3 month lag** | 5 |
| Food prices (markets) | IDW over 3 nearest markets + distance to nearest | seasonal amplitude, fieldwork-window anomaly and seasonal position, relative price of nutritious groups, staple USD/kg | 5-6 per group |

Rules that follow from our own results:

- **No "deviation from the national mean" column.** Within-country
  rank-normalisation (fix 3 of protocol v2) already removes the national
  level before pooling; an explicit deviation is the same information twice.
  Absolute values matter only for transporting *levels*, which do not
  transport anyway (LV-01, AR-01).
- **No within-buffer SD except terrain.** SD columns double the count and
  the ablation (DA-01) already shows capacity is a liability; the field
  does not use them.
- **Keep the two temporal roles apart.** Climatology, amplitude and phase
  are properties of a place and transport to an unsurveyed country.
  Fieldwork-window anomalies are properties of a survey; they exist only for
  surveyed clusters and belong in the target-adjustment step (remove the
  survey-timing artefact from the training outcome), or become predictors
  only under an explicit "the next survey will run in month M" scenario.
  Script 41 builds both kinds and labels them `sea_*` and `tmc_*` for that
  reason.
- **Seasonality as amplitude and phase, not twelve monthly columns.** The
  twelve `lst_night_m01..m12` columns already in the vocabulary are the
  expensive way to say "amplitude 6 degrees, peak in April".

At 2 km / 5 km buffers and this recipe, a cluster row carries roughly 40
static, 10 slow, 30 dynamic-climate and 15 price columns: about 100 columns
in eight domains, against 451 in the current district vocabulary.

## 4. Food prices specifically

What the literature codes, and the column each becomes (all built today in
script 38 for districts and for cluster points):

| Concept | Source | Column |
|---|---|---|
| Seasonal price gap for staples (peak-to-trough, log points) | Gilbert, Christiaensen & Kaminski 2017 | `fpt_staple_seas_amp` |
| Price anomaly at the time of the survey, z against the market's own seasonal norm | WFP ALPS; FAO IPA | `fpt_staple_anom_z` |
| Where in the seasonal cycle the survey fell | (as above) | `fpt_staple_seas_pos`, `fpt_months_to_peak` |
| Spatial price index: local vs national median, same commodity and unit, same months | standard | `fpt_*_rel_win` |
| Relative price of nutritious foods against the staple | Headey & Alderman 2019; Bai, Herforth & Masters 2022 | `fpt_animal_to_staple`, `fpt_pulses_to_staple`, `fpt_vegfruit_to_staple` |
| Absolute staple price, USD per kg | for cross-country level comparison only | `fpt_staple_usd_kg` |

Two cautions from building them. WFP series are thin in the survey years
(Ghana 2017 has about 20 observations per market), so many market-months are
empty and the anomaly falls back to the country profile; and script 07 took
Gambia's survey year as 2021 when the fieldwork was January-April 2018, so
the nine `fprice_*` columns already in the vocabulary are three years late
for Gambia. Script 39 scores the new block against them.

## 5. The model at the cluster level

- **Outcome.** Binomial: deficient count out of tested, with Kish effective
  n per cluster as the trials, or the survey-weighted cluster prevalence with
  precision weights; and the cluster mean of the negated log biomarker for
  the level target. Both as in protocol v2.
- **Structure.** Cluster nested in district nested in region; a district
  random intercept, or the spatial smoother, alongside the covariate index.
- **Folds.** Unchanged: the three estimands are still district in-fill,
  region extrapolation and country transport, so folds are cut by district
  and region, never by cluster. Clusters in the same district never split
  across folds.
- **Prediction to districts.** Two ways, and the difference matters. The
  cheap way predicts at the district's mean covariates, which is the
  ecological prediction we make now with a better-trained model. The right
  way, as in model-based geostatistics, predicts on a population grid and
  averages up, because a nonlinear model at the mean covariate is not the
  mean of the model over the district. Start cheap; the grid version needs
  gridded covariates the pipeline does not yet store.
- **What it does not fix.** The survey's sample size. AR-01 and G4-02 stand:
  a cluster-level model orders districts better; it does not replace
  respondents.

## 6. Order of work, if this is taken up

1. Buffer extraction for the 328 clusters at 2 km / 5 km through the
   existing `reduceRegions` path (one GEE run; the account and cache
   layout are in `gee_rgee_setup`); static and slow layers first.
2. Rebuild the climate block as climatology + amplitude + phase from the
   monthly layers already extracted, and the fieldwork anomalies from
   FW-01.
3. Fit the protocol-v2 arms at the cluster level under the same three
   estimands; report district-level scores next to the current ones.
4. Only then decide on grid prediction.

## 7. First pass (same day)

Steps 1-3 were run on 4 September from rasters already on disk
(`scripts/cluster_level/01-03`, no Earth Engine call): 323 clusters, 117
buffer covariates, the protocol-v2 arms under the three estimands, scored at
the cluster and after aggregation to districts. The cluster-fitted models did
not beat the district-fitted ones on any estimand and were worse for
transport (climate + soil index at Admin-2: 0.24 vs 0.37 on the level); the
fieldwork-window block added +0.02 in-country. The likely reasons and the
next things to try are in `SANDBOX_LOG_2026-09.md` (CL-01/02/03). The
resolution argument in section 1 stands as an argument; on this pass the
noise in a cluster of 8-19 respondents outweighed it.

A second pass (CL-04) weighted the index by Kish n and shrank cluster
outcomes toward their district by empirical Bayes. On prevalence the
estimated between-cluster variance within districts is zero, so the shrunk
outcome is the district mean itself; that variant still trails the district
fit, which means the 2-5 km covariates are not better predictors than
district means on this vocabulary. The repairs recovered a third to a half
of the gap for region extrapolation and Admin-1 transport, and nothing
in-fill. The cluster track is a sensitivity analysis until the vocabulary
is matched layer for layer and grid prediction is tried.
