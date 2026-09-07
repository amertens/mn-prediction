# Cluster-level analysis (parallel to the district and regional estimations)

Started 2026-09-04. Nothing here replaces the Admin-2 / Admin-1 work in
`scripts/protocol_v2/` or the targets pipeline; it is the same question one
rung down, at the survey cluster, where the outcome is measured. Design and
rationale: `docs/findings/CLUSTER_LEVEL_DESIGN_2026-09.md`. Not part of the
`{targets}` DAG.

| Script | Produces | Notes |
|---|---|---|
| `01_build_cluster_targets.R` | `results/tables/cluster_level/targets_cluster.csv` | 323 clusters, all with GPS; uniform outcome, negated log biomarker level, Kish n; fieldwork date from FW-01 |
| `02_extract_cluster_covariates.R` | `data/covariates/cluster/predictors_cluster.csv` (+ metadata) | Buffer means (2 km urban / 5 km rural, urban from GHSL SMOD) over rasters already on disk; no Earth Engine call. Dynamic layers get climatology, amplitude, peak month and fieldwork-window columns (`_fw`, `_fw_anom`, `_prev3`, role = fieldwork) |
| `03_cluster_benchmarks.R` | `results/tables/cluster_level/benchmarks_cluster*.csv` | Protocol-v2 arms fitted on clusters under the same three estimands (district in-fill, region extrapolation, country transport); scored at the cluster and after aggregation to districts, next to the district-level results |
| `04_weighted_shrunken.R` | `results/tables/cluster_level/benchmarks_cluster_ws_*.csv` | Kish-n-weighted index x empirical-Bayes shrinkage of cluster outcomes toward district/region (training fold only); recovers a third to a half of the gap to district fitting, never all of it |

Conventions carried over from protocol v2: within-country rank normalisation
before pooling, domain principal components, the zero-tuning domain index as
the estimator of record, precision weights from effective n, folds cut by
district and region (never by cluster).

Environment: `CL_URBAN_KM`, `CL_RURAL_KM` (buffer radii, default 2 / 5);
for `03`: `CL_REPS` (in-fill replicates, 10), `CL_SETS` (comma list of
`transportable`, `with_fieldwork`, `climate_soil`), `CL_ARMS`, `CL_TAG`
(output suffix for a partial run), `CL_SUMMARY_ONLY=1` (re-print the
comparison from existing tables).

## First pass, 2026-09-04

Cluster-fitted models aggregated to districts do not beat the district-fitted
protocol-v2 models on any estimand, and transport is worse (climate + soil
index at Admin-2: 0.244 vs 0.368 on the biomarker level). The
fieldwork-window block adds +0.02 in-country on the level target. Full table
and the list of things to try next (weighted index, shrunken cluster
outcomes, matched vocabulary, 10 km radius, grid prediction) in
`docs/findings/SANDBOX_LOG_2026-09.md`, entry CL-01/02/03.

Data notes: Sierra Leone's rasters are named both `Sierra_Leone` and
`Sierra Leone`; seven of its space-named soil GeoTIFFs crash GDAL part-way
through a read, so `02` uses the single-band `*_clean.tif` copies written by
`scratchpad/clean_sl_soil.R` (one process per file).


## Additions, 2026-09-07

| Script | Output | What it does |
|---|---|---|
| `00_export_gee_geoms.R` | `data/external_cache/gee_geoms/*.geojson` | Simplified Admin-2 polygons and cluster buffers (2 km urban / 5 km rural from the GHSL flag in `predictors_cluster.csv`) for the Earth Engine reducers |
| `05_merge_new_layers.R` | `predictors_cluster*.csv` (+ metadata) | Adds the IHME 5 km surfaces (validated columns), GLW4 livestock density, water / coast distance and the ESPEN helminth block to the cluster table |
| `06_matched_vocabulary.R` | `results/tables/cluster_level/matched_vocabulary_*.csv` | District-level protocol arms on the cluster vocabulary aggregated to districts, the full set, and both, next to the cluster fits aggregated to districts: separates the fitting unit from the layers |

Environment added: `CL_OUT_TAG` (suffix for an alternative extraction, e.g.
`_r10` for 2 km / 10 km buffers; honoured by `02`, `00`, `05`, script 48 and
the Earth Engine script), `CL_PRED_TAG` (which predictor table `03`, `04` and
`06` read), `CL_OUT_DIR` (where `03` and `04` write, so a tagged run does not
mix into the main summary, which globs every `benchmarks_cluster_cells*.csv`
in its folder), `CL_SKIP_ADMIN2=1` (Earth Engine script: clusters only), and
`GLW_YEAR` / `GLW_TAG` for script 48 (2020 FAO density rasters are the
default; 2015 Dataverse counts rank districts identically). The 10 km rural
sensitivity (CL-05) lives under `results/tables/cluster_level/r10/`.
