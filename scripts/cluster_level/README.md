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
