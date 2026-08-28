# External data sources: endpoints, credentials, refresh commands

Status as of the `analysis-updates-2026-08` branch. Every reachability claim
below was tested from this machine; every credential claim was checked against
the filesystem rather than assumed.

## Summary

| Source | Credential | Reachable | Used by | Status on this branch |
|---|---|---|---|---|
| Earth Engine | OAuth, present | yes | GEE Admin-2 covariates | working, not re-run |
| WHO VMNIS | none | see note | WS7 anchors | ingested, no machine endpoint |
| Stevens et al. 2022 | none | n/a | WS7 anchors (composite only) | ingested |
| GBD / GHDx | free account, **absent** | yes (site) | WS7 anchors, dashboard comparator | not retrieved |
| WorldPop | none | yes | population weighting | not retrieved |
| Malaria Atlas | none | yes | MAP_ covariates | cached, not refreshed |
| SoilGrids / ISRIC | none | yes | soil covariates | cached, not refreshed |
| geoBoundaries | none | yes | admin boundaries | not retrieved |
| GFDx | none | yes | fortification covariates | not retrieved |
| DHS Program | account, present | yes | hemoglobin track (flag-gated) | not retrieved, flag unset |

## Earth Engine

**Credential.** OAuth credentials at
`C:/Users/andre/.config/earthengine/credentials`, account directory
`amertens@berkeley.edu`. `rgee` 1.1.7, with the reticulate Python pinned at
`scripts/build_gee_legacy_parity.R:60`. Cloud project `mn-prediction-420517`
(`scripts/build_gee_legacy_parity.R:62`, overridable with `EE_PROJECT`).

**Verified.** `ee_Initialize(drive = FALSE)` succeeds and a live collection query
returns (`MODIS/061/MOD13Q1` size 609).

**Cost.** Extraction pulls results through `reduceRegions()$getInfo()` in chunks
(`sandbox_parsimony/python/extract_tanzania_gee.py:10`) and never calls
`Export.toDrive` or `Export.toCloudStorage`, so no Drive or Cloud Storage is
involved. Whether Earth Engine compute is billed depends on how project
`mn-prediction-420517` is registered: noncommercial and academic registrations
use Earth Engine free, subject to quota rather than charge. That registration
state cannot be read from this machine and should be confirmed in the Cloud
console before a long extraction.

**Refresh.**

```bash
Rscript scripts/refresh_gee_coverage.R
Rscript scripts/build_gee_legacy_parity.R
```

**Not run on this branch.** The planned displacement-integrated cluster-buffer
extraction was dropped on evidence; see `docs/findings/WS5_cluster_spatial_model.md`.

## WHO VMNIS

**Credential.** None.

**Endpoint.** There is no machine endpoint in use. `scripts/pull_vmnis_validation.R`
reads a local WHO long-format extract held outside this repository
(`mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta`). This is the reframe
branch: the loader validates a schema rather than fetching, and the provenance
row records that.

**Coverage limit, which is a finding rather than a caveat.** For the four
primary outcomes across the four countries, VMNIS supplies an anchor for 3 of 16
cells, has no iron entry for any country, and its three vitamin A entries
predate their surveys by 15 to 19 years. See
`docs/findings/WS7_anchored_transport.md`.

**Refresh.**

```bash
Rscript scripts/pull_vmnis_validation.R
Rscript scripts/build_anchors.R
```

## GBD / GHDx

**Credential.** A free GHDx account is required. None is present: `.Renviron`
contains no GBD, GHDx or healthdata key.

**Reachable.** `https://ghdx.healthdata.org/` returns HTTP 200.

**Status.** Not retrieved. The placeholder
`dashboard/data/gbd_estimates.rds` (`dashboard/data-raw/02_gbd_placeholder.R:72`)
is left in place. The dashboard roadmap's stated choice is to land the export or
remove the stub; removing it deletes a dashboard input, so it was flagged for
confirmation rather than done.

**To land it.** Export national prevalence with uncertainty for the matching
causes and risks from the GBD Results Tool, save to
`dashboard/data-raw/gbd_estimates.csv`, then run
`dashboard/data-raw/02_gbd_placeholder.R`, which already prefers a real CSV when
one is present (`:65`).

## WorldPop

**Credential.** None. `wpgpDownloadR` is **not installed**; use the REST API
directly.

**Reachable.** `https://hub.worldpop.org/` returns HTTP 200.

**Status.** Not retrieved. Consequence: no population-weighted aggregation
exists anywhere on this branch. WS5 and WS7 aggregate with survey n or cluster
count instead, and every affected row records that in an `aggregation` or
`aggregation_weight` column so the figures cannot be read as population
weighted.

## Malaria Atlas

**Credential.** None. `malariaAtlas` 1.6.4 installed.

**Reachable.** `https://data.malariaatlas.org/` returns HTTP 200.

**Status.** Cached under `data/external_cache/malaria_atlas`, not refreshed.
When refreshing, keep to the canonical 43-surface set (the `Malaria`,
`Interventions` and `Blood_Disorders` prefixes) rather than the full catalogue.

## SoilGrids / ISRIC

**Credential.** None.

**Reachable.** `https://rest.isric.org/soilgrids/v2.0/properties/layers` returns
HTTP 200.

**Refresh.**

```bash
Rscript scripts/build_soilgrids_admin2.R
```

## geoBoundaries

**Credential.** None.

**Reachable.** `https://www.geoboundaries.org/api/current/gbOpen/MWI/ADM2/`
returns HTTP 200.

**Status.** Not retrieved. The outstanding work it would serve is the
Admin1+Admin2 join-key migration. `R/admin2_key_hygiene.R` documents the Malawi
duplicate-territorial-authority problem and the GADM inland-water polygons, but
exports only `snap_water_to_land()` (`:116`) and is **not wired into
`_targets.R`**. That migration touches every join in the pipeline and was not
attempted on this branch.

## GFDx

**Credential.** None.

**Reachable.** `https://fortificationdata.org/` returns HTTP 200.

**Status.** Not retrieved. Intended use is national `gfdx_*` fortification
covariates, tested against the WS3 level-offset decomposition. These would be
national covariates and must never be treated as subnational variation.

## DHS Program

**Credential.** `rdhs` 0.8.1 installed, `~/.rdhs.json` present (HOME resolves to
`C:\Users\andre\OneDrive\Documents`).

**Reachable.** `https://api.dhsprogram.com/rest/dhs/surveys` returns HTTP 200.

**Status.** Not retrieved. The hemoglobin track is gated on `ENABLE_HB_TRACK`,
which is unset, so nothing was downloaded.

**Note on running `rdhs` non-interactively.** Set `R_USER` and `HOME` to
`OneDrive/Documents` so `~/.rdhs.json` loads, and set
`options(rappdir_permission = TRUE)`.

## Provenance

`metadata/external_provenance.csv` records, per artifact: source, access method,
URL, the script that retrieved it, a version or date, and licence. It currently
covers the VMNIS extract, the Stevens 2022 extracts, and the derived anchor
table. Any further acquisition should append a row there and add an md5 stamp
target in `_targets.R`, following the `gee_parity_stamp_*` and
`assay_lineage_stamp` pattern.
