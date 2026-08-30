# Admin-2 predictor covariates: download, assembly, harmonisation

Everything that produces a predictor for an admin-2 area lives here. Four
numbered stages, each with one job, run in order.

```
scripts/covariates/
  01_extract_alphaearth.R      STAGE 1  expensive downloads -> data/<source>/    (cached; run once)
  02_build_country_predictors.R STAGE 2  all sources -> one raw table per country
  03_harmonize.R                STAGE 3  raw vocabularies -> one canonical vocabulary
  04_document_and_qc.R          STAGE 4  data dictionary + cross-country QC
  run_all.R                              runs 1-4
```

```bash
Rscript scripts/covariates/run_all.R
```

## Why the stages are split this way

**Stage 1 is the only step that touches a network.** Each source writes a
cached file under `data/<source>/` and is skipped if that file exists, so
re-running the pipeline never re-downloads. Set `AEF_FORCE=1` (or delete the
file) to refresh deliberately.

**Stage 2 assembles, it does not transform.** Each country gets one table
holding every admin-2 source in that source's *native* column names. This table
is a faithful record of what each provider delivered, which means a
harmonisation decision can be revised and re-run without re-downloading.

It also fixes a real gap: the area-level model previously saw only the
Earth-Engine block. The DHS, SoilGrids, MapSPAM and ESPEN admin-2 aggregates
were merged into the *individual*-level data and never reached the primary
area-level estimator. Stage 2 puts them in one place so their use is a
modelling choice rather than an accident of plumbing.

**Stage 3 harmonises.** One canonical vocabulary, one set of units, and an
explicit record of every column kept and dropped.

**Stage 4 documents and checks.** A data dictionary you can hand to a
collaborator, plus a QC report that re-tests the thing harmonisation is
supposed to guarantee.

## What is in the harmonised dataset

**208 predictors × 725 admin-2 areas × 5 countries**, from **22 source families**
across **19 providers** in **11 conceptual domains**. Every predictor is present
in every country — that is what makes it the pooled/LOCO covariate set.

| domain | predictors | % |
|---|---:|---:|
| Land surface (learned) | 64 | 30.8 |
| Soil chemistry | 45 | 21.6 |
| Climate | 40 | 19.2 |
| Health services and demography | 18 | 8.7 |
| Vegetation | 13 | 6.2 |
| Cropping | 12 | 5.8 |
| Land cover | 6 | 2.9 |
| Built environment | 4 | 1.9 |
| Population | 4 | 1.9 |
| Access | 1 | 0.5 |
| Terrain | 1 | 0.5 |

| family | provider | n | | family | provider | n |
|---|---|---:|---|---|---|---:|
| alphaearth | Google DeepMind | 64 | | landcoverlayers | Copernicus | 6 |
| soil (iSDAsoil) | ISDA | 38 | | popdensity | CIESIN | 3 |
| dhs | DHS Program | 18 | | productivity | NASA MODIS | 2 |
| lst_night | Oxford MAP / MODIS | 17 | | accessibility | Oxford MAP | 1 |
| terraclimate | University of Idaho | 13 | | ccnl | BNU | 1 |
| mapspam | IFPRI / MapSPAM | 12 | | dailyevi | NASA MODIS | 1 |
| trmm | NASA | 10 | | elevation | USGS | 1 |
| ndvi | NOAA CDR | 8 | | ghsbuilts, ghspop | JRC | 1 + 1 |
| soilgrids | ISRIC | 7 | | globalhumanmodification | CSP | 1 |
| | | | | lai8days | NASA MODIS | 1 |
| | | | | wapor | FAO | 1 |
| | | | | wsf | DLR | 1 |

Areas per country: Ghana 260, Malawi 243, Tanzania 171, Gambia 37,
Sierra Leone 14.

Licensing: CC BY 4.0 (69 predictors), Earth Engine catalogue open (64), public
domain (23), academic-attribution (19), DHS Program terms — registration
required (18), CC0 (13), FAO open data (1), free-for-research (1). The DHS block
is the only part that cannot be redistributed without the recipient holding
their own DHS registration.

### Against the previous covariate set

| | before | after |
|---|---|---|
| Predictors shared by all five countries | 140 | **208** |
| ...of which iSDAsoil depth/dispersion columns | 99 (71%) | 38 (18%) |
| Distinct source families in the shared set | 10 | **22** |
| Conceptual domains | not tracked | **11** |
| Sources reaching the area-level model | Earth Engine only | + DHS, SoilGrids, MapSPAM, AlphaEarth |
| Variables passing the cross-country QC | not tested | 195 / 208 |
| Value-identical duplicate columns | not tested | 0 |

Recovered by renaming alone: `terraclimate` (13), `lst_night` (17),
`landcoverlayers` (6 cover fractions), `productivity`, `lai8days`, `dailyevi`,
`wapor`. Added: 64 AlphaEarth dimensions, 18 DHS admin-2 indicators, 7 SoilGrids
properties, 12 MapSPAM crop variables.

## Verification

```bash
Rscript scripts/covariates/07_verify_harmonization.R
```

18 checks, non-zero exit on any failure, so it can gate a rebuild. Current
build: **18 PASS, 0 WARN, 0 FAIL**
([results/tables/harmonization_verification.csv](../../results/tables/harmonization_verification.csv)).

The check that matters most is **10-values-reproduce**: every canonical value is
re-derived independently from the per-country raw table plus the declared unit
conversion and compared to what shipped. 1,040 of 1,040 (country, variable)
pairs match exactly. That is what catches a rule silently pointing at the wrong
source column, or a conversion applied twice.

The rest cover: no unclassified raw column; no canonical name fed by two
different source families in one country; `collapse='latest'` has a unique
winner and `collapse='none'` is genuinely 1:1; exclusions honoured; every
predictor documented with a family and a domain; one row per country × admin-2;
no all-NA, constant, or value-identical columns; NDVI within [-1, 1]; night LST
in a plausible Celsius range; cover fractions within [0, 100]; AlphaEarth
vectors unit-norm to 1.1e-16; no declared conversion leaving a >20× spread gap;
and coverage.csv agreeing with the shipped columns.

Two defects this suite caught after the first build, both now fixed:

- **37 predictors had no domain.** Rules had been added for DHS, SoilGrids,
  MapSPAM and ESPEN without matching rows in `source_registry.csv`, so those
  columns shipped undocumented.
- **11 of 14 population-density columns were exact duplicates.** GPW v4.11
  publishes 5-year epochs and the export replicated each epoch across the
  intervening years, so `popdens_y2010 == y2011 == y2012`, and so on. Only the
  published epochs (2010, 2015, 2020) are now kept.

## The problem stage 3 exists to solve

The pooled and LOCO models match covariates **by exact column name**. Measured
on the built store (2026-08):

- 1,515 distinct admin-2 covariate names exist across the five countries
- only **140** are shared by all five
- 99 of those 140 are the iSDAsoil block, leaving ~10 distinct measurement families

Yet **33 of 35 raster families are present in every country**. Almost nothing
was actually missing. Dynamic layers carry the survey year in the column name
(Gambia 2018 exported `_2017`/`_2018`; Sierra Leone 2013 exported
`_2012`/`_2013`), so the same measurement appeared under a different string in
each country and the name intersection threw it away. `terraclimate`,
`landcoverlayers`, `lst_night`, `productivity`, `lai8days`, `dailyevi` and
`wapor` were all lost this way. DHS admin-2 aggregates had the identical defect
in a different vocabulary: `dhs2019_AN_ANEM_W_ANY_adm2` vs `dhs2010_...`.

## Where the rules live

No harmonisation rule lives in code. Four CSVs under `metadata/covariates/`
drive the engine in `R/covariates/canonicalize.R`:

| file | one row per | what it decides |
|---|---|---|
| `source_registry.csv` | source family | what it is, provider, domain, causal pathway, unit, licence, status |
| `harmonization_rules.csv` | name pattern | raw name -> canonical name, keep or drop, and why |
| `unit_conversions.csv` | canonical pattern x country | multiply/add to reach the canonical unit, and why |
| `exclusions.csv` | canonical pattern | variables that cannot be reconciled and must not ship |
| `admin2_sources.csv` | source | where stage 2 finds it on disk |

A rule is `(order, regex, action, canonical template, collapse policy, reason)`.
First match by `order` wins; `canonical` may use backreferences. `collapse`
resolves the case where several raw columns in **one** country map to the same
canonical name:

- `none` — 1:1 (calendar-year series stay separate columns)
- `latest` — keep the highest embedded year (the 2-year exports)
- `mean` — average them (month climatologies spanning several years)

Anything matching no rule is reported as `unmatched` rather than silently
dropped, so a new export vintage shows up as a line in the QC summary.

## Unit defects this found and fixed

Same column name did not mean same scale. These were real, and three of them
were inside the covariate set the models were already using:

| variable | defect | fix |
|---|---|---|
| `lst_night_*` | Tanzania in **Kelvin** (~292), the other four in **Celsius** (~20) | subtract 273.15 for Tanzania |
| `npp_*` | Tanzania pre-scaled, raster exports raw DN (10⁴ apart) | multiply raster path by 1e-4 |
| `lai_*` | same defect, factor 10 | multiply raster path by 0.1 |
| `ndvi_*` | all countries raw DN | multiply by 1e-4 everywhere |
| `tclim_*` | TerraClimate scale factors never applied | per-band factor (0.1 / 0.01 / 0.001) |

Root cause: the Earth-Engine path applies the manifest `scale_factor` and the
legacy raster path never did, so any layer with a scale factor diverged between
the country extracted through Earth Engine (Tanzania) and the four extracted
from rasters.

`lst_night` additionally mixes two *products* — Oxford MAP for four countries, a
MODIS-derived source for Tanzania. After conversion the levels agree, but this
is flagged in the QC report as harmonised-across-products and should be
validated before it carries weight in a published model.

## What is deliberately excluded

Documented in `exclusions.csv` and in the `drop` rules, with reasons:

- **`soil_nitrogen_*`** — country medians of 742 / 1,090 / 2,455 / 19,051 /
  23,784 with maxima to 2.45e7. No constant factor reconciles them. This one was
  **inside the previous 140-column intersection** and was manufacturing a
  country effect in the pooled model.
- **Categorical class codes** (`landcovertype`, `ghsl_smod`,
  `discrete_classification`) — a zonal *mean* of class codes is not a quantity.
  Reinstate as per-class area fractions.
- **Non-commensurable cross-band summaries** (`fldas`, `atmosphere`,
  `gpw_demographic`, and every `*_annual_*` over a multivariate band set) — these
  average different physical quantities together; the FLDAS "climate" summary is
  dominated by surface pressure.
- **QC and assessment bands** (`*_qc`, `*_assessment`,
  `data_density_indicator`, `change_confidence`).
- **`esa_worldcereal`** — the band is masked outside cropland, so a zonal mean
  returns the class code rather than a cropland fraction, and tile ids are
  country-specific.

Several of these are recoverable with a better extraction (class fractions
instead of class-code means). They are dropped rather than fixed here because
fixing them means re-extracting, not renaming.

## AlphaEarth Foundations

`GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL` — 64 dimensions, 10 m native, annual
2017–2025. Verified against the live catalogue before wiring it in.

One fixed **2017** vintage is used for every country. Choosing each country's
nearest available year would make embedding vintage collinear with country — a
country fixed effect wearing a covariate label. Surveys run 2010–2018, so this
is a genuine temporal mismatch for Sierra Leone (−4 years) and Tanzania
(−7 years); measured drift over Sierra Leone's 14 districts between 2017 and
2024 gives between-district correlations of 0.70–0.97 by dimension. Treat AEF
as a **time-invariant land-surface descriptor**, not a covariate measured at
survey time.

Pixel vectors are unit-norm; a polygon mean is not, so the aggregate is
re-normalised to the unit hypersphere. Per-dimension z-scoring is deliberately
**not** applied — it destroys the angular geometry the embedding was trained on.
Tree learners are unaffected either way; `glmnet` standardises internally, which
is worth a sensitivity check.

## Adding a new predictor source

1. Write a stage-1 script that downloads it once to `data/<source>/` and skips
   when the file exists.
2. Add a row to `metadata/covariates/admin2_sources.csv` (where stage 2 finds it).
3. Add a row to `metadata/covariates/source_registry.csv` (what it is, unit,
   licence, causal pathway).
4. Add one rule to `metadata/covariates/harmonization_rules.csv`.
5. Re-run stages 2–4 and read `qc_summary.md`.

No code change is needed for any of this. If you skip step 4, the columns appear
in the QC summary under "matched by no rule" rather than vanishing.

## Adding a new country

Add it to the `COUNTRIES` list in stages 1 and 2 (GADM code, DHS survey name and
year), run stage 1 for that country, then stages 2–4. Read `coverage.csv`
afterwards: a new country with a thin export will shrink the shared vocabulary
for everyone, and `coverage.csv` names exactly which variables it cost.

## Using the harmonised set in the models

The area-level model reads whichever vocabulary `COVARIATE_VOCAB` selects:

```bash
COVARIATE_VOCAB=harmonized Rscript -e "targets::tar_make()"
```

Unset (or `legacy`) reproduces the published analysis of record exactly. The
switch lives in [`R/covariates/area_vocabulary.R`](../../R/covariates/area_vocabulary.R)
and is read inside `extract_gee_admin2()` and `extract_area_covariates()`, so
one setting moves every downstream consumer at once — `area_model_*`,
`run_area_comparison()`, `build_area_loco_dataset()`, the benchmark suite — and
the two vocabularies can never end up mixed in one design matrix. Under
`harmonized` the returned columns carry no `gee_` prefix, which makes a mix
detectable rather than silent.

After flipping it: `targets::tar_invalidate(matches("area_|loco|benchmark"))`.

Verify with `Rscript scripts/covariates/05_smoke_test_rewire.R [Country] [outcome]`,
which asserts that the harmonised set loads complete, that no legacy names leak
into it, and that the default path is unchanged.

The locked eight-feature soil comparator (`fit_predict_spatial_plus_soil`) now
lists both vocabularies and intersects with what is present. Seven of the eight
resolve under `harmonized`; `gee_soilaluminium_annual_min` has no canonical
equivalent because iSDAsoil `_annual_*` columns average depth means together
with their standard deviations and are dropped as non-commensurable.

## Outputs

```
data/covariates/
  country/<Country>_predictors_admin2_raw.csv   every source, native names
  country/<Country>_source_manifest.csv         what resolved, what was missing
  harmonized/predictors_admin2_harmonized.csv   pooled, canonical, analysis-ready
  harmonized/column_map.csv                     raw -> canonical, per country, with reasons
  harmonized/coverage.csv                       canonical variable x country
  harmonized/data_dictionary.csv                shareable variable documentation
  harmonized/qc_report.csv                      per-variable cross-country checks
  harmonized/qc_summary.md                      human-readable QC
```

`data_dictionary.csv` plus `predictors_admin2_harmonized.csv` are the two files
to share with a collaborator; the dictionary carries provider, licence, unit and
causal pathway for every column.
