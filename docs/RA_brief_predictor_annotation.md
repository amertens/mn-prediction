# RA brief — predictor annotation and taxonomy

**Status:** ready to start. **Estimated effort:** 3–5 days.
**Deliverable:** a completed `results/tables/protocol_v2/variable_sheet.csv`.

## Why this matters (read this first)

The project's best-transporting model is a **domain index**: each of 18
conceptual domains is collapsed to a single score, the scores are weighted by
their association with the outcome, and the weighted sum ranks districts. It has
no tuned hyperparameters and it beats an elastic net over all 373 raw predictors
by up to +0.357 rank correlation (16 of 18 cells).

That means **the domain taxonomy is doing more work than any estimator in the
project** — and the taxonomy is a hand-curated column that was never built for
this purpose. Two measurements show where it is weak:

| Domain | members | PC1 variance explained | PCs needed for 80% |
|:---|---:|---:|---:|
| Agricultural production, land use | 93 | **0.17** | 17 |
| Infant and child morbidity/mortality | 37 | **0.13** | 17 |
| Soil characteristics | 45 | 0.24 | 11 |
| Education, employment, SES | 9 | 0.56 | 3 |
| Ruralness, population density | 9 | **0.84** | 1 |

The big domains have no dominant factor, so collapsing them to one number
discards most of their information. Ruralness genuinely is one factor and loses
nothing. Better sub-domains would let the large domains be represented by
several coherent scores instead of one incoherent average — and a follow-up
experiment confirms the gain is real: using more components per domain raises
leave-one-country-out transport from 0.151 to 0.255 (level target).

## What exists already

`results/tables/protocol_v2/variable_sheet.csv` — 373 variables × 38 fields,
built by `scripts/protocol_v2/06_build_variable_sheet.R`. Every field is tagged
with its provenance so you always know what is documented and what is a guess:

- `definition_source` — `data_dictionary` (294 rows, trustworthy) or
  `MISSING_needs_RA` (**79 rows**).
- `subdomain_source` — `AI_PROPOSED_unverified` for **every** row. The proposed
  sub-domain is a naming regularity (strip depth bands, months, years and
  statistic suffixes, keep the first two tokens), not a semantic judgement.
- `unit_source` — same split as definitions: 294 documented, 79 missing.
- Empirical columns computed from the data: `var_type`, `n_levels`, `median`,
  `IQR`, `min`, `max`, `pct_missing_overall`, `worst_country_pct_missing`,
  `countries_present`, `country_specific`, `constant_in_some_country`.
- Six blank `ra_*` columns for you to fill.

`results/tables/protocol_v2/variable_sheet_gaps.csv` is the 79-row subset with
no documented definition or unit — start there.

## Tasks

### 1. Fill the 79 undocumented variables (highest priority)

These are the "extra domain" blocks added after the original harmonisation:
IHME, Malaria Atlas, DHS prior-round aggregates, MapSPAM, ESPEN, food security.
None has an entry in `data/covariates/harmonized/data_dictionary.csv`.

For each, from the **provider's own documentation** (not from the variable
name), record in the `ra_*` columns: what the variable measures, its unit, its
measurement year or window, its native spatial resolution, and how it was
aggregated to Admin-2. Cite the documentation URL or PDF page in `ra_notes`.
Where the provider's definition contradicts the variable name, say so
explicitly — that is a finding, not a nuisance.

### 2. Verify or replace the AI-proposed sub-domains

Fill `ra_verified_subdomain` for all 373. The proposed values are a starting
point, and their quality varies enormously by domain: Malaria has 23 proposed
sub-domains for 23 variables (useless — every variable is its own group), while
Household assets, Dietary diversity, Education and Healthcare access each got
exactly **one** proposed sub-domain for 7–12 variables (also useless — no
structure at all).

Aim for sub-domains that are **internally coherent and separately
interpretable**, with roughly 3–15 members each. The 93-variable agriculture
domain is the most valuable target: splitting it into, say, crop mix, livestock,
land cover, production intensity and market orientation would let each be scored
separately instead of averaged into noise.

### 3. Add multi-label mechanism tags

Fill `ra_mechanism_tag` with **all** that apply, semicolon-separated, from:
`intake`, `absorption`, `requirement`, `loss`, `status_measurement`,
`socioeconomic`, `environmental`, `confounder_only`.

Single-domain assignment is currently forced and it misclassifies genuinely
cross-cutting variables — `dhs_hh_cattle` is agriculture *and* household assets
*and* a plausible dietary (animal-source food) proxy, and it is one of the
replicated survivors, so getting it right matters.

### 4. Add a distal–proximal ordinal

Fill `ra_distal_proximal` with `1_proximal` (directly about what people eat or
absorb), `2_intermediate` (food availability, access, care), or `3_distal`
(climate, soil, geography). This lets the project test whether proximal
predictors outperform distal ones — a hypothesis it has asserted but never
tested with a clean variable partition.

### 5. Record measurement year

Fill `ra_measurement_year` with the actual year or window the variable
measures, **not** the vintage of the file. This is the input to the standing
temporal-alignment work (`docs/RA_tasks_temporal_alignment.md`). Note that
surveys span 2013–2021, so a predictor's lag relative to each survey differs by
country, and nobody has documented what lag each variable carries.

## Acceptance criteria

- All 373 rows have `ra_verified_subdomain`, `ra_mechanism_tag` and
  `ra_distal_proximal` filled.
- All 79 gap rows have `ra_verified_definition` and a documentation citation.
- No sub-domain has a single member unless that variable is genuinely unique.
- No sub-domain in a domain of 20+ variables has more than about 20 members.
- A short note listing every case where documentation contradicted the variable
  name, or where you could not find documentation at all.

## What NOT to do

Do not rename or delete variables, and do not edit
`predictors_admin2_shared.csv` or any file under `data/`. This task produces
**annotation only**. The modelling code reads the taxonomy from the sheet; the
data itself must not move underneath it.
