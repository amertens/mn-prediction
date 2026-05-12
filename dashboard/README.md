# Micronutrient Burden Dashboard

Interactive Shiny dashboard for exploring sub-national micronutrient
deficiency predictions across The Gambia, Ghana, Sierra Leone, and Malawi.

## Tabs

1. **Map explorer** — Interactive choropleth with layer toggles (predicted prevalence, observed prevalence, confidence interval width, population at risk, WHO classification). Click any district for a side-panel detail card.
2. **District profiles** — Pick a country and district to see all available micronutrient outcomes side by side, with country-average comparison line and 95% conformal CIs.
3. **National burden** — Population-weighted national prevalence with toggleable population year (survey year vs 2023 projection), "hidden burden" indicator, and pipeline-vs-survey comparison.
4. **Scenarios** — Two modes:
    - **Coverage scenario:** simulate an intervention applied to all districts, only above-average districts, or top-N highest prevalence; specify coverage and effect size; see before/after maps and cases averted.
    - **Sensitivity scenario:** illustrative tool for projecting how prevalence might shift under climate, food price, food security, or conflict perturbations.
5. **GBD comparison** — Compare pipeline national estimates against IHME Global Burden of Disease estimates. Currently using placeholder values; an RA task is documented inline for replacing with actual GBD Results Tool exports.
6. **Methods** — Plain-language description of modelling, validation, conformal CIs, limitations, and data sources.

Every table includes a 1–3 paragraph methods note below it.

## Structure

```
dashboard/
├── app.R                              # Entry point — sources global.R, builds UI/server
├── global.R                           # Shared data loading and config
├── deploy.R                           # shinyapps.io deployment script
├── README.md
├── R/
│   ├── fct_helpers.R                  # Reusable helpers (formatting, joining, aggregation)
│   ├── mod_map_explorer.R             # Tab 1
│   ├── mod_district_profile.R         # Tab 2
│   ├── mod_national_burden.R          # Tab 3
│   ├── mod_scenarios.R                # Tab 4 (coverage + sensitivity)
│   ├── mod_gbd_compare.R              # Tab 5
│   └── mod_methods.R                  # Tab 6
├── data/                              # Pre-built data (built by data-raw scripts)
│   ├── admin2_predictions.rds
│   ├── admin2_population.rds          (with 2023 projection)
│   ├── admin2_boundaries.rds
│   ├── admin1_boundaries.rds
│   ├── national_estimates.rds
│   ├── cv_performance.rds
│   ├── gbd_estimates.rds              (placeholder — see RA task in module)
│   └── metadata.rds
├── data-raw/                          # Build scripts and tests (not deployed)
│   ├── 01_prepare_dashboard_data.R    # Main data builder
│   ├── 02_gbd_placeholder.R           # GBD framework data
│   ├── smoke_test.R
│   ├── test_app_construction.R
│   ├── test_endpoints.R
│   ├── test_endpoints_v2.R
│   ├── test_server.R
│   └── test_deploy_ready.R
└── www/                               # Optional static assets
```

## Running locally

```r
# 1. Build dashboard data (only needed when pipeline outputs change)
Rscript dashboard/data-raw/01_prepare_dashboard_data.R
Rscript dashboard/data-raw/02_gbd_placeholder.R

# 2. Launch the app
setwd("dashboard")
shiny::runApp()
```

## Deploying to shinyapps.io

```bash
Rscript dashboard/deploy.R
```

Prerequisites and one-time setup are documented at the top of `deploy.R`.
On the current machine the `amertens` shinyapps.io account is already
configured.

## Required packages

- `shiny`, `bslib`, `bsicons`
- `dplyr`, `tidyr`, `sf`
- `leaflet`, `plotly`, `reactable`, `htmltools`
- `rsconnect` (for deployment only)

## Data refresh

Whenever the upstream pipeline (`_targets_full`) is rerun, regenerate the
dashboard data:

```bash
Rscript dashboard/data-raw/01_prepare_dashboard_data.R
Rscript dashboard/data-raw/02_gbd_placeholder.R
```

The dashboard never reads from `_targets_full` directly — only from the
curated RDS files in `dashboard/data/`. This keeps the dashboard fast and
makes it easy to deploy to a server that doesn't have the full pipeline.

## Outstanding RA tasks

- **GBD Results Tool data** — The GBD comparison tab uses placeholder values.
  An RA needs to download actual GBD prevalence estimates for the relevant
  countries, outcomes, and years from
  https://vizhub.healthdata.org/gbd-results/ and save as
  `dashboard/data-raw/gbd_estimates.csv`. Detailed instructions are in the
  GBD Comparison tab itself and at the top of
  `dashboard/data-raw/02_gbd_placeholder.R`.

- **Population projection refinement** — Current 2023 projections use
  uniform country-level annual growth rates from World Bank Population
  Estimates. For more accurate sub-national projections, this could be
  replaced with WorldPop 2023 raster data extracted to Admin-2 polygons.
