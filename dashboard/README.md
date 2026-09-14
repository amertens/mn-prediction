# Micronutrient Burden Dashboard

Shiny dashboard of district rankings of micronutrient deficiency for The
Gambia, Ghana, Sierra Leone and Malawi, built from public data and scored
against the four national biomarker surveys under the corrected protocol
(protocol v2), with a ranking for Cote d'Ivoire, which has no survey.

Rebuilt on 2026-09-13. Everything the app shows comes from the committed
protocol result tables plus one fit (the deployment ranking); the earlier
layers (person-level SuperLearner, area-level recipe, Fay-Herriot, BYM2, the
old leaderboard, the P1 to P8 comparison, the GBD placeholder, the Sierra
Leone chiefdom layer) were removed because each rested on an evaluation the
audit withdrew or on a model the protocol found no better than chance.

Live app: <https://amertens.shinyapps.io/micronutrient-burden/>

## Tabs

| Menu | Tab | What it shows |
|---|---|---|
| | Start here | What the dashboard is, how to read the map, a worked example, the checklist of what the models can and cannot do yet |
| Where is deficiency? | Map explorer | Priority score per district, chance of the worst fifth, planning prevalence anchored to the national survey, the survey's own estimate; click a district for its numbers and drivers |
| | District profiles | One district across every outcome, with the exact predictor decomposition of its score |
| | Cote d'Ivoire | 33 districts ranked from climate and soil alone, with rank uncertainty |
| What drives it? | What drives the estimate | Back-projected index weights per outcome, which data groups carry the model and which travel, the twenty-layer composite, the recurring gradient |
| | What tracks which nutrient | Cross-country sign replication of district associations (signal probes) |
| | Predictor catalogue | Every predictor with definition, source, coverage, weight per outcome, replication and a small map |
| Can we trust it? | How well it works | The three tests against matched comparators and the null, the reliability ceiling, the learning curve, the geostatistical comparator, what else was tried |
| | What the ranking buys | Burden reached by the worst fifth, calibration of the worst-fifth probability, WHO band accuracy |
| | Methods | The model, the protocol, performance, surveys, data, what changed, limits |
| | Plan a survey | The anchor-and-rank survey design against district and regional surveys |

## Structure

```
dashboard/
├── app.R                      # entry point
├── global.R                   # data loading, headline numbers (Q), caveats, about, glossary
├── deploy.R                   # shinyapps.io deployment
├── R/
│   ├── fct_helpers.R          # joins, decomposition, plot helpers
│   ├── mod_start_here.R
│   ├── mod_map_explorer.R
│   ├── mod_district.R
│   ├── mod_civ.R
│   ├── mod_importance.R
│   ├── mod_nutrient_signal.R
│   ├── mod_catalogue.R
│   ├── mod_trust.R
│   ├── mod_targeting.R
│   ├── mod_methods.R
│   └── mod_survey_design.R
├── data/                      # built bundles (gitignored)
│   ├── admin2_index.rds       # the deployment ranking, survey estimates, anchors, per-column weights
│   ├── civ_index.rds
│   ├── protocol_evidence.rds  # benchmark, targeting, ceiling, design, importance, comparator tables
│   ├── predictor_catalogue.rds
│   ├── nutrient_signal.rds
│   ├── admin2_population.rds, admin2_boundaries.rds, admin1_boundaries.rds, metadata.rds
│   └── oos_cote_divoire.rds   # kept for its Cote d'Ivoire boundaries
├── data-raw/
│   ├── 01_prepare_dashboard_data.R   # population, boundaries, metadata (targets-based; run rarely)
│   ├── 03_build_nutrient_signal.R
│   ├── 05_build_protocol_v2_bundles.R  # everything else
│   ├── smoke_test.R
│   └── test_server.R
└── report/                    # printable country briefs and the overview
```

## Building the data and running

```r
# from the repo root
Rscript dashboard/data-raw/05_build_protocol_v2_bundles.R   # a minute or two
Rscript dashboard/data-raw/03_build_nutrient_signal.R
Rscript dashboard/data-raw/smoke_test.R
Rscript dashboard/data-raw/test_server.R
# then
setwd("dashboard"); shiny::runApp()
```

Rebuild the bundles whenever the protocol tables under
`results/tables/protocol_v2/` or `results/tables/policy_deck/` change. The
builder fits the index once per country and outcome (24 cells) and reads
everything else from those tables.

## Deploying

```bash
Rscript dashboard/deploy.R
```

## Briefs

```bash
Rscript dashboard/report/render_reports.R            # all countries + overview
Rscript dashboard/report/render_reports.R ghana      # one
```

The briefs source `global.R`, so they cannot disagree with the screen.

## Known gaps

- The variable annotation sheet defines 308 of the 454 predictors; the
  catalogue shows the rest as "definition pending" and is the worklist.
- The chance of being in the worst fifth exists for surveyed districts in
  three countries; Sierra Leone's 14 districts cannot be cross-validated.
- Nothing in the app is linkable (no URL state), and the boundary files are
  most of the bundle size. Both are on the roadmap.
