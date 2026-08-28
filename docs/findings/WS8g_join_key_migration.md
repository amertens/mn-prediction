# WS8g. Admin-2 join-key migration

Migrate every Admin-2 join from the Admin-2 name alone to the
(Admin1, Admin2) pair, and add a regression test asserting unit counts per
country.

## A Stage 0 correction

The Stage 0 audit recorded `R/admin2_key_hygiene.R` as exporting only
`snap_water_to_land()` and being unwired. That was wrong. The module also
defines `is_water_admin2()`, `drop_water_admin2()`, `dedupe_admin2_key()`,
`clean_admin2_keys()` and `assert_unique_admin2()`, and `clean_admin2_keys()` is
called from seven places. The audit grepped for one function name and the
export list rather than for each function.

One part of the note was right: **`assert_unique_admin2()` was defined and never
called anywhere.** It exists precisely to catch the defect described below, and
it was not guarding anything.

## What was actually wrong

Measured before any change was made.

Admin-2 name collisions exist in exactly one country
(source: `results/tables/corrected/admin2_unit_counts.csv`):

| Country | GADM polygons | Duplicate land names | Duplicate (Admin1, Admin2) pairs |
|---|---|---|---|
| Gambia | 37 | 0 | 0 |
| Ghana | 260 | 0 | 0 |
| Sierra Leone | 14 | 0 | 0 |
| Malawi | 256 | 4 | 0 |

Malawi's four are genuine same-named districts in different regions: TA Lundu
(Blantyre and Chikwawa), TA Ngabu (Chikwawa and Nsanje), TA Pemba (Dedza and
Salima), TA Malemia (Nsanje and Zomba). Three of the four are surveyed, carrying
30, 38 and 45 respondents, 113 of Malawi's 3,099 rows.

The survey side is clean: all four countries have as many (Admin1, Admin2) pairs
as Admin-2 names, and every one of Malawi's 89 surveyed pairs matches a GADM
pair exactly. The collisions are entirely on the GADM-derived covariate side,
and joining a unique survey table to a non-unique GADM table on the name alone
is what produced two distinct failures.

### Failure 1: covariates averaged across regions

`dedupe_admin2_key()` collapsed same-named districts by averaging their numeric
columns. Malawi's covariate table therefore had 239 rows rather than 243, and
respondents in TA Lundu, Chikwawa were joined to covariates that were half from
TA Lundu, Blantyre.

Worse, the collapse kept whichever region came first. The stored covariate table
labels TA Lundu as **Blantyre** and TA Malemia as **Nsanje**, while the survey
observed them in **Chikwawa** and **Zomba**. The retained label is the region
the survey did not observe.

### Failure 2: silent row multiplication

`build_area_loco_dataset()` joined the pooled table to GADM centroids with
`left_join(by = "Admin2")` against the un-deduped 256-row polygon table.
Measured before the fix:

| Country | Surveyed districts | Pooled rows before | Pooled rows after |
|---|---|---|---|
| Gambia | 30 | 30 | 30 |
| Ghana | 75 | 75 | 75 |
| Sierra Leone | 14 | 14 | 14 |
| **Malawi** | **87** | **90** | **87** |
| Pooled total | | 209 | 206 |

TA Lundu, TA Malemia and TA Pemba each appeared twice, once per region's
centroid. Those three districts carried **double weight** in every area-level
fit, and one copy of each pair carried the **wrong region's coordinates**, which
corrupts the spatial comparators specifically. Every workstream on this branch
that used `build_area_loco_dataset()` inherited it: the WS2b hygiene comparison,
the WS3 uniform-BRINDA LOCO, the WS7 anchored transport, and the WS6
calibration learning curve.

## What changed

`R/admin2_key_hygiene.R`:

- `admin2_key()` builds the canonical key, Admin1 + Admin2 where Admin1 exists
  and Admin2 alone otherwise.
- `can_pair_join()` and `admin2_join_by()` decide **per join** whether both
  sides support the pair, so a partially migrated pipeline stays correct rather
  than requiring a flag day.
- `dedupe_admin2_key()` now collapses on the pair. What remains to collapse is
  genuine multi-part geometry within one region, which is what collapsing was
  for; two districts that merely share a name are no longer merged.
- `assert_unique_admin2()` is pair-aware, and is now actually called.
- `warn_if_join_multiplied()` and `report_pair_join_losses()` are new.

`R/admin2_analysis.R`: `compute_svy_admin2()` carries Admin1 through. It is
attached by lookup rather than by regrouping the `survey_mean()` call, so the
prevalences and standard errors are byte-identical and the only change is an
added column. If a survey ever does contain a within-country name collision the
lookup stops being one-to-one and it warns rather than guessing.

`R/area_level_comparison.R`: both joins in `build_area_loco_dataset()` use
`admin2_join_by()`, and the centroid table is collapsed first when the pair key
is unavailable. This is the fix for failure 2.

`R/transportability_area.R`: `assemble_area_transport()` cannot use the pair,
because `build_admin2_covariates()` does not carry Admin1. Both sides are unique
by name today so it cannot multiply, but that is a property of the inputs rather
than of the join, so it now asserts uniqueness instead of relying on it.

## A hazard the migration introduces, and its guard

A pair join is correct only when both sides agree on Admin1. They disagree while
the store is **partially rebuilt**: once `svy_admin2` gains Admin1 but
`gee_admin2` still carries the name-collapsed labels, TA Lundu (Chikwawa) no
longer matches TA Lundu (Blantyre) and the district silently disappears.

Dropping is the safe failure, since an absent district beats one wearing another
region's covariates, but it must not be silent. `report_pair_join_losses()`
detects exactly this case, distinguishing an Admin1 disagreement from a
genuinely absent district by checking whether the row would have matched on the
name alone. Verified by simulating the partial state:

```
[admin2_hygiene] Malawi survey-to-GEE join: 2 district(s) dropped by the
(Admin1, Admin2) join that WOULD match on the name alone: TA Lundu, TA Malemia.
The two sides disagree on Admin1, which usually means one of them predates the
join-key migration. Rebuild the covariate targets.
```

In practice the partial state does not arise: `tar_outdated()` confirms
`gee_admin2_*` and `svy_admin2_*` invalidate together, because both depend on
functions this change touched.

## What still has to happen

The row-multiplication fix is live now and needs no rebuild: it is entirely in
the join.

The covariate-averaging fix needs `gee_admin2_*` rebuilt. Verified against raw
GADM, the pair-aware dedupe recovers the four Malawi districts and changes
nothing elsewhere:

| Country | Raw polygons | Name-only dedupe | Pair-aware dedupe | Recovered |
|---|---|---|---|---|
| Gambia | 37 | 37 | 37 | 0 |
| Ghana | 260 | 260 | 260 | 0 |
| Sierra Leone | 14 | 14 | 14 | 0 |
| Malawi | 256 | 239 | 243 | 4 |

The store was not rebuilt here. Rebuilding `gee_admin2_malawi` re-runs zonal
statistics over 256 polygons, and the cascade it triggers is a pipeline run the
project should choose to start rather than something to slip into a migration
commit. Until then the stored covariate table keeps the averaged values, which
is the pre-existing behaviour rather than a regression.

**These numbers will change on the next full rebuild.** Malawi's pooled area
count drops from 90 to 87, three districts stop carrying double weight and gain
their own region's coordinates, and their covariates stop being cross-region
averages. That is a correctness fix, and the regression gate exists to make the
resulting diff explicit rather than silent.

## Regression test

`scripts/test_admin2_keys.R` asserts, per country: no duplicated key in the
survey table, none in the covariate table, no duplicated (Admin1, Admin2) pair
among GADM land polygons even where names repeat, and that pooling never
produces more rows than there are surveyed districts. It exits 1 on failure and
writes `results/tables/corrected/admin2_unit_counts.csv`.

It passes on all four countries.

## Reviewer passes

**Statistical reviewer.**

- *What changed in the estimates.* Nothing yet, because no target was rebuilt.
  On rebuild, three Malawi districts lose a duplicate row each. The affected
  quantity is the area-level fit's effective sample, not any individual-level
  estimate.
- *Survey estimates untouched.* Admin1 is attached to `svy_admin2` by lookup
  after the design-based computation, so no prevalence or standard error moves.
  Checked by construction: the `survey_mean()` call is unmodified.
- *Direction of the fix.* Removing duplicated rows reduces a district's weight
  from two to one, which is the correct weight. It is not a tuning choice.
- *Failure mode chosen deliberately.* Where the pair key cannot match, the
  district is dropped rather than rescued by name, because a name rescue would
  reattach the wrong region's covariates, which is the original defect.
- *No selection or leakage surface.* This workstream fits no model.

**Reproducibility reviewer.**

- `tar_manifest()` parses 845 targets, unchanged.
- `tar_outdated()` confirms `gee_admin2_*` and `svy_admin2_*` invalidate
  together, so the partial-rebuild hazard cannot arise from a normal `tar_make`.
- The regression gate reports all 52 frozen baselines unchanged, as expected:
  no table was regenerated.
- Every helper degrades to name-only behaviour when Admin1 is absent, so tables
  that never carried Admin1 are unaffected.
- Paths resolve through `here::here()`; no absolute path, no `setwd()`.
