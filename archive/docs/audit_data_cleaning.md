# Data-cleaning / harmonization audit

READ-ONLY audit of the data-cleaning and harmonization code (2026-06-23).
Scope: `R/data_prep.R`, `R/config.R`, `R/external_data.R`, `R/merge_external.R`,
`R/ingest_new_country.R`, `R/cluster_aggregation.R`, `R/feature_engineering*.R`,
`R/food_security.R`, `R/admin2_analysis.R` (uniform-cutoff path), and the
per-country cleaning in `src/{Gambia,Ghana,SierraLeone,malawi,IHME,dm}`.

Findings that touch the merged caches were verified against the on-disk RDS
files (`readRDS` only — the `_targets` store was not touched). Verified items
are marked **[verified]**; unverified inferences are marked **[flag]**.

---

## HIGH severity

### H1. Outcome leakage: `gw_exclude_patterns` is incomplete; iron/VitA/Hb outcome-derived columns leak into predictors **[verified, Sierra Leone; flag others]**
- `R/config.R:343-347` (SL), `:103-107` (Gambia), `:231-235` (Ghana), `:479-483` (Malawi); applied in `R/data_prep.R:445-453`.
- The GW leakage filter matches `"Ferr"` (two r's) but the SL iron outcome
  column is `gw_wFerAdjBR1` / `gw_LNwFerAdjBR1` — spelled `"FerAdj"` (one r) —
  so it does NOT match. The exact-column exclusion at `data_prep.R:451` removes
  the *active* outcome column (`gw_wFerAdjBR1`) per-outcome, but its log sibling
  **`gw_LNwFerAdjBR1` is never the named outcome and therefore always leaks** as
  a near-perfect linear proxy of the women-iron outcome.
- Also surviving the SL filter (verified present in the merged data as kept
  predictors): `gw_wVitADef`, `gw_cVitADef`, `gw_wVitADefAdjBR1`,
  `gw_wVitAInsuff` (pattern is `"VAD"`, not `"VitA"`), and `gw_wHb`, `gw_cHb`,
  `gw_HbCat` (hemoglobin, which directly determines anemia/IDA). When predicting
  one population, the OTHER population's deficiency column also leaks
  (e.g. `gw_cVitADef` present while modeling women_vitA).
- Why it's a bug: a column that is (a deterministic function of) the outcome is
  available as a predictor, inflating in-sample and CV performance and
  invalidating transportability claims.
- Severity: HIGH.
- Fix: broaden the exclude patterns to cover the actual spellings
  (`"Fer"`, `"VitA"`, `"Hb"`, `"Iod"`, `"Goiter"`, `"anemia|anaemia"`,
  `"_LN"`-prefixed log biomarkers), or — more robustly — switch to a
  positive allowlist of permitted gw_ predictor stems instead of a denylist.

### H2. Cross-country inflammation-adjustment inconsistency in the VitA uniform-transport outcome **[verified, Malawi]**
- `R/admin2_analysis.R:20` (`UNIFORM_TRANSPORT_TAGS` includes `child_vitA`,
  `women_vitA`) + `:317-327` (overwrites the binary with
  `apply_threshold(oc$continuous, 0.70)`).
- The continuous RBP column differs by country: Gambia/Ghana use
  Thurnham-adjusted RBP (`gw_*RBPAdjThurn`, `config.R:42,57,142,153`), SL uses
  `gw_*RBPAdj` (`config.R:269,280`), **Malawi uses RAW `rbp`** (`config.R:380,390`).
- Verified in the Malawi cache: the reported `vitA_def` for preschool children
  (0.199) ≠ raw `rbp<0.70` (0.239) → children's binary is BRINDA-adjusted; but
  for women and school-age children reported `vitA_def` == raw `rbp<0.70`
  exactly → no adjustment. Malawi also HAS `crp`/`agp`, so adjustment is
  feasible but was only applied to children (see H3).
- Net effect: the uniform-transport override applies the SAME cutoff (0.70) to
  DIFFERENTLY adjusted continuous columns across countries (Thurnham vs raw),
  so the "harmonized" VitA outcome is not actually harmonized — it reintroduces
  the cross-country level offset it is meant to remove. Matches the
  `fe_transport_level_offset` memo class.
- Severity: HIGH (undermines LOCO/transport comparability for vitamin A).
- Fix: compute a Thurnham/BRINDA-adjusted continuous RBP for ALL countries
  (Malawi has crp/agp; use `harmonize_rbp_thurnham()` in
  `R/ingest_new_country.R`) and point every country's VitA `continuous` at the
  adjusted column before applying the 0.70 cutoff.

### H3. Malawi VAD binary applies inflammation adjustment to children only, raw to women/school-age **[verified]**
- `src/malawi/1_DHS_mn_data.R` (define_mn_deficiency): children VAD uses
  `brinda_adjust(rbp)`, women/school-age use raw `rbp < 0.70`.
- Verified against the cache (see H2): children adjusted, women/school-age raw.
- Why it's a bug: within one survey the VAD definition is inconsistent across
  populations; raw RBP overstates deficiency under inflammation.
- Severity: HIGH (within-country); compounds H2.
- Fix: apply the same adjustment to all populations and save an adjusted
  continuous RBP column.

### H4. Iron continuous-outcome scale/adjustment differs by country (log-Brinda vs linear-Thurnham vs mixed vs regression-adjusted) **[verified]**
- Gambia: `gw_LogFerAdj` (LOG scale, Brinda), cutoffs `log(12)`/`log(15)`,
  `cutoff_scale="log"` (`config.R:68-83`).
- Ghana: `gw_cFerrAdjThurn`/`gw_wFerrAdjThurn` (LINEAR µg/L, Thurnham), cutoffs
  12/15 (`config.R:170-188`).
- Sierra Leone: `gw_cFerrAdj` (children) vs `gw_wFerAdjBR1` (women) — DIFFERENT
  variables/adjustments for the two populations (`config.R:291,302`). Verified:
  `gw_wFerrAdj` has only n=389 (mostly NA) while `gw_wFerAdjBR1` has n=774, so
  the workaround is necessary but means women use a different adjustment than
  children.
- Malawi: `sf_reg` (verified linear µg/L, median 34.6 — a regression-adjusted
  ferritin, distinct from raw `fer` median 46), cutoffs 12/15.
- `apply_threshold()` (`R/admin1_analysis.R:13`) ignores `cutoff_scale`; correct
  output depends entirely on the config keeping each cutoff on the same scale as
  its continuous column. This coupling is correct as currently configured
  (Gambia log cutoff + log column; others linear+linear) but is fragile: any
  future edit that points a log column at a linear cutoff (or vice versa) will
  silently produce wrong prevalences with no error.
- Severity: HIGH for cross-country pooling/transport; MED within-country.
- Fix: standardize one adjustment + scale for the continuous iron outcome across
  all countries; have `apply_threshold` (or the caller) assert/honor
  `cutoff_scale` rather than relying on hand-matched config values.

---

## MEDIUM severity

### M1. Ghana B12/folate outlier capping lives only in a sandbox script, not the production clean path **[flag]**
- `sandbox/20_upstream_data_fixes.R` documents Ghana `gw_B12_pmol_L`
  mean ≈ 12,139, sd ≈ 50,012 (implausible — likely a pg/mL vs pmol/L unit error
  in some rows) and caps B12 at 1500 and folate at 50 nmol/L. These caps are not
  in `src/Ghana/1_GW_Ghana_data_clean.R` or `R/data_prep.R`.
- If the sandbox cleaning is not wired into `Ghana_merged_dataset.rds`, the
  production Ghana B12/folate continuous outcomes contain order-of-magnitude
  outliers. (`R/ingest_new_country.R` provides `cap_b12_outliers()` /
  `cap_folate_outliers()` but they are only called in the new-country ingest
  path, not for Ghana.)
- Severity: MED-HIGH. FLAG: the Ghana merged RDS is not present at the configured
  path in this checkout, so this could not be verified on disk.
- Fix: move the caps into the Ghana clean script (or call the `cap_*` helpers in
  `load_merged_data`); investigate whether high B12 is a unit error (×0.738)
  rather than true outliers.

### M2. `prune_predictor_cols` IHME count-drop regex can over/under-match **[verified pattern]**
- `R/data_prep.R:344`: `grep("^ihme_.*(_counts?$|number_of_)", ..., ignore.case=TRUE)`.
- `number_of_` is matched anywhere (not anchored), and `_counts?$` only matches a
  trailing token. Columns like `ihme_..._count_per_1000` (rate-like but
  containing "count") would survive, while a genuinely count-encoded column not
  ending in `count(s)` and lacking `number_of_` would survive too. Low blast
  radius (only affects which predictors are dropped, not outcomes), but worth a
  spot-check of the actual ihme_ column inventory.
- Severity: MED. Fix: validate the regex against the real ihme_ column list;
  prefer an explicit list of count columns.

### M3. SL 1/2→0/1 recode is silently skipped if a column isn't EXACTLY coded {1,2} **[verified columns are {1,2}]**
- `R/data_prep.R:156-170` (`recode_12`) and the duplicate at `:397-404`
  (`build_outcome_dataset`) only fire when `all(vals %in% c(1,2))`.
- Verified: SL `gw_wFolDef`, `gw_wB12DefWHO`, `gw_cIDA`, `gw_wIDA` are all coded
  {1,2} in the cache, so the recode currently fires correctly for folate/B12,
  and iron is additionally protected by the UNIFORM_TRANSPORT_TAGS override
  (which rebuilds the binary from the continuous column). So this is presently
  correct — but it is brittle: if any future column contains a stray 0/9/NA-code
  the guard fails silently and a {1,2} (i.e. "2 = not deficient") column flows
  downstream as if 2 were a positive value.
- Severity: MED (latent). Fix: assert the expected value set and `stop()`/`warn`
  on mismatch instead of silently passing through; consolidate the two redundant
  recode code paths into one.

### M4. Gambia GADM lookup uses ISO2 "GM" instead of ISO3 **[flag]**
- `src/Gambia/2_GW_Gambia_data_merge.R` calls `load_gadm_cached("GM", ...)`;
  config `gadm_code = "GMB"` and SL uses `"SLE"`.
- `load_gadm_cached` (`R/data_prep.R:29`) passes the code to
  `geodata::gadm(country=...)`, which expects ISO3 or a full country name; ISO2
  "GM" is not accepted and the download would fail. Likely masked by an existing
  cache file.
- Severity: MED (reproducibility — breaks on a clean machine). Fix: use `"GMB"`.

### M5. Admin-2 fuzzy name matching always assigns a best match before thresholding **[verified pattern]**
- `R/merge_external.R:83-91` (agrep, `max.distance=0.2`, picks closest by nchar
  on multi-hit), the DHS block at `:187-189` (`agrep` then `fuzzy[1]` — takes the
  FIRST hit, not the closest), `R/food_security.R:403-442` (adist, normalized
  <0.3), and `R/external_data.R:1103-1110` (agrep then `fuzzy[1]`).
- Two concerns: (a) the DHS and external-data fuzzy paths take `fuzzy[1]`
  (arbitrary order) rather than the closest match; (b) survey districts with no
  true counterpart can still receive a spurious match, attaching the wrong
  district's covariates. Match rates are logged but individual mappings are not
  surfaced for review.
- Severity: MED. Fix: pick the minimum-distance candidate consistently across all
  fuzzy paths; log the full survey→source name map for hand-verification; raise
  the threshold or require a confirmed mapping table for the 4 known countries.

### M6. `merge_food_security` doc/behavior mismatch on phase 3-5 and Admin2 partial-merge **[verified]**
- `R/food_security.R:12` (module doc) says "phase 3+5" but the code computes
  phase 3+4+5 (`:183-197`) — code is correct (crisis+emergency+famine), the doc
  comment is wrong/misleading.
- `merge_by_admin` (`:334-356`): when Admin2 `match_rate >= 0.3`, it left-joins
  only the matched source rows and RETURNS immediately — it never falls back to
  Admin1 for the unmatched 70%, so up to ~2/3 of districts can silently get NA
  fsec_ values while the log reports a "successful" Admin2 merge.
- Severity: MED. Fix: fix the doc; for Admin2 matches below ~0.7, fill unmatched
  districts via the Admin1 fallback rather than returning early.

---

## LOW severity / fragility

### L1. Malawi zinc time-of-day fallback applies the fasting (higher) cutoff to unknown-time samples **[flag]**
- `src/malawi/1_DHS_mn_data.R`: zinc cutoff switches on `time` (morning/fasting
  65/66 vs afternoon 57/59); the `else`/NA branch uses the fasting (higher)
  cutoff, which the comment calls "conservative" but for a less-than rule a
  higher cutoff INCREASES measured deficiency. Verify intended direction.
- Severity: LOW. Fix: confirm the time coding and whether the higher-cutoff
  fallback is intended.

### L2. Aggressive silent row drops in cleaning **[flag]**
- `src/Gambia/1_GW_Gambia_data_clean.R` drops rows with missing lat/long and
  missing `cnum` (with a "temp drop, debug later" comment). Clusters lacking GPS
  are permanently removed; verify how many biomarker observations are lost.
- Severity: LOW-MED.

### L3. Long positional `subset(select = -c(...))` drop lists are fragile **[verified]**
- Gambia/Ghana/SL clean scripts hard-code large column-drop lists; these also
  drop the `_logcrpcoeff*` adjustment-coefficient columns, which is why the SL
  RBP/ferritin adjustment METHOD can no longer be determined from the data
  (relevant to H2/H4). If upstream `.dta` columns change, the wrong columns are
  dropped silently. Severity: LOW (fragility).

### L4. Ghana log-ferritin columns derived in clean script but absent from cache **[flag]**
- `config.R:170-173` notes `gw_cLogFerrAdjThurn`/`gw_wLogFerrAdjThurn` are
  created by the clean script but "absent from the current merged cache," so the
  config falls back to the LINEAR `gw_*FerrAdjThurn` columns. This means the
  Ghana cache is stale relative to the clean script, and Ghana iron is modeled on
  a linear scale while Gambia uses log — interacts with H4. Severity: LOW-MED.
  Fix: re-run the Ghana merge so the cache matches the clean script, or
  explicitly standardize the scale.

---

## Items checked and found CORRECT (no bug)

- **Malawi folate** (`R/data_prep.R:81-88`): the previously-reported
  `× 2.266` + `<3` double-bug is fixed — `fol_nmol <- fol` (no conversion) and a
  `< 10 nmol/L` WHO cutoff, matching Ghana/SL. Verified `fol`/`fol_nmol` handling.
- **Malawi B12** (`:94-100`): `< 148 pmol/L`, correct units/direction.
- **Malawi `sf_reg`** EXISTS in the cache (verified; linear ferritin, median 34.6)
  — an earlier inference that it was missing was a FALSE POSITIVE.
- **SL config columns** (`gw_cFerrAdj`, `gw_wFerAdjBR1`, `gw_wFolate`,
  `gw_wFolDef`, `gw_B12`, `gw_wB12DefWHO`, RBP/IDA) all VERIFIED present in cache.
- **Cutoff directions** are uniformly `"less"` for all deficiencies — correct
  (deficiency = below threshold) for ferritin/RBP/folate/B12/zinc.
- **`merge_external.R` dedup** (`:111-133`) correctly collapses duplicate-name
  GADM polygons (e.g. "Lake Malawi") by averaging before the join, and guards
  against row-count inflation (`:140-143`).
- **External `NA_real_` assignments** in `R/external_data.R` (lines 773, 1135-
  1156, 1648-1672, etc.) are legitimate "download/extract failed" fallbacks, not
  the VMNIS-class overwrite of a usable column.
