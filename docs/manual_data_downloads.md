# Manual data downloads — new predictor sources

Two of the new predictor domains cannot be pulled programmatically (access changed
in 2026) and need a one-time manual download before their build script can run.
Everything else (`spam_`, `vas_`, `fao_`) is fully automated and already built.

Countries in scope: **Gambia (GMB), Ghana (GHA), Sierra Leone (SLE), Malawi (MWI),
Tanzania (TZA)**. Companion spec with full rationale/schema:
`docs/ra_new_data_sources_specs.md`.

---

## 1. ESPEN — soil-transmitted helminths / schistosomiasis (`espen_`)

**Why manual:** WHO/AFRO has **paused issuing ESPEN API keys**, so the programmatic
API is closed. The public download page needs no key.

**Steps** (5 countries × 2 diseases = 10 small CSVs):

1. Open <https://espen.afro.who.int/maps-data/data-query-tools/download-data>
2. Set the form:
   - **Country:** Gambia (then Ghana, Sierra Leone, Malawi, Tanzania)
   - **Disease:** `STH` — then repeat with `SCH`
   - **Level:** `IU`
   - **From / To:** the full available year range
3. Click **Download Data**; save each CSV into `data/ESPEN/raw/`. Filenames don't
   matter (the script auto-detects country + disease from contents), but keeping
   the disease word — e.g. `Ghana_STH.csv` — avoids ambiguity.
4. Build:
   ```bash
   Rscript scripts/build_espen_admin2.R
   ```
5. Iron-signal test (baseline vs +helminth):
   ```bash
   Rscript sandbox_fe/21_espen_iron_test.R
   ```

*If keys reopen:* request via `ntd.espen@who.int` or
<https://espen.afro.who.int/espen-api-platform>, set `ESPEN_API_KEY` in
`~/.Renviron`, install `ESPENAPI` — the script then pulls automatically.

---

## 2. FPN — Food Prices for Nutrition affordability (`fpn_`)

**Why manual:** the World Bank retired the `source=88` **data** API (migrated to
Data360, which errors), and the old `CoNA_*` codes no longer exist (now `CoHD_*` /
`CoCA_*`).

**Steps** (one CSV covers all countries):

1. Open <https://databank.worldbank.org/source/food-prices-for-nutrition>
2. **Country:** Gambia, Ghana, Sierra Leone, Malawi, Tanzania
3. **Series** (current FPN 4.0 codes):
   - `CoHD_headcount`, `CoHD_PPP`, `CoHD_fexp` — healthy-diet cost/affordability
   - `CoCA_headcount` — energy-sufficient diet
   - `CoHD_lns_prop` (legumes/nuts/seeds cost share), `CoHD_asf_prop` (animal-source cost share)
4. **Time:** all available years
5. **Download → CSV** (default wide layout); save into `data/FPN/raw/`
6. Build:
   ```bash
   Rscript scripts/build_fpn_affordability.R
   ```

---

## Already automated — no action needed

| Domain | Prefix | Builder | Status |
|---|---|---|---|
| MapSPAM crop mix | `spam_` | `scripts/build_mapspam_admin2.R` | built, all 5 countries |
| Vitamin A supplementation | `vas_` | `scripts/build_vas_national.R` | built, all 5 countries |
| FAOSTAT food supply | `fao_` | `scripts/build_faostat_supply.R` | built, all 5 countries |

After any build script runs, the `file.exists()`-guarded joins in each
`src/<Country>/2_GW_<Country>_data_merge.R` pick the data up automatically on the
next merge run — no further wiring needed.

---

## Optional — the two spec-only sources (not built or wired yet)

- **GFDx** (fortification, `gfdx_`): download the dataset CSV from
  <https://fortificationdata.org/full-gfdx-datasets/> (no key).
- **GDD** (diet intake, `gdd_`): install the `nutriR` package for the subnational
  inadequacy distributions; the full estimate set is a gated bulk request at
  <https://globaldietarydatabase.org/>.
