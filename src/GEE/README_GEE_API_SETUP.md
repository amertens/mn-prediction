# Earth Engine API auto-extraction — setup

Replaces the manual EE Code-Editor `.tif` export with a scripted pull via the
Earth Engine API (`rgee`). Layers are defined once in
[`gee_layer_manifest.R`](gee_layer_manifest.R) (reconstructed from
`data/GEE/GEE_export_metadata.xlsx`) and extracted by
[`extract_gee_ee_api.R`](extract_gee_ee_api.R).

## Your existing setup (already done)

You configured rgee for the mn-proxies project — see
`../mn-proxies/src/archive/GEE-download.R`. Reusable as-is:

- **EE account:** `amertens@berkeley.edu` (already EE-registered).
- **Python env:** `C:/Users/andre/AppData/Local/r-miniconda/envs/rgee/python.exe`
  (the `rgee` conda env from `rgee::ee_install(py_env = "rgee")`).

`4_TZ_GEE_extract.R` already points at both, so on your machine the only step is
re-authenticating if the token has expired.

## Cloud project (resolved)

Modern EE init requires a **Google Cloud project**. Yours is
**`mn-prediction-420517`** — hard-set as the default in
`4_TZ_GEE_extract.R` and `00_ee_connect_test.R` (override via `EE_PROJECT`).
Confirm the **Earth Engine API is enabled** for it at
<https://console.cloud.google.com/apis/library/earthengine.googleapis.com>.

## Run order

1. **Connectivity test (once):**
   ```r
   rgee::ee_Authenticate()        # only if the token is stale (opens browser)
   source(here::here("src/GEE/00_ee_connect_test.R"))
   ```
   Confirms auth + project + a real reduction (SRTM elevation and iSDAsoil zinc
   over a Tanzania box). Must print "Connectivity test PASSED".
2. **Full Tanzania extraction:**
   ```r
   source(here::here("src/Tanzania/4_TZ_GEE_extract.R"))
   ```

I never need your credentials — you run `ee_Authenticate()` locally; the scripts
just use the resulting token.

> Headless/automated alternative: a service-account JSON key
> (`ee_Initialize(user, service_account = TRUE, ...)`). Not needed for
> interactive runs.

## Run (Tanzania)

```r
source(here::here("src/Tanzania/4_TZ_GEE_extract.R"))
```
Writes `data/GEE/Tanzania_2010_admin2_gee.csv` and
`data/GEE/TZ2010_buffers_<date>.csv`, which the merge step picks up
automatically.

## Validation checklist (first run)

- Watch the per-layer `ok:` log — any layer that fails to build/extract is
  skipped with a message (e.g. the **deprecated Oxford MAP LST asset** may need
  the MODIS substitute noted in the manifest).
- Check the printed **[parity]** line: how many admin-2 column names match the
  legacy Sierra Leone file. EE band names differ from the old terra-extracted
  names, so expect partial overlap — see the naming caveat below.
- Spot-check a soil column (e.g. `gee_a2_SoilZinc`) for non-NA, sane values.

## Known caveats

- **Column-name parity.** EE band names ≠ the legacy terra-derived names, so the
  auto-extracted columns won't be byte-identical to the existing four countries.
  Fine for within-Tanzania modelling; for the **pooled / LOCO common-vocabulary
  set**, either re-extract all countries through this script (recommended for
  consistency) or add a name crosswalk.
- **2010 coverage.** Time-varying layers cover 2010 (NDVI→2013, LST→2015,
  CCNL→2013; GPW demographics is natively 2010). Single-epoch layers
  (WorldCereal 2021, WSF 2015, accessibility 2019, human-modification 2016) use
  their fixed epoch — same for every country.
- **Non-EE layers** (SPAM crop, Köppen-Geiger, Agro-Ecological Zones, FEWS NDVI
  anomaly, GLW4 livestock) are global statics already in the repo and flow in
  through the existing local-raster path — not this API script.
- **Per-dataset band/epoch details** (which band, which epoch snap) are encoded
  in the manifest from the metadata sheet but were written without a live EE
  test; the first run's logs + parity check are how we confirm them.
