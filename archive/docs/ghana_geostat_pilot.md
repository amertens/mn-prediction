# Pilot: cluster-level geostatistics for Ghana (vs admin-2 aggregation)

**Question:** does modeling at the survey-cluster (EA) level with a spatial field
beat the current pipeline, which aggregates to admin-2 first (forcing effective
n = #areas)? Script: `sandbox/30_ghana_geostat_pilot.R`. Outcome: Ghana child
iron deficiency (`gw_cIDAdjThurn`, ferritin < 12).

## Setup
- **1,165 children in 90 clusters** (median 12/cluster), GPS from `lon`/`lat`,
  across 75 admin-2 areas. Overall ID prevalence 23.9%.
- Each cluster = a binomial observation `(n_def, n)` at its GPS centroid.
- **INLA-SPDE** Matérn spatial field (PC priors), 809-node mesh.
- Covariates: the existing `gee_a2_*` proxies — but these are **admin-2-constant**
  (verified), so they enter as area-level fixed effects; the SPDE field supplies
  the point-level resolution. (True point-level covariates would need GEE
  re-extraction at cluster GPS — a follow-on.)
- **Honest leave-area-out CV** (10 folds over admin-2): covariates prescreened
  *inside each fold* on training clusters only (top-6 |cor|, decorrelated),
  standardized with train stats. Cluster predictions aggregated (n-weighted) to
  admin-2 and compared to the direct survey prevalence on the common 75 areas.

## Results (leave-area-out CV, aggregated to admin-2)

| Spec | MAE (pp) | RMSE (pp) | r |
|---|---|---|---|
| **A** covariate-only (6 proxies, cluster binomial) | 12.0 | 15.3 | 0.571 |
| **B** spatial-only (SPDE field, no covariates) | 12.6 | 16.2 | 0.494 |
| **C** spatial + covariates (full) | 11.5 | 14.8 | **0.605** |
| current admin-2 area-SL (honest OOF, iron) | — | — | ~0.12 |

## Findings

1. **Big jump in honest admin-2 prediction: r ≈ 0.12 → 0.60.** Moving to the
   cluster level recovers signal the admin-2 aggregation throws away. MAE ~11–12 pp
   on a 24% base.
2. **Two distinct levers, both real:**
   - *Parsimony + cluster binomial* (spec A): 6 prescreened proxies in a
     sample-size-weighted cluster model already reach r = 0.57 — far above the
     437-covariate area-SL (consistent with the project's "reduce dimensions"
     finding). Much of the gain is simply *not* pre-collapsing to 75 area rows.
   - *Spatial field* (spec B): geography **alone**, no covariates, gives r = 0.49.
     Adding it on top of covariates (C vs A) is a modest but positive lift
     (r 0.571 → 0.605; MAE 12.0 → 11.5).
3. **Honest comparison caveats:** the ~0.12 baseline is the area-SL's
   cross-outcome median honest OOF, not a Ghana-iron-specific refit; and the
   direct survey admin-2 target is itself noisy (≈1 cluster/area in Ghana), which
   *attenuates* r — so 0.60 is achieved against a noisy yardstick.

## Recommendation / next steps
- **The lever works** — worth promoting from pilot to a proper track. In rough
  order of value:
  1. **Re-extract covariates at cluster GPS** (true point-level) — unlocks the
     part of the geostatistical benefit this pilot couldn't use.
  2. **Pool countries at the cluster level** with country effects + the shared
     spatial field — this is where n multiplies (the transport goal).
  3. **Apply to the other biomarkers**; consider a joint multi-outcome field.
- Spatial field's marginal value is modest in Ghana specifically (dense, small
  country); expect larger gains where survey coverage is sparser and covariates
  are re-extracted at points.
