# Area-Level Modelling & Honest Uncertainty — Fay-Herriot SAE prototype

Follow-up to the feature-engineering work (see `FINDINGS.md`). Goal: address the root
finding that **effective n = number of areas (14–87), not individuals**, by modelling at the
Admin2 level with a proper small-area-estimation (SAE) framework that (a) weights areas by
precision, (b) borrows strength from proxy covariates, and (c) yields honest uncertainty.

Prototype: `sandbox_fe/19_area_sae_prototype.R` (figure `FE_sae.png`). Outcomes/countries:
Ghana child_vitA, Malawi child_vitA, Ghana child_iron, Gambia women_iron.

## Method

Per Admin2 we form a **design-based direct estimate** `p_dir` (survey weights + cluster design)
with sampling variance `psi`, plus **bundle-reduced, supervised top-k proxy covariates**
(decorrelated, k≤5). We fit a **Fay-Herriot area model**

```
p_dir_i = x_i' b + u_i + e_i ,   e_i ~ N(0, psi_i) [known],  u_i ~ N(0, s2u)
```

with EBLUP + analytic MSE (`sae::eblupFH`/`mseFH`). The shrinkage weight
`gamma_i = s2u / (s2u + psi_i)` controls how much each area trusts its own direct estimate vs
the covariate model.

## Findings

**1. Admin2 is below the surveys' design resolution.** Every survey has only **1–3 clusters per
Admin2** (Gambia up to 13). Design-based per-area variances are therefore essentially
non-estimable (svyby SE ≈ 0 with a single PSU). **Direct Admin2 estimation is not viable; a
model-based approach (FH or multilevel) is mandatory for any area map.**

**2. There is a large optimism gap — in-sample fit overstates out-of-sample area skill.**
Correlation between area prevalence and the covariate prediction:

| cell | in-sample | leave-area-out CV | CV skill vs grand mean |
|---|---|---|---|
| Ghana child_iron | 0.70 | **0.62** | **+21%** |
| Gambia women_iron | 0.72 | 0.41 | +8% |
| Ghana child_vitA | 0.49 | 0.35 | +1% |
| Malawi child_vitA | 0.45 | **0.05** | **−13% (worse than the mean)** |

The pipeline's in-sample fit / importance numbers can overstate real skill substantially
(Malawi child_vitA: 0.45 → 0.05; the model is *worse* than just predicting the national mean).
**Report leave-area-out CV, not in-sample fit.**

**3. Proxy skill is real but outcome- and country-specific.** Ghana child_iron (+21%, CV cor
0.62) and Gambia women_iron (+8%) carry genuine area signal; Ghana child_vitA is ≈ the mean;
Malawi child_vitA has none. This matches the domain-biology result (iron has the
malaria/health signal that transfers; vitamin A is weaker and environment-driven).

**4. SAE borrows strength correctly.** Few-cluster areas are pulled toward the covariate model
(mean gamma for ≤median clusters 0.32–0.60 vs 0.67–0.85 for >median). This is the principled
way to handle the 1–3-cluster-per-area problem.

**5. Honest area uncertainty is wide; FH is the best achievable.** Direct-estimate 95% CI
half-widths are ±0.17–0.21 (uninformative for ranking areas); FH tightens them 19–33% to
±0.12–0.14 by borrowing strength. Still wide — an honest reflection of the data, and far more
defensible than an individual-level interval that implicitly credits ~1000 independent obs.

## Recommendations for the primary pipeline

- **Adopt an area-level SAE model as co-primary** for Admin2 maps (Fay-Herriot at minimum; a
  unit-level / multilevel model — e.g. `lme4`/`INLA` BYM with a spatial random effect — is the
  natural upgrade and can use the cluster-level proxies). Keep the individual-level SL as a
  sensitivity analysis.
- **Switch reported model skill and variable importance to leave-area-out CV.** In-sample
  numbers are optimistic and the gap is outcome-specific, so it cannot be applied as a flat
  correction.
- **Report area estimates with FH (or fully-Bayesian) intervals**, and state plainly that
  Admin2 is below the survey design resolution so estimates are model-dependent.
- **Where CV skill ≤ 0 (e.g. Malawi child_vitA), fall back to the direct/national estimate** —
  do not publish a proxy-model map that is worse than the mean.

## Caveats / next steps
- `psi` currently uses a binomial sampling variance (design effect ≈ 1 because 1 PSU/area makes
  the within-area DEFF non-estimable). This likely makes `psi` — and thus FH CIs — slightly
  optimistic. **Apply a design effect estimated at a level with multiple PSUs (Admin1/national)**
  before production use.
- k-covariate supervised selection inside CV is honest but small-n unstable; a fixed
  bundle + ridge (or the BYM spatial model) would be steadier.
- Spatial smoothing (neighbour structure) is not yet used; a BYM/SPDE random effect should
  recover additional area signal and is the recommended production form.

## Files
- `19_area_sae_prototype.R`, `20_sae_figure.R`
- `results_19_sae.rds` (per-area detail), `results_19_sae_summary.csv`, `FE_sae.png`
