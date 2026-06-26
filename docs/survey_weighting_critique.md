# Critical evaluation: survey weights & subnational (admin-2) prevalence

Scope: how `compute_svy_admin2()` / `compute_national_estimates()` use the
survey design (weights, PSU, strata) to produce admin-2 prevalences, whether
that is sound given the surveys' design, and whether the **right weights** are
used. Evidence below is from the merged country data + `R/config.R`.

---

## 1. Admin-2 is below the surveys' design resolution (the core problem)

Clusters (PSUs) per admin-2 area:

| Country | admin-2 areas | median clusters/area | % areas ≤2 clusters | n/area (median) |
|---|---|---|---|---|
| Gambia | 30 | **1** | 77% | 64 |
| Ghana | 75 | **1** | 97% | 29 |
| Sierra Leone | 14 | 4 | 14% | 96 |
| Malawi | 89 | **1** | 97% | 32 |

For **Gambia, Ghana, and Malawi the typical admin-2 area contains a single
sampled cluster.** Consequences:

- **Design variance is not estimable from one PSU.** The code sets
  `survey.lonely.psu = "adjust"`, which lets `survey_mean()` return *a* standard
  error, but with 1 PSU per domain that SE is a borrowed/centered surrogate, not
  a genuine design-based estimate. The admin-2 CIs are therefore not trustworthy
  for most areas.
- **The point estimate is one EA's (weighted) mean.** A single enumeration area
  is not representative of a district, so the "direct survey admin-2 prevalence"
  carries large cluster-level bias — it reflects one village, not the area.
- **Weights barely operate within a 1-cluster domain** (near-constant), so the
  weighted and unweighted area means almost coincide; weighting cannot rescue
  representativeness here.

**Implication:** the direct `svy_prev` at admin-2 should be treated as a *noisy
target*, not ground truth. This is exactly the justification for the SAE layer
(Fay-Herriot / BYM2) that borrows strength across areas, and for reporting
**leave-area-out CV** rather than in-sample agreement. Only **Sierra Leone**
(14 districts, ~4 clusters each) supports semi-direct admin-2 estimates; for the
other three, admin-2 numbers are model-dependent by necessity. (Consistent with
the project's `sae_area_level` note.)

---

## 2. Gambia uses the WRONG weight for biomarker outcomes (concrete, fixable)

The configured weight `gw_svy_weight` is **NA for 51 and zero for 28** records,
and crucially it does **not cover the blood (biomarker) subsample** — the raw
child VAD computed with `gw_svy_weight` returns **NA**, because biomarker-bearing
children lack a valid value on that weight.

Gambia carries dedicated **biomarker-subsample weights**:
`gw_c_blood_weight` (covers all 1358 children, 0 NA/0 zero) and
`gw_w_blood_weight` (covers all 1851 women). These are the correct weights for
blood-based outcomes (they incorporate differential blood-draw non-response).

- child raw VAD: `gw_svy_weight` → **NA**; `gw_c_blood_weight` → **26.6%**;
  unweighted → 30.8% (cor(weights)=0.87).

Current handling is crude and inconsistent: `compute_national_estimates()`
replaces NA/≤0 weights with **1** (`wts[is.na(wts)|wts<=0] <- 1`), while
`compute_svy_admin2()` *drops* NA-weight rows but keeps zero-weight rows. Either
way the biomarker sample is mis-weighted.

**Recommendation:** for biomarker outcomes, use `gw_c_blood_weight` /
`gw_w_blood_weight` (population-specific) for Gambia instead of `gw_svy_weight`.

---

## 3. One weight per country for BOTH children and women

Every country sets a single `weight_col` used for all outcomes. Children and
women are different analytic populations with different selection probabilities
(and, for biomarkers, different blood-draw response), so a shared weight is
loose:

- **Gambia** has population-specific blood weights (above) — clearly should be
  split by population.
- **Ghana**: configured `gw_sWeight` (range 0.47–1.35) is **nearly uncorrelated
  with the available `gw_child_weight`** (r = 0.12; mean abs diff = 1.09), and
  `gw_PSU_weight` / `gw_PSUStrat_weight` also exist. `gw_sWeight` looks like a
  generic person/household weight, not the child biomarker weight. **Needs a
  codebook check** to confirm which weight each population should use — flag, not
  yet a confirmed bug.
- **Sierra Leone**: only `gw_svy_weight` is present and it *does* vary by
  population (child 0.58–3.62, women 0.41–5.68) — probably fine.

---

## 4. Malawi ignores available stratification & uses one expansion weight

- `strata_col = NULL`, yet the data contains `mregion` (and urban/rural is
  derivable). The MNS was a stratified design; omitting strata mis-estimates
  variance (usually conservative) and discards design information.
- `svy_weight` is an **expansion weight** (range 2.6e4–5.8e6, ~225×) applied to
  **both** children and women from one column. Expansion vs normalized weights
  don't change a *mean*, but the single-column-for-both-populations pattern is
  the same concern as §3.

**Recommendation:** add `strata = mregion` (with `lonely.psu="adjust"`) for
Malawi; confirm whether separate child/women weights exist upstream.

---

## Priority of fixes

1. **Gambia blood weights** (§2) — concrete mis-weighting of every Gambia
   biomarker estimate; highest value, low effort.
2. **Frame admin-2 direct estimates as model-dependent** (§1) — ensure
   manuscript/dashboard never present 1-cluster admin-2 `svy_prev` (or its CI) as
   a reliable direct estimate; lean on SAE + leave-area-out CV. Mostly a
   reporting/》framing fix.
3. **Population-specific weights** (§3) — verify Ghana child weight vs `gw_sWeight`
   against the codebook before switching.
4. **Malawi strata** (§4) — add `mregion` strata; modest precision/correctness gain.

All four are *weighting/representativeness* issues independent of the BRINDA
work; none blocks the current regen.
