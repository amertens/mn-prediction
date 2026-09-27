# Abstract — external validation of sub-national micronutrient deficiency rankings (XV-01 / XV-02)

Drafted 22 September 2026 from `XV-01_EXTERNAL_VALIDATION_2026-09-22.md`.
All figures traceable to `results/tables/external_validation/xv_transport_pooled.csv`.

---

**Remotely sensed climate and soil rank sub-national micronutrient deficiency in
six countries absent from the training data**

**Background.** Sub-national estimates of micronutrient deficiency are needed
where biomarker surveys do not reach. Area-level models built on environmental
proxies can rank districts within surveyed countries, but whether those
rankings transport to an unsurveyed country has been assessed only by
leave-one-country-out cross-validation inside a single harmonised panel — a
test that shares data processing, inflammation adjustment and survey weighting
with the training data, and so cannot distinguish a transportable signal from a
shared analytic pipeline.

**Methods.** We fitted a parsimonious two-domain index — 59 climate variables
(TerraClimate 1991–2020 normals; MODIS land-surface temperature) and soil
properties (Africa-only iSDAsoil, or global SoilGrids v2.0) — to
survey-weighted first-administrative-level prevalence and biomarker
concentrations from four African micronutrient surveys (Gambia 2018, Ghana
2017, Sierra Leone 2013, Malawi 2015–16; 53 units). Predictors were
rank-normalised within country and reduced to domain principal components
oriented on training rows only; outcomes were z-scored within country. Held-out
rankings were scored against WHO Vitamin and Mineral Nutrition Information
System (VMNIS) sub-national deposits for six countries never used in training:
Zambia 2023, Ethiopia 2015, Sudan 2018, Nigeria 2021, Pakistan 2018–19 and
India 2016–18 (6–29 units each; 41 country–outcome–target cells). Because
cross-survey biomarker concentrations, deficiency cut-offs and inflammation
adjustments are not comparable, the estimand was the within-country Spearman
correlation, which is invariant to any country-constant offset. Inference used
a country-block permutation null that preserves correlation among a country's
outcomes.

**Results.** Across the four African countries the index ranked held-out units
better than chance on both targets (biomarker concentration: mean ρ = 0.402,
12 of 12 cells positive, p = 0.001; prevalence: ρ = 0.332, 17 of 21, p < 0.001),
above the internal leave-one-country-out estimates previously obtained at the
same tier (0.281 and 0.267). Substituting global SoilGrids — which contains no
plant-available micronutrients — for iSDAsoil cost nothing across 33 paired
cells (ρ = 0.374 vs 0.358; SoilGrids superior in 17 of 33). With that
substitution the index transported off-continent to South Asia (prevalence:
ρ = 0.394, 7 of 8 cells positive, p < 0.001), including the two best-powered
cells in the study (India, child vitamin A ρ = 0.497 across 27 states, p = 0.006;
child iron ρ = 0.366 across 28 states, p = 0.032). Iron, vitamin B12 and folate
carried the African signal; vitamin A did not, yet was the strongest outcome in
South Asia.

**Conclusions.** Rankings of sub-national micronutrient deficiency derived from
environmental proxies generalise beyond the surveys, the harmonisation pipeline
and the continent on which they were fitted. Soil micronutrient layers are not
the active ingredient, so the approach requires only globally available
rasters. Validation was restricted to first-administrative-level published
aggregates without design effects, and the off-continent arm to prevalence, so
absolute burden estimates remain untested; published sub-national VMNIS records
nonetheless constitute an under-used external validation resource obtainable at
no additional data-collection cost.

---

*Word count (Background–Conclusions): ~390.*

## Notes for whoever uses this

- The comparison to "0.281 and 0.267" is a **benchmark, not a paired test**:
  those internal figures come from `admin1_transport.csv`, fitted on the older
  `predictors_admin2_shared.csv` vocabulary, whereas the external figures use
  the re-extracted common vocabulary. The soil comparison (0.374 vs 0.358) *is*
  paired, on the same 33 cells and the same folds.
- The "41 cells" is 33 African (12 level + 21 prevalence, iSDA) plus 8
  off-continent (prevalence only).
- Zinc is absent throughout: Malawi is its only training country, so it cannot
  transport.
- If a word limit bites, the first sentence of Background and the final clause
  of Conclusions are the most compressible.
