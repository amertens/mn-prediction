# XV-01. The climate + soil index ranks sub-national units in four countries it has never seen

22 September 2026. Scripts `scripts/external_validation/01`–`04`; results in
`results/tables/external_validation/`.

## The question

Every transport number this project has published is leave-one-country-out
*inside* the four-country panel — our harmonisation, our BRINDA adjustment,
our weighting, our four surveys. The claim being made of it is that the
rankings carry to a country with no survey at all. That had never been scored
against a country outside the panel.

The one prior external check was
`scripts/policy_deck/24_civ_2007_survey_check.R` — women's B12 over nine Côte
d'Ivoire eco-regions, rank agreement 0.95, prompted by S. Hess at the 18
September check-in. One outcome, nine zones, a compass crosswalk, a 2007
survey. Encouraging, not evidence.

## The data that made it possible

WHO VMNIS deposits sub-national rows that the project's national-only pull
(`scripts/pull_vmnis_validation.R`, `Representativeness == "national"`)
discards: **6,803 rows at "1st administration level"** and 3,665 "regional
(within country)", against 9,761 national. The earlier plan doc
(`docs/brinda_vmnis_loco_validation_plan.qmd`) concluded sub-national VMNIS was
out of scope because region names are free text. That is wrong — the
per-indicator exports carry a structured `Representativeness.` field, a
`Representativeness name`, and a `SurveyId`. The prevalences were also missed
because each indicator names its own column (`Depleted iron stores
prevalence`, `Inadequacy prevalence`), not the generic
`Prevalenceofdeficiency`: 787 of 866 sub-national ferritin rows carry a
prevalence once you look in the right column.

Four countries clear the bar — enough units, populated prevalence, our
populations, and in Africa so the soil half (iSDAsoil) exists:

| Country | Survey | Units | Cells |
|:---|:---|---:|---:|
| Zambia | 2023 NFNC/TDRC | 9 provinces | 8 |
| Ethiopia | ENMS 2015 | 11 regions | 7 |
| Sudan | 2018 | 15 states | 4 |
| Nigeria | NFCMS 2021 | 6 zones | 8 |

Nigeria's NFCMS 2021 is deposited at the six geopolitical zones, not the 36
states, which is the single biggest missed opportunity in the set.

**BRINDA cannot do this.** Its pooled outputs are national by country × year,
and nothing is on disk (`archive/national_prediction/data/` holds a README).
Sub-national BRINDA needs a consortium request.

## Why a ranking test is the right test

Cross-survey biomarker *levels* are not comparable — raw ferritin spans six-fold
across our own four surveys (`fe_transport_level_offset`) — and VMNIS cut-offs
and adjustments vary between deposits. A **within-country Spearman** is
invariant to any country-constant offset and to any assay or cut-off choice
held fixed inside one survey. It is the one quantity that survives the
comparison, and it happens to be exactly the quantity the project claims.
Nothing here tests levels or absolute prevalence.

## Design

Training is the four panel countries at their admin-1 rung (53 units); the test
country is held out entirely. Following `16_admin1_transport.R`: predictors
area-averaged as a simple mean of districts, rank-normalised within country,
domain PCs to 80% variance oriented from the training rows only; outcome
z-scored within country.

Covariates are **re-extracted for both sides by one extractor**
(`03_extract_climate_soil.py`, 1,602 districts, 59 climate + 64 soil columns).
The existing `predictors_admin2_shared.csv` climate/soil columns come from an
older analyst export whose scaling is not fully documented, and the pooled
model matches covariates by name while treating a mis-scaled name as
comparable — the failure WS-G caught. Extracting both sides with the same code
makes them commensurate by construction. Parity was checked, not assumed:
`pr_ann_mean` reproduces the existing `clim_pr_ann_mean` **exactly** (Gambia
854.9, Ghana 1260.2) and iSDA zinc lands within 3% of `soil_zinc_mean_0_20`
(Gambia 1.16 vs 1.129, Ghana 2.29 vs 2.243), the residual being polygon
simplification.

## Result

| Target | Cells | Mean ρ | Positive | Country-block null 95th | p |
|:---|---:|---:|---:|---:|---:|
| **level** | 12 | **0.402** | **12 / 12** | 0.201 | **0.0010** |
| **prevalence** | 21 | **0.332** | 17 / 21 | 0.176 | **0.0005** |

The elastic net on the same domain axes is well behind on the level (0.194,
8/12, p = 0.046) and close on prevalence (0.290, p = 0.003) — the same ordering
the panel data shows.

**Two nulls are reported because they assume different things.** The
independent null draws one permutation per cell and is too generous, because a
country's outcomes are measured on the same units and are correlated. The
country-block null draws one relabelling of a country's units per replicate and
applies it to every outcome of that country. It is duly fatter (0.201 vs 0.150
on the level; 0.176 vs 0.129 on prevalence) — the correlation concern was
real — and the index clears it either way.

For scale, the panel's own internal LOCO at the same tier is 0.281 (level) and
0.267 (prevalence). **The external numbers are higher, not lower.**

### Which outcomes carry it

Iron, B12 and folate transport; vitamin A does not.

- Zambia women's B12 level **0.800**, women's iron **0.650**, child iron **0.583**
- Sudan women's iron level 0.554, women's vitamin A level 0.540, child iron 0.489
- Ethiopia women's folate prevalence **0.873**, women's B12 0.536
- All four negative cells are vitamin A or a thin Nigeria cell: Zambia child
  vitamin A prevalence −0.427, Ethiopia women's vitamin A prevalence −0.428,
  Sudan women's iron prevalence −0.064, Nigeria women's folate −0.086

This matches the panel's own finding that child vitamin A shows no reliable
geography on multi-cluster units in three of four countries (P7's basis), and
it matches DC-H2's warning that vitamin A prevalence sits near the floor for
women in several surveys. Nigeria's child-iron prevalence ρ = 1.000 is real but
rests on six units.

## What this does and does not settle

**Does.** The ranking claim survives contact with countries outside the panel
and with labels this project did not make, harmonise or weight. The central
falsification condition in
`PREREGISTRATION_NEW_COUNTRIES_2026-09.md` — "if P1 fails, the transport result
was a property of West/Southern Africa rather than of the proxies" — did not
trigger; Ethiopia, Sudan and Zambia are not West Africa and Sudan is a
different agro-ecology altogether.

**Does not.**

1. **This is not the pre-registered P1 scoring.** P1 is specified on a new
   country's *own* survey microdata at the district rung and the regional tier,
   with effective n from the measured design effect. This is admin-1 only, on
   published aggregates, with no design effect available and no precision
   weighting of the test country's units. It is an independent external test
   that is *consistent with* P1, not a substitute for it.
2. **P2 is untested here.** Only the climate + soil index was run; the full
   24-domain index cannot be built for these countries without the other
   domains, so the comparison P2 makes has no second arm.
3. **Boundary vintage is reconstructed, not observed.** Zambia's Muchinga
   (2011) is returned to Northern except Chama; Sudan's Central Darfur, East
   Darfur and West Kurdufan (2012–13) split back to pre-split parents; Abyei is
   dropped. Each is recorded in `gadm_to_vmnis_crosswalk.csv` with a note, but
   each is a judgement.
4. **Six units is not a test.** Nigeria cannot clear its own permutation null
   at any effect size below ρ ≈ 0.77. It contributes to the pooled result and
   should not be quoted alone.
5. **Survey vintage.** Sudan 2018–19 and Ethiopia 2015 are matched to
   survey-year climate, but the soil and climatology blocks are long-run
   normals, so vintage matters less here than it would for a level claim.

## Files

- `scripts/external_validation/01_build_vmnis_targets.py` → `data/external_validation/vmnis_admin1_targets.csv` (257 labels, 27 cells)
- `scripts/external_validation/02_build_crosswalk_geoms.R` → crosswalk + geometry
- `scripts/external_validation/03_extract_climate_soil.py` → 1,602 districts × 123 columns
- `scripts/external_validation/04_transport_test.R` → `xv_transport.csv`, `xv_transport_summary.csv`, `xv_transport_pooled.csv`
- Survey windows added to `metadata/survey_years.csv` with `in_protocol = FALSE`

---

# XV-02. The off-continent test: Pakistan and India, with a substituted soil block

Added 22 September 2026, same scripts.

## Why a substitution was needed, and why it needed a control

iSDAsoil is Africa-only, so the pre-registered soil half cannot be built for
Pakistan or India. The substitute is **global SoilGrids v2.0** (10 properties ×
4 depths = 40 columns: pH, organic carbon, nitrogen, CEC, clay, sand, silt,
bulk density, coarse fragments, organic carbon density).

That substitution is **not** obviously neutral: SoilGrids carries no
plant-available micronutrients — no Zn, Fe, Ca, Mg, P, K, S — which are
mechanistically the interesting half of iSDA for a micronutrient outcome. A
weak off-continent result would then be unreadable: transport failing abroad,
or just a worse soil block? So both blocks were run on the four **African**
countries, on the same cells, to price the substitution before using it.

## The substitution is free

Paired over the same 33 African cells:

| Soil block | Mean ρ | Better in |
|:---|---:|---:|
| iSDAsoil (Africa-only, with micronutrients) | 0.358 | 16 / 33 |
| SoilGrids v2.0 (global, no micronutrients) | **0.374** | 17 / 33 |

A dead heat, if anything favouring the global block (+0.016). Pooled, on the
index: level 0.433 vs 0.402, prevalence 0.340 vs 0.332 — SoilGrids is never
worse.

**This is a finding in its own right, and it is not a comfortable one for the
mechanistic story.** The plant-available soil micronutrient layers that only
iSDA has contribute nothing the generic physical and chemical properties do not
already carry. Soil zinc is not what makes the soil block work. That is
consistent with the project's existing reading that the transported axis is
rural-subsistence agro-ecology rather than a soil-to-diet nutrient pathway
(`transport_domains_climate_soil`), and it removes the Africa-only constraint
from the deployable recipe.

## The off-continent result

Prevalence only — neither deposit carries biomarker means, so there is no level
target. India's CNNS sampled preschool/school-age children and adolescents but
**no women**, so only child outcomes can be scored there; child zinc drops
because Malawi is its only training country.

| Set | Cells | Mean ρ | Positive | Country-block null 95th | p |
|:---|---:|---:|---:|---:|---:|
| **Off-continent, pooled** | 8 | **0.394** | 7 / 8 | 0.225 | **0.0005** |
| India | 2 | 0.431 | 2 / 2 | 0.221 | 0.0015 |
| Pakistan | 6 | 0.381 | 5 / 6 | 0.290 | 0.0065 |

For comparison, the African arm on prevalence is 0.340. **The off-continent
number is not lower.**

Per cell:

| Country | Outcome | Units | ρ | own-cell p |
|:---|:---|---:|---:|---:|
| India | child vitamin A | 27 | **0.497** | **0.006** |
| India | child iron | 28 | 0.366 | 0.032 |
| Pakistan | women's folate | 8 | 0.690 | 0.044 |
| Pakistan | child vitamin A | 8 | 0.619 | 0.058 |
| Pakistan | women's iron | 8 | 0.500 | 0.101 |
| Pakistan | women's vitamin A | 8 | 0.405 | 0.168 |
| Pakistan | child iron | 8 | 0.333 | 0.218 |
| Pakistan | women's B12 | 8 | **−0.262** | 0.746 |

India's two cells are the **only cells in the whole study that clear their own
individual permutation null**, because 27–29 states is the only place the
design has real per-cell power. Everything else rests on pooling.

**A reversal worth noticing.** Vitamin A is the outcome that *failed* in Africa
(three of the four negative cells there). Off-continent it is the *strongest*
(India 0.497, Pakistan 0.619). No explanation is offered here; it is recorded
because it cuts against reading the African vitamin A failure as a property of
vitamin A itself.

## What XV-02 does and does not settle

**Does.** The transported ranking is not a property of Africa. It survives to
South Asia, on a different continent, a different agro-ecology, different
surveys, different assays, and a soil block with no overlap in its
micronutrient content. And the Africa-only constraint on the recipe is lifted:
global SoilGrids performs as well as iSDA everywhere it was tested.

**Does not.**

1. **Prevalence only, 8 cells, 2 countries.** The level target — where the
   African arm was strongest — cannot be scored off-continent at all.
2. **Pakistan's 8 units are at the floor.** No Pakistani cell clears its own
   null at conventional levels, and its mean rests on the pooled test.
3. **India contributes 2 cells**, so its sign test is uninformative (0.25 is
   the smallest attainable p with two cells); its strength is unit count, not
   cell count.
4. **This is still not the pre-registered P1 scoring**, for the same reasons
   given for XV-01: admin-1, published aggregates, no design effect, no
   precision weighting.
5. **The soil comparison is on African cells only.** That the substitution is
   free in Africa is evidence, not proof, that it is free in South Asia — there
   is no iSDA there to check against, by construction.
