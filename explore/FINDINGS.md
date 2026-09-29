# `explore/` findings log

Running record of the hypothesis-generating probes in this folder. One entry
per probe: the question, the design, the result, the honest number next to the
loose one, and a verdict.

**Verdicts.** *candidate* = worth pre-registering for a country these models
have never seen · *dead end* = tried, does not work, recorded so nobody
re-tries it · *needs data* = the idea is untested because the data is not here.

**Everything here is post hoc on the same four countries** unless an entry says
otherwise. See `README.md` for the isolation contract and the reference numbers.

---

## GATE (2026-09-28) — the harness reproduces the record

**Question.** Can a generalised scorer, which accepts an arbitrary arm and an
arbitrary feature matrix, reproduce the numbers in
`results/tables/protocol_v2/benchmarks_v2_cells.csv`?

**Design.** `explore/R/harness.R` reuses the project's own `exp_cell`
construction, `prep_predictors_v2` (fix 3), `domain_representation_v2` (fix 4),
`make_folds_v2` seeds (fix 1) and `score_v2`. Run the protocol's own arms
(`null_train_mean`, `spatial`, `domain_index`) over all cells, both targets,
all three estimands, and compare cell by cell against the record.

**Result. CLEAN.** 98.8% of 376 cell × arm comparisons within 0.02; mean
|diff| 0.0009. In-fill and region reproduce to three decimals
(level 0.394 / 0.457 / −0.199; prev 0.298 / 0.307 / −0.195). LOCO transport
reproduces **exactly** (level 0.287, prev 0.192, max |diff| 0.0000). Residual
differences are confined to the constant-prediction null arm under the region
estimand on prevalence, where Spearman is tie-dominated.

**Two traps it caught**, both of which would have silently corrupted every
probe downstream:

1. **Predictor tiers.** The headline record runs
   `V2_PREDICTOR_TIERS = "open,survey_public"` (no DHS) for *all three*
   estimands, not just transport. A probe on the default (all tiers) is
   computed on a bigger predictor set than the record it is compared to.
   `exp_load()` now defaults to the headline set.
2. **Where the domain axes are built under LOCO.** Building them per country
   and intersecting the columns lets each country learn its own PC1 sign, so a
   domain score means the opposite thing in two countries. The first harness
   did this and lost up to **0.76 Spearman** on a single cell (Gambia
   women_iron prev: −0.236 against the record's +0.524). The axes must be
   rebuilt inside each fold on the pooled matrix with
   `sign_rows = training countries`, as `02b_merge_and_loco.R` does.

**Verdict.** Gate passed; probes may proceed. Trap 2 is worth carrying into
any future cross-country work in the main pipeline too — it is a silent,
sign-flipping failure that looks like a null result rather than a bug.

---

## RV-01 (2026-09-28) — the cross-survey level offset is not a predictor-side batch effect, and it is nutrient-specific

**Question.** The project's transport estimand is a rank claim only, because
biomarker levels carry large cross-survey offsets (raw ferritin 6× across
countries; AS-01 rules out the assay, and it is not the adjustment method).
Genomics calls this a batch effect and removes it by projecting out the
unwanted-variation subspace rather than discarding the signal that shares a
scale with it. Does that work here?

**Design.** `explore/scripts/03_ruv_transport.R`. Two predictor-side
projections, both learned on training rows and never using the outcome:
`ruv_bc_k` projects out the top-k between-country discriminant directions
(whitened between-class scatter, k ≤ 3 with four countries); `ruv_pc_k`
projects out the top-k principal components as the unsupervised comparator.
Scored on LOCO transport (22 cells) and, as the control, on within-country
in-fill — a genuine batch correction should help across countries and be
roughly neutral within one.

**Result 1 — predictor-side RUV is provably a no-op, not merely ineffective.**
`ruv_bc` reproduces the index to three decimals (transport level 0.297 vs
0.298; prevalence 0.193 vs 0.196). The reason is structural, and worth stating
because it forecloses a whole family of ideas: **within-country rank
normalisation (the project's fix 3) maps each column to
`qnorm((rank − 0.5)/n)` inside each country, so every column has mean ≈ 0 in
every country — and therefore so does every linear combination of columns.**
Checked directly on the pooled 206 × 329 matrix: the largest country-mean
deviation is 0.058, and country is unpredictable from PC1–PC4
(F = 0.0, η² = 0.000 for all four). The between-class scatter is identically
zero, so there is no linear direction to remove. Fix 3 already does everything
a linear predictor-side batch correction could do.

**Result 2 — removing principal components destroys signal, in both
directions.** `ruv_pc_k1` / `k3` fall to 0.164 / 0.207 on transport (from
0.298) and to 0.153 / 0.143 on within-country in-fill (from 0.394). Hurting
*both* estimands is the signature of removing signal rather than batch: the
dominant directions of the predictor matrix are the agro-ecological gradient
the whole approach rests on.

**Result 3 (the useful one) — the level offset is concentrated in iron and
folate, and is small for B12 and vitamin A.** Between-country share of the
level variance, per outcome:

| outcome | countries | between-country share | range of country means | mean within-country sd |
|---|---|---|---|---|
| child_iron | 4 | **0.804** | 1.607 | 0.327 |
| women_iron | 4 | **0.710** | 1.032 | 0.304 |
| women_folate | 3 | 0.577 | 0.791 | 0.345 |
| women_vitA | 4 | 0.349 | 0.221 | 0.131 |
| child_vitA | 4 | 0.280 | 0.138 | 0.103 |
| women_b12 | 3 | **0.143** | 0.282 | 0.362 |

**What this suggests.** The rank-only restriction on transport is a *global*
response to a problem that is largely an *iron and folate* problem. For B12 —
where between-country variance is 14% and within-country spread is the largest
of any outcome — transporting a LEVEL, not just a ranking, may be feasible.
That is a concrete pre-registrable hypothesis, and it lines up with B12
already being the strongest cell on the record (Malawi 0.70).

**Verdict.** *dead end* for predictor-side RUV — and worth keeping as a dead
end, because the reason is a mathematical identity that rules out the whole
family. *candidate* for outcome-side level transport restricted to B12 and
vitamin A. Caveat: three countries for B12 and folate, and these shares are
computed on the same four surveys that would be used to fit any correction.

## TM-01 phase 1 (2026-09-28) — seasonal phase is a coin flip; the lag hypothesis is untested

**Question.** Fieldwork windows are only 2–3 months per country, so calendar
spread is not the temporal axis available. What varies across clusters is
*seasonal phase*: two clusters visited in the same week sit at different points
in their own agricultural year because their rainy seasons peak in different
months. Since ferritin and retinol integrate months of intake, phase should
matter.

**Design.** `explore/scripts/09_temporal_phase1.R`, at the survey cluster
(323 clusters with GPS), folds cut by district so a cluster's own district
never trains on it. Phase was constructed — it does not exist in the cluster
table — as the circular months from each dynamic layer's peak month to the
cluster's fieldwork month, plus its sin/cos encoding (21 columns from 7 layers).
Arms are ridge on climatology, + the fieldwork-window block, + phase, and phase
alone. Scored at cluster level and aggregated to districts, 24 cells.

**Result. Phase adds nothing on average.** Mean gain over the fieldwork block
is **+0.011, improving 12 of 24 cells** — a coin flip. Phase alone scores 0.024
at cluster level and −0.003 at district level: it carries no signal by itself.
The fieldwork-window block does add (+0.026, 16 of 24), replicating FW-01's
+0.02 on independent folds.

**One pattern worth recording, with a confound attached.** The cells where
phase helps most are Malawi women_zinc (+0.095; phase alone scores **0.233**
there, against climatology's −0.128), Sierra Leone women_vitA (+0.075), Gambia
child_vitA (+0.059), Sierra Leone women_folate (+0.058) and Malawi child_zinc
(+0.048). Malawi zinc leading this list is suspicious rather than promising:
ZN-02 established a *collection artefact* in Malawi serum zinc (afternoon draw
−3.5%, date +0.7%/day), and phase is a deterministic function of the fieldwork
month. So the strongest "seasonal" effect in the table is the cell where a
date artefact is most plausible. Any follow-up must separate season from
collection date before reading this as biology.

**What this does and does not rule out.** It tests *phase*, one of the three
temporal ideas. It does not test the mechanistically stronger one — a **lag
stack** over the 0–24 months before the draw, so the *previous growing season*
is represented rather than the 3-month window the cluster table already
carries. That needs monthly rasters spanning 2013–2018; the on-disk rasters
cover only parts of 2014–15, so it is only reachable through the staged Earth
Engine extraction (`10_gee_cluster_monthly.py`).

**Verdict.** *dead end* for seasonal phase as a predictor. *needs data* for the
lag hypothesis. The phase result is mild evidence against a large temporal
signal at this resolution, which should temper expectations for the lag test
rather than cancel it.

## AE-01 (2026-09-28) — AlphaEarth as a kernel: the hypothesis was wrong

**Question.** The 64 `aef_*` dimensions are a learned embedding whose inner
product encodes similarity. Scored as "a domain" the record puts them at
±0.012, i.e. nothing — but the domain representation collapses them through
their leading principal components, which for a learned embedding ought to be
the one representation that destroys the information. Does the embedding carry
signal when used as a **kernel** instead?

**Design.** `explore/scripts/02_alphaearth_kernel.R`. BLUP on the linear and
cosine kernels over the 64 dimensions, against `aef_index` — the current
representation, the index run on the embedding's own domain axes — plus
climate+soil alone and embedding-plus-climate+soil, on identical folds.

**Result. No.** On medians the kernel looked better (in-fill level 0.405 vs the
PC representation's 0.378), but the **paired** comparison, which is the right
one because both arms see the same cells and the same folds, reverses it:

| comparison | wins | mean difference |
|---|---|---|
| kernel vs PC, in-fill level | **4 of 18** | — |
| kernel vs PC, transport level | 12 of 22 | +0.016 |
| kernel vs PC, transport prevalence | 9 of 22 | −0.023 |
| adding the embedding to climate+soil, transport level | **4 of 22** | −0.020 |
| adding the embedding to climate+soil, transport prevalence | 5 of 22 | −0.009 |

The representation is not what was holding AlphaEarth back. And adding the
embedding to climate+soil actively *hurts* transport, in 18 of 22 cells.

**A methodological note worth carrying.** The median across cells and the
paired per-cell comparison disagreed in direction here. With 18–22 cells of
very unequal difficulty, the median is the wrong summary for "does A beat B" —
every arm comparison in this folder should be read paired.

**One thing that survives.** `cs_linear` — a single BLUP kernel on climate and
soil — reached **0.469** on in-fill level, above the spatial smoother (0.457)
and well above the domain index (0.394), and transports at 0.360. The kernel
machinery is worth keeping; the embedding is not what it should be fed.

**Verdict.** *dead end* for the embedding, by a well-powered paired test rather
than a weak one. The record's ±0.012 was right and for the right reason.

---

## MX-01 (2026-09-28) — a literature-grounded, nutrient-specific feature set does not beat 383 generic columns, and shows no nutrient specificity

**Question.** The predictor set is nutrient-agnostic: the same columns are
offered to zinc, B12, folate, iron and vitamin A. The literature is not. Do
~10 mechanistic features built from the crop basket beat 383 generic ones?

**Design.** `explore/scripts/06a` + `06` + `07`. The raw 42-crop MapSPAM
production grid (the pipeline currently sees only 5 collapsed group shares)
crossed with a documented food-composition table
(`explore/data/food_composition.csv`) to derive, per district: the
**phytate:zinc molar ratio** of the production basket (the Wessells & Brown
mechanism behind national zinc-deficiency estimates), **provitamin-A carotenoid
density** and **oil-palm share** (crude red palm oil is the dominant West
African plant source), phytate load on non-haem iron, folate density, and
composition axes the 5-group version cannot express. 193 of 206 target
districts join; the 13 misses are urban wards and Malawi Boma towns with no
cropland, where a production basket is undefined.

The features pass an agronomy sanity check: carotenoid density peaks on Ghana's
oil-palm belt (Wassa East, Mpohor, Birim North, Akyemansa), phytate:zinc on
Gambia's groundnut-and-millet Foni districts, cassava share in Ghana's Central
region. They are measuring what they claim to measure.

**The control that makes this readable.** Every cell was scored twice: with the
features matched to its own nutrient, and with a **mismatched** set (the zinc
cell gets the vitamin A features, and so on). A generic agro-ecology proxy
would do equally well either way.

**Result. A clean negative, on both counts.** In-fill, level target, median
Spearman over cells:

| arm | median |
|---|---|
| spatial | 0.457 |
| index + mechanistic | 0.394 |
| domain index | 0.392 |
| all 16 mechanistic | 0.140 |
| **mechanistic, mismatched nutrient** | **0.075** |
| **mechanistic, matched nutrient** | **0.058** |

The matched set beats the mismatched set in **6 of 18 cells** — below chance.
There is no nutrient specificity. And adding the features to the index moves it
from 0.392 to 0.394, i.e. nothing.

The single exception is Ghana women_folate (matched 0.264 vs mismatched 0.147,
+0.117), where the mechanism is pulse share. One cell out of eighteen is what
this design produces by luck.

**Limitations that do not rescue it.** Production is not consumption; MapSPAM is
2010 against surveys from 2013–18; the composition table is indicative rather
than a country-specific food-composition table; and ~10 columns is low
capacity. But none of these explain the *specificity* failure: the wrong
nutrient's mechanism does as well as the right one, which is what the mismatch
control exists to detect.

**Verdict.** *dead end.* Worth recording as a strong negative rather than a
weak one, because the feature set was built carefully, validated against
agronomy, and tested with a specificity control. It is evidence that the
district-level signal is not dietary composition — consistent with the
project's own finding that the replicated covariate associations run opposite
to individual-level nutrition for legumes and cattle, a rural-subsistence axis
rather than a diet axis.

## MT-01 (2026-09-28) — multi-trait BLUP: the first real positive

**Question.** XO-01 found cross-nutrient borrowing null under transport but the
same nutrient in the other population beating the index. That is low-rank
structure across the 24 cells, and the project fits every cell on its own.
Quantitative genetics fits correlated traits jointly and gains most exactly
where single-trait estimation is noisiest — which at n = 14–87 is everywhere.

**Design.** `explore/scripts/05_multitrait.R`. A held-out district has no
biomarker measured at all, so the model cannot condition on its other
nutrients; the gain has to come from estimating the *shared spatial component*
on more data. On training districts only, the first eigenvector of the trait
correlation matrix gives loadings, and their weighted mean is a general
deficiency factor `g`. `mt_blup` = BLUP of `g` on the climate+soil and spatial
kernels, scaled by the trait's loading, plus a trait-specific BLUP of the
residual. `st_blup` is the same kernels, single trait — the controlled
contrast. Held-out districts are held out for every trait simultaneously, and
the loadings and `g` are computed on training rows only.

**Result. Multi-trait wins, and for the predicted reason.**

| | mean gain | wins | Wilcoxon p |
|---|---|---|---|
| level | +0.022 | 13/18 | 0.027 |
| prevalence | +0.056 | 13/18 | 0.005 |
| pooled | **+0.039** | 26/36 | 0.0004 |

Cells within a country share districts and the same general factor, so those 36
comparisons are not independent. At the **country × target block** level, which
is the honest unit: **6 of 6 blocks positive**, sign p = 0.031, mean of block
means +0.038, and every country agrees (Gambia +0.034, Ghana +0.037, Malawi
+0.043).

**The shrinkage signature is present**, which matters more than the p-value.
Theory says borrowing should help most where the single-trait fit is weakest,
and it does: Spearman(single-trait score, gain) = **−0.388, p = 0.020**; by
tercile of single-trait performance the mean gain is **+0.086 (weak)**, +0.013
(mid), +0.018 (strong). The largest gains are all rescues of cells that were
worse than useless alone — Malawi women_iron prevalence −0.137 → +0.042, Ghana
women_vitA prevalence −0.037 → +0.133, Malawi child_zinc prevalence −0.180 →
−0.020. Where the single-trait fit was already strong (the B12 cells) it
changes nothing, −0.004 to −0.006.

**Trait correlations, which are the mechanism.** Within country the outcomes
are strongly correlated, including *across* nutrients: Gambia women_vitA ~
child_vitA 0.837, women_iron ~ child_vitA 0.781, women_iron ~ child_iron 0.772;
Malawi women_zinc ~ child_zinc 0.552; Ghana women_b12 ~ child_iron 0.382. This
is not in tension with XO-01's null — that was cross-nutrient borrowing under
*transport*, where the country offset intervenes. Within a country, the
nutrients share district-level structure and it is usable.

**Limits.** Four countries, and Sierra Leone is absent because its 14 districts
cannot be folded for in-fill. The kernels were fixed to climate+soil plus
space, chosen from the record rather than tuned here, which is a virtue for
honesty and a constraint on the ceiling. This was not pre-registered.

**Verdict.** *candidate*, and the strongest so far. Concretely: when a survey
measures several biomarkers, fit them jointly with a shared latent factor
rather than one cell at a time — and expect the gain in the weak cells, which
are the ones currently reported as failures.

## EB-01 (2026-09-28) — cross-cell moderation does not improve the index, but the screen it produces is worth keeping

**Question.** The index sets its weights from one cell's 14–87 districts, and
`fe_effective_n` names exactly that as the root cause of unstable prediction.
There are 24 cells sharing one predictor vocabulary — the many-features ×
many-contrasts layout empirical-Bayes moderation was built for. Does shrinking
each cell's weights toward the cross-cell consensus help?

**Design.** `explore/scripts/04_moderated_meta_index.R`. The index's own
statistic `z = fisher-z(spearman) × sqrt(n−3)` has unit variance by
construction, so the hierarchy is textbook: `z_cj ~ N(θ_cj, 1)`,
`θ_cj ~ N(μ_j, τ²_j)`, posterior mean `μ_j + τ²/(τ²+1)·(z_cj − μ_j)`, with
μ and τ² estimated from the *other* cells only. Cells of the same country share
districts, so companion rows are always filtered to exclude the held-out
district keys — otherwise the weights would read the held-out district's own
survey through another biomarker. Under LOCO the whole held-out country is
dropped from the companions.

**Result 1 — the estimator is null.** Paired against the index:

| estimand | mean gain | wins | sign p |
|---|---|---|---|
| in-fill level | −0.017 | 6/18 | 0.24 |
| in-fill prevalence | −0.009 | 8/18 | 0.82 |
| transport level | +0.009 | 11/22 | 1.00 |
| transport prevalence | +0.003 | 12/22 | 0.83 |

Moderation neither helps nor hurts. The index's per-cell weights are evidently
not the binding constraint.

**Result 2 (the keeper) — a cross-cell replicated predictor screen.** 173 of
329 predictors clear FDR < 0.05 pooled over 24 cells. The domain composition
**independently corroborates the climate+soil transport finding by a completely
different route** — a meta-analysis of marginal associations rather than
leave-one-country-out ablation:

| domain | predictors at FDR < 0.05 |
|---|---|
| Climate and weather | 37 |
| Satellite embedding | 35 |
| Soil characteristics | 28 |
| Ecosystem productivity/greenness | 11 |
| everything else (15 domains) | 62 |

The single strongest predictor in the study is **`glw_ruminant_share`**
(pooled z = 6.27, FDR 1.2e-07, sign consistency 0.92). Positive z means *more*
deficiency, so more ruminants goes with worse status — the counter-intuitive
direction the project already documented as a rural-subsistence axis rather
than a diet axis. `glw_cattle_km2` is also in the top 22 (z = 4.64). Otherwise
the head of the list is temperature range and variability, vapour-pressure
deficit, soil magnesium, organic carbon, pH and nitrogen.

**The two results together make a methodological point worth carrying.** The
satellite embedding contributes **35 predictors at FDR < 0.05 with sign
consistency 0.79–0.88** — among the most replicated associations in the whole
set — and AE-01 showed it adds **nothing** predictively and actively hurts when
added to climate+soil. Replicated marginal association and incremental
predictive value are different things, and at p = 329 with heavy collinearity
they come apart sharply. Any prioritisation built from a signal scan (the RA
annotation work, for instance) should be read with that in mind.

**Verdict.** *dead end* for the moderated estimator. *candidate* for the screen
as a prioritisation tool, with the caveat above. `explore/out/04_moderated_axis_stats.csv`
carries the full ranked list with pooled z, τ², FDR and sign consistency.

## TM-01 phase 2 (2026-09-28) — the previous-growing-season hypothesis is rejected, and the lag profile says why

**Question.** The mechanistically strongest temporal idea: micronutrient stores
integrate months of intake, so the quality of the **previous growing season**
(lags 9–15 months before the blood draw) should predict status better than the
long-run climatology of the pixel and better than the 3-month window the
cluster table already carries. Phase 1 could not test this — the monthly
rasters on disk cover only parts of 2014–15 against surveys from 2013–18 — so
the Earth Engine extraction was run for it.

**The extraction.** `explore/scripts/10_gee_cluster_monthly.py`: 36 monthly
values for CHIRPS rainfall, MODIS NDVI and land-surface temperature, and FLDAS
soil moisture and evaporation, at all 323 cluster buffers (2 km urban / 5 km
rural, the cluster track's own convention), aligned so lag 0 is each cluster's
own fieldwork month. One `reduceRegions` per layer over a 120-band monthly
image rather than one call per month — five Earth Engine calls, 56,470 rows.
Values verified against expectation: CHIRPS 110 mm/month with 11% dry-season
zeros, LST 304.7 K, NDVI 0.45, soil moisture 0.28 m³/m³, evaporation
2.6e-05 kg m⁻² s⁻¹ over 10,326 distinct values.

Each value is an **anomaly** against that cluster's own mean for the same
calendar month, because raw monthly values are dominated by the season the
survey happened to run in, which is near-constant within a country.

**Result 1 — no representation of the lag stack beats climatology.** Paired,
district scoring, 24 cells:

| arm | columns | mean gain | wins | sign p |
|---|---|---|---|---|
| all 36 lags per layer | 180 | +0.037 | 14/24 | 0.54 |
| 8 interpretable windows | 40 | +0.021 | 15/24 | 0.31 |
| 5-df penalised distributed lag | 25 | +0.010 | 13/24 | 0.84 |
| **lags 9–15 only** | 5 | **+0.004** | 15/24 | 0.31 |

The narrow hypothesis, stated as sharply as it can be, adds **+0.004**. Alone
it scores 0.081 against climatology's 0.347.

**Result 2 — the lag profile is flat, which is the real finding.** Mean
|marginal correlation| by lag over all cells and layers:

- lags 9–15 (previous growing season): **0.160**
- every other lag: **0.151** (difference +0.009, p = 0.075)
- lags 24–35, two to three years before the draw: **0.158**
- profile range across all 36 lags: 0.117 to 0.186, sd 0.018
- **argmax is at lag 31**, not in 9–15
- no decay with lag (Spearman(lag, |r|) = 0.183, p = 0.29)

A real temporal mechanism predicts a peak at 9–15 and decay away from it. What
is there instead is a flat profile in which weather from *three years before
the blood draw* is as associated with status as weather from the last growing
season. That is the signature of a **spatial** association showing through
every lag — persistent differences between places — not a temporal one. It is
the same conclusion the domain work reached from the other direction: the
signal is agro-ecological regime, not recent conditions.

**Limitation.** The anomaly baseline is each cluster's own mean over only three
years, so it removes the cluster's level but leaves the regional anomaly field,
which is spatially correlated. A longer baseline (the extraction fetches 120
months and currently writes 36) would sharpen the anomalies but cannot change
a profile that peaks at lag 31.

**Verdict.** *dead end*, and a well-powered one: the hypothesis had a specific
prediction about *where* in the lag profile the association should sit, and the
data contradicts it rather than merely failing to confirm it. Taken with phase
1, the whole fine-temporal direction is closed at this resolution — which is
worth knowing, because it was expensive to reach and would otherwise stay open.

## KB-01 (2026-09-28) — kernel BLUP helps where the index is weakest: transport, not in-fill

**Question.** n = 14–87 with p = 383 is the regime quantitative genetics lives
in, and the field's answer is not variable selection but a relationship matrix
plus REML-estimated shrinkage — one *estimated* hyperparameter rather than a
tuned one, which is how "capacity is a liability" gets solved instead of
avoided. Does it beat the zero-tuning index?

**Design.** `explore/scripts/01_kernel_blup.R` with
`explore/R/methods_kernel.R`. Kernels: all predictors, climate+soil, five
broad blocks, a Matérn-style spatial kernel on centroids, and combinations.
Shrinkage by REML — the exact single-kernel eigendecomposition solver, or
direct optimisation of the profiled likelihood for several kernels. All 144
multi-kernel fits converged.

**Result 1 — transport: a consistent win.** Mean Spearman over 22 held-out
cells, level target: `blup_cs_spatial` **0.386** (21/22 positive),
`blup_cs` 0.360, `blup_all` 0.345, `blup_5k` 0.334, against the domain index's
0.298 and the spatial kernel alone at 0.199 (geography does not transport, as
expected).

Paired per cell the gain is +0.088 but only 14 of 22 cells (p = 0.29). Cells
within a held-out country share districts and predictors, so the country is the
honest unit — and at that level it is unambiguous:

| held-out country | gain, level | gain, prevalence |
|---|---|---|
| Gambia | +0.079 | +0.059 |
| Ghana | +0.041 | +0.009 |
| Malawi | +0.078 | +0.019 |
| Sierra Leone | +0.151 | +0.174 |
| **all 4 positive, both targets** | **+0.087** | **+0.065** |

Eight of eight country × target blocks positive, sign p = 0.008. All four
kernel arms beat the index in the mean on both targets — seven of seven
arm × target comparisons in the same direction.

**Result 2 — in-country: the kernel LOSES.** Paired against the index on
in-fill level, `blup_cs` is **−0.060, winning 5 of 18 cells**; `blup_all`
−0.037 (6/18); `blup_cs_spatial` −0.042 (6/18). The medians say the opposite
(`blup_cs` 0.469 vs the index's 0.394) — this is the AE-01 lesson again, and it
is now the second time in this folder that the median and the paired test have
disagreed in **direction**. Read the paired column.

**Result 3 — the nested test confirms the standing conclusion rather than
challenging it.** Comparing a covariate BLUP to the GAM smoother confounds
"covariates help" with "the kernel is a better smoother". The clean contrast is
identical machinery with covariates added as a second kernel:

| contrast | mean gain | wins | sign p |
|---|---|---|---|
| blup_cs_spatial vs blup_spatial | +0.015 | 9/18 | 1.00 |
| blup_5k vs blup_spatial | +0.026 | 9/18 | 1.00 |
| blup_cs vs blup_spatial | −0.003 | 8/18 | 0.81 |

Exactly chance. **Within a surveyed country, covariates add nothing on top of a
spatial smoother even when the smoother and the covariate model are the same
estimator** — standing conclusion 2 survives a test built specifically to
break it.

**Result 4 — the variance decomposition does not support a domain ranking.**
All 144 fits converged, but with collinear kernels REML returns a *sparse*
solution: one block takes the variance and the rest go to zero, and which one
varies by cell — space dominates 10 of 24 cells, climate 6, embedding 4, soil 2,
agriculture 2, with residual a median 0.52. Reported as instability, not as an
attribution: the median share of every covariate block is 0.000 and that means
"usually not selected", not "contributes nothing".

**Result 5 (added after the synthesis) — the transport win is NOT the kernel.**
The comparison above is against the *full* domain index. The project already
has a better transport arm on the record: the climate+soil index (`index_cs`,
0.368 on the record, reproduced here at 0.374). Head to head on the same 44
transport cells:

| | mean Spearman | vs index_cs | cells | blocks |
|---|---|---|---|---|
| index_cs (the pre-registered candidate) | 0.324 | — | — | — |
| blup_cs_spatial | 0.324 | **−0.0003** | 23/44 | 3/8 (p = 0.73) |
| blup_cs | 0.315 | −0.0096 | 21/44 | 2/8 (p = 0.29) |

**Statistically indistinguishable.** The entire apparent gain was "restrict to
climate and soil instead of using all domains" — which is the project's own
existing pre-registered candidate. The kernel machinery and the REML-estimated
shrinkage add nothing on top of it.

**Verdict (corrected).** *Not an independent candidate.* KB-01 re-derives the
known climate+soil transport result by a completely different estimator —
worth having as independent corroboration that the finding is about the
*predictors*, not about the elastic-net-and-index machinery that produced it —
but it is not a new method to adopt. *Dead end* for in-country use, where it is
worse than the index and adds nothing over geography.

The corollary is the useful part: **two unrelated estimators, plus EB-01's
cross-cell meta-analysis of marginal associations, now independently converge
on climate + soil.** That is stronger evidence for the pre-registered recipe
than any one of them alone.

## NP-01 (2026-09-28) — no untried n≪p estimator beats the index, and XO-01's defect turns out to be costless

**Question.** The project has tried elastic net, the zero-tuning index,
HAL/PCHAL and a SuperLearner. The standard chemometrics and genomics answers to
n ≪ p are untried here. Do any of them help?

**Design.** `explore/scripts/08_np_estimators.R`: partial least squares
(2 components, no tuning), principal components regression (its unsupervised
twin, to separate "low rank helps" from "supervision helps"), supervised PCA
(Bair–Tibshirani, screened on training rows only), MCP, stability selection
over 25 subsamples, CV-tuned ridge on all predictors, and `index_std` — the
domain index with each axis standardised before summing, which is the fix for
the defect XO-01 identified (the index sums *un-standardised* axes, so any axis
with small spread is silently under-weighted).

**Result 1 — nothing beats the index, and the selection methods are
significantly worse.** Block-level, both targets pooled:

| arm | in-fill gain | blocks | transport gain | blocks |
|---|---|---|---|---|
| index_std | −0.003 | 4/6 | −0.012 | 3/8 |
| ridge_cv | −0.066 | 1/6 | **+0.027** | 6/8 (p = 0.29) |
| spca | −0.033 | **0/6** | −0.023 | 3/8 |
| pls2 | −0.080 | 1/6 | −0.037 | 2/8 |
| pcr2 | −0.078 | **0/6** | −0.040 | 1/8 |
| stabsel | −0.133 | **0/6** | −0.115 | 3/8 |
| mcp | −0.213 | **0/6** | −0.103 | 1/8 |

MCP, stability selection, PCR and supervised PCA are harmful in-country with
every block agreeing (p = 0.031 each). **Four more estimator families confirm
"capacity is a liability at n = 14–87"**, and specifically that *selection*
loses: the two arms that pick a subset (MCP, stability selection) are the two
worst. Plain CV-ridge on all 383 predictors is the best of the new arms for
transport (+0.027, 6 of 8 blocks) but does not reach significance.

**Result 2 — XO-01's defect is real but costs nothing.** Standardising the
index's axes before summing changes the result by **−0.003 in-fill (4 of 6
blocks) and −0.012 on transport (3 of 8)**. The under-weighting XO-01
identified is genuine arithmetic, but the axes evidently do not differ enough
in spread for it to matter. That closes an open item: the index does not need
fixing on this account.

**Verdict.** *dead end* for all seven arms. Useful as a negative because it is
broad — the standard n ≪ p toolkit, applied honestly, does not beat a
zero-tuning index at this sample size.

---

# Synthesis (2026-09-28)

Nine probes, all scored on the project's own folds and metrics, with the gate
confirming the harness reproduces `benchmarks_v2_cells.csv` (LOCO exact).
Full table: `explore/out/leads.csv`.

## The one new lead

**MT-01, multi-trait BLUP.** When a survey measures several biomarkers, fit
them jointly through a shared latent factor instead of one cell at a time:
+0.038 over single-trait on identical kernels and folds, **6 of 6 country ×
target blocks positive** (p = 0.031). It carries the signature that
distinguishes a real effect from a lucky one — the gain is concentrated exactly
where single-trait fitting is weakest (weak tercile +0.086, strong +0.018;
Spearman(baseline, gain) = −0.388, p = 0.020), and the large gains are rescues
of cells currently reported as failures (Malawi women_iron prevalence
−0.137 → +0.042; Ghana women_vitA prevalence −0.037 → +0.133).

Not pre-registered, four countries, Sierra Leone absent. Pre-register before
quoting.

## The corroboration, which is worth as much

**Three unrelated routes now converge on climate + soil**, where before there
was one:

1. the existing leave-one-country-out domain ablation (0.368 on the record);
2. **KB-01**, a REML-shrunk kernel BLUP — an estimator sharing no machinery
   with the index — which reaches the same place and is *indistinguishable*
   from the climate+soil index head to head (−0.0003 over 44 transport cells,
   3/8 blocks, p = 0.73);
3. **EB-01**, a cross-cell empirical-Bayes meta-analysis of marginal
   associations, in which climate (37), satellite embedding (35), soil (28) and
   greenness (11) supply 111 of the 173 predictors clearing FDR < 0.05.

The pre-registered recipe is better supported than it was, and now specifically
supported as a claim about the *predictors* rather than about the
elastic-net-and-index machinery that first produced it.

## What was closed

| direction | why it is closed |
|---|---|
| predictor-side batch correction (RUV/SVA) | a *mathematical identity*: rank-normalising within country forces every linear combination to mean ≈ 0 in every country, so there is nothing to project out (verified: country unpredictable from PC1–4, F = 0.0, η² = 0.000) |
| AlphaEarth as a kernel | not a representation problem; the kernel loses to the PC version (0/6 blocks) and adding the embedding to climate+soil hurts (4/22 cells) |
| nutrient-specific crop-basket mechanism | features validated against agronomy, then the *wrong* nutrient's features did as well as the right one's (6/18) — no specificity |
| seasonal phase at the cluster | +0.011, 12 of 24 cells — chance |
| previous growing season (lags 9–15) | +0.004; and the lag profile is **flat**, peaking at lag 31, with weather three years before the draw as associated as the last harvest — a spatial association, not a temporal one |
| the standard n ≪ p toolkit | PLS, PCR, supervised PCA, MCP, stability selection, CV-ridge: none beats the index; the two *selection* methods are the two worst |
| XO-01's un-standardised-axis defect | real arithmetic, costless in practice (−0.003, 4/6 blocks) |

## What the closures say together

Every attempt to extract more from the *predictor side* failed — finer temporal
resolution, finer representation, mechanistic specificity, batch correction,
better estimators. The one thing that worked borrowed strength across
**outcomes**, not predictors. Combined with the nested test (covariates add
nothing on top of a spatial smoother even with identical machinery: 9/18) and
with the flat lag profile, the picture is consistent: at Admin-2 the
predictable signal is a smooth agro-ecological surface, it is close to
saturated by climate + soil, and the remaining headroom is in how the
*outcomes* are pooled and how much survey noise sits under them — not in more
or better covariates.

## Two methodological cautions that generalise

1. **Read paired, not median.** Twice — AE-01 and KB-01 — the median across
   cells and the paired per-cell comparison disagreed in *direction*. With
   18–22 cells of very unequal difficulty the median is not a comparison.
2. **The country is the unit, not the cell.** Cells of one country share
   districts and predictors. The first synthesis reported *zero* candidates
   because it sign-tested over cells; both real positives are unambiguous over
   country blocks. Blocks must pool the two targets — with four countries,
   binom.test(4, 4) = 0.125 can never clear 0.05.

## If this is taken further

- Pre-register MT-01 for the next country, and score it on the weak cells
  specifically, since that is where the mechanism says the gain lives.
- RV-01's by-product is untested and cheap: the between-country share of level
  variance is 0.80 for child iron but **0.14 for women B12**, so the rank-only
  restriction on transport may be over-general. Transporting a *level* for B12
  and vitamin A is a concrete, falsifiable next probe.
- `explore/out/04_moderated_axis_stats.csv` ranks all 329 predictors by
  cross-cell replicated association — usable for annotation prioritisation,
  but read with EB-01's caveat: the satellite embedding supplies 35 of the most
  replicated associations in the set and adds nothing predictively.

---

## LT-01 (2026-09-28) — transporting a LEVEL: my own suggestion was wrong as stated, and the fix is spread calibration

**Question.** RV-01 found the cross-survey level offset is not a global problem
— between-country share of level variance is 0.80 for child iron but 0.14 for
women B12 — and I suggested on that basis that transporting a *level* rather
than a ranking might be feasible for B12. This probe tests it.

**Design.** `explore/scripts/12_level_transport.R`. Leave-one-country-out on
the **raw, un-standardised** level (the within-country standardisation is
exactly what the rank-only protocol does, and is what is being tested).
Predictors still rank-normalised within country. Two arms are **oracles used as
measuring devices, not methods**: `null_true_mean` predicts every district at
the held-out country's true mean (a perfect anchor, no ranking), and the
anchored arms shift the prediction to that mean. MAE is reported relative to
the held-out country's own sd, so it is comparable across biomarkers.

**Result 1 — RV-01's ordering is confirmed, precisely.** Share of squared
transport error that is pure country offset: women_b12 0.159, women_vitA 0.274,
child_vitA 0.339, women_iron 0.419, child_iron 0.511, women_folate 0.665. That
reproduces RV-01's between-country variance shares (0.14 / 0.35 / 0.28 / 0.71 /
0.80 / 0.58) from a completely separate computation.

**Result 2 — but the suggestion was wrong. An unanchored transported level is
unusable for every outcome, B12 included.** `index_cs` MAE is **1.117× the
held-out country's own sd** for B12 (child vitA 1.120, women vitA 1.173, iron
1.93–2.34). Its skill against simply predicting the training countries' mean is
**negative** for every outcome (−0.043 for B12). A prediction whose average
error exceeds one standard deviation is not a usable level.

**Result 3 — the binding constraint was not the offset, it was spread.** With a
perfect anchor the ranking is good (Spearman 0.47 for B12, 0.44 vitamin A) yet
`index_cs_anchored` still scores 0.854 MAE/sd — *worse* than the 0.798 a
constant at the true mean achieves for a normal outcome (√(2/π); the
`null_true_mean` arm lands at 0.762–0.828, confirming the metric behaves).
The index is scaled to the training countries' sd, i.e. as if its ranking were
perfect. A prediction correlating ρ with truth should carry spread ρ × sd.

**Result 4 — shrinking the spread by a nested ρ fixes it, for every outcome.**
ρ estimated by leave-one-*training*-country-out, so the held-out country never
informs it. MAE/sd, anchor + calibrated spread vs the flat anchor:

| outcome | flat anchor | anchor + calibrated ranking | nested ρ |
|---|---|---|---|
| **women_b12** | 0.762 | **0.707** | **0.456** |
| child_vitA | 0.791 | 0.726 | 0.126 |
| women_vitA | 0.798 | 0.747 | 0.150 |
| women_iron | 0.772 | 0.761 | 0.058 |
| child_iron | 0.805 | 0.793 | 0.031 |
| women_folate | 0.828 | 0.828 | 0.050 |

Better in all six. Across the 22 outcome × held-out-country units the mean
reduction is +0.033 MAE/sd; **8 units are exact ties** because ρ shrank to ~0
and the arm correctly collapses to the constant. Of the 14 non-tied units the
calibrated arm wins **11** (sign p = 0.057), the best gains are +0.185 and
+0.184, and the worst loss is −0.050.

**B12 is genuinely the exception, and now there is a number for why.** Its
nested ρ is **0.456**, three to fifteen times every other outcome (0.031–0.150).
That is the same ordering RV-01 gave, arrived at independently.

**Verdict (superseded — see the reconciliation below).**

**Limits.** The anchor is an oracle here; in deployment it comes from a small
survey with its own error, which will eat some of a gain this size. B12 and
folate have three countries, so training is on two. Not pre-registered.

### LT-01 reconciliation against script 78 (LV-02) — the "refinement" was already the production design, and the result reproduces LV-02's null

Checked at the user's request against `scripts/protocol_v2/78_flat_national_comparator.R`
and `docs/findings/LV-02_NZ-01_LEVELS_2026-09-28.md`.

**1. The spread calibration is not new. It is already in the design.** LV-02's
A1 arm is

```
A1_u = expit(logit(anchor) + rho_train * sd_train * z_u)
```

with `rho_train` the mean nested leave-one-training-country-out Spearman,
floored at 0 — the same quantity LT-01 computes and the same place it is
applied. So LT-01's "refined candidate" is a description of what the project
already does, not a change to it. What LT-01 adds is the *diagnostic* for why
that shrinkage is necessary: without it, an anchored ranking is over-dispersed
and MAE punishes over-dispersion hard enough that a real ranking
(Spearman 0.47) scores worse than a flat constant (0.854 vs 0.798). That
explains a design choice; it does not improve it.

**2. My headline number is LV-02's number.** Counted the way LV-02 counts —
ties are a "no", because the question is whether to add the ranking at all —
LT-01 gives **11 of 22 outcome × country units better, 3 worse, 8 tied**.
LV-02 internal, climate+soil: **11 of 22 cells**. The same split. I reported
11 of 14 by dropping the ties; LV-02's accounting is the right one for the
deployment question, and under it my probe *reproduces its null* rather than
refining it. Dropping ties answers a different and narrower question ("when the
method does act, does it help?") and should not have been presented as the
headline.

**3. LV-02 is the better test and it was pre-registered.** Its reading was
fixed before any result ("A1 beats A0 on average AND in a majority of cells"),
mine was not. Its anchor is a realistic draw from a survey fraction
f = 0.05/0.25/1.00, carrying its own sampling error; **mine was an oracle**,
the exact held-out mean. And it reproduces all 14,080 published A1 draws to
5e-14, which mine has no equivalent of.

**4. One genuine tension, which is a scope difference rather than a
contradiction.** LT-01 finds B12's nested ρ = **0.456**, the highest of any
outcome; LV-02 reports that its 8 external ties "are the B12 and folate cells:
their nested training rho is negative or zero". These are different
measurements — LT-01 is Admin-2, the **level** target, the four study countries;
LV-02's external arm is admin-1, **prevalence**, WHO VMNIS deposits with two
training countries. So "B12 carries an unusually transportable ranking" holds
in the internal Admin-2 level setting and **does not replicate** at admin-1 on
prevalence in the external deposits. Recorded as rung- and target-specific, not
as a general property of B12.

**What survives from LT-01.** Result 1, the offset share of squared transport
error by outcome (B12 0.159 → folate 0.665), independently reproducing RV-01's
between-country variance shares — new and unaffected. Result 2, that an
*unanchored* transported level is unusable for every outcome including B12
(MAE 1.117 × the held-out sd, negative skill against the training-country mean)
— stands, and agrees with LV-02's "almost all of the level accuracy comes from
the national figure". Result 3, the over-dispersion diagnostic — stands as an
explanation.

**Verdict (corrected).** *Not a candidate.* LT-01 reproduces LV-02's
pre-registered null on the same 11-of-22 split, using a weaker design with an
oracle anchor, and its proposed "fix" is the production design's existing
behaviour. LV-02's reading stands unchanged: at the shrinkage the training
countries justify, the transported ranking neither helps nor hurts district
levels on average, and its value is in the ordering, not the level. The one
piece worth keeping is the B12 ρ = 0.456 at the Admin-2 level rung, flagged as
not replicating at admin-1 on prevalence.

**Process note for this folder.** LT-01 was written up before checking an
existing script that tested the same question. The probe was cheap and the
error was caught, but the order should be reversed: search `docs/findings/` and
`scripts/protocol_v2/` for the question before running a probe, not after.

---

## Prior-art check on the other probes (2026-09-28)

Prompted by the LT-01 failure, each probe's question was searched against
`docs/findings/`, `docs/` and `scripts/`. Results:

| probe | closest prior art | status |
|---|---|---|
| **MT-01** multi-trait BLUP | **XO-01** (`63_cross_outcome_borrowing.R`) | **distinct.** XO-01 borrows a held-out *country's* other biomarkers as predictors under transport; MT-01 fits a country's outcomes jointly for a held-out *district* that has no biomarker measured at all. Different estimand (C vs A) and different mechanism (observed covariate vs shared spatial component). MT-01's trait correlations reproduce XO-01's — same-nutrient-across-population strong, cross-nutrient weak and sign-unstable — so it builds on XO-01 rather than repeating it. |
| EB-01 screen | `scripts/covariates/16_bivariate_fdr.R` | **partly prior art.** Script 16 runs one bivariate test per predictor × country × outcome with BH FDR — per cell. EB-01 pools *across* the 24 cells with a moderated random-effects statistic, which is the new part; the FDR-screen idea is not. |
| KB-01 kernel/REML shrinkage | Fay-Herriot / EBLUP in the main pipeline; `WSB_ESTIMATOR_TOURNAMENT` | **partly prior art.** Area-level random-effects shrinkage with REML-estimated variance is already the pipeline's primary SAE estimator; KB-01's novelty was the kernel form over covariates, and its result was in any case reduced to corroboration. |
| AE-01 AlphaEarth | script 42, `addon_alphaearth_*` | prior art for the *domain* arm (±0.012 on the record); the kernel representation was new and failed. |
| RV-01 RUV | `WS3_measurement_harmonization`, `AS-01` assay lineage | question previously approached through assay/adjustment harmonisation, not through subspace removal. New, and closed structurally. |
| MX-01 crop mechanism | `WSF_NUTRITION_PROXIMAL`, `build_mapspam_admin2.R` (5 group shares) | new at the 42-crop / composition-table level; the 5-group version is prior art and was already weak. |
| TM-01 phase, TM-01b lags | `FW-01`, `FP-02`, `41_climate_timematched.R`, `CN-01` | fieldwork windows and time-matched climate are prior art (+0.02); phase and the 0–35 month lag stack were new. |
| NP-01 estimators | `SL-06`, `HP-01/02/03`, `WSB_ESTIMATOR_TOURNAMENT` | the tournament covered SuperLearner/HAL/enet; PLS, PCR, supervised PCA, MCP and stability selection were new. |
| **LT-01 level transport** | **LV-02** (`78_flat_national_comparator.R`) | **superseded** — see the reconciliation above. |

**Standing rule for this folder, from LT-01.** Search `docs/findings/` and
`scripts/protocol_v2/` for the question *before* running a probe. One grep would
have caught LT-01 before it was written up.

---

# Seasonality and collection timing (2026-09-28) — probes CG-01, MA-01, PT-01, SO-01

## The exploration that gated all four

Individual collection timing exists for two countries: Malawi
(`date_interview`, `time_blood_draw`, `fast` in `clean_malawi_mn_data.RDS`) and
Gambia (`gw_cIntDate`, `gw_wIntDate`, and **`gw_cVASDate`** — the vitamin A
supplementation date). Ghana is month-resolution only; Sierra Leone's dates are
reconstructed from date of birth plus age.

Three structural facts, and they decide most of what follows:

1. **A survey team visits a cluster in about one day** — median within-cluster
   date spread **1 day**, maximum 6.
2. **Date is therefore almost the same variable as place**:
   R²(date ~ cluster) = **0.999** Malawi, **0.989** Gambia;
   R²(date ~ district) = **0.908** Malawi.
3. **The raw date slopes are small but not zero**: Malawi log RBP
   −0.0009/day (t = −2.3), AGP +0.0013/day (t = 2.0) — −6% and +9% across the
   70-day window.

So collection date cannot serve as an individual predictor net of geography,
and cannot be adjusted out of district estimates, because it is collinear with
them. Fasting is unusable (65 of 3,099 fasted). Gambia's days-since-VAS is
mechanistically right — VAS within 60 days raises log RBP by 0.038, about 3.9%
— but n = 276 gives t = 1.0.

## CG-01 — the fieldwork calendar is NOT a smooth spatial surface, and the standing conclusion survives

**The hypothesis (mine, and it was wrong).** Fieldwork teams move through space
contiguously, so the survey calendar should itself be a smooth spatial surface.
If so, the project's spatial smoother might be fitting the *survey schedule*
rather than geography — which would reinterpret standing conclusion 2.

**Result — refuted on its own premise.** The calendar is not spatially smooth.
R²(fieldwork date ~ thin-plate spline in lon/lat), per country:

| country | R² |
|---|---|
| Sierra Leone | 0.572 |
| Ghana | **−0.010** |
| Malawi | **−0.023** |
| Gambia | **−0.054** |

Negative adjusted R² means the spline fits worse than a constant. In three of
four countries a spatial smoother **cannot** represent the fieldwork calendar at
all: teams are evidently assigned in a way that scrambles date across geography
rather than sweeping through it. The premise fails, so the threat does not
arise.

Consistent with that, the calendar explains very little district outcome
variance — R²(outcome ~ date) = 0.011 Ghana, 0.025 Malawi, 0.030 Gambia (0.112
Sierra Leone, on 14 districts) — and residualising the outcome on fieldwork date
inside the training fold moves the spatial smoother by **+0.005** (level) and
**−0.008** (prevalence). Date alone predicts *negatively* (−0.108 level, −0.146
prevalence).

**This is a reassuring negative.** The district targets the whole project
predicts are not a survey calendar in disguise, and they cannot be, because the
calendar has no spatial structure to be confused with geography.

**A third median-versus-paired reversal, noted in passing.** On these folds the
index beats the spatial smoother in **13 of 18** cells paired, while the record's
median ordering has spatial ahead (0.457 vs 0.394). That is the same reversal
AE-01 and KB-01 found. It does not disturb standing conclusion 2, which rests on
the *nested* test (covariates on top of the smoother), and KB-01 confirmed that
at 9/18 — chance. The two arms are about equally good and neither adds to the
other.

## MA-01 — removing collection artefacts does not improve district reliability

Split-half reliability of district means (200 random within-district splits),
raw against artefact-adjusted (artefacts removed *within cluster*, so the
adjustment cannot absorb between-district signal):

| country | marker | n | artefacts | raw | adjusted | gain |
|---|---|---|---|---|---|---|
| Malawi | rbp | 2,886 | time of draw + fasting | 0.350 | 0.349 | −0.001 |
| Malawi | vitb12 | 810 | time of draw + fasting | 0.696 | 0.694 | −0.001 |
| Malawi | zn_gdl | 2,885 | time of draw + fasting | 0.703 | 0.700 | −0.004 |
| Gambia | gw_cRBP | 1,020 | days since VAS | 0.453 | 0.457 | +0.004 |

Mean gain **−0.0006**, 1 of 4 markers improved. The artefacts are real at the
individual level and irrelevant at the district level, which is what the
exploration predicted (adjusting Malawi district means for time of draw left the
ranking at Spearman 0.999).

**The by-product is more useful than the result.** Those raw split-half numbers
are a *measurement* explanation for a pattern the project has been trying to
explain with covariates: **vitamin A has far lower district reliability
(RBP 0.350 Malawi, 0.453 Gambia) than B12 (0.696) or zinc (0.703)**. Vitamin A
cells are the weak ones and B12 the strongest cell on the record — and this says
that ordering is set at the blood draw, before any model sees the data. A
predictor cannot recover district signal the biomarker does not reliably carry.

## PT-01 — collection timing does not help individual prediction, except zinc

Individual ridge on the deficiency indicator, cluster-blocked 5-fold CV, timing
block against no covariates:

| country | marker | no covariates | + timing | Brier skill |
|---|---|---|---|---|
| Malawi | rbp | 0.457 | 0.453 | −0.003 |
| Malawi | vitb12 | 0.439 | 0.440 | −0.005 |
| **Malawi** | **zn_gdl** | 0.453 | **0.552** | **+0.008** |
| Gambia | gw_cRBP | 0.448 | 0.441 | −0.004 |

Only zinc moves, AUC 0.453 to **0.552**, and it is the one biomarker ZN-02
already identified as carrying a collection artefact (afternoon draw −3.5%,
+0.7%/day). So this reproduces ZN-02 at the individual level through a different
route. It is a measurement artefact, not nutrition signal, and MA-01 shows it
does not propagate to district estimates.

*Reading note:* the no-covariate arm predicts each fold's own training
prevalence rather than one constant, so its AUC is near 0.45 rather than exactly
0.50. Compare the arms to each other, not to 0.5.

## SO-01 — the cross-country level offset against season: suggestive, one country, not a test

Seasonal position of each survey, from the Earth Engine monthly NDVI stack:

| country | fieldwork month | months since NDVI peak | NDVI at fieldwork (share of annual range) | mean offset |
|---|---|---|---|---|
| Sierra Leone | Nov | 1 | 0.88 | −0.159 |
| Gambia | Mar | 6 | **0.05** | **+0.448** |
| Ghana | May | 7 | 0.67 | −0.038 |
| Malawi | Jan | 10 | 0.54 | −0.102 |

Spearman(months since peak, offset) = +0.20 over four countries — nothing. But
the *greenness* ordering is suggestive: the one country surveyed at its seasonal
trough (Gambia, at 5% of its annual NDVI range) carries by far the largest
deficiency offset, and the correlation of NDVI-at-fieldwork with the offset is
about −0.6.

**That is one country doing all the work**, and Gambia differs from the others
in many ways besides season. Recorded as consistent with a seasonal contribution
to the offset, and as nothing more. Four countries cannot test this; it needs
either a country surveyed twice in different seasons, or multi-round surveys the
design does not have.

## Verdicts

*Dead end* for CG-01 (premise refuted — and reassuringly so), MA-01 and PT-01.
*Needs data* for SO-01. **Keep** the district split-half reliabilities from
MA-01: vitamin A 0.35–0.45 against B12 and zinc at about 0.70 is a
measurement-side explanation for which cells work, and it is cheap to extend to
the remaining biomarkers and countries.

---

## RL-01 (2026-09-28) — district reliability for every cell, and what dichotomising costs

**Extension of MA-01's by-product** to all 27 country × outcome cells the
pipeline defines, using the project's own config, loader and population masks so
the populations match the targets exactly. Respondents are split at random
within district (300 splits), district means correlated, Spearman-Brown
corrected to full sample.

### Reliability by nutrient (continuous level, Spearman)

| nutrient | reliability | cells |
|---|---|---|
| selenium | **0.939** | 2 |
| folate | 0.720 | 3 |
| B12 | 0.704 | 3 |
| iodine | 0.682 | 1 |
| zinc | 0.681 | 2 |
| iron | 0.532 | 8 |
| **vitamin A** | **0.518** | 8 |

**The project's two headline nutrients are its two least reliably measured**, and
they are the only two carried by all four countries. By country: Gambia 0.724,
Malawi 0.655, Ghana 0.594, Sierra Leone 0.489.

### Reliability predicts model accuracy

Against the record's in-fill level accuracy for the domain index:
**Spearman = 0.521, p = 0.028 over 18 cells.** Roughly half the cross-cell
variation in how well the model does is explained by how reliably the biomarker
measures districts at all — before any covariate is involved. Malawi women_b12
(reliability 0.843, accuracy 0.694) and Gambia women_vitA (0.805, 0.697) sit at
the top; Malawi women_vitA (0.301, 0.302) and child_vitA (0.377, 0.230) at the
bottom.

**Zinc is the exception that confirms the reading**: reliably measured (0.72
women, 0.64 children) and unpredictable (index −0.007 and −0.105). That
independently reproduces ZN-01 — "a reliable target the proxies do not touch" —
from a different direction.

### Prior art, and the discrepancy that turned into the finding

Applying this folder's standing rule: **WS1a already computes empirical
split-half reliability** (`R/reliability_empirical.R`,
`results/tables/reliability_empirical.csv`, 200 splits, Spearman-Brown, two
schemes). My numbers correlate with it at only 0.22, with a mean absolute
difference of 0.281, and WS1a reports **0.000** for Malawi and Sierra Leone B12
where I get 0.84 and 0.80.

That is not an error in either. **WS1a computes the reliability of district
PREVALENCE from the binary outcome with Pearson; the table above is the
continuous LEVEL with Spearman.** Different quantities. Recomputing both in one
framework prices the difference:

**Dichotomising costs +0.203 reliability on average, and the continuous level is
more reliable in 25 of 27 cells.** The cost is concentrated where the cutoff is
extreme — exactly as theory predicts, since a binary at a rare threshold has
little variance and large sampling noise:

| cell | prevalence | level reliability | prevalence reliability | cost |
|---|---|---|---|---|
| Sierra Leone women_b12 | 0.005 | 0.797 | −0.350 | **1.147** |
| Sierra Leone women_vitA | 0.023 | 0.607 | −0.027 | 0.634 |
| Malawi women_vitA | 0.019 | 0.282 | −0.127 | 0.409 |
| Ghana women_vitA | 0.013 | 0.634 | 0.234 | 0.401 |
| Malawi women_b12 | 0.107 | 0.840 | 0.493 | 0.347 |
| … | | | | |
| Malawi child_selenium | 0.855 | 0.935 | 0.917 | 0.017 |
| Sierra Leone child_vitA | 0.183 | 0.122 | 0.198 | −0.076 |
| Sierra Leone women_iron | — | −0.006 | 0.484 | −0.490 |

**And the reliability gap tracks the accuracy gap**: Spearman(cost of
dichotomising, level-minus-prevalence accuracy) = **0.507, p = 0.034** over 18
cells, with a mean accuracy gap of +0.098.

**That explains a pattern the project has carried on the record without a
mechanism.** The level target beats prevalence everywhere — in-fill 0.394 vs
0.298, transport 0.287 vs 0.192 — and this says why: dichotomising at a clinical
cutoff destroys district signal, most severely for rare deficiencies, and the
cells where it destroys most are the cells where the level target wins by most.
It is a measurement explanation for a modelling observation.

### Caveats

- **The cluster-split column is not per-cell interpretable.** Only **13
  districts** per country have two or more clusters (Gambia 57%, Ghana 83%,
  Malawi 85% single-cluster), so that estimate rests on 13 units and swings
  wildly (−0.65 to +0.82). In aggregate it suggests the respondent split is
  optimistic by about +0.23, consistent with CE-01, but no single cell's value
  should be read. Sierra Leone is the only country whose design supports it
  (0% single-cluster) and there it is low: 0.174 mean.
- Two Sierra Leone binaries are coded 1/2 rather than 0/1, so their `prev_mean`
  is not a prevalence. The split-half correlation is invariant to that linear
  recoding, so their reliability is unaffected.
- Reliability here is unweighted; the pipeline's targets are survey-weighted.

### Verdict

*Candidate, and cheap to act on.* Two concrete uses. **(1)** Judge a cell's model
accuracy against its own reliability rather than against 1.0 — Malawi women_vitA
at reliability 0.30 cannot be predicted at 0.6 by anything, and reporting it as
a failure of the covariates misattributes the cause. **(2)** Prefer the
continuous level target wherever a decision permits it, and treat prevalence at
a rare cutoff as the expensive choice it is: Sierra Leone women B12 at 0.5%
prevalence has no usable district signal on the binary and good signal on the
level.

---

## GN-01 (2026-09-28) — measured grain and soil chemistry (GeoNutrition Malawi): the mechanism is real, the increment is zero

**Why this source.** MX-01 built nutrient-specific features by multiplying crop
*production* by a generic food-composition table — assuming maize everywhere has
the same zinc — and failed with no nutrient specificity. I diagnosed that as the
composition proxy being wrong. GeoNutrition tests the diagnosis directly: it
*measures* the mineral concentration of grain grown at 1,812 georeferenced
Malawi sites plus the soil chemistry beneath. It also fills a gap this session
verified — the 575-column store contains **no environmental selenium at all**
(iSDA carries Al, Ca, CEC, Fe, Mg, P, K, S, C, Zn; no Se) and no environmental
iodine, only two fortification-*programme* fields.

**Source.** Gashu et al. 2022, *Scientific Data*; figshare
`10.6084/m9.figshare.15911973`, **CC BY 4.0**. Malawi national sampling
April–June 2018, 1,900 sites of which **820 were drawn from the 2015/16 Malawi
DHS frame** — the same frame as this project's biomarker survey.

**Build.** 1,809 of 1,812 sites joined to GADM Admin-2; **86 of Malawi's 87
target districts** covered; 28 variables carried (9 grain, 19 soil) as district
medians of log concentration. Malawi selenium and iodine district targets were
built here because they are configured in `R/config.R` but not yet in
`targets_v2.csv`. Caveat: grain is the 2018 harvest against a Dec 2015 – Feb
2016 biomarker survey, so grain concentration is a proxy for the district's
typical grain, not for what was eaten.

**Result — a clean null.** Adding the block to the index, Malawi in-fill, level:
mean gain **−0.005, better in 5 of 11 cells, sign p = 1.00**. The cells where
mechanism predicts the most are the ones that get *worse*:

| outcome | index | index + GeoNutrition | GeoNutrition alone | gain |
|---|---|---|---|---|
| women_iodine | 0.348 | 0.385 | 0.368 | **+0.037** |
| women_folate | 0.416 | 0.426 | 0.290 | +0.010 |
| child_iron | 0.362 | 0.368 | 0.165 | +0.006 |
| **child_selenium** | **0.469** | **0.452** | 0.276 | **−0.018** |
| women_zinc | −0.021 | −0.040 | −0.192 | −0.019 |
| **women_selenium** | **0.430** | **0.389** | 0.130 | **−0.042** |

**The largest gain is women_iodine — the one nutrient GeoNutrition does not
measure at all** (no iodine in either the grain or the soil panel). That is the
signature of noise, and it is the same pattern MX-01's mismatch control showed.

**But the mechanism is real, and that is the interesting part.** Marginal
district-level associations with selenium status (y_level is *negated* log serum
Se, so negative = more selenium goes with better status — the correct direction):

| variable | child_selenium | women_selenium |
|---|---|---|
| **grain Se** | **−0.351** (n=81) | **−0.362** |
| soil Se, total | −0.112 | −0.085 |
| soil Se, adsorbed | −0.127 | −0.106 |
| soil Se, organic | −0.029 | −0.040 |
| soil Se, soluble | −0.035 | −0.031 |
| grain S | −0.233 | −0.195 |
| soil pH | −0.218 | −0.134 |

So the soil→grain→serum chain the GeoNutrition project is built on **is visible
in this project's own district data** — but it is the **grain** step that carries
it (|r| ≈ 0.35), not the soil step (|r| 0.03–0.13).

**Why the increment is nevertheless zero.** The index alone already reaches
**0.469** on child selenium — better than grain Se's marginal |r| of 0.35. The
remotely sensed store evidently already captures as much of the soil-to-grain
selenium gradient as measured grain chemistry supplies, presumably through the
climate, terrain and soil variables that drive both. Ground-truth grain
chemistry is redundant with it rather than additional to it.

*(An in-sample R² of 0.69 for the store predicting soil Se is not evidence of
this — n = 84 against 383 columns will fit anything — and is excluded from the
reading.)*

**What it corrects.** My MX-01 diagnosis was wrong. The composition proxy was
not the problem: supplying measured, spatially varying composition changes
nothing, on the outcome where the mechanism is most direct *and* where
measurement is most reliable (RL-01: Malawi selenium district reliability 0.94,
the highest in the study). A measurement-noise explanation is therefore
excluded. The district-level signal is not dietary mineral supply.

**What it changes for the data-source shortlist.** A global *soil* selenium map
(Jones et al. 2017) is now a much weaker prospect than it looked: soil Se is the
weak link here (|r| 0.09–0.13) and grain Se the strong one, and grain-composition
surveys do not exist outside Malawi and Ethiopia.

**Verdict.** *Dead end* as a predictor block, on a well-powered test. **Keep the
data** — it is CC BY 4.0, already built to Admin-2 in
`explore/out/17_geonutrition_admin2.csv`, and the grain-Se association is a
publishable external corroboration of the soil-to-serum chain in this project's
own units, independent of the GeoNutrition team's own analysis.

---

## Methods from the Bayesian SAE preprint that this project has NOT applied

Read: *Mapping Subnational Vulnerability to Inadequate Micronutrient Intake
using a Bayesian Small Area Estimation Framework* (arXiv 2604.14971). Rwanda
(validation), Senegal, Nigeria; ADM2; outcome is household inadequate apparent
intake from HCES, not biomarkers.

Most of its machinery the project already has: BYM2 with PC priors (SL→BYM2 and
the surveyPrev cluster model, DS-01), beta-binomial cluster models, design-based
direct estimates with Taylor linearisation, population-weighted aggregation to
ADM1, posterior credible intervals (CP-01). **Three things it does that this
project does not:**

**1. Joint variance–mean smoothing — the strongest candidate.** The paper models
the sampling variance rather than plugging it in:

```
log(V_ℓ) = γ₀ + γ₁ log(p_ℓ(1−p_ℓ)) + γ₂ log(n_ℓ) + τ_ℓ,   τ_ℓ ~ N(0, σ²_τ)
```

This project's Fay-Herriot (`R/benchmark_models.R:229–270`) instead treats
sampling variances as **known**, as `p(1−p)/n` with the **design effect fixed at
1.5** — and the project's own evidence says that constant is wrong in both
directions (WS1: the reconciling value has median 0.969; protocol v2: deff is
~2.4). Modelling the variance removes the constant entirely and lets the data
set it per area. Directly relevant to NZ-01's finding that about half of
prevalence error is survey noise.

**2. Phantom clusters for single-cluster areas.** The paper augments ADM2 areas
with one cluster using synthetic ADM1-level observations so a variance can be
estimated at all. This project's FH does the opposite: `sv <- pmax(sv, 1e-8)`
**floors** degenerate one-cluster areas, which gives them a near-zero sampling
variance and therefore near-maximal weight on their own noisy direct estimate.
Given that **85% of Malawi and 83% of Ghana districts are single-cluster**, this
is the project's single most exposed design assumption, and the paper offers a
principled alternative to a floor.

**3. Mean Interval Score (Winkler) and CV reliability thresholds.** The paper
scores intervals with a proper scoring rule that combines width and calibration,
and flags estimates by coefficient of variation (<16.6% unrestricted,
16.6–33.3% caution, >33.3% unreliable). CP-01 currently reports coverage and
width separately, so a narrow-but-miscalibrated arm and a wide-but-honest one
are not directly comparable. MIS makes them comparable in one number, and the
CV bands are a ready-made, externally recognised way to grey out cells the
dashboard should not present.

**Recommended order.** (2) is the highest value and the cheapest to test —
replace the `1e-8` floor with either a modelled variance or the paper's phantom-
cluster augmentation and re-run the FH arm. (1) is the principled version of the
same fix. (3) is presentational but cheap and would improve the dashboard's
honesty about which cells to show.

---

## VF-01 (2026-09-28) — the Fay-Herriot variance assumption is self-correcting, so the preprint's fix does not bite

**The concern, and it was mine.** The Bayesian SAE preprint models the sampling
variance and augments single-cluster areas; this project's FH
(`R/benchmark_models.R:229–270`) instead fixes the design effect at **1.5** and
**floors** single-cluster areas at `1e-8`. Both assumptions are measurably wrong
on this data — the district design effect already in `targets_v2` as
`deff_binary` has **median 2.57** (IQR 1.56–3.46, max 6.1), and **85% of Malawi,
83% of Ghana and 57% of Gambia district-rows are single-cluster**. I predicted
that understating the sampling variance would make FH trust noisy districts too
much when fitting β, and that fixing it would move numbers.

**Design.** FH implemented directly so that *only* the variance rule differs
between arms: production (`deff` 1.5, floored), each district's own measured
design effect, the preprint's joint variance–mean smoothing
(`log v = γ₀ + γ₁log(p(1−p)) + γ₂log(n)` fitted on training districts), and the
preprint's phantom-cluster idea (a single-cluster district borrows the pooled
within-Admin1 variance of the multi-cluster districts rather than being
floored). In-fill CV on the prevalence target, 18 cells.

**Result — all four are identical.** Median Spearman **0.227** for every arm;
median MAE 10.64–10.65 points. Paired against the production rule:

| arm | mean gain | cells | sign p |
|---|---|---|---|
| measured design effect | −0.0003 | 8/18 | 0.82 |
| joint variance–mean smoothing | −0.0001 | 8/18 | 0.82 |
| phantom clusters | −0.0003 | 8/18 | 0.82 |

**Why — and this is the finding.** Not because the sampling variances are small:
γ = A/(A+D) runs from 0.00 to 0.78, so shrinkage is substantial. It is because
**A and D are estimated jointly from the same data**. The moment estimator is
A = Σw(r²−D)/Σw, so raising D lowers A by a compensating amount and the *total*
A+D — which is what the weights 1/(A+D) depend on — barely moves. Mean γ under
the production rule against the measured one, per cell: 0.579/0.574,
0.674/0.701, 0.469/0.496, 0.443/0.434. The variance-component split is weakly
identified; the total is not. **Misspecifying the design effect is absorbed by
the between-area variance, and prediction is untouched.**

**What that means for the production code.** The `deff = 1.5` constant and the
`1e-8` floor are **harmless for prediction accuracy**, which is what the
benchmark measures — so this is not a defect to fix, and my recommendation to
fix it was wrong. The split does still matter for anything that *reports* γ as
"how much this district's own survey is trusted", because with D understated by
about 1.7× the reported γ is correspondingly too high. That affects published
EBLUP estimates for surveyed districts, not cross-validated accuracy.

**Two things worth noticing in passing.**

*FH is not competitive on ranking but is better on level.* Against the domain
index it is −0.066 mean Spearman, better in only **3 of 18** cells — but its MAE
is **10.64 points against the index's 13.07**. That is the project's
level-versus-ranking split showing up inside a single estimator.

*In about a quarter of cells FH estimates A = 0 exactly* — Gambia child_iron,
Malawi child_iron, Sierra Leone child_vitA, child_iron, women_vitA and
women_b12. A = 0 means γ = 0: FH concludes the observed district differences are
**entirely sampling noise** and collapses to the synthetic estimate. That is an
independent, model-based version of the same verdict RL-01 reaches by split-half
reliability, and the two agree on the weakest cell in the study (Sierra Leone
child vitamin A: reliability 0.117, A = 0).

**Verdict.** *Dead end*, and specifically a **withdrawal of my own
recommendation**. The preprint's variance machinery is well motivated in its own
setting — where it is used to report uncertainty for surveyed areas — but it
cannot improve out-of-sample prediction here, for a structural reason that
applies to any Fay-Herriot fitted this way.

### Also closed: the `r_share` scale concern

I raised twice the possibility that `r_share = pearson_r / r_max` divides a
level-target correlation by a prevalence-scale ceiling, since
`admin2_reliability()` builds the ceiling from `svy_prev` and `n_svy` only.
Checked all three call sites — `R/area_level_comparison.R:387`,
`R/cluster_mbg.R:293`, `R/corrected/p12_distributional.R:296`. All three are
prevalence-scale throughout (`Y <- train_df$svy_prev`, `mae_pp`/`rmse_pp` in
percentage points, `ind$def`). **The scales match; there is no defect.**
Recorded because I raised it.

---

## RL-01 robustness — the reliabilities survive survey weighting

RL-01's split-half reliabilities are unweighted; the pipeline's targets are
survey-weighted. Recomputed with weighted district means (same 300 splits):

- mean |difference| **0.040**
- Spearman between the two orderings **0.988** over 27 cells
- the headline is **unchanged**: reliability predicts index in-fill accuracy at
  **0.521, p = 0.028** under both weightings, to three decimals

Only Sierra Leone moves materially (child_vitA 0.121 → −0.306; women_iron
−0.001 → −0.208), and both were already the weakest cells in the table.
Weighting makes the weak cells look worse, not better. The nutrient ordering —
selenium, folate, B12 high; iron and vitamin A low — is unaffected.

## CV-01 — the standard CV bands, applied to the published district estimates

The preprint gates estimates by coefficient of variation on the conventional
survey-statistics bands: **< 16.6%** publish unrestricted, **16.6–33.3%**
publish with caution, **> 33.3%** unreliable. This project reports interval
width and coverage but has no such gate. Applied to the **direct** survey
district prevalences using the measured effective n already in `targets_v2`:

| | share of 1,350 district × outcome estimates |
|---|---|
| unrestricted (CV < 16.6%) | **5%** |
| caution (16.6–33.3%) | 9% |
| **unreliable (CV > 33.3%)** | **86%** |

Median CV by country: Ghana **100%**, Malawi **96%**, Sierra Leone **77%**,
Gambia **38%**. Nine cells have **100%** of their districts in the unreliable
band, including all three vitamin A cells outside Gambia.

**Read this correctly, in two ways.**

*First, it is about the DIRECT estimates, not the model.* These are the survey's
own district figures. The modelled estimates should be better — that is what
small-area estimation is for. So the honest reading is not "the project's
outputs are unreliable" but **a quantification of why the modelling is
necessary at all**: 86% of the raw district figures are too imprecise to publish
by a standard external criterion.

*Second, CV is scale-dependent and structurally penalises rare outcomes.*
CV = SE/p, so a low-prevalence outcome is penalised however precisely it is
measured. The cells that pass are simply the high-prevalence ones — Sierra Leone
folate (79% prevalent, 93% of districts unrestricted), Malawi zinc (~55–60%),
Gambia iron. This is the *same* structural fact RL-01 found from the other
direction: rare deficiencies lose most from dichotomisation, and they are also
the ones whose prevalence CV is worst. One cause, two symptoms — a binary at a
rare threshold carries little information.

**Verdict.** *Candidate, presentational.* The CV bands are a cheap, externally
recognised gate the dashboard could apply to any district figure it shows raw,
and they make the case for modelling concrete. They should **not** be applied to
the modelled estimates without recomputing the CV from the model's own posterior
or conformal interval, and they should always be reported next to prevalence,
because a "good" CV here mostly means "common deficiency".

---

## SL-07 (2026-09-28) — a SuperLearner of many *simple* models: tested, and it loses

**The question.** The project has tested a SuperLearner over its protocol *arms*
(SL-06) and over tuned learners (SL-01/02/03), and the zero-tuning index beat
all of them. The other library design was untested: very many very *simple*
candidates — one variable each, pairs, single principal components. That is a
different bias/variance trade from a handful of flexible learners.

**Prediction, recorded before running.** A non-negative convex combination of
univariate least-squares fits is algebraically close to a shrunken ridge on the
same variables, and WS-01 already found the index *is* a max-shrinkage ridge
that nothing beats in-country. So this should land near the index rather than
above it.

**Design.** `explore/scripts/20_simple_library_sl.R`. A proper SuperLearner —
inner 5-fold CV over the training rows for honest candidate predictions, NNLS
meta-weights, candidates refitted on the full training fold. Libraries: every
domain axis alone (~66 learners); plus all pairs among the top-10 screened;
plus single principal components; plus an intercept-only candidate. A rank
meta-learner variant weights candidates by inner-CV Spearman instead of NNLS.

**Result — worse than the index, significantly.** Blocks are country × target:

| estimand | arm | mean gain | cells | blocks | block p |
|---|---|---|---|---|---|
| in-fill | sl_uni | −0.104 | 4/36 | **0/6** | 0.031 |
| in-fill | sl_uni_pair | −0.086 | 5/36 | **0/6** | 0.031 |
| in-fill | sl_pc | −0.086 | 6/36 | **0/6** | 0.031 |
| in-fill | sl_rank | −0.060 | 5/36 | **0/6** | 0.031 |
| transport | sl_uni | −0.072 | 18/44 | 1/8 | 0.070 |
| transport | sl_pc | −0.059 | 19/44 | 1/8 | 0.070 |
| **transport** | **sl_rank** | **+0.009** | 18/44 | 4/8 | 1.000 |

Every block agrees the simple-library SuperLearner is worse in-country. On
transport it ties at best, and only with the rank meta-learner.

**And it reproduces SL-02 by a different route.** SL-02 found that a
rank-aligned meta-learner "recovers the index's performance but does not exceed
it". Here, on a completely different library — hundreds of simple OLS fits
rather than the protocol arms — the same thing happens: `sl_rank` is the best
variant, lands at +0.009 on transport and −0.060 in-country, and does not beat
the index. Two unrelated libraries, same destination.

**Why, and the theory says so in advance.** The SuperLearner oracle inequality
(van der Laan & Dudoit 2003; van der Vaart, Dudoit & van der Laan 2006) bounds
the CV-selector's risk by the oracle's plus a term of order **(1 + log K)/n**,
with the library permitted to grow polynomially in n. That is an *asymptotic*
licence. At n = 14–87 the penalty is not negligible: with K ≈ 120 and n = 30,
log(120)/30 ≈ 0.16 before the leading constant. The promise that "adding
candidates is nearly free" is precisely the promise that fails at this sample
size — and it fails in the measured direction and roughly the measured
magnitude.

**Does PC-HAL supersede this? No — they are different objects.** HAL is a single
penalised regression over a large indicator basis with an L1 constraint on total
variation, fitted jointly, with a rate guarantee (n^(−1/3) under càdlàg and
bounded variation). A SuperLearner of simple models is a cross-validated convex
combination of separately fitted models, with an oracle-inequality guarantee.
Neither contains the other. Empirically they arrive at the same place: HP-01,
HP-02 and HP-03 found hapc/PCHAL ties the index, and SL-07 finds the simple
library ties it at best. Same destination, different routes — which is itself
the informative part.

**Verdict.** *Dead end*, well powered (0 of 6 blocks in-country). The estimator
space is now mapped thoroughly enough to stop: index, elastic net, SuperLearner
over arms, SuperLearner over simple models, HAL, PCHAL, kernel BLUP/REML, PLS,
PCR, supervised PCA, MCP, stability selection, CV-ridge, BYM2, Fay-Herriot and
model-based geostatistics all land at or below a zero-tuning index. **The
binding constraint is not the estimator.**

### What the literature review turned up that IS untried

Three things, in descending order of value.

**1. Data-adaptive target parameters (van der Laan, Hubbard, Pfeiffer).** This
addresses the project's *actual* binding constraint rather than its estimator.
Every finding in this folder and much of the main record is post hoc on the same
four countries, which is why nothing can be quoted without a pre-registration
caveat. The method partitions the sample into a parameter-generating split and
an estimation split: the target parameter is *defined as* whatever the selection
algorithm picks on the first, and inference is done on the second, averaged over
V splits. That makes "the best domain subset chosen by looking at the data" a
legitimate estimand with valid inference, instead of a caveat. With four
countries it is thin, but it is the right framework for exactly the problem the
project keeps hitting.

**2. Moderated variance for semiparametric estimators (Hejazi, Boileau, van der
Laan & Hubbard 2023).** Generalises limma's empirical-Bayes moderation from
t-statistics to the *variance estimators* of asymptotically linear
semiparametric estimators, aimed explicitly at "modest sample sizes" in
high-dimensional biology. EB-01 in this folder did an ad hoc version of exactly
this (moderating per-predictor statistics across 24 cells) and found the
estimator null but the screen useful; this is the rigorous form, and it would
give the screen valid inference rather than a nominal FDR.

**3. Random-projection ensembles (Cannings & Samworth 2017, JRSS-B).** The
closest formal match to the "ensemble of simple models" intuition: apply a base
learner to many random low-dimensional projections, keep the best within
disjoint groups, aggregate with a data-driven vote threshold. Its error bound
does not depend on the ambient dimension under a low-dimensional structure
assumption, and there is an R package (`RPEnsemble`). Given SL-07 and KB-01 both
land at the index, the prior on this helping is low — but it is the one member
of the family with a dimension-free guarantee, and it is cheap.

---

# The binding constraint (2026-09-28)

## What it is not

Sixteen estimator families now land at or below a zero-tuning index: index,
elastic net, SuperLearner over protocol arms (SL-06), SuperLearner over hundreds
of simple models (SL-07), HAL, PCHAL (HP-01/02/03), kernel BLUP with REML
shrinkage (KB-01), PLS, PCR, supervised PCA, MCP, stability selection, CV-ridge
(NP-01), BYM2, Fay-Herriot (VF-01) and model-based geostatistics (MB-01). It is
also not predictor representation (AE-01), not temporal resolution (TM-01a/b),
not dietary mechanism (MX-01, and GN-01 with *measured* grain chemistry), not a
predictor-side batch effect (RV-01, a mathematical no-op), not the Fay-Herriot
variance assumption (VF-01, self-correcting), not the index's weight estimation
(EB-01), and not the un-standardised-axis defect (NP-01, costless).

## The gap, quantified

With `r_max = sqrt(split-half reliability)` as the most any predictor could
reach, over the 18 scorable in-fill cells (level target):

- median attainable **r_max = 0.79**
- median achieved (domain index) **0.39**
- median **r_share = 0.54**, median headroom **0.353**
- 2 of 18 cells are above 80% of attainable; **8 of 18 are below 50%**

So **measurement error is not the whole story.** Measurement explains *which*
cells do well — RL-01 found reliability predicts accuracy at ρ = 0.52,
p = 0.028 — but it does not explain why nearly all of them fall about half
short of what their own measurement permits. Roughly half the attainable signal
is being left on the table, and nothing tested this session recovers any of it.

## What the constraint is

**The number of independent (area, outcome) observations available to learn a
mapping from ~300 collinear predictors — not the number of people surveyed, and
not the choice of learner.**

Four independent lines of evidence point at it, and they are the only things in
the project's record with a consistent direction:

1. **The only monotone dose–response in the whole project is countries.** Each
   added training country buys about +0.05 of transported rank accuracy,
   monotone across all four specifications tested, 1→2→3 countries (TC-02).
2. **The only intervention that worked this session adds observations per area.**
   MT-01 (multi-trait BLUP, borrowing across a survey's biomarkers) gains +0.038,
   6 of 6 country × target blocks, and the gain concentrates exactly where
   single-cell estimation is weakest (weak tercile +0.086 against strong +0.018,
   ρ = −0.388, p = 0.020) — the shrinkage signature of an n-limited problem.
3. **Pooling that does *not* add observations fails.** EB-01 pooled the
   *parameters* across all 24 cells by empirical-Bayes moderation and was null
   (6/18, 8/18 in-country). So it is not pooling per se; it is pooling that
   supplies new observations.
4. **The theory says the penalty bites at this n, and it does.** The SuperLearner
   oracle inequality carries a term of order (1 + log K)/n; at n = 30 with
   K ≈ 120 that is ≈ 0.16 before the leading constant, and SL-07 measured a loss
   of the predicted sign and rough magnitude.

## How to loosen it, in order of evidential support

**1. More areas, by adding countries.** The strongest evidence in the project
(+0.05 per country, monotone), and XV-01/02 already showed the recipe transports
to six countries never used in training, including off-continent, with global
SoilGrids substituting for Africa-only iSDA at no cost. WHO VMNIS carries
sub-national deposits that served as a held-out label set. This is a
data-acquisition task, not a modelling one, and it is the only lever with a
demonstrated dose–response.

**2. More outcomes per area, modelled jointly.** MT-01, pre-registered before
the next country. Expect the gain in the weak cells specifically, which is where
the mechanism says it lives and where the project currently reports failures.

**3. Treat the aggregation rung as a lever, because the project has already
measured the trade.** Transport is 0.50 at the first sub-national tier, 0.31 at
Admin-2 (19/22 cells), and 0.24 at the cluster (CL-01). Coarser units are fewer
but individually more reliable, and transport improves monotonically as they get
coarser. Finer is not better: the cluster track has 1.6× more units and lost.
This states an explicit tension rather than a fix — policy needs Admin-2, the
statistics prefer Admin-1 — and it means an Admin-1 model benchmarked down to
Admin-2 deserves a fair test as a *deployment* option, not just as a diagnostic.

**4. Free reliability, from a reporting choice.** RL-01: dichotomising at a
clinical cutoff costs **0.203** of district reliability (level more reliable in
25 of 27 cells), the cost is worst for rare deficiencies (Sierra Leone women's
B12 at 0.5% prevalence loses 1.147), and the reliability gap **tracks** the
accuracy gap (ρ = 0.507, p = 0.034). Preferring the continuous level wherever a
decision permits it buys accuracy for nothing.

**5. Stop modelling the cells whose own data says there is nothing there.**
VF-01 found Fay-Herriot estimates the between-district variance as exactly zero
in about a quarter of cells, agreeing with RL-01's split-half on the weakest.
CV-01 found 86% of direct district estimates exceed a 33.3% coefficient of
variation. Reporting a model failure in a cell with no measurable district
signal misattributes the cause.

**6. Make post-hoc selection legitimate, rather than loosening the constraint.**
Every result in this folder is chosen on the same four countries, which is why
none can be quoted without a caveat. Data-adaptive target parameters (van der
Laan, Hubbard & Pfeiffer) define the estimand *as* whatever the selection
algorithm picks on a parameter-generating split and do inference on a held-out
split. Thin at four countries, but it converts a standing caveat into a defined
parameter.

## The single most concentrated opportunity: zinc

Malawi zinc is the largest block of unexplained headroom in the study, by more
than double:

| cell | reliability | r_max | achieved | headroom |
|---|---|---|---|---|
| Malawi child_zinc | 0.642 | 0.801 | **−0.105** | **0.906** |
| Malawi women_zinc | 0.717 | 0.847 | **−0.007** | **0.853** |

It is **well measured** (reliability 0.64–0.72, above iron and vitamin A) and
**completely unpredicted** (r_share ≈ 0). Measurement noise is therefore
excluded as the explanation. That is a precise quantification of ZN-01 — "a
reliable target the proxies do not touch" — and it makes zinc the highest-value
single target in the project: a reliably measured outcome with ~0.9 of
attainable correlation unexplained. GN-01 already ruled out measured grain zinc
and isotopically exchangeable soil zinc as the answer.

## One caution about a tempting number

`r_share` falls with the number of districts (Gambia 30 districts 0.71, Ghana 75
0.57, Malawi 87 0.36; ρ = −0.561, p = 0.015 over 18 cells). **Do not read that
as "more areas is worse".** District count is *perfectly confounded* with
country here — three countries, three counts — and the cell-level test treats
nested cells as independent. It is a country effect (consistent with the
project's own "Ghana teaches, Gambia learnable") and the design cannot separate
the two.

---

## ZW-01 (2026-09-28) — why zinc is unpredicted, and a correction to HR-02

**Answer: there is no Admin-2 geography in zinc to predict.** Its reproducible
between-unit variation is *within-district cluster* variation, which Admin-2
covariates cannot reach by construction.

Nested variance components of the continuous biomarker (district / cluster-
within-district / residual, moment-matched), against what the index achieves:

| cell | district share | cluster share | achieved |
|---|---|---|---|
| Malawi women_selenium | **0.669** | 0.046 | — |
| Malawi child_selenium | **0.544** | 0.117 | — |
| Malawi **women_b12** | **0.417** | 0.040 | **0.694** |
| Malawi women_folate | 0.341 | 0.000 | 0.412 |
| Malawi women_zinc | 0.113 | 0.124 | **−0.007** |
| Malawi **child_zinc** | **0.000** | **0.172** | **−0.105** |
| Malawi women_vitA | 0.000 | 0.053 | 0.302 |

The ordering is the story: the well-predicted cells (B12 0.417, folate 0.341,
selenium 0.544–0.669) are the ones with real district geography; zinc has the
**largest cluster share in the study (0.172)** and a district share
indistinguishable from zero. This confirms ZN-02 — "district 0 | TA 0 | cluster
0.013 … almost nothing separates districts" — by an independent route and on
women as well as children.

### The correction: HR-02's zinc headroom was an artefact of my own metric

HR-02 called zinc "the single most concentrated opportunity", with r_max 0.80
and headroom 0.906. That was wrong, and wrong for a reason I should have caught:
**`r_max = sqrt(split-half reliability)` counts cluster effects as district
signal wherever a district is a single cluster** — and 85% of Malawi districts
are. The split-half is an upper bound, not a ceiling.

Worse, **the project already had the right number and I did not use it.**
`results/tables/protocol_v2/variance_components_ceiling.csv` (VC-01, script 34)
carries `ceiling_vc`, a cluster-adjusted honest ceiling, alongside
`ceiling_within_vc`, the split-half bound, and `cluster_share_of_ceiling`. For
Malawi child zinc on the level it reports:

- honest ceiling `ceiling_vc` = **0.000**
- split-half ceiling `ceiling_within_vc` = 0.770
- `cluster_share_of_ceiling` = **1.000**

So the project's own table says 100% of child zinc's apparent ceiling is cluster
effect and the attainable district correlation is zero. My 0.906 of "headroom"
does not exist. Zinc is not the largest opportunity in the study; it is the
clearest case of *no* opportunity, and the record said so before I started.

### What that does to the binding-constraint numbers

HR-02's headroom used the inflated split-half bound throughout. On the honest
ceiling, averaged over the 27 admin-2 cells:

| | prevalence | level |
|---|---|---|
| honest ceiling (`ceiling_vc`) | 0.474 | **0.591** |
| split-half bound (`ceiling_within_vc`) | 0.592 | 0.735 |
| achieved (index) | 0.301 | **0.399** |
| **r_share against the honest ceiling** | **0.64** | **0.68** |
| index at or above the honest ceiling | 2 of 27 | 1 of 27 |
| cluster share of the split-half bound | 0.217 | 0.214 |

So the model is getting about **68% of what is honestly attainable on the
level**, not the 54% HR-02 reported. **The remaining headroom is roughly
0.19, not 0.35 — about half what I claimed.** That tightens the
binding-constraint conclusion rather than changing its direction: the
observation-count argument stands (it rests on TC-02's monotone dose–response,
MT-01's shrinkage signature and EB-01's null, none of which used r_max), but
there is materially less room left than I said.

**Process note, second occurrence.** After LT-01 I wrote a standing rule for this
folder: search `docs/findings/` and `scripts/protocol_v2/` for the question
before running a probe. I did not apply it to HR-02, and an existing results
table contradicted my headline. The rule needs to extend to *results tables*,
not just findings documents and scripts.
