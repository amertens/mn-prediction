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
