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
