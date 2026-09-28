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
