# Two readings of the same results — 1 September 2026 (revision c)

> **SUPERSEDED by `TWO_READINGS_2026-09d.md`.** Revision c reported the methods
> audit. Revision d reports what happened when its five recommendations were
> implemented and measured, when the unused data sources were wired in, and when
> two further defects were found in the process — a Ghana survey weight that was
> step 1 of a documented 3-step construction, and a modelled surface joined at
> the wrong administrative level. Revision c's `median r 0.058 / 0.206` is
> **superseded**; under the corrected protocol the figures are 0.286 and 0.398.
> A status table for every revision-c claim is at the end of revision d.


This supersedes `TWO_READINGS_2026-09b.md`, which is preserved beside it and now
carries a supersession banner. Both abstracts below describe the **same**
measurements. Neither contains a claim the other contradicts.

**What forced the revision.** A 31-agent methods audit and a six-probe signal
battery (report: <https://claude.ai/code/artifact/d65faecb-1196-4c71-863c-57105860ee90>;
scripts in `scripts/signal_probes/`, outputs in `results/tables/signal_probes/`)
established two things the previous revision could not have known. First, three
of the numbers that carried revision b's central negative are **artifacts of the
evaluation design rather than measurements of the data** — every one of them
verified by an independent adversarial check (24 of 24 confirmed, none refuted).
Second, the covariates carry **cross-country-replicated, transportable
between-region signal** that every significance test in the project was
structurally unable to detect, because all of them permute outcomes within
Admin-1 region and therefore keep that signal in every null draw.

The previous revision opened by saying the protocol matters more than any
modelling choice. That remains true, and is now the point: the protocol was not
merely consequential, it was in three places **aimed away from where the signal
is**, and in one place scored against a baseline the project had already
withdrawn.

---

## Reading A: the optimistic abstract

**Geospatial and demographic proxies carry replicated, transportable signal for
micronutrient deficiency, and earlier negative reports of this were measurement
artifacts**

*Background.* Sub-national targeting of micronutrient programmes is limited by
survey cost, and the standing hope has been that geospatial covariates could
extend surveys to unsurveyed districts and countries. Recent negative reports —
including our own — have concluded that they cannot. Those conclusions rest on
evaluation designs whose power against the relevant alternative was never
quantified.

*Methods.* Across 23–24 country–outcome cells in four national biomarker surveys
we tested covariate association with instruments matched to the alternative that
matters: association pooled across countries by random-effects meta-analysis and
calibrated by a **region-permutation family-wise null**; an omnibus dense-signal
statistic (mean r², the linear-kernel SKAT/globaltest statistic) against nulls of
graded spatial strictness; leave-one-country-out transport with a
**zero-hyperparameter** index; and, within each country, mass univariate testing
with Benjamini–Hochberg q-values and a Storey π₀ estimate over each country's
native ~590-predictor vocabulary. Skill was scored on both the deficiency
prevalence and the underlying biomarker concentration, under both random and
region-blocked folds.

*Results.* **Covariate signal is real, replicated, and transportable.** Seven
predictors on the prevalence scale and fourteen on the biomarker scale survive
family-wise permutation control with **sign agreement in all four countries** —
night land-surface temperature, soil chemistry, land-use embeddings, prior-round
child wasting and women's education — and 8 of 18 domain scores survive. An
index over those domain scores, with signs and weights taken **only from the
other three countries and nothing fitted on the held-out one**, ranks the
held-out country's regions at mean Spearman **0.309** across 21 cells, positive
in 18, global permutation p = **0.0005**. Within countries, hundreds of
associations survive BH control in three of four countries, with Storey π₀
falling as low as **0.26** — implying that a majority of the vocabulary carries
some association in the strongest cells — and in Malawi — the only country with genuine district-level resolution — the
covariate block carries information **beyond any smooth spatial surface** in 4 of
8 cells. Two corrections raise measured skill: the published within-country
median r of 0.058 rests on a single unreplicated fold draw, and re-randomising
that same protocol gives **0.217** (leave-one-region-out: **0.235**); and
modelling the biomarker concentration rather than the dichotomised indicator
takes cross-validated skill from **−0.03 to +0.24** under the project's own
strict folds. Against the only information-symmetric covariate-free baseline —
the jackknifed regional mean — the covariate model **wins** (r 0.156 vs 0.076),
and under fold-matched comparison covariates and a spatial smoother are
statistically tied. Ranking transports: targeting the worst-ranked fifth of
districts reaches **1.23×** the burden of untargeted allocation in a country held
out entirely, above chance in 12 of 16 cells.

*Conclusion.* The deliverable is a **ranked priority list** for countries without
a survey, and it rests on measured, replicated, cross-country signal rather than
on hope. The validated reliability ceiling (0.47–0.61) leaves substantial
headroom that no model has yet claimed.

---

## Reading B: the pessimistic abstract

**The signal is real, small, and mostly geography; nothing here is yet a
deployable district map**

*Background.* A substantial literature reports that environmental and remotely
sensed covariates predict micronutrient deficiency sub-nationally, usually
validated against the same survey's estimates under folds permitting spatial
leakage. Our own attempt to falsify that literature over-corrected: several of
our negative headlines were design artifacts. Correcting them does not vindicate
the literature.

*Methods.* As above, with the emphasis on what survives strict spatial control
and on the size, not merely the existence, of the effects.

*Results.* Every effect that survives is **modest**. The best transport result is
a rank correlation near **0.31** from an untuned index; the corrected
within-country skill is **0.22–0.24**, not the 0.058 previously published but far
from usable. Residualising on region centroids collapses the within-country
domain associations entirely (max pooled |z| **0.73**), and under a
Moran-eigenvector null the covariate block adds nothing beyond a smooth spatial
surface in **Gambia (0 of 4 cells) and Ghana (0 of 6)** — only Malawi, with 87
districts, retains it. Sierra Leone, at 14 areas, shows **nothing at all** on
either outcome scale under any null, and is uninformative rather than negative.
The malaria domain is null in every probe. Women's vitamin A prevalence is
exactly zero in **47–87%** of districts and is not a usable district-level
modelling target. Transport rests on **four countries**, with the B12 and folate
cells drawing on only one or two training countries, so transportability is
asserted on an n that cannot support it. The surviving predictors are existence
proofs from a low-power design, not a ranked shortlist: the family-wise threshold
implies roughly 4% power at plausible effect sizes. And several mechanical
defects remain in the production path — a name-only Admin-2 join that fans
Malawi's 87 surveyed districts to 90 training rows, corrupted DHS water-access
variables, a pooled cross-country covariate set with 22% of columns
scale-incomparable and uncentered, and Ghana's entire 153-column DHS block
silently removed by complete-case filters.

*Conclusion.* The honest claim is **information, not accuracy**: proxies carry
replicated between-region signal that transports, and they do not yet support a
district prevalence map. Published sub-national accuracies should still be
treated as unreliable unless they report the reliability of their validation
target, use administratively blocked folds, calibrate multiplicity against a
null that preserves the structure they claim to exploit, and state whether their
baseline saw the held-out unit. Applied here, those conditions shrink the effect
to a ranking — they no longer remove it.

---

## What moved since revision b

Every item below was verified independently before being written down.

**1. The within-country headline was a fold-draw accident.** `median r 0.058`
comes from one random 3-fold assignment of Admin-1 regions
(`wsl1_within_country.R:117-119`, seed re-executed per cell so one partition is
shared by all of a country's outcomes and draw errors correlate). That draw sits
in the bottom decile of its own protocol's distribution. Re-randomised 20×:
median **0.217**. Leave-one-region-out, same estimator and predictors:
**0.235**. Gambia child vitamin A: committed 0.058 against a re-draw median of
0.449. Two auditors reproduced the committed table exactly before varying it, so
the comparison is like-for-like. **Every number downstream of that single draw
moves, including the continuous-vs-binary gain and the beats-the-null tallies in
both abstracts of revision b.**

**2. "0 of 294 predictors survive FDR" was unattainable by construction.** With
999 permutations the smallest attainable p is 0.001, so a single true predictor
of *arbitrary* strength receives BH q ≥ 0.296; a within-cell survivor would
require ≥ 6 predictors simultaneously tied at the permutation floor, and the
maximum observed in any cell is 2. The result was guaranteed before the data were
seen. The same table's conventional-null columns carry **476** within-cell and
**310** global BH survivors. Median power of the screen at latent ρ = 0.3 is
**0.005**. The claim is withdrawn as evidence of absence.

**3. The covariate-free "number to beat" had already been withdrawn.** The flat
regional mean (r 0.516) assigns each held-out district its region's survey
estimate *computed from that district's own respondents*. The project ran the
jackknife control on 2026-08-31 and withdrew the claim; the jackknifed arm scores
**r 0.076**, losing to the covariate model at 0.156. `EVALUATION_PROTOCOL.md`
rule 3.6 was written two hours earlier and still names 0.516 "the number to
beat", and the submission harness was sized so that this non-compliant arm could
be scored at all.

**4. The spatial-smoother comparison was protocol-mismatched.** The "spatial
beats every covariate arm" rows score a leave-one-**district**-out smoother
against leave-one-**region**-out covariate arms, on 18 cells against 24 (the
smoother silently skips all six Sierra Leone cells). Under the one fold-symmetric
comparison in the repository the two are tied: r **0.156 vs 0.137**, covariates
winning **10 of 24** cells — and revision b's "6 of 24" is a miscount against its
own source table.

**5. Every "no signal" test conditioned away the signal.** The per-cell FDR
screen, the WS-A pooled scan and the WS-F confirmatory test all permute outcomes
within Admin-1 only, so each null draw retains the predictor's full between-region
association. Against purely between-region signal these tests have power exactly
α. Under that handicap WS-A still produced two family-wise survivors at ~4%
power, which is better read as positive evidence than as a near-miss.

**6. The mechanistic predictor set is enriched, not clean.** WS-F reported "of 31
chosen in advance, none survive". The set in fact shows **7 of 29 predictors
nominally significant against 1.45 expected** (4.8×, binomial p ≈ 5×10⁻⁴) with
mechanistically coherent directions — cereal share +, soil zinc −, pulses +,
soil organic carbon −, soil pH +. It is the clearest distributed signal in the
project and was reported as its cleanest null.

**7. The domain-ablation claim was overstated and partly mechanical.** "Every
single domain performs better alone than all 373 together" is not supported by
its own table: a domain alone beats the full set in **132 of 300** rows (44%),
and per domain in 9 of 17. The supportable version is that the *best* single
domain beats the full set in 17 of 18 cells. "Deleting a domain helps in 72% of
cases" counts 100 exact-zero rows as help; strictly it is **117 of 300** (39%).
Both figures are also inflated by the fixed top-20 screen, under which at most
~5% of predictors ever enter a model, so a domain's drop cost is mechanically
zero whenever it is never selected.

**8. "The choice of null moves results by 2.80 pp" does not reproduce.** Direct
measurement on the same cells gives **0.2–0.4 pp** (max 1.78). No committed table
computing 2.80 was found.

**9. The noise ceiling in the outward documents is the invalidated one.** The
manuscript still frames the negative through the analytic ceiling (median r_max
0.098, "92% of the ceiling captured") that the project's own WS1 proved ~4.7×
biased low against simulated truth. The binding ceiling is the empirical one,
**0.47–0.61**, which means most cells retain substantial recoverable
between-district signal that no model has captured. Cells whose empirical ceiling
truncates to zero are *unestimable*, not signal-free: Malawi women's B12 has a
ceiling of 0.000 and beyond-geography covariate association at p = 0.002.

**10. The withdrawn anchoring gain is still the headline recommendation
downstream.** The regional-anchoring result (mean r 0.164 → 0.413) was withdrawn
as circular — jackknifed, the hard anchor beats no-anchor in 8 of 24 cells at
r 0.147 vs 0.156 — yet it remains in the manuscript abstract, the slides, and the
deployed dashboard.

**11. New positive result: signal that replicates and transports.** Detailed in
Reading A. The single most important number is the leave-one-country-out index at
mean Spearman 0.309 over 21 cells (p = 0.0005), because it is the only test in
the project that measures what a country with no survey would actually receive,
and it involves no tuning of any kind.

**12. Harmonisation is exonerated.** Every within-country hit at q < 0.05 in the
native pre-harmonisation vocabularies has its base name inside the shared
373-predictor set. Signal was not lost in harmonisation, consistent with the
project's earlier finding that doubling Gambia's raw extraction changed
downstream results bit-for-bit.

**13. Linkage is exonerated, by direct geometry.** Every survey cluster's GPS
coordinates were point-in-polygon tested against the covariate polygons: 100%
agreement in Gambia, Ghana and Sierra Leone, and Malawi's only three
disagreements are the documented lake-displacement clusters. All headline joins
are pair-keyed with full coverage and no fan-out. **The negative was manufactured
in the evaluation and inference design, not in the data.**

---

## Limitations that bias the comparison

**Against Reading A.** The pooled scans use random-effects meta-analysis over
four country clusters, which the audit criticised elsewhere in this project for
producing a heavy-tailed null; here the null is calibrated by permutation on the
identical statistic, so the *level* is controlled by construction, but power is
low and the survivors should be read as existence proofs rather than as a ranked
shortlist. Between-region association cannot distinguish covariate information
from a smooth gradient that happens to track deficiency, and within Gambia and
Ghana it demonstrably cannot be separated. The within-country mass-univariate
results use the exchangeable null, which counts spatially structured association
as signal; they are strong evidence against a linkage or cleaning fault, not
evidence of hundreds of independently useful predictors. Covariates are
aggregated to Admin-1 with unweighted means, and population weighting was not
swept.

**Against Reading B.** Three of its predecessors' central negatives are now known
to be design artifacts, and the remaining negative statements are conditional on
a within-region framing that discards roughly 58% of covariate variance — of
which perhaps 30–40 percentage points is structure in excess of noise. Its
strongest surviving negative, the collapse under spatial partialling, is
irrelevant to the deployment case it is most often quoted against: no spatial
smoother can be fitted for a country with no survey, so signal that is
collinear with geography is still the only signal available there. And Sierra
Leone, which contributes negative cells throughout, is a 14-area country in which
nothing is detectable at any effect size.

**Against both.** Four countries remains a small n for any cross-country claim.
The audit was a code-and-statistics audit: it did not re-extract rasters, re-run
GEE, or re-derive biomarker cut-offs from source documents, so the mechanical
defects it reports are verified but the absence of others is not proven. And
nothing measured here yet constitutes a deployable district map — the next round
has to earn that separately.

---

## Status of the numbers in revision b

| Revision b claim | Status |
|:---|:---|
| Median r 0.058 (prevalence), 0.206 (level) | **Withdrawn** — single unlucky fold draw; 0.217/0.235 on replication |
| "0 of 294 survive FDR", "2 of 294 pooled" | **Reframed** — the zero was unattainable by design; pooled survivors stand |
| "31 mechanistic predictors, none survive" | **Withdrawn** — set is enriched 4.8× |
| "No domain load-bearing; every domain better alone" | **Corrected** — 44% of rows, 9 of 17 domains; partly a screen artifact |
| "Unshrunken survey estimate outranks every model" | **Qualified** — reverses above ρ ≈ 0.6, a regime Gambia plausibly occupies |
| "Fold construction moves results by 0.262" | **Stands**, and is now the central finding rather than a caveat |
| "Choice of null moves results by 2.80 pp" | **Withdrawn** — does not reproduce; 0.2–0.4 pp |
| WHO band verification (item 5) | **Stands** |
| Individual-specific cut-offs (item 6) | **Stands** |
| Two withdrawn cells (item 7) | **Stands**, but both are still present in machine-readable aggregates and the dashboard bundle |
| Predictor-breadth claim (item 8) | **Stands**, and is reinforced by the native-vocabulary scan |
