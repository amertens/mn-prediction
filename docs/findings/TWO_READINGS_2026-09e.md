# Two readings of the same results — 3 September 2026 (revision e)

> **Superseded on 2026-09-04 by `TWO_READINGS_2026-09f.md`.** The regional
> transport figure below (0.50–0.56, 12 of 12) was computed on three countries
> because of a population-file join bug (BUG-01); on all 22 cells it is 0.29,
> 17 positive. Revision f carries the corrected numbers and the 4 September
> results. This file is kept as the record of what was believed on 3 September.

Supersedes `TWO_READINGS_2026-09d.md`, which is preserved beside it. Revision d
reported the corrected protocol's numbers. This revision folds in the fifteen
sandbox experiments of 2–3 September (`docs/findings/SANDBOX_LOG_2026-09.md`,
scripts `protocol_v2/19`–`33`), the directional correction to the association
scan, and the null calibration of the transport counts. Both abstracts below
describe the same evidence; they differ in what they choose to emphasise, and
each is written so that a reader who checks the sources finds nothing
overstated.

Every number comes from `results/tables/protocol_v2/` under one protocol:
replicated folds, three separately scored estimands, information-matched
baselines, precision-weighted scoring, a 451-predictor vocabulary reduced to
domain principal components, and — new in this revision — permutation nulls
for the transport counts and cluster-split reliability ceilings.

---

## Reading A: the optimistic abstract

**Remotely sensed climate and soil rank districts for micronutrient targeting
across borders, from as little as one training country, and the simplest
estimator is the best one**

*Background.* Sub-national targeting of micronutrient programmes is limited by
the cost of biomarker surveys. Earlier negative findings — including our own —
rested on evaluation designs whose power and information symmetry had not been
quantified; correcting them left a positive result whose *mechanism* was
unexamined. This revision identifies it.

*Methods.* Across 24 country–outcome cells in four national biomarker surveys
(Gambia, Ghana, Malawi, Sierra Leone; 146 districts on a common administrative
rung) we scored a 451-predictor vocabulary, reduced to principal components
within 18 domains, under three separately defined estimands with
information-matched baselines: in-fill within a surveyed country, extrapolation
to an unvisited region, and leave-one-country-out transport. Random folds were
replicated ten times; districts were weighted by effective sample size from a
measured design effect. Transport counts were calibrated against a
country-block permutation null; domain contributions were measured by
drop-one and only-one ablation; the reliability ceiling was recomputed
splitting survey clusters rather than respondents.

*Results.* **Covariates beat the honest survey baseline within a country**:
against a jackknifed regional mean, rank correlation 0.398 vs 0.320 on the
biomarker level and 0.286 vs 0.193 on prevalence, better in 15 of 18 cells.
**Rankings transport to a country never seen in training**: district-level
mean Spearman 0.28–0.31 (19–21 of 22 cells positive), regional-level 0.50–0.56
(12 of 12), and the mean correlation exceeds every one of 500 country-block
permutation replicates at both tiers (null 95th percentiles 0.16 and 0.08).
**Two remotely sensed domains carry it.** Removing soil or climate costs
0.03–0.04 of transported accuracy; removing any of ten survey-derived domains
costs nothing or helps. A two-domain climate + soil index trained on a *single*
country ranks a second country's districts positively in 98 percent of fits,
transports better than the full index at the district rung (0.368 vs 0.252,
positive in 22 of 22) and equally at the regional tier, and still gains with
each added training country (+0.03–0.04, against +0.05 for the full index).
**Capacity is a liability, and the ensemble confirms it.** A zero-tuning
index over domain summaries beats a penalised fit over all 451 columns by up to
0.38, beats every tuned learner, and beats a SuperLearner that contains it
under four meta-learner losses; the ensemble at best recovers the index. The
associations it rests on replicate in sign across all four countries and are
nutrient-specific — legume cultivation and cattle ownership for child vitamin A
and iron, soil chemistry and night-time temperature for women's iron, child
wasting for women's vitamin A — describing a rural agro-ecological gradient
that is the same thing in every country because it is measured from orbit.

*Conclusion.* The deliverable is a **ranked district priority list** for
countries with or without a survey, carried by layers that exist for every
country on earth. We have pre-registered seven quantitative predictions for
the next two countries, including that the two-domain index will transport at
least as well as the full model; if they hold, the recipe is a handful of
freely available rasters and no tuning.

---

## Reading B: the pessimistic abstract

**Geography does the work; the covariates' contribution is a ranking, mostly
regional, and several of the numbers that look strongest are upper bounds**

*Background.* Correcting an over-stated negative does not establish a positive,
and identifying a mechanism does not license the numbers found while
identifying it. This reading states what survives when each result is scored
against the cheapest alternative, the fairest baseline, and its own null.

*Methods.* As above, with attention to what each arm adds over the cheapest
comparator and to how each headline behaves under a null.

*Results.* **Within a surveyed country geography is sufficient.** A
covariate-free spatial smoother scores 0.391 where the full model scores
0.398; covariates fitted to the smoother's residuals improve it in 5–7 of 18
cells with a median gain of zero. **The covariate advantage over survey
averages is a ranking advantage, not a reach advantage.** The model's
worst-ranked fifth captures 24 percent of national deficiency burden against
21 percent for a jackknifed regional mean, better in only 8–11 of 18 cells,
and against 30 percent for a regional mean that includes the district's own
respondents; the margin concentrates where regions hold three to five
districts, exactly where the jackknifed baseline is most penalised. **The
model does not substitute for survey sample**: fitted symmetrically on a
reduced survey, its district error degrades in lockstep with the survey's
(+3.36 vs +3.33 pp from full to 15 percent of sample). **Transport is a
ranking at the regional tier and nothing more.** Absolute levels do not
transport; anchoring to a country's own national value halves error, anchoring
to a published external estimate triples it. Burden captured under transport
is at chance (lift 1.01). "Twelve of twelve" overstates the evidence:
with correlated outcomes inside a country, no transport still yields 15–16 of
22 positive cells five percent of the time; only the mean correlation is
unambiguous. **The climate + soil result was found on the data it is scored
on** and is reported as a prediction, not a result; honest nested selection
recovers +0.03 on the biomarker level and nothing on prevalence. **The
reliability ceiling — and so the headroom — is an upper bound.** In three of
four countries most districts are a single survey cluster (57, 83 and 85
percent of units), so the split-half ceiling counts everything a cluster
shares as geography; on multi-cluster units 16–27 percent of it is cluster
effect, and at Malawi's district rung only child iron keeps reliable geography
once clusters are split. Malawi, where nearly all apparent headroom sat, was
mostly artefact; zinc, the one outcome with a high ceiling and no proxy
signal, cannot be distinguished from a collection artefact at this rung.
**The replicated associations are not dietary effects**: legumes and cattle
track *more* deficiency, the reverse of the household-level relationship —
markers of where subsistence agriculture is, useful for deciding where to
look, uninformative about what to change. An earlier version of this document
had those signs inverted.

*Conclusion.* The honest product is a **ranked list, strongest at the first
sub-national tier, for countries without a survey** — the one setting in which
no spatial smoother or survey average can be fitted at all. Within a surveyed
country a smoother delivers most of what the covariates do. Every claim about
prevalence levels, burden reached, survey savings or headroom should be read
as an upper bound or withdrawn.

---

## What moved since revision d

| Revision d claim | Status in e |
|:---|:---|
| Transport 0.288, positive in 21 of 22; "targeting lift 1.7×" | **Restated.** Ranking stands (0.28–0.31 district, 0.50–0.56 regional; mean ρ beats 500/500 null replicates). The 1.7× was top-*quartile overlap*, not burden; burden capture under transport is at chance (1.01). Cell counts have a fat null and are now secondary. |
| "Simplest estimator wins" | **Strengthened.** Beats a SuperLearner containing it under MSE, rank, weighted-rank and burden-capture meta-learners (SL-01…04). |
| "Reliability ceiling 0.47–0.61 leaves substantial headroom" | **Downgraded to upper bound.** Most units are single clusters; 16–27% of the ceiling is cluster effect on multi-cluster units; Malawi's headroom was mostly artefact (CE-01, MW-01). |
| Covariates vs spatial smoother, 0.398 vs 0.391 | **Stands, sharpened:** nested residual test null in all four combos (5–7 of 18). |
| Women's vitamin A "not a usable district target" | **Qualified.** True at TA/district rung where most cells are zero; at Admin-1 it transports in every cell and tracks child wasting (4/4). |
| Malaria null | **Stands, qualified:** no burden indicator replicates; IRS *coverage* reaches z 3.66 (FWER 0.14). |
| Directional claims (legumes, cattle protective) | **Inverted and corrected 2026-09-02.** Both track more deficiency; agro-ecology axis. |
| — | **New:** transport carried by climate + soil (DA-01/02/03, TC-02), pre-registered for Ethiopia/Pakistan. |
| — | **New:** model does not substitute for survey sample (G4-02). |
| — | **New:** burden margin over survey averages is small and partly jackknife artefact (CF-01, R6-02). |
| — | **New:** production area-level SuperLearner had no weights and no blocking; fixed (MAE 18.1 → 15.0 pp). |

## Paper outline

Unchanged from revision d in structure; Section 2 (what travels with
deficiency) now leads with the agro-ecology axis, Section 4 (prediction) adds
the domain-ablation and null-calibration results, and Section 5 (protocol) adds
the cluster-split ceiling. See `SLIDE_OUTLINE_2026-09.md` for the talk version.
