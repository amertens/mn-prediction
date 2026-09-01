# Two readings of the same results — 1 September 2026 (revision d)

Supersedes `TWO_READINGS_2026-09c.md`, which is preserved beside it. Revision c
reported the methods audit. This revision reports what happened when the audit's
five recommendations were **implemented and measured**, when the unused data
sources were **wired in**, and when two further defects were **found in the
process of doing so**.

Every number below comes from `results/tables/protocol_v2/`, produced under one
protocol: replicated folds, three separately-scored estimands, information-matched
baselines, precision-weighted scoring, and a 451-predictor vocabulary. Details and
a plain-language guide to every arm are in `docs/findings/PROTOCOL_V2.md`.

---

## Reading A: the optimistic abstract

**Geospatial and demographic proxies rank districts for micronutrient targeting,
they beat the honest survey-only baseline, and they transport to countries that
have never been surveyed**

*Background.* Sub-national targeting of micronutrient programmes is limited by
survey cost. Reports that geospatial covariates can extend surveys — including
our own negative ones — have rested on evaluation designs whose power and
information symmetry were never quantified.

*Methods.* Across 24 country–outcome cells in four national biomarker surveys
(Gambia 2021, Ghana 2017, Malawi 2015–16, Sierra Leone 2013) we scored a
451-predictor vocabulary under three separately-defined estimands, each against
a baseline that saw the same information: in-fill within a surveyed country,
extrapolation to an unvisited region, and leave-one-country-out transport. Folds
with any random component were replicated ten times; districts were weighted by
an effective sample size derived from a measured design effect; and both the
deficiency prevalence and the underlying biomarker concentration were modelled.

*Results.* **Covariates beat the honest covariate-free baseline in every
within-country panel.** Against the jackknifed regional survey mean — the
comparator that does not read the held-out district's own respondents — a
spatial-plus-covariate model reaches rank correlation **0.398 against 0.320** on
the biomarker level, and a covariate-only index reaches **0.286 against 0.193**
on prevalence. **Ranking transports to an unsurveyed country**: an index built
on three countries ranks districts in the fourth at mean Spearman **0.288,
positive in 21 of 22 cells**, capturing **42–44 percent** of the truly
worst-quarter districts against a chance level of 25 percent — a targeting lift
of about 1.7×. **The simplest estimator wins.** A zero-tuning index over domain
summaries — nothing selected, nothing fitted — is the best arm on prevalence in
both within-country panels, beating a penalised regression over all 451 raw
columns by 0.19 and 0.38. And the validated reliability ceiling of 0.47–0.61
leaves substantial headroom above every number here.

*Conclusion.* The deliverable is a **ranked district priority list**, available
now, for countries with or without a survey. It rests on replicated folds and
information-matched baselines rather than on a favourable protocol.

---

## Reading B: the pessimistic abstract

**Within a surveyed country, geography does the work; the covariates' distinct
contribution is narrow, and no result here is a prevalence map**

*Background.* Correcting an over-stated negative does not establish a positive.
Revision c withdrew three headline negatives as artifacts of the evaluation
design. What survives correction is smaller and more specific than either the
prior negative or the prior optimism.

*Methods.* As above, with attention to what each arm adds over the cheapest
alternative.

*Results.* **A covariate-free spatial smoother, using no predictors at all,
scores 0.391 where the full covariate model scores 0.398.** Adding 451
predictors to geography buys **0.007**. On the prevalence target the covariate
index does lead (0.286 against 0.268), but the margin is inside the resolvable
range at these sample sizes. **Effect sizes remain modest**: the best transport
result is a rank correlation near 0.29, and prevalence transport is 0.15.
**Nothing here is a level**: transport outcomes are standardised within country
before pooling because cross-survey biomarker offsets are large, so no
prevalence is claimed and none should be reported. **More predictors made
prediction worse**: a penalised fit over the full vocabulary is *negative*
(−0.107) under region extrapolation, positive in only 6 of 24 cells, and adding
61 predictors today made it worse, not better. Sierra Leone, at 14 districts,
is uninformative for any within-country claim; women's vitamin A is exactly zero
in 47–87 percent of districts and is not a usable district target; the malaria
domain is null in every probe. Four countries cannot support a transportability
claim, and the B12 and folate transport cells rest on one or two training
countries.

*Conclusion.* The honest product is a **ranking, not a map**, and within a
surveyed country a spatial smoother delivers most of it. The covariates earn
their place in exactly two situations: when the comparator is a survey-derived
baseline that must be jackknifed to be fair, and when there is no survey at all
and no smoother can be fitted.

---

## What moved since revision c

**1. The five fixes were implemented, and each was measured in isolation.**
Replication: the median range of one cell's score across ten fold draws is
**0.174**, maximum 0.705. Effective sample size: measured design effects are
**1.0–5.6, median ≈ 2.4**, not the assumed 1.5, so 77 percent of Malawi's
districts have n_eff < 5. Continuous target: **+0.10 to +0.11** in all three
estimands. Domain representation: the zero-tuning index beats the 451-column
fit by up to **+0.38**. Estimand split: the same arm scores 0.32 under in-fill
and 0.19 against a baseline that previously read its own answer.

**2. A survey weight was wrong for four years, and the survey's own report says
so.** Ghana's configured `gw_sWeight` takes three distinct values, one per
stratum. Appendix 11 of the Ghana Micronutrient Survey 2017 report documents a
three-step construction and states that "the final survey weights for each PSU
were calculated by multiplying the standardized PSU weights by the standardized
stratum weights". Verified in the data: `gw_PSU_weight * gw_sWeight ==
gw_PSUStrat_weight` to 2e-06. The pipeline was using step 1 of 3. Corrected,
national prevalence moves 1.3–1.5 pp and individual districts up to 4.8 pp.
Gambia, Malawi and Sierra Leone were audited at the same time and are correct.

**3. A modelled surface was being joined at the wrong administrative level.**
IHME's `adm2_name` is not Admin-2 in every country: for Malawi it is the 27
districts, which are this project's **Admin-1**, while Admin-2 is the ~239
Traditional Authorities. Measured, 27 of 30 Malawi names match Admin-1 exactly
and **0 match Admin-2** — yet the existing harmonisation fuzzy-matched district
names onto TA names and "recovered" 23 of 30, every one a silent mis-linkage.

**4. Three substantial sources reached zero predictors and now reach 61.**
IHME (93 GB, 13 topics) had been logged as "kept" and then deleted by the
all-four-countries filter. HFID is not health-facility data but food-insecurity
data carrying Food Consumption Score and reduced Coping Strategies Index at
Admin-2, of which only two IPC phase columns had ever been taken. Neither the
food-price series nor FAOSTAT supply reached the vocabulary at all. Added with
per-country coverage recorded instead of the filter, plus WorldPop age-sex,
GHS-SMOD urbanisation, Köppen and AEZ zone classes, and the Meta Relative
Wealth Index. Vocabulary 373 → 451.

**5. A caution of mine was wrong and is withdrawn.** Ten predictors were briefly
excluded by default as "near-outcome", on my own judgement. Checked against
source documentation: GFDx anaemia is WHO's 2011 estimate, which pre-dates every
survey here; GFDx zinc is derived from FAO food-balance-sheet availability with
no biomarker survey involved; IHME's anaemia surface is defined on haemoglobin,
a different biomarker from the ferritin, RBP, folate, B12 and zinc outcomes
modelled here. All are now included by default, flagged as modelled surfaces
with their provenance recorded.

**6. An earlier claim of mine was too broad and is narrowed.** "Pooling raw
predictors across countries destroys transport" was partly an artifact of
setting `standardize = FALSE` in the comparison. Re-run fairly: a penalised fit
barely notices the scaling scheme (0.212 vs 0.246), while any **averaging**
composite collapses without it (0.014 vs 0.174). The fix is load-bearing here
because the best-transporting arm is an average — but the general claim needed
its scope.

**7. More components per domain beat one.** The single sign-aligned mean is the
worst of seven representations tested; retaining principal components to 80
percent of each domain's variance raises transport from 0.151 to 0.255.
Supervision adds nothing, so the gain comes from representing more of each
domain's variance rather than from aiming it at the outcome.

---

## Status of revision c's numbers

| Revision c claim | Status in d |
|:---|:---|
| Median r 0.058 / 0.206 | **Superseded.** Under the corrected protocol: 0.286 (prev) and 0.398 (level) in-fill. |
| "0 of 294 survive FDR" was unattainable | **Stands.** |
| Number-to-beat 0.516 withdrawn | **Stands**, and `EVALUATION_PROTOCOL.md` rule 3.6 is now corrected to the jackknifed arm (0.076). |
| Covariates tie the spatial smoother | **Stands**, and sharpens: 0.398 vs 0.391 on level, covariates ahead on prevalence. |
| Transport 0.309 at Admin-1 (probe P3a) | **Complemented** by 0.288 at Admin-2 under the full protocol, positive in 21 of 22. |
| Ceiling 0.47–0.61 binding | **Stands.** |
| "6 of 24" miscount → 10 of 24 | **Stands**, corrected in three documents. |
| Malaria null; Sierra Leone uninformative | **Stands.** |

---

## Paper outline

A single paper, structured so the negative and the positive are the same
argument rather than competing ones.

**Title.** *What geospatial covariates can and cannot do for sub-national
micronutrient estimation: a protocol-first re-analysis of four national
biomarker surveys.*

**1. Introduction.** The targeting problem and the survey-cost constraint. The
literature's claim. The observation that motivates the paper: published
sub-national accuracies are rarely accompanied by the reliability of their
validation target, the information content of their baseline, or the
replication of their folds.

**2. Data.** Four surveys, 24 country–outcome cells, the 451-predictor
vocabulary and its 19 domains. A table of what each source contributes and at
what administrative level it natively lives — with the IHME admin-level finding
as a worked example of why that column matters.

**3. Protocol (the methodological core).** Three estimands with
information-matched baselines. Replicated folds, and the measured cost of not
replicating (median range 0.174). Effective sample size from measured design
effects (median 2.4, not 1.5). Continuous versus dichotomised targets. Within-
country normalisation before pooling. This section is the paper's contribution
even if every effect size were null.

**4. Results.** The leaderboard, by estimand. Covariates versus the jackknifed
baseline. Covariates versus geography. Transport, with top-quartile capture as
the decision-relevant metric. The dimension result: simplest estimator wins.

**5. What the protocol changes.** Side-by-side of each headline under the old
and corrected protocols — the paper's most transferable content, because the
same errors are available to any group doing this work.

**6. Limitations.** Four countries. Ranking not level. Geography does most of
the within-country work. Sierra Leone. Degenerate low-prevalence targets.

**7. Recommendations for the field.** The four reporting conditions, plus the
two added here: replicate any random fold and report the spread, and state the
information your baseline saw.

**Supplement.** Data-source provenance and assumptions; the variable annotation
sheet; the effect of each fix in isolation; full per-cell tables.
