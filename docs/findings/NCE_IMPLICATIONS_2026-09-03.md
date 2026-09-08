# What the 3 September sandbox results mean for the NCE draft

Read this before the next pass on `NCE_Proxy_2026-09-0x.docx`. Each item names
the sentence in the current draft, what the new evidence says, and the
wording I would use. Sources are entries in
[SANDBOX_LOG_2026-09.md](SANDBOX_LOG_2026-09.md); numbers are reproducible
from `scripts/protocol_v2/19`–`32`.

## Sentences to change

**1. The burden sentence.** Draft: *"directing an intervention to the fifth of
districts the model ranks worst reaches 24% of a country's total deficiency
burden, compared with 19% using regional survey averages … regional averages
are no better than chance … while the model is."*
Evidence (CF-01, R6-02): 24% vs 19% is correct as computed, but per cell the
model beats the jackknifed regional mean in only 8–11 of 18, the pooled margin
comes from countries whose regions hold 3–5 districts, and a regional mean that
*includes* the district reaches 30%. The ranking claim (0.29 vs 0.19; 15 of 18)
is the robust one.
Wording: keep the numbers, scope the claim — *"for districts a survey has not
reached, regional averages are no better than chance, while the model is"* —
and do not repeat the 26% as a headline elsewhere. Drop the Malawi "12%"
clause if it is still present.

**2. The survey-size sentence.** Draft (already softened): *"Whether the model
can substitute for part of that sample … is an additional objective we propose
to quantify."*
Evidence (G4-02): tested symmetrically, it does not substitute. The model's
district error degrades in lockstep with the survey's as the sample shrinks
(+3.36 vs +3.33 pp from 100% to 15%). What it adds is roughly constant at any
sample size: ~0.6 pp lower district error and a much better district ordering.
Wording: *"At any given survey size the model improves the ordering of
districts and modestly reduces district error; it does not reduce the sample a
survey needs. The extension will quantify the regional-versus-district
trade-off in survey design and what the model adds at each sample size."* Do
not frame savings as coming from a smaller sample.

**3. The headroom sentence (Edit 6, if used).** Draft: *"the models currently
capture about two-thirds of the accuracy that is attainable given survey
sampling noise."*
Evidence (CE-01, MW-01): 57–85% of Admin-2 units in three countries are single
survey clusters, so the ceiling counts cluster-level effects as geography; on
multi-cluster units 16–27% of the ceiling is cluster effect, and at Malawi's
district rung only child iron retains reliable geography once clusters are
split. Malawi — where almost all the apparent headroom sat — was largely
artefact.
Wording: *"By a split-half reliability estimate, the models capture roughly
two-thirds of the between-district variation the surveys can resolve; that
estimate is an upper bound on remaining headroom, because most districts are
represented by a single survey cluster."* Or omit the sentence — it now buys
less than it costs to defend.

## Sentences that got stronger

**4. Transport.** *Regional rankings transport in 12 of 12; district in 19 of
22* — **corrected 2026-09-04 (BUG-01).** The "12 of 12" was four outcomes in
three countries: a population-file spelling dropped Sierra Leone from every
population-joined run. On all 22 regional cells the mean rank correlation is
0.30 with 17 positive (0.43 and 12 of 12 in the two countries with eight or
more regions), against a permutation null whose 95th percentile is 0.16; with
the climate + soil index it is 0.45 (20 of 22). The district figure (mean 0.28 with 17 of 22 positive for the
zero-tuning index, 0.28 with 21 of 22 for the penalised fit; 0.37 and 22 of
22 with climate + soil) is unaffected. Use: *"district
and regional rankings transport to a country never used in training, with a
mean rank correlation of about 0.3 at both tiers, far outside a permutation
null (95th percentile 0.16); a two-domain index built from remotely sensed
climate and soil layers reaches 0.37 at district level and 0.46 at regional
level."* Do not quote 12 of 12 or 0.50–0.56 anywhere.

**4b. A design the NCE can name (AR-01, new).** A survey that measures only
the national prevalence — about five percent of a full survey's sample,
fifty to seventy respondents — combined with the transported ranking gives a
district error of 10 percentage points, against 20 for a district survey of
the same size and 12 for a regional one; the district survey needs about 40
percent of the full sample before it does better. It does not beat the survey
on burden captured. Suggested sentence: *"For a country without a
micronutrient survey, a small national biomarker sample combined with the
transported district ranking gives district-level estimates that a
conventional survey only matches at roughly eight times the sample."*

**5. Adding countries.** *Each added training country buys ~0.05 of
transported accuracy* — unchanged, and it holds for a two-domain index too
(TC-02: +0.03–0.04 per country, from a much higher start). Keep as written.

**6. What to collect (new, if wanted).** The transport result is carried by
remotely sensed climate and soil layers; ten survey-derived domains are dead
weight across borders (DA-01/02/03). A two-domain index trained on a *single*
country ranks a second country's districts positively in 98% of tests. This
is the concrete reason the extension's data-source work should prioritise the
remotely sensed layers — but it was identified on these four countries, so it
belongs in the NCE as a **pre-registered prediction** for Ethiopia and
Pakistan ([PREREGISTRATION_NEW_COUNTRIES_2026-09.md](PREREGISTRATION_NEW_COUNTRIES_2026-09.md)),
not as a result. Suggested sentence for the "add countries" activity: *"We
have pre-registered seven quantitative predictions for the new countries,
including that a parsimonious index built only from remotely sensed climate
and soil layers will transport district rankings at least as well as the full
model."*

## Sentences that are unaffected

The 15-of-18 ranking claim, the 0.29 vs 0.19 correlations, the Malawi 43% /
49% example, the 20% national figure (harmonised definition), and paragraph 6
as redlined. The SuperLearner work (SL-01…04) changes nothing in the NCE: it
cites the index directly.

## One thing not to say

That the ensemble or the SuperLearner is the model. Four meta-learner losses
were tested; none beats the zero-tuning index at these sample sizes. If the
"12-learner ensemble" language from the original proposal is still anywhere
in the draft, replace it with *"a parsimonious index over domain summaries,
which under our corrected evaluation outperformed more complex ensembles."*
