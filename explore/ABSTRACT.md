# Exploratory feature engineering and modelling for sub-national micronutrient prediction: ten probes, one lead

Draft abstract, 2026-09-28. Everything below is hypothesis-generating and
post hoc on four countries; nothing is a production claim.

---

**Background.** Sub-national prediction of micronutrient deficiency from
remotely sensed and administrative proxies is limited by a weak predictor
signal at Admin-2: across 24 country × outcome cells (14–87 districts each,
383 predictors), a zero-tuning domain index beats every tuned learner, and
within a surveyed country covariates add nothing over a spatial smoother. We
asked whether feature engineering or estimators drawn from quantitative
genetics, chemometrics and environmental epidemiology could find signal that
the existing pipeline misses, for specific nutrients, specific countries, or
broadly.

**Methods.** Ten probes were run in an isolated sandbox against a scorer that
reuses the parent project's own cell construction, within-country rank
normalisation, domain representation, fold seeds and metrics; the scorer
reproduced the published benchmark table to a mean absolute difference of
0.0009 over 376 cell × arm comparisons, with leave-one-country-out transport
exact. Probes covered (i) multi-kernel REML-BLUP and AlphaEarth satellite
embeddings used as kernels rather than principal components; (ii) removal of
unwanted variation (RUV/SVA) for the cross-survey biomarker offset;
(iii) empirical-Bayes moderation of index weights pooled across cells;
(iv) multi-trait BLUP through a shared latent factor; (v) nutrient-specific
mechanistic features derived from the 42-crop MapSPAM grid and a food
composition table (phytate:zinc molar ratio, provitamin-A carotenoid density,
oil-palm share); (vi) seasonal phase and a 0–35 month lag stack extracted from
Earth Engine at 323 survey-cluster buffers; and (vii) seven n ≪ p estimators
(PLS, PCR, supervised PCA, MCP, stability selection, CV-ridge, standardised
index). Comparisons were paired on identical folds and tested over country
blocks.

**Results.** One arm survived: **multi-trait BLUP**, fitting a survey's
biomarkers jointly rather than one cell at a time, improved district in-fill by
0.038 Spearman over single-trait fitting on identical kernels (6 of 6
country × target blocks, p = 0.031), with the gain concentrated where
single-trait fitting was weakest (weak tercile +0.086, strong +0.018;
ρ(baseline, gain) = −0.388, p = 0.020) and rescuing cells currently scored as
failures. Transport results converged on climate and soil from three
independent routes — the existing domain ablation, a REML-shrunk kernel
statistically indistinguishable from the pre-registered climate+soil index
(−0.0003 over 44 cells), and a cross-cell empirical-Bayes screen in which
climate, satellite embedding, soil and greenness supplied 111 of 173
predictors surviving FDR < 0.05 — strengthening that recipe as a claim about
predictors rather than machinery. Eight directions were closed, several
structurally: predictor-side batch correction is a mathematical no-op, because
within-country rank normalisation forces every linear projection to have zero
country separation (F = 0.0, η² = 0.000); the satellite embedding's failure is
not a representation artefact (kernel loses to PCs in 0 of 6 blocks);
nutrient-matched mechanistic features did no better than deliberately
mismatched ones (6 of 18 cells), despite validating against agronomy; seasonal
phase was a coin flip (12 of 24); and the previous-growing-season hypothesis
was contradicted rather than merely unsupported, the lag-association profile
being flat and peaking at 31 months, with weather three years before the blood
draw as associated as the last harvest. No n ≪ p estimator beat the index, and
the two performing explicit variable selection were the worst.

**Conclusions.** Every attempt to extract more from the predictor side failed;
the only arm that helped borrowed strength across outcomes. Combined with a
nested test showing covariates add nothing over a spatial smoother even under
identical machinery, and with the flat lag profile, the picture is that the
predictable Admin-2 signal is a smooth agro-ecological surface close to
saturated by climate and soil, and that remaining headroom lies in how outcomes
are pooled and how much survey noise sits beneath them rather than in more or
better covariates. Two methodological cautions generalise: medians across
heterogeneous cells disagreed in direction with paired comparisons twice, and
sign tests must block by country, since cells within a country share districts.

**Limitations.** Four countries, post hoc selection, no pre-registration; the
single positive rests on six blocks with one country excluded for having too
few districts to fold. A separate probe suggesting that transported *levels*
might be recoverable for vitamin B12 was retracted on checking against an
existing pre-registered analysis, which it reproduced (11 of 22 units) rather
than refined, and whose proposed correction it already implemented.
