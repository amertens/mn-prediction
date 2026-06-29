# Project status: micronutrient deficiency prediction (June 2026)

A snapshot of where the analysis and dashboard stand, what is outstanding, and
the key scientific insights from the current results. Companion to the
methodological roadmap (`REVIEW_AND_ROADMAP_2026-06.md`), the manuscript
(`manuscript_mcn.qmd` + `manuscript_mcn_appendix.qmd`), and the exploratory
findings (`EXPLORATORY_ANALYSES.md`).

## Analysis: where it stands

- **One pipeline.** The production pipeline and the parallel P1-P8 corrected
  pipeline have been consolidated into a single `_targets.R` (498 production
  targets plus 241 folded-in corrected targets, 739 total), built into
  `_targets_full/`, which is the analysis of record. The standalone
  `_targets_corrected.R` is archived.
- **Corrected methods are canonical.** The leakage-free corrected methods
  (in-fold preprocessing, cluster and spatial-block cross-validation, out-of-fold
  calibration, sampling-error-aware admin-2 error, decision-value scoring,
  split-conformal and design-based intervals, out-of-support trust flags,
  partial-pooling small-area estimation) are the headline; the original pipeline
  is retained as a sensitivity comparison.
- **Manuscript.** A high-level Methods section in the main text plus a detailed
  supplementary appendix (data sources, cleaning and linkage, outcome
  harmonization, modeling, accuracy assessment, sensitivity analyses,
  reproducibility). Both rendered to Word.
- **Supporting artifacts.** A multi-outcome "simplified subset" (data now
  committed) with Andy Kim's project plan; exploratory findings preserved in
  `docs/exploratory_*.md`; repo cleaned (sandbox workspaces retired, stale
  `_targets/` and `_targets_fast_backup/` stores removed). All committed to main.

## Dashboard: where it stands

- **Internal app is live and current:** https://amertens.shinyapps.io/micronutrient-burden/
  It includes the area-level SAE map (HAL, BYM2, Fay-Herriot layers), the Sierra
  Leone Admin-3 (chiefdom) layer, the populated District Factors panel
  (model-agnostic SHAP), corrected confidence-interval and calibration language,
  the misclassification-risk layer, the Start-here guide, and the Decision-value
  and Scenarios tabs. It carries an "In development" banner.
- **Public app** (`app_public.R`) is built but not deployed, blocked on the
  shinyapps.io five-app account limit. The internal URL is shareable meanwhile.

## Outstanding to-dos

Methodological (highest value first; see `REVIEW_AND_ROADMAP_2026-06.md`):

1. Decision-focused ensembling / top-N targeting: tune the ensemble to correctly
   rank the highest-burden districts.
2. National-anchor calibration: shift transported predictions so their mean
   matches a known national prevalence, fixing the level bias.
3. Sequential multi-outcome modeling; proper binomial loss for prevalence;
   aggregate uncertainty plus a value-of-information ("where to sample next") map.
4. Anemia validation track: rerun the method on DHS haemoglobin (large samples,
   nearly every sub-Saharan country) to validate the methodology on a trusted
   outcome and expand the cross-country panel.

Project and operations:

- Andy Kim to start on the simplified subset (top-N ranking plus national-anchor
  calibration) toward the September Ghana presentation.
- GBD comparison: source real GBD Results Tool data (RA task) and refine the
  predictor-family covariate-comparison table.
- Public app: deploy once a shinyapps slot is freed or the tier upgraded; address
  reconnect/stability.
- Code hygiene: a reviewed deduplication pass on the legacy `src/` tree; clarify
  whether `national_prediction/` is active or superseded.
- Manuscript: fill the repository-URL placeholder and co-author list; decide
  whether to push the corrected metrics fully into the headline result tables.
- Whole-region sub-Saharan expansion: deferred until the core models are stronger.

## Key scientific insights from current results

1. **The binding constraint is the data, not the model.** Four surveys, ~206
   admin-2 areas, individually weak proxies, and biomarker comparability that is
   itself disputed in the field.
2. **Effective sample size is the number of areas (14 to 87), not individuals.**
   The predictors are area-resolved, so almost all of their variance is
   between-area; this drives the weak and unstable transfer.
3. **Honest evaluation changes the story.** Optimistic and in-sample evaluation
   overstated performance. Under cluster and spatial-block cross-validation,
   within-country discrimination falls toward chance (for example Gambia child
   iron AUC 0.57 to 0.48), and the directly-fit area ensemble's apparent
   excellence was an in-sample artifact.
4. **Within-country ranking is modest; national estimates are good.** The
   individual-level SuperLearner aggregated to admin-2 ranks districts at a median
   Pearson r near 0.53 (the area ensemble near 0.12 under honest evaluation), but
   national prevalence is accurate to within about 5 percentage points. National
   estimation is the product that already works.
5. **Cross-country transport fails mainly on level, not ranking.** Predictions
   are pulled toward the training-country mean. Iron transports best; vitamin A
   essentially does not; Malawi, the only non-West-African country, anti-correlates
   when held out. National-anchor calibration targets exactly this failure.
6. **The ensemble does not beat simple small-area methods** on most
   leave-one-country-out splits. The best transportable model is the simplest: a
   thin-plate spatial spline on district centroids plus a few soil-micronutrient
   features (mean leave-one-country-out r near 0.33), beating all 152
   earth-observation covariates. Soil micronutrients are the closest covariate to
   a causal pathway (soil to crop to intake to biomarker).
7. **Calibration and uncertainty are now honest.** Many binary models had
   negative Brier skill when raw, repaired by out-of-fold recalibration;
   split-conformal intervals cover roughly 87 to 100 percent against a 90 percent
   target; honest area intervals remain wide (Fay-Herriot tightens them 19 to
   33 percent).
8. **Strategic reframe.** Because ranking survives where absolute level does not,
   the fundable deliverables are (a) correctly flagging the highest-burden
   districts and (b) value-of-information sampling, rather than precise absolute
   prevalence.
