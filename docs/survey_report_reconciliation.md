# Survey-report reconciliation and critical review of the biomarker outcomes

*2026-09-15. Source material: the four national micronutrient survey final
reports (GMNS 2018, GMS 2017, SLMS 2013, MNS 2015-16) and the 21 publications
in `Downloads/Micronutrient survey reports/` that cite or re-analyse them.
Every report was read in full for design, laboratory methods, case
definitions, national and sub-national tables and appendices; every paper for
its definitions, sample sizes and headline numbers. The RA's index lists a
22nd paper (Malawi: "Malaria infection confounds inflammation-adjusted
micronutrient biomarker concentrations", 2020) that is not in the folder.*

## 1. Bottom line

1. **Two data-provenance problems** that no amount of modelling fixes:
   * The **Sierra Leone child file is the anaemic subset only** (n = 532 =
     the report's anaemia numerator; every child has Hb <= 10.9 g/dL; the
     survey assayed 654 children). All SL child results in the pipeline are
     estimates for anaemic children. **Request the full child file from
     GroundWork.**
   * The **Malawi RBP -> retinol calibration is not usable** (r = 0.36 in
     children, 0.19 in women; the survey's own re-analysis, Likoswe 2021, and
     the survey team, Williams 2021, both say so). The pipeline's vitamin A
     rule until today (`VITA_RULE=retinol_equiv`, VA-01) built on that line and
     yielded 0.4% child VAD for Malawi and 27% for Ghana, reordering the
     countries relative to every published estimate. **The default is now
     `rbp070`** (uniform BRINDA-adjusted RBP < 0.70); `retinol_equiv` is the
     sensitivity analysis.
2. **Five coding defects, now fixed** (section 5): IDA-labelled columns used as
   the iron-deficiency binary in three countries; a home-grown BRINDA that
   under-counted Malawi iron deficiency by half; a name-only Admin-2 join that
   double-counted 113 Malawian respondents; pregnant women and out-of-range
   children retained in every outcome; zinc sentinel values (-100, 0) counted
   as deficient. After the fixes every iron, zinc, folate and B12 national
   estimate reproduces its report (section 3).
3. **The analysis approach is sound where the reports are informative, and
   the reports point to three changes** (section 6): make uniform
   BRINDA-adjusted RBP < 0.70 the primary vitamin A outcome and keep the
   retinol-equivalent rule as a sensitivity analysis; add the Admin-1
   report tables as a standing reconciliation test; and add plasma selenium
   (already in the Malawi file) as a positive control with a known
   geospatial driver.

## 2. What the four surveys actually did (the parts that matter here)

| | Gambia GMNS 2018 | Ghana GMS 2017 | Sierra Leone SLMS 2013 | Malawi MNS 2015-16 |
|---|---|---|---|---|
| Fieldwork / season | 13 Mar - 4 May 2018, **dry** (malaria RDT < 1%; inflammation 28% children) | 27 Apr - 9 Jun 2017, onset of rains (malaria 21%; inflammation 46%) | 11 Nov - 2 Dec 2013, **end of rains** (malaria 53%; inflammation 73%) | Dec 2015 - Feb 2016, **rainy / lean** (malaria 28%; inflammation 57%) |
| Design, representativeness | 14 strata (8 LGAs x urban/rural), 71 EAs; representative by **LGA** | 3 agro-ecological belts x 30 EAs; representative by **belt** | urban/rural x 30 EAs; representative by **stratum** only | 105 of the 850 DHS clusters, 35 per region; representative by **region**; MDHS household weights |
| Blood | capillary (microtainer) | capillary in children (venous in MRDR subsample), venous in women | venous EDTA plasma; **shipment thawed** (sTfR dropped) | venous serum + trace-element tube |
| Assays | VitMin ELISA (ferritin, sTfR, RBP, CRP, AGP); retinol HPLC at VitMin, n = 14 | VitMin ELISA; retinol + MRDR HPLC Wisconsin, n = 300; folate/B12 Cobas e411 (Davis) | VitMin ELISA; retinol HPLC Davis, n = 33 women only; folate/B12 Cobas e411 | VitMin ELISA; retinol + MRDR HPLC INCAP; folate **microbiologic assay** (CDC); B12 Cobas 6000; zinc AES (CHORI) |
| Ferritin adjustment | **BRINDA** (both groups) | **Thurnham** (both) | **Thurnham** (internal correction factors) | **BRINDA** internal regression |
| RBP adjustment | BRINDA in children; none in women | Thurnham in children; none in women | Thurnham in both | **none** |
| VAD cut-off | RBP < 0.70 | RBP < 0.70 | RBP < 0.70 | **RBP < 0.46** (calibrated) |
| Women's iron cut-off | ferritin < 15 | < 15 | states < 15, **tabulated < 12** (see 3.3) | < 15 |
| Population | NPW; children 6-59 mo | NPW; 6-59 mo | NPW; 6-59 mo | NPW; 6-59 mo |

Two consequences run through everything below. (a) The surveys were fielded
in different seasons with inflammation prevalence from 28% to 76%, so any
residual (un-adjusted) inflammation effect on ferritin and RBP differs by
country and is confounded with the national level offsets the LOCO analysis
tries to transport. (b) The reports themselves mix two inflammation
adjustments and, for Ghana, publish both: child VAD is 28.9% unadjusted,
20.8% Thurnham, 13.1% BRINDA; child ID 16.5 / 21.5 / 29.5%; women ID 13.0 /
13.8 / 20.5% (GMS Appendix 14-15). Petry 2019 give the same for Gambia (ID
59.0% BRINDA vs 54.2% Thurnham). The choice of adjustment moves a national
prevalence by as much as the gap between countries.

## 3. Report vs pipeline, national prevalence (weighted %)

"Before" is the pipeline as found (`_targets_full`, protocol-v2 definitions);
"after" is with the fixes in section 5. Vitamin A rows show both rules; `rbp070`
became the default the same day (section 5).

| Country | Outcome | Report | Before | After (`retinol_equiv`, now the sensitivity rule) | After (`rbp070`, now the default) |
|---|---|---|---|---|---|
| Gambia | child VAD | 18.3 (BRINDA) | 17.3 | 16.7 | 16.7 |
| Gambia | women VAD | 1.8 | 1.7 | 1.6 | 1.6 |
| Gambia | child ID | **59.0** | 57.3 (binary col was IDA: 37.5) | **59.0** | |
| Gambia | women ID | **41.4** | 41.5 (IDA: 28.2) | **41.4** | |
| Ghana | child VAD | 20.8 Thurnham / 13.1 BRINDA | 27.2 | 27.2 | 14.7 |
| Ghana | women VAD | 1.5 | 2.9 | 2.9 | 1.7 |
| Ghana | child ID | 21.5 | 21.5 | 21.5 | |
| Ghana | women ID | **13.7** | 14.0 (IDA: 8.9) | **13.7** | |
| Ghana | folate def | 53.8 | 53.8 | 53.8 | |
| Ghana | B12 def | 6.9 | 6.9 | 6.9 | |
| Sierra Leone | child VAD | 17.4 (all children) | 6.5 | 6.5 (anaemic subset) | 12.0 |
| Sierra Leone | women VAD | 1.8 | 1.0 | 1.0 | 1.4 |
| Sierra Leone | child ID | 5.2 | 5.1 | 5.1 (anaemic subset) | |
| Sierra Leone | women ID | 8.3 (**< 12**; 12.8 at < 15) | 18.0 (BRINDA col) | **12.8** (Thurnham < 15) | |
| Sierra Leone | folate / B12 | 79.2 / 0.5 | 79.2 / 0.5 | 79.2 / 0.5 | |
| Malawi | child VAD | 3.6 (RBP < 0.46, unadj.); 0.4 adj. | 0.4 | 0.4 | 10.3 |
| Malawi | women VAD | 0.3 | 0.3 | 0.3 | 1.6 |
| Malawi | child ID | **21.7** | 21.5 (binary col: 10.3) | **21.7** | |
| Malawi | women ID | **15.1** | 15.9 (binary col: 8.1; pregnant incl.) | **15.1** | |
| Malawi | B12 def | 12.9 (< 150); Qi 2024: 11.8 (< 148) | 11.4 | **11.8** | |
| Malawi | child / women zinc | **60.4 / 62.5** | 61.1 / 64.1 | **60.4 / 62.5** | |

Iron, zinc, folate and B12 are now exact to the reports (the survey's own
flags are used where they exist). Every remaining gap is vitamin A.

A second check the reports allow and the pipeline did not use: Admin-1
tables. Mapping GADM districts to Gambia's eight LGAs, the pipeline's LGA
prevalences track GMNS Tables 18/21/32 with r = 0.99 (child VAD), 0.98
(women ID), 0.87 (child ID; the miss is Kanifing, whose boundary my district
map only approximates), MAE 1.6-3 pp. The individual -> area machinery is
doing what it should; the same test is available for Ghana's belts, Sierra
Leone's regions and Malawi's regions and should become a standing target
(section 6).

## 4. Findings in detail

### 4.1 Sierra Leone: the child file is the anaemic subset (data provenance)

`Sierra Leone_Child data_Davis-Berkeley.dta` has 532 rows; `cHb` runs 4.2-10.9
g/dL; `cAnemiaSev` = 27 severe / 315 moderate / 190 mild / 0 none;
`cAnemiaYN` is labelled "Not anemic" for all 532 (a mislabel). SLMS Table 23:
anaemia **532** of 710 children (76.3%); Table A8-10: 190 / 315 / 27. Ferritin
and RBP were assayed in 654 children (Table 26 footnote; Wirth 2018 "654
children and 774 women who provided blood samples"), so 168 non-anaemic
children with assays are absent. In the file `cIDA == cFeDefAdj` for every
child, as it must be when all are anaemic. The women's file is complete (945;
871 with Hb, 44.8% anaemic = report). Consequences: SL child VAD/ID
prevalences and their district pattern are for anaemic children; the
weights (`cStatWt`, 10 distinct values) were built for the full sample; n is
486 not 654. Not previously noticed (config and docs cite 532 as the sample).
`R/config.R` now carries the warning. Nothing in code can fix this.

### 4.2 The RBP -> retinol calibrations are not comparable (vitamin A)

Each survey calibrated the same VitMin RBP ELISA against HPLC retinol on a
subsample, in **different regression directions and with very different
sample sizes and fits**:

| Survey | n | Direction fitted | Line | R^2 | RBP at retinol 0.70 |
|---|---|---|---|---|---|
| Gambia | **14** | RBP on retinol | RBP = 0.978 ret + 0.0153 | 0.91 | 0.70 |
| Ghana | 300 | RBP on retinol | RBP = 1.1486 ret - 0.016 | 0.90 | 0.79 (the data's `cVADAdjSherryThurn`) |
| Sierra Leone | 33 women | retinol on RBP | ret = 0.196 + 0.788 RBP | 0.82 | 0.64 |
| Malawi PSC | 76 | retinol on RBP | ret = 0.379 + 0.755 RBP | **0.13** (refit; Likoswe: 0.20) | 0.43 |
| Malawi WRA | 95 | retinol on RBP | ret = 0.73 + 0.68 RBP (refit) | **0.03** | n/a |

`metadata/rbp_retinol_calibration.csv` inverts the Gambia and Ghana lines and
uses the SL and Malawi lines as published, so the four "retinol-equivalent"
scales are not the same estimand: an inverted RBP-on-retinol line has slope
1/b, a retinol-on-RBP line has slope r^2/b, and the two differ by the factor
r^2 (0.9 in Ghana, 0.13 in Malawi). Refitting Malawi PSC in the Ghana/Gambia
direction gives an RBP cut-off of 0.87, not 0.43; Deming regression 0.85.

Malawi is the decisive case. In the 74 children with both assays, **20.3%
have HPLC retinol < 0.70 and 23.0% have RBP < 0.70, but 1.4% have RBP < 0.46**.
The inverse regression predicts everyone toward the mean retinol (0.95) and
by construction produces near-zero deficiency; `retinol_equiv` then applies
BRINDA on top (0.4%). Likoswe et al. 2021 (Nutrients 13:849) re-analysed the
same subsample, found no inflammation adjustment rescues the fit, note the
SR:RBP > 1 ratio is biologically implausible and that the INCAP retinol had
no QC material, and recommend **BRINDA-adjusted RBP < 0.70 -> VAD = 10%** in
PSC (24% unadjusted) - exactly the pipeline's `rbp070` figure (10.2%).
Williams et al. 2021 (AJCN 113:854, the survey team) add that the 2009 Malawi
survey's calibrated cut-off was **0.78** with the same RBP laboratory and a
different retinol laboratory, that "analytic error if high retinol
concentrations were extrapolated above the standard curve" is a candidate
explanation, and that later UPLC retinol on n = 171 gave 1.8%
inflammation-adjusted VAD; MRDR found no depleted child. So Malawi's true
VAD is low, RBP < 0.70 (10%) overstates it, 0.46 (0.4%) understates it, and
RBP alone cannot resolve which - but the retinol-equivalent rule adds
calibration noise to precisely the between-country level differences the
LOCO analysis studies, and it does so with a 14-sample line for Gambia and a
0.13-R^2 line for Malawi. Under `rbp070` the countries read Gambia 17,
Ghana 15, SL 12 (anaemic subset), Malawi 10; under `retinol_equiv` Ghana 27,
Gambia 17, SL 6.5, Malawi 0.4. Ghana's MRDR (6.7% depleted) and Malawi's (0%)
both say RBP < 0.70 is a "low RBP" indicator loosely related to liver stores;
that is true of every country equally under `rbp070`, which is the point.

Recommendation (section 6): `VITA_RULE=rbp070` primary (same assay, same
laboratory, one adjustment, the estimand Likoswe recommend and the three
GroundWork reports used), `retinol_equiv` as a labelled sensitivity analysis,
and if a calibration is retained, fit all four in one direction (Deming) and
drop Malawi's.

### 4.3 Sierra Leone women's iron: the report used the child cut-off

The file has 97 women with Thurnham-adjusted ferritin < 15 (12.8% weighted)
but the report gives 65 (8.3%). With **< 12** the file gives 66 (8.4%), and
IDA (< 12 & Hb < 120) 53 (6.7%) against the report's 52 (6.1%). SLMS Table 3
and the Table 32 footnote say < 15, but Wirth 2016 Fig 1's caption says "iron
deficiency, ferritin <12 µg/L" for women. The published 8.3% is a report
error; the pipeline keeps the WHO < 15. Anaemia (390/871), VAD (18), folate
(608) and B12 (4) reproduce the report exactly, so the file is otherwise
consistent with it.

### 4.4 Iron-deficiency binaries were iron-deficiency ANAEMIA in three countries

`R/config.R` pointed `child_iron`/`women_iron` at `gw_cIDA_Brinda` /
`gw_wIDA_Brinda` (Gambia; labels "Iron Deficiency Anemia with ID ..."; 42% /
31% vs the report's ID 59.0 / 41.4), `gw_wIDA_Thurn` (Ghana women; 7.7% vs ID
13.7) and `gw_cIDA` / `gw_wIDA` (Sierra Leone; "Both iron deficiency and
anemia"). The continuous targets were iron deficiency, so binary and level
described different conditions. `resolve_uniform_outcome()` masked this on
the area path by re-deriving the binary from the continuous column, but the
row filter and every direct reader of `oc$binary` (`national_estimates.R`,
the individual-level SL, `corrected/p1`, `cluster_aggregation.R`, `mrp.R`,
`leakage_report.R`) used IDA. Fixed: the survey's own ID flags are configured
(`gw_cID_Brinda`, `gw_wID_Brinda`, `gw_wIDAdjThurn`, `gw_cFeDefAdj`,
`gw_wFeDefAdjThurn`), and SL women's continuous is now the survey's Thurnham
ferritin (`gw_wFerrAdjThurn`, which exists; the old note that "gw_wFerrAdj
has 0 values for women" was looking at the child file's mothers' column).

### 4.5 Malawi: home-grown BRINDA under-counted iron deficiency by half

`src/malawi/1_DHS_mn_data.R` regressed log ferritin on log CRP + log AGP and
took `exp(residual + intercept)` for everyone, i.e. re-centred every child to
CRP = 1 mg/L and AGP = 1 g/L. For the majority with CRP below 1 mg/L that
*raised* ferritin. Result: `iron_def` 107/1102 = 9.7% (women 74/752) against
the survey's own BRINDA flag `sf_c1` (= `sf_reg` < 12/15 exactly) of 222 =
20.1% (women 131; report 21.7 / 15.1). The configured continuous was already
`sf_reg`, so the two targets disagreed within Malawi as well. Fixed: binaries
point at `sf_c1`; the script's function now implements the reference-decile
BRINDA (only observations above the 10th-percentile reference are adjusted,
coefficients clamped by the direction of the acute-phase response) and uses
`sf_reg` when present; the SAC ferritin cut-off is corrected to 15.

### 4.6 Malawi: duplicated Admin-2 names doubled 113 respondents

GADM level 2 for Malawi repeats four TA names across districts (TA Lundu:
Blantyre/Chikwawa; TA Malemia: Nsanje/Zomba; TA Ngabu: Chikwawa/Nsanje; TA
Pemba: Dedza/Salima). The individual-level GEE merge in `_targets.R` joined
on the Admin-2 name alone, so `merged_malawi` had 3212 rows for 3099
respondents: every person in the surveyed TA Lundu, TA Malemia and TA Pemba
appeared twice, once with the other district's covariates, and
`targets_v2.csv` carried n_raw = 46 for TA Malemia (23 children). The
`[merge_gee] Row count changed` warning fired and was ignored.
`R/admin2_key_hygiene.R` (2026-08-28) had fixed the area-level joins but not
this one. Fixed: the merge uses `admin2_join_by()` (Admin1 + Admin2), the
covariate table is de-duplicated on that key, and a changed row count is now
an error.

### 4.7 Population definitions: pregnant women and children outside 6-59 months

All four reports tabulate biomarkers for non-pregnant women and children
6-59 months. The pipeline filtered on the child/women flag only. Gambia's
women's file has 158 self-reported pregnant women, 32 with ferritin/RBP (15
iron deficient); Malawi's 34, 31 with RBP/B12/zinc (the survey's own `sf_reg`
is already NA for them); Ghana's 153, only 5 with assays. The Gambia child
file has 21 children aged 60-64 months and 5 aged < 6 months with assays. With
both restrictions the Gambia ID numerators reproduce the report exactly (632
women; children 59.0%) and the VAD numerator drops from 222 to the report's
218. Fixed via one shared `outcome_population_mask()` used by
`build_outcome_dataset()`, the cluster track, the LOCO pooling and
`corrected/p12`; `preg_col`, `child_age_col`, `child_age_range` are declared
per country in `R/config.R`.

### 4.8 Zinc

The survey's `low_zn` (IZiNCG cut-offs by age, draw time and fasting) and
the local `zinc_def` agreed for 1085/1086 children; the exception was a
`zn_gdl` of -100 counted as deficient, and one woman had 0. Sentinels are now
NA and the binaries use `low_zn` (report: 60.4 / 62.5%). Two of the Malawi
papers (Likoswe 2020; Gebremedhin 2020) show inflammation adjustment lowers
the child prevalence by 2-10 pp; the report and the pipeline leave zinc
unadjusted, which is defensible but should be stated. Likoswe also drop
values > 125 µg/dL as contamination; the file has such values (max 378) and
they enter the LEVEL target as -log(zinc). Consider capping for the level.

### 4.9 Folate is not harmonisable across assays as configured

Ghana and Sierra Leone measured serum folate on the Roche Cobas e411
immunoassay; Malawi used the CDC microbiologic assay, for which the WHO < 10
nmol/L cut-off was set (the Malawi report itself uses < 6.8 and < 14 and RBC
folate < 748). Immunoassays read systematically lower than the MBA, so the
Ghana 54% / SL 79% / Malawi 19% contrast is part assay, part fortification.
`women_folate` should be flagged as within-country only, or a published
cross-calibration applied. B12 was Cobas in all three (comparable); Malawi's
`vitb12` is pg/mL (the 1.355 conversion is confirmed by Qi 2024: 11.8% <
148 pmol/L, pipeline 11.8%).

### 4.10 Smaller items

* Gambia and Ghana children gave **capillary** blood; `metadata/assay_sources.csv`
  said venous (corrected).
* Sierra Leone's Thurnham factors are internal correction factors, not the
  published ones (children 0.38 / 0.49 / 0.74; women 0.70 / 0.58 / **1.02**
  for AGP-only, which raises ferritin). This is the survey's choice; note it.
* The Ghana report's Table 27 footnote says "BRINDA" while the methods say
  Thurnham (a copy error; the numbers are Thurnham). GroundWork's own later
  paper on the same data (Donkor 2021) uses BRINDA for both ferritin and RBP.
* Malawi's `strata_col = NULL` is right for the point estimates (the report
  used household weights only) but the design has region x residence strata.
* SL child `cRBPAdj` is adjusted only for children >= 6 months (label); all
  children in the file are >= 6 months.

## 5. What was changed in the code (2026-09-15)

| File | Change |
|---|---|
| `_targets.R` (merged_<country>) | GEE covariates joined on Admin1 + Admin2 via `admin2_join_by()`; covariate table de-duplicated on the key; row-count change is an error. |
| `R/config.R` | Iron binaries -> survey ID flags (Gambia `gw_cID_Brinda`/`gw_wID_Brinda`; Ghana `gw_wIDAdjThurn`; SL `gw_cFeDefAdj`/`gw_wFeDefAdjThurn`; Malawi `sf_c1`); SL women's continuous -> `gw_wFerrAdjThurn`; zinc binaries -> `low_zn`; `preg_col`, `child_age_col`, `child_age_range` per country; SL anaemic-subset warning. |
| `R/data_prep.R` | `outcome_population_mask()` (flag + non-pregnant + 6-59 months), used by `build_outcome_dataset()`; Malawi zinc sentinels -> NA; `sf_c1` logged. |
| `R/cluster_aggregation.R`, `R/cluster_mbg.R`, `R/transportability.R`, `R/corrected/p12_distributional.R` | use `outcome_population_mask()`. |
| `src/malawi/1_DHS_mn_data.R` | reference-decile BRINDA with direction clamp; `sf_reg` preferred; SAC cut-off 15; zinc sentinels. (Raw `MW_*.DTA` are not in the repo; the fix applies on the next re-run.) |
| `metadata/assay_sources.csv`, `metadata/adjustment_inventory.csv` | capillary specimens; inventory regenerated. |

Every cached target downstream of `merged_*` and `outcome_data_*` is now
outdated; `targets::tar_make()` (full mode for the headline tables) and the
protocol-v2 scripts 01-02 need to be re-run, then the dashboard bundles.

| `R/brinda_adjustment.R` | (PI decision, same day) `VITA_RULE` default switched from `retinol_equiv` to **`rbp070`**; `VITA_RULE=retinol_equiv` runs the sensitivity analysis. `docs/manuscript_mcn_v2.qmd` methods paragraph updated to match. |

Still open (analytical decisions for the PI):
* Whether to cap serum zinc at 125 µg/dL for the level target.
* Whether to keep `women_folate` in the cross-country (LOCO) analyses.

## 6. Should the analysis approach change?

**Keep.** The ecological, area-level framing; survey-weighted district
targets with the design effect estimated where PSUs are plentiful; uniform
inflammation adjustment re-derived from the raw assays (the reports show
why: the adjustment alone moves Ghana child VAD 13-29%); one shared outcome
resolver. The Gambia LGA reconciliation (r = 0.98-0.99) shows the
aggregation is right.

**Change.**

1. *Vitamin A estimand* - `rbp070` primary, `retinol_equiv` sensitivity, per
   4.2. Label the outcome "low RBP (< 0.70 µmol/L, BRINDA-adjusted)" rather
   than VAD in the manuscript; note Ghana MRDR 6.7% and Malawi MRDR 0%.
2. *Report reconciliation as a target.* Encode the published national and
   Admin-1 tables (`metadata/mn surveys/`) as a small CSV and add a target
   that compares `svy_admin1_*` / national estimates to them under the
   survey's own definitions, failing on > 3 pp. It would have caught 4.1,
   4.4, 4.5 and 4.6 immediately. The Gambia LGA map used here is in
   `scratchpad/lga_check2.R`; Ghana belts, SL regions and Malawi regions are
   direct Admin-1 matches.
3. *Selenium as a positive control.* `sel` (plasma Se, µg/L) is in the Malawi
   file for 990 children / 802 women; district alone explains 42% of its
   log-variance, and its geology/soil-pH driver is established (Phiri 2019;
   Phiri 2020 in this folder; Gashu 2021). If the covariate stack recovers
   Se's spatial structure it validates the machinery; if it does not, that is
   diagnostic. (A parallel session added `child_selenium`/`women_selenium`/
   `women_iodine` outcomes to `R/config.R` today, with `sel_def` = plasma Se
   < 84.6 µg/L: 86% of children, in line with Phiri 2019.)
4. *Sierra Leone children.* Obtain the full file; until then present SL child
   results as anaemic-subgroup estimates and exclude them from LOCO
   vitamin A/iron summaries, or drop the cell.
5. *Season and inflammation as design covariates.* The four surveys span dry
   to peak-malaria seasons. The country-level residual inflammation effect is
   confounded with the level offset; carry fieldwork-window malaria/inflammation
   prevalence as a country-level covariate in the transport model (the
   pipeline has the windows in `metadata/survey_years.csv`) and report
   sensitivity to CRP-only vs CRP+AGP adjustment, which the Malawi papers show
   matters when malaria is high.
6. *Heterogeneity is the finding, not the noise.* Petry 2021 (Ghana) show the
   anaemia-risk-factor relationships differ by belt; the Malawi papers show
   malaria drives ferritin, sTfR, RBP and zinc. This argues for the current
   within-country-centred transport model and against pooling raw levels; it
   also argues for reporting the zone x covariate interaction rather than
   a single national coefficient.

## 7. The 21 publications in one line each

*Malawi (12 in folder).* **Likoswe 2021** (Nutrients 13:849): RBP-retinol r ~ 0.2, calibration unusable, recommends BRINDA RBP < 0.70 (10%). **Williams 2021** (AJCN): VAD near-eliminated by triangulation (UPLC retinol, MRDR, retinyl esters 18% > 5%); 2009 cut-off was 0.78. **Likoswe 2020** (Nutrients 12:1563): zinc 62 -> 59 (ICF) -> 52% (BRINDA); drops > 125 µg/dL. **Gebremedhin 2020** (Nutrition): zinc adjustment changes prevalence 1-3 pp. **Ntenda 2025** (Int Health): PSC ID 19.6% (BRINDA), FID 50%; rural and South higher. **Ntenda 2022** (Malaria J): SAC malaria -> anaemia and sTfR, not ferritin. **Ntenda 2026** (Int Health): SAC inflammation lowers RBP/zinc, raises ferritin; Central/South lower RBP. **McGann 2018** (Blood Adv): alpha-thal 43%, G6PD 20% of boys, sickle trait 9%; ID 21.5%. **Ntenda 2019** (IDP): clinical malaria aOR 4.6 for anaemia. **Rhodes 2020** (J Nutr): NPW only; VAD = RBP < 0.46. **Qi 2024** (Birth Defects Res): NPW n = 778; serum folate median 18 nmol/L (MBA), B12 < 148 = 11.8% - confirms units. **Phiri 2020** (Environ Int): urine Se tracks plasma Se between clusters. *Missing:* the 2020 "malaria confounds inflammation-adjusted biomarkers" paper.

*Gambia (3).* **Petry 2019** (Nutrients 11:2275): the survey's primary paper; BRINDA; Thurnham would give ID 54.2 / 36.8; retinol HPLC n = 14 "confirmed" 0.70; women's VAI cut-off 1.05. **Akindutire 2025** (Front Public Health): uses **DHS 2019-20**, not the GMNS. **AJOL 2020**: uses **MICS 2018**; irrelevant.

*Sierra Leone (3).* **Wirth 2016** (PLoS One): primary paper; 654 children with assays; women's ID caption "< 12". **Rohner 2016** (Nutrients 8:74): iodine; NPW 817. **Wirth 2018** (BMC Res Notes): sickle/thal in a non-random 388-child subsample; confirms 654.

*Ghana (3).* **Petry 2021** (MCN): Thurnham; risk factors differ by belt. **Donkor 2021** (Life): 6-23 mo, **BRINDA** ferritin and RBP (ID 45%, VAD 10%). **Christian 2022** (Nutrients 14:1427): NPW n = 1063; B12 < 148.
