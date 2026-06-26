# Biomarker survey acquisition guide

**Purpose:** Concrete acquisition instructions for the candidate micronutrient
*biomarker* surveys to add to the current 4-country panel (Gambia GMNS 2018,
Ghana GMS 2017, Sierra Leone SLMS 2013, Malawi MNS 2015–16). Every survey below
has **serum/blood biomarkers** — ferritin (iron), RBP or serum retinol
(vitamin A), serum zinc, folate, B12, and/or urinary iodine — **not** hemoglobin/
anemia alone.

**Harmonization shortcut:** any survey marked **[BRINDA]** was pooled by the
BRINDA project with standardized ferritin + retinol/RBP + CRP/AGP inflammation
adjustment — request via BRINDA to get comparable biomarker definitions out of
the box. BRINDA overview: https://pmc.ncbi.nlm.nih.gov/articles/PMC4785469/

**Subnational caveat:** all candidates resolve to **admin-1 (region/province)**
except **Pakistan NNS 2018 (district)**. Adding countries broadens cross-country
range but does not fix the admin-2 resolution gap.

---

## Tier 1 — priority additions (all-African, full panel)

### 1. Nigeria — NFCMS 2021 (National Food Consumption & Micronutrient Survey)
- **Country-year:** Nigeria, 2021
- **Biomarkers:** serum ferritin, sTfR, serum retinol + RBP, serum zinc, folate, vitamin B12, urinary iodine, CRP, AGP
- **Populations:** children 6–59 mo; women of reproductive age; adolescent girls
- **Subnational:** national + geopolitical zone / selected states
- **Access:** request microdata via GHDx; final report public
  - GHDx record: https://ghdx.healthdata.org/record/nigeria-food-consumption-and-micronutrient-survey-2021
  - Final report (PDF): https://www.iita.org/wp-content/uploads/2024/05/NFCMS-2021-Final-Report.pdf
- **Notes:** most recent + largest full panel; top priority.

### 2. Ethiopia — ENMS 2015 (Ethiopian National Micronutrient Survey, EPHI)
- **Country-year:** Ethiopia, 2015
- **Biomarkers:** ferritin, sTfR, RBP, serum zinc, folate, vitamin B12, urinary iodine, CRP, AGP
- **Populations:** preschool children, school-age children, women of reproductive age
- **Subnational:** national; all 9 regions + 2 city administrations (admin-1)
- **Access:** request to Ethiopian Public Health Institute (EPHI); GHDx record exists
  - GHDx record: https://ghdx.healthdata.org/record/ethiopia-national-micronutrient-survey-2015
  - EPHI: https://ephi.gov.et/ (Food Science & Nutrition Research Directorate)

### 3. Kenya — KNMS 2011 (Kenya National Micronutrient Survey)
- **Country-year:** Kenya, 2011
- **Biomarkers:** ferritin, RBP, serum zinc, folate, CRP/AGP (urinary iodine in some groups)
- **Populations:** preschool children, school-age children, WRA, pregnant women, men 15–54
- **Subnational:** national; 8 former provinces (not powered for the 47 counties)
- **Access:** report public; microdata via KNBS NADA catalog
  - Report (PDF): https://www.gainhealth.org/sites/default/files/publications/documents/kenya-national-micronutrient-survey-report-2011.pdf
  - KNBS NADA microdata catalog: https://statistics.knbs.or.ke/nada/index.php/catalog/4
- **Notes:** national summary already hand-extracted into `results/external/knms2011_national_summary.csv`. **[BRINDA]** (Kenya 2007/2010 child surveys are in BRINDA).

### 4. Cameroon — NMS 2009 (HarvestPlus / IRD)
- **Country-year:** Cameroon, 2009
- **Biomarkers:** ferritin, sTfR, serum retinol, plasma zinc, folate, vitamin B12
- **Populations:** children + women
- **Subnational:** 3 macro-regions
- **Access:** direct request to PI **Reina Engle-Stone (UC Davis)**; well-documented in published papers. **[BRINDA]**
  - Contact via UC Davis Dept. of Nutrition: https://nutrition.ucdavis.edu/

### 5. Tanzania — TDHS 2010 micronutrient subsample (NUT5)
- **Country-year:** Tanzania, 2010
- **Biomarkers:** vitamin A (RBP), iron (serum ferritin), urinary iodine
- **Populations:** children 6–59 mo + women of reproductive age (national subsample)
- **Subnational:** national subsample (not full admin-1)
- **Access:** **free** via the DHS Program (registration + project request); biomarkers are in a separate micronutrient report/recode
  - NUT5 report: https://dhsprogram.com/pubs/pdf/NUT5/NUT5.pdf
  - DHS data request: https://dhsprogram.com/data/
- **Notes:** easiest access (DHS) but a narrower panel (no zinc/folate/B12).

---

## DHS surveys with genuine biomarker subsamples (the short list)

DHS biomarker rounds beyond Hb/anemia are rare — two flavors: a vitamin-A
(retinol → RBP-on-dried-blood-spot) wave, and a few full-panel micronutrient
surveys co-designed on the DHS frame.

| Country-year | Biomarkers | Population | Access + URL |
|---|---|---|---|
| **Tanzania 2010 (NUT5)** | RBP, ferritin, iodine | children 6–59mo, WRA | DHS Program — https://dhsprogram.com/pubs/pdf/NUT5/NUT5.pdf |
| **Cambodia 2014 (CDHS-linked MNS)** | ferritin, sTfR, folate, B12, RBP, CRP/AGP | 792 children + 720 women | DHS Program; method paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4924189/ |
| **Malawi 2015–16 (MDHS / MNS, FR319)** | ferritin, sTfR, RBP, zinc, B12, folate, iodine | PSC, SAC, WRA, men | *already in panel* — https://dhsprogram.com/pubs/pdf/FR319/FR319.m.final.pdf |
| **Uganda 2000–01 (FR128)** | serum retinol (HPLC) | children 6–59mo, WRA | DHS Program — https://www.dhsprogram.com/pubs/pdf/FR128/FR128.pdf |
| **Uganda 2006 / 2011 / 2016** | RBP (DBS) — vitamin A only | children, WRA | DHS Program (vitamin-A subsample) |
| **Zimbabwe DHS (~2010–11)** | RBP (DBS) — vitamin A only | children, WRA | DHS Program *(exact round to confirm)* |
| **Rwanda RDHS 2019–20** | micronutrient module (confirm biomarker list) | — | World Bank microdata https://microdata.worldbank.org/index.php/catalog/4065 |

**DHS access mechanics:** free account at https://dhsprogram.com → per-project
dataset request. Biomarker values are usually in a **separate biomarker recode
file**, not the standard recode. Reference: https://dhsprogram.com/Methodology/Biomarkers.cfm.
Malawi MNS lab inquiries: immpact@cdc.gov.

---

## Non-African comparators (full panel; useful for harmonization range)

| Country-year | Biomarkers | Subnational | Access + URL |
|---|---|---|---|
| **Pakistan NNS 2018–19** | ferritin, retinol, vit D, B12, folate, zinc, iodine, CRP/AGP | **district** (rare depth) | GHDx https://ghdx.healthdata.org/record/pakistan-national-nutrition-survey-2018-2019 ; HDX https://data.humdata.org/dataset/pakistan-national-nutrition-survey |
| **Bangladesh NMS 2011–12** | ferritin, sTfR, retinol/RBP, zinc, folate, B12, iodine | national; urban/rural/slum | Report https://www.gainhealth.org/resources/reports-and-publications/national-micronutrients-status-survey-2011-12-final-report ; GHDx https://ghdx.healthdata.org/record/bangladesh-national-micronutrients-status-survey-2011 **[BRINDA]** |
| **Côte d'Ivoire NMS 2007** | ferritin, retinol, folate, B12 | 9 eco-regions | request (GAIN/PI) **[BRINDA]** — relevant: CIV is our in-flight OOS test country |
| **South Africa SANHANES-1 2011–12** | ferritin, sTfR, RBP | provincial (9) | HSRC DUA https://hsrc.ac.za/news/research-outputs/curated-data-the-south-african-national-health-and-nutrition-examination-survey-sanhanes-1/ |
| **Zambia FCMSS 2020** | ferritin, retinol, zinc | provincial | request (NFNC / Intake) |

---

## MICS — do NOT use for biomarker outcomes

Standard MICS rounds (incl. **MICS6**) collect **no serum micronutrient
biomarkers** — no ferritin/RBP/retinol/zinc/folate/B12, and no urinary iodine.
What MICS has is a **household salt-iodine commodity test** (rapid kit, not a
biomarker), vitamin A **supplementation-coverage** questions (not serum
retinol), anthropometry, and—only in some early rounds—Hb/anemia. Use MICS for
**household covariates only**, never as a deficiency outcome.

- MICS FAQ: https://mics.unicef.org/faq
- MICS microdata (for covariates): https://mics.unicef.org/surveys

---

## Access-route cheat sheet
| Route | Surveys | How |
|---|---|---|
| **DHS Program** (free) | Tanzania 2010, Cambodia 2014, Uganda, Malawi MNS | account → project request → biomarker recode file |
| **GHDx request** | Nigeria 2021, Ethiopia 2015, Pakistan 2018, Bangladesh 2011 | record page → "request access" |
| **National stats office** | Kenya (KNBS NADA) | NADA catalog registration |
| **Institute / PI request** | Ethiopia (EPHI), Cameroon (Engle-Stone, UC Davis), Zambia (NFNC) | email PI/institute |
| **DUA** | South Africa SANHANES (HSRC) | data-use agreement |
| **BRINDA** | Kenya, Cameroon, Côte d'Ivoire, Liberia, Bangladesh, Malawi | request pooled, inflation-harmonized biomarkers |
