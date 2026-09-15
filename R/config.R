# =============================================================================
# R/config.R
#
# Pipeline configuration: country registries, outcome definitions, domains.
# Returns a list-of-lists; each country has its own config block.
# Adding a new country = adding a new entry to get_country_configs().
# =============================================================================

#' Return all country configurations
#'
#' Each entry contains: data_path, gadm_code, survey design columns,
#' outcome definitions, domain prefixes, and pipeline parameters.
#' @return named list of country config lists
get_country_configs <- function() {
  list(

    # ── Gambia (GMNS 2018) ────────────────────────────────────────────────
    # Gambia National Micronutrient Survey. Fieldwork: 13 Mar – 4 May 2018.
    # DHS year: 2019 (Gambia DHS 2019-20). No 2018 DHS round exists for
    # Gambia; 2019 is the closest available (+1 year post-survey). The prior
    # round (2013) has a 5-year gap and is significantly less representative.
    Gambia = list(
      country     = "Gambia",
      gadm_code   = "GMB",
      survey_year = 2018L,
      dhs_year    = 2019L,    # closest available; no 2018 DHS round exists
      data_path   = here::here("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"),
      raster_dir  = here::here("data", "Gambia_GEE_rasters"),

      # Survey design columns
      cluster_id    = "gw_cnum",
      admin1_col    = "Admin1",
      admin2_col    = "Admin2",
      psu_col       = "gw_cnum",
      weight_col    = "gw_svy_weight",
      strata_col    = "gw_strata",
      child_flag    = "gw_child_flag",
      # 2026-09-15 (survey-report reconciliation): the GMNS report and Petry
      # et al. 2019 estimate every biomarker in NON-PREGNANT women and in
      # children 6-59 months. The women's file carries 158 self-reported
      # pregnant women, 32 of them with ferritin/RBP, and the child file 21
      # children aged 60-64 months (aged out since the MICS listing) and 5
      # aged <6 months with assays. build_outcome_dataset() now applies both
      # restrictions; with them the ID count among women reproduces the
      # report's numerator (632) exactly.
      preg_col      = "gw_wPreg",              # 1 = pregnant by self-report
      child_age_col = "gw_cAgeMonths_GMNS",    # completed months; keep 6-59
      child_age_range = c(6, 59),

      # Outcomes (each is population x micronutrient)
      outcomes = list(
        child_vitA = list(
          tag            = "child_vitA",
          label          = "Vitamin A (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cRBPAdjThurn",
          binary         = "gw_cVAD_Thurn",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_vitA = list(
          tag            = "women_vitA",
          label          = "Vitamin A (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wRBPAdjThurn",
          binary         = "gw_wVAD_Thurn",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # 2026-09-15 fix: the binary was gw_cIDA_Brinda / gw_wIDA_Brinda, whose
        # Stata labels read "Iron Deficiency ANEMIA with ID from adjusted
        # ferritin" (ID + Hb < 110/120): 42% / 31% of the sample, against the
        # report's ID of 59.0% / 41.4%. The continuous target (BRINDA log
        # ferritin < log 12/15) IS iron deficiency, so the two targets described
        # different conditions. resolve_uniform_outcome() masked this on the
        # area path by re-deriving the binary from the continuous, but every
        # consumer that reads oc$binary directly (national_estimates, the
        # individual-level SL, corrected/p1, cluster_aggregation, mrp, the
        # !is.na(binary) row filter) got IDA. gw_cID_Brinda / gw_wID_Brinda are
        # the survey's own BRINDA iron-deficiency flags (GMNS Table 18 / 32).
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_LogFerAdj",
          binary         = "gw_cID_Brinda",
          cutoff         = log(12),
          cutoff_dir     = "less",
          cutoff_scale   = "log"
        ),
        women_iron = list(
          tag            = "women_iron",
          label          = "Iron deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_LogFerAdj",
          binary         = "gw_wID_Brinda",
          cutoff         = log(15),
          cutoff_dir     = "less",
          cutoff_scale   = "log"
        )
      ),

      # Predictor domain prefixes (for variable selection & ablation)
      domains = list(
        GW     = list(prefix = "gw_"),
        DHS    = list(prefix = "dhs2019_"),
        MICS   = list(prefix = "mics_"),
        IHME   = list(prefix = "ihme_"),
        LSMS   = list(prefix = "lsms_"),
        MAP    = list(prefix = "MAP_"),
        WFP    = list(prefix = "wfp_", extra = "nearest_market_id"),
        FLUNET = list(prefix = "flunet_"),
        GEE    = list(prefix = "gee_"),
        FSEC   = list(prefix = "fsec_"),
        SOIL   = list(prefix = "soil_"),
        ESPEN  = list(prefix = "espen_"),
        SPAM   = list(prefix = "spam_"),
        VAS    = list(prefix = "vas_"),
        FAO    = list(prefix = "fao_"),
        FPN    = list(prefix = "fpn_")
      ),

      # Patterns for outcome-leakage removal in gw_ predictors
      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "LogFer", "logfer",
        "IDA", "Brinda", "Thurn",
        "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
      )

      # ── Future outcomes (Gambia has no folate/B12/zinc biomarker data) ──
    ),

    # ── Ghana (GMS 2017) ────────────────────────────────────────────────
    # Ghana Micronutrient Survey. Fieldwork: 27 Apr – 9 Jun 2017.
    # 2450 individuals (1165 children, ~992 women, some overlap).
    # No gw_child_flag in merged data — derived in load_merged_data().
    # DHS available: 2014, 2016, 2017. Currently using 2014 (3-year gap).
    # TODO: Switch to DHS 2017 once Admin-2 aggregation is run for that year
    #       (src/DHS/DHS_admin2_aggregation.R already configured for Ghana 2017).
    Ghana = list(
      country     = "Ghana",
      gadm_code   = "GHA",
      survey_year = 2017L,
      dhs_year    = 2014L,    # TODO: upgrade to 2017 when pre-computed files available
      data_path   = here::here("data", "IPD", "Ghana", "Ghana_merged_dataset.rds"),
      raster_dir  = here::here("data", "Ghana_GEE rasters"),

      # Survey design columns
      cluster_id    = "gw_cnum",
      admin1_col    = "Admin1",
      admin2_col    = "Admin2",
      psu_col       = "gw_cnum",
      # 2026-09-01 CORRECTION. gw_sWeight is only the STRATUM component of the
      # design weight - it takes 3 distinct values, one per stratum. The Ghana
      # Micronutrient Survey 2017 report, Appendix 11 ("Description of Survey
      # Weights"), documents a three-step construction and states that "the
      # final survey weights for each PSU were calculated by multiplying the
      # standardized PSU weights by the standardized stratum weights".
      # Verified in the data: gw_PSU_weight * gw_sWeight == gw_PSUStrat_weight
      # to 2e-06. Using gw_sWeight alone ignores PSU selection probability and
      # shifts national prevalence by 1.3-1.5 pp (child iron 0.202 -> 0.215,
      # women iron 0.124 -> 0.140) and individual districts by up to 4.8 pp.
      weight_col    = "gw_PSUStrat_weight",
      strata_col    = "gw_strata",
      child_flag    = "gw_child_flag",   # derived in load_merged_data()
      # 2026-09-15: report estimates are for non-pregnant women and children
      # 6-59 months. Pregnant women in the GMS gave finger-prick Hb/malaria
      # only, so only 5 carry ferritin/RBP -- the filter is near no-op here
      # but keeps the population definition identical across countries.
      preg_col      = "gw_wPreg",
      child_age_col = "gw_cAgeMonths",
      child_age_range = c(6, 59),

      outcomes = list(
        child_vitA = list(
          tag            = "child_vitA",
          label          = "Vitamin A (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cRBPAdjThurn",
          binary         = "gw_cVADAdjThurn",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_vitA = list(
          tag            = "women_vitA",
          label          = "Vitamin A (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wRBPAdjThurn",
          binary         = "gw_wVADAdjThurn",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # 2026-05-13 fix (review N-2): switch from un-adjusted gw_cLogFerr /
        # gw_wLogFerr to Thurnham-adjusted log-ferritin (derived in
        # 1_GW_Ghana_data_clean.R as gw_cLogFerrAdjThurn / gw_wLogFerrAdjThurn).
        # The pre-2026-05 config paired un-adjusted continuous ferritin with
        # adjusted binary outcomes — inconsistent across countries (Gambia/SL
        # use adjusted) and within Ghana (the binary was already adjusted).
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cFerrAdjThurn",   # Thurnham-adjusted ferritin (linear µg/L).
                                                  # The log-transformed column is created by
                                                  # 1_GW_Ghana_data_clean.R but is absent from the
                                                  # current merged cache; use the linear column.
          binary         = "gw_cIDAdjThurn",
          cutoff         = 12,                    # WHO: serum ferritin <12 µg/L (children) = ID
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_iron = list(
          tag            = "women_iron",
          label          = "Iron deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wFerrAdjThurn",   # Thurnham-adjusted ferritin (linear µg/L)
          # 2026-09-15 fix: was gw_wIDA_Thurn, labelled "Iron Deficiency Anemia
          # with ID from adjusted ferritin" (ID + Hb < 120; 7.7% of the sample,
          # report IDA 8.9%). gw_wIDAdjThurn is the survey's iron-deficiency
          # flag (Thurnham ferritin < 15; report ID 13.7%, GMS Table 39) and
          # matches the continuous target and the children's gw_cIDAdjThurn.
          binary         = "gw_wIDAdjThurn",
          cutoff         = 15,                    # WHO: serum ferritin <15 µg/L (women) = ID
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),

        women_folate = list(
          tag            = "women_folate",
          label          = "Folate deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wFolate",      # serum folate, nmol/L
          binary         = "gw_wFolateDef",   # pre-computed, 0/1
          cutoff         = 10,                 # WHO: <10 nmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_b12 = list(
          tag            = "women_b12",
          label          = "Vitamin B12 deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wB12",          # serum B12, pmol/L
          binary         = "gw_wB12Def",       # pre-computed, 0/1
          cutoff         = 148,                 # WHO: <148 pmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        )
      ),

      domains = list(
        GW      = list(prefix = "gw_"),
        DHS2014 = list(prefix = "dhs2014_"),
        DHS2016 = list(prefix = "dhs2016_"),
        DHS2017 = list(prefix = "dhs2017_"),
        MICS    = list(prefix = "mics_"),
        IHME    = list(prefix = "ihme_"),
        LSMS    = list(prefix = "lsms_"),
        MAP     = list(prefix = "MAP_"),
        WFP     = list(prefix = "wfp_", extra = "nearest_market_id"),
        FLUNET  = list(prefix = "flunet_"),
        GEE     = list(prefix = "gee_"),
        FSEC    = list(prefix = "fsec_"),
        SOIL    = list(prefix = "soil_"),
        ESPEN   = list(prefix = "espen_"),
        SPAM    = list(prefix = "spam_"),
        VAS     = list(prefix = "vas_"),
        FAO     = list(prefix = "fao_"),
        FPN     = list(prefix = "fpn_")
      ),

      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "LogFer", "logfer", "LogRBP",
        "IDA", "Brinda", "Thurn", "Ferr", "ferr",
        "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
      )
    ),

    # ── Sierra Leone (SLMS 2013) ────────────────────────────────────────
    # Sierra Leone Micronutrient Survey. Fieldwork: 11 Nov – 2 Dec 2013.
    # 1477 individuals (532 children, 945 women).
    #
    # *** DATA-PROVENANCE WARNING (found 2026-09-15) ***
    # The child file GroundWork supplied ("Sierra Leone_Child data_Davis-
    # Berkeley.dta", n = 532) contains ONLY THE ANAEMIC CHILDREN of the survey.
    # Every one of its 532 children has Hb <= 10.9 g/dL; the SLMS report's
    # Table 23 gives anaemia = 532 cases of 710 children (76.3%), and Table
    # A8-10's severity counts (190 mild / 315 moderate / 27 severe) match the
    # file's cAnemiaSev tabulation exactly. The report's ferritin and RBP
    # denominators are 654 children, so 168 non-anaemic children with assays
    # are missing from this file, and cIDA == cFeDefAdj for every child (all
    # are anaemic, so ID and IDA coincide). Every Sierra Leone child outcome
    # in this pipeline is therefore estimated on a selected subgroup, with
    # weights that were built for the full sample. The full child file must
    # be requested from GroundWork; until then treat SL child results as
    # anaemic-subgroup estimates (see docs/survey_report_reconciliation.md).
    #
    # Biomarker columns use different naming (no Thurn/Brinda suffix):
    #   Vitamin A: gw_cRBPAdj / gw_wRBPAdj (continuous), no binary VAD column
    #   Iron: gw_cFerrAdj / gw_wFerrAdjThurn (continuous, Thurnham), binary
    #         gw_cFeDefAdj / gw_wFeDefAdjThurn (Thurnham ferritin < 12 / < 15)
    # Binary VAD derived in load_merged_data() from continuous RBP < 0.70.
    # No MICS, LSMS, FluNet, or WFP data. DHS 2013 aligned.
    # The women's file is non-pregnant women only (wPregNow == 2 throughout),
    # so no preg_col is needed.
    SierraLeone = list(
      country     = "Sierra Leone",
      gadm_code   = "SLE",
      survey_year = 2013L,
      dhs_year    = 2013L,    # aligned with survey year
      data_path   = here::here("data", "IPD", "Sierra Leone",
                                "SierraLeone_merged_dataset.rds"),
      raster_dir  = here::here("data", "Sierra_Leone_GEE_rasters"),

      cluster_id    = "gw_cnum",
      admin1_col    = "Admin1",
      admin2_col    = "Admin2",
      psu_col       = "gw_cnum",
      weight_col    = "gw_svy_weight",
      strata_col    = "gw_hStratum",
      child_flag    = "gw_child_flag",
      child_age_col = "gw_cAge",          # months (6.0-59.6 in the file)
      child_age_range = c(6, 59),

      outcomes = list(
        child_vitA = list(
          tag            = "child_vitA",
          label          = "Vitamin A (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cRBPAdj",
          binary         = "gw_cVAD",          # derived in load_merged_data()
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_vitA = list(
          tag            = "women_vitA",
          label          = "Vitamin A (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wRBPAdj",
          binary         = "gw_wVAD",           # derived in load_merged_data()
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # 2026-09-15 fix: the binaries were gw_cIDA / gw_wIDA, labelled "Both
        # iron deficiency (Thurnham adjusted ferritin) and anemia". The survey's
        # iron-deficiency flags are gw_cFeDefAdj ("Iron deficiency using Thurnham
        # adjusted ferritin") and gw_wFeDefAdjThurn ("... (<15 mcg/L)"), coded
        # 1 = deficient / 2 = not; build_outcome_dataset() recodes {1,2}.
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cFerrAdj",
          binary         = "gw_cFeDefAdj",
          cutoff         = 12,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_iron = list(
          tag            = "women_iron",
          label          = "Iron deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          # 2026-09-15: was gw_wFerAdjBR1 (BRINDA-adjusted, 134 women < 15 =
          # 18% weighted), chosen because "gw_wFerrAdj has 0 values for women"
          # -- that column is the CHILD file's mothers' variable; the women's
          # file carries gw_wFerrAdjThurn (n = 774). Thurnham is the survey's
          # own method, is what the children's gw_cFerrAdj uses, and gives 97
          # women < 15 (12.8% weighted). NOTE the SLMS report's published 8.3%
          # (65 cases) is ferritin < 12, the CHILD cut-off, applied to women
          # (66 cases < 12 in this file; Wirth et al. 2016 Fig 1 caption says
          # "<12"); < 15 is the WHO definition and is kept here.
          continuous     = "gw_wFerrAdjThurn",
          binary         = "gw_wFeDefAdjThurn",
          cutoff         = 15,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),

        women_folate = list(
          tag            = "women_folate",
          label          = "Folate deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wFolate",      # serum folate, nmol/L
          binary         = "gw_wFolDef",      # recoded from 1/2 to 1/0
          cutoff         = 10,                 # WHO: <10 nmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_b12 = list(
          tag            = "women_b12",
          label          = "Vitamin B12 deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_B12",           # serum B12, pmol/L
          binary         = "gw_wB12DefWHO",   # recoded from 1/2 to 1/0
          cutoff         = 148,                # WHO: <148 pmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        )
      ),

      domains = list(
        GW     = list(prefix = "gw_"),
        DHS    = list(prefix = "dhs2013_"),
        IHME   = list(prefix = "ihme_"),
        MAP    = list(prefix = "MAP_"),
        GEE    = list(prefix = "gee_"),
        FSEC   = list(prefix = "fsec_"),
        SOIL   = list(prefix = "soil_"),
        ESPEN  = list(prefix = "espen_"),
        SPAM   = list(prefix = "spam_"),
        VAS    = list(prefix = "vas_"),
        FAO    = list(prefix = "fao_"),
        FPN    = list(prefix = "fpn_")
      ),

      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "Ferr", "ferr", "FeD", "IDA",
        "Brinda", "Thurn",
        "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
      )
    )

    # ── Malawi (MNS 2015-16, within DHS round) ──────────────────────────
    # Malawi Micronutrient Survey. Fieldwork: 8 Dec 2015 – 16 Feb 2016.
    # Biomarker columns use raw names (rbp, fer, vitA_def, iron_def)
    # rather than gw_ prefix. Population filtering uses a text `population`
    # column, not a binary child_flag.
    # DHS 2015 used (−1yr); 78% of interviews in Jan–Feb 2016.
    , Malawi = list(
      country     = "Malawi",
      gadm_code   = "MWI",
      survey_year = 2016L,    # majority of fieldwork in Jan-Feb 2016
      dhs_year    = 2015L,    # closest available DHS (−1 year)
      data_path   = here::here("data", "IPD", "Malawi",
                                "Malawi_merged_dataset.rds"),
      raster_dir  = here::here("data", "Malawi_GEE_rasters"),
      cluster_id  = "gw_cnum",
      admin1_col  = "Admin1",
      admin2_col  = "Admin2",
      psu_col     = "gw_cnum",
      weight_col  = "svy_weight",
      strata_col  = NULL,          # no strata in Malawi MNS
      child_flag  = "population",  # text column, not binary
      # 2026-09-15: the MNS report tabulates biomarkers for NON-PREGNANT WRA
      # (34 of 838 women were pregnant; 31 of them have RBP/B12/zinc assays,
      # and the survey's own sf_reg is already NA for them). Filter them out
      # of every women's outcome, as the report does. age_month is 6-59 in
      # the file, so the child age filter is a guard only.
      preg_col        = "preg",           # 1 = pregnant; NA = unknown (kept)
      child_age_col   = "age_month",
      child_age_range = c(6, 59),

      outcomes = list(
        child_vitA = list(
          tag            = "child_vitA",
          label          = "Vitamin A (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "rbp",
          binary         = "vitA_def",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_vitA = list(
          tag            = "women_vitA",
          label          = "Vitamin A (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "rbp",
          binary         = "vitA_def",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # 2026-09-15 fix: the binary was `iron_def`, derived in
        # src/malawi/1_DHS_mn_data.R by a home-grown "BRINDA" that regressed
        # log ferritin on log CRP + log AGP and took exp(residual + intercept)
        # for EVERY child -- i.e. it re-centred everyone to CRP = 1 mg/L and
        # AGP = 1 g/L, raising ferritin for the (majority) low-inflammation
        # children instead of lowering it only above the BRINDA reference
        # deciles. It flagged 107 of 1102 children (9.7%) against the survey's
        # own inflammation-corrected flag sf_c1 (222 = 20.1%; report Table 7.1:
        # 21.7% weighted) and 74 vs 131 women (report 15.1%). sf_c1 is exactly
        # sf_reg < 12 (PSC) / < 15 (WRA), the survey's BRINDA internal-regression
        # ferritin (report section 2.9), so binary and continuous now agree.
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "sf_reg",
          binary         = "sf_c1",
          cutoff         = 12,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_iron = list(
          tag            = "women_iron",
          label          = "Iron deficiency (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "sf_reg",
          binary         = "sf_c1",
          cutoff         = 15,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),

        women_folate = list(
          tag            = "women_folate",
          label          = "Folate deficiency (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "fol_nmol",         # serum folate, nmol/L (raw `fol` is already nmol/L)
          binary         = "folate_def",       # derived: fol_nmol < 10 nmol/L (WHO)
          cutoff         = 10,                 # WHO serum-folate deficiency <10 nmol/L,
                                               # harmonized with Ghana/SL (gw_wFolate < 10 nmol/L).
                                               # 2026-06-16 fix: was 3 (a <3 ng/mL threshold applied
                                               # to a nmol/L column) -> Malawi folate def ~0% vs
                                               # 55-80% elsewhere; a cross-survey bug, not biology.
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_b12 = list(
          tag            = "women_b12",
          label          = "Vitamin B12 deficiency (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "vitb12",            # serum B12 in pg/mL (RP-01, 2026-09-08: the report's cut-offs are pmol/L;
          unit_factor    = 1 / 1.355,           #   1 pmol/L = 1.355 pg/mL; converted before the cut-off is applied)
          binary         = "b12_def",           # derived: vitb12 / 1.355 < 148 pmol/L
          cutoff         = 148,                 # WHO: <148 pmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # 2026-09-15: binary switched from the locally recomputed `zinc_def` to
        # the survey's own `low_zn` (IZiNCG cut-offs by age, time of draw and
        # fasting status; report Tables 9.1 / 9.3: 60.4% PSC, 62.5% WRA). The two
        # agreed for 1085/1086 children (the exception is a -100 sentinel in
        # zn_gdl that zinc_def coded as deficient; load_merged_data() now sets
        # zn_gdl <= 0 to NA) and differed for 13 women. Zinc is NOT
        # inflammation-adjusted, matching the report; Likoswe 2020 / Gebremedhin
        # 2020 show adjustment lowers the PSC prevalence by 2-10 pp.
        child_zinc = list(
          tag            = "child_zinc",
          label          = "Zinc deficiency (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "zn_gdl",            # serum zinc, ug/dL
          binary         = "low_zn",             # survey flag, IZiNCG cutoffs
          cutoff         = 65,                   # IZiNCG morning/fasting (conservative)
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_zinc = list(
          tag            = "women_zinc",
          label          = "Zinc deficiency (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "zn_gdl",
          binary         = "low_zn",
          cutoff         = 66,                   # IZiNCG morning/fasting (conservative)
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        # Added 2026-09-15 (MW-SE / MW-IO). The MNS 2015-16 also measured plasma
        # selenium (all groups; Phiri et al. 2019) and urinary iodine (women,
        # school-age children). Malawi-only outcomes: no cross-border test, in-
        # country arms only (scripts/protocol_v2/61_malawi_selenium_iodine.R).
        # `sel_def` / `iod_def` are derived in load_merged_data() from the
        # cut-offs below. Selenium is the one outcome in the panel with a known
        # geochemical driver (soil pH / soil Se), so it is the cleanest
        # soil -> biomarker test case the project has.
        child_selenium = list(
          tag            = "child_selenium",
          label          = "Selenium deficiency (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "sel",               # plasma selenium, ug/L
          binary         = "sel_def",           # derived: sel < 84.6 ug/L
          cutoff         = 84.6,                # 1.07 umol/L, optimal GPX3 activity (Thomson 2004);
                                                # the threshold the MNS analyses used (Phiri et al. 2019)
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_selenium = list(
          tag            = "women_selenium",
          label          = "Selenium deficiency (women)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "sel",
          binary         = "sel_def",
          cutoff         = 84.6,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_iodine = list(
          tag            = "women_iodine",
          label          = "Iodine insufficiency (women, UIC < 100 ug/L)",
          population     = "women",
          child_flag_val = "women",
          continuous     = "iod",               # urinary iodine concentration, ug/L
          binary         = "iod_def",           # derived: iod < 100 ug/L (WHO insufficient, non-pregnant)
          cutoff         = 100,                 # WHO classifies POPULATION iodine status by the median UIC;
                                                # the individual < 100 share is the convention IO-01 used
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        )
      ),

      domains = list(
        GW     = list(prefix = "gw_"),
        DHS    = list(prefix = "dhs2015_"),
        IHME   = list(prefix = "ihme_"),
        MAP    = list(prefix = "MAP_"),
        GEE    = list(prefix = "gee_"),
        WFP    = list(prefix = "wfp_"),
        FSEC   = list(prefix = "fsec_"),
        SOIL   = list(prefix = "soil_"),
        ESPEN  = list(prefix = "espen_"),
        SPAM   = list(prefix = "spam_"),
        VAS    = list(prefix = "vas_"),
        FAO    = list(prefix = "fao_"),
        FPN    = list(prefix = "fpn_")
      ),

      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "Ferr", "ferr", "FeD", "IDA",
        "Brinda", "Thurn",
        "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
      )
    )

    # Tanzania (TDHS 2009-10) is DROPPED as of 2026-08-26 — see
    # get_country_config_tanzania_archived_2010() below for why, and for the
    # preserved config to restore once a replacement dataset is available.
  )
}

#' ARCHIVED — Tanzania (TDHS 2009-10) country config. NOT called anywhere.
#'
#' Dropped from get_country_configs() on 2026-08-26 after a project
#' collaborator (Omar), on the 2026-08 UC Davis / BMGF bi-weekly call,
#' identified the TDHS 2010 dried-blood-spot RBP (vitamin-A) measurements as
#' unreliable ("completely wrong... totally incorrect... not even internally
#' consistent"), consistent with independent data-quality anomalies found in
#' this pipeline (a unit-label error in `TZ61BIOMARKER.DOC`, and catastrophic
#' national-level bias — 77-92 pp — when Tanzania was held out under
#' leave-one-country-out evaluation). His recommendation was to wait for the
#' Tanzania 2023 DHS round instead. This function preserves the exact config
#' used for the 2010 round so it is trivial to compare against or adapt once
#' 2023 data is accessible: copy the `Tanzania = list(...)` value below back
#' into get_country_configs()'s returned list, updating survey_year/dhs_year,
#' data_path, and any changed variable names for the new round. The upstream
#' cleaning scripts (`src/Tanzania/`) are similarly archived, not deleted —
#' see `archive/src/Tanzania/`.
#' @return named list with one element, `Tanzania` (matching the shape of
#'   every other entry in get_country_configs())
get_country_config_tanzania_archived_2010 <- function() {
  list(
    # ── Tanzania (TDHS 2009-10, DHS micronutrient module) ────────────────
    # SCAFFOLD — outcome half only. Unlike the other four countries (which
    # are standalone micronutrient surveys with analyst-prepped biomarker
    # files), Tanzania is a DHS round: the biomarker data live in the OB
    # ("other biomarkers") recode keyed by blood-sample bar code and must be
    # linked to the PR recode (cluster, weight) and the GE GPS shapefile.
    # That linkage is built in archive/src/Tanzania/1_GW_Tanzania_data_clean.R.
    #
    # VITAMIN A ONLY. The TDHS 2010 OB module measured RBP (vit A), sTfR +
    # CRP (iron), and urinary/salt iodine. It did NOT measure serum ferritin,
    # so iron is NOT comparable to the ferritin-based definition used by the
    # other countries and is intentionally omitted. Iodine does not fit the
    # individual→admin-2 prevalence framework (population-median metric,
    # single country, no LOCO) and is omitted. See README_TANZANIA_TODO.md.
    Tanzania = list(
      country     = "Tanzania",
      gadm_code   = "TZA",
      survey_year = 2010L,
      dhs_year    = 2010L,    # same survey (DHS micronutrient module)
      data_path   = here::here("data", "IPD", "Tanzania 2010",
                                "Tanzania_merged_dataset.rds"),
      raster_dir  = here::here("data", "Tanzania_GEE_rasters"),

      # Survey design columns (emitted by 1_GW_Tanzania_data_clean.R)
      # NOTE: no micronutrient-subsample weight exists in TDHS 2010; gw_svy_weight
      # is the household weight (HV005/1e6). The MN biomarkers were a subsample,
      # so this is the known subsample-weight caveat (cf. Gambia). Document in
      # the manuscript; revisit if a subsample weight is recoverable.
      cluster_id    = "gw_cnum",
      admin1_col    = "Admin1",
      admin2_col    = "Admin2",
      psu_col       = "gw_cnum",
      weight_col    = "gw_svy_weight",
      strata_col    = "gw_strata",
      child_flag    = "gw_child_flag",

      outcomes = list(
        child_vitA = list(
          tag            = "child_vitA",
          label          = "Vitamin A (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cRBPAdjCRP",   # CRP-adjusted RBP, µmol/L (OB: rbpadcrp)
          binary         = "gw_cVAD",          # derived in clean: RBP < 0.70
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_vitA = list(
          tag            = "women_vitA",
          label          = "Vitamin A (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wRBPAdjCRP",
          binary         = "gw_wVAD",
          cutoff         = 0.70,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        )
        # No iron (sTfR not ferritin-comparable), folate, B12, zinc, or iodine.
      ),

      # Tanzania-specific DHS year prefix. Other domains follow the standard
      # extraction scripts (none run yet — see README_TANZANIA_TODO.md).
      domains = list(
        GW     = list(prefix = "gw_"),
        DHS    = list(prefix = "dhs2010_"),
        IHME   = list(prefix = "ihme_"),
        MAP    = list(prefix = "MAP_"),
        GEE    = list(prefix = "gee_"),
        WFP    = list(prefix = "wfp_", extra = "nearest_market_id"),
        FSEC   = list(prefix = "fsec_"),
        SOIL   = list(prefix = "soil_"),
        LSMS   = list(prefix = "lsms_"),
        ESPEN  = list(prefix = "espen_"),
        SPAM   = list(prefix = "spam_"),
        VAS    = list(prefix = "vas_"),
        FAO    = list(prefix = "fao_"),
        FPN    = list(prefix = "fpn_")
      ),

      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "stfr", "sTfR", "STFR", "CRP", "crp",
        "iodine", "Iodine", "UIC", "uic",
        "Ferr", "ferr", "FeD", "IDA", "Brinda", "Thurn",
        "Folate", "folate", "Fol", "B12", "b12", "Zinc", "zinc"
      )
    )
  )
}


#' Return global pipeline parameters
#'
#' @param mode "fast" for quick beta testing (~minutes), "full" for publication
#'   quality (~hours). Fast mode uses a minimal 3-learner SL stack (mean +
#'   glmnet + ranger) and 10-20 bootstrap replicates. Full mode uses the
#'   complete stack (5 learners + screener pipelines) and 200 bootstrap reps.
#' @return named list of pipeline parameters
get_pipeline_params <- function(mode = Sys.getenv("PIPELINE_MODE", "fast")) {

  mode <- match.arg(mode, c("fast", "full"))

  base <- list(
    seed = 12345L,
    mode = mode,

    # ── Feature-engineering options (validated in sandbox_fe/; see FINDINGS.md) ──
    # All default to LEGACY behaviour so existing runs are unchanged. Toggle via
    # environment variables for a clean A/B against the current pipeline.
    #
    #   FE_NORMALIZE=rank   -> rank/quantile scaling (step_percentile) instead of
    #                          z-score. +0.003–0.012 AUC where it helps; never a
    #                          material loss; transfers better. Safe to adopt.
    #   FE_BUNDLES=true     -> restrict predictors to an outcome-specific,
    #                          biology-driven bundle (vitamin A -> environment;
    #                          iron/folate/B12/zinc -> health). Matches/beats the
    #                          full set with 1/3–1/8 the features (parsimony +
    #                          interpretability); improves women_iron transfer.
    #
    # NOTE: unsupervised per-domain PCA is intentionally NOT offered here. The
    # metadata evaluation (sandbox_fe/17_pca_metadata_eval.R) showed prefix
    # "domains" are incoherent for PCA and climate-construct loadings are not
    # transportable across countries. Experimental construct-level reduction
    # lives in R/feature_engineering_constructs.R, unwired by design.
    normalize_method    = tolower(Sys.getenv("FE_NORMALIZE", "zscore")),
    use_outcome_bundles = tolower(Sys.getenv("FE_BUNDLES", "false")) %in%
                            c("1", "true", "yes"),

    #   GEE_COVARIATE_HYGIENE=true -> drop the cross-band `_annual_*` summaries
    #                          taken over NON-commensurable bands (soil depth
    #                          means averaged with their stdevs, FLDAS averaging
    #                          pressure with wind, land-cover class codes, ...)
    #                          plus static layers' identical per-year copies.
    #                          See R/gee_band_semantics.R. Recorded here so the
    #                          setting is captured in `pipeline_params` and thus
    #                          in the target hash for everything downstream of
    #                          it -- an env var alone is invisible to {targets}.
    #                          Area-level targets read the flag directly, so
    #                          after changing it run:
    #                            targets::tar_invalidate(matches("area_|loco"))
    gee_covariate_hygiene = gee_hygiene_enabled(),

    #   ANALYSIS_PROFILE=smoke -> restrict the 2026-08 analysis-extension
    #                          workstreams (nested LOCO, level decomposition,
    #                          distributional default, cluster MBG, subsample)
    #                          to one representative outcome and one country
    #                          pair, for development. "full" runs the whole
    #                          grid. This does NOT change PIPELINE_MODE, which
    #                          still controls the SL stack and bootstrap counts;
    #                          the two are independent and both are recorded
    #                          here so each enters the target hash.
    analysis_profile = tolower(Sys.getenv("ANALYSIS_PROFILE", "full"))
  )

  if (!base$analysis_profile %in% c("smoke", "full"))
    stop("ANALYSIS_PROFILE must be 'smoke' or 'full', got '",
         base$analysis_profile, "'", call. = FALSE)

  if (mode == "fast") {
    # ── Fast / beta-testing mode ────────────────────────────────────────
    # ~10-30 min total.  Good enough to verify the pipeline runs end-to-end,
    # check relative performance across outcomes, and debug visualisations.
    c(base, list(
      K        = 5L,        # CV folds (same — fewer folds barely saves time)
      B_boot   = 20L,       # LEGACY, see note below
      B_area   = 10L,       # area-level model bootstrap (already fast)
      sl_stack = "fast",    # 5 learners: mean, lasso, enet, ranger, xgboost
      sl_with_gp = FALSE,   # fast mode never included GP
      prescreen_pval = 0.2  # same prescreening threshold
    ))
  } else {
    # ── Full / publication mode ─────────────────────────────────────────
    # ~2-4 hours.  Full learner stack with screener pipelines, production
    # bootstrap replicates, and all ablation domains.
    c(base, list(
      K        = 5L,
      B_boot   = 200L,      # LEGACY, see note below
      B_area   = 500L,      # area-level model bootstrap
      sl_stack = "full",    # 12 learners: + ridge, 2 ranger, xgb_deep, 3 BART
      sl_with_gp = FALSE,   # Gaussian process: see note below
      prescreen_pval = 0.2
    ))
  }
}

# ── What actually differs between fast and full ────────────────────────────
# Verified against the code, not the comments, 2026-08:
#
#   sl_stack   5 learners (mean, lasso, elastic_net, ranger_main,
#              xgb_conservative) vs 12 (adds ridge, ranger_low_mtry,
#              ranger_deep, xgb_deep, bart_small/100/200).
#              BART matters here: it is the learner that holds up on the rare
#              outcomes (women's vitamin A and B12 sit at 1-3% prevalence)
#              where the others collapse to the mean.
#   sl_with_gp FALSE in both modes as of 2026-08. gaussianprocess used to be
#              the 13th full-mode learner. kernlab::gausspr is O(n^3) in rows:
#              measured at p = 140 it takes 0.2 s at n = 400, 1.0 s at n = 800
#              and 9.4 s at n = 1600. The area-level targets (n = 30-370) never
#              felt it, but the individual-level LOCO targets pool n = 10,011
#              and 13,107, where that curve implies 8-18 h per target. Three of
#              them ran 4.9 CPU-hours apiece without finishing before this was
#              tracked down. The gp_sensitivity target refits one small
#              country x outcome WITH the GP so the cost of dropping it is
#              measured rather than assumed.
#   B_area     10 vs 500 area-level bootstrap replicates
#              (the only other live difference; R/admin2_analysis.R)
#   K, prescreen_pval, conformal intervals, ablation domains: IDENTICAL.
#
# B_boot is retained only because R/bootstrap.R reads it, and that file is
# preserved for reference and is NOT in the targets graph (see _targets.R:83) --
# run_bootstrap_ci() has no call sites. Changing B_boot therefore has no effect
# on a tar_make() run.
#
# B_admin2 was removed in 2026-08: it was documented as "admin2 individual SL
# bootstrap, 4 vs 50" but had ZERO readers anywhere in the repo, so full mode
# never ran the 50 replicates the setting implied.


#' Outcome-specific predictor bundles (biology-driven feature selection)
#'
#' Derived from the domain-transfer analysis (sandbox_fe/12_domain_signal.R,
#' 15_feature_bundles.R): vitamin-A deficiency transfers best on environmental
#' / vegetation surfaces; iron (and other haematological / dietary) deficiencies
#' on malaria burden + modelled health surfaces. Bundles are expressed as column
#' name prefixes; only proxy-domain prefixes are used (never gw_).
#'
#' @param tag Outcome tag (e.g., "child_vitA", "women_iron").
#' @return character vector of column-name prefixes, or NULL for "use all proxy".
bundle_prefixes_for_outcome <- function(tag) {
  if (is.null(tag)) return(NULL)
  env_bundle    <- c("gee_", "soil_", "fsec_")               # environment-forward
  health_bundle <- c("MAP_", "ihme_", "fsec_")               # malaria + health burden
  if (grepl("vitA", tag, ignore.case = TRUE)) return(env_bundle)
  if (grepl("iron|folate|b12|zinc|anemia|anaemia", tag, ignore.case = TRUE))
    return(health_bundle)
  NULL  # unknown outcome -> fall back to full proxy set
}
