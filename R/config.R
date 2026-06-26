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
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_LogFerAdj",
          binary         = "gw_cIDA_Brinda",
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
          binary         = "gw_wIDA_Brinda",
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
        SOIL   = list(prefix = "soil_")
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
      weight_col    = "gw_sWeight",
      strata_col    = "gw_strata",
      child_flag    = "gw_child_flag",   # derived in load_merged_data()

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
          binary         = "gw_wIDA_Thurn",
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
        SOIL    = list(prefix = "soil_")
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
    # Biomarker columns use different naming (no Thurn/Brinda suffix):
    #   Vitamin A: gw_cRBPAdj / gw_wRBPAdj (continuous), no binary VAD column
    #   Iron: gw_cFerrAdj / gw_wFerrAdj (continuous), gw_cIDA / gw_wIDA (binary)
    # Binary VAD derived in load_merged_data() from continuous RBP < 0.70.
    # No MICS, LSMS, FluNet, or WFP data. DHS 2013 aligned.
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
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = 1L,
          continuous     = "gw_cFerrAdj",
          binary         = "gw_cIDA",
          cutoff         = 12,
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        women_iron = list(
          tag            = "women_iron",
          label          = "Iron deficiency (women)",
          population     = "women",
          child_flag_val = 0L,
          continuous     = "gw_wFerAdjBR1",   # gw_wFerrAdj has 0 values for women; wFerAdjBR1 has 774
          binary         = "gw_wIDA",
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
        SOIL   = list(prefix = "soil_")
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
        child_iron = list(
          tag            = "child_iron",
          label          = "Iron deficiency (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "sf_reg",
          binary         = "iron_def",
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
          binary         = "iron_def",
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
          continuous     = "vitb12",            # serum B12, pmol/L
          binary         = "b12_def",           # derived: vitb12 < 148 pmol/L
          cutoff         = 148,                 # WHO: <148 pmol/L
          cutoff_dir     = "less",
          cutoff_scale   = "original"
        ),
        child_zinc = list(
          tag            = "child_zinc",
          label          = "Zinc deficiency (children)",
          population     = "children",
          child_flag_val = "preschool children",
          continuous     = "zn_gdl",            # serum zinc, ug/dL
          binary         = "zinc_def",           # pre-computed with IZiNCG cutoffs
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
          binary         = "zinc_def",
          cutoff         = 66,                   # IZiNCG morning/fasting (conservative)
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
        SOIL   = list(prefix = "soil_")
      ),

      gw_exclude_patterns = c(
        "RBP", "rbp", "VAD", "Ferr", "ferr", "FeD", "IDA",
        "Brinda", "Thurn",
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
                            c("1", "true", "yes")
  )

  if (mode == "fast") {
    # ── Fast / beta-testing mode ────────────────────────────────────────
    # ~10-30 min total.  Good enough to verify the pipeline runs end-to-end,
    # check relative performance across outcomes, and debug visualisations.
    c(base, list(
      K        = 5L,        # CV folds (same — fewer folds barely saves time)
      B_boot   = 20L,       # admin1 bootstrap
      B_admin2 = 4L,        # admin2 individual SL bootstrap
      B_area   = 10L,       # area-level model bootstrap (already fast)
      sl_stack = "fast",    # 3-learner stack: mean + glmnet + ranger
      prescreen_pval = 0.2  # same prescreening threshold
    ))
  } else {
    # ── Full / publication mode ─────────────────────────────────────────
    # ~2-4 hours.  Full learner stack with screener pipelines, production
    # bootstrap replicates, and all ablation domains.
    c(base, list(
      K        = 5L,
      B_boot   = 200L,      # admin1 bootstrap
      B_admin2 = 50L,       # admin2 individual SL bootstrap
      B_area   = 500L,      # area-level model bootstrap
      sl_stack = "full",    # full stack with screener pipelines
      prescreen_pval = 0.2
    ))
  }
}


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
