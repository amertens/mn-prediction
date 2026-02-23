# =============================================================================
# src/analysis/config.R
#
# Central configuration for the Gambia micronutrient prediction pipeline.
# All tunable parameters live here; downstream scripts read from `cfg`.
# =============================================================================

cfg <- list(

  # ---- Countries -------------------------------------------------------
  countries = "Gambia",

  # ---- Raw data path (used with here::here()) --------------------------
  data_path = file.path("data", "IPD", "Gambia", "Gambia_merged_dataset.rds"),

  # ---- Key column names in the merged dataset --------------------------
  cluster_id  = "gw_cnum",
  admin1      = "Admin1",
  month       = "gw_month",
  child_flag  = "gw_child_flag",   # 1 = child, 0 = women

  # ---- Outcomes per (population x micronutrient) ----------------------
  #
  # continuous : column name of the continuous biomarker
  # binary     : pre-computed binary deficiency indicator (if present)
  #              Used only to compare against modelled prevalence; the SL
  #              is always fitted on the continuous scale and predictions
  #              are thresholded deterministically (see note below).
  # cutoff     : threshold applied to the *continuous* prediction to
  #              classify an individual as deficient.
  # cutoff_dir : "less"    → deficient if continuous < cutoff
  #              "greater" → deficient if continuous > cutoff
  # cutoff_scale: "original" | "log"  (documents the scale, not a transform)
  # model_tag  : short tag used in output filenames
  # population : "children" | "women"
  # child_flag_val : value of gw_child_flag that selects this population
  #
  # NOTE on modeling strategy: the task requires direct binary modeling
  # where feasible.  Pre-computed binary indicators exist for all four
  # outcomes (gw_cVAD_Thurn, gw_wVAD_Thurn, gw_cIDA_Brinda, gw_wIDA_Brinda).
  # Scripts 02–05 fit BOTH a continuous SL (deliverable A, using `slmod`)
  # and a binary SL (deliverable B/C/D, using `slmod2_bin`) and prefer the
  # binary model's predicted probabilities for prevalence estimation.
  # Thresholding of continuous predictions is implemented as a fallback and
  # is used for RMSE/MAE evaluation.

  outcomes = list(

    child_vitA = list(
      tag           = "child_vitA",
      label         = "Vitamin A (children)",
      population    = "children",
      child_flag_val = 1L,
      continuous    = "gw_cRBPAdjThurn",   # RBP (µmol/L), Thurnham-adjusted
      binary        = "gw_cVAD_Thurn",      # 1 = VAD (RBP < 0.70 µmol/L)
      cutoff        = 0.70,                  # µmol/L
      cutoff_dir    = "less",
      cutoff_scale  = "original",
      model_file_cont = "res_GW_Gambia_SL_child_vitA.rds",
      model_file_bin  = "res_bin_GW_Gambia_SL_child_vitA.rds",
      table_tag     = "children_vitA",
      scatter_tag   = "children_vitA"
    ),

    women_vitA = list(
      tag           = "women_vitA",
      label         = "Vitamin A (women)",
      population    = "women",
      child_flag_val = 0L,
      continuous    = "gw_wRBPAdjThurn",
      binary        = "gw_wVAD_Thurn",
      cutoff        = 0.70,
      cutoff_dir    = "less",
      cutoff_scale  = "original",
      model_file_cont = "res_GW_Gambia_SL_women_vitA.rds",
      model_file_bin  = "res_bin_GW_Gambia_SL_women_vitA.rds",
      table_tag     = "women_vitA",
      scatter_tag   = "women_vitA"
    ),

    child_iron = list(
      tag           = "child_iron",
      label         = "Iron deficiency (children)",
      population    = "children",
      child_flag_val = 1L,
      continuous    = "gw_LogFerAdj",      # log(ferritin µg/L), Brinda-adjusted
      binary        = "gw_cIDA_Brinda",    # 1 = IDA (ferritin < 12 µg/L)
      cutoff        = log(12),             # log scale; ferritin < 12 µg/L in children
      cutoff_dir    = "less",
      cutoff_scale  = "log",
      model_file_cont = "res_GW_Gambia_SL_child_iron.rds",
      model_file_bin  = "res_bin_GW_Gambia_SL_child_iron.rds",
      table_tag     = "children_iron",
      scatter_tag   = "children_iron"
    ),

    women_iron = list(
      tag           = "women_iron",
      label         = "Iron deficiency (women)",
      population    = "women",
      child_flag_val = 0L,
      continuous    = "gw_LogFerAdj",
      binary        = "gw_wIDA_Brinda",   # 1 = IDA (ferritin < 15 µg/L)
      cutoff        = log(15),            # log scale; ferritin < 15 µg/L in women
      cutoff_dir    = "less",
      cutoff_scale  = "log",
      model_file_cont = "res_GW_Gambia_SL_mom_iron.rds",
      model_file_bin  = "res_bin_GW_Gambia_SL_mom_iron.rds",
      table_tag     = "women_iron",
      scatter_tag   = "women_iron"
    )
  ),

  # ---- Predictor domain prefixes --------------------------------------
  # Each element: list(prefix = <string>, extra = <optional char vec>)
  # "GW" covers non-outcome gw_ predictors (leaking vars removed in 01_)
  domains = list(
    GW     = list(prefix = "gw_"),
    DHS    = list(prefix = "dhs2019_"),
    MICS   = list(prefix = "mics_"),
    IHME   = list(prefix = "ihme_"),
    LSMS   = list(prefix = "lsms_"),
    MAP    = list(prefix = "MAP_"),
    WFP    = list(prefix = "wfp_", extra = "nearest_market_id"),
    FLUNET = list(prefix = "flunet_"),
    GEE    = list(prefix = "gee_")
  ),

  # gw_ variable exclusion patterns (outcome leakage prevention)
  # Any gw_ column matching any of these patterns (grepl) is dropped from
  # predictors.  Patterns cover all four continuous/binary outcome vars.
  gw_exclude_patterns = c(
    "RBP", "rbp",        # retinol-binding protein (Vitamin A)
    "VAD",               # vitamin A deficiency indicator
    "LogFer", "logfer",  # log ferritin (Iron)
    "IDA",               # iron deficiency anaemia indicator
    "Brinda", "Thurn"    # adjustment-method suffixes used in outcome names
  ),

  # ---- Cross-validation (cluster-blocked) ----------------------------
  K = 5L,   # number of CV folds

  # ---- Bootstrap ------------------------------------------------------
  B = 200L,  # bootstrap replicates (increase for publication; 200 for dev)

  # ---- Output paths ---------------------------------------------------
  out_models  = here::here("results", "models"),
  out_tables  = here::here("results", "tables"),
  out_figures = here::here("results", "figures"),

  # ---- Reproducibility seed ------------------------------------------
  seed = 12345L
)
