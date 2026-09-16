#!/usr/bin/env bash
# =============================================================================
# scripts/protocol_v2/rerun_downstream.sh   [RR-11, 2026-09-15]
#
# Every protocol script, policy-deck script and dashboard bundle that the
# September decks and the online dashboard read, run AT THE HEADLINE TIERS
# (V2_PREDICTOR_TIERS, default open,survey_public - the deployable proxy set;
# DOWNSTREAM_TIERS overrides). Run after rerun_benchmarks.sh, which writes the
# benchmarks_v2_*.csv these scripts read. Each step is logged to
# logs/rr11_<name>.log and the run continues on error.
#
#   nohup bash scripts/protocol_v2/rerun_downstream.sh > logs/rr11_downstream.log 2>&1 &
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
export V2_PREDICTOR_TIERS="${DOWNSTREAM_TIERS:-open,survey_public}"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
mkdir -p logs
step() { echo "[$(date +%H:%M:%S)] $*"; }
run() { local name="$1" script="$2"; shift 2; local t0=$(date +%s); env "$@" "$RS" -e "source('$script')" > "logs/rr11_${name}.log" 2>&1; local rc=$?; step "   $name exit=$rc ($(( $(date +%s) - t0 )) s)"; [ $rc -ne 0 ] && tail -3 "logs/rr11_${name}.log"; return 0; }
[ -s results/tables/protocol_v2/benchmarks_v2_summary.csv ] || { step "ABORT: benchmarks_v2_summary.csv missing - run rerun_benchmarks.sh first"; exit 1; }
step "downstream at tiers $V2_PREDICTOR_TIERS"

step "1. protocol scripts"
for f in 05_domain_representation 12_nce_targeting_metrics 13_admin1_nested_spatial 15_training_country_curve 16_admin1_transport \
         23_domain_ablation_loco 25_nested_domain_selection 28_climate_soil_admin1 33_transport_null_calibration 35_anchor_and_rank \
         43_source_ablation_loco 44_risk_category_accuracy 47_iodine_in_country 56_weight_sources 57_index_importance 61_malawi_selenium_iodine \
         06_build_variable_sheet 19_sl_with_domain_index 20_sl_rank_loss 30_training_curve_climate_soil 34_variance_components_ceiling \
         36_urbanicity_conditioning 45_lsms_flunet_block 46_individual_level_models; do
  s="scripts/protocol_v2/${f}.R"; [ -f "$s" ] || { step "   $f: not found, skipped"; continue; }
  run "$f" "$s"
done
run 19_sl_with_domain_index_hapc scripts/protocol_v2/19_sl_with_domain_index.R SL_HAPC=1
# 39 is the parameterised add-on harness (needs ADDON_FILE); the Ghana deck reads its LSMS / FluNet run (AD-lsms), which
# scores the block script 45 has just rebuilt against the current base vocabulary. Run unparameterised it stops at once.
run 39_addon_lsms_flunet scripts/protocol_v2/39_addon_feature_test.R ADDON_FILE=data/covariates/harmonized/predictors_admin2_lsms_flunet.csv ADDON_TAG=lsms_flunet "ADDON_DOMAIN=LSMS and FluNet"

step "2. policy-deck figures and tables"
for f in 01_figures_main 02_figure_ghana_map 03_figure_geostat 04_civ_climate_soil_prediction 05_civ_map 06_civ_rank_uncertainty 07_worst_fifth_probability; do
  run "deck_${f}" "scripts/policy_deck/${f}.R"
done

step "3. dashboard bundles"
run dashboard_bundles dashboard/data-raw/05_build_protocol_v2_bundles.R
step "done - render docs/slides/*.qmd and run dashboard/deploy.R once the numbers are checked"
