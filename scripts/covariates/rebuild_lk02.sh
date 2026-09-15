#!/usr/bin/env bash
# =============================================================================
# scripts/covariates/rebuild_lk02.sh
#
# Orchestrates the LK-02 / AB-01 predictor-set rebuild AFTER the Malawi DHS
# files have been regenerated and the first shared-set build has run
# (logs/lk02_builder_pass1.log). Run from the project root under Git Bash:
#
#   bash scripts/covariates/rebuild_lk02.sh > logs/lk02_orchestrator.log 2>&1 &
#
# Steps
#   1. wait for builder pass 1 (stage 2-4 + verify + builder chain)
#   2. cluster-model BYM2 for the DHS columns new to the shared set, three
#      countries in parallel (DHS_CM_ONLY_NEW=1)
#   3. wait for the Malawi full refit (MNS clusters excluded), then its
#      incremental pass for the new columns
#   4. builder pass 2 -> 07 -> 08 -> 59 -> [60 if data/ACLED/*.csv] -> 53
#   5. audit tables
# Models are NOT run.
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
ts() { date +"%H:%M:%S"; }
step() { echo "[$(ts)] $*"; }

step "1. waiting for builder pass 1"
until grep -q "predictors_admin2_shared_metadata.csv" logs/lk02_builder_pass1.log 2>/dev/null; do sleep 30; done
step "   builder pass 1 done: $(grep -c . logs/lk02_builder_pass1.log) log lines"

step "2. incremental cluster-model fits for Gambia / Ghana / SierraLeone (parallel, staggered)"
# Staggered by 90 s: three INLA namespaces attaching at the same instant raced
# on a shared temp file (libloc_*.rds) and two of three died on 2026-09-15.
# One retry per country after a failed exit.
run_incremental() {
  local cn=$1
  DHS_COUNTRY=$cn DHS_CM_ONLY_NEW=1 "$RS" -e "source('scripts/covariates/build_dhs_admin2_clustermodel.R')" > "logs/lk02_clustermodel_${cn}_new.log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "exit=$rc (retrying once)" >> "logs/lk02_clustermodel_${cn}_new.log"; sleep 60
    DHS_COUNTRY=$cn DHS_CM_ONLY_NEW=1 "$RS" -e "source('scripts/covariates/build_dhs_admin2_clustermodel.R')" >> "logs/lk02_clustermodel_${cn}_new.log" 2>&1
    rc=$?
  fi
  echo "exit=$rc" >> "logs/lk02_clustermodel_${cn}_new.log"
}
for cn in Gambia Ghana SierraLeone; do run_incremental $cn & sleep 90; done
wait
for cn in Gambia Ghana SierraLeone; do step "   $cn: $(grep '^exit=' logs/lk02_clustermodel_${cn}_new.log | tail -1) $(grep -c ' areas ' logs/lk02_clustermodel_${cn}_new.log) fits"; done

step "3. waiting for the Malawi full refit"
until grep -q "^exit=" logs/lk02_clustermodel_Malawi.log 2>/dev/null; do sleep 60; done
step "   Malawi full refit: $(grep '^exit=' logs/lk02_clustermodel_Malawi.log)"
DHS_COUNTRY=Malawi DHS_CM_ONLY_NEW=1 "$RS" -e "source('scripts/covariates/build_dhs_admin2_clustermodel.R')" > logs/lk02_clustermodel_Malawi_new.log 2>&1; echo "exit=$?" >> logs/lk02_clustermodel_Malawi_new.log
step "   Malawi incremental: $(grep '^exit=' logs/lk02_clustermodel_Malawi_new.log)"

step "4. builder pass 2 -> 07 -> 08 -> 59 -> [60] -> 62 -> 53 -> stamp -> spine"
"$RS" -e "source('scripts/covariates/build_shared_predictor_set.R')" > logs/lk02_builder_pass2.log 2>&1; step "   builder exit=$?"
"$RS" -e "source('scripts/protocol_v2/07_build_food_environment.R')" > logs/lk02_07_food.log 2>&1; step "   07 exit=$?"
"$RS" -e "source('scripts/protocol_v2/08_build_extra_sources.R')" > logs/lk02_08_extrasrc.log 2>&1; step "   08 exit=$?"
"$RS" -e "source('scripts/protocol_v2/59_build_addback_sources.R')" > logs/lk02_59_addback.log 2>&1; step "   59 exit=$?"
if ls data/ACLED/*.csv >/dev/null 2>&1; then
  "$RS" -e "source('scripts/protocol_v2/60_build_acled_conflict.R')" > logs/lk02_60_acled.log 2>&1; step "   60 exit=$?"
else step "   60 skipped (no ACLED export under data/ACLED/)"; fi
"$RS" -e "source('scripts/protocol_v2/62_append_blocks.R')" > logs/lk02_62_blocks.log 2>&1; step "   62 exit=$? (HCES, RTFP, MODIS NDVI blocks; each builder runs separately)"
"$RS" -e "source('scripts/protocol_v2/53_apply_exclusion_policy.R')" > logs/lk02_53_policy.log 2>&1; step "   53 exit=$?"
"$RS" -e "source('scripts/covariates/stamp_predictor_metadata.R')" > logs/lk02_stamp.log 2>&1; step "   stamp exit=$? (tier, subnational, year_used)"
"$RS" -e "source('scripts/covariates/build_admin2_spine.R')" > logs/lk02_spine.log 2>&1; step "   spine exit=$?"

step "5. audit"
"$RS" -e "source('scripts/covariates/audit_predictor_set.R')" > logs/lk02_audit.log 2>&1; step "   audit exit=$?"
step "shared set columns: $(head -1 data/covariates/harmonized/predictors_admin2_shared.csv | tr ',' '\n' | wc -l) (incl. 3 keys); metadata rows: $(($(wc -l < data/covariates/harmonized/predictors_admin2_shared_metadata.csv) - 1))"
echo "[$(ts)] REBUILD COMPLETE" > logs/lk02_REBUILD_COMPLETE
step "done"
