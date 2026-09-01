#!/usr/bin/env bash
# =============================================================================
# scripts/accuracy_impact/run_overnight_rebuild.sh
#
# Rebuild the predictor vocabulary from Gambia's re-extracted GEE output, then
# re-run every analysis that depends on it, in dependency order.
#
# WHY A DRIVER AND NOT A SEQUENCE OF MANUAL CALLS
# -----------------------------------------------
# Gambia's Admin-2 GEE extraction went from 222 to 435 columns, so it is no
# longer the binding constraint on the harmonised intersection. Everything
# downstream of that file is stale: the harmonised set, the shared set, and
# every result computed on either. Running them by hand invites a half-updated
# state where some tables reflect the new vocabulary and some do not, which is
# the failure this project has spent two days cataloguing.
#
# Each stage logs to its own file, and a stage that fails does NOT stop the
# ones that do not depend on it. The summary at the end reports what ran, what
# failed, and what was skipped, so a partial night is legible in the morning.
#
#   bash scripts/accuracy_impact/run_overnight_rebuild.sh
# =============================================================================
set -u
cd "$(dirname "$0")/../.." || exit 1
export R_USER="C:/Users/andre/OneDrive/Documents"
export HOME="C:/Users/andre/OneDrive/Documents"
export COVARIATE_VOCAB=harmonized
RS="/c/Program Files/R/R-4.4.2/bin/Rscript.exe"
LOG=logs/overnight_$(date +%Y%m%d_%H%M%S)
mkdir -p "$LOG"
echo "logs -> $LOG"

STAGES=()
run () {                      # run <name> <script> [args...]
  local name="$1"; shift
  local t0=$(date +%s)
  echo "[$(date +%H:%M:%S)] START  $name"
  "$RS" "$@" > "$LOG/$name.log" 2>&1
  local rc=$?
  local dt=$(( $(date +%s) - t0 ))
  if [ $rc -eq 0 ]; then echo "[$(date +%H:%M:%S)] OK     $name  (${dt}s)"
  else echo "[$(date +%H:%M:%S)] FAIL   $name  (rc=$rc, ${dt}s)"; fi
  STAGES+=("$name:$rc:${dt}s")
  return $rc
}

# --- 1. vocabulary ---------------------------------------------------------
# 02 rebuilds each country's raw predictor table from its source files, which
# is where Gambia's new GEE columns enter. 03 re-derives the harmonised
# intersection from those tables.
run build_country_predictors scripts/covariates/02_build_country_predictors.R
run harmonize                scripts/covariates/03_harmonize.R
run verify_harmonization     scripts/covariates/07_verify_harmonization.R

# --- 2. extra domains and the shared set -----------------------------------
run harmonize_extra_domains  scripts/covariates/harmonize_extra_domains.R
run build_shared             scripts/covariates/build_shared_predictor_set.R

# --- 3. analyses that consume the shared set -------------------------------
run within_country           scripts/accuracy_impact/wsl1_within_country.R
run loco_decomposition       scripts/accuracy_impact/wsi1_loco_decomposition.R
run targeting_lift           scripts/accuracy_impact/wsj1_targeting_lift.R
run predictor_consistency    scripts/accuracy_impact/wsa_predictor_consistency.R

# --- 4. the test suite, last, so a broken guard is visible in the morning ---
run tests                    tests/testthat.R

echo
echo "================ SUMMARY ================"
printf '%s\n' "${STAGES[@]}" | while IFS=: read -r n rc dt; do
  if [ "$rc" = "0" ]; then printf "  OK    %-26s %s\n" "$n" "$dt"
  else printf "  FAIL  %-26s %s (rc=%s)\n" "$n" "$dt" "$rc"; fi
done
echo
echo "predictor counts:"
for f in data/covariates/harmonized/predictors_admin2_harmonized.csv \
         data/covariates/harmonized/predictors_admin2_shared.csv; do
  [ -f "$f" ] && echo "  $(basename "$f"): $(head -1 "$f" | tr ',' '\n' | wc -l) columns"
done
echo "logs in $LOG"
