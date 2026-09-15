#!/usr/bin/env bash
# =============================================================================
# scripts/covariates/rebuild_shared_tail.sh
#
# The tail of the shared-set rebuild, from the builder to the audit, with the
# DHS cluster models and the stand-alone block builders left as they are:
#   builder -> 07 -> 08 -> 59 -> [60 if data/ACLED/*.csv] -> 62 -> 53 -> stamp -> spine -> audit
# Run after any change to a block file, an exclusion rule, the tier / domain /
# alignment tables, or scripts 07/08/59. The stand-alone blocks (HCES, RTFP,
# MODIS NDVI, MICS) are built by their own scripts first when their inputs
# change. About five minutes. From the project root under Git Bash:
#   bash scripts/covariates/rebuild_shared_tail.sh
# Logs: logs/tail_<step>.log; the last line of each step reports its exit code.
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
mkdir -p logs
step() { echo "[$(date +%H:%M:%S)] $*"; }
run() { local name="$1" script="$2"; "$RS" -e "source('$script')" > "logs/tail_${name}.log" 2>&1; local rc=$?; step "   $name exit=$rc"; [ $rc -ne 0 ] && { echo "   see logs/tail_${name}.log"; tail -5 "logs/tail_${name}.log"; }; return $rc; }
step "builder -> 07 -> 08 -> 59 -> [60] -> 62 -> 53 -> stamp -> spine -> audit"
[ "${SKIP_BUILDER:-0}" = "1" ] || run builder scripts/covariates/build_shared_predictor_set.R || exit 1
[ "${SKIP_BUILDER:-0}" = "1" ] || run 07 scripts/protocol_v2/07_build_food_environment.R || exit 1
run 08 scripts/protocol_v2/08_build_extra_sources.R || exit 1
run 59 scripts/protocol_v2/59_build_addback_sources.R || exit 1
if ls data/ACLED/*.csv >/dev/null 2>&1; then run 60 scripts/protocol_v2/60_build_acled_conflict.R; else step "   60 skipped (no ACLED export under data/ACLED/)"; fi
run 62 scripts/protocol_v2/62_append_blocks.R || exit 1
run 53 scripts/protocol_v2/53_apply_exclusion_policy.R || exit 1
run stamp scripts/covariates/stamp_predictor_metadata.R || exit 1
run spine scripts/covariates/build_admin2_spine.R || exit 1
run audit scripts/covariates/audit_predictor_set.R || exit 1
step "done: $(grep -c . data/covariates/harmonized/predictors_admin2_shared_metadata.csv) metadata lines (incl. header)"
