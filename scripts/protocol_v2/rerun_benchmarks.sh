#!/usr/bin/env bash
# =============================================================================
# scripts/protocol_v2/rerun_benchmarks.sh   [RR-11 launcher, 2026-09-15]
#
# The protocol-v2 re-run on the finalised 522-column set.
#   0. refuse to start unless the shared metadata has exactly EXPECTED rows
#      (rebuild trap 3: a stale set runs silently) and the join lint passes;
#      copy the previous results to results/tables/protocol_v2_pre_RR11_<date>/
#      so the before/after tables can be written
#   1. two shard sets of 02 (in-country estimands), four countries in parallel:
#        default   every tier (DHS in)                      benchmarks_v2_raw_<country>.csv
#        no DHS    V2_PREDICTOR_TIERS=open,survey_public    benchmarks_v2_raw_<country>_nodhs.csv
#   2. 02b (merge + transport) four times:
#        headline   default shards, transport open+survey_public        benchmarks_v2_*.csv
#        with DHS   default shards, transport all tiers                 benchmarks_v2_*_withdhs.csv
#        open only  default shards, transport open                      benchmarks_v2_*_open.csv
#        no DHS     no-DHS shards, transport open+survey_public         benchmarks_v2_*_nodhs.csv
#      (in-country rows are identical across the first three; every row
#      carries `predictor_tiers`)
#   3. RUN_DOWNSTREAM=1: the headline-figure scripts, sequentially, continuing
#      on error (each logged): 05 12 13 15 16 23 25 28 33 35 43 44 47 56 57 61
# From the project root under Git Bash:
#   RUN_DOWNSTREAM=1 nohup bash scripts/protocol_v2/rerun_benchmarks.sh > logs/rr11_launcher.log 2>&1 &
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
EXPECTED="${EXPECTED_PREDICTORS:-522}"
OUT=results/tables/protocol_v2
mkdir -p logs
step() { echo "[$(date +%H:%M:%S)] $*"; }

# ── 0. guards and backup ──────────────────────────────────────────────────────
n=$(($(grep -c . data/covariates/harmonized/predictors_admin2_shared_metadata.csv) - 1))
if [ "$n" -ne "$EXPECTED" ]; then step "ABORT: shared metadata has $n predictors, expected $EXPECTED (set EXPECTED_PREDICTORS to override)"; exit 1; fi
step "shared set: $n predictors"
"$RS" -e 'suppressPackageStartupMessages(library(testthat)); r <- as.data.frame(test_file("tests/testthat/test-admin2-join-lint.R", reporter = "silent")); quit(status = if (sum(r$failed)) 1 else 0)' > logs/rr11_lint.log 2>&1 \
  || { step "ABORT: the Admin-2 join lint fails (logs/rr11_lint.log)"; exit 1; }
step "join lint clean"
BK="results/tables/protocol_v2_pre_RR11_$(date +%Y%m%d)"
if [ ! -d "$BK" ]; then mkdir -p "$BK" && cp "$OUT"/*.csv "$BK"/ 2>/dev/null; step "previous results copied to $BK ($(ls "$BK" | wc -l) files)"; else step "backup $BK already exists, kept"; fi

# ── 1. shards ─────────────────────────────────────────────────────────────────
run_shards() {   # $1 = shard suffix ("" or "_nodhs"), $2 = tiers ("" = default)
  local suf="$1" tiers="$2"
  for cn in gambia ghana malawi sierraleone; do
    V2_COUNTRY=$cn V2_SHARD_SUFFIX="$suf" V2_PREDICTOR_TIERS="$tiers" "$RS" -e "source('scripts/protocol_v2/02_run_benchmarks_v2.R')" > "logs/rr11_02_${cn}${suf}.log" 2>&1 &
    sleep 20
  done
  wait
  for cn in gambia ghana malawi sierraleone; do
    [ -s "$OUT/benchmarks_v2_raw_${cn}${suf}.csv" ] || { step "ABORT: shard ${cn}${suf} produced no output (logs/rr11_02_${cn}${suf}.log)"; exit 1; }
  done
}
step "1a. in-country shards, every tier"
run_shards "" ""
step "1b. in-country shards, no DHS (open,survey_public)"
run_shards "_nodhs" "open,survey_public"
step "   shards done"

# ── 2. transport arms ─────────────────────────────────────────────────────────
step "2. transport arms"
"$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_headline.log 2>&1; step "   headline (default shards; transport open+survey_public) exit=$?"
V2_PREDICTOR_TIERS="open,survey_public,survey_dhs" V2_OUT_SUFFIX="_withdhs" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_withdhs.log 2>&1; step "   with DHS exit=$?"
V2_PREDICTOR_TIERS="open" V2_OUT_SUFFIX="_open" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_open.log 2>&1; step "   open only exit=$?"
V2_SHARD_SUFFIX="_nodhs" V2_PREDICTOR_TIERS="open,survey_public" V2_OUT_SUFFIX="_nodhs" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_nodhs.log 2>&1; step "   no DHS anywhere exit=$?"
step "headline in $OUT/benchmarks_v2_summary.csv (arms: _withdhs, _open, _nodhs)"

# ── 3. downstream headline-figure scripts (optional) ─────────────────────────
if [ "${RUN_DOWNSTREAM:-0}" = "1" ]; then
  step "3. downstream scripts (headline figures; each logged, continuing on error)"
  for f in 05_domain_representation 12_nce_targeting_metrics 13_admin1_nested_spatial 15_training_country_curve 16_admin1_transport \
           23_domain_ablation_loco 25_nested_domain_selection 28_climate_soil_admin1 33_transport_null_calibration 35_anchor_and_rank \
           43_source_ablation_loco 44_risk_category_accuracy 47_iodine_in_country 56_weight_sources 57_index_importance 61_malawi_selenium_iodine; do
    s="scripts/protocol_v2/${f}.R"; [ -f "$s" ] || { step "   $f: script not found, skipped"; continue; }
    t0=$(date +%s); "$RS" -e "source('$s')" > "logs/rr11_${f}.log" 2>&1; rc=$?
    step "   $f exit=$rc ($(( $(date +%s) - t0 )) s)"
  done
fi
step "done"
