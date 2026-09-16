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
#   2. 02b (merge + transport) three times:
#        headline   no-DHS shards, transport open+survey_public      benchmarks_v2_*.csv
#                   (the deployable proxy set everywhere; in-fill loses nothing
#                   without DHS - RR-11: index 0.402 either way - and transport gains)
#        with DHS   default shards, transport all tiers               benchmarks_v2_*_withdhs.csv
#        open only  no-DHS shards, transport open                     benchmarks_v2_*_open.csv
#      Every row carries `predictor_tiers`.
#   3. RUN_DOWNSTREAM=1: scripts/protocol_v2/rerun_downstream.sh - every
#      protocol, policy-deck and dashboard script the decks read, at the
#      headline tiers (DOWNSTREAM_TIERS overrides)
# From the project root under Git Bash:
#   RUN_DOWNSTREAM=1 nohup bash scripts/protocol_v2/rerun_benchmarks.sh > logs/rr11_launcher.log 2>&1 &
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
EXPECTED="${EXPECTED_PREDICTORS:-542}"
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
    f="$OUT/benchmarks_v2_raw_${cn}${suf}.csv"
    { [ -s "$f" ] && [ "$f" -nt logs/rr11_lint.log ]; } || { step "ABORT: shard ${cn}${suf} produced no fresh output (logs/rr11_02_${cn}${suf}.log)"; exit 1; }
  done
}
step "1a. in-country shards, every tier"
run_shards "" ""
step "1b. in-country shards, no DHS (open,survey_public)"
run_shards "_nodhs" "open,survey_public"
step "   shards done"

# ── 2. transport arms ─────────────────────────────────────────────────────────
step "2. transport arms"
V2_SHARD_SUFFIX="_nodhs" V2_PREDICTOR_TIERS="open,survey_public" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_headline.log 2>&1; step "   headline (no-DHS shards; transport open+survey_public) exit=$?"
V2_PREDICTOR_TIERS="open,survey_public,survey_dhs" V2_OUT_SUFFIX="_withdhs" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_withdhs.log 2>&1; step "   with DHS (default shards; transport all tiers) exit=$?"
V2_SHARD_SUFFIX="_nodhs" V2_PREDICTOR_TIERS="open" V2_OUT_SUFFIX="_open" "$RS" -e "source('scripts/protocol_v2/02b_merge_and_loco.R')" > logs/rr11_02b_open.log 2>&1; step "   open only (no-DHS shards; transport open) exit=$?"
step "headline in $OUT/benchmarks_v2_summary.csv (arms: _withdhs, _open)"

# ── 3. downstream (optional) ──────────────────────────────────────────────────
if [ "${RUN_DOWNSTREAM:-0}" = "1" ]; then bash scripts/protocol_v2/rerun_downstream.sh; fi
step "done"
