#!/usr/bin/env bash
# =============================================================================
# scripts/policy_deck/rr14b_deck_tables.sh   [2026-09-28]
# The deck tables the v5 main slides read, rebuilt on the post-fix benchmarks
# (RR-14) ahead of the rest of the downstream run: Ghana held out (25), the
# per-nutrient district pairs (28), the measurability screen (20), and the
# external check (external_validation/04). The RR-14 wrapper runs the same
# scripts again at its step 4; running them twice is harmless.
#   bash scripts/policy_deck/rr14b_deck_tables.sh > logs/rr14b_deck_tables.log 2>&1
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
export V2_PREDICTOR_TIERS="open,survey_public"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
step() { echo "[$(date '+%a %H:%M:%S')] $*"; }
run() { local name="$1" script="$2"; local t0=$(date +%s); "$RS" -e "source('$script')" > "logs/rr14b_${name}.log" 2>&1; local rc=$?; step "   $name exit=$rc ($(( $(date +%s) - t0 )) s)"; [ $rc -ne 0 ] && tail -3 "logs/rr14b_${name}.log"; return 0; }
step "post-fix deck tables"
run 25_ghana_heldout scripts/policy_deck/25_mnf15_v2_ghana_heldout.R
run 28_percell_pairs scripts/policy_deck/28_mnf15_v3_percell_pairs.R
run 20_mnf15_figures scripts/policy_deck/20_mnf15_figures.R
run xv04_transport_test scripts/external_validation/04_transport_test.R
step "done"
