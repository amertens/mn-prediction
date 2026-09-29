#!/usr/bin/env bash
# =============================================================================
# scripts/protocol_v2/rerun_survey_fixes_2026-09-27.sh   [RR-14, 2026-09-27]
#
# Propagates the 15 September survey-report corrections
# (docs/survey_report_reconciliation.md: vitamin A on the RBP < 0.70 rule,
# non-pregnant women and children 6-59 months only, Sierra Leone women's iron
# on Thurnham, Malawi repeated Traditional Authority names, Malawi zinc on the
# survey's own flag) into every result table the MNF15 decks read. The
# modelling targets (targets_v2.csv) were last built on 8 September.
#
#   0. snapshot the current results OUTSIDE OneDrive, and stop unless the copy
#      is complete: results/tables, results/figures, results/models,
#      data/covariates/harmonized, _targets_full, logs
#   1. the pipeline's outcome datasets for the four countries (full-mode store)
#   2. protocol-v2 script 01 (targets_v2.csv), stop unless it was rewritten
#   3. benchmarks + transport + downstream (rerun_benchmarks.sh, RUN_DOWNSTREAM=1,
#      SKIP_DASHBOARD=1: the dashboard's bundles are left to the session that owns it)
#   4. tables the decks read that the standard downstream does not rebuild
#      (external validation 04, conformal bands 66, survey planner 64, the
#      policy-deck tables 20, 24, 25, 28, and the map tables of 10)
# Figure scripts with pinned numbers (26, 29, 31, the v4 deck) are NOT run here:
# they stop on purpose when a number moves and are re-run after review.
#
# Not included: re-cleaning Malawi from the raw MNS files (src/Malawi): the
# fixed pipeline reads the survey's own flags (sf_reg, sf_c1, low_zn), which
# the current Malawi dataset already carries.
#
# Restore: copy the snapshot folders back over the project (see SNAP below).
#   nohup bash scripts/protocol_v2/rerun_survey_fixes_2026-09-27.sh > logs/rr14_survey_fixes.log 2>&1 &
# =============================================================================
set -u
cd "C:/Users/andre/OneDrive/Documents/mn-prediction"
export R_USER="C:/Users/andre/OneDrive/Documents" HOME="C:/Users/andre/OneDrive/Documents"
RS="C:/Program Files/R/R-4.4.2/bin/Rscript.exe"
SNAP="C:/Users/andre/mn-prediction-snapshots/pre_rebuild_2026-09-27"
step() { echo "[$(date '+%a %H:%M:%S')] $*"; }
run() { local name="$1" script="$2"; shift 2; local t0=$(date +%s); env "$@" "$RS" -e "source('$script')" > "logs/rr14_${name}.log" 2>&1; local rc=$?; step "   $name exit=$rc ($(( $(date +%s) - t0 )) s)"; [ $rc -ne 0 ] && tail -3 "logs/rr14_${name}.log"; return 0; }

# ── 0. snapshot ───────────────────────────────────────────────────────────────
step "0. snapshot to $SNAP"
if [ -f "$SNAP/SNAPSHOT_COMPLETE" ]; then
  step "   complete snapshot already there (taken $(cat "$SNAP/SNAPSHOT_COMPLETE")); kept, not overwritten"
else
if [ -e "$SNAP" ]; then step "ABORT: $SNAP exists but is not marked complete; check it by hand"; exit 1; fi
mkdir -p "$SNAP/results" "$SNAP/data/covariates"
for d in results/tables results/figures results/models data/covariates/harmonized _targets_full logs; do
  cp -a "$d" "$SNAP/$d" || { step "ABORT: copying $d failed"; exit 1; }
  a=$(find "$d" -type f | wc -l); b=$(find "$SNAP/$d" -type f | wc -l)
  [ "$a" -eq "$b" ] || { step "ABORT: $d has $a files, snapshot $b"; exit 1; }
  step "   $d: $a files copied"
done
cp -a docs/slides/MNF15-talk-2026-09-v4.pptx "$SNAP/" 2>/dev/null
echo "Snapshot of the results before the RR-14 rebuild ($(date)). Restore a folder by copying it back into C:/Users/andre/OneDrive/Documents/mn-prediction/." > "$SNAP/README.txt"
date > "$SNAP/SNAPSHOT_COMPLETE"
step "   snapshot complete ($(du -sh "$SNAP" | cut -f1))"
fi

# ── 1. outcome datasets ──────────────────────────────────────────────────────
step "1. pipeline outcome datasets (full-mode store)"
T0=$(date +%s)
PIPELINE_MODE=full "$RS" scripts/protocol_v2/rr14_outcome_datasets.R > logs/rr14_tar_make.log 2>&1 || { step "ABORT: tar_make failed (logs/rr14_tar_make.log)"; tail -20 logs/rr14_tar_make.log; exit 1; }
step "   outcome datasets rebuilt ($(( $(date +%s) - T0 )) s)"

# ── 2. modelling targets ─────────────────────────────────────────────────────
step "2. protocol-v2 targets (script 01)"
touch logs/rr14_stamp
"$RS" -e "source('scripts/protocol_v2/01_build_targets_v2.R')" > logs/rr14_01_build_targets.log 2>&1 \
  || { step "ABORT: script 01 failed (logs/rr14_01_build_targets.log)"; tail -20 logs/rr14_01_build_targets.log; exit 1; }
[ results/tables/protocol_v2/targets_v2.csv -nt logs/rr14_stamp ] || { step "ABORT: targets_v2.csv was not rewritten"; exit 1; }
step "   targets_v2.csv rewritten ($(($(wc -l < results/tables/protocol_v2/targets_v2.csv) - 1)) rows)"

# ── 3. benchmarks, transport, downstream ─────────────────────────────────────
step "3. benchmarks + downstream (dashboard skipped)"
RUN_DOWNSTREAM=1 SKIP_DASHBOARD=1 bash scripts/protocol_v2/rerun_benchmarks.sh
[ results/tables/protocol_v2/benchmarks_v2_cells.csv -nt logs/rr14_stamp ] || { step "ABORT: benchmarks were not rewritten"; exit 1; }

# ── 4. tables the decks read outside the standard downstream ─────────────────
step "4. extra tables"
export V2_PREDICTOR_TIERS="open,survey_public"
run xv04_transport_test scripts/external_validation/04_transport_test.R
run 66_conformal scripts/protocol_v2/66_conformal_prevalence_bands.R
run 64_survey_planner scripts/protocol_v2/64_survey_planner_validation.R
run deck_20_mnf15_figures scripts/policy_deck/20_mnf15_figures.R
run deck_24_civ_2007 scripts/policy_deck/24_civ_2007_survey_check.R
run deck_25_ghana_heldout scripts/policy_deck/25_mnf15_v2_ghana_heldout.R
run deck_28_percell_pairs scripts/policy_deck/28_mnf15_v3_percell_pairs.R
t0=$(date +%s); "$RS" scripts/policy_deck/10_viz_tables.R A,B,B2 > logs/rr14_deck_10_viz_default.log 2>&1; step "   deck_10_viz_default exit=$? ($(( $(date +%s) - t0 )) s)"
t0=$(date +%s); VZ_A_OUTCOMES=women_b12 VZ_A_COUNTRIES=Ghana,Malawi VZ_B2_PAIRS="Malawi:women_b12:deploy_malawi_women_b12_cal.csv" "$RS" scripts/policy_deck/10_viz_tables.R A,B2 > logs/rr14_deck_10_viz_b12.log 2>&1; step "   deck_10_viz_b12 exit=$? ($(( $(date +%s) - t0 )) s)"
step "done - review the numbers, then re-run the figure scripts (26 FIG_VER=v4, 29, 31, 01/03 with FIG_OUT) and the deck"
