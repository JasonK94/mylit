#!/bin/bash
# ==============================================================================
# MultiNicheNet re-run at anno1 level (26 cell types)
# Previous run at anno2 (8 compartments) had empty group_prioritization
# because DE was too coarse → no significant LR pairs
# ==============================================================================

set -e

BASE="/data/user3/sobj/stroke_hc_v8_2"
MNN_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R"

# ---- Layer 1: HC vs IS at anno1 ----
echo ">>> [1/2] MNN Layer 1: HC vs IS at anno1"
mkdir -p "$BASE/cci/mnn/L1_cohort_anno1"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_1_hc_is.qs" \
  -o "$BASE/cci/mnn/L1_cohort_anno1" \
  -g cohort \
  -s patient_name \
  -c anno1 \
  -f "Stroke-HC" \
  --min_cells 10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L1_cohort_anno1.log"

# ---- Layer 2: g3 at anno1 ----
echo ">>> [2/2] MNN Layer 2: g3 at anno1"
mkdir -p "$BASE/cci/mnn/L2_g3_anno1"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_2_is_g3.qs" \
  -o "$BASE/cci/mnn/L2_g3_anno1" \
  -g g3 \
  -s patient_name \
  -c anno1 \
  -f "X2-X1" \
  --min_cells 10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L2_g3_anno1.log"

echo ">>> MNN anno1 complete!"
echo "Output: $BASE/cci/mnn/L1_cohort_anno1/ and L2_g3_anno1/"
