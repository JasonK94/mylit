#!/bin/bash
# ==============================================================================
# CCI Analysis: CellChat + MultiNicheNet (Dual-layer)
# Layer 1: HC vs IS (5_1_hc_is.qs, cohort)
# Layer 2: g3==2 vs g3==1 (5_2_is_g3.qs, g3)
# ==============================================================================
set -e

BASE="/data/user3/sobj/stroke_hc_v8_2"
LOG_DIR="${BASE}/logs"
CELLCHAT_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R"
MNN_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R"

INPUT_L1="${BASE}/5_1_hc_is.qs"
INPUT_L2="${BASE}/5_2_is_g3.qs"

mkdir -p "$LOG_DIR"

echo "=== CCI Analysis Started: $(date) ==="

# ======================================================================
# CELLCHAT - Layer 2: g3==2 vs g3==1 (smaller dataset, run first)
# ======================================================================
echo ""
echo ">>> [1/4] CellChat Layer 2: g3 comparison (anno2)"
echo "    Started: $(date)"

Rscript "$CELLCHAT_SCRIPT" \
  -i "$INPUT_L2" \
  -g anno2 \
  -s patient_name \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  -p 0.05 \
  -m 5 \
  -c 4 \
  -o "${BASE}/cci/cellchat/L2_g3_anno2" \
  2>&1 | tee "${LOG_DIR}/cellchat_L2_g3_anno2.log"

echo "    >>> CellChat L2 anno2 completed: $(date)"

# ======================================================================
# CELLCHAT - Layer 1: HC vs IS
# ======================================================================
echo ""
echo ">>> [2/4] CellChat Layer 1: HC vs IS (anno2)"
echo "    Started: $(date)"

Rscript "$CELLCHAT_SCRIPT" \
  -i "$INPUT_L1" \
  -g anno2 \
  -s patient_name \
  -a cohort \
  --subset_aggregate "Stroke,HC" \
  -d "Secreted Signaling" \
  -p 0.05 \
  -m 5 \
  -c 4 \
  -o "${BASE}/cci/cellchat/L1_cohort_anno2" \
  2>&1 | tee "${LOG_DIR}/cellchat_L1_cohort_anno2.log"

echo "    >>> CellChat L1 anno2 completed: $(date)"

# ======================================================================
# MULTINICHENET - Layer 2: g3==2 vs g3==1
# ======================================================================
echo ""
echo ">>> [3/4] MultiNicheNet Layer 2: g3 comparison (anno2)"
echo "    Started: $(date)"

Rscript "$MNN_SCRIPT" \
  -i "$INPUT_L2" \
  -s patient_name \
  -g g3 \
  -c anno2 \
  -f "X2-X1" \
  --logfc_threshold 0.10 \
  --min_cells 5 \
  -o "${BASE}/cci/mnn/L2_g3_anno2" \
  --cores 4 \
  2>&1 | tee "${LOG_DIR}/mnn_L2_g3_anno2.log"

echo "    >>> MNN L2 anno2 completed: $(date)"

# ======================================================================
# MULTINICHENET - Layer 1: HC vs IS
# ======================================================================
echo ""
echo ">>> [4/4] MultiNicheNet Layer 1: HC vs IS (anno2)"
echo "    Started: $(date)"

Rscript "$MNN_SCRIPT" \
  -i "$INPUT_L1" \
  -s patient_name \
  -g cohort \
  -c anno2 \
  -f "Stroke-HC" \
  --logfc_threshold 0.10 \
  --min_cells 5 \
  -o "${BASE}/cci/mnn/L1_cohort_anno2" \
  --cores 4 \
  2>&1 | tee "${LOG_DIR}/mnn_L1_cohort_anno2.log"

echo "    >>> MNN L1 anno2 completed: $(date)"

echo ""
echo "=== CCI Analysis Finished: $(date) ==="
