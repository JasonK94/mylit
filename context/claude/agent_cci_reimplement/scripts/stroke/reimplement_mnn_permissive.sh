#!/bin/bash
# ==============================================================================
# MNN Re-run with permissive thresholds
# ==============================================================================
# Problem: group_prioritization_tbl = 0 rows for all previous runs
# Root cause: lr_target_prior_cor is NULL (insufficient samples per cell type)
#
# Strategy: Relax all thresholds significantly
#   - min_cells: 10 → 5 (allow sparser cell types)
#   - logfc_threshold: 0.10 → 0.05 (weaker DE included)
#   - fraction_cutoff: 0.05 → 0.01 (express in fewer cells)
#   - p_val_threshold: 0.05 → 0.10 (more lenient significance)
#
# Run at anno2 level (8 types) where more samples available per cell type
# ==============================================================================

set -e

BASE="/data/user3/sobj/stroke_hc_v8_2"
MNN_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R"

echo "=== MNN Permissive Re-run ==="
echo "Started: $(date)"

# ---- Layer 1: HC vs Stroke at anno2 (permissive) ----
echo ""
echo ">>> [1/4] MNN L1 anno2 (permissive)"
mkdir -p "$BASE/cci/mnn/L1_cohort_anno2_permissive"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_1_hc_is.qs" \
  -o "$BASE/cci/mnn/L1_cohort_anno2_permissive" \
  -g cohort \
  -s patient_name \
  -c anno2 \
  -f "Stroke-HC" \
  --min_cells 5 \
  --logfc_threshold 0.05 \
  --fraction_cutoff 0.01 \
  --p_val_threshold 0.10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L1_anno2_permissive.log"

echo "    L1 anno2 permissive done: $(date)"

# ---- Layer 2: g3 at anno2 (permissive) ----
echo ""
echo ">>> [2/4] MNN L2 anno2 (permissive)"
mkdir -p "$BASE/cci/mnn/L2_g3_anno2_permissive"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_2_is_g3.qs" \
  -o "$BASE/cci/mnn/L2_g3_anno2_permissive" \
  -g g3 \
  -s patient_name \
  -c anno2 \
  -f "X2-X1" \
  --min_cells 5 \
  --logfc_threshold 0.05 \
  --fraction_cutoff 0.01 \
  --p_val_threshold 0.10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L2_anno2_permissive.log"

echo "    L2 anno2 permissive done: $(date)"

# ---- Layer 1: HC vs Stroke at anno1 (permissive) ----
echo ""
echo ">>> [3/4] MNN L1 anno1 (permissive)"
mkdir -p "$BASE/cci/mnn/L1_cohort_anno1_permissive"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_1_hc_is.qs" \
  -o "$BASE/cci/mnn/L1_cohort_anno1_permissive" \
  -g cohort \
  -s patient_name \
  -c anno1 \
  -f "Stroke-HC" \
  --min_cells 5 \
  --logfc_threshold 0.05 \
  --fraction_cutoff 0.01 \
  --p_val_threshold 0.10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L1_anno1_permissive.log"

echo "    L1 anno1 permissive done: $(date)"

# ---- Layer 2: g3 at anno1 (permissive) ----
echo ""
echo ">>> [4/4] MNN L2 anno1 (permissive)"
mkdir -p "$BASE/cci/mnn/L2_g3_anno1_permissive"
taskset -c 16-31 Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_2_is_g3.qs" \
  -o "$BASE/cci/mnn/L2_g3_anno1_permissive" \
  -g g3 \
  -s patient_name \
  -c anno1 \
  -f "X2-X1" \
  --min_cells 5 \
  --logfc_threshold 0.05 \
  --fraction_cutoff 0.01 \
  --p_val_threshold 0.10 \
  --cores 16 \
  2>&1 | tee "$BASE/logs/mnn_L2_anno1_permissive.log"

echo "    L2 anno1 permissive done: $(date)"

echo ""
echo "=== MNN Permissive Re-run Complete: $(date) ==="

# Check results
for d in L1_cohort_anno2_permissive L2_g3_anno2_permissive L1_cohort_anno1_permissive L2_g3_anno1_permissive; do
  rds="$BASE/cci/mnn/$d/multinichenet_results.rds"
  qs="$BASE/cci/mnn/$d/multinichenet_results.qs"
  if [ -f "$rds" ] || [ -f "$qs" ]; then
    echo "  $d: results file exists"
  else
    echo "  $d: NO results found"
  fi
done
