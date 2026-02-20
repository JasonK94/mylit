#!/bin/bash
# ==============================================================================
# FGS & TML Sweep: Run for multiple n_features values
# Input: 5_2_is_g3.qs (IS patients with valid g3)
# ==============================================================================
set -e

INPUT="/data/user3/sobj/stroke_hc_v8_2/5_2_is_g3.qs"
OUTPUT_BASE="/data/user3/sobj/stroke_hc_v8_2/fgs"
FGS_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R"
LOG_DIR="/data/user3/sobj/stroke_hc_v8_2/logs"

mkdir -p "$OUTPUT_BASE" "$LOG_DIR"

# Verify input exists
if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

echo "=== FGS Sweep Started: $(date) ==="
echo "Input: $INPUT"
echo "Output base: $OUTPUT_BASE"

# n_features values to sweep
FEATURES=(50 100 200)

for N in "${FEATURES[@]}"; do
    OUTPUT="${OUTPUT_BASE}/whole_IS_${N}"
    LOG="${LOG_DIR}/fgs_whole_IS_${N}_$(date +%Y%m%d_%H%M%S).log"
    
    echo ""
    echo ">>> Running FGS n_features=${N}"
    echo "    Output: ${OUTPUT}"
    echo "    Log: ${LOG}"
    echo "    Started: $(date)"
    
    # CPU limit to prevent resource exhaustion (use 16 cores)
    taskset -c 0-15 Rscript "$FGS_SCRIPT" \
        --input "$INPUT" \
        --output "$OUTPUT" \
        --target_var g3 \
        --control_vars "sex,age,GEM" \
        --n_features "$N" \
        --cores 16 \
        2>&1 | tee "$LOG"
    
    EXIT_CODE=${PIPESTATUS[0]}
    if [ $EXIT_CODE -eq 0 ]; then
        echo "    >>> COMPLETED n_features=${N}: $(date)"
    else
        echo "    >>> FAILED n_features=${N} (exit code: $EXIT_CODE): $(date)"
        echo "    Check log: $LOG"
    fi
done

echo ""
echo "=== FGS Sweep Finished: $(date) ==="
