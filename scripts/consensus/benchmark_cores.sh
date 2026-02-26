#!/bin/bash
# scripts/consensus/benchmark_cores.sh

LOG_DIR="benchmark_logs"
mkdir -p "$LOG_DIR"

INPUT="/data/user3/sobj/is2_IS_3_clustered.qs"
BASE_OUT="/data/user3/sobj/consensus/AG_benchmark"
SCRIPT="scripts/consensus/run_deg_consensus_cli.R"

# Define cores to test
CORES=(1 2 4 8)

echo "Starting Benchmark: $(date)"

for c in "${CORES[@]}"; do
    OUT_DIR="${BASE_OUT}/cores_${c}"
    LOG_FILE="${LOG_DIR}/run_core_${c}.log"
    
    echo "Launching run with ${c} cores -> ${OUT_DIR}"
    
    # Run in background
    (
        /usr/bin/time -v Rscript "$SCRIPT" \
            --input "$INPUT" \
            --output "$OUT_DIR" \
            --cluster "anno3big" \
            --group "g3" \
            --contrast "2 - 1" \
            --covariates "sex" \
            --methods "edgeR-LRT,DESeq2-Wald,limma-voom,nebula" \
            --cores "$c" \
            > "$LOG_FILE" 2>&1
    ) &
done

wait
echo "All benchmark runs submitted/finished. Check output timestamps."
