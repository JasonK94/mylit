#!/bin/bash
# Test all L2 methods individually using taskset
# This script tests each L2 method one by one to avoid CPU saturation

# CPU 코어 제한 (taskset 사용)
# 기본값: 8코어 (0-7)
CORES="${1:-0-7}"

echo "========================================"
echo "Testing all L2 methods individually"
echo "Using CPU cores: $CORES"
echo "========================================"
echo ""

# 지원되는 모든 methods
METHODS=("glm" "ranger" "xgbTree" "glmnet" "svmRadial" "mlp" "mlpKerasDropout" "nnet" "earth")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_SCRIPT="$SCRIPT_DIR/test_l2_method_individual.R"

if [ ! -f "$TEST_SCRIPT" ]; then
    echo "Error: Test script not found: $TEST_SCRIPT"
    exit 1
fi

# 결과 저장 디렉토리
RESULTS_DIR="/data/user3/sobj"
SUMMARY_FILE="$RESULTS_DIR/tml2_v7_test_summary.txt"

echo "Results will be saved to: $RESULTS_DIR"
echo "Summary will be saved to: $SUMMARY_FILE"
echo ""

# Summary 파일 초기화
echo "L2 Method Test Summary" > "$SUMMARY_FILE"
echo "Generated: $(date)" >> "$SUMMARY_FILE"
echo "CPU cores used: $CORES" >> "$SUMMARY_FILE"
echo "========================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

SUCCESS_COUNT=0
FAIL_COUNT=0

for METHOD in "${METHODS[@]}"; do
    echo "----------------------------------------"
    echo "Testing: $METHOD"
    echo "----------------------------------------"
    
    # taskset을 사용하여 CPU 코어 제한
    if taskset -c "$CORES" Rscript "$TEST_SCRIPT" "$METHOD"; then
        echo "✓ $METHOD: SUCCESS"
        echo "✓ $METHOD: SUCCESS" >> "$SUMMARY_FILE"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ $METHOD: FAILED"
        echo "✗ $METHOD: FAILED" >> "$SUMMARY_FILE"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
    sleep 2  # 각 테스트 사이에 짧은 대기
done

# Summary 작성
echo "" >> "$SUMMARY_FILE"
echo "========================================" >> "$SUMMARY_FILE"
echo "Total tested: ${#METHODS[@]}"
echo "Successful: $SUCCESS_COUNT"
echo "Failed: $FAIL_COUNT"
echo "" >> "$SUMMARY_FILE"
echo "Total tested: ${#METHODS[@]}" >> "$SUMMARY_FILE"
echo "Successful: $SUCCESS_COUNT" >> "$SUMMARY_FILE"
echo "Failed: $FAIL_COUNT" >> "$SUMMARY_FILE"

echo ""
echo "========================================"
echo "All tests completed"
echo "Summary saved to: $SUMMARY_FILE"
echo "========================================"

