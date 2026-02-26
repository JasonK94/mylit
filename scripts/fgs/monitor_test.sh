#!/bin/bash
LOG_FILE="logs/fgs/test_real_data_v14.log"
CMD="nohup taskset -c 0-7 Rscript scripts/fgs/test_real_data.R > $LOG_FILE 2>&1 &"

echo "Starting test v14..."
eval $CMD
PID=$!
echo "Test started with PID $PID"

while kill -0 $PID 2>/dev/null; do
    echo "Test running... $(date)"
    sleep 1800 # 30분 대기
done

echo "Test finished at $(date)"
# 결과 요약 저장
tail -n 100 $LOG_FILE > logs/fgs/latest_summary.txt
