#!/bin/bash
# Formula 1 전체 분석 상태 확인 스크립트

echo "=========================================="
echo "Formula 1 전체 분석 상태 확인"
echo "=========================================="
echo ""

# 프로세스 확인
echo "1. 실행 중인 프로세스:"
ps aux | grep -E "Rscript.*run_formula1_analysis" | grep -v grep || echo "  실행 중인 프로세스 없음"
echo ""

# 로그 파일 확인
LOG_FILE="/tmp/formula1_full_analysis.log"
if [ -f "$LOG_FILE" ]; then
    echo "2. 로그 파일: $LOG_FILE"
    echo "   크기: $(du -h $LOG_FILE | cut -f1)"
    echo "   줄 수: $(wc -l < $LOG_FILE)"
    echo ""
    
    echo "3. 최근 로그 (마지막 30줄):"
    tail -30 "$LOG_FILE"
    echo ""
    
    echo "4. 주요 이벤트:"
    grep -E "데이터 로드|메타데이터|Formula|분석 시작|runNEBULA|완료|실패|Error|error|결과 저장|분석 완료" "$LOG_FILE" | tail -20 || echo "  이벤트 없음"
    echo ""
else
    echo "2. 로그 파일이 아직 생성되지 않았습니다."
    echo ""
fi

# 결과 파일 확인
echo "5. 생성된 결과 파일:"
ls -lht /data/user3/sobj/IS6_formula1*.qs 2>/dev/null | head -5 || echo "  결과 파일 아직 생성되지 않음"
echo ""

echo "=========================================="
echo "상태 확인 완료"
echo "=========================================="

