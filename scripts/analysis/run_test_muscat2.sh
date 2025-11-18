#!/bin/bash
# ============================================================================
# runMUSCAT 테스트 실행 스크립트
# ============================================================================

# 작업 디렉터리로 이동
cd /home/user3/GJC_KDW_250721

# R 스크립트 실행
Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_run_muscat2.R 2>&1 | tee /tmp/test_muscat2.log

# 결과 확인
if [ $? -eq 0 ]; then
    echo "========================================"
    echo "테스트 성공"
    echo "========================================"
    echo "로그 파일: /tmp/test_muscat2.log"
    echo "결과 파일: /data/user3/sobj/test_muscat2_v1_result.qs"
else
    echo "========================================"
    echo "테스트 실패"
    echo "========================================"
    echo "로그 파일: /tmp/test_muscat2.log"
    echo "오류 내용을 확인하세요."
fi

