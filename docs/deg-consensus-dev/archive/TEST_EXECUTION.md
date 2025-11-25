# Phase 2 테스트 실행 가이드

## 테스트 실행 방법

### 1. R 세션 시작
```bash
cd /home/user3/GJC_KDW_250721
R
```

### 2. R 세션에서 실행
```r
# 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")

# 테스트 스크립트 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_phase2_limma_interactive.R")
```

## 예상 결과

테스트가 성공하면 다음 파일들이 생성됩니다:
- `/data/user3/sobj/test_limma_voom_v1_result.qs`
- `/data/user3/sobj/test_limma_trend_v1_result.qs`
- `/data/user3/sobj/test_deg_consensus_phase2_limma.qs`

## 디버깅

오류가 발생하면:
1. traceback()으로 스택 트레이스 확인
2. 각 단계별로 데이터 구조 확인
3. 메타데이터 컬럼명 확인

## 알려진 이슈

1. **pb_sub$group_id 참조 문제**: 수정 완료 - 각 클러스터별로 샘플의 group_id를 올바르게 추출하도록 수정

