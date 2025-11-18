# 테스트 및 분석 실행 가이드

## 빠른 시작

### 1. R 세션 시작
```bash
cd /home/user3/GJC_KDW_250721
R
```

### 2. 테스트 실행
```r
# 대화형 테스트 스크립트 실행 (권장)
source("/home/user3/data_user3/git_repo/_wt/analysis/test_interactive.R")
```

## 개발된 함수

1. **runMUSCAT2_v1**: MUSCAT 분석 함수 (결측치 처리 개선)
2. **runNEBULA2_v1**: NEBULA 분석 함수 (결측치 처리 개선)
3. **runNEBULA2_v1_with_pseudobulk**: Pseudobulk와 NEBULA 결합 함수

## 테스트 스크립트

- `test_interactive.R`: 대화형 테스트 (권장)
- `test_simple.R`: 기본 데이터 확인
- `test_run_muscat2.R`: runMUSCAT2_v1 전용 테스트
- `test_functions.R`: 전체 함수 테스트

## 문서

- `QUICK_START.md`: 빠른 시작 가이드
- `TEST_INSTRUCTIONS.md`: 상세한 테스트 지침
- `TEST_SUMMARY.md`: 테스트 요약
- `DEVELOPMENT_SUMMARY.md`: 개발 요약

## 다음 단계

1. R 세션에서 `test_interactive.R` 실행
2. 결과 확인 및 분석
3. 전체 데이터로 확장 테스트
