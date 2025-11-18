# 개발 컨텍스트 요약 (2025-11-12)

## 프로젝트 개요
- **프로젝트 경로**: `/home/user3/data_user3/git_repo/mylit`
- **R 패키지**: `myR` (로컬 개발 패키지)
- **작업 디렉토리**: `/home/user3/GJC_KDW_250721`

## 최근 완료된 작업

### 1. Find Gene Signature (FGS) 함수 개발
- **`find_gene_signature_v5.2`**: 완성됨
  - 모든 메서드 구현 완료 (random_forest, ranger, lasso, ridge, elastic_net, xgboost, pca_loadings, nmf_loadings, gam, limma, wilcoxon)
  - `mgcv::bam` 사용으로 속도 개선
  - `limma` 레벨 sanitization 추가 (make.names)
  - 위치: `myR/R/test.R` (line ~3106-3885)

- **`find_gene_signature_v5.3`**: 완성됨
  - v5.2 기반, GAM에 동적 k 적용
  - `gam.k = NULL`이면 `k_dynamic = max(3, min(10, floor(n_unique_vals / gam.k_dynamic_factor)))`
  - `gam.min_unique` threshold로 유전자 필터링
  - 위치: `myR/R/test.R` (line ~3888-4119)

- **FGS_v5.2, FGS_v5.3 별칭 추가**
  - `myR/R/signature.R`에 별칭 함수 추가 (line ~1784-1820)
  - NAMESPACE에 export 추가

### 2. 테스트 결과
- **v5.3 GAM 동적 k**: 정상 작동 확인 (다운샘플링 데이터)
- **pseudobulk 호환성**: v5.2, v5.3 모두 pseudobulk 데이터와 호환 확인
- **테스트 데이터**: `/data/user3/sobj/IS_scvi_251107_ds2500.qs` (다운샘플링 버전)

### 3. Differential Expression 분석 함수
- **`runMUSCAT`**: 정상 작동 확인 ✅
  - 위치: `myR/R/test_analysis.R`
  - pseudobulking 자동 수행

- **`runMAST_v1`**: MAST 패키지 호환성 문제 ⚠️
  - 위치: `myR/R/test_analysis.R` (line ~40-105)
  - 문제: MAST 1.28.0에서 `FromSeurat` deprecated
  - 시도: `MAST::SceToSingleCellAssay()` 사용, 하지만 내부적으로 `FromSeurat` 호출
  - 상태: 수정 중, `FromMatrix` fallback 추가했으나 여전히 오류

- **`runNEBULA_v1`**: 수정 완료, 최종 테스트 필요
  - 위치: `myR/R/test_analysis.R` (line ~107-275)
  - 수정: HL method 추가, fallback 구현
  - 상태: 코드 수정 완료, 테스트 필요

## 현재 상태

### 완료된 것
1. ✅ FGS_v5.2, v5.3 구현 및 테스트
2. ✅ FGS_v5.2, v5.3을 signature.R에 추가
3. ✅ runMUSCAT 정상 작동 확인
4. ✅ v5.3 GAM 동적 k 기능 확인
5. ✅ pseudobulk 호환성 확인

### 진행 중/이슈
1. ⚠️ **runMAST_v1**: MAST 패키지 호환성 문제
   - MAST 1.28.0에서 `FromSeurat` 제거됨
   - `SceToSingleCellAssay`가 내부적으로 `FromSeurat` 호출
   - 해결 방안 필요

2. 🔄 **runNEBULA_v1**: 코드 수정 완료, 최종 테스트 필요

## 주요 파일 위치

### 핵심 함수
- `myR/R/test.R`: `find_gene_signature_v5.2`, `find_gene_signature_v5.3`, `TML6`
- `myR/R/signature.R`: `FGS_v5.2`, `FGS_v5.3` (별칭)
- `myR/R/test_analysis.R`: `runMAST_v1`, `runNEBULA_v1`, `runMUSCAT`

### 헬퍼 함수
- `fgs_preprocess_data_v5.2`: `myR/R/test.R` (line ~2939)

### 설정 파일
- `myR/NAMESPACE`: export 목록
- `st/start.R`: R 세션 초기화 스크립트

## 테스트 데이터
- **다운샘플링 데이터**: `/data/user3/sobj/IS_scvi_251107_ds2500.qs`
- **원본 데이터**: `/data/user3/sobj/IS_scvi_251107.qs`

## 다음 작업 제안

1. **runMAST_v1 수정**
   - MAST 패키지 최신 버전 문서 확인
   - `FromMatrix` 직접 사용 또는 다른 변환 방법 시도
   - 또는 MAST 대신 다른 DE 방법 제안

2. **runNEBULA_v1 최종 테스트**
   - 다운샘플링 데이터로 테스트
   - 성공 시 full 데이터로 테스트

3. **전체 파이프라인 통합 테스트**
   - FGS → TML6 → DE 분석 전체 플로우 확인

## 참고사항
- 모든 테스트는 `qs::qsave()`로 결과 저장
- 실행 시간 기록 권장
- 다운샘플링 데이터로 먼저 테스트 후 full 데이터 사용

