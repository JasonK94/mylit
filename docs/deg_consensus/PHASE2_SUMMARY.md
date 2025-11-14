# Phase 2 진행 상황

## 완료된 작업

### 1. limma 계열 방법론 구현
- **파일**: `myR/R/deg_consensus/deg_methods_limma.R`
- **구현된 함수**:
  - `runLIMMA_voom_v1()`: limma-voom 분석
  - `runLIMMA_trend_v1()`: limma-trend 분석

### 2. 주요 특징
- runMUSCAT2_v1과 동일한 인터페이스
- 동일한 전처리 파이프라인 (NA 처리, pseudobulk 생성 등)
- 클러스터별 분석 지원
- 표준화된 출력 형식

### 3. 테스트 스크립트
- **파일**: `scripts/test_phase2_limma.R`
- limma-voom, limma-trend 직접 테스트
- run_deg_consensus 통합 테스트

## 다음 단계

### 1. 테스트 실행
```r
# R 세션에서 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_phase2_limma.R")
```

### 2. 추가 구현 필요
- `runLIMMA_wt_v1()`: limma with weights
- `runDREAM_v1()`: dream 분석
- edgeR 독립 구현 (LRT, QLF, robust)
- DESeq2 독립 구현 (Wald, LRT)

## 테스트 결과 확인

테스트 실행 후 다음 파일들이 생성되어야 함:
- `/data/user3/sobj/test_limma_voom_v1_result.qs`
- `/data/user3/sobj/test_limma_trend_v1_result.qs`
- `/data/user3/sobj/test_deg_consensus_phase2_limma.qs`

## 주의사항

1. **함수 로드**: 테스트 전에 필요한 함수들이 로드되어 있어야 함
2. **패키지**: limma, edgeR, muscat 등 필수 패키지 설치 확인
3. **메모리**: 큰 데이터셋의 경우 메모리 사용량 주의

