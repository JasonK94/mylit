# Pseudotime Analysis 테스트 결과

## 테스트 환경
- 작업 디렉토리: `/home/user3/GJC_KDW_250721`
- R 세션: renv 및 start.R 자동 로드
- 테스트 데이터: `/data/user3/sobj/IS_scvi_251107_ds2500.qs` (다운샘플)
- 원본 데이터: `/data/user3/sobj/IS_scvi_251107.qs` (전체)

## 완료된 테스트

### 1. 데이터 전처리 ✅
- **스크립트**: `test_pseudotime_preprocessing.R`
- **결과**:
  - 다운샘플 데이터 전처리 완료
  - 원본 데이터 전처리 완료
  - sex 변수: M, F로 변환 (인코딩 문제 해결)
  - g3 변수: factor로 변환
  - 시간 변수 파싱 완료
  - 시간 차이 계산 완료

**생성된 파일**:
- `/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs` (321 MB)
- `/data/user3/sobj/IS_scvi_251107_preprocessed.qs` (3.7 GB)

### 2. 기본 Pseudotime 분석 ✅
- **스크립트**: `test_pseudotime_basic.R`
- **결과**:
  - Monocle3 trajectory inference: ✅ 성공
    - 2500 cells, 51795 genes
    - Pseudotime 계산 완료 (root cell 지정 후)
  - Slingshot trajectory inference: ✅ 성공 (UMAP 사용)
    - UMAP reduction 사용 시 정상 작동
    - PCA 사용 시 singular matrix 오류 (예상됨)

**생성된 파일**:
- `/data/user3/sobj/IS_scvi_251107_ds2500_monocle3.qs` (18 MB)
- `/data/user3/sobj/IS_scvi_251107_ds2500_slingshot.qs`

### 3. 고급 분석 ✅
- **스크립트**: `test_pseudotime_advanced.R`
- **결과**:
  - Metadata 상관관계 동적 분석: ✅ 성공
    - DDIT4, UTY, S100B 유전자 분석
    - nih_change, age, arrival_gcs_score와의 상관관계 계산
    - 구간별 상관관계 변화 분석 완료
  - 패턴 탐지: ✅ 성공
    - Bimodal 패턴 탐지 성공
    - 증가/감소, oscillatory 패턴 분석 가능

**생성된 파일**:
- `/data/user3/sobj/IS_scvi_251107_ds2500_correlation_dynamics.qs`
- `/data/user3/sobj/IS_scvi_251107_ds2500_pattern_detection.qs`

## 발견된 이슈 및 해결

### 1. 인코딩 문제 ✅ 해결
- **문제**: sex 변수의 인코딩이 깨져서 표시됨
- **해결**: iconv를 이용한 안전한 변환 로직 추가
- **상태**: 해결 완료

### 2. Slingshot na 인자 호환성 ✅ 해결
- **문제**: `slingPseudotime(na = FALSE)` 인자가 일부 버전에서 지원되지 않음
- **해결**: tryCatch로 fallback 처리
- **상태**: 해결 완료

### 3. Monocle3 Pseudotime 계산 ✅ 해결
- **문제**: root cell이 지정되지 않아 pseudotime이 계산되지 않음
- **해결**: Seurat 객체의 클러스터 정보를 사용하여 root cell 자동 선택
- **상태**: 해결 완료

### 4. PCA vs UMAP
- **발견**: Slingshot은 PCA보다 UMAP에서 더 안정적으로 작동
- **권장**: UMAP reduction 사용

## 주요 Feature 확인

모든 주요 유전자가 데이터에 존재함:
- ✅ DDIT4, UTY, S100B, XIST, HLA-B, CCL4, HLA-C, TXNIP (8/8)

## Metadata 확인

중요한 metadata 모두 확인됨:
- ✅ nih_change (numeric)
- ✅ sex (factor: M, F) - 전처리 완료
- ✅ g3 (factor: 1, 2) - 전처리 완료
- ✅ arrival_gcs_score, age, ini_nih, nih_end1 (numeric)
- ✅ icu_adm_dt, ia_start, ia_end (datetime) - 전처리 완료
- ✅ 시간 차이 변수들 (icu_to_ia_start_hours 등) - 계산 완료

## 다음 단계

1. ✅ 전처리 완료
2. ✅ 기본 분석 완료
3. ✅ 고급 분석 완료
4. ⏳ 원본 데이터로 최종 테스트 (선택적)
5. ⏳ 추가 함수 개발 (Branching Analysis 등)

## 실행 방법

```r
# 작업 디렉토리에서 실행
cd /home/user3/GJC_KDW_250721

# 1. 전처리
Rscript test_pseudotime_preprocessing.R

# 2. 기본 분석
Rscript test_pseudotime_basic.R

# 3. 고급 분석
Rscript test_pseudotime_advanced.R

# 4. 결과 요약
Rscript test_results_summary.R
```

