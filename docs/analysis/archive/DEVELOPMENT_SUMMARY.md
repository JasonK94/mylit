# 개발 요약: MUSCAT2_v1, NEBULA2_v1, runNEBULA (pseudobulk mode)

## 개요
main2 워크트리에서 runMAST 시리즈 함수 개발을 진행했습니다. MAST2_v1은 문제가 많아 제외하고, MUSCAT2_v1, NEBULA2_v1, runNEBULA (pseudobulk mode) 함수를 개발했습니다.

## 개발된 함수

### 1. runMUSCAT
**위치**: `myR/R/test_analysis.R` (line 528-745)

**기능**:
- MUSCAT 기반 pseudobulk differential expression 분석
- runMUSCAT (legacy mode) 개선 버전
- 결측치 처리 강화 (g3 등)
- NA 필터링 개선
- 오류 메시지 개선

**주요 개선사항**:
- `remove_na_groups` 파라미터 추가 (기본값: TRUE)
- 필수 컬럼 검증
- NA 값이 있는 셀 자동 제거
- 그룹 수 검증 (최소 2개 그룹 필요)
- 단계별 상세 메시지 출력

**파라미터**:
- `sobj`: Seurat 객체
- `cluster_id`: 클러스터 컬럼명 (기본: "seurat_clusters")
- `sample_id`: 샘플 ID 컬럼명 (기본: "hos_no")
- `group_id`: 그룹 컬럼명 (기본: "type")
- `batch_id`: 배치 컬럼명 (선택)
- `contrast`: 비교할 그룹 (예: "IS - SAH")
- `method`: DE 방법 ("edgeR", "DESeq2", "limma-trend", "limma-voom")
- `remove_na_groups`: NA 값이 있는 셀 제거 여부 (기본: TRUE)

### 2. runNEBULA
**위치**: `myR/R/test_analysis.R` (line 779-934)

**기능**:
- NEBULA 기반 single-cell level differential expression 분석
- runNEBULA_v1 개선 버전
- 결측치 처리 강화 (g3 등)
- 변수별 NA 개수 리포트
- factor 레벨 검증

**주요 개선사항**:
- `remove_na_cells` 파라미터 추가 (기본값: TRUE)
- 변수별 NA 개수 리포트
- factor 변수 레벨 수 확인
- offset 값 검증 (양수 확인)
- 더 상세한 오류 메시지

**파라미터**:
- `sobj`: Seurat 객체
- `layer`: Assay layer (기본: "counts")
- `fixed_effects`: 고정 효과 변수 (기본: c("g3"))
- `covar_effects`: 공변량 변수 (선택)
- `patient_col`: 환자 ID 컬럼명 (기본: "hos_no")
- `offset`: Offset 컬럼명 (기본: "nCount_RNA")
- `min_count`: 최소 발현 세포 수 (기본: 10)
- `remove_na_cells`: NA 값이 있는 셀 제거 여부 (기본: TRUE)

### 3. runNEBULA (pseudobulk mode)
**위치**: `myR/R/test_analysis.R` (line 977-1206)

**기능**:
- Pseudobulk와 NEBULA 분석 수직 결합
- 클러스터-샘플별 집계 후 NEBULA 실행
- Pseudobulk 메타데이터 및 count matrix 반환

**주요 특징**:
- Pseudobulking: cluster_id와 sample_id로 집계
- Offset 계산 방법: sum, mean, n_cells
- Pseudobulk 메타데이터 자동 생성
- NEBULA 분석 결과와 함께 pseudobulk 데이터 반환

**파라미터**:
- `sobj`: Seurat 객체
- `layer`: Assay layer (기본: "counts")
- `cluster_id`: 클러스터 컬럼명 (기본: "seurat_clusters")
- `sample_id`: 샘플 ID 컬럼명 (기본: "hos_no")
- `group_id`: 그룹 컬럼명 (기본: "type")
- `fixed_effects`: 고정 효과 변수 (기본: c("g3"))
- `covar_effects`: 공변량 변수 (선택)
- `patient_col`: 환자 ID 컬럼명 (기본: "hos_no")
- `offset_method`: Offset 계산 방법 ("sum", "mean", "n_cells")
- `min_count`: 최소 발현 세포 수 (기본: 10)
- `min_cells_per_pb`: Pseudobulk 샘플당 최소 세포 수 (기본: 3)
- `remove_na_cells`: NA 값이 있는 셀 제거 여부 (기본: TRUE)
- `keep_clusters`: 특정 클러스터만 사용 (선택)

**반환값**:
- `nebula_result`: NEBULA 분석 결과
- `pseudobulk_meta`: Pseudobulk 메타데이터
- `pseudobulk_counts`: Pseudobulk count matrix

## 파일 변경사항

### 수정된 파일
1. **myR/R/test_analysis.R**
   - 725줄 추가
   - runMUSCAT 함수 추가 (217줄)
   - runNEBULA 함수 추가 (155줄)
   - runNEBULA (pseudobulk mode) 함수 추가 (229줄)

2. **myR/NAMESPACE**
   - export 추가:
     - `runMUSCAT`
     - `runNEBULA`
     - `runNEBULA (pseudobulk mode)`

### 새로 생성된 파일
1. **test_functions.R**
   - 실제 데이터로 테스트하는 스크립트
   - 각 함수별 테스트 플래그 포함

2. **test_functions_quick.R**
   - 함수 존재 확인 및 기본 검증 스크립트
   - 실제 데이터 없이 함수 시그니처 확인

## 주요 개선사항

### 1. 결측치 처리
- g3 등 메타데이터의 NA 값 자동 감지 및 제거
- 변수별 NA 개수 리포트
- 제거 여부 선택 가능 (`remove_na_cells`, `remove_na_groups`)

### 2. 데이터 검증
- 필수 컬럼 존재 확인
- 그룹 수 검증 (최소 2개 그룹)
- factor 레벨 수 확인
- offset 값 검증 (양수 확인)

### 3. 오류 처리
- 상세한 오류 메시지
- 단계별 진행 메시지
- 변수별 NA 개수 리포트

### 4. 로깅
- 단계별 상세 메시지 출력
- 진행 상황 표시
- 결과 요약 정보

## 테스트 방법

### 1. 빠른 테스트 (함수 존재 확인)
```r
source("test_functions_quick.R")
```

### 2. 실제 데이터 테스트
```r
# test_functions.R 수정
test_muscat2 <- TRUE  # 테스트 플래그 활성화
test_nebula2 <- TRUE
test_nebula2_pb <- TRUE

# 스크립트 실행
source("test_functions.R")
```

### 3. 개별 함수 테스트
```r
# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# runMUSCAT 테스트
res_muscat2 <- runMUSCAT(
  sobj = sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  contrast = "IS - SAH",
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만
)

# runNEBULA 테스트
res_nebula2 <- runNEBULA(
  sobj = sobj,
  layer = "counts",
  fixed_effects = "g3",
  patient_col = "hos_no",
  offset = "nCount_RNA",
  remove_na_cells = TRUE
)

# runNEBULA (pseudobulk mode) 테스트
res_nebula2_pb <- runNEBULA (pseudobulk mode)(
  sobj = sobj,
  layer = "counts",
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  fixed_effects = "g3",
  patient_col = "hos_no",
  offset_method = "sum",
  remove_na_cells = TRUE,
  keep_clusters = c("0", "1", "2")
)
```

## Git 커밋 정보

- **브랜치**: main2
- **커밋 해시**: f51d825
- **커밋 메시지**: "Add MUSCAT2_v1, NEBULA2_v1, and runNEBULA (pseudobulk mode) functions"

## 다음 단계

1. **실제 데이터로 테스트**
   - 다운샘플링 데이터로 먼저 테스트
   - 성공 시 full 데이터로 확장

2. **성능 최적화**
   - 큰 데이터셋에서의 성능 확인
   - 필요시 최적화

3. **문서화**
   - 함수 도움말 완성
   - 사용 예제 추가

4. **통합 테스트**
   - 다른 함수와의 통합 테스트
   - 전체 파이프라인 테스트

## 주의사항

1. **g3 결측치**: g3 변수에 결측치가 있을 수 있으므로 `remove_na_cells=TRUE` 사용 권장
2. **메모리**: 큰 데이터셋의 경우 메모리 사용량 주의
3. **실행 시간**: NEBULA 분석은 시간이 오래 걸릴 수 있음
4. **패키지 의존성**: 필요한 패키지들이 설치되어 있어야 함
   - muscat
   - nebula
   - SingleCellExperiment
   - SummarizedExperiment
   - Seurat
   - 기타

## 참고사항

- MAST2_v1 함수는 개발하지 않음 (문제가 많음)
- 모든 함수는 main2 워크트리에서 개발됨
- 테스트는 다운샘플링 데이터로 먼저 수행 권장
- 결과는 qs::qsave()로 저장 가능

