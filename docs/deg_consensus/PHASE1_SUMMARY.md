# Phase 1 완료 요약

## 완료된 작업

### 1. 메인 통합 함수 구현
- **파일**: `myR/R/deg_consensus/run_deg_consensus.R`
- **함수**: `run_deg_consensus()`

### 2. 주요 기능

#### 2.1 지원하는 방법론
- **muscat 기반** (기존 함수 활용):
  - `muscat-edgeR`: `runMUSCAT2_v1(method="edgeR")`
  - `muscat-DESeq2`: `runMUSCAT2_v1(method="DESeq2")`
  - `muscat-limma-voom`: `runMUSCAT2_v1(method="limma-voom")`
  - `muscat-limma-trend`: `runMUSCAT2_v1(method="limma-trend")`
- **nebula** (기존 함수 활용):
  - `nebula`: `runNEBULA2_v1()`
- **Phase 2에서 구현 예정**:
  - `limma-voom`, `limma-trend`, `limma-wt`, `dream`
  - `edgeR-LRT`, `edgeR-QLF`, `edgeR-robust`
  - `DESeq2-Wald`, `DESeq2-LRT`

#### 2.2 함수 인터페이스
```r
run_deg_consensus(
  sobj,
  contrast = "2 - 1",
  methods = c("muscat-edgeR", "muscat-DESeq2", "nebula"),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM",
  pb_min_cells = 3,
  filter_genes = "edgeR",
  keep_clusters = NULL,
  remove_na_groups = TRUE,
  parallel = FALSE,
  n_cores = 4,
  verbose = TRUE
)
```

#### 2.3 반환값
```r
list(
  results = list(),        # 각 방법론별 결과
  methods_run = character(), # 성공한 방법론
  methods_failed = character(), # 실패한 방법론
  errors = list(),         # 에러 메시지
  contrast = character(), # 사용된 contrast
  metadata = list()        # 메타데이터 정보
)
```

### 3. 구현 세부사항

#### 3.1 방법론 핸들러 매핑
- 각 방법론을 해당 함수로 매핑하는 핸들러 시스템
- Phase 2에서 구현될 함수들도 미리 등록되어 있음

#### 3.2 인자 자동 조정
- **muscat 방법론**: `cluster_id` 필수 확인
- **nebula 방법론**: 
  - `fixed_effects` = `group_id`로 자동 설정
  - `patient_col` = `sample_id`로 자동 설정
  - muscat 전용 인자 제거

#### 3.3 에러 처리
- 각 방법론별 try-catch로 독립 실행
- 실패한 방법론은 건너뛰고 계속 진행
- 에러 메시지 수집 및 반환

### 4. 다음 단계 (Phase 2)

#### 4.1 limma 계열 표준화
- `runLIMMA_voom_v1()`: runMUSCAT2_v1 스타일로 구현
- `runLIMMA_trend_v1()`: runMUSCAT2_v1 스타일로 구현
- `runLIMMA_wt_v1()`: runMUSCAT2_v1 스타일로 구현
- `runDREAM_v1()`: runMUSCAT2_v1 스타일로 구현

#### 4.2 edgeR 독립 구현
- `runEDGER_LRT_v1()`: LRT 옵션
- `runEDGER_QLF_v1()`: QLF 옵션
- `runEDGER_robust_v1()`: robust 옵션

#### 4.3 DESeq2 독립 구현
- `runDESEQ2_Wald_v1()`: Wald test
- `runDESEQ2_LRT_v1()`: LRT

### 5. 테스트 방법

```r
# 데이터 로드
library(qs)
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR")
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")  # runMUSCAT2_v1, runNEBULA2_v1 로드

# 실행
result <- run_deg_consensus(
  sobj = is5s,
  contrast = "2 - 1",
  methods = c("muscat-edgeR", "muscat-DESeq2", "nebula"),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM"
)

# 결과 확인
str(result)
result$methods_run
result$methods_failed
```

## 주의사항

1. **함수 의존성**: `runMUSCAT2_v1`과 `runNEBULA2_v1`이 환경에 로드되어 있어야 함
2. **nebula 인자**: nebula는 contrast 문자열을 사용하지 않으므로 자동 제거됨
3. **Phase 2 함수들**: 아직 구현되지 않았으므로 호출 시 에러 발생

