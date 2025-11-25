# 테스트 함수 가이드

## 주요 테스트 함수 요약표

| 함수명 | 입력 객체 | 출력 객체 | 저장 위치 | 상태 | 비고 |
|--------|----------|----------|----------|------|------|
| **find_gene_signature_v4.2** | Seurat 객체 또는 expression matrix | List (method별 결과) | 수동 저장 필요 | ✅ 작동 | v4.2는 안정적 |
| **find_gene_signature_v5.2** | Seurat 객체 또는 expression matrix | List (method별 결과) | 수동 저장 필요 | ✅ 작동 | 모든 메서드 구현 완료 |
| **find_gene_signature_v5.3** | Seurat 객체 또는 expression matrix | List (method별 결과) | 수동 저장 필요 | ✅ 작동 | v5.2 + 동적 k |
| **FGS_v5.2** | find_gene_signature_v5.2와 동일 | find_gene_signature_v5.2와 동일 | 수동 저장 필요 | ✅ 작동 | signature.R의 별칭 |
| **FGS_v5.3** | find_gene_signature_v5.3와 동일 | find_gene_signature_v5.3와 동일 | 수동 저장 필요 | ✅ 작동 | signature.R의 별칭 |
| **runMAST_v1** | Seurat 객체 | data.frame (DE 결과) | 수동 저장 필요 | ⚠️ 수정 중 | summary_res 추출 문제 |
| **runNEBULA_v1** | Seurat 객체 | NEBULA result 객체 | 수동 저장 필요 | 🔄 테스트 필요 | 코드 수정 완료 |
| **runMUSCAT_v5** | Seurat 객체 | data.frame (DE 결과) | 수동 저장 필요 | ✅ 작동 | Pseudobulk 자동 수행 |
| **TML6** | L1 signatures + holdout data | List (best_model, l2_train 등) | 수동 저장 필요 | ✅ 작동 | Meta-learner |

---

## 상세 사용법

### 1. Find Gene Signature 함수들

#### find_gene_signature_v4.2
```r
# 입력
data <- is4  # Seurat 객체 또는 expression matrix
target_var <- "g3"
control_vars <- "hos_no"
test_n <- 100

# 실행
result <- find_gene_signature_v4.2(
  data = is4,
  target_var = "g3",
  control_vars = "hos_no",
  test_n = 100
)

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/FGS_v4.2_result.qs")
```

#### find_gene_signature_v5.2
```r
# 입력
data <- is4  # Seurat 객체 또는 expression matrix
target_var <- "g3"
control_vars <- "hos_no"
method <- c("random_forest", "limma", "gam")  # 선택 가능한 메서드

# 실행
result <- find_gene_signature_v5.2(
  data = is4,
  target_var = "g3",
  control_vars = "hos_no",
  method = c("random_forest", "limma", "gam")
)

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/FGS_v5.2_result.qs")
```

#### find_gene_signature_v5.3
```r
# 입력
data <- is4  # Seurat 객체 또는 expression matrix
target_var <- "g3"
control_vars <- "hos_no"
gam.k <- NULL  # NULL이면 동적 k 계산
gam.k_dynamic_factor <- 5
gam.min_unique <- 10

# 실행
result <- find_gene_signature_v5.3(
  data = is4,
  target_var = "g3",
  control_vars = "hos_no",
  gam.k = NULL,
  gam.k_dynamic_factor = 5,
  gam.min_unique = 10
)

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/FGS_v5.3_result.qs")
```

---

### 2. Differential Expression 분석 함수들

#### runMAST_v1
```r
# 입력
sobj <- is4_light  # Seurat 객체
formula <- "~ g3"
lrt_variable <- "g3"
min_cells_expr <- 5

# 실행
result <- runMAST_v1(
  sobj = is4_light,
  formula = "~ g3",
  min_cells_expr = 5,
  n_cores = 1,
  lrt_variable = "g3"
)

# 출력: data.frame
# - primerid: Gene identifier
# - p_value_hurdle: P-value
# - coef: Coefficient
# - ci.hi, ci.lo: Confidence intervals

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/runMAST_v1_result.qs")
```

#### runNEBULA_v1
```r
# 입력
sobj <- is4  # Seurat 객체
fixed_effects <- c("g3", "celltype")
covar_effects <- c("batch")
patient_col <- "hos_no"
offset <- "nCount_RNA"

# 실행
result <- runNEBULA_v1(
  sobj = is4,
  fixed_effects = c("g3", "celltype"),
  covar_effects = c("batch"),
  patient_col = "hos_no",
  offset = "nCount_RNA"
)

# 출력: NEBULA result 객체
# - result$summary: 요약 통계
# - result$random: Random effects

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/runNEBULA_v1_result.qs")
```

#### runMUSCAT_v5
```r
# 입력
sobj <- is4  # Seurat 객체
cluster_id <- "seurat_clusters"
sample_id <- "hos_no"
group_id <- "g3"
contrast <- "2 - 1"  # g3=2 vs g3=1

# 실행
result <- runMUSCAT_v5(
  sobj = is4,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "g3",
  contrast = "2 - 1",
  method = "edgeR"
)

# 출력: data.frame (클러스터별 DE 결과)
# - cluster_id: 클러스터 ID
# - gene: Gene identifier
# - logFC, p_val, p_adj 등

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/runMUSCAT_v5_result.qs")
```

---

### 3. Meta-Learner 함수

#### TML6 (train_meta_learner_v6)
```r
# 입력
l1_signatures <- fgs_result  # find_gene_signature 결과
holdout_data <- is5_30  # Seurat 객체 (holdout set)
target_var <- "g3"
l2_methods <- c("glm", "ranger", "xgbTree")

# 실행
result <- TML6(
  l1_signatures = fgs_result,
  holdout_data = is5_30,
  target_var = "g3",
  l2_methods = c("glm", "ranger", "xgbTree")
)

# 출력: List
# - best_model: 최고 성능 모델
# - best_model_name: 모델 이름
# - l2_train: L2 학습 데이터
# - l1_signatures: L1 signatures
# - model_comparison: 모델 비교 결과

# 저장 (수동)
qs::qsave(result, "/data/user3/sobj/TML6_result.qs")
```

---

## 테스트 데이터 위치

| 데이터명 | 경로 | 용도 | 크기 |
|---------|------|------|------|
| **원본 데이터** | `/data/user3/sobj/IS_scvi_251107.qs` | Full 분석 | ~30K cells |
| **다운샘플링** | `/data/user3/sobj/IS_scvi_251107_ds2500.qs` | 빠른 테스트 | ~2.5K cells |
| **경량 데이터** | `/data/user3/sobj/IS_scvi_251107_lightweight.qs` | 디버깅 | 500 cells |
| **is5 데이터** | `/data/user3/sobj/IS5_g3NA_removal_251110.qs` | TML6 holdout | Full |

---

## 권장 테스트 순서

### 1단계: 빠른 테스트 (경량 데이터)
```r
# 데이터 로드
is4_light <- qs::qread("/data/user3/sobj/IS_scvi_251107_lightweight.qs")

# FGS 테스트
fgs_result <- find_gene_signature_v5.2(
  data = is4_light,
  target_var = "g3",
  control_vars = "hos_no",
  method = c("limma", "wilcoxon")  # 빠른 메서드만
)

# 저장
qs::qsave(fgs_result, "/data/user3/sobj/FGS_v5.2_lightweight_test.qs")
```

### 2단계: 중간 테스트 (다운샘플링 데이터)
```r
# 데이터 로드
is4_ds <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# FGS 전체 메서드 테스트
fgs_result <- find_gene_signature_v5.3(
  data = is4_ds,
  target_var = "g3",
  control_vars = "hos_no",
  method = c("random_forest", "limma", "gam", "nmf_loadings")
)

# 저장
qs::qsave(fgs_result, "/data/user3/sobj/FGS_v5.3_ds2500_test.qs")
```

### 3단계: Full 분석 (원본 데이터)
```r
# 데이터 로드
is4 <- qs::qread("/data/user3/sobj/IS_scvi_251107.qs")

# FGS 전체 메서드 실행
fgs_result <- find_gene_signature_v5.3(
  data = is4,
  target_var = "g3",
  control_vars = "hos_no"
)

# 저장
qs::qsave(fgs_result, "/data/user3/sobj/FGS_v5.3_full.qs")
```

---

## 저장 파일 네이밍 규칙

권장 형식: `{함수명}_{버전}_{데이터}_{기타정보}.qs`

예시:
- `FGS_v5.2_ds2500_test.qs`
- `FGS_v5.3_full_seed42.qs`
- `runMAST_v1_lightweight_test.qs`
- `runNEBULA_v1_ds2500_test.qs`
- `runMUSCAT_v5_full_edgeR.qs`
- `TML6_is5_30_result.qs`

---

## 주의사항

1. **모든 함수는 수동 저장 필요**: 함수 내부에서 자동 저장하지 않음
2. **실행 시간 기록 권장**: 큰 데이터셋의 경우 시간이 오래 걸릴 수 있음
3. **메모리 관리**: Full 데이터셋 사용 시 메모리 부족 주의
4. **Seed 고정**: 재현성을 위해 `set.seed(42)` 사용 권장

---

## 문제 해결

### runMAST_v1 오류
- 현재 `summary_res` datatable 추출 문제 발생 중
- 해결 전까지는 다른 DE 방법 사용 권장 (runMUSCAT_v5, runNEBULA_v1)

### 메모리 부족
- 다운샘플링 데이터로 먼저 테스트
- 메서드별로 나눠서 실행

### 실행 시간이 너무 오래 걸림
- `gam`, `nmf_loadings` 메서드는 특히 느림
- 다운샘플링 데이터 사용 또는 해당 메서드 제외

