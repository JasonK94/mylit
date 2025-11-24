# Differential Expression 함수들의 데이터 프로세싱 절차 및 사용 함수

## 개요

본 문서는 `runMAST_v1`, `runNEBULA_v1`, `runMUSCAT` 함수들의 데이터 프로세싱 절차와 사용 함수를 상세히 설명합니다.

## 1. runMAST_v1

### 1.1 데이터 프로세싱 절차

#### 단계 1: Seurat → SingleCellAssay (SCA) 변환
```r
sce <- Seurat::as.SingleCellExperiment(sobj)
sca <- MAST::SceToSingleCellAssay(sce)
# Fallback: MAST::FromMatrix() 사용
```

**사용 함수:**
- `Seurat::as.SingleCellExperiment()`: Seurat 객체를 SingleCellExperiment로 변환
- `MAST::SceToSingleCellAssay()`: SCE를 SingleCellAssay로 변환
- `MAST::FromMatrix()`: Fallback - Matrix에서 직접 SCA 생성

#### 단계 2: 유전자 필터링
```r
keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
sca_filtered <- sca[keep_genes, ]
```

**사용 함수:**
- `MAST::freq()`: 각 유전자의 발현 비율 (0~1) 계산

#### 단계 3: 정규화
```r
SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
```

**사용 함수:**
- `MAST::cpm()`: Counts Per Million 계산 (log 변환 옵션)

#### 단계 4: MAST::zlm 실행
```r
zfit <- MAST::zlm(formula_obj, sca = sca_filtered, method = "glmer", parallel = TRUE, nCores = n_cores)
```

**사용 함수:**
- `MAST::zlm()`: Zero-inflated hurdle model 피팅

#### 단계 5: LRT 결과 요약
```r
summary_res <- summary(zfit, doLRT = lrt_variable)
```

**사용 함수:**
- `summary()`: MAST zlm 객체의 요약 및 LRT 검정

### 1.2 데이터 흐름

```
Seurat object
  ↓ [Seurat::as.SingleCellExperiment]
SingleCellExperiment (SCE)
  ↓ [MAST::SceToSingleCellAssay]
SingleCellAssay (SCA)
  ↓ [MAST::freq() + 필터링]
Filtered SCA (발현 세포 수 기준)
  ↓ [MAST::cpm()]
Log2(CPM+1) 정규화된 SCA
  ↓ [MAST::zlm()]
MAST zlm 객체
  ↓ [summary() with LRT]
결과 데이터프레임
```

## 2. runNEBULA_v1

### 2.1 데이터 프로세싱 절차

#### 단계 0: 데이터 추출
```r
meta <- sobj@meta.data
counts <- GetAssayData(sobj, layer = layer)
```

**사용 함수:**
- `GetAssayData()`: Seurat 객체에서 count matrix 추출

#### 단계 1: 유전자 필터링
```r
keep_genes <- rowSums(counts > 0) >= min_count
counts_filtered <- counts[keep_genes, ]
```

**사용 함수:**
- `rowSums()`: 각 유전자의 발현 세포 수 계산

#### 단계 2: NA 값 확인 및 제거
```r
keep_cells_idx <- complete.cases(meta[, vars_to_check])
meta_clean <- meta[keep_cells_idx, ]
counts_clean <- counts_filtered[, keep_cells_idx]
```

**사용 함수:**
- `complete.cases()`: NA가 없는 행 확인

#### 단계 3: 디자인 행렬 및 벡터 생성
```r
meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)
design_matrix <- model.matrix(as.formula(formula_str), data = meta_clean)
```

**사용 함수:**
- `as.factor()`: 범주형 변수 변환
- `model.matrix()`: 설계 행렬 생성

#### 단계 4: group_cell()로 데이터 정렬
```r
data_grouped <- nebula::group_cell(
  count = counts_clean,
  id = id_vector,
  pred = design_matrix,
  offset = offset_vector
)
```

**사용 함수:**
- `nebula::group_cell()`: NEBULA를 위한 데이터 그룹화 및 정렬

#### 단계 5: NEBULA 실행
```r
re_nebula <- nebula::nebula(
  count = data_grouped$count,
  id = data_grouped$id,
  pred = data_grouped$pred,
  offset = data_grouped$offset,
  model = "NBLMM",
  method = "HL"
)
```

**사용 함수:**
- `nebula::nebula()`: Negative Binomial mixed-effects model 피팅

### 2.2 데이터 흐름

```
Seurat object
  ↓ [GetAssayData()]
Count matrix (dgCMatrix)
  ↓ [rowSums() + 필터링]
Filtered count matrix (발현 세포 수 기준)
  ↓ [complete.cases()]
NA 제거된 데이터
  ↓ [as.factor() + model.matrix()]
설계 행렬 + factor 변환된 메타데이터
  ↓ [nebula::group_cell()]
Grouped data (patient별로 정렬)
  ↓ [nebula::nebula()]
NEBULA 결과 객체
```

## 3. runMUSCAT

### 3.1 데이터 프로세싱 절차

#### 단계 1: Seurat → SCE, prepSCE
```r
sce <- Seurat::as.SingleCellExperiment(sobj)
sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)
```

**사용 함수:**
- `Seurat::as.SingleCellExperiment()`: Seurat 객체를 SCE로 변환
- `muscat::prepSCE()`: SCE를 MUSCAT 분석에 맞게 준비

#### 단계 2: Pseudobulking
```r
pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))
```

**사용 함수:**
- `muscat::aggregateData()`: 클러스터와 샘플별로 pseudobulk 생성

#### 단계 2-1: Pseudobulk 메타데이터 보강
```r
pb_meta <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
```

**사용 함수:**
- `dplyr::left_join()`: 메타데이터 병합

#### 단계 3: Contrast 그룹 필터링
```r
keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
pb_sub <- pb[, keep_idx]
```

**사용 함수:**
- `SummarizedExperiment::colData()`: SCE의 메타데이터 접근

#### 단계 4: Design matrix 생성 및 DE 분석
```r
design <- stats::model.matrix(~ 0 + group + batch, data = colData(pb_sub))
```

**사용 함수:**
- `stats::model.matrix()`: 설계 행렬 생성
- `limma::voom()`: edgeR의 voom 변환 (limma-voom 사용 시)
- `edgeR::glmQLFit()`, `edgeR::glmQLFTest()`: edgeR 분석
- `DESeq2::DESeq()`: DESeq2 분석
- `limma::lmFit()`, `limma::eBayes()`, `limma::topTable()`: limma 분석

#### 단계 5: muscat::pbDS() 또는 직접 DE 분석
```r
resDS <- muscat::pbDS(pb_sub, method = method, ...)
```

**사용 함수:**
- `muscat::pbDS()`: Pseudobulk differential expression 분석

### 3.2 데이터 흐름

```
Seurat object
  ↓ [Seurat::as.SingleCellExperiment]
SingleCellExperiment (SCE)
  ↓ [muscat::prepSCE()]
Prepared SCE (cluster_id, sample_id, group_id 포함)
  ↓ [muscat::aggregateData()]
Pseudobulk SCE (클러스터별로 assays 분리)
  ↓ [dplyr::left_join()]
메타데이터 보강된 Pseudobulk
  ↓ [그룹 필터링]
Contrast 그룹만 포함된 Pseudobulk
  ↓ [stats::model.matrix() + DE 분석]
DE 분석 결과 (edgeR/DESeq2/limma)
```

## 4. 공통 사용 부분

### 4.1 공통 데이터 변환

#### 1. Seurat → SingleCellExperiment 변환
**사용 함수:**
- `Seurat::as.SingleCellExperiment()`

**사용 함수:**
- `runMAST_v1`: ✓ 사용 (SCE를 SCA로 추가 변환)
- `runMUSCAT`: ✓ 사용 (직접 사용)
- `runNEBULA_v1`: ✗ 사용 안 함 (직접 GetAssayData 사용)

#### 2. 유전자 필터링
**공통 패턴:**
- 발현 세포 수 기준 필터링

**차이점:**
- `runMAST_v1`: `MAST::freq() * ncol(sca) >= min_cells_expr`
- `runNEBULA_v1`: `rowSums(counts > 0) >= min_count`
- `runMUSCAT`: `muscat::pbDS()` 내부에서 처리

#### 3. 메타데이터 처리
**공통 패턴:**
- NA 값 확인 및 제거
- Factor 변환

**차이점:**
- `runMAST_v1`: SCA 내부에 메타데이터 포함
- `runNEBULA_v1`: `complete.cases()` 사용, `as.factor()` 적용
- `runMUSCAT`: `muscat::prepSCE()` 내부 처리, `droplevels()` 적용

#### 4. 설계 행렬 생성
**공통 패턴:**
- `model.matrix()` 사용

**차이점:**
- `runMAST_v1`: Formula 기반 (사용자 제공)
- `runNEBULA_v1`: Fixed effects 기반 (자동 생성)
- `runMUSCAT`: `stats::model.matrix()` 사용 (batch 포함 가능)

### 4.2 공통 R 패키지

#### Seurat 패키지
- **사용:** 모든 함수
- **용도:** 입력 데이터 처리

#### SummarizedExperiment 패키지
- **사용:** `runMAST_v1`, `runMUSCAT`
- **용도:** SCE 객체 처리

#### SingleCellExperiment 패키지
- **사용:** `runMAST_v1`, `runMUSCAT`
- **용도:** SCE 객체 생성 및 조작

#### S4Vectors 패키지
- **사용:** `runMUSCAT`
- **용도:** DataFrame 생성

### 4.3 공통 처리 패턴

#### 1. NA 값 처리
```r
# 패턴 1: complete.cases() 사용 (runNEBULA_v1)
keep_cells_idx <- complete.cases(meta[, vars_to_check])

# 패턴 2: is.na() 체크 (runMUSCAT)
na_mask <- is.na(meta[[group_id]]) | is.na(meta[[cluster_id]])
```

#### 2. Factor 변환
```r
# 패턴 1: lapply + as.factor (runNEBULA_v1)
meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)

# 패턴 2: 직접 변환 (runMUSCAT)
sce$cluster_id <- droplevels(factor(colData(sce)$cluster_id))
```

#### 3. 데이터 필터링
```r
# 패턴: boolean indexing 사용
filtered_data <- data[keep_idx, ]
```

## 5. 주요 차이점

### 5.1 데이터 레벨

| 함수 | 데이터 레벨 | 특징 |
|------|------------|------|
| `runMAST_v1` | Single-cell | 각 세포를 개별적으로 분석 |
| `runNEBULA_v1` | Single-cell | 각 세포를 개별적으로 분석 (random effects 포함) |
| `runMUSCAT` | Pseudobulk | 클러스터와 샘플별로 집계 후 분석 |

### 5.2 정규화 방법

| 함수 | 정규화 방법 | 함수 |
|------|------------|------|
| `runMAST_v1` | Log2(CPM+1) | `MAST::cpm(log = TRUE)` |
| `runNEBULA_v1` | None (offset 사용) | Offset 변수로 정규화 |
| `runMUSCAT` | DE 방법에 따라 다름 | edgeR/DESeq2/limma 내부 정규화 |

### 5.3 모델 타입

| 함수 | 모델 타입 | 특징 |
|------|----------|------|
| `runMAST_v1` | Hurdle model | Zero-inflated + continuous |
| `runNEBULA_v1` | Negative Binomial mixed-effects | Random effects 포함 |
| `runMUSCAT` | edgeR/DESeq2/limma | Pseudobulk 기반 |

## 6. 공통 개선 가능 사항

### 6.1 NA 처리 통일
- 모든 함수에서 동일한 NA 처리 로직 사용
- `complete.cases()` 기반 통일

### 6.2 유전자 필터링 통일
- 동일한 필터링 기준 및 함수 사용
- `rowSums(counts > 0) >= min_count` 패턴 통일

### 6.3 메타데이터 검증 통일
- 필수 컬럼 존재 여부 확인
- Factor 변환 로직 통일

### 6.4 오류 처리 통일
- Try-catch 패턴 통일
- 오류 메시지 형식 통일

## 7. 요약

### 공통 부분
1. **Seurat → 다른 객체 변환**: `runMAST_v1`, `runMUSCAT`
2. **유전자 필터링**: 모든 함수 (방법만 다름)
3. **메타데이터 처리**: 모든 함수 (NA 제거, factor 변환)
4. **설계 행렬 생성**: `runNEBULA_v1`, `runMUSCAT`

### 고유 부분
1. **runMAST_v1**: SCA 변환, MAST::zlm, Hurdle model
2. **runNEBULA_v1**: group_cell, nebula::nebula, Random effects
3. **runMUSCAT**: Pseudobulking, muscat::pbDS, edgeR/DESeq2/limma

### 공통 사용 패키지
- `Seurat`: 모든 함수
- `SummarizedExperiment`: `runMAST_v1`, `runMUSCAT`
- `SingleCellExperiment`: `runMAST_v1`, `runMUSCAT`
- `dplyr`: `runMUSCAT`
- `stats`: `runNEBULA_v1`, `runMUSCAT`

