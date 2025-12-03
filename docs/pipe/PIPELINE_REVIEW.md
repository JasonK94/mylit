# Pipeline 스크립트 리뷰 및 배치 효과 보정 분석

## 문제: GEM1~4와 GEM5~8이 UMAP에서 분리됨

이는 **배치 효과 보정이 제대로 되지 않았다**는 신호입니다.

## 각 Step별 배치 정보 전달 확인

### Step 1: Read & Demultiplexing (`pipe1_read_demulti.R`)

✅ **GEM 정보 추가됨:**
```r
obj$GEM <- gem_name
```

✅ **문제 없음** - GEM 정보가 metadata에 저장됨

### Step 2: Normalization & Clustering (`pipe2_nmz_clustering.R`)

⚠️ **확인 필요:** GEM 정보가 그대로 전달되는지 확인

### Step 3: SoupX (`pipe3_ambient_removal.R`)

⚠️ **확인 필요:** SoupX 후에도 GEM 정보가 유지되는지

### Step 4: SCTransform (`pipe4_sctransform.R`)

⚠️ **확인 필요:** SCTransform 후에도 GEM 정보가 유지되는지

### Step 5: Doublet Detection (`pipe5_doubletfinder.R`)

⚠️ **확인 필요:** DoubletFinder 후에도 GEM 정보가 유지되는지

### Step 6: Integration (`pipe6_integration.R`)

#### RPCA Integration

**현재 코드:**
```r
batch_col <- get_param("metadata_gem_col", config_list, "GEM")
obj.list <- SplitObject(merged, split.by = batch_col)
```

✅ **문제 없음** - GEM으로 split하여 각 batch를 독립적으로 SCTransform

**하지만:**
- `PrepSCTIntegration`이 제대로 실행되는지 확인 필요
- `FindIntegrationAnchors`의 `reduction = "rpca"`가 제대로 작동하는지 확인

#### scVI Integration

**현재 코드:**
```r
batch_col <- get_param("scvi_batch_col", config_list, "GEM")
merged <- IntegrateLayers(
  object = merged,
  method = SeuratWrappers::scVIIntegration,
  batch = batch_col,
  ...
)
```

⚠️ **잠재적 문제:**
1. **JoinLayers 문제:** scVI 전에 `JoinLayers`를 호출하는데, 이게 맞는지 확인 필요
2. **Layer 문제:** `layers = "counts"`가 맞는지 확인 필요
3. **PCA 문제:** scVI 전에 PCA를 실행하는데, 이게 batch-aware한지 확인 필요

## 가능한 원인

### 1. scVI Integration 전 데이터 준비 문제

**현재:**
```r
# Join layers
merged <- JoinLayers(merged, assay = "RNA")
DefaultAssay(merged) <- "RNA"

# Normalize and PCA
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = nfeatures, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = npcs, verbose = FALSE)
```

**문제:** 
- `NormalizeData`, `ScaleData`, `RunPCA`가 **전체 데이터에 대해** 실행됨
- 이는 batch-specific normalization이 아님
- scVI는 자체적으로 normalization을 하지만, PCA는 batch-aware하지 않을 수 있음

### 2. scVI Integration 파라미터 문제

**현재:**
```r
merged <- IntegrateLayers(
  object = merged,
  method = SeuratWrappers::scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.scvi",
  batch = batch_col,
  layers = "counts",
  ...
)
```

**확인 필요:**
- `layers = "counts"`가 맞는지 (Seurat v5에서는 "counts" layer 사용)
- `orig.reduction = "pca"`가 필요한지 (scVI는 자체적으로 PCA를 할 수 있음)

### 3. RPCA Integration 문제

**현재:**
```r
# Split by GEM
obj.list <- SplitObject(merged, split.by = batch_col)

# SCTransform each batch
obj.list <- lapply(obj.list, function(x) {
  x <- SCTransform(x, ...)
  x <- RunPCA(x, ...)
  x
})

# PrepSCTIntegration
obj.list <- PrepSCTIntegration(object.list = obj.list, ...)

# FindIntegrationAnchors
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  normalization.method = "SCT",
  reduction = "rpca",
  ...
)
```

**확인 필요:**
- `PrepSCTIntegration`이 제대로 실행되는지
- 각 batch의 SCTransform이 독립적으로 실행되는지
- `FindIntegrationAnchors`의 `reduction = "rpca"`가 제대로 작동하는지

## 해결 방안

### 1. scVI Integration 개선

**제안:**
```r
# Join layers (필요한 경우만)
if (packageVersion("Seurat") >= "5.0.0") {
  # Check if layers exist
  if (length(Layers(merged, assay = "RNA")) > 1) {
    merged <- JoinLayers(merged, assay = "RNA")
  }
}

# scVI는 자체적으로 normalization을 하므로, PCA는 선택적
# 하지만 orig.reduction이 필요하면 PCA를 실행
if (!"pca" %in% names(merged@reductions)) {
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = nfeatures, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, npcs = npcs, verbose = FALSE)
}

# scVI integration
merged <- IntegrateLayers(
  object = merged,
  method = SeuratWrappers::scVIIntegration,
  new.reduction = "integrated.scvi",
  batch = batch_col,
  layers = "counts",
  # orig.reduction은 선택적 (scVI가 자체적으로 할 수 있음)
  verbose = TRUE
)
```

### 2. 배치 정보 확인

**각 step에서 GEM 정보가 유지되는지 확인:**
```r
# Step 1 후
table(sl[[1]]$GEM)

# Step 6 전
table(merged$GEM)

# Step 6 후
table(integrated$GEM)
```

### 3. Integration 결과 검증

**UMAP에서 배치 효과 확인:**
```r
DimPlot(integrated, group.by = "GEM", reduction = "umap")
```

**배치 효과가 보정되었다면:**
- GEM1~4와 GEM5~8이 UMAP에서 섞여 있어야 함
- 각 GEM이 특정 영역에 뭉치지 않아야 함

## 다음 단계

1. 각 step에서 GEM 정보가 유지되는지 확인
2. scVI integration 파라미터 조정
3. RPCA integration 검증
4. Integration 결과 시각화 및 검증

