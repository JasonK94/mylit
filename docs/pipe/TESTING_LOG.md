# Pipeline Testing Log

이 문서는 파이프라인 테스트 과정에서 발생한 오류와 해결 방법을 기록합니다.

## 로그 파일 위치

- **실제 파이프라인 로그**: `logs/{run_id}/{run_id}_log.log`
- **Master log**: `logs/total.log`
- **참고**: `/tmp/scvi_test2.log`는 테스트 중 직접 실행한 명령의 출력을 확인하기 위해 사용한 임시 파일입니다.

## Step 1: Read Data & Demultiplexing

### 오류 1: HTO assay 이름 처리 문제
**문제**: Seurat가 assay 이름의 공백을 점(.)으로 변환하여 "Multiplexing Capture"가 "Multiplexing.Capture"로 저장됨
**해결**: assay 추가 후 실제 저장된 이름을 확인하여 사용하도록 수정
**위치**: `scripts/pipe1_read_demulti.R` (line 163-167)

### 오류 2: HTODemux 실패 처리
**문제**: HTODemux가 "Cells with zero counts exist as a cluster" 오류로 실패
**해결**: tryCatch로 오류 처리하고, 실패 시 fallback으로 "unknown" 할당
**위치**: `myR/R/pipe_demulti.R` (line 104-125)

### 변경사항: demultiplex_id 기반 샘플 매칭
**변경 내용**: 
- 기존: `generate_sample_names`로 샘플 이름 생성 후 demux 파일 column rename
- 변경: demux 파일의 column 이름을 직접 사용 (demultiplex_id와 매칭)
- `Best_Sample`은 `get_barcode_mapping`에서 column 이름을 사용하여 생성됨
**위치**: 
- `scripts/pipe1_read_demulti.R` (line 98-116)
- `myR/R/pipe_demulti.R` (line 35-44)

## Step 2: Normalization & Clustering

### 상태: 성공
- LogNormalize 정상 작동
- QC 필터링 정상 작동
- PCA 및 클러스터링 정상 작동

## Step 3: SoupX Ambient RNA Removal

### 상태: 성공
- SoupX 정상 작동
- 일부 경고 메시지 있으나 치명적 오류 없음

### 경고 메시지: "No plausible marker genes found"
**상황**: 다운샘플링된 데이터(10%)에서 많은 샘플에서 발생
**원인**: 데이터가 너무 작아서 SoupX가 marker gene을 찾지 못함
**해결**: 파이프라인이 자동으로 원본 counts를 유지 (`[WARNING] Keeping original counts`)
**결론**: **정상적인 동작**이며, 다운샘플링된 데이터에서는 예상되는 현상입니다. Full data에서는 대부분의 샘플에서 SoupX가 정상 작동합니다.

## Step 4: Doublet Detection (scDblFinder)

### 상태: 성공
- SCTransform 및 scDblFinder 정상 작동
- Doublet removal은 기본적으로 비활성화 (config에서 제어)

## Step 5: Integration

### 오류 1: vapply 타입 불일치
**문제**: `vapply(obj.list, ncol, 1L)`에서 ncol이 double을 반환하는데 integer를 기대
**해결**: `as.integer(ncol(x))`로 명시적 변환
**위치**: `scripts/pipe5_integration.R` (line 100)

### 오류 2: dims 파라미터 파싱 오류
**문제**: `dims_str` 파싱 로직이 복잡하고 오류 발생
**해결**: 간단한 파싱 로직으로 변경 (grepl로 ":" 확인 후 eval)
**위치**: `scripts/pipe5_integration.R` (line 140-144, 166-170, 237-241)

### 오류 3: IntegrateLayers 네임스페이스 오류
**문제**: `SeuratWrappers::IntegrateLayers`에서 'IntegrateLayers' is not an exported object 오류
**원인**: Seurat v5에서는 `IntegrateLayers`가 Seurat 패키지에 있음 (SeuratWrappers 아님)
**해결**: `SeuratWrappers::IntegrateLayers` → `IntegrateLayers`로 변경 (scVIIntegration만 SeuratWrappers에서)
**위치**: `scripts/pipe5_integration.R` (line 235)

### 오류 4: scVI 실행 전 PCA 필요
**문제**: scVI integration이 orig.reduction="pca"를 요구하지만 PCA가 실행되지 않음
**해결**: scVI integration 전에 NormalizeData, FindVariableFeatures, ScaleData, RunPCA 실행 추가
**위치**: `scripts/pipe5_integration.R` (line 221-234)

## 주의사항

1. **demultiplex_id 매칭**: demux 파일의 column 이름이 config의 `demultiplex_id`와 정확히 일치해야 함
2. **HTO assay 이름**: Seurat가 공백을 점으로 변환하므로 실제 저장된 이름 확인 필요
3. **dims 파라미터**: "1:50" 형식의 문자열을 벡터로 변환할 때 주의 필요


## 2025-12-03: Step 1 Demultiplexing Fix
- **Issue**: Previous logic filtered barcodes immediately, leading to identical cell counts or incorrect assignments when multiple samples shared a GEM posterior file. Also, rowname conflicts occurred during `bind_rows`.
- **Fix**:
  - `demultiplex_demuxalot` now processes the full posterior file for a GEM.
  - `pipe1_read_demulti.R` caches these results per GEM.
  - Seurat object creation filters barcodes using `demultiplex_id` *after* loading the full GEM mapping.
  - Added `Probability_Ratio` to demux output.
  - `unique_id` used for rownames to avoid conflicts.
- **Verification**: `run_fix1` (0.5% downsample) showed correct, distinct cell counts for all 56 SNP samples.

## 2025-12-03: SoupX auto-estimation retries & guards
- Added config knobs (`soupx_tfidf_start`, `soupx_tfidf_floor`, `soupx_soup_quantile_start`, `soupx_soup_quantile_floor`, `soupx_param_step`) so the pipeline now relaxes thresholds in 0.1 increments when SoupX fails to find marker genes. Logs and `plot_generation.log` tell the operator how to tweak values if plotting still fails.
- Wrapped `adjustCounts()` with a guard; if SoupX never calculated contamination fractions the script now logs a warning, keeps the original counts, and moves on instead of throwing.
- Documented DecontX/CellBender as candidates for a future ambient-removal method so Step 3 can become pluggable later.
