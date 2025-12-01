# Pipeline Testing Log

이 문서는 파이프라인 테스트 과정에서 발생한 오류와 해결 방법을 기록합니다.

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

## 주의사항

1. **demultiplex_id 매칭**: demux 파일의 column 이름이 config의 `demultiplex_id`와 정확히 일치해야 함
2. **HTO assay 이름**: Seurat가 공백을 점으로 변환하므로 실제 저장된 이름 확인 필요
3. **dims 파라미터**: "1:50" 형식의 문자열을 벡터로 변환할 때 주의 필요

