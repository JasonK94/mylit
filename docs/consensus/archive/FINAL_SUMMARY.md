# Multi-Model DEG Consensus Engine - 최종 완료 요약

## 완료된 작업

### Phase 1: 통합 함수 구현 ✅
- `run_deg_consensus()`: 여러 방법론을 통합 실행하는 메인 함수
- 파일: `myR/R/deg_consensus/run_deg_consensus.R`

### Phase 2: 각 방법론 구현 ✅

#### 2.1 limma 계열
- `runLIMMA_voom_v1()`: limma-voom 분석
- `runLIMMA_trend_v1()`: limma-trend 분석
- 파일: `myR/R/deg_consensus/deg_methods_limma.R`

#### 2.2 edgeR 독립 구현
- `runEDGER_LRT_v1()`: edgeR Likelihood Ratio Test
- `runEDGER_QLF_v1()`: edgeR Quasi-Likelihood F-test
- 파일: `myR/R/deg_consensus/deg_methods_edger.R`

#### 2.3 DESeq2 독립 구현
- `runDESEQ2_Wald_v1()`: DESeq2 Wald test
- `runDESEQ2_LRT_v1()`: DESeq2 Likelihood Ratio Test
- 파일: `myR/R/deg_consensus/deg_methods_deseq2.R`

#### 2.4 기존 함수 활용
- `runMUSCAT2_v1()`: muscat 기반 (edgeR, DESeq2, limma-voom, limma-trend 지원)
- `runNEBULA2_v1()`: nebula 기반

### Phase 3: 결과 형식 통합 ✅
- `standardize_deg_results()`: 다양한 형식의 결과를 표준화
- `build_deg_matrices()`: 유전자 × 방법론 행렬 구성
- 파일: `myR/R/deg_consensus/deg_standardize.R`

### Phase 4-5: Consensus 분석 ✅
- `compute_agreement_scores()`: 방법론 간 일치도 계산
- `perform_deg_pca()`: 방법론 수준 PCA
- `cluster_deg_methods()`: 방법론 클러스터링
- `compute_consensus_scores()`: Consensus 점수 계산
- `generate_consensus_deg_list()`: 최종 Consensus DEG 리스트 생성
- 파일: `myR/R/deg_consensus/deg_consensus_analysis.R`

## 지원하는 방법론

### Pseudobulk 기반
1. **muscat-edgeR**: `runMUSCAT2_v1(method="edgeR")`
2. **muscat-DESeq2**: `runMUSCAT2_v1(method="DESeq2")`
3. **muscat-limma-voom**: `runMUSCAT2_v1(method="limma-voom")`
4. **muscat-limma-trend**: `runMUSCAT2_v1(method="limma-trend")`
5. **limma-voom**: `runLIMMA_voom_v1()`
6. **limma-trend**: `runLIMMA_trend_v1()`
7. **edgeR-LRT**: `runEDGER_LRT_v1()`
8. **edgeR-QLF**: `runEDGER_QLF_v1()`
9. **DESeq2-Wald**: `runDESEQ2_Wald_v1()`
10. **DESeq2-LRT**: `runDESEQ2_LRT_v1()`

### Single-cell 기반
11. **nebula**: `runNEBULA2_v1()`

## 사용 방법

### 1. 전체 파이프라인 실행
```r
# 함수 로드
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 데이터 로드
is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 전체 파이프라인 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_full_pipeline.R")
```

### 2. 개별 단계 실행
```r
# Phase 1-2: 방법론 실행
result <- run_deg_consensus(
  sobj = is5s,
  contrast = "2 - 1",
  methods = c("muscat-edgeR", "limma-voom", "edgeR-LRT"),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM"
)

# Phase 3: 표준화 및 행렬 구성
standardized <- lapply(result$methods_run, function(m) {
  standardize_deg_results(result$results[[m]], m)
})
names(standardized) <- result$methods_run

matrices <- build_deg_matrices(standardized)

# Phase 4-5: Consensus 분석
agreement <- compute_agreement_scores(matrices$significance)
consensus_scores <- compute_consensus_scores(matrices, agreement)
consensus_deg <- generate_consensus_deg_list(consensus_scores)
```

## 출력 파일

테스트 실행 후 생성되는 파일:
- `/data/user3/sobj/test_deg_consensus_full_pipeline.qs`: Phase 1-2 결과
- `/data/user3/sobj/test_deg_consensus_phase3.qs`: Phase 3 결과
- `/data/user3/sobj/test_deg_consensus_final_result.qs`: 최종 결과

## 주요 특징

1. **표준화된 인터페이스**: 모든 방법론이 동일한 인터페이스 사용
2. **자동 에러 처리**: 실패한 방법론은 건너뛰고 계속 진행
3. **유연한 결과 통합**: 다양한 결과 형식을 자동으로 표준화
4. **Consensus 분석**: 방법론 간 일치도를 기반으로 robust DEG signature 생성

## 다음 단계 (선택)

1. **시각화 함수 추가** (Phase 6)
   - PCA plot
   - Clustering dendrogram
   - Agreement score distribution
   - Consensus DEG volcano plot

2. **추가 방법론**
   - limma-wt
   - dream
   - edgeR-robust

3. **성능 최적화**
   - 병렬 처리
   - 메모리 최적화

