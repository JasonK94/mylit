# Multi-Model DEG Consensus Engine - 개발 계획서

## 프로젝트 개요

여러 DEG 방법론을 동일한 데이터셋에 적용하여 유전자 수준의 통계량을 추출하고, 방법론 수준에서의 클러스터링을 통해 consensus DEG signature를 생성하는 시스템을 구축합니다.

## 목표

1. **g3==2 vs g3==1 그룹 간 DEG 분석**을 여러 방법론으로 수행
2. 각 방법론에서 **유전자 수준 통계량 추출** (effect size, p-value, FDR 등)
3. **방법론 수준 클러스터링** 수행 (PCA, NMF, UMAP, hierarchical clustering 등)
4. **Consensus DEG signature** 생성

## 데이터 사양

- **테스트 데이터**: `/data/user3/sobj/IS_scvi_251107_ds2500.qs` (2500 cells downsampled)
- **전체 데이터**: `/data/user3/sobj/IS_scvi_251107.qs` (full dataset)
- **포맷**: Seurat object
- **메타데이터**:
  - `sample`: `hos_no`
  - `target variable`: `g3` (1, 2, "NA")
  - `batch`: `GEM`
  - `biological factor`: `sex`
  - `cluster`: `anno3.scvi`

## 개발 단계

### Phase 1: 통합 함수 구현 (메인 함수)

#### 1.1 메인 통합 함수 설계 및 구현
```r
run_deg_consensus(
  sobj,
  contrast = "2 - 1",  # g3==2 vs g3==1
  methods = c("muscat-edgeR", "muscat-DESeq2", "muscat-limma-voom", 
              "muscat-limma-trend", "nebula", "limma-voom", 
              "limma-trend", "limma-wt", "dream",
              "edgeR-LRT", "edgeR-QLF", "edgeR-robust",
              "DESeq2-Wald", "DESeq2-LRT"),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM",
  ...
)
```

**기능:**
- 각 방법론을 순차/병렬로 실행
- 결과를 리스트로 수집
- Phase 3의 통합 함수로 전달

### Phase 2: limma 계열 방법론 표준화 + 각 방법론 옵션 다양화

#### 2.1 muscat 기반 방법론 (기존 함수 활용)
- [x] **muscat-edgeR**: `runMUSCAT2_v1(method="edgeR")` 활용
- [x] **muscat-DESeq2**: `runMUSCAT2_v1(method="DESeq2")` 활용
- [x] **muscat-limma-voom**: `runMUSCAT2_v1(method="limma-voom")` 활용
- [x] **muscat-limma-trend**: `runMUSCAT2_v1(method="limma-trend")` 활용
- **참고**: muscat 내부에서 이미 edgeR, DESeq2, limma를 지원하므로 별도 구현 불필요

#### 2.2 limma 계열 독립 구현 (runMUSCAT2_v1 스타일로 표준화)
- [ ] **limma-voom** (`runLIMMA_voom_v1`)
  - Pseudobulk 생성 (muscat 스타일)
  - voom 변환
  - limma 분석
  - runMUSCAT2_v1과 동일한 인터페이스 및 출력 형식
- [ ] **limma-trend** (`runLIMMA_trend_v1`)
  - Pseudobulk 생성
  - logCPM 변환
  - limma-trend 분석
- [ ] **limma-wt** (`runLIMMA_wt_v1`)
  - Pseudobulk 생성
  - limma with weights 분석
- [ ] **dream** (`runDREAM_v1`)
  - Pseudobulk 생성
  - Random effects 포함
  - dream 분석

#### 2.3 edgeR 독립 구현 (옵션 다양화)
- [ ] **edgeR-LRT** (`runEDGER_LRT_v1`)
  - Pseudobulk 생성
  - LRT (Likelihood Ratio Test)
  - runMUSCAT2_v1과 동일한 인터페이스
- [ ] **edgeR-QLF** (`runEDGER_QLF_v1`)
  - Pseudobulk 생성
  - QLF (Quasi-Likelihood F-test)
- [ ] **edgeR-robust** (`runEDGER_robust_v1`)
  - Pseudobulk 생성
  - Robust estimation

#### 2.4 DESeq2 독립 구현 (옵션 다양화)
- [ ] **DESeq2-Wald** (`runDESEQ2_Wald_v1`)
  - Pseudobulk 생성
  - DESeq2 Wald test
  - runMUSCAT2_v1과 동일한 인터페이스
- [ ] **DESeq2-LRT** (`runDESEQ2_LRT_v1`)
  - Pseudobulk 생성
  - DESeq2 LRT

#### 2.5 Single-cell 기반 방법론
- [x] **nebula**: `runNEBULA2_v1` 활용 (기존 함수)
- [ ] **MAST**: 제외 (버그 많음)

### Phase 3: 결과 형식 통합 함수

#### 3.1 통계량 추출 및 표준화 함수
```r
standardize_deg_results(
  deg_result,      # 각 방법론의 결과
  method_name,     # 방법론 이름 (예: "muscat-edgeR", "limma-voom")
  source_function, # 원본 함수 이름 (예: "runMUSCAT2_v1", "runLIMMA_voom_v1")
  ...
)
```

**표준화 항목:**
- `gene`: 유전자 ID (표준화)
- `logFC` 또는 `beta`: Effect size
- `pvalue`: P-value
- `pvalue_adj` 또는 `FDR`: 조정된 p-value (FDR)
- `statistic`: Test statistic
- `se`: Standard error (가능한 경우)
- `method`: 방법론 이름
- `cluster_id`: 클러스터 ID (해당하는 경우)

**처리 사항:**
- 컬럼명 통일 (p_val_adj vs adj_p_val vs FDR 등)
- NA 처리
- 유전자 ID 매칭

#### 3.2 행렬 구성 함수
```r
build_deg_matrices(
  standardized_results_list,  # 표준화된 결과 리스트
  genes = NULL,              # 공통 유전자 집합 (NULL이면 교집합)
  ...
)
```

**출력:**
- `B[g, m]`: Effect size matrix (genes × methods)
- `P[g, m]`: P-value matrix
- `L[g, m]`: -log10(p) matrix
- `S[g, m]`: Significance indicator matrix (FDR < 0.05)
- `SE[g, m]`: Standard error matrix (가능한 경우)
- `T[g, m]`: Test statistic matrix

### Phase 3: 통계량 추출 및 행렬 구성 (Week 3-4)

#### 3.1 통계량 추출 함수
```r
extract_deg_statistics(
  deg_result,  # 각 방법론의 결과
  method_name, # 방법론 이름
  gene_col = "gene",  # 유전자 컬럼명
  ...
)
```

**추출 항목:**
- Effect size (beta/logFC)
- Standard error
- Test statistic
- P-value
- FDR (Benjamini-Hochberg)
- Optional: weights, dispersion

#### 3.2 행렬 구성 함수
```r
build_deg_matrices(
  deg_results_list,  # 모든 방법론 결과 리스트
  genes = NULL,  # 공통 유전자 집합 (NULL이면 교집합)
  ...
)
```

**출력:**
- `B[g, m]`: Effect size matrix (genes × methods)
- `P[g, m]`: P-value matrix
- `L[g, m]`: -log10(p) matrix
- `S[g, m]`: Significance indicator matrix (FDR < 0.05)
- `SE[g, m]`: Standard error matrix
- `T[g, m]`: Test statistic matrix

### Phase 4: Consensus 분석 (Week 4-5)

#### 4.1 Agreement Score 계산
```r
compute_agreement_scores(
  significance_matrix,  # S[g, m]
  ...
)
```

- `A_g = mean(S[g, m] across models)`: 각 유전자에 대한 방법론 간 일치도

#### 4.2 Unsupervised Integration

##### 4.2.1 PCA
```r
perform_deg_pca(
  beta_matrix,      # B[g, m]
  logp_matrix,      # L[g, m]
  significance_matrix,  # S[g, m]
  ...
)
```

##### 4.2.2 NMF
```r
perform_deg_nmf(
  beta_matrix,
  logp_matrix,
  significance_matrix,
  rank = NULL,  # 자동 선택
  ...
)
```

##### 4.2.3 UMAP / Diffusion Map
```r
perform_deg_umap(
  beta_matrix,
  logp_matrix,
  significance_matrix,
  ...
)
```

##### 4.2.4 Clustering
```r
cluster_deg_methods(
  deg_matrices,  # 통합된 행렬들
  method = c("kmeans", "hierarchical", "hdbscan"),
  ...
)
```

### Phase 5: Consensus DEG Signature 생성 (Week 5-6)

#### 5.1 Consensus 점수 계산
```r
compute_consensus_scores(
  deg_matrices,
  agreement_scores,
  clustering_results,
  weights = NULL,  # 방법론별 가중치
  ...
)
```

#### 5.2 최종 DEG 리스트 생성
```r
generate_consensus_deg_list(
  consensus_scores,
  fdr_threshold = 0.05,
  agreement_threshold = 0.5,  # 최소 50% 방법론에서 유의
  ...
)
```

### Phase 6: 시각화 및 결과 출력 (Week 6)

#### 6.1 시각화 함수
- [ ] PCA plot (methods in PC space)
- [ ] NMF heatmap
- [ ] UMAP plot
- [ ] Hierarchical clustering dendrogram
- [ ] Agreement score distribution
- [ ] Consensus DEG volcano plot
- [ ] Method comparison heatmap

#### 6.2 결과 출력
- [ ] Per-method DEG tables
- [ ] Gene × method matrices
- [ ] Consensus DEG list
- [ ] Stability metrics
- [ ] Clustering results

## 파일 구조

```
_wt/deg-consensus/
├── myR/
│   └── R/
│       ├── deg_consensus.R          # 메인 함수
│       ├── deg_methods/
│       │   ├── deg_muscat.R         # muscat wrapper
│       │   ├── deg_limma.R          # limma/voom/dream
│       │   ├── deg_edger.R          # edgeR variants
│       │   ├── deg_deseq2.R         # DESeq2 variants
│       │   ├── deg_nebula.R         # nebula wrapper
│       │   └── deg_mast.R           # MAST (optional)
│       ├── deg_extraction.R         # 통계량 추출
│       ├── deg_matrices.R           # 행렬 구성
│       ├── deg_consensus_analysis.R # Consensus 분석
│       ├── deg_clustering.R         # 방법론 클러스터링
│       └── deg_visualization.R      # 시각화
├── docs/
│   └── deg_consensus/
│       ├── DEVELOPMENT_PLAN.md     # 이 문서
│       ├── context_deg.md           # 스펙 문서
│       ├── TEST_INSTRUCTIONS_example.md
│       └── API_REFERENCE.md         # API 문서
└── scripts/
    └── test_deg_consensus.R         # 테스트 스크립트
```

## 구현 우선순위

### Priority 1 (필수)
1. `runMUSCAT2_v1` 통합
2. `runNEBULA2_v1` 통합
3. limma/voom 구현
4. edgeR-LRT 구현
5. 통계량 추출 함수
6. 행렬 구성 함수
7. Agreement score 계산
8. PCA 기반 클러스터링

### Priority 2 (중요)
1. DESeq2-Wald 구현
2. dream 구현
3. NMF 기반 클러스터링
4. UMAP 기반 클러스터링
5. Consensus DEG 리스트 생성

### Priority 3 (선택)
1. edgeR-QLF, edgeR-robust
2. DESeq2-LRT
3. limma-trend, limma-wt
4. MAST
5. HDBSCAN 클러스터링

## 테스트 계획

### 단위 테스트
- 각 방법론별 독립 테스트
- 통계량 추출 정확도 검증
- 행렬 구성 정확도 검증

### 통합 테스트
- 전체 파이프라인 테스트 (is5s 데이터)
- 결과 일관성 검증
- 메모리/성능 테스트

### 검증 테스트
- 기존 방법론 결과와 비교
- Consensus 결과의 생물학적 타당성 검증

## 성능 고려사항

1. **메모리 관리**
   - 큰 행렬의 경우 희소 행렬 활용
   - 중간 결과 디스크 저장 옵션

2. **병렬 처리**
   - 방법론별 병렬 실행
   - 클러스터링 알고리즘 병렬화

3. **캐싱**
   - Pseudobulk 결과 캐싱
   - 중간 계산 결과 저장

## 의존성

### 필수 패키지
- Seurat
- muscat
- nebula
- limma
- edgeR
- DESeq2
- SingleCellExperiment
- SummarizedExperiment

### 선택 패키지
- MAST
- NMF
- uwot (UMAP)
- dbscan (HDBSCAN)
- cluster

## 다음 단계

1. Phase 1 시작: 기존 함수 통합 및 인터페이스 설계
2. 각 Phase별 마일스톤 설정
3. 코드 리뷰 및 테스트 계획 수립
4. 문서화 및 예제 작성

## 참고 자료

- `docs/deg_consensus/context_deg.md`: 상세 스펙
- `docs/deg_consensus/TEST_INSTRUCTIONS_example.md`: 테스트 가이드
- `myR/R/test_analysis.R`: 기존 DEG 함수들
- muscat documentation
- nebula documentation

