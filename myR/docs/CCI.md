# Cell-to-Cell Interaction (CCI) Analysis Tool

## 프로젝트 의의

### 배경 및 목적

단일세포 RNA 시퀀싱(scRNAseq) 데이터에서 세포 간 상호작용(Cell-to-Cell Interaction, CCI)을 분석하는 것은 면역학, 발달생물학, 암 연구 등 다양한 분야에서 중요한 과제입니다. 특히, 특정 조건에서 receiver cell type의 유전자 발현 변화에 기여하는 sender cell type을 식별하는 것은 질병 메커니즘 이해와 치료 표적 발굴에 핵심적입니다.

본 프로젝트는 **NicheNet**을 중심으로 하여, **사전 계산된 DEG(Differential Expression Genes) 리스트를 직접 입력받아** CCI 분석을 수행할 수 있는 통합 도구를 개발하는 것을 목표로 합니다. 이는 기존의 DEG 분석 파이프라인과의 원활한 통합을 가능하게 하며, 다양한 DEG 분석 엔진(예: MUSCAT, limma-dream, edgeR 등)의 결과를 직접 활용할 수 있게 합니다.

### 주요 특징

1. **DEG 리스트 직접 입력**: 사전 계산된 DEG 결과를 데이터프레임으로 입력받아 분석
2. **자동 Sender 식별**: Receiver를 제외한 모든 클러스터를 자동으로 sender로 식별하거나, 사용자가 직접 지정 가능
3. **모듈화된 구조**: 데이터 준비, 검증, 분석, 결과 저장이 독립적인 모듈로 구성
4. **논문 수준의 시각화**: Circos plot을 포함한 고품질 플롯 생성
5. **완전한 문서화**: 사용 가이드, 테스트 스크립트, 개발 로그 포함

### 기존 도구와의 차별점

기존 `run_nichenet_analysis()` 함수는 Seurat 객체 내부에서 `FindMarkers()`를 호출하여 DEG를 계산합니다. 본 도구는:
- 외부에서 계산된 DEG 리스트를 직접 사용 가능
- 다양한 DEG 분석 엔진의 결과를 통합
- DEG 분석과 CCI 분석의 분리로 유연성 향상

## 작업 내용

### 1. 워크트리 및 브랜치 구조

```
_wt/cci/
  docs/cci/              # 문서
    cci.md              # 메인 문서
    TEST_INSTRUCTIONS.md # 테스트 가이드
    devlog.md           # 개발 로그
    SUCCESS_SUMMARY.md  # 성공 요약
    PLOT_OUTPUT_SUMMARY.md # 플롯 출력 요약
  scripts/cci/          # 테스트 스크립트
    test_cci.R          # 전체 테스트
    test_cci_actual.R   # 실제 데이터 테스트
    test_cci_interactive.R # 대화형 테스트
    test_cci_with_plots.R # 플롯 생성 테스트
  myR/R/cci/            # 함수 모듈
    run_cci_analysis.R  # 메인 함수
    prepare_cci_data.R  # 데이터 준비
    utils_cci.R         # 유틸리티
    save_cci_results.R # 결과 저장
    nichenet_wrapper.R # NicheNet 래퍼
```

### 2. 핵심 함수 모듈

#### `prepare_cci_data.R`
- `validate_cci_inputs()`: 입력 데이터 검증
- `extract_receiver_degs()`: Receiver 클러스터의 DEG 추출
- `identify_sender_clusters()`: Sender 클러스터 자동 식별 또는 사용자 지정
- `filter_expressed_genes()`: 발현 유전자 필터링 (NicheNet 호환)
- `prepare_cci_summary()`: 분석 준비 데이터 요약

#### `run_cci_analysis.R`
메인 CCI 분석 함수로 다음 단계를 수행:
1. 입력 검증
2. Receiver DEG 추출
3. Sender 클러스터 식별
4. 발현 유전자 필터링
5. 데이터 요약 생성
6. NicheNet 분석 실행 (`run_nichenet_analysis()` 호출)
7. 결과 컴파일 및 저장

#### `save_cci_results.R`
- `save_cci_intermediate()`: 중간 결과 저장
- `save_cci_final()`: 최종 결과 저장
- `load_cci_results()`: 저장된 결과 로드

#### `utils_cci.R`
- `format_deg_summary()`: DEG 요약 테이블 생성
- `identify_top_senders()`: Top sender 식별
- `create_sender_receiver_map()`: Sender-receiver 매핑 생성

### 3. NicheNet 통합

기존 `run_nichenet_analysis()` 함수(`myR/R/CCI.R`)를 활용하여:
- Ligand activity 예측
- Ligand-target 네트워크 분석
- Ligand-receptor 상호작용 시각화
- Circos plot 생성

#### 최신 최적화
- Receiver 별 DEG를 여러 번 계산하지 않도록 `receiver_de_table` 인자를 통해 사전 계산된 DEG를 재사용할 수 있습니다.
- DEG 테이블 컬럼명이 `avg_log2FC`, `p_val_adj`와 달라도 `receiver_gene_col`, `receiver_logfc_col`, `receiver_pval_col` 파라미터로 매핑할 수 있습니다.
- NicheNet 실행 전 Sender/Receiver 개수, DEG 수, 예상 소요 시간을 로그로 출력하여 진행 상황을 명확히 제공합니다.
- 실행 중간 결과는 자동으로 `.qs` 체크포인트(예: `nichenet_results.qs`)로 저장되어, 중단 시 해당 파일로부터 재현 가능합니다.
- 대용량 객체는 체크포인트 저장 후 `rm()` + `gc()`를 통해 정리하여 메모리 사용량을 최소화합니다.

### 4. 시각화 기능

다음 플롯들이 자동 생성됩니다:
- **Ligand-Target Heatmap**: 우선순위가 높은 ligand와 예측된 target 유전자 간의 상호작용
- **Ligand-Receptor Heatmap**: Ligand와 receptor 간의 상호작용
- **Ligand Activity Histogram**: Ligand activity 분포
- **Ligand AUPR Heatmap**: Ligand activity (AUPR) 히트맵
- **Circos Plot**: Ligand-receptor 상호작용을 원형 다이어그램으로 시각화 (논문 수준)

## 개발 로그

### 2025-11-14: 프로젝트 시작 및 초기 구현

#### 목표
- CCI 분석 도구 개발
- DEG 리스트 직접 입력 지원
- 모듈화된 구조 구축
- 완전한 문서화

#### 완료된 작업

**1. 워크트리 및 브랜치 생성**
- `cci` 브랜치 생성
- `_wt/cci` 워크트리 생성
- 표준 디렉터리 구조 생성 (`docs/cci/`, `scripts/cci/`, `myR/R/cci/`)

**2. 문서 작성**
- `cci.md`: 메인 기능 문서
- `TEST_INSTRUCTIONS.md`: 테스트 가이드
- `devlog.md`: 개발 로그
- `SUCCESS_SUMMARY.md`: 성공 요약
- `PLOT_OUTPUT_SUMMARY.md`: 플롯 출력 요약

**3. 함수 구현**
- 5개 모듈 파일 구현 (총 ~1500줄)
- 입력 검증, DEG 추출, Sender 식별, 발현 유전자 필터링, 결과 저장 기능

**4. 테스트 스크립트 작성**
- 5개 테스트 스크립트 작성
- 실제 데이터로 테스트 성공

### 2025-11-14: 테스트 및 디버깅

#### 발견된 문제 및 해결

**1. 발현 유전자 필터링 문제**
- **문제**: `filter_expressed_genes()` 함수에서 0개 유전자 반환
- **원인**: `Idents(sobj)` 설정이 `nichenetr::get_expressed_genes()` 호출 전에 필요
- **해결**: `cluster_col` 파라미터 추가 및 Idents 설정 로직 개선

**2. NicheNet 데이터 경로 문제**
- **문제**: NicheNet 데이터 다운로드/읽기 실패
- **해결**: `/data/user3/git_repo/human/` 경로 사용

**3. FindMarkers DEG 필터링 문제**
- **문제**: FindMarkers는 DEG를 찾지만 필터링 후 0개
- **원인**: `p_val_adj`가 모두 1.0 이상 (조정된 p-value가 1.0을 초과할 수 있음)
- **해결**: `p_val_adj_cutoff=1.1`, `logfc_cutoff=0.05`로 완화

#### 테스트 결과

**성공한 테스트**:
- ✅ 데이터 로드 및 검증
- ✅ DEG 추출 (3-10개 DEG)
- ✅ Sender 자동 식별 (23개 sender)
- ✅ 발현 유전자 필터링 (Sender: 9459-11838개, Receiver: 3607-4559개)
- ✅ NicheNet 분석 완료

**최종 분석 결과**:
- **Receiver**: CD4+ T-cells
- **Senders**: 5개 (Monocyte, NK cell, T cell 중심)
- **DEGs**: 2121개 (NicheNet ligand_target_matrix와 교집합)
- **Potential ligands**: 149개
- **Top 10 ligands**: MIF, OSM, CXCL2, SPON2, FURIN, ADAM10, S100A9, TNFSF12, CTSD, NAMPT

### 2025-11-14: 플롯 생성 및 Circos plot 개선

#### 완료된 작업

**1. 모든 플롯 생성 성공**
- Ligand-Target Heatmap (329KB)
- Ligand-Receptor Heatmap (126KB)
- Ligand Activity Histogram (36KB)
- Ligand AUPR Heatmap (84KB)
- **Circos Plot** (PDF: 58KB, PNG: 72KB)

**2. Circos plot 개선 계획**
- Plot title 컨트롤 추가 (NULL일 시 condition, donor/receiver 자동 명시)
- 색깔 범례 명확히 표기 (우상단 등)
- 논문 수준의 퀄리티 향상

## 사용 방법

### 기본 사용법

```r
# 함수 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# DEG 리스트 준비
deg_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = "CD4+ T-cells",
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01)
)

# CCI 분석 실행
results <- run_cci_analysis(
  sobj = sobj,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "CD4+ T-cells",
  sender_clusters = c("Monocytes / Macrophages", "NK Cells"),  # 선택적
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  nichenet_data_dir = "/data/user3/git_repo/human",  # 중요!
  output_dir = "/data/user3/sobj/cci_plots_output",  # 플롯 저장
  run_circos = TRUE,  # Circos plot 생성
  p_val_adj_cutoff = 1.1,
  logfc_cutoff = 0.05
)
```

### 파라미터 설명

#### 필수 파라미터
- `sobj`: Seurat 객체
- `cluster_col`: 클러스터 정보 컬럼명
- `deg_df`: DEG 결과 데이터프레임 (gene, cluster, avg_log2FC, p_val_adj 컬럼 필수)
- `receiver_cluster`: Receiver 클러스터 ID
- `condition_col`: 조건 컬럼명
- `condition_oi`: 관심 조건 값
- `condition_ref`: 참조 조건 값

#### 선택적 파라미터
- `sender_clusters`: Sender 클러스터 벡터 (NULL이면 자동 식별)
- `species`: "human" 또는 "mouse" (기본: "human")
- `output_dir`: 결과 저장 디렉터리 (기본: NULL)
- `run_circos`: Circos plot 생성 여부 (기본: TRUE)
- `nichenet_data_dir`: NicheNet 데이터 디렉터리 (기본: NULL, 자동 다운로드)
- `p_val_adj_cutoff`: DEG p-value 임계값 (기본: 0.05, 완화: 1.1)
- `logfc_cutoff`: logFC 임계값 (기본: 0.25, 완화: 0.05)
- `top_n_ligands`: 분석할 상위 ligand 개수 (기본: 20)

## 출력 구조

### 결과 리스트

```r
results <- list(
  nichenet_results = list(
    ligand_activities = ...,           # Ligand activity 점수
    best_upstream_ligands = ...,       # 우선순위가 높은 ligand 목록
    active_ligand_target_links_df = ..., # Ligand-target 연결
    ligand_receptor_network_df = ...,  # Ligand-receptor 상호작용
    plot_ligand_target_network = ...,  # Ligand-target 히트맵
    plot_ligand_receptor_network = ..., # Ligand-receptor 히트맵
    plot_ligand_activity_hist = ...,   # Activity 히스토그램
    plot_ligand_aupr_heatmap = ...,    # AUPR 히트맵
    plot_circos = ...,                 # Circos plot (recorded plot)
    DE_table_receiver = ...            # Receiver DEG 테이블
  ),
  sender_clusters = ...,              # 사용된 sender 클러스터
  receiver_cluster = ...,              # Receiver 클러스터
  deg_summary = ...,                   # DEG 요약
  saved_path = ...                     # 저장 경로
)
```

### 저장된 파일

`output_dir`을 지정하면 다음 파일들이 저장됩니다:
- `NicheNet_Ligand_Target_Heatmap.png`
- `NicheNet_Ligand_Receptor_Heatmap.png`
- `NicheNet_Ligand_Activity_Histogram.png`
- `NicheNet_Ligand_AUPR_Heatmap.png`
- `NicheNet_Circos_LR.pdf` (고해상도)
- `NicheNet_Circos_LR.png` (웹용)
- `NicheNet_Circos_LR_legend.txt` (범례 정보)

## 주의사항

1. **DEG 데이터프레임 형식**
   - 필수 컬럼: `gene`, `cluster`, `avg_log2FC` (또는 `logFC`), `p_val_adj`
   - `cluster` 컬럼의 값이 `receiver_cluster`와 정확히 일치해야 함

2. **NicheNet 데이터**
   - 첫 실행 시 Zenodo에서 자동 다운로드 (시간 소요)
   - `/data/user3/git_repo/human/` 경로 사용 권장

3. **메모리 사용량**
   - 큰 데이터셋의 경우 메모리 사용량이 많을 수 있음
   - 다운샘플 데이터로 먼저 테스트 권장

4. **FindMarkers 파라미터**
   - `p_val_adj_cutoff`와 `logfc_cutoff`는 데이터에 따라 조정 필요
   - 조정된 p-value가 1.0을 초과할 수 있으므로 `p_val_adj_cutoff=1.1` 권장

## 문제 해결

### 오류: "No DEGs found in receiver"
- `deg_df`에서 `receiver_cluster`와 일치하는 DEG가 있는지 확인
- `cluster` 컬럼 값이 정확히 일치하는지 확인

### 오류: "No potential ligands found"
- Sender 클러스터에 ligand가 발현되는지 확인
- Receiver 클러스터에 receptor가 발현되는지 확인
- `min_pct_expressed` 값을 낮춰보기

### NicheNet 데이터 다운로드 실패
- 인터넷 연결 확인
- 수동으로 데이터 다운로드 후 `nichenet_data_dir` 지정

## 참고 자료

- NicheNet 공식 문서: https://github.com/saeyslab/nichenetr
- 기존 CCI.R 함수: `myR/R/CCI.R`
- 워크트리 표준화 제안: `docs/lds/STANDARDIZATION_PROPOSAL.md`

## 향후 개선 계획

1. **Circos plot 개선**
   - Plot title 자동 생성 (condition, donor/receiver 정보 포함)
   - 색깔 범례 명확히 표기 (우상단 등)
   - 논문 수준의 퀄리티 향상

2. **기능 개선**
   - 조건 없이 DEG 리스트만으로 분석 가능하도록 개선
   - NicheNet 함수 직접 호출 옵션 추가
   - 추가 시각화 기능

3. **문서 보완**
   - 사용 예시 추가
   - 문제 해결 가이드 보완
   - 성능 최적화 가이드

