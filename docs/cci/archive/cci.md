# Cell-to-Cell Interaction (CCI) Analysis Tool

## 개요
scRNAseq 데이터에서 Cell-to-Cell Interaction 분석을 수행하는 통합 도구입니다. NicheNet을 중심으로 하여 ligand-receptor 상호작용을 분석하고, DEG 리스트를 직접 입력받아 receiver cell type의 변화에 기여하는 sender cell type을 식별합니다.

## 핵심 개념

### 1. Sender-Receiver 모델
- **Sender**: 신호를 보내는 세포 유형 (ligand 발현)
- **Receiver**: 신호를 받는 세포 유형 (receptor 발현 및 DEG 발생)
- **Ligand-Receptor Pair**: NicheNet 데이터베이스 기반 상호작용

### 2. DEG 기반 분석
- Receiver cell type에서 조건별로 발생한 DEG를 분석
- DEG 리스트를 직접 입력받아 NicheNet 분석에 활용
- Sender cell type의 ligand activity를 예측하여 우선순위 결정

## 주요 함수

### `run_cci_analysis()`
CCI 분석의 메인 함수입니다.

**입력**:
- `sobj`: Seurat 객체
- `cluster_col`: 클러스터 컬럼명 (예: "anno3.scvi")
- `deg_df`: DEG 결과 데이터프레임
- `receiver_cluster`: receiver 클러스터 ID
- `sender_clusters`: sender 클러스터 벡터 (선택적, NULL이면 자동 식별)
- `condition_col`: 조건 컬럼 (예: "g3")
- `condition_oi`: 관심 조건 (예: "2")
- `condition_ref`: 참조 조건 (예: "1")
- `species`: "human" 또는 "mouse" (기본: "human")
- `output_dir`: 출력 디렉터리 (기본: NULL)

**출력**:
- NicheNet 분석 결과 리스트
- Sender-receiver 매핑 정보
- DEG 요약 정보
- 저장 경로

## 사용 예시

### 기본 사용법
```r
# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# DEG 리스트 준비 (예시)
deg_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = c("Cluster1", "Cluster1", "Cluster1"),
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01)
)

# CCI 분석 실행
results <- run_cci_analysis(
  sobj = sobj,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Cluster1",
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human"
)
```

### 다중 Sender 분석
```r
# 여러 sender 클러스터 지정
results <- run_cci_analysis(
  sobj = sobj,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Cluster1",
  sender_clusters = c("Cluster2", "Cluster3", "Cluster4"),
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1"
)
```

## 파라미터 설명

### 필수 파라미터
- **sobj**: 분석할 Seurat 객체
- **cluster_col**: 클러스터 정보가 담긴 메타데이터 컬럼명
- **deg_df**: DEG 결과 데이터프레임 (gene, cluster, logFC, p_val_adj 컬럼 필수)
- **receiver_cluster**: 분석할 receiver 클러스터 ID
- **condition_col**: 비교할 조건이 담긴 메타데이터 컬럼명
- **condition_oi**: 관심 조건 값
- **condition_ref**: 참조 조건 값

### 선택적 파라미터
- **sender_clusters**: sender 클러스터 벡터. NULL이면 자동 식별
- **species**: "human" 또는 "mouse" (기본: "human")
- **output_dir**: 결과 저장 디렉터리 (기본: NULL, 저장 안 함)
- **min_pct_expressed**: 발현 비율 임계값 (기본: 0.10)
- **p_val_adj_cutoff**: DEG p-value 임계값 (기본: 0.05)
- **logfc_cutoff**: logFC 임계값 (기본: 0.25)
- **top_n_ligands**: 분석할 상위 ligand 개수 (기본: 20)

## 출력 구조

### 결과 리스트 구성
```r
results <- list(
  nichenet_results = ...,      # NicheNet 분석 결과
  sender_receiver_map = ...,    # Sender-receiver 매핑
  deg_summary = ...,           # 사용된 DEG 요약
  output_path = ...            # 저장 경로
)
```

### NicheNet 결과 구조
- `ligand_activities`: Ligand activity 점수
- `best_upstream_ligands`: 우선순위가 높은 ligand 목록
- `active_ligand_target_links_df`: Ligand-target 연결 정보
- `ligand_receptor_network_df`: Ligand-receptor 상호작용
- `plot_ligand_target_network`: Ligand-target 네트워크 히트맵
- `plot_ligand_receptor_network`: Ligand-receptor 네트워크 히트맵
- `plot_circos`: Circos plot (선택적)

## 주의사항

1. **DEG 데이터프레임 형식**: 
   - 필수 컬럼: `gene`, `cluster`, `logFC` (또는 `avg_log2FC`), `p_val_adj`
   - `cluster` 컬럼의 값이 `receiver_cluster`와 일치해야 함

2. **메타데이터 요구사항**:
   - `cluster_col`에 지정한 컬럼이 존재해야 함
   - `condition_col`에 지정한 컬럼이 존재해야 함
   - 조건 값에 NA가 있으면 제외됨

3. **NicheNet 데이터**:
   - 첫 실행 시 Zenodo에서 자동 다운로드
   - 이후에는 로컬 캐시 또는 전역 환경에서 로드

4. **메모리 사용량**:
   - 큰 데이터셋의 경우 메모리 사용량이 많을 수 있음
   - 다운샘플 데이터로 먼저 테스트 권장

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
- 기존 nichenet_analysis.R: `myR/R/analysis/nichenet_analysis.R`

