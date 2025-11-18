# CCI (Cell-to-Cell Interaction) Analysis Tool

## 개요
scRNAseq 데이터에서 Cell-to-Cell Interaction 분석을 수행하는 통합 도구입니다. NicheNet을 중심으로 하여 ligand-receptor 상호작용을 분석하고, DEG 리스트를 직접 입력받아 receiver cell type의 변화에 기여하는 sender cell type을 식별합니다.

## 빠른 시작

### 1. 환경 설정
```r
# 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/cci/myR")

# 또는 함수 소스 직접 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")

# run_nichenet_analysis가 들어있는 CCI.R은 워크트리 파일을 우선 사용
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  stop("CCI.R not found in worktree or main repository.")
}
```

### 2. 기본 사용법
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

## 프로젝트 구조

```
_wt/cci/
  docs/
    cci/
      cci.md                    # 메인 기능 문서
      TEST_INSTRUCTIONS.md      # 테스트 가이드
      devlog.md                 # 개발 로그
      CCI_DEVELOPMENT_PLAN.md   # 개발 계획
  scripts/
    cci/
      test_cci.R                # 테스트 스크립트
  myR/
    R/
      cci/
        run_cci_analysis.R      # 메인 CCI 분석 함수
        prepare_cci_data.R      # 데이터 준비 및 검증
        nichenet_wrapper.R      # NicheNet 분석 래퍼
        save_cci_results.R      # 결과 저장 유틸리티
        utils_cci.R             # CCI 유틸리티 함수
```

> `guide.md` 기준: CCI 워크트리에서 사용하는 스크립트는 `scripts/cci/`, 문서는 `docs/cci/`에만 추가합니다.

## 주요 함수

### `run_cci_analysis()`
CCI 분석의 메인 함수입니다. Seurat 객체, 클러스터 정보, DEG 리스트를 입력받아 NicheNet 분석을 수행합니다.

**필수 파라미터**:
- `sobj`: Seurat 객체
- `cluster_col`: 클러스터 컬럼명
- `deg_df`: DEG 데이터프레임
- `receiver_cluster`: Receiver 클러스터 ID
- `condition_col`: 조건 컬럼명
- `condition_oi`: 관심 조건

**선택적 파라미터**:
- `sender_clusters`: Sender 클러스터 벡터 (NULL이면 자동 식별)
- `species`: "human" 또는 "mouse" (기본: "human")
- `auto_save`: 자동 저장 여부 (기본: TRUE)

## 데이터 요구사항

### DEG 데이터프레임 형식
필수 컬럼:
- `gene`: 유전자명
- `cluster`: 클러스터 ID
- `avg_log2FC` 또는 `logFC`: 로그 폴드 체인지
- `p_val_adj` 또는 `FDR`: 조정된 p-value

### Seurat 객체 메타데이터
필수 컬럼:
- 클러스터 정보 컬럼 (예: `anno3.scvi`)
- 조건 정보 컬럼 (예: `g3`)

### Receiver DEG 재사용
- `run_nichenet_analysis()`는 `receiver_de_table`을 직접 전달받아 `FindMarkers()` 재실행 없이 바로 NicheNet을 돌릴 수 있습니다.
- 컬럼명이 다르면 `receiver_gene_col`, `receiver_logfc_col`, `receiver_pval_col`로 매핑하세요.
- `run_cci_analysis()`는 자체적으로 receiver DEG를 추출하여 동일한 방식으로 전달하므로, 동일 receiver를 반복 실행하더라도 중복 계산이 발생하지 않습니다.

## 출력

결과는 리스트 형태로 반환되며, 다음을 포함합니다:
- `nichenet_results`: NicheNet 분석 결과
- `sender_receiver_map`: Sender-receiver 매핑
- `deg_summary`: DEG 요약 정보
- `saved_path`: 저장된 파일 경로 (auto_save=TRUE인 경우)

## 테스트

테스트 스크립트 실행:
```r
source("/home/user3/data_user3/git_repo/_wt/cci/scripts/cci/test_cci.R")
```

## 문서

- **상세 문서**: `docs/cci/cci.md`
- **모듈별 상세 가이드**: `docs/cci/CCI_module.md`
- **테스트 가이드**: `docs/cci/TEST_INSTRUCTIONS.md`
- **개발 로그**: `docs/cci/devlog.md`

## 참고 자료

- NicheNet 공식 문서: https://github.com/saeyslab/nichenetr
- CCI.R 함수 (워크트리 우선): `/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R` (없으면 `/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R`)

