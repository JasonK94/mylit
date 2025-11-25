# myR Package: Single-Cell RNA-Seq Analysis Integrated Guide

이 문서는 `myR` 패키지의 전체 분석 모듈 통합 가이드입니다. 단일세포 RNA 시퀀싱(scRNAseq) 데이터 분석 파이프라인에서 각 모듈의 역할, 입력/출력, 그리고 모듈 간 관계를 설명합니다.

## 1. Overview (개요)

`myR`는 단일세포 RNA 시퀀싱 데이터를 위한 통합 분석 패키지입니다. 차등 발현 분석(DEG), 세포 간 상호작용(CCI), 차등 풍부도 분석, 환자 수준 분석, trajectory 분석 등 다양한 분석 모듈을 제공합니다.

### 주요 특징
- **통합 파이프라인**: Seurat 객체를 중심으로 한 일관된 데이터 구조
- **다양한 분석 방법론**: Pseudobulk, Mixed Model, Single-cell 등 다양한 DEG 방법론 지원
- **Consensus 분석**: 여러 방법론의 결과를 통합하여 신뢰도 높은 결과 도출
- **시각화**: 논문 수준의 플롯 자동 생성

## 2. Module Overview & Workflow (모듈 개요 및 워크플로우)

### 2.1 전체 분석 파이프라인

```mermaid
%%{init: {'theme':'base', 'themeVariables': {'primaryEdgeColor':'#000000', 'primaryEdgeThickness':4, 'primaryTextColor':'#000000', 'primaryBorderColor':'#000000', 'edgeLabelBackground':'#ffffff', 'tertiaryColor':'#000000'}}}%%
flowchart TD
    Start([scRNAseq Data<br/>Seurat Object])
    
    Start --> Preprocess[Preprocessing<br/>QC, Normalization,<br/>Dimensionality Reduction]
    
    Preprocess --> AnalysisBranch
    Preprocess --> AbundanceBranch
    Preprocess --> TrajectoryBranch
    Preprocess --> PatientLevelBranch
    
    subgraph AnalysisBranch["DEG Analysis Branch"]
        direction TB
        Analysis[analysis Module<br/>NEBULA Mixed Model]
        DegConsensus[deg-consensus Module<br/>Multi-model Consensus]
        LDS[lds Module<br/>Limma-Dream-SVA]
        
        Analysis --> DegConsensus
        LDS --> DegConsensus
        DegConsensus --> DEGResults[("DEG Results<br/>Gene Lists")]
    end
    
    subgraph AbundanceBranch["Differential Abundance Branch"]
        direction TB
        Milo[milo Module<br/>Neighborhood-based DA]
        Milo --> AbundanceResults[("DA Results<br/>Neighborhood Changes")]
    end
    
    subgraph TrajectoryBranch["Trajectory Analysis Branch"]
        direction TB
        Pseudotime[pseudotime Module<br/>Monocle3, Slingshot]
        Pseudotime --> TrajectoryResults[("Trajectory Results<br/>Pseudotime, Lineages")]
    end
    
    subgraph PatientLevelBranch["Patient-Level Analysis Branch"]
        direction TB
        PtUmap[pt.umap Module<br/>Patient Dimensionality Reduction]
        PtUmap --> PatientResults[("Patient-Level Features<br/>PCA, UMAP, PERMANOVA")]
    end
    
    DEGResults --> CCI
    DEGResults --> FGS
    
    subgraph CCI["Cell-Cell Interaction"]
        direction TB
        CCIAnalysis[cci Module<br/>NicheNet Analysis]
        CCIAnalysis --> CCIResults[("CCI Results<br/>Ligand-Receptor Networks<br/>Circos Plots")]
    end
    
    subgraph FGS["Gene Signature Discovery"]
        direction TB
        FGSModule[fgs Module<br/>Find Gene Signature<br/>Random Forest, Lasso, GAM, etc.]
        FGSModule --> SignatureResults[("Gene Signatures<br/>L1 Signatures<br/>TML Meta-Learner")]
    end
    
    DEGResults --> Output[("Final Outputs<br/>.qs files, Plots")]
    AbundanceResults --> Output
    TrajectoryResults --> Output
    PatientResults --> Output
    CCIResults --> Output
    SignatureResults --> Output
    
    style Start fill:#e1f5ff,stroke:#0066cc,stroke-width:3px
    style Preprocess fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style AnalysisBranch fill:#d1ecf1,stroke:#17a2b8,stroke-width:2px
    style AbundanceBranch fill:#cfe2ff,stroke:#0d6efd,stroke-width:2px
    style TrajectoryBranch fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style PatientLevelBranch fill:#e2e3e5,stroke:#6c757d,stroke-width:2px
    style CCI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style FGS fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style Output fill:#d4edda,stroke:#28a745,stroke-width:3px
```

### 2.2 모듈별 상세 정보

| Module | Purpose | Input | Output | Key Methods |
|--------|---------|-------|--------|-------------|
| **analysis** | Mixed-Effects Model DEG | Seurat Object, Formula | DEG Results (.qs) | NEBULA, muscat |
| **deg-consensus** | Multi-model DEG Consensus | Seurat Object, Multiple DEG Results | Consensus DEG List | limma, edgeR, DESeq2, muscat, nebula, dream |
| **lds** | SVA-corrected LMM DEG | Seurat/DGEList, Formula | MArrayLM Object, SV Matrix | limma-dream, SVA |
| **milo** | Differential Abundance | Seurat Object | DA Results (.qs) | MiloR, GLM |
| **cci** | Cell-Cell Interaction | Seurat Object, DEG List | Ligand-Receptor Networks | NicheNet |
| **fgs** | Gene Signature Discovery | Seurat Object, Target Variable | Gene Signatures, TML Model | Random Forest, Lasso, GAM, PCA, NMF |
| **pt.umap** | Patient-Level Analysis | Seurat Object | Patient Features, PCA/UMAP | CLR Transform, Signature Scores |
| **pseudotime** | Trajectory Inference | Seurat Object | Pseudotime, Lineages | Monocle3, Slingshot, GAM |

## 3. Module Details (모듈 상세)

### 3.1 DEG Analysis Modules

#### analysis Module
- **목적**: 환자 수준의 변동성을 고려한 Mixed-Effects Model 기반 DEG 분석
- **입력**: Seurat Object, Formula (예: `~ g3 + sex + anno3.scvi + (1|GEM/hos_no)`)
- **출력**: NEBULA 결과 요약 (.qs), 통계 테이블
- **특징**: Nested Random Effect 지원, 복잡한 실험 설계 지원
- **상세 문서**: `docs/analysis/ANALYSIS_INTEGRATED_GUIDE.md`

#### deg-consensus Module
- **목적**: 여러 DEG 방법론의 결과를 통합하여 Consensus DEG 리스트 생성
- **입력**: Seurat Object 또는 사전 계산된 DEG 리스트들
- **출력**: Consensus DEG List, Agreement Scores, Visualization Plots
- **특징**: 10개 이상의 방법론 지원, Agreement Score 기반 필터링
- **상세 문서**: `docs/deg-consensus-dev/DEG_CONSENSUS_INTEGRATED_GUIDE.md`

#### lds Module
- **목적**: SVA를 이용한 숨겨진 공변량 탐지 및 Limma-Dream 기반 LMM DEG 분석
- **입력**: Seurat/DGEList/Matrix, Formula, 메타데이터
- **출력**: MArrayLM 객체, SV 행렬, 상관관계 히트맵
- **특징**: GeoMx 데이터에 적합, 잔차 분산 기반 SV 탐지
- **상세 문서**: `docs/lds/LDS_INTEGRATED_GUIDE.md`

### 3.2 Differential Abundance Module

#### milo Module
- **목적**: Neighborhood 기반 차등 풍부도(Differential Abundance) 분석
- **입력**: Seurat Object, Patient/Sample ID, Cluster ID, Target Variable
- **출력**: DA Results (.qs), UMAP Plots, Beeswarm Plots
- **특징**: kNN 그래프 기반 Neighborhood 생성, Block Method 지원
- **상세 문서**: `docs/milo/MILO_INTEGRATED_GUIDE.md`

### 3.3 Trajectory Analysis Module

#### pseudotime Module
- **목적**: 세포 분화 경로 추론 및 Pseudotime 기반 유전자 동역학 분석
- **입력**: Seurat Object, 메타데이터 (시간 변수 등)
- **출력**: Pseudotime, Lineages, Gene Dynamics Results
- **특징**: Monocle3, Slingshot 지원, GAM 기반 동역학 모델링
- **상세 문서**: `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE.md`

### 3.4 Patient-Level Analysis Module

#### pt.umap Module
- **목적**: 셀 레벨 데이터를 환자 수준으로 집계하여 차원 축소 및 이상 탐지
- **입력**: Seurat Object, Cluster ID, Patient ID
- **출력**: Patient-Level Features, PCA/UMAP, PERMANOVA Results
- **특징**: Frequency View, Signature View, Latent View 지원
- **상세 문서**: `docs/pt.umap/PT_UMAP_INTEGRATED_GUIDE.md`

### 3.5 Cell-Cell Interaction Module

#### cci Module
- **목적**: NicheNet 기반 세포 간 상호작용 분석 (Ligand-Receptor-Target 네트워크)
- **입력**: Seurat Object, DEG List (Receiver 클러스터)
- **출력**: Ligand-Receptor Networks, Circos Plots, Heatmaps
- **특징**: DEG 리스트 직접 입력 지원, 자동 Sender 식별
- **상세 문서**: `docs/cci/CCI_INTEGRATED_GUIDE.md`

### 3.6 Gene Signature Discovery Module

#### fgs Module
- **목적**: 그룹을 구분하는 Gene Signature 발견 및 Meta-Learner 구축
- **입력**: Seurat Object, Target Variable, Control Variables
- **출력**: L1 Signatures (각 방법론별), TML Model (L2), Gene Importance
- **특징**: Random Forest, Lasso, GAM, PCA, NMF 등 다양한 방법론 지원
- **상세 문서**: `docs/fgs/FGS_TML_INTEGRATED_GUIDE.md`

## 4. Typical Analysis Workflows (일반적인 분석 워크플로우)

### 4.1 DEG Analysis Workflow
```
1. Preprocessing (QC, Normalization)
   ↓
2. analysis Module (NEBULA) 또는 lds Module (Limma-Dream-SVA)
   ↓
3. deg-consensus Module (여러 방법론 통합)
   ↓
4. DEG 리스트 활용:
   - cci Module (Cell-Cell Interaction 분석)
   - fgs Module (Signature Discovery)
```

### 4.2 Comprehensive Analysis Workflow
```
1. Preprocessing
   ↓
2. Parallel Analysis:
   - DEG Analysis (analysis, deg-consensus)
   - Differential Abundance (milo)
   - Trajectory Analysis (pseudotime)
   - Patient-Level Analysis (pt.umap)
   ↓
3. Integration:
   - CCI Analysis (cci) using DEG results
   - Signature Discovery (fgs) using DEG results
   ↓
4. Final Outputs: .qs files, Plots
```

## 5. Data Format & Requirements (데이터 형식 및 요구사항)

### 공통 입력 형식
- **Seurat Object**: `.qs` 포맷으로 저장된 Seurat 객체
- **메타데이터 컬럼**:
  - `hos_no` 또는 `patient_id`: 환자/샘플 식별자
  - `anno3.scvi` 또는 `cluster_id`: 클러스터 정보
  - `g3` 또는 `target_var`: 주요 그룹 변수
  - 기타 공변량: `sex`, `age`, `batch` 등

### 출력 형식
- **결과 객체**: `.qs` 포맷 (R 객체 저장)
- **플롯**: `.png`, `.pdf` 포맷
- **통계 테이블**: CSV 또는 R 데이터프레임

## 6. Quick Start (빠른 시작)

### 기본 DEG 분석
```r
# 1. Load package
devtools::load_all("/path/to/myR")

# 2. Load data
sobj <- qs::qread("path/to/seurat_object.qs")

# 3. Run analysis (NEBULA)
source("scripts/analysis/run_formula1_analysis.R")

# 4. Run consensus
source("scripts/deg-consensus-dev/run_consensus_simple.R")
```

### CCI Analysis
```r
# 1. Prepare DEG list (from analysis or deg-consensus)
deg_df <- read.csv("deg_results.csv")

# 2. Run CCI analysis
source("myR/R/cci/run_cci_analysis.R")
results <- run_cci_analysis(
  sobj = sobj,
  deg_df = deg_df,
  receiver_cluster = "CD4+ T-cells"
)
```

## 7. Module Documentation (모듈 문서)

각 모듈의 상세한 사용법, 파라미터, 예제는 다음 문서를 참조하세요:

| Module | English Guide | Korean Guide |
|--------|---------------|--------------|
| **analysis** | `docs/analysis/ANALYSIS_INTEGRATED_GUIDE.md` | `docs/analysis/ANALYSIS_INTEGRATED_GUIDE_KR.md` |
| **deg-consensus** | `docs/deg-consensus-dev/DEG_CONSENSUS_INTEGRATED_GUIDE.md` | `docs/deg-consensus-dev/DEG_CONSENSUS_INTEGRATED_GUIDE_KR.md` |
| **lds** | `docs/lds/LDS_INTEGRATED_GUIDE.md` | `docs/lds/LDS_INTEGRATED_GUIDE_KR.md` |
| **milo** | `docs/milo/MILO_INTEGRATED_GUIDE.md` | `docs/milo/MILO_INTEGRATED_GUIDE_KR.md` |
| **cci** | `docs/cci/CCI_INTEGRATED_GUIDE.md` | `docs/cci/CCI_INTEGRATED_GUIDE_KR.md` |
| **fgs** | `docs/fgs/FGS_TML_INTEGRATED_GUIDE.md` | `docs/fgs/FGS_TML_INTEGRATED_GUIDE_KR.md` |
| **pt.umap** | `docs/pt.umap/PT_UMAP_INTEGRATED_GUIDE.md` | `docs/pt.umap/PT_UMAP_INTEGRATED_GUIDE_KR.md` |
| **pseudotime** | `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE.md` | `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE_KR.md` |

## 8. Best Practices (모범 사례)

1. **데이터 검증**: 각 모듈 실행 전 필수 메타데이터 컬럼 확인
2. **결과 저장**: 중간 결과는 `.qs` 포맷으로 저장하여 재사용
3. **메모리 관리**: 대규모 데이터셋 분석 시 다운샘플링으로 먼저 테스트
4. **Consensus 활용**: 단일 방법론보다 `deg-consensus`를 통한 통합 결과 사용 권장
5. **문서 참조**: 각 모듈의 통합 가이드에서 상세 파라미터 및 주의사항 확인

## 9. Troubleshooting (문제 해결)

### 공통 이슈
- **메모리 부족**: 데이터 다운샘플링 또는 배치 처리
- **의존성 오류**: `renv::activate()` 후 패키지 재설치
- **결과 불일치**: 각 모듈의 버전 및 파라미터 확인

### 모듈별 이슈
각 모듈의 통합 가이드에서 "Critical Warnings" 섹션을 참조하세요.

## 10. References (참고 자료)

- **Seurat**: https://satijalab.org/seurat/
- **NEBULA**: https://github.com/lhe17/nebula
- **NicheNet**: https://github.com/saeyslab/nichenetr
- **MiloR**: https://bioconductor.org/packages/MiloR/
- **Monocle3**: https://cole-trapnell-lab.github.io/monocle3/
- **Slingshot**: https://bioconductor.org/packages/Slingshot/

---

**Last Updated**: 2025-01-XX  
**Package Version**: See `myR/DESCRIPTION`


