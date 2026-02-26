# Cell-to-Cell Interaction (CCI) Analysis 통합 가이드

이 문서는 Cell-to-Cell Interaction (CCI) 분석 모듈의 통합 가이드입니다. NicheNet을 기반으로 하여 ligand-receptor 상호작용을 분석하고, 사전 계산된 DEG 리스트를 활용하여 유연한 분석을 지원합니다.

## 1. 소개 (Introduction)

### 목적
단일세포 RNA 시퀀싱(scRNAseq) 데이터에서 세포 간 상호작용을 분석하여, 특정 조건(예: 질병 vs 정상)에서 **Receiver 세포**의 유전자 발현 변화(DEG)를 유도하는 **Sender 세포**와 **Ligand**를 식별합니다.

### 핵심 기능
1.  **DEG 리스트 직접 입력**: `FindMarkers`, `limma`, `edgeR` 등 다양한 방법으로 계산된 DEG 리스트를 직접 입력받아 분석합니다.
2.  **자동 Sender 식별**: Receiver를 제외한 모든 클러스터를 자동으로 Sender로 설정하거나, 사용자가 지정할 수 있습니다.
3.  **NicheNet 통합**: NicheNet 데이터베이스를 활용하여 Ligand-Receptor-Target 유전자 네트워크를 구축합니다.
4.  **논문 수준 시각화**: Circos plot, 히트맵 등 다양한 시각화 결과물을 자동 생성합니다.

### 주요 용어
*   **Sender**: 신호(Ligand)를 보내는 세포 유형.
*   **Receiver**: 신호(Ligand)를 받아 반응하는 세포 유형. Target Gene의 발현 변화(DEG)가 관찰되는 세포입니다.
*   **Ligand**: Sender 세포에서 분비되어 Receiver 세포의 Receptor와 결합하는 물질.
*   **Receptor**: Receiver 세포 표면에 존재하며 Ligand와 결합하는 단백질.
*   **Target Gene**: Ligand-Receptor 상호작용에 의해 발현이 조절되는 하위 유전자 (DEG).

## 2. 워크플로우 시각화 (Workflow Visualization)

### 2.1. Basic NicheNet Workflow (Standard/Single-Condition)

이 워크플로우는 **단일 조건** 또는 **단순 비교**(예: Normal vs Tumor)에서 특정 Receiver 세포의 DEG를 설명하는 Ligand를 찾을 때 사용됩니다.

```mermaid
%%{init: {'theme':'base', 'themeVariables': {'primaryEdgeColor':'#333', 'primaryTextColor':'#000', 'lineColor':'#333'}}}%%
flowchart TD
    %% Nodes
    Start([Start: Seurat Object])
    
    subgraph DataPrep ["1. Data Preparation"]
        direction TB
        CalcDEG["DEG Calculation<br/>Method: Seurat FindMarkers (Wilcoxon/MAST)<br/>Comparison: Condition vs Control"]
        DefineReceiver["Define Receiver & Sender<br/>Receiver: Cell type of interest<br/>Sender: All or selected types"]
        ExtractExpressed["Filter Expressed Genes<br/>Threshold: min_pct > 0.10"]
    end

    subgraph NicheNetCore ["2. NicheNet Core Algorithm"]
        direction TB
        LoadModel[("NicheNet Model<br/>Ligand-Target Matrix<br/>(Weighted Regulatory Potential)")]
        
        CalcActivity["Calculate Ligand Activity<br/>Method: Pearson Correlation<br/>(Target Potential vs. DEG Status)"]
        
        PredictTargets["Infer Target Genes<br/>Identify DEGs with high<br/>regulatory potential from top ligands"]
        
        PrioritizePair["Ligand-Receptor Scoring<br/>Combine Activity + Expression"]
    end
    
    subgraph Output ["3. Outputs"]
        Heatmaps["Heatmaps<br/>1. Ligand Activity<br/>2. Ligand-Target Links"]
        Circos["Circos Plot<br/>(Simple Linkage)"]
    end

    %% Edges
    Start --> DataPrep
    CalcDEG -->|Receiver DEGs| CalcActivity
    DefineReceiver --> ExtractExpressed
    ExtractExpressed -->|Expressed Ligands/Receptors| PrioritizePair
    
    LoadModel --> CalcActivity
    CalcActivity -->|Activity Score| PredictTargets
    PredictTargets --> PrioritizePair
    
    PrioritizePair --> Output
    
    style Start fill:#f9f,stroke:#333
    style LoadModel fill:#eee,stroke:#999,stroke-dasharray: 5 5
    style CalcActivity fill:#d1ecf1,stroke:#17a2b8
    style Output fill:#d4edda,stroke:#28a745
```

### 2.2. MultiNicheNet Workflow (Complex/Multi-Sample)

이 워크플로우는 **Multi-Sample/Multi-Group** 데이터셋(IS2 vs IS3 등)에서 샘플 간 변동성을 고려한 통계적 비교 분석을 수행합니다. 현재 개발 중인 고도화 파이프라인입니다.

```mermaid
%%{init: {'theme':'base', 'themeVariables': {'primaryEdgeColor':'#333', 'primaryTextColor':'#000', 'lineColor':'#333'}}}%%
flowchart TD
    %% Nodes
    StartMNN(["Start: Seurat Object<br/>(Multi-Sample)"])
    
    subgraph MNNPrep ["1. Pre-processing"]
        ToSCE["Convert to SingleCellExperiment (SCE)"]
        Pseudobulk["Pseudobulk Aggregation<br/>Sum counts per Sample-CellType"]
    end
    
    subgraph DEAnalysis ["2. Differential Expression (muscat)"]
        MuscatRun["Run muscat (edgeR/limma-voom)<br/>Correction for Sample ID"]
        SaveDE[("Save: celltype_de.rds<br/>(LogFC, p-val per contrast)")]
    PrepData --> Sender
    
    subgraph Signal[Signaling Analysis]
        MNN[MultiNicheNet]
        CC[CellChat]
    end

    Prep --> MNN
    Prep --> CC
    
    subgraph Sender["3. Sender 식별"]
        direction TB
        CheckSender{sender_clusters<br/>지정됨?}
        CheckSender -->|Yes| UseSpecified["지정된 Sender 사용"]
        CheckSender -->|No| AutoIdentify["자동 식별<br/>Receiver 제외<br/>모든 클러스터"]
        UseSpecified --> SenderList[("Sender 클러스터 목록")]
        AutoIdentify --> SenderList
    end
    
    subgraph MNNAlg ["3. MultiNicheNet Core"]
        Loadprior["Load Ligand-Target Matrix"]
        PredictAct["Predict Ligand Activities<br/>Per Contrast & Receiver"]
        SaveAct[("Save: ligand_activities.rds")]
        
        ConstructNetwork["Construct Network"]
    end
    
    subgraph Scoring ["4. Prioritization & Scoring"]
        IntegrateScore["Calculate Prioritization Score<br/>Weighted Sum of:<br/>- Activity (Scaled)<br/>- Ligand/Receptor Expression<br/>- DE LogFC (Concordance)"]
        Dedup["Deduplication & Sorting<br/>1. Sort by Score (Desc)<br/>2. Distinct(L-R-S-R)"]
        SaveTable[("Save: prioritization_tables.rds")]
    end
    
    subgraph Viz ["5. Advanced Visualization"]
        CircosComp["Comparative Circos Plot (V10)<br/>- Inner Hint Tracks<br/>- Deterministic Layout<br/>- Deduped Interactions"]
    end

    %% Edges
    StartMNN --> MNNPrep
    ToSCE --> Pseudobulk
    Pseudobulk --> MuscatRun
    MuscatRun --> SaveDE
    SaveDE --> PredictAct
    
    Loadprior --> PredictAct
    PredictAct --> SaveAct
    SaveAct --> ConstructNetwork
    ConstructNetwork --> IntegrateScore
    
    IntegrateScore --> SaveTable
    SaveTable --> Dedup
    Dedup --> CircosComp


    %% 1. 스타일 클래스 정의 (따옴표 없이 작성)
    classDef mnnStyle fill:#f9f,stroke:#333
    classDef orangeStyle fill:#ffedcc,stroke:#d69e2e
    classDef greenStyle fill:#d4edda,stroke:#28a745

    %% 2. 노드에 클래스 적용
    class StartMNN mnnStyle
    class SaveDE,SaveAct,SaveTable orangeStyle
    class CircosComp greenStyle

```

## 3. CellChat 분석 워크플로우 (CellChat Workflow)

CellChat 분석은 다음 두 가지 전략 중 하나를 선택하여 수행할 수 있습니다. 각 단계별 상세 로직은 아래와 같습니다.

```mermaid
%%{init: {'theme':'base', 'themeVariables': {'primaryEdgeColor':'#000000', 'primaryEdgeThickness':4, 'primaryTextColor':'#000000', 'primaryBorderColor':'#000000', 'edgeLabelBackground':'#ffffff', 'tertiaryColor':'#000000'}}}%%
graph TD
    %% 노드 스타일 정의
    classDef step fill:#f9f,stroke:#333,stroke-width:2px;
    classDef logic fill:#e1f5fe,stroke:#0277bd,stroke-width:2px,stroke-dasharray: 5 5;
    classDef cluster fill:#d7ccc8,stroke:#5d4037,stroke-width:2px;
    classDef error fill:#ffcdd2,stroke:#c62828,stroke-width:2px;

    %% 시작
    Start([Seurat Object Input]) --> ChooseStrat{Strategy Selection}
    
    %% 전략 선택
    ChooseStrat --"Method 1: Pooled (Recommended)<br/>Power 우선"--> Pool[Split by Condition]
    ChooseStrat --"Method 2: Sample-wise (Strict)<br/>Rigor 우선"--> Split[Split by Sample]
    
    %% Pooled Analysis
    subgraph Pooled["Method 1: Pooled Analysis (High Power)"]
        direction TB
        Pool --> |"g3='2' (Stroke)"| Rep1[Rep1: All Stroke Cells]
        Pool --> |"g3='1' (Control)"| Rep2[Rep2: All Control Cells]
        
        Rep1 --> Step1P["Step 1: Preprocessing"]
        Step1P --> DEG_P["Identify DEG & Interactions"]
        class DEG_P logic
        DEG_P --"Wilcoxon Rank Sum Test<br/>(P < 0.05)"--> FilterP["Filter Interactions"]
        
        FilterP --> Step2P["Step 2: Probability"]
        Step2P --> Prob_P["computeCommunProb"]
        class Prob_P logic
        Prob_P --"Type: triMean<br/>Permutation Test (P < 0.05)"--> MatP[Comm Probability Matrix]
    end

    %% Sample-wise Analysis
    subgraph SampleWise["Method 2: Sample-wise Analysis (Rigorous)"]
        direction TB
        Split --> |Sample 1...N| Indiv[Individual Analysis]
        Indiv --> Step1S["Step 1: Preprocessing"]
        Step1S --> DEG_S["Identify DEG"]
        
        DEG_S --> Step2S["Step 2: Probability"]
        Step2S --> Prob_S["computeCommunProb"]
        
        Prob_S --> Merge[Merge by Condition]
        class Merge logic
        Merge --"mergeCellChat()"--> MergedObj[Merged Object]
    end

    %% 비교 분석
    MatP --> Compare["Step 3: Comparison Analysis"]
    MergedObj --> Compare
    
    subgraph Comparison["Step 3: Differential Analysis"]
        direction TB
        Compare --> DiffNet["Differential Network"]
        DiffNet --"Circle Plot (Red/Blue)<br/>netVisual_diffInteraction"--> Vis1[Changed Interactions]
        
        Compare --> InfoFlow["Information Flow"]
        InfoFlow --"Stacked Bar<br/>rankNet"--> Vis2[Signaling Pathway Changes]
        
        Compare --> DiffStrength["Interaction Strength"]
    end
    
    class Start,ChooseStrat step
    class DiffNet,InfoFlow,DiffStrength logic
```

### 상세 방법론 (Methodology Details)

1.  **전략 선택 (Strategy Selection)**
    *   **Pooled Analysis**: 조건별로 모든 샘플의 세포를 합쳐서 분석. Power가 높으나 False Positive 위험이 있어 `min.cells`를 높여(50+) 제어. **권장**.
    *   **Sample-wise Analysis**: 각 환자별로 독립 분석 후 병합. 엄밀하지만 세포 수가 적을 경우 Interaction이 검출되지 않을 수 있음(Zero Interaction).

2.  **DEG & Interaction 식별 (`identifyOverExpressedGenes`)**
    *   **Method**: Wilcoxon Rank Sum Test (기본값)
    *   이 함수는 단순한 발현량 Thresholding(임계값 설정)만이 아니라, 통계적 검정(Wilcoxon Test)을 수행하여 유의미한 과발현 유전자를 식별합니다.
    *   **Thresholding**: 계산된 P-value(`thresh.p`)와 Fold Change(`thresh.fc`), 발현 비율(`thresh.pc`) 등을 기준으로 필터링합니다.
    *   **P-value 이슈**: 대규모 데이터셋에서 Adjusted P-value (Bonferroni 등) 계산 시, $P_{adj} = P_{raw} \times N$ 공식에 의해 1.0을 초과하는 값이 나올 수 있습니다. 이는 수학적으로 가능한 값이며, 보통 1.0으로 간주하면 됩니다.

3.  **확률 계산 (`computeCommunProb`)**
    *   **Expression Value**: Trimean 사용 ($Avg = (Q1 + 2Q2 + Q3) / 4$)하여 이상치 영향 최소화.
    *   **Hill Function**: 리간드-수용체 결합 확률을 모델링.
        $$ P_{ij} = \frac{L_i \cdot R_j}{K_h + L_i \cdot R_j} $$
        ($L_i$: Sender의 리간드 발현량, $R_j$: Receiver의 수용체 발현량, $K_h$: 해리 상수)
    *   **Permutation Test**: 통계적 유의성 검증.
        *   **과정**: Sender와 Receiver의 세포 라벨(Identity)을 무작위로 섞음(Shuffle).
        *   **검정**: 무작위 상태에서의 확률 분포(Null Distribution) 생성 후, 관측된 확률이 상위 5% 안에 드는지 확인 ($P < 0.05$).

4.  **비교 분석 (`netVisual_diffInteraction`)**
    *   **Count**: 상호작용의 개수 (유의미한 L-R 쌍의 수).
    *   **Weight (Interaction Strength)**: 상호작용 확률(Communication Probability)의 합 또는 강도입니다. 단순한 확률(Probability) 자체가 아니라, 해당 그룹 간의 모든 유의미한 L-R Pair의 통신 확률을 누적(Accumulate)한 값으로 "상호작용 강도"를 의미합니다.
    *   **Difference Calculation**: 두 조건 간의 Count 또는 Weight 차이를 계산.
        $$ \Delta_{weight} = Weight_{Condition2} - Weight_{Condition1} $$
    *   **Red Edge**: $\Delta > 0$ (증가) / **Blue Edge**: $\Delta < 0$ (감소).

## 4. 개발 로그 및 개선사항 (Development Log & Improvements)

### 주요 변경 사항
*   **v1.0 (2025-11-14)**: 초기 구현. DEG 리스트 직접 입력 지원, 모듈화된 구조(준비, 분석, 저장) 구축.
*   **v1.1 (2025-12-16)**: 파일 구조 리팩토링 및 CellChat v2 호환성 업데이트. NicheNet 관련 파일 명확화.
*   **최적화**:
    *   `receiver_de_table` 재사용: 동일 Receiver에 대해 DEG 계산 반복 방지.
    *   메모리 관리: 중간 결과 저장 후 대용량 객체 정리(`rm`, `gc`).
    *   컬럼 매핑: 다양한 DEG 테이블 포맷(`avg_log2FC` vs `logFC` 등) 호환 지원.

### 해결된 이슈
*   **발현 유전자 필터링 오류**: `Idents(sobj)` 설정 누락으로 인한 0개 유전자 반환 문제 해결.
*   **DEG 필터링 문제**: 조정된 p-value가 1.0을 초과하는 경우에 대비하여 cutoff 완화(`1.1`).
*   **NicheNet 데이터 다운로드**: 로컬 경로(`/data/user3/git_repo/human`) 우선 사용 및 자동 다운로드 지원.

## 4. 사용자 가이드 및 주의사항 (User Guide & Warnings)

### 기본 사용법

```r
source("myR/R/cci_nichenet_run.R")
library(qs)

# 1. 데이터 로드
sobj <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# 2. DEG 리스트 준비 (예시)
deg_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = "CD4+ T-cells",
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01)
)

# 3. 분석 실행
results <- run_cci_analysis(
  sobj = sobj,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "CD4+ T-cells",
  sender_clusters = c("Monocytes", "NK Cells"), # NULL이면 자동 식별
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  output_dir = "/data/user3/sobj/cci_output"
)
```

### Critical Warnings (주의사항)
1.  **메모리 사용량**: NicheNet 분석은 메모리를 많이 사용하므로, 대규모 데이터셋 분석 시 주의하십시오. 다운샘플링 데이터로 먼저 테스트하는 것을 권장합니다.
2.  **DEG 매칭**: `deg_df`의 `cluster` 컬럼 값은 `receiver_cluster` 파라미터와 정확히 일치해야 합니다.
3.  **데이터 경로**: NicheNet 데이터가 없는 경우 최초 실행 시 다운로드를 시도하며 시간이 소요됩니다. 사내 서버 경로(`/data/user3/git_repo/human`)를 활용하세요.

## 5. 방법론 (Methodology)

### 분석 로직
1.  **입력 검증**: Seurat 객체와 DEG 리스트의 정합성 확인.
2.  **Receiver DEG 추출**: 입력된 `deg_df`에서 Receiver 클러스터에 해당하는 DEG 추출.
3.  **Sender 식별**: 지정된 Sender 또는 전체 클러스터에서 Sender 후보 식별.
4.  **발현 필터링**: Sender와 Receiver에서 일정 비율(`min_pct_expressed`) 이상 발현되는 유전자만 선별.
5.  **NicheNet 분석**:
    *   Ligand-Target Potential 점수 계산.
    *   Ligand Activity 예측 및 우선순위화.
    *   Top Ligand의 Target Gene 네트워크 추론.
6.  **결과 저장 및 시각화**: 결과 객체(.qs) 및 Plot 생성.

## 6. 부록 (Appendix)

### 출력 파일 구조
*   `nichenet_results.qs`: 전체 분석 결과 객체.
*   `NicheNet_Ligand_Target_Heatmap.png`: 주요 Ligand와 Target 유전자 간의 관계.
*   `NicheNet_Circos_LR.pdf`: Ligand-Receptor 상호작용 Circos plot.
*   `analysis_summary.qs`: 분석 메타데이터 요약.

### 관련 스크립트
*   `myR/R/cci_nichenet_run.R`: 메인 분석 함수 (구 `run_cci_analysis.R`).
*   `myR/R/cci_nichenet_wrapper.R`: NicheNet 래퍼 (구 `CCI.R`)
*   `scripts/cci/test_is5_downsample.R`: 다운샘플링 데이터 테스트 스크립트.
