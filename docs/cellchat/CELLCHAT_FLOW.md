# CellChat Data Transformation Flow

This document outlines the step-by-step transformation of data within the CellChat pipeline, focusing on input/output logic rather than implementation details.

## 1. Initialization
**Input**: Seurat Object (`sobj`)
**Function**: `createCellChat(object = sobj, group.by = "...")`
**Transformation**:
- Extracts Norm Data -> `cc@data` (LogNormalized matrix)
- Extracts Metadata -> `cc@meta`
- Extracts Identities -> `cc@idents`

## 2. Ligand-Receptor Identification (Preprocessing)
### A. Identify Over-expressed Genes
**Function**: `identifyOverExpressedGenes(cc)`
**Logic**: 
- Validates "expressed" genes (e.g., detected in >10% of cells in a group).
- Performs **Wilcoxon Rank Sum Test** (One-vs-Rest) for each cell group.
- P-value threshold ($P < 0.05$).
**Output**: `cc@var.features$features` (List of significant DEGs per group).

### B. Identify Over-expressed Interactions
**Function**: `identifyOverExpressedInteractions(cc)`
**Logic**:
- Mapping: Matches identified DEGs (from A) to Ligand-Receptor pairs in `cc@DB`.
- Constraint: Both Ligand AND Receptor must be DEGs (or "expressed" depending on parameters).
**Output**: `cc@LR$LRsig` (Subset of DB containing only valid LR pairs for this dataset).
- *This is the candidate list for probability calculation.*

## 3. Probability Calculation (Core Logic)
### A. Compute Communication Probability
**Function**: `computeCommunProb(cc, type = "triMean")`
**Logic**:
1. **Expression Value**: Uses *Trimean* (Q1, Q2, Q3 weighted avg) of expression per cell group.
   - $AvgExp = (Q1 + 2 \cdot Q2 + Q3) / 4$
2. **Law of Mass Action**:
   - $Prob_{i,j} = \frac{L_i \cdot R_j}{K_h + L_i \cdot R_j}$ (Hill function)
   - Where $L_i$ = Ligand exp in Group i, $R_j$ = Receptor exp in Group j.
3. **Significance Test (Permutation)**:
   - Randomly permute cell labels (Sender/Receiver identities).
   - Re-calculate Prob many times (e.g., 100).
   - Count how many times random Prob > observed Prob -> P-value.
   - Set Prob to 0 if $P > 0.05$.
**Output**: `cc@net$prob` (3D Array: `Sources` x `Targets` x `LR_Pairs`).

### B. Filter Communication
**Function**: `filterCommunication(cc, min.cells = 10)`
**Logic**:
- Set interactions to 0 if the cell group has fewer than `min.cells`.
**Output**: Updates `cc@net$prob`.

## 4. Aggregation (Post-Processing)
### A. Aggregate Network
**Function**: `aggregateNet(cc)`
**Logic**:
- Sums/Averages probabilities across all LR pairs to get Group-to-Group strength.
**Output**: 
- `cc@net$count`: Number of significant interactions ($P < 0.05$).
- `cc@net$weight`: Sum of interaction probabilities.

### B. Merge CellChats (For Comparison)
**Function**: `mergeCellChat(list(cc1, cc2))`
**Logic**:
- Merges `net` slots from multiple objects into a list or array.
- Enables differential analysis (e.g., `cc2 - cc1`).

---
**Summary**:
`Seurat` -> `DEGs` -> `Valid LRs` -> `Probabilities (Hill + Permutation)` -> `Significant Network` -> `Visualizations`

---

## 5. Metric Interpretation Guide

### A. `identifyCommunicationPatterns` Heatmap: **Contribution**
- **Cell Patterns (행)**: 각 cell type이 해당 pattern에 얼마나 기여하는지 (NMF W matrix)
- **Communication Patterns (열)**: 각 signaling pathway가 해당 pattern에 기여하는 정도 (NMF H matrix)
- **값 범위**: 0~1, 높을수록 해당 pattern의 핵심 요소

### B. `netAnalysis_river` (Sankey Diagram)
- Ribbon **두께** = Cell type ↔ Pattern ↔ Signaling 간 연결 강도 (contribution)
- 두꺼울수록 해당 조합이 특정 pattern을 특징짓는 핵심 요소

### C. `netAnalysis_dot` Contribution
- Dot **크기** = contribution 값
- Dot **색상** = contribution 강도 (연속 스케일)

### D. `netVisual_diffInteraction` 화살표
- **빨간색**: 두 번째 그룹 (Stroke) 방향으로 **증가**
- **파란색**: 첫 번째 그룹 (Control) 방향으로 **증가** (= Stroke에서 감소)
- **두께**: 차이의 절대값 (count 또는 weight)
- **방향**: Sender → Receiver

### E. `rankNet` RIF 및 색상
- **RIF (Relative Information Flow)**: 각 signaling pathway의 총 통신 확률 합
- **그룹별 색** (파랑/빨강): 해당 그룹에서 **유의하게 활성화**된 pathway
- **검정색**: 두 그룹 간 통계적 차이 없음 (`do.stat = TRUE` 시 Wilcoxon test)

### F. `netVisual_bubble` Communication Prob
- **점 크기**: Communication probability (Hill function 결과)
- **X축**: Sender|Receiver cell pair
- **Y축**: L-R pair (pathway 내)
- **`thresh`**: P-value 필터 (단, merged object에서 무효할 수 있음)

### G. `netVisual_heatmap` Relative Values
- **Relative value**: (Group2 - Group1) / max(|diff|) 정규화
- **음수**: Control 방향 증가 (Stroke에서 감소)
- **양수**: Stroke 방향 증가
- **축 Bar graph**: 누적 incoming/outgoing strength (음수 가능)

### H. `netVisual_aggregate` 주의사항
- **Merged object에서 사용 불가**: `mergeCellChat` 결과는 `@netP$Control`, `@netP$Stroke` 구조
- **해결책**: 개별 object (`cc1` 또는 `cc2`)에서 호출
```r
netVisual_aggregate(cc2, signaling = "MIF", layout = "circle")  # Stroke만
```



```mermaid
%%{init: {'theme':'base', 'themeVariables': {'primaryEdgeColor':'#333333', 'primaryEdgeThickness':3, 'primaryTextColor':'#000000', 'primaryBorderColor':'#666666', 'edgeLabelBackground':'#ffffff', 'tertiaryColor':'#999999'}}}%%
flowchart TD

  classDef step fill:#f3f0ff,stroke:#333,stroke-width:2px;
  classDef logic fill:#e1f5fe,stroke:#0277bd,stroke-width:2px,stroke-dasharray: 5 5;
  classDef io fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px;
  classDef warn fill:#fff3e0,stroke:#ef6c00,stroke-width:2px;
  classDef ext fill:#ffebee,stroke:#c62828,stroke-width:2px;

  Start(["Seurat Object / Cell labels / meta: sample, condition"]):::io
  Start --> Prep0["Pre-checks: celltype labels fixed; sample_id/condition present; min.cells policy 결정"]:::step
  Prep0 --> Choose{"Strategy?"}:::step

  %% -------------------------
  %% Strategy 1: Pooled
  %% -------------------------
  Choose -- "Method 1: Pooled by condition<br/>High power, descriptive" --> PoolSplit["Split object by condition"]:::step

  subgraph Pooled["Method 1: Pooled CellChat per condition"]
    direction TB

    PoolSplit --> CondA["Condition A pooled cells"]:::io
    PoolSplit --> CondB["Condition B pooled cells"]:::io

    CondA --> A0["Create CellChat object"]:::step
    CondB --> B0["Create CellChat object"]:::step

    A0 --> A1["subsetData: signaling genes by DB"]:::step
    B0 --> B1["subsetData: signaling genes by DB"]:::step

    A1 --> A2["identifyOverExpressedGenes<br/>cell-group specificity; Wilcoxon etc"]:::logic
    B1 --> B2["identifyOverExpressedGenes<br/>cell-group specificity; Wilcoxon etc"]:::logic

    A2 --> A3["identifyOverExpressedInteractions<br/>prune DB pairs using OE genes"]:::logic
    B2 --> B3["identifyOverExpressedInteractions<br/>prune DB pairs using OE genes"]:::logic

    A3 --> A4{"Optional? smoothData / projectData<br/>PPI-based smoothing"}:::warn
    B3 --> B4{"Optional? smoothData / projectData<br/>PPI-based smoothing"}:::warn

    A4 --> A5["computeCommunProb<br/>Hill/mass-action; robust mean; permutation pval"]:::logic
    B4 --> B5["computeCommunProb<br/>Hill/mass-action; robust mean; permutation pval"]:::logic

    A5 --> A6["filterCommunication<br/>(min.cells 등 희소군 제거)"]:::step
    B5 --> B6["filterCommunication<br/>(min.cells 등 희소군 제거)"]:::step

    A6 --> A7["computeCommunProbPathway<br/>pathway sum"]:::step
    B6 --> B7["computeCommunProbPathway<br/>pathway sum"]:::step

    A7 --> A8["aggregateNet<br/>count / weight matrices"]:::step
    B7 --> B8["aggregateNet<br/>count / weight matrices"]:::step
  end

  %% Merge & compare (descriptive)
  A8 --> MergePool["mergeCellChat A,B"]:::logic
  B8 --> MergePool

  MergePool --> ComparePool["compareInteractions<br/>count diff, weight diff"]:::logic
  ComparePool --> VisDiff1["netVisual_diffInteraction<br/>edge Δcount/Δweight"]:::step
  ComparePool --> VisHeat1["netVisual_heatmap<br/>matrix diff"]:::step
  ComparePool --> Rank1["rankNet / information flow<br/>pathway-level shifts"]:::step

  MergePool --> Curation1{"Optional: DE overlay / curation"}:::warn
  Curation1 --> MapDE1["netMappingDEG<br/>condition DE -> edges annotation"]:::logic
  MapDE1 --> Subset1["subsetCommunication<br/>filter by ligand/receptor logFC etc"]:::step
  Subset1 --> Bubble1["netVisual_bubble / chord<br/>focused LR/pathways"]:::step

  MergePool --> Patterns1{"Optional: patterns / manifold"}:::warn
  Patterns1 --> NMF1["identifyCommunicationPatterns<br/>outgoing/incoming summary -> NMF"]:::logic
  NMF1 --> River1["netAnalysis_river<br/>pattern↔cell↔pathway alluvial"]:::step
  NMF1 --> Dot1["netAnalysis_dot<br/>loadings dot plot"]:::step

  Patterns1 --> Sim1["computeNetSimilarity<br/>pathway networks similarity"]:::logic
  Sim1 --> Emb1["netEmbedding<br/>UMAP/etc on similarity"]:::logic
  Emb1 --> Clus1["netClustering<br/>community / KNN graph"]:::logic
  Clus1 --> VisEmb1["Embedding scatter + clusters<br/>visualize conserved vs specific signaling"]:::step

  %% -------------------------
  %% Strategy 2: Sample-wise
  %% -------------------------
  Choose -- "Method 2: Sample-wise CellChat<br/>Rigor via replicate variance" --> SplitS["Split object by sample_id"]:::step

  subgraph SampleWise["Method 2: Per-sample CellChat replicate-aware"]
    direction TB

    SplitS --> S1["For each sample s: create CellChat_s"]:::step
    S1 --> S2["Run same CC core pipeline per sample:<br/>subsetData -> OE genes/interactions -> computeCommunProb -> filterCommunication -> aggregateNet"]:::logic
    S2 --> S3["Output per-sample matrices:<br/>weight_s i,j, count_s i,j<br/>and/or pathway matrices"]:::io

    %% Two branches after per-sample outputs
    S3 --> BranchS{"Goal after per-sample run?"}:::step

    BranchS -- "A) Descriptive per-condition visualization" --> MergeCond["Group samples by condition<br/>mergeCellChat list"]:::logic
    MergeCond --> CompareDesc["compareInteractions + plots<br/>descriptive only"]:::warn

    BranchS -- "B) Inferential stats across samples<br/>recommended for rigor claims" --> StatsExt["External inference step<br/>edge-wise stats"]:::ext
    StatsExt --> Model1["Example:<br/>weight_ij,s ~ condition + 1|sample<br/>or Wilcoxon/permutation on sample summaries"]:::ext
    Model1 --> FDR1["Multiple testing control<br/>FDR on edges/pathways"]:::ext
    FDR1 --> SigEdges["Significant differential edges/pathways<br/>replicate-aware"]:::io
    SigEdges --> VisRig["Visualize significant set<br/>(overlay on CC networks / heatmap / bubble)"]:::step
  end

  %% Connect sample-wise descriptive into common visualization lane
  CompareDesc --> VisDiff2["netVisual_diffInteraction<br/>descriptive"]:::step
  CompareDesc --> Rank2["rankNet / info flow<br/>descriptive"]:::step

  %% Styling for main nodes
  class Start,Prep0,Choose step
```
