# Cell-to-Cell Interaction (CCI) Analysis Integrated Guide

This document is the integrated guide for the Cell-to-Cell Interaction (CCI) analysis module. It analyzes ligand-receptor interactions based on NicheNet and supports flexible analysis by utilizing pre-computed DEG lists.

## 1. Introduction

### Purpose
Analyze cell-to-cell interactions in single-cell RNA sequencing (scRNAseq) data to identify **Sender cells** and **Ligands** that induce gene expression changes (DEG) in **Receiver cells** under specific conditions (e.g., disease vs. normal).

### Key Features
1.  **Direct DEG List Input**: Accepts DEG lists computed by various methods such as `FindMarkers`, `limma`, `edgeR`, etc.
2.  **Automatic Sender Identification**: Automatically sets all clusters except Receiver as Senders, or allows user specification.
3.  **NicheNet Integration**: Constructs Ligand-Receptor-Target gene networks using the NicheNet database.
4.  **Publication-Quality Visualization**: Automatically generates various visualization outputs including Circos plots and heatmaps.

### Key Terms
*   **Sender**: Cell type that sends signals (Ligands).
*   **Receiver**: Cell type that receives signals (Ligands) and responds. Cells where expression changes (DEG) of Target Genes are observed.
*   **Ligand**: Substance secreted by Sender cells that binds to Receptors on Receiver cells.
*   **Receptor**: Protein present on the surface of Receiver cells that binds to Ligands.
*   **Target Gene**: Downstream genes (DEG) whose expression is regulated by Ligand-Receptor interactions.

## 2. Workflow Visualization (시각화)

```mermaid
flowchart TD
    Input[Input Data: Seurat Object, DEG List] --> Prep[Data Preparation]
    Prep --> Sender[Identify Sender Clusters]
    Prep --> Filter[Filter Expressed Genes]
    
    Sender --> NicheNet[Run NicheNet Analysis]
    Filter --> NicheNet
    
    NicheNet --> LigandAct[Predict Ligand Activities]
    LigandAct --> Prioritize[Prioritize Ligands]
    Prioritize --> Network[Infer Ligand-Target Network]
    
    Network --> Vis[Visualization]
    Vis --> Heatmap[Heatmaps: Ligand-Target, Ligand-Receptor]
    Vis --> Circos[Circos Plot]
    
    Network --> Output[Save Results (.qs)]
```

## 3. Development Log & Improvements

### Major Changes
*   **v1.0 (2025-11-14)**: Initial implementation. Direct DEG list input support, modular structure (preparation, analysis, saving) established.
*   **Optimizations**:
    *   `receiver_de_table` reuse: Prevents repeated DEG calculations for the same Receiver.
    *   Memory management: Cleans up large objects after saving intermediate results (`rm`, `gc`).
    *   Column mapping: Supports various DEG table formats (`avg_log2FC` vs `logFC`, etc.).

### Resolved Issues
*   **Expressed Gene Filtering Error**: Fixed issue where 0 genes were returned due to missing `Idents(sobj)` setting.
*   **DEG Filtering Problem**: Relaxed cutoff (`1.1`) to handle cases where adjusted p-value exceeds 1.0.
*   **NicheNet Data Download**: Prioritizes local path (`/data/user3/git_repo/human`) and supports automatic download.

## 4. User Guide & Warnings

### Basic Usage

```r
source("myR/R/cci/run_cci_analysis.R")
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

### 주요 파라미터
*   `sobj`: Seurat 객체.
*   `deg_df`: DEG 데이터프레임. 필수 컬럼: `gene`, `cluster`, `avg_log2FC`, `p_val_adj`.
*   `receiver_cluster`: 분석 대상 Receiver 클러스터 이름.
*   `nichenet_data_dir`: NicheNet 데이터 경로 (기본값: `/data/user3/git_repo/human`).

### Critical Warnings
1.  **Memory Usage**: NicheNet analysis uses significant memory. Exercise caution when analyzing large datasets. It is recommended to test with downsampled data first.
2.  **DEG Matching**: The `cluster` column value in `deg_df` must exactly match the `receiver_cluster` parameter.
3.  **Data Path**: If NicheNet data is not available, the first run will attempt to download it, which may take time. Utilize the internal server path (`/data/user3/git_repo/human`).

## 5. Methodology

### Analysis Logic
1.  **Input Validation**: Verify consistency between Seurat object and DEG list.
2.  **Receiver DEG Extraction**: Extract DEGs corresponding to the Receiver cluster from the input `deg_df`.
3.  **Sender Identification**: Identify Sender candidates from specified Senders or all clusters.
4.  **Expression Filtering**: Select only genes expressed above a certain percentage (`min_pct_expressed`) in both Sender and Receiver.
5.  **NicheNet Analysis**:
    *   Calculate Ligand-Target Potential scores.
    *   Predict and prioritize Ligand Activities.
    *   Infer Target Gene networks for top Ligands.
6.  **Result Saving and Visualization**: Generate result objects (.qs) and plots.

## 6. Appendix

### Output File Structure
*   `nichenet_results.qs`: Complete analysis result object.
*   `NicheNet_Ligand_Target_Heatmap.png`: Relationships between major Ligands and Target genes.
*   `NicheNet_Circos_LR.pdf`: Ligand-Receptor interaction Circos plot.
*   `analysis_summary.qs`: Analysis metadata summary.

### Related Scripts
*   `myR/R/cci/run_cci_analysis.R`: Main analysis function.
*   `scripts/cci/test_is5_downsample.R`: Downsampled data test script.

