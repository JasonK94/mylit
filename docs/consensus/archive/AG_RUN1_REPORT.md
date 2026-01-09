# AG_run1 Analysis Report

## Overview
- **Dataset**: `/data/user3/sobj/is2_IS_3_clustered.qs`
- **Levels**:
  - `Broad_anno3big` (5-7 clusters)
  - `Fine_anno3` (22-25 clusters)
- **Methods**: 12 methods including standard PB, Dream, Nebula, and Muscat variants.
- **Output Directory**: `/data/user3/sobj/consensus/AG_run1`

## Implementation Details
1. **Methods**:
   - Standalone: `edgeR-LRT`, `edgeR-QLF`, `DESeq2-Wald`, `DESeq2-LRT`, `limma-voom`, `limma-trend`, `dream`.
   - Single-Cell Mixed Models: `nebula` (NBLMM).
   - Wrapper: `muscat` (edgeR/DESeq2/limma variants).
   
2. **Consensus & Similarity**:
   - **Feature**: Test Statistics (`statistic` column: t, F, Z).
   - **Similarity Metric**: Spearman correlation of statistic vectors (pairwise complete).
   - **Clustering**: Hierarchical clustering (Ward.D2) of the correlation matrix.
   - **Outputs**:
     - `AG_run1_method_similarity_heatmap.png`: Heatmap of method correlations.
     - `AG_run1_method_clustering_dendrogram.png`: Dendrogram of methods.

3. **Covariates**:
   - `sex`, `age` included in design matrices where supported.

## Script
- `scripts/consensus/run_AG_run1.R`
