# Comparison & Consensus Plan

## 1. Comparison Metrics
*   **Correlation**: LogFC correlation between methods (Heatmap).
*   **Overlap**: Jaccard Index / UpSet Plot of DEGs (FDR < 0.05).
*   **Rank Consistency**: Rank-biased Overlap (RBO).

## 2. Consensus Logic
*   **Current**: Majority Vote / Average Rank.
*   **Future**: Precision-weighted Averaging (based on p-value/SE).

## 3. Implementation
*   Target: `aggregate_cluster_deg_consensus.R`.
*   Timeline: After method stabilization (Phase 2 completion).
