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
