# FGS & pt.umap Integration Guide

## Overview
This guide explains how to use gene signatures derived from FGS/TML to perform patient-level dimensionality reduction using `pt.umap`.

## Workflow

1.  **Generate Signatures**: Use FGS/TML to identify key gene signatures and their importance (CMGI).
2.  **Compute Patient Scores**: Calculate signature scores for each patient.
3.  **Run pt.umap**: Use the patient score matrix as input for UMAP.

## Step-by-Step

### 1. Select Signatures
You can use:
*   **All L1 Signatures**: Use all signatures found by FGS.
*   **Top Signatures**: Select signatures with high importance in TML (CMGI).
*   **Top Genes**: Use the top genes identified by CMGI to form a new "Meta Signature".

### 2. Compute Patient Scores
Instead of cell-type frequency, we use the aggregated expression of signatures per patient.

```r
# Assuming 'sobj' is your Seurat object
# 'fgsa' is FGS result, 'tmla' is TML result

# Option A: Use L1 Signatures directly
sigs <- fgsa

# Option B: Use Top Genes from CMGI
cmgi <- compute_meta_gene_importance(tmla)
top_genes <- head(cmgi$gene_summary$gene, 50)
# Create a list for AddModuleScore
sigs <- list(MetaSig = top_genes)

# Calculate scores per cell
# Seurat's AddModuleScore adds columns to meta.data
sobj <- Seurat::AddModuleScore(sobj, features = sigs, name = "MetaSig_")
score_col <- paste0("MetaSig_", 1) # Seurat appends index

# Aggregate per patient (e.g., mean score)
# Assuming 'patient_id' is the column for patients
meta <- sobj@meta.data
pt_scores <- stats::aggregate(meta[[score_col]], by = list(Patient = meta$patient_id), FUN = mean)
rownames(pt_scores) <- pt_scores$Patient
pt_matrix <- as.matrix(pt_scores[, 2, drop = FALSE])
colnames(pt_matrix) <- "MetaSig"
```

### 3. Run pt.umap
The `pt_matrix` (Patients x Signatures) replaces the Frequency Matrix.

```r
# Run UMAP on the score matrix
# You might need to adapt pt.umap input format or use standard UMAP
library(umap)
umap_res <- umap(pt_matrix)

# Visualize
plot_df <- data.frame(
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2],
  Patient = rownames(pt_matrix)
)

ggplot(plot_df, aes(x = UMAP1, y = UMAP2, label = Patient)) +
  geom_point() +
  geom_text(vjust = 1.5) +
  theme_minimal() +
  labs(title = "Patient UMAP based on FGS/TML Signatures")
```

## Why this works
*   **Biological Relevance**: Instead of just cell abundance, we use functional states (signatures) to group patients.
*   **Resolution**: Captures subtle differences in gene expression programs that frequency analysis might miss.
