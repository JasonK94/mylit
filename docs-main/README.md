# myR: Single-Cell RNA-Seq Analysis Package

## Overview
`myR` is an integrated analysis package for single-cell RNA sequencing (scRNAseq) data. It provides various analysis modules including differential expression analysis (DEG), cell-cell interaction (CCI), differential abundance analysis, patient-level analysis, and trajectory analysis.

## 📚 Integrated Guide

**For a complete module overview and workflow, please refer to:**

- **English**: [`docs/INTEGRATED_GUIDE.md`](docs/INTEGRATED_GUIDE.md)
- **Korean**: [`docs/INTEGRATED_GUIDE_KR.md`](docs/INTEGRATED_GUIDE_KR.md)

The integrated guide includes:
- Complete analysis pipeline visualization (Mermaid diagrams)
- Role, input/output, and key methodologies for each module
- Common analysis workflow examples
- Links to detailed module documentation

## Main Modules

| Module | Purpose | Documentation |
|--------|---------|---------------|
| **analysis** | Mixed-Effects Model DEG (NEBULA) | `docs/analysis/` |
| **deg-consensus** | Multi-model DEG Consensus | `docs/deg-consensus-dev/` |
| **lds** | Limma-Dream-SVA | `docs/lds/` |
| **milo** | Differential Abundance | `docs/milo/` |
| **cci** | Cell-Cell Interaction (NicheNet) | `docs/cci/` |
| **fgs** | Gene Signature Discovery | `docs/fgs/` |
| **pt.umap** | Patient-Level Analysis | `docs/pt.umap/` |
| **pseudotime** | Trajectory Inference | `docs/pseudotime-dev/` |

For detailed usage of each module, please refer to the "Module Documentation" section in the integrated guide.

---

## CCI (Cell-to-Cell Interaction) Analysis Tool

### Overview
An integrated tool for performing Cell-to-Cell Interaction analysis on scRNAseq data. It focuses on NicheNet to analyze ligand-receptor interactions and directly accepts DEG lists to identify sender cell types that contribute to changes in receiver cell types.

## Quick Start

### 1. Environment Setup
```r
# Load package
devtools::load_all("path/to/myR")

# Or load function source directly
source("path/to/myR/R/cci/run_cci_analysis.R")

# CCI.R containing run_nichenet_analysis - prioritize worktree file
cci_core_worktree <- "path/to/worktree/myR/R/CCI.R"
cci_core_mainrepo <- "path/to/mainrepo/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  stop("CCI.R not found in worktree or main repository.")
}
```

### 2. Basic Usage
```r
# Load data
library(qs)
sobj <- qs::qread("path/to/your/data.qs")

# Prepare DEG list (example)
deg_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = c("Cluster1", "Cluster1", "Cluster1"),
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01)
)

# Run CCI analysis
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

## Project Structure

```
worktree/cci/
  docs/
    cci/
      cci.md                    # Main feature documentation
      TEST_INSTRUCTIONS.md      # Test guide
      devlog.md                 # Development log
      CCI_DEVELOPMENT_PLAN.md   # Development plan
  scripts/
    cci/
      test_cci.R                # Test script
  myR/
    R/
      cci/
        run_cci_analysis.R      # Main CCI analysis function
        prepare_cci_data.R      # Data preparation and validation
        nichenet_wrapper.R      # NicheNet analysis wrapper
        save_cci_results.R      # Result saving utility
        utils_cci.R             # CCI utility functions
```

> According to `guide.md`: Scripts used in the CCI worktree should be added only to `scripts/cci/`, and documentation only to `docs/cci/`.

## Main Functions

### `run_cci_analysis()`
The main function for CCI analysis. It accepts a Seurat object, cluster information, and DEG list, and performs NicheNet analysis.

**Required Parameters**:
- `sobj`: Seurat object
- `cluster_col`: Cluster column name
- `deg_df`: DEG dataframe
- `receiver_cluster`: Receiver cluster ID
- `condition_col`: Condition column name
- `condition_oi`: Condition of interest

**Optional Parameters**:
- `sender_clusters`: Sender cluster vector (NULL for automatic identification)
- `species`: "human" or "mouse" (default: "human")
- `auto_save`: Auto-save option (default: TRUE)

## Data Requirements

### DEG Dataframe Format
Required columns:
- `gene`: Gene name
- `cluster`: Cluster ID
- `avg_log2FC` or `logFC`: Log fold change
- `p_val_adj` or `FDR`: Adjusted p-value

### Seurat Object Metadata
Required columns:
- Cluster information column (e.g., `anno3.scvi`)
- Condition information column (e.g., `g3`)

### Receiver DEG Reuse
- `run_nichenet_analysis()` can directly accept `receiver_de_table` to run NicheNet without re-executing `FindMarkers()`.
- If column names differ, map them using `receiver_gene_col`, `receiver_logfc_col`, `receiver_pval_col`.
- `run_cci_analysis()` extracts receiver DEG internally and passes it in the same way, so no duplicate computation occurs even when running the same receiver repeatedly.

## Output

Results are returned as a list containing:
- `nichenet_results`: NicheNet analysis results
- `sender_receiver_map`: Sender-receiver mapping
- `deg_summary`: DEG summary information
- `saved_path`: Saved file path (when auto_save=TRUE)

## Testing

Run test script:
```r
source("path/to/test_cci.R")
```

## Documentation

- **Detailed Documentation**: `docs/cci/cci.md`
- **Module-Specific Guide**: `docs/cci/CCI_module.md`
- **Test Guide**: `docs/cci/TEST_INSTRUCTIONS.md`
- **Development Log**: `docs/cci/devlog.md`

## References

- NicheNet Official Documentation: https://github.com/saeyslab/nichenetr
- CCI.R function (worktree priority): `path/to/worktree/myR/R/CCI.R` (or `path/to/mainrepo/myR/R/CCI.R` if not found)
