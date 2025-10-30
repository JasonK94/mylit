# myR Package Refactoring Plan

## Executive Summary

This document outlines a comprehensive reorganization and refactoring plan for the `myR` scRNA-seq analysis package. The current codebase has grown organically and would benefit from improved organization, reduced redundancy, and better separation of concerns.

## Current State Analysis

### Existing Files (17 total)
1. **CCI.R** - Cell-cell interaction analysis (NicheNet)
2. **cluster_frequency.R** - Cluster composition and statistical analysis
3. **demulti_utils.R** - Demultiplexing utilities
4. **documents.R** - Documentation helpers
5. **gene_list.R** - Gene signature collections
6. **GeoMx.R** - GeoMx spatial profiling analysis
7. **legacy.R** - Deprecated/legacy functions (2500+ lines!)
8. **markers.R** - Marker gene analysis utilities
9. **pathway_analysis.R** - GO/GSEA pathway analysis
10. **plots.R** - Visualization functions (2400+ lines!)
11. **pseudobulk_deg.R** - Pseudobulk differential expression (very large file!)
12. **pseudotime.R** - Trajectory inference (Slingshot, tradeSeq, Monocle3)
13. **seurat_to_scanpy.R** - Seurat to Python/AnnData conversion
14. **signature.R** - Module scoring and signature analysis
15. **test_claude.R** - Test file (should be removed)
16. **test.R** - Test file (should be removed)
17. **utils.R** - General utility functions

### Key Issues Identified

#### 1. **Organizational Issues**
- **plots.R is too large** (2400+ lines) with mixed responsibilities:
  - Basic plots (histograms, density, line plots)
  - Complex statistical plots (boxplots with tests)
  - Specialized plots (CDF, cumulative line, scatter-smooth)
  - Heatmaps (gene sets, individual genes)
  - Composition plots (bar graphs)
  
- **legacy.R contains outdated functions** that should either be:
  - Updated and integrated
  - Removed if truly obsolete
  - Documented as deprecated

- **pseudobulk_deg.R is extremely large** and should be split

- **Inconsistent naming conventions**:
  - Some functions use `my` prefix (myhm_genesets, mybox, mydensity)
  - Others don't follow a clear convention
  - Mix of snake_case and camelCase

#### 2. **Redundancy Issues**
- **Multiple heatmap functions** with overlapping functionality:
  - `myhm_genesets2()`, `myhm_genesets3()`, `myhm_genesets4()` in legacy.R and plots.R
  - `myhm_genes2()`, `myhm_genes4()` in legacy.R and plots.R
  - Should be consolidated into a single, flexible implementation

- **Module scoring methods duplicated**:
  - Simple averaging in `myhm_genesets*()`
  - AddModuleScore wrapper in signature.R
  - Need unified approach

- **Similar plotting functions**:
  - Multiple scatter plot functions with different features
  - Box plot variations (mybox for Seurat, mybox_df for data frames, mybox_geomx for GeoMx)

- **Statistical testing duplicated**:
  - cluster_frequency.R has post-hoc analysis
  - plots.R has similar functionality in mybox functions

#### 3. **Efficiency Issues**
- **No caching mechanisms** for expensive operations
- **Limited parallelization** in some functions
- **Redundant data transformations**
- **No input validation helpers** (repeated validation code)

## Proposed New Structure

### Core Principles
1. **Separation of Concerns**: Group related functions together
2. **Thin Files**: No single file > 800 lines
3. **Clear Hierarchy**: Core → Analysis → Visualization → Utilities
4. **Consistent Naming**: Adopt clear naming conventions
5. **Reduce Redundancy**: Consolidate similar functions

### New File Organization

```
R/
├── core/
│   ├── data_loading.R           # Seurat object I/O, conversion
│   ├── preprocessing.R          # QC, normalization helpers
│   ├── validation.R             # Input validation utilities
│   └── conversion.R             # Format conversions (Seurat ↔ SCE ↔ AnnData)
│
├── analysis/
│   ├── markers/
│   │   ├── finding.R            # FindMarkers wrappers
│   │   ├── filtering.R          # marker_trim, marker_filter
│   │   └── synthesis.R          # Meta-analysis functions
│   │
│   ├── differential_expression/
│   │   ├── pseudobulk_core.R   # Core pseudobulk functions
│   │   ├── pseudobulk_edger.R  # edgeR-based analysis
│   │   ├── pseudobulk_limma.R  # limma-based analysis (if needed)
│   │   └── linear_models.R     # linear_seurat, regression analysis
│   │
│   ├── cell_communication/
│   │   ├── nichenet.R           # NicheNet analysis
│   │   └── circos.R             # Circos plot utilities
│   │
│   ├── spatial/
│   │   └── geomx.R              # GeoMx analysis pipeline
│   │
│   ├── trajectory/
│   │   ├── slingshot.R          # Slingshot wrapper
│   │   ├── monocle.R            # Monocle3 analysis
│   │   ├── tradeseq.R           # tradeSeq analysis
│   │   └── dynamics.R           # Gene dynamics (GAM)
│   │
│   ├── pathway/
│   │   ├── go_analysis.R        # GO enrichment
│   │   ├── gsea.R               # GSEA analysis
│   │   └── pathway_utils.R      # Conversion, formatting
│   │
│   └── composition/
│       ├── cluster_frequency.R  # Statistical tests for composition
│       ├── cluster_stats.R      # Summary statistics
│       └── confounders.R        # Confounder adjustment
│
├── visualization/
│   ├── basic_plots.R            # Basic univariate plots (histogram, density, line)
│   ├── scatter_plots.R          # Scatter, smooth, colored variants
│   ├── distribution_plots.R     # CDF, violin, box plots
│   ├── composition_plots.R      # Bar graphs (cmb, acmb, cml)
│   ├── heatmaps.R               # Unified heatmap functions
│   ├── feature_plots.R          # Gene/module expression plots
│   ├── statistical_plots.R      # Plots with statistical tests
│   └── plot_themes.R            # Common themes and styling
│
├── signatures/
│   ├── gene_lists.R             # Curated gene signature collections
│   ├── module_scoring.R         # Module score calculation
│   ├── enrichit.R               # enrichIt wrapper
│   └── signature_viz.R          # Signature visualization
│
├── utilities/
│   ├── demultiplexing.R         # Demux utilities
│   ├── general_utils.R          # General helper functions
│   ├── sample_utils.R           # Sample sorting, naming
│   └── plot_utils.R             # Plot saving, conflict resolution
│
└── deprecated/
    └── legacy_functions.R       # Clearly marked deprecated functions
```

## Specific Refactoring Tasks

### Phase 1: Foundation (Core & Utilities)

#### Task 1.1: Create Core Module
- **New file**: `R/core/validation.R`
  - Extract common validation patterns
  - Create reusable validation functions
  - Add informative error messages

- **New file**: `R/core/conversion.R`
  - Move `save_seurat_to_h5ad()` here
  - Add reverse conversion if needed
  - Document format requirements

- **New file**: `R/core/data_loading.R`
  - Functions for loading various formats
  - Seurat object initialization helpers

#### Task 1.2: Create Utilities Module
- **New file**: `R/utilities/general_utils.R`
  - Move `downsample_sobj()`, `printmy()`, `printMy()`
  - Add new utility functions as needed

- **New file**: `R/utilities/sample_utils.R`
  - Move `sort_samples()` (deduplicate - it's defined twice!)
  - Add sample naming utilities

- **New file**: `R/utilities/plot_utils.R`
  - Move `save_plot_with_conflict_resolution()`
  - Add common plot formatting utilities

### Phase 2: Visualization Consolidation

#### Task 2.1: Consolidate Heatmap Functions
**Target**: `R/visualization/heatmaps.R`

Create a unified heatmap interface:
```r
plot_heatmap <- function(
    data,                    # Seurat object or matrix
    features,                # Genes or gene sets
    feature_type = c("genes", "gene_sets"),
    group.by = "seurat_clusters",
    value = c("average", "sum"),
    assay = "SCT",
    scale_method = c("feature", "group", "none"),
    show_stats = FALSE,
    ...
) {
    # Unified implementation
}
```

**Consolidate**:
- `myhm_genesets2()`, `myhm_genesets3()`, `myhm_genesets4()` → `plot_heatmap(..., feature_type="gene_sets")`
- `myhm_genes2()`, `myhm_genes4()` → `plot_heatmap(..., feature_type="genes")`

#### Task 2.2: Reorganize Plot Types
- **`basic_plots.R`**: `mybar()`, `mydensity()`, `myline()`, `mylines()`
- **`distribution_plots.R`**: `cdf()`, `cdf_multi()`, box plot functions
- **`scatter_plots.R`**: All scatter/smooth variants, correlation plots
- **`composition_plots.R`**: `cmb()`, `acmb()`, `cml()`
- **`statistical_plots.R`**: `mybox()`, `mybox_df()`, `mybox_geomx()`, `vln_p()`

#### Task 2.3: Create `plot_themes.R`
- Common ggplot2 themes
- Color palettes
- Styling functions

### Phase 3: Analysis Module Refactoring

#### Task 3.1: Split pseudobulk_deg.R
Currently this file is very large. Split into:

- **`R/analysis/differential_expression/pseudobulk_core.R`**:
  - `prepare_pseudobulk()`
  - Data aggregation functions
  - Common utilities

- **`R/analysis/differential_expression/pseudobulk_edger.R`**:
  - `prepare_pseudobulk_edgeR()`
  - `run_pseudobulk_deg()`
  - edgeR-specific functions

- **`R/analysis/differential_expression/linear_models.R`**:
  - `pseudobulk_linear_fit()` (from legacy.R)
  - `post_hoc_slope_comparison()` (from legacy.R)
  - New signature: `linear_seurat()` (from signature.R)

#### Task 3.2: Consolidate Trajectory Analysis
Move from `pseudotime.R` to separate files:

- **`R/analysis/trajectory/slingshot.R`**: `run_slingshot_from_seurat()`
- **`R/analysis/trajectory/tradeseq.R`**: `analyze_gene_dynamics_tradeSeq()`
- **`R/analysis/trajectory/dynamics.R`**: `analyze_gene_dynamics()`, `process_gene_list_dynamics()`

#### Task 3.3: Refactor Marker Analysis
- **`R/analysis/markers/finding.R`**: Wrappers for FindMarkers
- **`R/analysis/markers/filtering.R`**: `marker_trim()`, `marker_filter()`, `lrf()`
- **`R/analysis/markers/synthesis.R`**: `synthesize_markers()`, `synthesize_ranks()`
- **`R/analysis/markers/utilities.R`**: `all_markers_to_list()`, `marker_print_all()`, `marker_print()`

#### Task 3.4: Module Scoring Consolidation
**Target**: `R/signatures/module_scoring.R`

Consolidate:
- `AddMultipleModuleScores()` from signature.R
- `PlotModuleScoreHeatmap()` from signature.R  
- `CompareModuleScoringMethods()` from signature.R
- Integration with heatmap visualization

### Phase 4: Efficiency Improvements

#### Task 4.1: Add Caching
- Implement caching for expensive operations
- Add cache invalidation logic
- Document caching behavior

#### Task 4.2: Parallelization
- Review all major functions for parallelization opportunities
- Use `future` or `parallel` consistently
- Add progress bars for long operations

#### Task 4.3: Input Validation
- Create validation helper functions
- Standardize error messages
- Add informative warnings

### Phase 5: Documentation & Testing

#### Task 5.1: Update Documentation
- Roxygen2 documentation for all functions
- Add examples for main workflows
- Create vignettes for common analyses

#### Task 5.2: Create Tests
- Unit tests for core functions
- Integration tests for pipelines
- Test coverage reports

## Implementation Strategy

### Step-by-Step Approach

1. **Week 1-2: Foundation**
   - Create new directory structure
   - Implement core and utility modules
   - Test basic functionality

2. **Week 3-4: Visualization**
   - Consolidate heatmap functions
   - Reorganize plot types
   - Create unified interfaces

3. **Week 5-6: Analysis Modules**
   - Split pseudobulk_deg.R
   - Refactor marker analysis
   - Consolidate module scoring

4. **Week 7: Polish & Documentation**
   - Update all documentation
   - Add examples
   - Create migration guide

5. **Week 8: Testing & Review**
   - Write tests
   - Code review
   - Performance benchmarks

### Breaking Changes & Migration

#### Deprecated Functions
Mark as deprecated (with `.Deprecated()`) but keep for backward compatibility:
- Old heatmap functions (`myhm_*`)
- Legacy functions in legacy.R
- Redundant plot functions

#### Migration Guide
Create `MIGRATION.md` documenting:
- Function name changes
- Parameter changes
- New recommended workflows
- Examples of old vs new code

## Success Metrics

1. **Code Quality**
   - ✓ No file > 800 lines
   - ✓ < 10% code duplication
   - ✓ All functions documented
   - ✓ Consistent naming conventions

2. **Performance**
   - ✓ Major functions 10-20% faster (through efficiency improvements)
   - ✓ Parallelization for computationally intensive tasks
   - ✓ Reduced memory footprint

3. **Usability**
   - ✓ Clear function interfaces
   - ✓ Consistent parameter names
   - ✓ Helpful error messages
   - ✓ Examples for all major functions

4. **Maintainability**
   - ✓ Clear file organization
   - ✓ Modular design
   - ✓ Easy to add new features
   - ✓ Test coverage > 70%

## Next Steps

1. Review and approve this plan
2. Start with Phase 1 (Foundation)
3. Regular check-ins and code reviews
4. Iterative implementation with testing
5. Final review and merge

## Questions for Discussion

1. Are there specific workflows that must be preserved exactly?
2. Which functions are most critical for backward compatibility?
3. Are there any performance bottlenecks to prioritize?
4. Preferred naming conventions?
5. Testing framework preferences?

---

**Document Version**: 1.0  
**Date**: 2025-10-10  
**Author**: Code Refactoring Initiative


