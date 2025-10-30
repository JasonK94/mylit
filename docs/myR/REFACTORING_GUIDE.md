# myR Package Refactoring Guide

## Overview

This document describes the refactoring of the myR package for improved organization, maintainability, and code reusability.

## New Directory Structure

```
myR/R/
├── core/                              # Core data processing
│   ├── validation.R                   # Input validation functions
│   ├── data_preparation.R             # Data extraction and preparation
│   └── seurat_conversion.R            # Seurat object conversions
│
├── utilities/                         # General utilities
│   ├── general_utils.R                # General helper functions
│   ├── sample_utils.R                 # Sample manipulation utilities
│   ├── plot_utils.R                   # Plot saving and styling utilities
│   └── demulti_utils.R                # Demultiplexing utilities
│
├── analysis/                          # Analysis modules
│   ├── markers/
│   │   └── marker_processing.R        # Marker gene processing and filtering
│   ├── differential_expression/
│   │   ├── pseudobulk_deg.R           # Pseudobulk DE analysis
│   │   └── linear_regression.R        # Linear regression for gene expression
│   ├── cell_communication/
│   │   └── nichenet_analysis.R        # NicheNet CCI analysis
│   ├── spatial/
│   │   └── geomx_analysis.R           # GeoMx spatial profiling
│   ├── trajectory/
│   │   └── trajectory_inference.R     # Slingshot and GAM analysis
│   ├── pathway/
│   │   └── pathway_enrichment.R       # GO, KEGG, GSEA analysis
│   └── composition/
│       └── cluster_frequency.R        # Cluster composition analysis
│
├── visualization/                     # Visualization modules
│   ├── basic_plots.R                  # Basic plotting functions
│   ├── composition_plots.R            # Composition and proportion plots
│   ├── heatmaps.R                     # Heatmap functions
│   └── scatter_plots.R                # Scatter and correlation plots
│
├── signatures/                        # Gene signature modules
│   ├── signature_scoring.R            # Module score calculation
│   └── signature_discovery.R          # Signature discovery methods
│
└── deprecated/                        # Legacy functions (to be removed)
    └── legacy.R                       # Old implementations

```

## File Migration Map

### From `CCI.R`
- `ligand_to_target()` → `analysis/cell_communication/nichenet_analysis.R`
- `run_nichenet_analysis()` → `analysis/cell_communication/nichenet_analysis.R`
- `prepare_nichenet_circos_data()` → `analysis/cell_communication/nichenet_analysis.R`
- `draw_nichenet_circos_plot()` → `analysis/cell_communication/nichenet_analysis.R`

### From `cluster_frequency.R`
- `seurat_group_stats()` → `analysis/composition/cluster_frequency.R`
- `seurat_posthoc_analysis()` → `analysis/composition/cluster_frequency.R`
- `extract_sample_metadata()` → `analysis/composition/cluster_frequency.R`
- `perform_prior_tests()` → `analysis/composition/cluster_frequency.R`
- `plot_cluster_fractions()` → `analysis/composition/cluster_frequency.R`
- `calculate_sample_signatures()` → `analysis/composition/cluster_frequency.R`
- `adjust_confounders()` → `analysis/composition/cluster_frequency.R`

### From `demulti_utils.R`
- `get_best_two()` → `utilities/demulti_utils.R`
- `get_barcode_mapping()` → `utilities/demulti_utils.R`
- `is_doublet()` → `utilities/demulti_utils.R`
- `generate_sample_values()` → `utilities/demulti_utils.R`
- `generate_sample_names()` → `utilities/demulti_utils.R`
- Additional helper functions added

### From `GeoMx.R`
- `prepare_geomx_data()` → `analysis/spatial/geomx_analysis.R`
- `q3_normalize()` → `analysis/spatial/geomx_analysis.R`
- `perform_qc()` → `analysis/spatial/geomx_analysis.R`
- `find_deg_geomx()` → `analysis/spatial/geomx_analysis.R`
- `plot_deg_volcano()` → `analysis/spatial/geomx_analysis.R`
- `plot_deg_heatmap()` → `analysis/spatial/geomx_analysis.R`
- `run_geomx_analysis()` → `analysis/spatial/geomx_analysis.R`

### From `markers.R`
- `marker_trim()` → `analysis/markers/marker_processing.R`
- `marker_filter()` → `analysis/markers/marker_processing.R`
- `lrf()` → `analysis/markers/marker_processing.R`
- `all_markers_to_list()` → `analysis/markers/marker_processing.R`
- `marker_print_all()` → `analysis/markers/marker_processing.R`
- `marker_print()` → `analysis/markers/marker_processing.R`
- `synthesize_markers()` → `analysis/markers/marker_processing.R`
- `synthesize_ranks()` → `analysis/markers/marker_processing.R`

### From `pathway_analysis.R`
- `convert_gene_ids()` → `analysis/pathway/pathway_enrichment.R`
- `prepare_gene_lists()` → `analysis/pathway/pathway_enrichment.R`
- `get_pathway_sets()` → `analysis/pathway/pathway_enrichment.R`
- `run_go_analysis()` → `analysis/pathway/pathway_enrichment.R`
- `run_kegg_analysis()` → `analysis/pathway/pathway_enrichment.R`
- `run_gsea_analysis()` → `analysis/pathway/pathway_enrichment.R`
- `format_results()` → `analysis/pathway/pathway_enrichment.R`
- `myGO()` → `analysis/pathway/pathway_enrichment.R`

### From `plots.R`
- `mybar()` → `visualization/basic_plots.R`
- `myline()` → `visualization/basic_plots.R`
- `mylines()` → `visualization/basic_plots.R`
- `mydensity()` → `visualization/basic_plots.R`
- `cdf()` → `visualization/basic_plots.R`
- `cdf_multi()` → `visualization/basic_plots.R`
- `cmb()` → `visualization/composition_plots.R`
- `acmb()` → `visualization/composition_plots.R`
- `cml()` → `visualization/composition_plots.R`
- `upset_gene_lists()` → `visualization/composition_plots.R`
- `vln_p()` → `visualization/composition_plots.R`
- `mybox()` → `visualization/scatter_plots.R` (to be created)
- `mybox_df()` → `visualization/scatter_plots.R` (to be created)
- `scatter_smooth_colored()` → `visualization/scatter_plots.R` (to be created)
- `myhm_genesets2()` → `visualization/heatmaps.R` (to be created)
- `myhm_genes2()` → `visualization/heatmaps.R` (to be created)

### From `pseudobulk_deg.R`
- `prepare_pseudobulk_edgeR()` → `analysis/differential_expression/pseudobulk_deg.R`
- `run_pseudobulk_deg()` → `analysis/differential_expression/pseudobulk_deg.R`

### From `pseudotime.R`
- `run_slingshot_from_seurat()` → `analysis/trajectory/trajectory_inference.R`
- `save_plot_with_conflict_resolution()` → `utilities/plot_utils.R`
- `analyze_gene_dynamics()` → `analysis/trajectory/trajectory_inference.R`
- `process_gene_list_dynamics()` → `analysis/trajectory/trajectory_inference.R`
- `analyze_gene_dynamics_tradeSeq()` → `analysis/trajectory/trajectory_inference.R`

### From `seurat_to_scanpy.R`
- `save_seurat_to_h5ad()` → `core/seurat_conversion.R` (to be created)

### From `signature.R`
- `AddMultipleModuleScores()` → `signatures/signature_scoring.R`
- `PlotModuleScoreHeatmap()` → `signatures/signature_visualization.R` (to be created)
- `CompareModuleScoringMethods()` → `signatures/signature_visualization.R` (to be created)
- `add_signature_enrichit()` → `signatures/signature_scoring.R`
- `add_progeny_scores()` → `signatures/signature_scoring.R`
- `linear_seurat()` → `analysis/differential_expression/linear_regression.R` (to be created)
- `plot_top_genes()` → `analysis/differential_expression/linear_regression.R` (to be created)
- **NEW:** `find_gene_signature()` → `signatures/signature_discovery.R`
- **NEW:** `score_signature()` → `signatures/signature_scoring.R`
- **NEW:** `print.gene_signature()` → `signatures/signature_discovery.R`

### From `test_claude.R`
- `downsample_sobj()` → `utilities/seurat_utils.R` (to be created)
- `printmy()` → `analysis/markers/marker_processing.R` (variant of marker_print)
- `printMy()` → `analysis/markers/marker_processing.R` (variant of marker_print_all)
- `sort_samples()` → `utilities/sample_utils.R`
- `print_gene_combinations()` → `utilities/gene_list_utils.R` (to be created)

### From `gene_list.R`
- All gene lists → Keep as-is, no changes needed
- These are data objects, not functions

### From `legacy.R`
- Legacy functions → `deprecated/legacy.R`
- Mark for deprecation, provide migration paths

## Key Improvements

### 1. **Modular Organization**
- Functions grouped by purpose and application
- Clear separation of concerns
- Easier to find and maintain functions

### 2. **Consistent Naming**
- Core functions use clear, descriptive names
- Internal helpers prefixed with `.` (e.g., `.run_edger_analysis`)
- Exported functions documented with `@export`

### 3. **Better Documentation**
- Each module has a header describing its purpose
- Functions have comprehensive `@param` and `@return` documentation
- Examples provided where appropriate

### 4. **Reduced Redundancy**
- Similar functions consolidated (e.g., marker processing)
- Common patterns extracted into utilities
- Duplicate implementations removed

### 5. **Enhanced Functionality**
- Validation functions ensure robust input handling
- Error messages are clear and actionable
- Optional parameters with sensible defaults

### 6. **Improved Testability**
- Smaller, focused functions easier to test
- Clear input/output contracts
- Modular design facilitates unit testing

## Migration Strategy

### Phase 1: Create New Structure (COMPLETE)
✅ Created core utilities
✅ Created modular analysis files
✅ Created visualization modules
✅ Created signature modules

### Phase 2: Update Exports and Documentation (IN PROGRESS)
- [ ] Update NAMESPACE with new exports
- [ ] Update DESCRIPTION with dependencies
- [ ] Create roxygen2 documentation
- [ ] Update man/ pages

### Phase 3: Deprecation (PLANNED)
- [ ] Mark old functions as deprecated
- [ ] Add `.Deprecated()` calls with new function names
- [ ] Update existing code to use new functions
- [ ] Add deprecation warnings to NAMESPACE

### Phase 4: Testing (PLANNED)
- [ ] Create unit tests for core functions
- [ ] Create integration tests for workflows
- [ ] Test backward compatibility
- [ ] Validate examples in documentation

### Phase 5: Cleanup (PLANNED)
- [ ] Remove deprecated functions
- [ ] Clean up old files
- [ ] Finalize documentation
- [ ] Update README with new structure

## Usage Examples

### Old Way
```r
# Old scattered function calls
source("myR/R/CCI.R")
source("myR/R/markers.R")
source("myR/R/plots.R")

markers <- FindAllMarkers(sobj)
markers_filtered <- marker_filter(marker_trim(markers))
nichenet_results <- run_nichenet_analysis(...)
```

### New Way
```r
# Clean, organized imports
library(myR)

# All functions available through package namespace
markers <- FindAllMarkers(sobj)
markers_filtered <- marker_filter(marker_trim(markers))
nichenet_results <- run_nichenet_analysis(...)

# Or with explicit namespacing
markers_filtered <- myR::marker_filter(myR::marker_trim(markers))
```

## Backward Compatibility

During the transition period:
1. Old function names will continue to work with deprecation warnings
2. Documentation will show both old and new usage
3. Migration helpers provided for complex changes
4. Gradual deprecation over multiple versions

## Contributing

When adding new functions:
1. Place in appropriate module based on functionality
2. Follow naming conventions
3. Add comprehensive documentation
4. Include examples
5. Write unit tests
6. Update this guide

## Questions or Issues?

Contact the package maintainer or open an issue on the repository.

---

**Last Updated:** 2025-10-10  
**Status:** Phase 1 Complete, Phase 2 In Progress


