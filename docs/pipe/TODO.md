## Metadata & Quality (High Priority)
- [ ] **Robust Metadata Handling**: Refactor metadata joining to be more configuration-driven and robust to missing keys. Implement logic to handle duplicate columns by comparing values and keeping the most informative ones, or tagging conflicts.
- [ ] **QC Metrics**: Add more comprehensive QC metrics and visualizations.

## Configuration & Architecture (New)
- [ ] **Config Separation**: Split the current configuration into three distinct layers:
    1.  **Data Manifest** (e.g., `samples.csv`): Contains sample-specific paths (FASTQ, BAM, Cell Ranger output) and sample-level metadata (GEM ID, Patient ID). This replaces the current `config_complete.csv` role.
    2.  **Pipeline Config** (e.g., `pipeline_config.yaml`): Defines pipeline parameters (filtering thresholds, integration methods, resolution, scVI settings). This replaces `config_default.csv`.
    3.  **Execution Config** (e.g., `run_config.json`): Specifies run-specific settings (output directories, specific steps to run, resource allocation).
- [ ] **Dynamic Parameter Overrides**: Implement a mechanism to override pipeline config parameters via command-line arguments (e.g., `--set integration_method=scVI --set remove_doublets=TRUE`) without modifying the file.
- [ ] **Config Validation**: Add a validation step to ensure all required paths and parameters exist and are valid before starting the pipeline.

## Integration
- [ ] **scANVI**: Implement scANVI (Single-cell ANnotation using Variational Inference) in Step 6 for semi-supervised integration using cell type labels.
- [ ] **MrVI**: Explore MrVI (Multi-resolution Variational Inference) for multi-resolution integration.
- [ ] **Subset Integration**: Add functionality to integrate subsets of data (e.g., specific disease conditions) to verify batch correction without biological confounding.

## Preprocessing
- [ ] **Ambient RNA Removal**: Evaluate and implement newer/better ambient RNA removal tools beyond SoupX (e.g., CellBender).
- [ ] **Doublet Detection**: Integrate alternative doublet finders (e.g., DoubletFinder, Solo) or ensemble methods.
- [ ] **Demultiplexing**: Create a dedicated workflow/branch for demultiplexing (Step 0/1) to support various tools and improve robustness.

## General
- [ ] **Python Interoperability**: Improve `qs` <-> `h5ad` conversion to seamlessly leverage Python-based tools.

## Documentation & Directory structure updates.
- [ ] **Update Documentation**: Update documentation to reflect new configuration structure and workflow. 
.....

