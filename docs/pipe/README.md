# Pipeline Scripts & Documentation Index

This document serves as a comprehensive index to all pipeline-related scripts and documentation.

## Core Pipeline Scripts

The main pipeline execution scripts are located in `scripts/pipe/`:
- `pipe1_read_demulti.R` - Step 1: Data loading and demultiplexing (SNP/HTO)
- `pipe2_nmz_clustering.R` - Step 2: Normalization and clustering
- `pipe3_ambient_removal.R` - Step 3: Ambient RNA removal (SoupX)
- `pipe4_sctransform.R` - Step 4: SCTransform normalization
- `pipe5_doubletfinder.R` - Step 5: Doublet detection
- `pipe6_integration.R` - Step 6: Integration (RPCA/scVI)
- `pipe_validate.R` - Validation and QC script
- `start_pipe.R` - Lightweight startup script for package loading

## Utility Scripts

Utility scripts are located in `myR/R/`:
- `pipe_utils.R` - Configuration loading, logging, and path management
- `pipe_demulti.R` - Demultiplexing helper functions (SNP/HTO)
- `utils_demulti.R` - Low-level demultiplexing utilities

## Helper Scripts (Archive)

Historical helper scripts have been moved to `docs/pipe/archive/old_scripts/`:
- `prep_metadata.R` - Metadata preprocessing for HTO samples
- `fix_manifest_hto.R` - Manifest correction for HTO-only analysis
- `update_manifest.R` - General manifest update utility
- `update_config_columns.R` - Configuration column synchronization
- `verify_meta.R` - Metadata verification
- `create_downsampled_data.R` - Data downsampling for testing
- `vars_config.R` - Configuration variable management
- `load_renv.R` - renv loading utility

## Documentation

### Main Documentation (docs/pipe/)
- **[PIPE_INTEGRATED_GUIDE_KR.md](PIPE_INTEGRATED_GUIDE_KR.md)** - 파이프라인 통합 가이드 (국문, 메인)
- **[PIPE_INTEGRATED_GUIDE.md](PIPE_INTEGRATED_GUIDE.md)** - Pipeline Integrated Guide (English)
- **[CONTEXT.md](CONTEXT.md)** - Current pipeline status and known issues
- **[COMMANDS.md](COMMANDS.md)** - Common commands and usage examples
- **[DATA_FLOW.md](DATA_FLOW.md)** - Data flow and directory structure
- **[SCRIPT_STEPS_EXPLANATION.md](SCRIPT_STEPS_EXPLANATION.md)** - Detailed step-by-step explanations
- **[TODO.md](TODO.md)** - Open tasks and future improvements
- **[TESTING_LOG.md](TESTING_LOG.md)** - Testing history and results
- **[PIPELINE_REVIEW.md](PIPELINE_REVIEW.md)** - Code review and refactoring notes

### Quick Reference
- **[commands_little.md](commands_little.md)** - Quick command reference
- **[commands_from_10m.md](commands_from_10m.md)** - Commands from 10m conversation history
- **[quest_251204_1.md](quest_251204_1.md)** - Specific troubleshooting session (2024-12-04)

### Archive
- **[archive/](archive/)** - Historical logs, test results, and deprecated documents

## Configuration Examples

Example configuration files are in `config/`:
- `manifest_stroke.csv` - Main manifest for stroke dataset (SNP + HTO)
- `manifest_hto_fixed.csv` - HTO-only manifest with unique GEM names
- `manifest_hto_test.csv` - Minimal HTO test manifest
- `run_config_stroke.json` - Run configuration for stroke dataset
- `run_config_hto.json` - Run configuration for HTO-only analysis
- `meta_data_prep_251204.csv` - Pre-processed clinical metadata with HTO tags

## Execution Wrapper

- `scripts/pipe_wrapper.sh` - Shell wrapper for pipeline execution (if needed)

## Cross-References

- Related to FGS module: See `docs/fgs/FGS_INTEGRATED_GUIDE_KR.md`
- Related to pseudotime analysis: See `docs/pseudotime-dev/`
- Related to CCI analysis: See `docs/cci/`
- myR package documentation: See `myR/README_Korean.md`

---

**Last Updated**: 2025-12-09
